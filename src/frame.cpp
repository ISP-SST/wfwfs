/***************************************************************************
 *                                                                         *
 *   WFWFS (Wide-field wavefront sensor) software                          *
 *                                                                         *
 *   Copyright (C) 2018-2019 by Tomas Hillberg  <hillberg@astro.su.se>     *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 *                                                                         *
 ***************************************************************************/
#include "frame.hpp"

#include <iostream>
#include <math.h>

using namespace wfwfs;
using namespace std;

namespace bpx = boost::posix_time;

FrameQueue::FrameQueue() : frameSize(0), blockSize(0), counter(1), width(0), height(0), depth(0) {
    
}


void FrameQueue::resize( size_t w, size_t h, size_t nF, size_t bitdepth ) {
    
    lock_guard<mutex> lock( mtx );
    
    size_t elementSize = (bitdepth>8)?2:1;
    size_t newFrameSize = elementSize*w*h;
    size_t newBlockSize = newFrameSize*nF;
    
    width = w;
    height = h;
    depth = bitdepth;
    
    if( (newFrameSize == frameSize) && (newBlockSize == blockSize) ) return;
    
    next_ = current_ = frames.end();
    
    if( newBlockSize == 0 ) {
        data.reset();
        hist.reset();
        frames.clear();
        return;
    }

    frameSize = newFrameSize;
    blockSize = newBlockSize;

    size_t histMax = (1<<depth);
    hist.reset( new uint32_t[ nF*histMax ] );
    
    data.reset( new uint8_t[blockSize] );
    frames.resize( nF );
    uint8_t* ptr = data.get();
    uint32_t* histPtr = hist.get();
    size_t o(0);
    for( auto &f: frames ) {
        f.data = ptr;
        ptr += frameSize;
        f.hist = histPtr;
        histPtr += histMax;
        f.offset = o++;
    }
    
    next_ = frames.begin();
    
}


Frame& FrameQueue::getEmpty( void ) {
    
    lock_guard<mutex> lock( mtx );
    
    if( frames.empty() || next_ == frames.end() ) {
        throw runtime_error("The FrameQueue is empty, can not get Frame!");
    }
    
    current_ = next_++;
    if( next_ == frames.end() ) {
        next_ = frames.begin();
    }

    frames_by_id.erase( current_->id );
    current_->id = counter++;
    
    return *current_;
    
}


void FrameQueue::queue( Frame& f ) {
    
    lock_guard<mutex> lock( mtx );
    frames_by_id.emplace( f.id, f );
    
    cond.notify_all();

}


size_t FrameQueue::getNearestID( const bpx::ptime& timestamp ) {
    
    size_t nearestID(0);                                // 0 will translate to getting latest frame in the queue.
    bpx::time_duration smallest_difference(24,0,0,0);   // initialize to something large (24-hours)
    
    unique_lock<mutex> lock( mtx );
    
    for( const auto& fp: frames_by_id ) {
        bpx::time_duration diff( timestamp - fp.second.timestamp );
        if( diff.is_negative() ) diff.invert_sign();
        if( nearestID && (diff > smallest_difference) ) { // we passed the minimum since frames_by_id is sequential
            break;
        } else if( diff < smallest_difference ) {
            smallest_difference = diff;
            nearestID = fp.first;
        }
    }
    
    return nearestID;
    
}
    

Frame& FrameQueue::getFrame( size_t id, bool wait ) {

    unique_lock<mutex> lock( mtx );

    if( id == 0 ) { // get latest, wait only if queue is empty
        while( wait && frames_by_id.empty() ) {
            cond.wait( lock );
        }
        if( frames_by_id.empty() ) {
            throw runtime_error("The FrameQueue is empty, can not get Frame!");
        }
        return frames_by_id.rbegin()->second;
    }
    
    // if we get here: id > 0
    while( wait && (frames_by_id.empty() || id > frames_by_id.rbegin()->second.id ) ) { // wait for future id
        cond.wait( lock );
    }
    
    if( frames_by_id.count(id) ) {
        return frames_by_id.at(id);
    } else if( !frames_by_id.empty() ) {   // if frame with id has already been overwritten, just return the latest and let the caller decide what to do
        return frames_by_id.rbegin()->second;
    }
    
    throw runtime_error("The FrameQueue is empty, can not get Frame!");

    
}


void FrameQueue::dropFrame( size_t id ) {
    
    lock_guard<mutex> lock( mtx );
    frames_by_id.erase( id );
    
    cond.notify_all();
    
}


void FrameQueue::resetCounter( void ) {
    lock_guard<mutex> lock( mtx );
    counter = 1;
}


void FrameQueue::getStats( Frame& frame, size_t stride ) {

    uint32_t* histPtr = frame.hist;
    if( !histPtr || !stride ) return;
    
    size_t histMax = (1<<depth);
    
    std::fill_n( histPtr, histMax, 0 );
    //memset(frame.df_hist, 0, maxval * sizeof *frame.df_hist);
    int npixels = 0;

    size_t y1 = height/8;
    size_t y2 = 7*(height/8);
    if( depth > 8 ) {
        uint16_t* data = reinterpret_cast<uint16_t*>( frame.data );
        uint16_t bitmask = histMax-1;
        for( size_t y=y1; y<y2; y += stride ) {
            size_t i = y * width;
            size_t o1 = i + width/8;
            size_t o2 = i + 7*(width/8);
            for( size_t o=o1; o<o2; ++o ) {
                frame.hist[ data[o]&bitmask ]++;
                //uint16_t c = df_correct(image, o);
                //frame.df_hist[c]++;
                npixels++;
            }
        }
    } else {
        uint8_t bitmask = histMax-1;
        for( size_t y=y1; y<y2; y += stride ) {
            size_t i = y * width;
            size_t o1 = i + width/8;
            size_t o2 = i + 7*(width/8);
            for( size_t o=o1; o<o2; ++o ) {
                frame.hist[ frame.data[o] & bitmask]++;
                //uint8_t c = df_correct(frame.data, o);
                //frame.df_hist[c]++;
                npixels++;
            }
        }
    } 

    double sum = 0;
    double sumsquared = 0;
    int min = -1;
    int max = -1;
    //int df_max = -1;
    //int df_min = -1;

    for( size_t i = 0; i<histMax; i++ ) {
        double di = i;
        uint32_t val = frame.hist[i];
        if( val ) {
            sum += di * frame.hist[i];
            sumsquared += di * di * frame.hist[i];
            if( min < 0 ) {
                min = i;
            }
            max = i;
        }
//         if( frame.df_hist[i] ) {
//             if(df_min < 0)
//                 df_min = i;
//             df_max = i;
//         }
    }

    sum /= npixels;
    sumsquared /= npixels;

    frame.stats[0].avg = sum;
    frame.stats[0].rms = sqrt(sumsquared - sum * sum) / sum;
    frame.stats[0].min = min;
    frame.stats[0].max = max;
    //frame.df_max = df_max;
    //frame.df_min = df_min;
    
}


void FrameQueue::getStats( Frame& frame, const std::vector<Cell>& cells ){

    uint32_t* histPtr = frame.hist;
    if( !histPtr ) return;
    
    frame.stats.resize( cells.size() + 1 );
    
    size_t histMax = (1<<depth);
    
    std::fill_n( histPtr, histMax, 0 );
    //memset(frame.df_hist, 0, maxval * sizeof *frame.df_hist);
    int npixels = 0;

    size_t y1 = height/8;
    size_t y2 = 7*(height/8);
    if( depth > 8 ) {
        uint16_t* data = reinterpret_cast<uint16_t*>( frame.data );
        uint16_t bitmask = histMax-1;
        for( size_t y=y1; y<y2; ++y ) {
            size_t i = y * width;
            size_t o1 = i + width/8;
            size_t o2 = i + 7*(width/8);
            for( size_t o=o1; o<o2; ++o ) {
                frame.hist[ data[o]&bitmask ]++;
                //uint16_t c = df_correct(image, o);
                //frame.df_hist[c]++;
                npixels++;
            }
        }
    } else {
        uint8_t bitmask = histMax-1;
        for( size_t y=y1; y<y2; ++y ) {
            size_t i = y * width;
            size_t o1 = i + width/8;
            size_t o2 = i + 7*(width/8);
            for( size_t o=o1; o<o2; ++o ) {
                frame.hist[ frame.data[o] & bitmask]++;
                //uint8_t c = df_correct(frame.data, o);
                //frame.df_hist[c]++;
                npixels++;
            }
        }
    } 

    double sum = 0;
    double sumsquared = 0;
    int min = -1;
    int max = -1;
    //int df_max = -1;
    //int df_min = -1;

    for( size_t i = 0; i<histMax; i++ ) {
        double di = i;
        uint32_t val = frame.hist[i];
        if( val ) {
            sum += di * frame.hist[i];
            sumsquared += di * di * frame.hist[i];
            if( min < 0 ) {
                min = i;
            }
            max = i;
        }
//         if( frame.df_hist[i] ) {
//             if(df_min < 0)
//                 df_min = i;
//             df_max = i;
//         }
    }

    sum /= npixels;
    sumsquared /= npixels;

    frame.stats[0].avg = sum;
    frame.stats[0].rms = sqrt(sumsquared - sum * sum) / sum;
    frame.stats[0].min = min;
    frame.stats[0].max = max;
    //frame.df_max = df_max;
    //frame.df_min = df_min;
    
}
