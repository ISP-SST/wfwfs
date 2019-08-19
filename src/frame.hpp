#ifndef WFWFS_FRAME_HPP
#define WFWFS_FRAME_HPP
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

#include "cell.hpp"
#include "functions.hpp"

#include <atomic>
#include <condition_variable>
#include <memory>
#include <map>
#include <mutex>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace wfwfs {

    struct Frame {
        
        Frame() : data(nullptr), id(0), offset(0), uc(0), hist(nullptr), stats(1,{0}) {}
        Frame( Frame&& f ) : timestamp(f.timestamp), data(f.data), id(f.id), offset(f.offset),
            uc(f.uc.load()), hist(f.hist), stats(f.stats) {}
        
        boost::posix_time::ptime timestamp;
        uint8_t* data;
        size_t id;
        size_t offset;
        std::atomic<int> uc;    // Usage count, for preventing data-mangling.
        
        // stats
        uint32_t* hist;
        struct Stats {
            int min, max;
            double avg, rms;
        };
        std::vector<Stats> stats;

    };

    struct LockedFrame {
        LockedFrame( Frame& f ) : frame(f) { frame.uc++; };
        LockedFrame( const LockedFrame& lf ) : frame(lf.frame) { frame.uc++; };
        LockedFrame( LockedFrame&& ) = delete;
        ~LockedFrame() { frame.uc--; };
        Frame& frame;
    };
    
    struct FrameQueue {
        
        size_t frameSize;
        size_t blockSize;
        size_t counter;
        size_t width;
        size_t height;
        size_t depth;
        size_t currentFrame;
        
        std::unique_ptr<uint8_t[]> data;
        std::unique_ptr<uint32_t[]> hist;
        std::vector<Frame> frames;
        std::map<size_t, Frame&> frames_by_id;
        std::mutex mtx;
        std::condition_variable cond;
        
        FrameQueue();
        void resize( size_t w, size_t h, size_t nF, size_t d=8 );
        PointI get_size( void ) const { return PointI( height, width ); };
        
        Frame& getEmpty( void );
        void queue( Frame& );
        
        size_t getNearestID( const boost::posix_time::ptime& timestamp );
        Frame& getFrame( size_t id, bool wait=true );
        void dropFrame( size_t id );
        
        void resetCounter( void );
        void getStats( Frame& f, size_t stride=1 );
        void getStats( Frame& f, const std::vector<Cell>& );
        
        template <typename T>
        void addFrame( const Frame& f, T* accumulator ) {
            size_t nEl = width*height;
            if( depth > 8 ) {
                const uint16_t* data = reinterpret_cast<const uint16_t*>( f.data );
                std::transform( accumulator, accumulator+nEl, data, accumulator, add<T,uint16_t>() );
            } else {
                std::transform( accumulator, accumulator+nEl, f.data, accumulator, add<T,uint8_t>() );
            } 

        }
        
        template <typename T,typename U>
        void addFrame( const T* d, U* accumulator ) {
            size_t nEl = width*height;
            if( depth > 8 ) {
                const uint16_t* data = reinterpret_cast<const uint16_t*>( d );
                std::transform( accumulator, accumulator+nEl, data, accumulator, add<U,uint16_t>() );
            } else {
                const uint8_t* data = reinterpret_cast<const uint8_t*>( d );
                std::transform( accumulator, accumulator+nEl, data, accumulator, add<U,uint8_t>() );
            } 
        }

        
        template <typename T>
        void copyFrame( const Frame& f, T* out ) {
            size_t nEl = width*height;
            if( depth > 8 ) {
                const uint16_t* data = reinterpret_cast<const uint16_t*>( f.data );
                std::copy_n( data, nEl, out );
            } else {
                std::copy_n( f.data, nEl, out );
            } 

        }
        
        uint8_t* getData(void) { return data.get(); }
        uint8_t* getEnd(void) { return data.get()+blockSize; }
        
    };
    
}   // wfwfs

#endif // WFWFS_FRAME_HPP
