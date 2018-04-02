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
#include "filecam.hpp"
#include "filefits.hpp"

#include <chrono>
#include <memory>
#include <thread>

#include <boost/range/iterator_range.hpp>
#include <boost/thread/thread.hpp>

using namespace wfwfs;
using namespace std;

namespace {
    
    struct FileData {
        FileData() : m_nFrames(0), m_Current(0), m_stride(0), m_nElements(0), exposure(0), gain(0), interval(0), maxval(0) {}
        const uint8_t* get_next( bpx::ptime& time ) const {
            if( !m_data ) return nullptr;
            if( m_Current >= m_nFrames ) {
                m_Current = 0;
                //exit(0);
            }
            if( false && m_Current < frameTimes.size() ) time = frameTimes[m_Current];
            else time = bpx::microsec_clock::universal_time();
            return m_data.get() + m_Current++ * m_stride;
        }
        void load( string fn, shared_ptr<FileMeta>& meta ) {
            if( !meta ) meta = getMeta( fn, true );
            size_t dS = meta->dataSize();
            if( dS ) {
                shared_ptr<Fits> hdr = static_pointer_cast<Fits>(meta);
                if( hdr ) {
                    maxval = Fits::getValue<int>( hdr->primaryHDU.cards, "DATAMAX" );
                    gain = Fits::getValue<float>( hdr->primaryHDU.cards, "DETGAIN" );
                    interval = Fits::getValue<float>( hdr->primaryHDU.cards, "CADENCE" );
                    exposure = Fits::getValue<float>( hdr->primaryHDU.cards, "TEXPOSUR" );
                    if( exposure == 0.0 ) {
                        exposure = Fits::getValue<float>( hdr->primaryHDU.cards, "XPOSURE" );
                    }
                }
                frameTimes =  hdr->getStartTimes();
                m_nFrames = meta->getNumberOfFrames();
                m_nElements = meta->nElements();
                m_stride = dS/m_nFrames;
                m_data.reset( new uint8_t[dS] );
                readFile( fn, reinterpret_cast<char*>(m_data.get()), meta );
            } else {
                m_data.reset();
                m_nFrames = m_stride = 0;
            }
            m_Current = 0;
        }
        unique_ptr<uint8_t[]> m_data;
        std::vector<bpx::ptime> frameTimes;
        size_t m_nFrames;
        mutable size_t m_Current;
        size_t m_stride;
        size_t m_nElements;
        float exposure, gain, interval;
        int maxval;
    } images;
    
}

FileCam::FileCam( boost::property_tree::ptree& pt ) : Camera(pt) {

    string path = pt.get<string>( "camera_type", "" );
    set_path( path );
    
    // FIXME: only using a single file for now
    if( files.size() ) {
        auto it = files.begin();
        images.load( it->first, it->second );
        size_t w_ind = 0;
        size_t h_ind = 1;
        if( it->second->nDims() > 2 ) {
            w_ind++;
            h_ind++;
        }
        
        if( images.maxval ) {
            int maxval = images.maxval;
            cfg.depth = 1;
            while(maxval>>=1) cfg.depth++;
        } else {
            cfg.depth = it->second->elementSize()*8;
        }

        cfg.width = it->second->dimSize(w_ind);
        cfg.height = it->second->dimSize(h_ind);
        cfg.n_buf = images.m_nFrames;
        
        state.exposure = images.exposure;
        state.gain = images.gain;
//        state.interval = images.interval;
        
    }
    
    info.model = "FITS camera";
    info.firmware_version = path;
    
//     cout << "Using input from: " << path << endl;
//     cout << "   depth: " << cfg.depth << "  width: " << cfg.width << "  height: " << cfg.height << "  n_buf: " << cfg.n_buf << endl;
    
}


FileCam::~FileCam( void ) {

}


void FileCam::add_file( const bfs::path& p ) {
    
    if( bfs::is_regular_file( p ) ) {
        string fn = p.string();
        try {
            shared_ptr<FileMeta> meta = getMeta( fn, true );
            if( meta ) {
                files[fn] = meta;
            }
        } catch ( ... ) {
            // ignore
        }
    }

}


void FileCam::set_path( const bfs::path& p ) {
    
    if( bfs::is_regular_file( p ) ) {
        add_file( p );
    } else if( bfs::is_directory( p )) {
        for( auto& entry : boost::make_iterator_range( bfs::directory_iterator(p), {} ) ) {
            set_path( entry );
        }
    }

}


void FileCam::init( void ) {


}


void FileCam::cleanup( void ) {
    
}


void FileCam::run( frame_cb_t cb ) {
    
    using namespace std::chrono;

    running = true;
    if( cb ) qcb = cb;

    system_clock::time_point next_frame_at = system_clock::now();
    
    while( running ) {
        boost::this_thread::interruption_point();
        this_thread::sleep_until( next_frame_at );
        bpx::ptime ts;
        const uint8_t* f = images.get_next( ts );
        next_frame_at += duration<long,std::micro>( static_cast<long>(state.interval*1E6) );
        lock_guard<mutex> lock(mtx);
        if( qcb ) {
            qcb(f, ts);
        }
    }

}


void FileCam::stop( void ) {
    
    running = false;
    
}


void FileCam::copyRawData( const uint8_t* in, uint8_t* out, uint8_t outDepth ) const {
    
    int shift = 0;
    size_t nEl = images.m_stride;
    if( cfg.depth > 8 ) nEl /= 2;
    if( cfg.depth == outDepth ) {
        memcpy( out, in, images.m_stride );
    } else if( outDepth > cfg.depth ) {
        shift = (outDepth - cfg.depth);
        if( cfg.depth > 8 ) {
            const uint16_t* in16 = reinterpret_cast<const uint16_t*>(in);
            std::transform( in16, in16+nEl, reinterpret_cast<uint16_t*>(out), [shift]( uint16_t a ) { return a<<shift; } );
        } else {
            if( outDepth > 8 ) {
                std::transform( in, in+nEl, reinterpret_cast<uint16_t*>(out), [shift]( uint8_t a ) { return static_cast<uint16_t>(a)<<shift; } );
            } else {
                std::transform( in, in+nEl, out, [shift]( uint8_t a ) { return a<<shift; } );
            }
        }
    } else {
        shift = (cfg.depth - outDepth);
        if( cfg.depth > 8 ) {
            const uint16_t* in16 = reinterpret_cast<const uint16_t*>(in);
            if( outDepth > 8 ) {
                std::transform( in16, in16+nEl, reinterpret_cast<uint16_t*>(out), [shift]( uint16_t a ) { return a>>shift; } );
            } else {
                std::transform( in16, in16+nEl, out, [shift]( uint16_t a ) { return a>>shift; } );
            }
        } else {
            std::transform( in, in+nEl, out, [shift]( uint8_t a ) { return a>>shift; } );
        }
    }

}
