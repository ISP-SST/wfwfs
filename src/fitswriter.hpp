#ifndef WFWFS_FITSWRITER_HPP
#define WFWFS_FITSWRITER_HPP
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

#include "filefits.hpp"
#include "frame.hpp"

#include <condition_variable>
// #include <cstdint>
#include <list>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace wfwfs {

    class FitsWriter {

    public:
        FitsWriter( FrameQueue&, const std::string&, int nThreads=1, bool compress=false );
        ~FitsWriter();

        void reset( int, int, bool );
        void makeHdr( void );
        
        void save( const std::string& filename, size_t& next_frame, int nframes );
        void wait( void );
        
        void save_meta( const std::vector<std::string>& m ) { extra_meta = m; };

        void updateCard( std::string key, std::string card ) { if( hdr ) Fits::updateCard( hdr->primaryHDU.cards, key, card ); };
        void updateCard( std::string card ) { if( hdr ) Fits::updateCard( hdr->primaryHDU.cards, card ); };
        
        static void clear_gmeta(void) { std::lock_guard<std::mutex> lock(globalMtx); globalMeta.clear(); };
        static void update_gmeta( std::string card ) { std::lock_guard<std::mutex> lock(globalMtx); Fits::updateCard( globalMeta, card ); };
        template <typename T>
        static void update_gmeta( std::string key, T value, std::string comment="" ) {
            std::lock_guard<std::mutex> lock(globalMtx); Fits::updateCard( globalMeta, key, Fits::makeCard(key, value, comment) );
        }

    private:
        
        static std::shared_ptr<uint8_t> get_buf( size_t N );
        static void return_buf( std::shared_ptr<uint8_t>& );
        static void clear_bufs( void );
        
        void open_file( const std::string& filename, size_t sz=0 );
        void close_file( void );
        void write_exptime_table( void );
        void write_acc( void );
       
        void thread_run(void);
        void push( const Frame& f );
        void push(void*, bpx::ptime);
        std::shared_ptr<uint8_t> pop();

        bool running;
        bool do_fsync;
        bool do_compress;
        bool do_write;
        int index;
        int nframes;
        int nthreads;
        int bytes_per_pixel;
        int npixels;
        int fd;

        size_t pcount, maxRowSize;
        size_t hdrEnd, dataStart;
        std::vector<std::string> extra_meta;
        std::shared_ptr<Fits> hdr;
        std::list<std::shared_ptr<uint8_t>> frames;
        std::list<boost::posix_time::ptime> times;
        std::vector<std::thread> threads;
        std::vector<int> offsets;
        std::mutex queueMtx;//, writeMtx;
        std::condition_variable cond;
        
        bool do_acc;
        std::string acc_filename;
        Array<uint32_t> acc;
        size_t frame_count;
        boost::posix_time::ptime first_ts, last_ts;
        

        static int activeCount, totalCount;
        static std::list<std::shared_ptr<uint8_t>> buffers;
        static std::vector<std::string> globalMeta;
        static std::mutex globalMtx, bufMtx;
        
        FrameQueue& fq;
        
    };

}   // wfwfs


#endif      // WFWFS_FITSWRITER_HPP
