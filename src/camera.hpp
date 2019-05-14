#ifndef WFWFS_CAMERA_HPP
#define WFWFS_CAMERA_HPP
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

#include <condition_variable>
#include <mutex>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>

namespace wfwfs {

    class Camera {
        
    public:
        
        typedef std::function<void(const uint8_t*, boost::posix_time::ptime)> frame_cb_t;
        
        struct Config {
            int depth;
            int width;
            int height;
            int n_buf;                 // size of the camera-side buffer
        } cfg;

        struct Info {
            std::string id;
            std::string model;
            std::string firmware_version;
            std::string serial;
        } info;

        struct State {
            double exposure;
            double gain;
            double interval;
            double offset;
            double temp;
        } state;

        explicit Camera( boost::property_tree::ptree& );
        virtual ~Camera( void );

        virtual double get_exposure( void ) const { std::lock_guard<std::mutex> lock(mtx); return state.exposure; };
        virtual double get_gain( void ) const { std::lock_guard<std::mutex> lock(mtx); return state.gain; };
        virtual double get_interval( void ) const { std::lock_guard<std::mutex> lock(mtx); return state.interval; };
        virtual double get_offset( void ) const { std::lock_guard<std::mutex> lock(mtx); return state.offset; };
        virtual double get_temp( void ) const { std::lock_guard<std::mutex> lock(mtx); return state.temp; };
        virtual PointI get_size( void ) const { std::lock_guard<std::mutex> lock(mtx); return PointI( cfg.height, cfg.width ); };

        std::vector<std::string> get_messages( void ) {
            std::lock_guard<std::mutex> lock(mtx);
            std::vector<std::string> ret;
            std::swap( msgList, ret );		// swap to also clear the message list
            return ret;
        }
        
        virtual void set_exposure( double e ) { std::lock_guard<std::mutex> lock(mtx); state.exposure = e; };
        virtual void set_gain( double g ) { std::lock_guard<std::mutex> lock(mtx); state.gain = g; };
        virtual void set_interval( double i ) { std::lock_guard<std::mutex> lock(mtx); state.interval = i; };
        virtual void set_offset( double o ) { std::lock_guard<std::mutex> lock(mtx); state.offset = o; };
        
        void set_queue_callback( frame_cb_t cb ) { std::lock_guard<std::mutex> lock(mtx); qcb = cb; };

        // allow camera-specific parsing of incoming telnet commands 
        virtual bool parse_cmd( const std::string& cmd, const std::string& what, std::string& line, std::string& reply );

        virtual void init( void ) {};
        virtual void cleanup( void ) {};

        virtual void run( frame_cb_t cb ) = 0;
        virtual void stop( void ) {};
        
        virtual void copyRawData( const uint8_t* in, uint8_t* out, uint8_t outDepth ) const {};
        
    protected:
        
        mutable std::mutex mtx;
        std::condition_variable cond;
        mutable std::vector<std::string> msgList;
        boost::property_tree::ptree& settings;
        bool running;
        frame_cb_t qcb;
        
    };

}   // wfwfs

#endif // WFWFS_CAMERA_HPP
