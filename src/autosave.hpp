#ifndef WFWFS_AUTOSAVE_HPP
#define WFWFS_AUTOSAVE_HPP
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

#include "dimmset.hpp"

#include "tcpconnection.hpp"

#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>


namespace wfwfs {
    
    class AutoSave {
    public:
        AutoSave();
        AutoSave( AutoSave&& );
        
        void parsePropertyTree( const boost::property_tree::ptree&, const std::string& );
        void setDir( std::string s ) { dir = s; };
        
        void start( const std::vector<DimmSet>& );
        void stop( void );
        
    private:
        
        enum AS_TYPE { TP_NONE=0, TP_TIME, TP_TELNET, TP_DIMM };
        
        bool wait_for_time( void );
        bool wait_for_telnet( void );
        bool wait_for_dimm( const std::vector<DimmSet>& );
        
        // common
        std::string dir;
        std::string name;
        std::string acc_name;
        int frames_per_file, nframes;
        int as_type, delay, min_interval, poll_interval;
        size_t repeats;
        bool compress;

        std::atomic<bool> running;
        std::thread trd;
        
        // timer
        std::string time_string;
        
        // telnet
        TcpConnection::Ptr conn;
        std::string telnet_host, telnet_port;
        std::string telnet_cmd, telnet_reply;
        
        // DIMM
        std::string dimm_name;
        int dimm_duration;
        float dimm_limit;
        

    };
    
}   // wfwfs


#endif // WFWFS_AUTOSAVE_HPP
