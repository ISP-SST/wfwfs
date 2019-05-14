#ifndef WFWFS_SEEINGLOG_HPP
#define WFWFS_SEEINGLOG_HPP
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

#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>


namespace wfwfs {
    
    class SeeingLog : public std::ofstream {
    public:
        SeeingLog( std::string n, int c );
        SeeingLog( SeeingLog&& );
        
        void setInterval( int i ) { interval = i; };
        void setName( std::string n ) { name = n; };
        void setDir( std::string s ) { dir = s; };
        void open( void );
        void close( void );
        void reopen( void ) { close(); open(); };
        void check( void );

        
        std::string get_name( void ) { return name; };
        void set_min_lock( float ml ) { min_lock = ml; };

        bool hasColumns( void ) { return !columns.empty(); }
        void addColumns( std::string );
        void printHeader( void );
        
        void run( std::vector<DimmSet>& );
        void stop( void );
        
    private:
        std::string dir;
        std::string name;
        int interval;
        std::vector< std::pair<std::string,int> > columns;
        bool running;
        boost::posix_time::ptime timestamp;
        std::thread trd;
        std::mutex mtx;

    };
    
}   // wfwfs


#endif // WFWFS_SEEINGLOG_HPP
