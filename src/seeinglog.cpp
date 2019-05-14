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
#include "seeinglog.hpp"

#include "seeing.hpp"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>


using namespace wfwfs;
using namespace std;

namespace bpx = boost::posix_time;
namespace bfs = boost::filesystem;


SeeingLog::SeeingLog( void ) : ofstream(), name(), min_lock(-1.0), interval(1), running(false),
                                                           timestamp(bpx::not_a_date_time) {
                                                               
}


SeeingLog::SeeingLog( std::string n, int i ) : ofstream(), name(n), min_lock(0.0), interval(i), running(false),
                                                           timestamp(bpx::not_a_date_time) {
    
}


SeeingLog::SeeingLog( SeeingLog&& rhs ) : dir( std::move(rhs.dir) ), name( std::move(rhs.name) ), min_lock(std::move(rhs.min_lock)), 
                                          interval( std::move(rhs.interval) ), columns( std::move(rhs.columns) ),
                                          running(false), timestamp(std::move(rhs.timestamp)), mtx() {
            
}


void SeeingLog::parsePropertyTree( const boost::property_tree::ptree& cfg_ptree, const std::string& base_dir ) {
    
    dir = cfg_ptree.get<string>( "outputdir", base_dir+"/logs/" );
    name = cfg_ptree.get<string>( "name", "log" );
    interval = cfg_ptree.get<int>( "interval", 1 );
    min_lock = cfg_ptree.get<float>( "min_lock", 0.0 );
    string cols = cfg_ptree.get<string>( "columns", "" );
    if( !cols.empty() ) {
        addColumns( cols );
    }
    
}


void SeeingLog::open( void ) {
    
    timestamp = bpx::second_clock::universal_time();
    
    bfs::path filePath( dir );
    filePath /= name + "_" + to_iso_extended_string( timestamp.date() ) + ".log";
    
    bfs::path dirPath = filePath.parent_path();
    if( !dirPath.empty() && !bfs::exists(dirPath) ) {
        if( !bfs::create_directories(dirPath) ) {
            cerr << boost::format( "failed to create log directory: %s" ) % dirPath << endl;
            return;
        }
    }
    
    lock_guard<mutex> lock(mtx);
    if( ofstream::is_open() ) {
        ofstream::close();
    }
    
    ofstream::open( filePath.string(), std::ofstream::out | std::ofstream::app );
    
    if( !ofstream::tellp() ) {      // print header if we opened a new file
        printHeader();
    }
    
    bfs::path linkPath( dir );
    linkPath /= name + ".log";

    if( bfs::is_symlink(linkPath) ) {
        bfs::remove(linkPath);
    }
    
    bfs::create_symlink( bfs::relative(filePath, dirPath), linkPath );


}


void SeeingLog::close( void ) {
    
    lock_guard<mutex> lock(mtx);
    ofstream::close();
    
}


void SeeingLog::check( void ) {
    
    bpx::ptime now = bpx::second_clock::universal_time();  
    if( now.date() != timestamp.date() ) {
        reopen();
    } else if( !ofstream::is_open() ) {
        open();
    }
}


void SeeingLog::addColumns( std::string line ) {
    
    while( !line.empty() ) {
        string n = popword( line );
        if( !name.empty() ) {
            int c = interval;
            try {
                size_t pos = line.find_first_of(" ");
                c = boost::lexical_cast<int>( line.substr(0,pos) );
                popword( line );            // cast ok, so strip it from line
            } catch( ... ) { }
            if( c < 1 ) c = interval;
            columns.push_back( {n,c} );
        }
    }
  
}


void SeeingLog::printHeader( void ) {
    
    *this << "#  Seeing measurements" << endl;
    *this << "#  Date: " << to_iso_extended_string( timestamp.date() ) << endl;
    *this << "#  Log interval: " << interval << "s" << endl;
    *this << "#  Columns: timestamp";
    for( auto& c: columns ) *this << "  \"" << c.first << "\" (" << c.second << "s average)";
    *this << endl;
    
}


void SeeingLog::run( std::vector<DimmSet>& dimm_sets ) {
    
    using namespace std::chrono;

    static const bpx::ptime epoch_time( boost::gregorian::date(1970,1,1) ); 
    running = true;
    
    check();
    
    trd = std::thread([&](){
        
        high_resolution_clock::time_point next = high_resolution_clock::now();
        
        while( running ) {
            
            next += seconds(interval);
            
            bpx::ptime now = bpx::microsec_clock::universal_time();
            bool has_finite_values(false);
            vector<float> values( columns.size(), 0.0 );
            
            for( size_t i(0); i<columns.size(); ++i ) {
                for( auto& ds: dimm_sets ) {
                    if( ds.get_name() == columns[i].first ) {
                        PointF v = ds.calculate_r0( now, columns[i].second, min_lock );
                        if( v.max_abs() > 0.0 ) {
                            has_finite_values = true;
                            values[i] = (v.x+v.y)/2;
                        }
                    }
                }
            }
            if( has_finite_values && !now.is_not_a_date_time() ) {   
                bpx::time_duration since_epoch = now - epoch_time;
                size_t nSecs = since_epoch.total_seconds();
                size_t nMicros = now.time_of_day().total_microseconds() - now.time_of_day().total_seconds()*1000000L;
                string log_line = boost::str( boost::format("%ld.%06ld") % nSecs % nMicros );

                for( auto& v: values ) {
                    log_line += " " + to_string( v );
                }
                
                lock_guard<mutex> lock(mtx);
                if( !ofstream::is_open() ) {
                    running = false;
                    break;
                }

                *this << log_line << endl;
                ofstream::flush();
            }
            
            this_thread::sleep_until( next );
            
        }
    });

    
}


void SeeingLog::stop( void ) {
    
    running = false;
    if( trd.joinable() ) {
        trd.join();
    }

}

