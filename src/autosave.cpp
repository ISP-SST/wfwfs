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
#include "autosave.hpp"

#include "daemon.hpp"
#include "seeing.hpp"
#include "translators.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>


using boost::algorithm::iequals;

using namespace wfwfs;
using namespace std;
using namespace std::chrono;

namespace bpx = boost::posix_time;
namespace bfs = boost::filesystem;

namespace {
    
    system_clock::time_point get_next_time( string time_str ) {
        typedef duration<int,std::ratio<60*60*24>> days_type;
        system_clock::time_point ret = system_clock::now();
        if( !time_str.empty() ) {
            if( time_str[0] == '+' ) {
                time_str = time_str.substr(1);
            } else {
                ret = time_point_cast<days_type>(ret);
            }
            int tmp = pop<int>( time_str, ":. " );
            if( tmp ) ret += hours(tmp);
            tmp = pop<int>( time_str, ":. " );
            if( tmp ) ret += minutes(tmp);
            tmp = pop<int>( time_str, ":. " );
            if( tmp ) ret += seconds(tmp);
        }
        if( ret < system_clock::now() ) {
            ret += days_type(1);
        }
        return ret;
    }
    
}


AutoSave::AutoSave() : frames_per_file(20), nframes(1), as_type(TP_NONE), delay(0), min_interval(0), poll_interval(5), repeats(-1),
                       compress(false), running(false), dimm_duration(10), max_consecutive(0), dimm_limit(0.1), consecutive_limit(0), trailing(false) {
    
}


AutoSave::AutoSave(AutoSave&& rhs) : dir(rhs.dir), name(rhs.name), acc_name(rhs.acc_name),
        frames_per_file(rhs.frames_per_file), nframes(rhs.nframes), as_type(rhs.as_type), delay(rhs.delay),
        min_interval(rhs.min_interval), poll_interval(rhs.poll_interval), repeats(rhs.repeats), compress(rhs.compress),
        running(false), time_string(rhs.time_string),
        telnet_host(rhs.telnet_host), telnet_port(rhs.telnet_port), telnet_cmd(rhs.telnet_cmd), telnet_reply(rhs.telnet_reply),
        dimm_name(rhs.dimm_name), dimm_duration(rhs.dimm_duration), max_consecutive(rhs.max_consecutive), dimm_limit(rhs.dimm_limit), consecutive_limit(rhs.consecutive_limit), trailing(rhs.trailing) {


}


void AutoSave::parsePropertyTree( const boost::property_tree::ptree& cfg_ptree, const std::string& base_dir ) {
    
    bfs::path tmpDir( cfg_ptree.get<string>( "outputdir", base_dir ) );
    string subdir = cfg_ptree.get<string>( "subdir", "" );
    if( !subdir.empty() ) {
        tmpDir /= subdir;
    }
    dir = tmpDir.string();
    name = cfg_ptree.get<string>( "name", "" );
    acc_name = cfg_ptree.get<string>( "acc_name", "" );
    string line = cfg_ptree.get<string>( "trigger", "" );
    string trigger_type = popword( line );
    if( iequals( trigger_type, "TIME" ) ) {
        as_type = TP_TIME;
        time_string = popword( line );
    } else if( iequals( trigger_type, "TELNET" ) ) {
        as_type = TP_TELNET;
        telnet_host = popword( line );
        telnet_port = popword( line );
        telnet_cmd = cfg_ptree.get<string>( "telnet_cmd", "" );
        telnet_reply = cfg_ptree.get<string>( "telnet_reply", "" );
    } else {    // trigger is name of a DIMM
        as_type = TP_DIMM;
        dimm_name = trigger_type;
        dimm_duration = max(1,pop<int>(line));
        dimm_limit = cfg_ptree.get<float>( "limit", dimm_limit );
        consecutive_limit = cfg_ptree.get<float>( "consecutive_limit", consecutive_limit );
        max_consecutive = cfg_ptree.get<int>( "max_consecutive", max_consecutive );
        trailing = cfg_ptree.get<bool>( "trailing", trailing );
    }
    
    poll_interval = cfg_ptree.get<int>( "poll_interval", 5 );
    frames_per_file = cfg_ptree.get<int>( "frames_per_file", frames_per_file );
    nframes = cfg_ptree.get<int>( "nframes", nframes );
    repeats = cfg_ptree.get<int>( "repeat", -1 );
    compress = cfg_ptree.get<int>( "compress", compress );
    delay = cfg_ptree.get<int>( "delay", delay );
    min_interval = cfg_ptree.get<int>( "min_interval", min_interval );

}



void AutoSave::start( const vector<DimmSet>& dimm_sets ) {
    
    if( running.exchange(true) ) {
        return;
    }

    static const bpx::ptime epoch_time( boost::gregorian::date(1970,1,1) );
    
    switch( as_type ) {
        case TP_TIME: ;
        case TP_TELNET: ;
        case TP_DIMM: break;
        default: stop(); return;
    }
    
    trd = std::thread([&](){
        size_t cnt(0);
        bool go(false);
        int consecutive_count(0);
        size_t first_frame(0);
        while( running && (cnt++ < repeats) ) {
            
            float this_climit(0.0);
            
            if( !go ) {
                consecutive_count = 0;
                first_frame = 0;
                switch( as_type ) {
                    case TP_TIME: go = wait_for_time(); break;
                    case TP_TELNET: go = wait_for_telnet(); break;
                    case TP_DIMM: go = wait_for_dimm(dimm_sets); break;
                    default: ;
                }
            } else {
                consecutive_count++;
            }
            
            if( running && go ) {
                Daemon& d = Daemon::get();
                high_resolution_clock::time_point next = high_resolution_clock::now();      // Save the timestamp when burst was started
                      
                if( first_frame == 0 ) {
                    if( delay > 0 ) {
                        this_thread::sleep_for( seconds(delay) );
                    } else if( delay < 0 ) {
                        bpx::ptime first_frametime = bpx::microsec_clock::universal_time() - bpx::time_duration( 0, 0, abs(delay), 0 );
                        first_frame = d.getNearestID( first_frametime );
                    }
                }
                
                if( as_type == TP_DIMM ) {
                    this_climit = consecutive_limit*get_r0(dimm_sets);
                }
                d.save_fits(  d.make_filename(name,-1,-1), nframes, frames_per_file,
                              compress, first_frame, d.make_filename(acc_name,-1,-1) );
                
                first_frame += nframes;
                go = false;
                
                if( (as_type == TP_DIMM) && !min_interval && (consecutive_count < max_consecutive) && (this_climit > 0.0) ) {
                    this_thread::sleep_for( seconds(dimm_duration) );       // Wait for a full period before checking r0
                    float r0 = get_r0(dimm_sets);
                    if( r0 >= this_climit ) {
                        go = true;
                        continue;
                    }
                }
                
                if( running && min_interval ) {
                    next += seconds(min_interval);              // sleep until min_interval has passed 
                    this_thread::sleep_until( next );
                }

            } else {
                usleep(100000);     // just to avoid a busy-loop.
            }
            
        }
    });

    
}


void AutoSave::stop( void ) {
    
    if( running.exchange(false) ) {
        if( trd.joinable() ) {
            trd.join();
        }
    }

}


bool AutoSave::wait_for_time( void ) {
    
    system_clock::time_point next = get_next_time( time_string );
    system_clock::time_point now = system_clock::now();
    
    while( now < next ) {
        this_thread::sleep_for( seconds(poll_interval) );
        if( !running ) return false;
        now = system_clock::now();
    }
    return true;
    
}


bool AutoSave::wait_for_telnet( void ) {
    
    bool reply_ok(false);
    Daemon& d = Daemon::get();

    while( running && !reply_ok ) {
        d.connect( conn, telnet_host, telnet_port );
        if( !conn ) break;
        conn->writeline( telnet_cmd );
        string reply;
        size_t n = conn->readline( reply );
        if( !running || (n == 0) ) break;
        if( iequals( reply, telnet_reply ) ) {
            reply_ok = true;
            break;
        }
        this_thread::sleep_for( seconds(poll_interval) );
    }
    
    return reply_ok;
}


bool AutoSave::wait_for_dimm( const vector<DimmSet>& dimm_sets ) {

    float previous_r0;
    bool reply_ok(false);
    bool r0_reached(false);
    for( auto& ds: dimm_sets ) {
        if( iequals( dimm_name, ds.get_name() ) ) {
            while( running && !reply_ok ) {
                PointF v = ds.get_r0( bpx::microsec_clock::universal_time(), dimm_duration );
                float r0 = (v.x+v.y)/2;
                if( !r0_reached && (r0 >= dimm_limit) ) {
                    r0_reached = true;
                    previous_r0 = r0;
                }
                if( r0_reached ) {
                    if( trailing && (r0 >= previous_r0) ) {
                        previous_r0 = r0;
                    } else {
                        reply_ok = true;
                        break;
                    }
                }
                this_thread::sleep_for( seconds(poll_interval) );
            }
        }
    }
    
    return reply_ok;
    
}


float AutoSave::get_r0( const vector<DimmSet>& dimm_sets ) {

    float r0(0.0);
    for( auto& ds: dimm_sets ) {
        if( iequals( dimm_name, ds.get_name() ) ) {
            PointF v = ds.get_r0( bpx::microsec_clock::universal_time(), dimm_duration );
            r0 = (v.x+v.y)/2;
        }
    }
    return r0;
    
}

