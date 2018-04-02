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
#include "camera.hpp"

#include "stringutil.hpp"

#include <boost/format.hpp>

using namespace wfwfs;
using namespace std;


Camera::Camera( boost::property_tree::ptree& pt ) : settings(pt), running(false) {

    cfg.height = settings.get<int>( "height", 1024 );
    cfg.width = settings.get<int>( "width", 1024 );
    cfg.depth = settings.get<int>( "depth", 8 );
    cfg.n_buf = settings.get<int>( "nframes", 8 );
 
    state.exposure = settings.get<double>( "exposure", 0.01 );
    state.interval = settings.get<double>( "interval", 0.01 );
 
    info.id = settings.get<string>( "camera_id", "" );

}


Camera::~Camera( void ) {

}


bool Camera::parse_cmd( const std::string& cmd, const std::string& what, std::string& extra, std::string& reply ) {

    reply.clear();
    if( cmd == "set" ) {
        double value = pop<double>( extra );
        if( what == "exposure" ) {
            set_exposure( value );
            reply = boost::str( boost::format("OK exposure %lf") % get_exposure() );
        } else if( what == "gain" ) {
            set_gain( value );
            reply = boost::str( boost::format("OK gain %lf") % get_gain() );
        } else if( what == "interval" ) {
            set_interval( value );
            reply = boost::str( boost::format("OK interval %lf") % get_interval() );
        } else if( what == "offset" ) {
            set_offset( value );
            reply = boost::str( boost::format("OK offset %lf") % get_offset() );
        }
    } else if( cmd == "get" ) {
        if( what == "exposure" ) {
            reply = boost::str( boost::format("OK exposure %lf") % get_exposure() );
        } else if( what == "gain" ) {
            reply = boost::str( boost::format("OK gain %lf") % get_gain() );
        } else if( what == "interval" ) {
            reply = boost::str( boost::format("OK interval %lf") % get_interval() );
        } else if( what == "offset" ) {
            reply = boost::str( boost::format("OK offset %lf") % get_offset() );
        } else if( what == "cfg" ) {
            reply = boost::str( boost::format("OK cfg %dx%d, %d-bits, %d buffers") % cfg.width % cfg.height % cfg.depth % cfg.n_buf );
        } else if( what == "info" ) {
            reply = boost::str( boost::format("OK info ID: \"%s\"  MODEL: \"%s\"  firmware: \"%s\"  serial: \"%s\"") % info.id % info.model % info.firmware_version % info.serial );
        } else if( what == "settings" ) {
            reply = boost::str( boost::format("OK settings exposure: %lf  gain: %lf  interval: %lf  offset: %lf")
            % state.exposure % state.gain % state.interval % state.offset );
        }
    }
    return !reply.empty();
}
