#ifndef WFWFS_PLEORACAM_HPP
#define WFWFS_PLEORACAM_HPP
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

#ifdef WFWFS_WITH_PLEORA

#include "camera.hpp"

// #include <PvSystem.h>
#include <PvDevice.h>
// #include <PvDeviceGEV.h>
// #include <PvDeviceU3V.h>
#include <PvStream.h>
// #include <PvStreamGEV.h>
// #include <PvStreamU3V.h>
#include <PvBuffer.h>


namespace wfwfs {

    class PleoraCam : public Camera {
        
    public:

        PleoraCam( boost::property_tree::ptree& );
        ~PleoraCam( void );

        void init( void );
        void cleanup( void );

        void run( frame_cb_t cb );
        void stop( void );
        
        void copyRawData( const uint8_t* in, uint8_t* out, uint8_t outDepth ) const;
        
        double get_exposure( void ) const;
        //double get_gain( void ) const;
        //double get_interval( void ) const;
        //double get_temp( void ) const;
        
        void set_exposure( double );
        void set_gain( double );
        void set_interval( double );

        // allow camera-specific parsing of incoming telnet commands 
        //bool parse_cmd( const std::string& cmd, const std::string& what, std::string& line, std::string& reply );
        
    private:
        
        void update_temp( void );
        bool maybeLog( const PvResult& res, std::string descr="" ) const;
        void CreateStreamBuffers( void );

        PvDevice* lDevice;
        PvGenParameterArray* lDeviceParams;
        PvStream *lStream;
    
        std::vector< std::shared_ptr<PvBuffer> > lBuffers;

        
    };

}   // wfwfs


#endif // WFWFS_WITH_PLEORA


#endif // WFWFS_PLEORACAM_HPP
