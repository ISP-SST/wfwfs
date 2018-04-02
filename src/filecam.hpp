#ifndef WFWFS_FILECAM_HPP
#define WFWFS_FILECAM_HPP
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
#include "fileio.hpp"

#include <map>

#include <boost/filesystem.hpp>
namespace bfs = boost::filesystem;


namespace wfwfs {

    class FileCam : public Camera {
        
    public:

        FileCam( boost::property_tree::ptree& );
        ~FileCam( void );

        void init( void );
        void cleanup( void );

        void run( frame_cb_t cb );
        void stop( void );
        
        void copyRawData( const uint8_t* in, uint8_t* out, uint8_t outDepth ) const;
        
        void set_exposure( double e ) {  };     // no changes to exposure/gain allowed, we keep the settings from the file
        void set_gain( double g ) {  };
        
    private:
        void add_file( const bfs::path& );
        void set_path( const bfs::path& );
        std::map< std::string, std::shared_ptr<FileMeta> > files;
        
    };

}   // wfwfs

#endif // WFWFS_FILECAM_HPP
