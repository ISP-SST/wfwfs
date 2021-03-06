#ifndef WFWFS_FILE_FILEMETA_HPP
#define WFWFS_FILE_FILEMETA_HPP
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

#include <exception>
#include <functional>
#include <iostream>
#include <map>

#include <boost/date_time/posix_time/posix_time.hpp>

namespace bpx = boost::posix_time;

namespace wfwfs {



    /*! @ingroup file FileIO
        *  @{
        */

    struct FileMeta {

        FileMeta ( void ) {}
        virtual ~FileMeta ( void ) = default;
        
        virtual std::vector<std::string> getText( bool raw=false ) { return std::vector<std::string>(1,""); }
        virtual size_t getNumberOfFrames(void) { return 1; }
        virtual bpx::ptime getStartTime(void) { return bpx::ptime(); };
        virtual bpx::ptime getEndTime(void) { return bpx::ptime(); };
        virtual bpx::ptime getAverageTime(void) { return bpx::ptime(); };
        virtual bpx::time_duration getExposureTime(void) { return bpx::time_duration(); };
        
        virtual std::vector<bpx::ptime> getStartTimes(void) { return std::vector<bpx::ptime>(); };
        virtual std::vector<size_t> getFrameNumbers(void) { return std::vector<size_t>(1,0); };
        
        virtual size_t dataSize(void) { return 0; };
        virtual size_t dimSize(size_t) { return 0; };
        virtual uint8_t elementSize(void) { return 0; };
        virtual uint8_t nDims(void) { return 0; }
        virtual size_t nElements(void) { return 0; };
        virtual int getIDLType(void) { return 0; };
        virtual int getFormat(void) { return 0; };
        
    };



    template <typename T, typename U>
    inline size_t readOrThrow ( T& strm, U* out, size_t nElements, const std::string& msg = "" ) {
        size_t nBytes = nElements*sizeof ( U );
        strm.read ( reinterpret_cast<char*> ( out ), nBytes );
        if ( !strm.good() ) {
            throw std::ios_base::failure ( "Read failed: " + msg );
        }

        return nBytes;
    }



    template <typename T, typename U>
    inline size_t writeOrThrow ( T& strm, const U* data, size_t nElements, const std::string& msg = "" ) {
        size_t nBytes = nElements*sizeof ( U );
        strm.write ( reinterpret_cast<const char*> ( data ), nBytes );
        if ( !strm.good() ) {
            throw std::ios_base::failure ( "Write failed: " + msg );
        }

        return nBytes;
    }




    /*! @} */

}

#endif // WFWFS_FILE_FILEMETA_HPP
