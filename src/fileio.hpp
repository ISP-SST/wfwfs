#ifndef WFWFS_FILEIO_HPP
#define WFWFS_FILEIO_HPP
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

#include "filemeta.hpp"
#include "array.hpp"

#include <string>
#include <memory>
#include <thread>
#include <vector>

namespace wfwfs {

    /*! @defgroup file FileIO
        *  @{
        */

    enum Format : uint8_t { FMT_NONE = 0,
                            FMT_ANA,
                            FMT_FITS,
                            FMT_NCDF,
                            FMT_MOMFBD,
                            FMT_PLAIN
                            };


    Format readFmt(const std::string& );
    Format guessFmt(const std::string& );

    template <typename T>
    void getOrRead( const std::string& fn, std::shared_ptr<T>& data );

    std::shared_ptr<wfwfs::FileMeta> getMeta(const std::string& fn, bool size_only=false);

    void readFile( const std::string& fn, char* data, std::shared_ptr<wfwfs::FileMeta>& meta );
    template <typename T>
    void readFile( const std::string& fn, wfwfs::Array<T>& data );
//         template <typename T>
//         void readFile( const std::string& fn, wfwfs::image::Image<T>& data, bool metaOnly=false );

    template <typename T>
    void writeFile( const std::string& fn, wfwfs::Array<T>& data );
//         template <typename T>
//         void writeFile( const std::string& fn, wfwfs::image::Image<T>& data );

    typedef std::function<void(char*,size_t,std::shared_ptr<wfwfs::FileMeta>&)> postLoadCallback;
    typedef std::function<void(char*,size_t,std::shared_ptr<wfwfs::FileMeta>&)> preSumCallback;
    
//         void loadFiles( const std::vector<std::string>& fn, char* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
//                        double* averages=nullptr, double* times=nullptr, std::string progressMsg="" );
    void loadFiles( const std::vector<std::string>& fn, char* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
                    postLoadCallback postLoad = postLoadCallback() );
    void sumFiles( const std::vector<std::string>& fn, double* out, size_t frameSize, uint8_t nThreads=std::thread::hardware_concurrency(),
                    preSumCallback preSum = preSumCallback() );
    
    enum ErrorHandling { EH_PRINT=1, EH_THROW };
    extern ErrorHandling errorHandling;         //<! Specify if routines should throw or print errors. Default is EH_PRINT;
    inline void setErrorHandling( ErrorHandling eh ) { errorHandling = eh; }


//std::shared_ptr<Image>
    /*! @} */


}

#endif // WFWFS_FILEIO_HPP
