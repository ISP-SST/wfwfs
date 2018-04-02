#ifndef WFWFS_FILE_FILEFITS_HPP
#define WFWFS_FILE_FILEFITS_HPP
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

#ifdef WFWFS_WITH_FITS

#include "fileio.hpp"
#include "filemeta.hpp"
// #include "image.hpp"
#include "array.hpp"
// #include "arrayutil.hpp"

#include <fitsio.h>

namespace wfwfs {

    /*! @ingroup file
     *  @{
     */

    /*! Container for reading/writing FITS files.
     *
     */
    struct Fits : public wfwfs::FileMeta {

        enum Magic { MAGIC_FITS = 0x504d4953 }; // = "PMIS"
        enum TypeIndex { FITS_NOTYPE = 0,
                            FITS_BYTE,          // uint8
                            FITS_WORD,          // int16
                            FITS_INT,           // int32
                            FITS_FLOAT,
                            FITS_DOUBLE,
                            FITS_COMPLEX,
                            FITS_STRING,
                            FITS_DCOMPLEX=9,
                            FITS_UWORD=12,
                            FITS_UINT,
                            FITS_LONG,
                            FITS_ULONG };
        static const uint8_t typeSizes[];   // = { 0, 1, 2, 4, 4, 8, 8, 0, 0, 16 };
        
        Fits( void );
        Fits( const std::string& );
        ~Fits();

        void close( void );
        void read( const std::string& );

        void write( std::ofstream& );

        std::vector<std::string> getText( bool );
        
        template <typename T>
        static std::string makeCard( std::string key, T value, std::string comment="" );
        static void insertCard( std::vector<std::string>& hdr, std::string card, size_t location=std::string::npos );
        static void insertCardAfter( std::vector<std::string>& hdr, std::string card, std::string after );
        static void insertCardBefore( std::vector<std::string>& hdr, std::string card, std::string before );
        static bool updateCard( std::vector<std::string>& hdr, size_t location, std::string card );
        static bool updateCard( std::vector<std::string>& hdr, std::string key, std::string card );
        static bool updateCard( std::vector<std::string>& hdr, std::string card );
        template <typename T>
        static T getValue( const std::vector<std::string>& hdr, std::string key);
        template <typename T>
        std::vector<T> getTableArray( std::string key );
        
        size_t getNumberOfFrames(void);
        bpx::ptime getStartTime(void);
        bpx::ptime getEndTime(void);
        bpx::ptime getAverageTime(void);
        bpx::time_duration getExposureTime(void);
        std::vector<bpx::ptime> getStartTimes(void);
        std::vector<size_t> getFrameNumbers(void);
        
        size_t dataSize(void);
        size_t dimSize(size_t);
        uint8_t elementSize(void);
        uint8_t nDims(void) { return primaryHDU.nDims; }
        size_t nElements(void);
        int getIDLType(void);
        
        double getMinMaxMean( const char* data, double* Min=nullptr, double* Max=nullptr );
        int getFormat(void) { return FMT_FITS; };

        struct hdu {
            int bitpix;
            int nDims;
            int dataType;               // data type as defined in cfitsio
            size_t elementSize;         // element size (in bytes) = abs(bitpix/8)
            size_t nElements;
            std::vector<int> dims;
            std::vector<std::string> cards;
            virtual void dummy(void)=0;
        };
        
        struct image_hdu : public hdu {
            image_hdu() : dHDU(0) {}
            void dummy(void){};
            int dHDU;       // index to hdu containing data, e.g. compressed tile image.
        };
        
        struct ascii_hdu : public hdu {
            void dummy(void){};
            struct table_info_t {
                int columnStart;            // = TBCOL, offset where this column starts
                std::string columnName;     // = TTYPEn, name of this data-column
                std::string columnFormat;   // = TFORM, Fortran ISO 2004 format string
                std::string columnUnit;     // = TUNIT, physical unit of the data
            };
            uint16_t nColumns;              // = TFIELDS, number of columns in this table
            std::string name;               // = EXTNAME
            std::vector<table_info_t> table_info;
            wfwfs::Array<char> data;
        };
        
        struct binary_hdu : public hdu {
            void dummy(void){};
            std::vector<std::string> data;
        };
        
        struct image_hdu primaryHDU;
        std::vector<std::shared_ptr<struct hdu>> extHDUs;
        
        fitsfile* fitsPtr_;
        int status_;

        /*! @name Read
            *  @brief Load a FITS file into a data block
            */
        //@{
        static void read( std::shared_ptr<wfwfs::Fits>& hdr, char* data );
        template <typename T>
        static void read( const std::string& filename, wfwfs::Array<T>& data, std::shared_ptr<wfwfs::Fits>& hdr );
//             template <typename T>
//             static void read( const std::string& filename, wfwfs::image::Image<T>& data, bool metaOnly=false );
        //@}
        
        /*! @name Write
            *  @brief Write data block into an FITS file.
            */
        //@{
        static void write( const std::string& filename, const char* data, std::shared_ptr<wfwfs::Fits> hdr, bool compress = false, int slice=5 );
        template <typename T>
        static void write( const std::string& filename, const wfwfs::Array<T>& data, std::shared_ptr<wfwfs::Fits> hdr=0, int sliceSize=0 );
//             template <typename T>
//             static void write( const std::string& filename, const wfwfs::image::Image<T>& image, int sliceSize=0 );
        template <typename T>
        static void write( const std::string& filename, const T* data, size_t n=1 );
        template <typename T>
        static void write( const std::string& filename, const std::vector<T>& v ) { write(filename,v.data(),v.size()); }
        //@}
        

    };

    /*! @} */

} // end namespace wfwfs

#endif  // WFWFS_WITH_FITS

#endif // WFWFS_FILE_FILEFITS_HPP
