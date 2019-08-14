#ifndef WFWFS_ARRAYSTATS_HPP
#define WFWFS_ARRAYSTATS_HPP
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

#include "array.hpp"


namespace wfwfs {

    enum StatType { ST_VALUES=1, ST_RMS, ST_NOISE=4, ST_ALL=7 };

    struct ArrayStats {
        typedef std::shared_ptr<ArrayStats> Ptr;

        ArrayStats() : clip(-1), cutoff(-1), min(0), max(0), median(0), med25(0), med75(0), sum(0), sqr_sum(0), norm(0),
                        mean(0), rms(0), stddev(0), noise(0), noiseRMS(0), nElements(0), hasInfinity(false) {}
        
        template <typename T>
        void getMinMaxMean( const T* data, size_t sizeY, size_t sizeX, size_t stride, const uint8_t* mask=nullptr );
        template <typename T>
        void getMinMaxMean( const T* data, size_t sizeY, size_t sizeX, const uint8_t* mask=nullptr ) {
            getMinMaxMean( data, sizeY, sizeX, sizeX, mask );
        }
        template <typename T> void getMinMaxMean(const T* data, size_t count, const uint8_t* mask=nullptr) {
            getMinMaxMean( data, 1, count, count, mask );
        }
        template <typename T> void getMinMaxMean(const wfwfs::Array<T>& data) {
            if( data.dense() ) getMinMaxMean( data.ptr(), data.nElements() );
            else {
                wfwfs::Array<T> tmp;
                data.copy(tmp);
                getMinMaxMean(tmp.get(),data.nElements());
            }
        }
        
        template <typename T>
        void getMedian( const T* data, size_t sizeY, size_t sizeX, size_t stride, const uint8_t* mask=nullptr );
        template <typename T>
        void getMedian( const T* data, size_t sizeY, size_t sizeX, const uint8_t* mask=nullptr ) {
            getMedian( data, sizeY, sizeX, sizeX, mask );
        }
        template <typename T> void getMedian( const T* data, size_t count, const uint8_t* mask=nullptr ) {
            getMedian( data, 1, count, count, mask );
        }
        template <typename T> void getMedian(const wfwfs::Array<T>& data) {
            if( data.dense() ) getMedian( data.ptr(), data.nElements() );
            else {
                wfwfs::Array<T> tmp;
                data.copy(tmp);
                getMedian( tmp.get(), data.nElements() );
            }
        }

        template <typename T> void getRmsStddev(const T* data, size_t count, const uint8_t* mask=nullptr);
        template <typename T> void getRmsStddev(const wfwfs::Array<T>& data) {
            if (data.dense()) getRmsStddev(data.ptr(),data.nElements());
            else {
                wfwfs::Array<T> tmp;
                data.copy(tmp);
                getRmsStddev(tmp.get(),data.nElements());
            }
        }
        
        template <typename T> void getStats(const T* data, size_t count, int flags=ST_ALL);
        template <typename T> void getStats(const wfwfs::Array<T>& data, int flags=ST_ALL);
        template <typename T> void getStats(uint32_t borderClip, const wfwfs::Array<T>& data, int flags=ST_ALL);
        
        size_t size( void ) const;
        uint64_t pack( char* ) const;
        uint64_t unpack( const char*, bool );

        int clip;
        double cutoff;
        double min, max, median, med25, med75;
        double sum, sqr_sum, norm, mean, rms, stddev;
        double noise, noiseRMS;
        size_t nElements;
        bool hasInfinity;

    };

}   // wfwfs


#endif  // WFWFS_ARRAYSTATS_HPP
