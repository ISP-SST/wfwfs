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
#include "arraystats.hpp"

#include <limits>

using namespace wfwfs;
using namespace std;

template <typename T>
void wfwfs::ArrayStats::getMinMaxMean(const T* data, size_t n) {
    
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    sum = 0;
    sqr_sum = 0;
    norm = 0;
    hasInfinity = false;
    const T* ptr = data + n;
    size_t count(0);
    while( ptr-- > data ) {
        double tmp = static_cast<double>( *ptr );
        if( tmp > max ) {
            max = tmp;
        } else if( tmp < min ) {
            min = tmp;
        }
        if ( isfinite( tmp ) ) { // will exclude NaN
            sum += tmp;
            norm += abs(tmp);
            sqr_sum += tmp*tmp;
            count++;
        } else hasInfinity = true;
    }
    mean = sum;
    if(count) mean /= count;
    
}
template void wfwfs::ArrayStats::getMinMaxMean( const char*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const uint8_t*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const int16_t*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const uint16_t*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const int32_t*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const uint32_t*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const float*, size_t );
template void wfwfs::ArrayStats::getMinMaxMean( const double*, size_t );

           
template <typename T>
void wfwfs::ArrayStats::getRmsStddev(const T* data, size_t n) {
    if(n > 1) {
        size_t count(0);
        rms = stddev = 0;
        const T* ptr = data + n;
        while( ptr-- > data ) {
            double tmp = static_cast<double>( *ptr );
            if( tmp == tmp ) { // will exclude NaN
                rms += tmp*tmp;
                count++;
            }
        }
        if (count) {
            rms /= count;
            stddev = sqrt(rms - mean*mean);
            rms = sqrt(rms);
        }
    } else if (n == 1) {
        double tmp = static_cast<double>( *data );
        rms = isfinite(tmp)?tmp:0;
        stddev = 0;
    }
}
template void wfwfs::ArrayStats::getRmsStddev( const float*, size_t );
template void wfwfs::ArrayStats::getRmsStddev( const double*, size_t );
template void wfwfs::ArrayStats::getRmsStddev( const int16_t*, size_t );
template void wfwfs::ArrayStats::getRmsStddev( const uint16_t*, size_t );
template void wfwfs::ArrayStats::getRmsStddev( const int32_t*, size_t );
template void wfwfs::ArrayStats::getRmsStddev( const uint32_t*, size_t );


template <typename T>
void wfwfs::ArrayStats::getStats( const Array<T>& data, int flags ) {

    getMinMaxMean(data);
    
    // TODO: median
    
    if( flags & ST_RMS ) {
        getRmsStddev(data);
    }
}
template void wfwfs::ArrayStats::getStats( const Array<float>&, int );
template void wfwfs::ArrayStats::getStats( const Array<double>&, int );
template void wfwfs::ArrayStats::getStats( const Array<int16_t>&, int );
template void wfwfs::ArrayStats::getStats( const Array<uint16_t>&, int );
template void wfwfs::ArrayStats::getStats( const Array<int32_t>&, int );
template void wfwfs::ArrayStats::getStats( const Array<uint32_t>&, int );


template <typename T>
void wfwfs::ArrayStats::getStats( const T* data, size_t count, int flags ) {

    getMinMaxMean(data, count);
    
    // TODO: median
    
    if( flags & ST_RMS ) {
        getRmsStddev(data,count);
    }
/*    
    if( flags & ST_NOISE ) {
        getNoise(data,count);
    }*/

}
template void wfwfs::ArrayStats::getStats( const float*, size_t, int );
template void wfwfs::ArrayStats::getStats( const double*, size_t, int );
template void wfwfs::ArrayStats::getStats( const int16_t*, size_t, int );
template void wfwfs::ArrayStats::getStats( const uint16_t*, size_t, int );
template void wfwfs::ArrayStats::getStats( const int32_t*, size_t, int );
template void wfwfs::ArrayStats::getStats( const uint32_t*, size_t, int );



template <typename T>
void wfwfs::ArrayStats::getStats( uint32_t borderClip, const Array<T>& data, int flags ) {

    if( !borderClip ) {
        getStats( data, flags );
        return;
    }

    std::vector<int64_t> first, last;
    const std::vector<size_t>& dims = data.dimensions();
    size_t nDims = dims.size();
    for( size_t i = 0; i < nDims; ++i ) {
        if( dims[i] == 1 ) {            // if dimSize == 1 we don't clip it.
            first.push_back( data.first()[i] );
            last.push_back( data.last()[i] );
        }
        else {
            if( 2 * borderClip > dims[i] ) {    // all data clipped, nothing to do
                return;
            }
            first.push_back( data.first()[i] + borderClip );
            last.push_back( data.last()[i] - borderClip );
        }
    }
    Array<T> clippedImage( data, first, last );
    getStats(clippedImage,flags);
    
}
template void wfwfs::ArrayStats::getStats( uint32_t, const Array<float>&, int );
template void wfwfs::ArrayStats::getStats( uint32_t, const Array<double>&, int );
template void wfwfs::ArrayStats::getStats( uint32_t, const Array<int16_t>&, int );
template void wfwfs::ArrayStats::getStats( uint32_t, const Array<uint16_t>&, int );
template void wfwfs::ArrayStats::getStats( uint32_t, const Array<int32_t>&, int );
template void wfwfs::ArrayStats::getStats( uint32_t, const Array<uint32_t>&, int );

/*
            int clip;
            double cutoff;
            double min, max, median;
            double sum, mean, rms, stddev;
            double noise, noiseRMS;
            uint8_t noiseType;              // flag indicating noise statistics (not used atm.)
*/
size_t wfwfs::ArrayStats::size( void ) const {
    return sizeof(int) + 12*sizeof(double);
}



uint64_t wfwfs::ArrayStats::pack( char* ptr ) const {
    using wfwfs::pack;
    uint64_t count = pack(ptr, clip);
    count += pack(ptr+count, cutoff);
    count += pack(ptr+count, min);
    count += pack(ptr+count, max);
    count += pack(ptr+count, median);
    count += pack(ptr+count, sum);
    count += pack(ptr+count, sqr_sum);
    count += pack(ptr+count, norm);
    count += pack(ptr+count, mean);
    count += pack(ptr+count, rms);
    count += pack(ptr+count, stddev);
    count += pack(ptr+count, noise);
    count += pack(ptr+count, noiseRMS);
    return count;
}


uint64_t wfwfs::ArrayStats::unpack( const char* ptr, bool swap_endian ) {
    using wfwfs::unpack;
    uint64_t count = unpack(ptr, clip, swap_endian);
    count += unpack(ptr+count, cutoff, swap_endian);
    count += unpack(ptr+count, min, swap_endian);
    count += unpack(ptr+count, max, swap_endian);
    count += unpack(ptr+count, median, swap_endian);
    count += unpack(ptr+count, sum, swap_endian);
    count += unpack(ptr+count, sqr_sum, swap_endian);
    count += unpack(ptr+count, norm, swap_endian);
    count += unpack(ptr+count, mean, swap_endian);
    count += unpack(ptr+count, rms, swap_endian);
    count += unpack(ptr+count, stddev, swap_endian);
    count += unpack(ptr+count, noise, swap_endian);
    count += unpack(ptr+count, noiseRMS, swap_endian);
    return count;

}
