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
#include "sads.hpp"

#include <numeric>

using namespace wfwfs;
using namespace std;


template <typename T>
bool wfwfs::sads( const T* __restrict__ ref_img, size_t ref_size, size_t ref_stride,
                  const T* __restrict__ img, size_t img_size, size_t img_stride,
                  PointF& shift, uint64_t* __restrict__ tmp ) {

    
    bool ret(false);
    bool swapped(false);
    if( ref_size > img_size ) {     // make sure the reference is the smallest image.
        swapped = true;
        std::swap( ref_img, img );
        std::swap( ref_size, img_size );
        std::swap( ref_stride, img_stride );
    }
    
    int n_positions = (img_size-ref_size)+1;
    float midpoint = (img_size-ref_size)/2;
    memset( tmp, 0, n_positions*n_positions*sizeof(uint64_t) );
    
    shift = -1.0;
    uint64_t min = numeric_limits<uint64_t>::max();

    for( int ys(0); ys<n_positions; ++ys ) {
        for( int xs(0); xs<n_positions; ++xs ) {
            int so = ys*n_positions+xs;
            uint64_t res = 0;
            const T* __restrict__ dPtr = img+ys*img_stride+xs;
            const T* __restrict__ rPtr = ref_img;
            for( size_t y(0); y<ref_size; ++y ) {
                for( size_t x(0); x<ref_size; ++x ) {           // Sum absolute differences for current row
                    res += std::abs( (int64_t)dPtr[x] - (int64_t)rPtr[x] );
                }
                rPtr += ref_stride;                             // Move pointers to next row.
                dPtr += img_stride;
            }
            if( res < min ) {                                   // save best position
                shift.x = xs;
                shift.y = ys;
                min = res;
            }
            tmp[so] = res*res;
        }
    }


    if( (shift.min() > 0) && (shift.max()+1 < n_positions) ) {
        PointD delta;
        subpixel( tmp, n_positions, shift.x, shift.y, delta.x, delta.y );
        if( delta.max_abs() < 1.5 ) {   // only accept corrections < 1.5
            shift -= delta;
        }
        ret = true;
    }
    
    shift -= midpoint;

    if( swapped ) {
        std::swap( shift.x, shift.y );
    }
    return ret;


}
template bool wfwfs::sads( const uint16_t* __restrict__, size_t, size_t,
                           const uint16_t* __restrict__, size_t, size_t,
                           PointF&, uint64_t* __restrict__ );




