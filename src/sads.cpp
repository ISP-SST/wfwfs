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
bool wfwfs::sads( const T* __restrict__ c1, size_t c1_size,
                  const T* __restrict__ c2, size_t c2_size,
                  PointF& shift, uint64_t* __restrict__ result ) {
    
    bool ret(false);
    int ss = (c2_size-c1_size)+1;
    int max_shift = (c2_size-c1_size)/2;
    memset( result, 0, ss*ss*sizeof(uint64_t) );
    
    shift = -1;
    uint64_t min = numeric_limits<uint64_t>::max();

    for( int ys(0); ys<ss; ++ys ) {
        for( int xs(0); xs<ss; ++xs ) {
            int so = ys*ss+xs;
            uint64_t res = 0;
            const T* __restrict__ dPtr = c2+ys*c2_size+xs;
            const T* __restrict__ rPtr = c1+max_shift*(c2_size+1);
            for( size_t y(0); y<c1_size; ++y ) {
                for( size_t x(0); x<c1_size; ++x ) {
                    int a = (int)dPtr[x] - (int)rPtr[x];
                    if( a < 0 ) res -= a;
                    else res += a;
                }
                rPtr += c2_size;
                dPtr += c2_size;
            }
            if( res < min ) {
                shift.x = xs;
                shift.y = ys;
                min = res;
            }
            result[so] = res*res;
        }
    }

    if( (shift.min() > 0) && (shift.max()+1 < ss) ) {
        PointD delta;
        subpixel( result, ss, shift.x, shift.y, delta.x, delta.y );
        if( delta.max_abs() < 1 ) {  // only accept corrections < 1 px
            shift += delta;
        }
        ret = true;
    }
    
    shift -= max_shift;
    
    return ret;


}
template bool wfwfs::sads( const uint16_t* __restrict__ c1, size_t c1_size,
                           const uint16_t* __restrict__ c2, size_t c2_size,
                           PointF& shift, uint64_t* __restrict__ result );
