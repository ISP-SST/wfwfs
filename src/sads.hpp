#ifndef WFWFS_SADS_HPP
#define WFWFS_SADS_HPP
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

#include "point.hpp"


namespace wfwfs {


    template <typename T>
    void subpixel( const T* data, size_t stride, int x, int y, double& x_delta, double& y_delta ) {
        
        if( (y-1)*stride+x < 1 ) {          // first element will be at negative offset relative to data
            x_delta = y_delta = NAN;
            return;
        }
        
        // Subpixel interpolation using a quadratic fit
        const T* row2 = data + y*stride + x - 1;
        const T* row3 = row2 + stride;
        const T* row1 = row2 - stride;
        
        const double row10s = static_cast<double>( row1[0] * row1[0] );
        const double row11s = static_cast<double>( row1[1] * row1[1] );
        const double row12s = static_cast<double>( row1[2] * row1[2] );
        const double row20s = static_cast<double>( row2[0] * row2[0] );
        const double row21s = static_cast<double>( row2[1] * row2[1] );
        const double row22s = static_cast<double>( row2[2] * row2[2] );
        const double row30s = static_cast<double>( row3[0] * row3[0] );
        const double row31s = static_cast<double>( row3[1] * row3[1] );
        const double row32s = static_cast<double>( row3[2] * row3[2] );
        
        const double a2 = (row22s - row20s) / 2.0;
        const double a3 = (row22s - 2*row21s + row20s);
        const double a4 = (row31s - row11s) / 2.0;
        const double a5 = (row31s - 2*row21s + row11s);
        const double a6 = (row32s - row30s - row12s + row10s) / 4.0;
        const double a7 = 1.0 / (a6*a6 - a3 * a5);

        x_delta = (a2 * a5 - a4 * a6) * a7;
        y_delta = (a3 * a4 - a2 * a6) * a7;

    }

    template <typename T>
    bool sads( const T* __restrict__ ref_img, size_t ref_size, size_t ref_stride,
               const T* __restrict__ img, size_t img_size, size_t img_stride,
               PointF& shift, uint64_t* __restrict__ tmp );
    
    template <typename T>
    bool sads( const T* __restrict__ ref_img, size_t ref_size,
               const T* __restrict__ img, size_t img_size,
               PointF& shift, uint64_t* __restrict__ tmp ) {
        return sads( ref_img, ref_size, ref_size, img, img_size, img_size, shift, tmp );
    }


}  // namespace wfwfs


#endif // WFWFS_SADS_HPP
