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
        const T* row1 = data + y*stride + x - 1;
        const T* row2 = row1 + stride;
        const T* row0 = row1 - stride;
        
        const double d00 = static_cast<double>( row0[0] );
        const double d01 = static_cast<double>( row0[1] );
        const double d02 = static_cast<double>( row0[2] );
        const double d10 = static_cast<double>( row1[0] );
        const double d11 = static_cast<double>( row1[1] );
        const double d12 = static_cast<double>( row1[2] );
        const double d20 = static_cast<double>( row2[0] );
        const double d21 = static_cast<double>( row2[1] );
        const double d22 = static_cast<double>( row2[2] );
        
        const double a0 = (d21 - d01) * 0.5;
        const double a1 = (d21 - 2*d11 + d01);
        const double a2 = (d12 - d10) * 0.5;
        const double a3 = (d12 - 2*d11 + d10);
        const double a4 = (d22 - d20 - d02 + d00) * 0.25;
        const double den = 1.0 / (a4*a4 - a3*a1);

        x_delta = (a2*a1 - a0*a4)*den;
        y_delta = (a3*a0 - a2*a4)*den;
        
        if( fabs(x_delta) > 0.5 || fabs(y_delta) > 0.5 ) {
            x_delta = -a2/a3;
            y_delta = -a0/a1;
        }

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
