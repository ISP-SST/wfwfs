#ifndef WFWFS_TYPES_HPP
#define WFWFS_TYPES_HPP
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

#include <cstdint>
#include <complex>

#ifdef WFWFS_WITH_FFTW3
#   include <fftw3.h>
#else
    typedef std::complex<double> fftw_complex;
#endif

namespace wfwfs {


        /*!  @ingroup util
         *  @{
         */
        
        /*!  @file      types.hpp
         *   @brief     Generic type definitions, commonly used all over the library
         *   @name      Types
         * 
         *   @author    Tomas Hillberg (hillberg@astro.su.se)
         *   @date      2018
         */
        
        
        typedef std::complex<double> complex_t;
        
#ifdef WFWFS_WITH_FFTW3
        inline const complex_t& operator*=( complex_t& lhs, const fftw_complex& rhs ) { return lhs *= reinterpret_cast<const complex_t&>(rhs); }
        inline const fftw_complex& operator*=( fftw_complex& lhs, const complex_t& rhs ) { reinterpret_cast<complex_t&>(lhs) *= rhs; return lhs; }
#endif
        
        /* @} */

}



#endif // WFWFS_TYPES_HPP


