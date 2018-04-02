#ifndef WFWFS_POINT_HPP
#define WFWFS_POINT_HPP
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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <sstream>
#include <iostream>

#include "datautil.hpp"

namespace wfwfs {


    template<typename T> struct PointType {
        PointType ( T yy=0, T xx=0 ) : x(xx), y(yy) {}
        PointType ( PointType<T>&& rhs ) : x(std::move(rhs.x)), y(std::move(rhs.y)) {}
        PointType ( const PointType<T>& rhs ) : x(rhs.x), y(rhs.y) {}
        template <typename U> PointType ( const PointType<U>& rhs ) : x(rhs.x), y(rhs.y) {}
        static inline uint64_t size(void) { return 2*sizeof(T); };
        uint64_t pack(char* ptr) const {
            uint64_t count = wfwfs::pack(ptr,x);
            count += wfwfs::pack(ptr+count,y);
            return count;
        }
        uint64_t unpack(const char* ptr, bool swap_endian=false) {
            uint64_t count = wfwfs::unpack(ptr,x,swap_endian);
            count += wfwfs::unpack(ptr+count,y,swap_endian);
            return count;
        }
        template <typename U> PointType<T> max(const PointType<U>& rhs) {
            PointType<T> tmp(*this);
            tmp.x=std::max<T>(tmp.x,rhs.x);
            tmp.y=std::max<T>(tmp.y,rhs.y);
            return std::move(tmp); }
        template <typename U> PointType<T> min(const PointType<U>& rhs) {
            PointType<T> tmp(*this);
            tmp.x=std::min<T>(tmp.x,rhs.x);
            tmp.y=std::min<T>(tmp.y,rhs.y);
            return std::move(tmp); }
        template <typename U> PointType<float> projectOnto( const PointType<U>& base, bool preserve_norm=true ) {
            PointType<float> tmp(base), result;
            double d = sqrt( tmp.norm() );
            tmp /= d;                                   // normalize base to simplify projection
            result.x = (tmp.x*x + tmp.y*y);             // longitudinal part (along "base", stored as "x" in result
            result.y = (-tmp.y*x + tmp.x*y);            // transverse part (defined as base + 90 degrees) , stored as "y" in result
            if( !preserve_norm ) {                      // preserve_norm =>  ||result|| == ||this||, else express result in units of ||base||
                result /= d;
            }
            return result;
        }
        inline bool isfinite( void ) const { return std::isfinite(x) && std::isfinite(y); }
        inline T max( void ) const { return std::max( x, y ); }
        inline T min( void ) const { return std::min( x, y ); }
        inline T max_abs( void ) const { return std::max( abs(x), abs(y) ); }
        inline T min_abs( void ) const { return std::min( abs(x), abs(y) ); }
        inline T norm( void ) const { return (x*x+y*y); }
        inline void round( void ) { x=std::round(x); y=std::round(y); }
        operator std::string() const { std::ostringstream out; out << "[" << y << "," << x << "]"; return out.str(); }
        //operator std::string() const { return "(" + std::to_string(y) + "," + std::to_string(x) + ")"; }
        template <typename U> PointType<T>& operator+=(const PointType<U>& rhs) { x += rhs.x; y += rhs.y; return *this; }
        template <typename U> PointType<T>& operator-=(const PointType<U>& rhs) { x -= rhs.x; y -= rhs.y; return *this; }
        template <typename U> PointType<T>& operator*=(const PointType<U>& rhs) { x *= rhs.x; y *= rhs.y; return *this; }
        template <typename U> PointType<T>& operator/=(const PointType<U>& rhs) { x /= rhs.x; y /= rhs.y; return *this; }
        template <typename U> PointType<T> operator+(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res+=rhs); }
        template <typename U> PointType<T> operator-(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res-=rhs); }
        template <typename U> PointType<T> operator*(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res*=rhs); }
        template <typename U> PointType<T> operator/(const PointType<U>& rhs) const { PointType<T> res(*this); return std::move(res/=rhs); }
        PointType<T>& operator=(const PointType<T>& rhs) { x = rhs.x; y = rhs.y; return *this; }
        template <typename U> PointType<T>& operator=(const PointType<U>& rhs) { x = rhs.x; y = rhs.y; return *this; }
        PointType<T>& operator=(const T& rhs) { x = rhs; y = rhs; return *this; }
        PointType<T>& operator+=(const T& rhs) { x += rhs; y += rhs; return *this; }
        PointType<T>& operator-=(const T& rhs) { x -= rhs; y -= rhs; return *this; }
        PointType<T>& operator*=(const T& rhs) { x *= rhs; y *= rhs; return *this; }
        PointType<T>& operator/=(const T& rhs) { x /= rhs; y /= rhs; return *this; }
        PointType<T> operator+(const T& rhs) const { PointType<T> tmp(*this); tmp += rhs; return std::move(tmp); }
        PointType<T> operator-(const T& rhs) const { PointType<T> tmp(*this); tmp -= rhs; return std::move(tmp); }
        PointType<T> operator-(void) const { PointType<T> tmp(-y,-x); return std::move(tmp); }
        PointType<T> operator*(const T& rhs) const { PointType<T> tmp(*this); tmp *= rhs; return std::move(tmp); }
        PointType<T> operator/(const T& rhs) const { PointType<T> tmp(*this); tmp /= rhs; return std::move(tmp); }
        template <typename U> bool operator==(const PointType<U>& rhs) const { return (x == rhs.x && y == rhs.y); }
        bool operator==(T rhs) const { return (x == rhs && y == rhs); }
        template <typename U> bool operator!=(const PointType<U>& rhs) const { return !(*this == rhs); }
        bool operator!=(T rhs) const { return !(*this == rhs); }
        template <typename U> bool operator<(const PointType<U>& rhs) const { if (y==rhs.y) return (x < rhs.x); return (y < rhs.y); }
        T x,y;
    };
    typedef PointType<double> PointD;
    typedef PointType<float> PointF;
    typedef PointType<int> PointI;
    typedef PointType<uint16_t> Point16;
    typedef PointType<uint32_t> Point;
    template<typename T>
    std::ostream& operator<<(std::ostream& os, const PointType<T>& pt) {
        os << (std::string)pt;
        return os;
    }


}

#endif // WFWFS_POINT_HPP
