#ifndef WFWFS_CELL_HPP
#define WFWFS_CELL_HPP
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

#include <memory>
#include <string>

namespace wfwfs {

    struct Cell {

        Cell( PointI p, int id=-1 ) : ID(id), pos(p) {};
        
        int ID;
        PointI pos;                             // pixel-coordinates for the center of the cell
        PointI grid_pos;                        // Nominal grid coordinates, the central subfield is defined to have coordinates (0,0)
        
        inline operator std::string() const { return (std::string)pos; }
        
    };
    inline std::ostream& operator<<( std::ostream& os, const Cell& c ) {
        os << (std::string)c;
        return os;
    }
    

}   // wfwfs


#endif // WFWFS_CELL_HPP
