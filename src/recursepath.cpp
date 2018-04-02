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
#include "recursepath.hpp"

#include <iostream>

using namespace wfwfs;
using namespace std;

namespace fs = boost::filesystem;


RecursePath::RecursePath( boost::filesystem::path& p, functype f, int sublevels ) : nSubLevels( sublevels ), callBack(f) {
    this->operator()( p );
}


RecursePath::RecursePath( const RecursePath& rhs, int sublevels ) : nSubLevels( sublevels ) {
    
    callBack = rhs.callBack;

}


RecursePath::RecursePath( const RecursePath& rhs ) {
    
    nSubLevels = rhs.nSubLevels-1;
    callBack = rhs.callBack;

}


void RecursePath::operator()( fs::directory_entry& p ) const {
    this->operator()( p.path() );
}


void RecursePath::operator()( const fs::path& p ) const {

    if( !callBack( p ) || ( nSubLevels < 0 ) ) {
        return;
    }

    try {
        fs::directory_iterator it( p );
        fs::directory_iterator end;
        for_each( it, end, RecursePath( *this, nSubLevels - 1 ) );
    }
    catch( fs::filesystem_error& fex ) {
        cerr << "Error in RecursePath at path \"" << p.string() << "\"" << endl;
        cerr << fex.what() << endl;
    }

}
