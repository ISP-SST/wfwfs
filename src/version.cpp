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
#include "version.hpp"

#include "revision.hpp"

using namespace wfwfs;
using namespace std;

uint64_t wfwfs::getVersionNumber ( void ) {
    return (uint64_t(versionMajor)<<32) + (versionMinor<<24) + (versionPatch<<16) + versionCommit;
}

string wfwfs::getVersionString ( void ) {
    return to_string(versionMajor) + "." + to_string(versionMinor) + "." + to_string(versionPatch) + "-" + to_string(versionCommit);
}

string wfwfs::getVersionString ( uint64_t v ) {
    return to_string((v>>32)&0xFF) + "." + to_string((v>>24)&0xFF) + "." + to_string((v>>16)&0xFF) + "-" + to_string(v&0xFFFF);
}

string wfwfs::getLongVersionString ( bool includeMessage ) {
    
    string ret = to_string(versionMajor) + "." + to_string(versionMinor)
         + "." + to_string(versionPatch)+ "-" + to_string(versionCommit);
#ifdef DEBUG_
    ret += " dbg";
#endif
    if( includeMessage ) ret += " (" + string(commitMessage) + ")";
    return ret;
}
