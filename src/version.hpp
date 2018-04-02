#ifndef WFWFS_VERSION_HPP
#define WFWFS_VERSION_HPP
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

#include <string>

namespace wfwfs {
    
    /*! @defgroup WFWFS
     *  @{
     */

    /*!  @file      version.hpp
     *   @brief     Functions for retrieving revision numbers/strings in some organized form.
     *   @name      Version/Revision information
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2018
     */

    /*!
     *   @fn         extern int getVersionNumber (void)
     *   @brief      Returns an integer representing the version
     *   @details    The different version numbers are packed into a 64-bit integer that can be used
     *               to compare versions. E.g. v. 1.2.3-9 would be returned as the integer 4328718345
     *               (= (1<<32)+(2<<24)+(3<<16)+9)
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2018
     */
    extern uint64_t getVersionNumber( void );


    /*!
     *   @fn         extern std::string getVersionString(void)
     *   @brief      Returns a human readable version string, e.g. "1.2.3"
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2018
     */
    extern std::string getVersionString( void );
    extern std::string getVersionString( uint64_t );


    /*!
     *   @fn         extern std::string getLongVersionString ( bool )
     *   @brief      Returns a string with version info and commit message.
     *   @details    E.g. If the build was made 7 commits after the last version-tag, the
     *               result would be something like: "1.2.3 (commit: 7 - Added implementation of getLongVersionString)"
     *
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2018
     */
    extern std::string getLongVersionString( bool includeMessage=true );

    /*! @} */

}

#endif // WFWFS_VERSION_HPP
