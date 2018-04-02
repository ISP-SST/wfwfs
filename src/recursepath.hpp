#ifndef WFWFS_RECURSEPATH_HPP
#define WFWFS_RECURSEPATH_HPP
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

#include <limits>

#include <boost/filesystem.hpp>


namespace wfwfs {

    /*! @enum EntryType
    *  @brief Identifiers for filesystem entries.s
    */
    enum EntryType { ET_UNDEF=0, ET_DIRECTORY=1, ET_FILE=2, ET_SYMLINK=4 };

    /*!  @class     RecursePath
     *   @brief     Wrapper class that calls a specified function for all files/folders in a tree.
     *   @author    Tomas Hillberg (hillberg@astro.su.se)
     *   @date      2018
     */
    class RecursePath  {
        
        typedef bool (*functype)( const boost::filesystem::path& );

    public:

        RecursePath( boost::filesystem::path& p, functype f, int nSubs=std::numeric_limits<int>::max() );
        RecursePath( const RecursePath& rhs, int sublevels );
        RecursePath( const RecursePath& rhs );

        void operator()( boost::filesystem::directory_entry& p ) const;
        void operator()( const boost::filesystem::path& p ) const;

    private:

        int nSubLevels;
        functype callBack;

    }; // end RecursePath

} // end namespace wfwfs



#endif // WFWFS_RECURSEPATH_HPP
