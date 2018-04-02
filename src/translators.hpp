#ifndef WFWFS_TRANSLATORS_HPP
#define WFWFS_TRANSLATORS_HPP
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

#include "stringutil.hpp"

#include <algorithm>
#include <limits>
#include <iostream>

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>

namespace wfwfs {

    struct BoolTranslator {

        typedef std::string internal_type;

        // Converts a string to bool (i.e. when calling get<bool>(key) on a ptree)
        boost::optional<bool> get_value(const internal_type& str) const {
            if(!str.empty()) {
                using boost::algorithm::iequals;

                if(iequals(str, "true") || iequals(str, "yes") || str == "1") {
                    return boost::optional<bool>(true);
                }
                else {
                    return boost::optional<bool>(false);
                }
            }
            else {
                return boost::optional<bool>(true);      // empty value => flag, so should default to true.
            }
        }

        // Converts a bool to string (i.e. when calling put<bool>(key) on a ptree)
        internal_type put_value(const bool& b) {
            if(b) return "1";
            return "0";
        }
    };


    // property_tree custom translator for vectors
    template <typename T>
    class VectorTranslator {

        typedef std::string internal_type;
        typedef typename std::vector<T> external_type;

        void parseSegment(external_type& out, const internal_type& elem) {
            size_t n1 = std::count(elem.begin(), elem.end(), '-');
            size_t n2 = std::count(elem.begin(), elem.end(), ':');
            if( (n1+n2) == 0 ) {
                out.push_back(boost::lexical_cast<T>(elem));
                return;
            } else if( (n1+n2) == 1 ) {
                size_t n = elem.find_first_of(":-");
                T first = boost::lexical_cast<T>(elem.substr(0, n));
                T last = boost::lexical_cast<T>(elem.substr(n + 1));
                bool rev(false);
                if( first > last ) {
                    std::swap(first,last);
                    rev = true;
                }
                while(first <= last) out.push_back(first++);
                if( rev ) std::reverse( out.begin(), out.end() );
            } else if( n2 == 2 ) {
                size_t n = elem.find_first_of(":");
                size_t nn = elem.find_first_of(":", n + 1);
                T first = boost::lexical_cast<T>(elem.substr(0, n));
                int step = boost::lexical_cast<T>(elem.substr(n + 1, nn - n - 1));
                T last = boost::lexical_cast<T>(elem.substr(nn + 1));
                bool rev(false);
                if( first > last ) {
                    std::swap(first,last);
                    if( step < 0 ) step = abs(step);
                    rev = true;
                }
                while(first <= last) {
                    out.push_back(first);
                    first += step;
                }
                if( rev ) std::reverse( out.begin(), out.end() );
            }
        }


    public:

        // Converts a string to a vector
        external_type get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                std::vector<std::string> tokens;
                boost::split(tokens, str, boost::is_any_of(", "));
                if(std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_signed) {    // only expand '-' & ':' for unsigned integers
                    for(auto & token : tokens) {
                        external_type res;
                        parseSegment( res, token );
                        result.insert( result.end(), res.begin(), res.end() );
                    }
                }
                else {
                    std::transform(tokens.begin(), tokens.end(), back_inserter(result), boost::lexical_cast<T, internal_type>);
                }
            }
            return result;
        }

        // Converts a vector to a comma-separated string
        internal_type put_value(const external_type& vec) {
            std::ostringstream oss;
            if(!vec.empty()) {
                if( std::numeric_limits<T>::is_integer && !std::numeric_limits<T>::is_signed ) {    // only expand '-' & ':' for unsigned integers
                    auto first = vec.begin();
                    while(first != vec.end()) {
                        
                        auto last = first + 1;
                        if(last == vec.end()) { // last element in vector
                            oss << *first;
                            break;
                        }

                        int step = *last - *first;
                        T next = *last + step;
                        size_t count = 2;
                        while(++last != vec.end() && *last == next) {
                            next += step;
                            count++;
                        }
                        if(count > 2) {
                            oss << *first;
                            if( step == 1 ) oss << "-";
                            else  oss << ":";
                            if( step != 1 ) oss << step << ":";
                            oss << *(last-1);
                        } else {
                            oss << *first;
                            last--;
                        }
                        if(last != vec.end()) oss << ",";
                        first = last;

                    }
                }
                else {
                    std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<T>(oss, ","));
                    oss << vec.back();
                }
            }
            return oss.str();
        }
    };


    // property_tree custom translator for vectors
    template <>
    class VectorTranslator<std::string> {

        typedef std::string internal_type;
        typedef typename std::vector<std::string> external_type;

    public:

        // Converts a string to a vector
        external_type get_value(const internal_type& str) {
            external_type result;
            if(!str.empty()) {
                boost::split(result, str, boost::is_any_of(","));
            }
            return result;
        }

        // Converts a vector to a comma-separated string
        internal_type put_value(const external_type& vec) {
            std::ostringstream oss;
            if(!vec.empty()) {
                std::copy(vec.begin(), vec.end() - 1, std::ostream_iterator<std::string>(oss, ","));
                oss << vec.back();
            }
            return oss.str();
        }
    };


}


namespace boost {

    namespace property_tree {

        template <> struct translator_between<std::string, bool> { typedef wfwfs::BoolTranslator type; };

        template<typename T>
        struct translator_between<std::string, std::vector<T>> {
            typedef wfwfs::VectorTranslator<T> type;
        };


    } // namespace property_tree
} // namespace boost

#endif  // WFWFS_TRANSLATORS_HPP
