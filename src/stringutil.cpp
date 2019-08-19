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

#include "translators.hpp"
#include "datautil.hpp"

#include <cstdlib>
#include <mutex>
#include <pwd.h>
#include <unistd.h>

#ifdef __GNUG__
#   include <cxxabi.h>
#endif

#include <boost/algorithm/string.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>


namespace bfs = boost::filesystem;
namespace bpt = boost::property_tree;
namespace bpo = boost::program_options;


using namespace wfwfs;
using namespace std;


bool wfwfs::onlyDigits(const string &s) {
    static const boost::regex e("[0-9]+");
    return boost::regex_match(s, e);
}

bool wfwfs::onlyAlpha(const string &s) {
    static const boost::regex e("^[[:alpha:]]*$");
    return boost::regex_match(s, e);
}

bool wfwfs::onlyAlnum(const string &s) {
    static const boost::regex e("^[[:alnum:]]*$");
    return boost::regex_match(s, e);
}

bool wfwfs::onlyHex(const string &s) {
    static const boost::regex e("^[[:xdigit:]]*$");
    return boost::regex_match(s, e);
}

bool wfwfs::isInteger(const string &s) {
    static const boost::regex e("^\\s*(-|0x|0X)?[[:xdigit:]]+\\s*$");
    return boost::regex_match(s, e);
}


bool wfwfs::isHex(const string &s) {
    static const boost::regex e("^\\s*[[:xdigit:]]+\\s*$");
    return boost::regex_match(s, e);
}


bool wfwfs::contains ( const string & haystack, const string & needle, bool ignoreCase ) {

    auto it = std::search (
                  haystack.begin(), haystack.end(),
                  needle.begin(),   needle.end(),
    [ignoreCase] ( char ch1, char ch2 ) {
        if ( ignoreCase ) return std::toupper ( ch1 ) == std::toupper ( ch2 );
        return ch1 == ch2;
    }
              );
    if ( it != haystack.end() ) return true;
    return false;

}


string wfwfs::replace_n( std::string input, const std::string& loc, const std::string& replace, size_t n) {
    
    if(n == 0) return input;
    
    size_t count(0);
    for( size_t offset = input.find(loc); offset != std::string::npos; ) {
        count++;
        offset = input.find(loc, offset+loc.length());
    }

    count = std::min(n, count);
    
    while( count-- ) {
        boost::replace_first( input, loc, replace );
    }
    
    return input;
    
}


string wfwfs::popword( std::string &line, const char *separator ) {
    
    string result;
    size_t pos = line.find_first_not_of( separator );
    if( pos == string::npos ) line.clear();
    else if( pos ) line.erase(0, pos);
    
    pos = line.find_last_not_of( separator );
    if( pos != string::npos ) line.erase(pos+1);
    else line.clear();
    
    if( !line.empty() && line[0] == ':' ) {     // TBD: is this good? doing it now for complying with Guus' code.
        result = line.substr( 1, string::npos );
        line.clear();
    } else {
        pos = line.find_first_of( separator );
        result = string( line, 0, pos );
        line.erase( 0, pos );
        pos = line.find_first_not_of( separator );
        if( pos == string::npos ) line.clear();
        else if( pos ) line.erase(0, pos);
        
        pos = line.find_last_not_of( separator );
        if( pos != string::npos ) line.erase(pos+1);
        else line.clear();
    }
    
//    boost::trim_left( line );
//    boost::trim_left( result );
    
    return result;
    
}


bool wfwfs::nocaseLess(const string& lhs, const string& rhs) {

          return std::lexicographical_compare( lhs.begin (), lhs.end (), rhs.begin (), rhs.end (),
                                               [] ( char ch1, char ch2 ) { return std::toupper( ch1 ) < std::toupper( ch2 ); }
                                             );

}


bool wfwfs::isRelative( const std::string &s ) {
    return (!s.empty() && s[0] != '/');
}


void wfwfs::maybeCreateDir( const bfs::path& p ) {
    
    if( !p.empty() && !bfs::exists(p) ) {
        if( !bfs::create_directories(p) ) {
            cout << boost::format( "failed to create directory: %s" ) % p << endl;
        }
    }

}


vector<set<string>> wfwfs::make_template( const vector<string>& list, string& out, string split_chars ) {
    
    vector<string> segments;
    vector< set<string> > seg_list;

    if( split_chars.empty() ) return seg_list;
    
    size_t n_segments(0);
    for( auto& str: list ) {
        boost::split(segments, str, boost::is_any_of(split_chars));
        if( !n_segments && (segments.size() >= 1) ) {    // TBD: how to treat size=1 (no split), return tpl="" or tpl="%1"
            n_segments = segments.size();
            seg_list.resize( n_segments );
        }
        if( segments.size() == n_segments ) {
            for( size_t j=0; j<n_segments; ++j ) {
                seg_list[j].insert(segments[j]);
            }
        } else if(n_segments) {
            //cout << "Multiple segment sizes!!  nS=" << n_segments << "  ss=" << segments.size() << endl;
        }
    }
    
    int arg_cnt(0);
    out = "";
    for( size_t j=0; j<n_segments; ++j ) {
        size_t nArgs = seg_list[j].size();
        if( nArgs > 1 ) {
            out += "%"+to_string(++arg_cnt);     // add a placeholder for this segment
        } else if( nArgs == 1 ) {
            out += *(seg_list[j].begin());       // a unique item, add it to template
        }
        if( j < n_segments-1 ) out += split_chars[0];       // separator    TBD: should this be a parameter
    }

    return std::move(seg_list);
    
}


string wfwfs::alignCenter(const string& s, size_t n, unsigned char c) {

    string ret(s);
    if( s.length() > n ) {
        ret.resize( n, c );
    } else {
        size_t nn = (n-s.length())/2;
        // if "n-s" is odd, the extra character goes on the left side
        ret.insert( ret.end(), nn, c );
        ret.insert( ret.begin(), (nn+nn%2), c );
    }
    return ret;

}


string wfwfs::alignLeft(const string& s, size_t n, unsigned char c) {

    string ret(s);
    ret.resize( n, c );
    return ret;

}


string wfwfs::alignRight(const string& s, size_t n, unsigned char c) {
    
    string ret(s);
    if( s.length() > n ) {
        ret.resize( n, c );
    } else {
        size_t nn = n-s.length();
        ret.insert( ret.begin(), nn, c );
    }
    return ret;

}


string wfwfs::getUname(__uid_t id) {
    
    static map<__uid_t,string> users;
    static mutex mtx;
    lock_guard<mutex> lock(mtx);
    
    if(!id) id = geteuid();
    
    auto it = users.find(id);
    if( it != users.end() ) return it->second;

    string tmp;
    struct passwd pwent;
    struct passwd *pwentp;
    char buf[1024];
    if( !getpwuid_r( id, &pwent, buf, 1024, &pwentp ) ) {
        tmp = pwent.pw_name;
    }
    else {
        tmp = std::to_string((int)id);
    }
    
    users.emplace( id, tmp );
    
    return tmp;
}


string expandTilde(string in) {
    if(in.empty() || in[0] != '~') return in;
    string tmp;
    size_t cut;
    struct passwd pwent;
    struct passwd *pwentp;
    char buf[1024];
    if(in.length() == 1 || in[1] == '/') {
        cut = 1;
        tmp = getenv("HOME");
        if(tmp.empty()) {
            if( !getpwuid_r( geteuid(), &pwent, buf, 1024, &pwentp ) ) {
                tmp = pwent.pw_dir;
            }
        }
    }
    else {
        cut = in.find_first_of('/');
        string user = in.substr(1, cut - 1);
        if( !getpwnam_r(user.c_str(), &pwent, buf, 1024, &pwentp) ) {
            tmp = pwent.pw_dir;
        }
    }
    if(tmp.empty()) return in;
    else return tmp + in.substr(cut);
}


string wfwfs::cleanPath(string in, string base) {

    if(in.empty()) return in;
    bfs::path fn, result, ain(in);
    if(!base.empty() && base[0] == '~') base = expandTilde(base);
    if(!base.empty() && base[0] != '/') result = bfs::current_path() / bfs::path(base);
    if(in[0] == '~') in = expandTilde(in);

    if(bfs::is_regular_file(ain)) {
        fn = ain.filename();
        ain = ain.parent_path();
    }
    auto it = ain.begin();
    if(in[0] != '/' && !base.empty()) {
        if(!bfs::is_directory(result)) return in;
    }
    else result = *it++;

    bool docanonical WFWFS_UNUSED = (result.string()[0] == '/');         // don't canonicalize relative paths
    for(; it != ain.end(); ++it) {
        if(*it == "..") result = result.parent_path();
        else if(*it != ".") {
#if BOOST_VERSION > 104800
            if(!exists(result / *it) && docanonical) {      // canonicalize the existing part (boost >= 1.48)
                result = canonical(result);
                docanonical = false;
            }
#else
            docanonical = false;
#endif
            result /= *it;
        }
    }

    return (result / fn).string();

}

void wfwfs::printProgress( const string& text, float progress ) {
    static mutex mtx;
    unique_lock<mutex> lock(mtx);
    if( progress >= 0 ) {
        printf( "\r%s (%.1f%%)", text.c_str(), progress );
    } else printf( "\r%s", text.c_str() );
    fflush(stdout);
}


template <typename T>
vector<T> wfwfs::stringToUInts(const string& str) {
    
    bpt::ptree tmpTree;                         // just to be able to use the VectorTranslator
    tmpTree.put( "tmp", str );
    return tmpTree.get<vector<T>>( "tmp", vector<T>() );

}
template vector<uint8_t> wfwfs::stringToUInts(const string&);
template vector<uint16_t> wfwfs::stringToUInts(const string&);
template vector<uint32_t> wfwfs::stringToUInts(const string&);
template vector<uint64_t> wfwfs::stringToUInts(const string&);


template <typename T>
std::string wfwfs::uIntsToString(const std::vector<T>& ints) {
    
    bpt::ptree tmpTree;
    tmpTree.put("tmp", ints);
    return tmpTree.get<string>( "tmp", "" );
    
}
template string wfwfs::uIntsToString(const vector<uint8_t>& );
template string wfwfs::uIntsToString(const vector<uint16_t>& );
template string wfwfs::uIntsToString(const vector<uint32_t>& );
template string wfwfs::uIntsToString(const vector<uint64_t>& );


std::string wfwfs::colorString( const std::string& in, StringColor col ) {

    static std::string a("\033[");
    static std::string b("m");
    static std::string c("\033[0m");
    return a + to_string(col) + b + in + c;

}


string wfwfs::tvToString( const timeval& a, bool millis ) {


    char tmp[15];
    strftime( tmp, 14, "%H:%M:%S", gmtime( &a.tv_sec ) );
    string ret( tmp );

    if ( millis ) {
        sprintf( tmp, ".%.3u", ( uint )( a.tv_usec / 1000 ) );
        ret += tmp;
    }

    return ret;

}


string wfwfs::tsToString( const timespec& a, bool millis ) {


    char tmp[15];
    strftime( tmp, 14, "%H:%M:%S", gmtime( &a.tv_sec ) );
    string ret( tmp );

    if ( millis ) {
        sprintf( tmp, ".%.6u", ( uint )( a.tv_nsec / 1000 ) );
        ret += tmp;
    }

    return ret;

}


#ifdef __GNUG__
    std::string wfwfs::demangle_name( const string& name ) {
        int status(0);
        std::unique_ptr<char, void(*)(void*)> res {
            abi::__cxa_demangle( name.c_str(), 0, 0, &status ),
            std::free
        };
        return (status==0) ? res.get() : name ;
    }
#else
    std::string wfwfs::demangle_name( const string& name ) { return name; }    // TODO: implement
#endif


string wfwfs::demangle_symbol( const string& sym ) {
    // Example symbol: ./module(func_name+0x15c) [0x8048a6d]
    static const boost::regex sym_re("([[:alnum:]_]+)[+ ]+([0x[:xdigit:]]+)");      // match func_name & offset
    boost::regex re( "(\\d+)+" );
    boost::smatch match;
    string ret = sym;
    if( boost::regex_search( sym, match, sym_re ) ) {
        string sym_name = string(match[1]);
        ret = demangle_name( sym_name );
    }
    return ret;

}


string wfwfs::printVariableMap( const bpo::variables_map vm ) {
    string ret;
    vm.begin();
    for( auto& it: vm ) {
        ret += it.first;
        const boost::any& value = it.second.value();
        if( !value.empty() ) {
            ret += "=";
            if( value.type() == typeid(int) ) ret += to_string( it.second.as<int>() );
            else if( value.type() == typeid(uint16_t) ) ret += to_string( it.second.as<uint16_t>() );
            else if( value.type() == typeid(uint32_t) ) ret += to_string( it.second.as<uint32_t>() );
            else if( value.type() == typeid(string) ) ret += it.second.as<string>();
            else if( value.type() == typeid(vector<string>) ) ret += printArray(it.second.as<vector<string>>(),"");
            else if( value.type() == typeid(vector<uint32_t>) ) ret += printArray(it.second.as<vector<uint32_t>>(),"");
            else {
                ret += " printVariableMap needs printing for type (" + demangle_name( value.type().name() ) + ")";
            }
        }
        if( it.second.defaulted() ) {
            ret += " (default)";
        }
        ret += "\n";
    }
    return ret;
}
