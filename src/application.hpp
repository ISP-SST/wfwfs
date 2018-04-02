#ifndef WFWFS_APPLICATION_HPP
#define WFWFS_APPLICATION_HPP
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

#include <boost/noncopyable.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>


namespace wfwfs {

    class Application : private boost::noncopyable {
        
    public:
        enum RunMode { LOOP=0, EXIT, RESET, RESTART };

        struct KillException {
            KillException( void ) {}
            virtual const char *what( void ) const throw() {
                return "application killed";
            }
        };

        struct ResetException {
            ResetException( void ) {}
            virtual const char *what( void ) const throw() {
                return "application reset";
            }
        };

        struct RestartException {
            RestartException( void ) {}
            virtual const char *what( void ) const throw() {
                return "application restart";
            }
        };

        struct ThreadExit {
            ThreadExit( void ) {}
            virtual const char *what( void ) const throw() {
                return "exiting thread";
            }
        };


        Application( boost::program_options::variables_map& vm, RunMode=EXIT );
        virtual ~Application( void );

        static void getOptions( boost::program_options::options_description& options, const std::string& );
        static boost::program_options::options_description& getOptionsDescription( const std::string& name = "application" );

        typedef void ( parserFunction )( boost::program_options::options_description&, boost::program_options::variables_map& );

        static std::pair<std::string, std::string> customParser( const std::string& s );
        static boost::program_options::options_description& parseCmdLine( int argc, char* argv[], boost::program_options::variables_map& vm,
                                                      boost::program_options::options_description* programOptions = nullptr,
                                                      boost::program_options::positional_options_description * positionalOptions = nullptr,
                                                      parserFunction customParser = nullptr );

        static void checkGeneralOptions( boost::program_options::options_description& desc, boost::program_options::variables_map& vm );

        int run( void );                    //!< application entry-point, basically just a loop that calls \c doWork()
        
        virtual void reset( void ) { runMode = RESET; };
        virtual void stop( void ) { runMode = EXIT; };

        std::string getName( void ) const;
        static std::string executableName;

    protected:
        /*! Application main method, default application returns false immediately (i.e. exits).
         *   Overload doWork with something interesting. 
         */
        virtual bool doWork( void ) { return false; };
        
        volatile RunMode runMode;
        
        int returnValue;

        static uint8_t logLevel;
        std::string applicationName;
        std::string settingsFile;

        boost::property_tree::ptree propTree;
        
        
    };


}

#endif // WFWFS_APPLICATION_HPP
