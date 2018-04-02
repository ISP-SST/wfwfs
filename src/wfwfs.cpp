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
#include "daemon.hpp"

#include "stringutil.hpp"

#include <syslog.h>

using namespace wfwfs;
using namespace std;

namespace bpo = boost::program_options;


namespace {
    
    // define options specific to this binary
    bpo::options_description getOptions( void ) {

        bpo::options_description options( "WFWFS Options" );
        options.add_options()
            ( "config,c", bpo::value<string>()->default_value( "/etc/ccd/config" ), "Configuration file to use." )
            ( "foreground,F", "Do not detach/background process.")
            ( "port,p", bpo::value<uint16_t>()->default_value( 15000 ), "Port to listen on, or connect to."
            " The environment variable WFWFS_PORT will be used as default if it is defined." )
            ( "threads,t", bpo::value<uint16_t>()->default_value( 0 ), "max number of threads to use.")
        ;

        return options;
    }

    // define environment variables to use as defaults if the corresponding command-line option is not specified
    string environmentMap( const string &envName ) {

        static const map<string, string> vmap = {
            { "WFWFS_CFG", "config" },      // If WFWFS_CFG is defined in the environment, it will override the default value for config
            { "WFWFS_PORT", "port" }        // If WFWFS_PORT is defined in the environment, it will override the default value for port
        };

        auto ci = vmap.find( envName );
        if( ci == vmap.end() ) {
            return "";
        } else {
            return ci->second;
        }
    }


}


int main( int argc, char *argv[] ) {

    bpo::variables_map vm, tmp_vm;
    bpo::options_description programOptions = getOptions();

    try {
        
        bpo::options_description& allOptions = Application::parseCmdLine( argc, argv, vm, &programOptions );

        // load matched environment variables according to the getOptionName() above.
        bpo::store( bpo::parse_environment( allOptions, environmentMap ), vm );
        
#if BOOST_VERSION > 104800  // FIXME verify exactly which version notify appears in
        vm.notify();
#endif
    
        bool detach = (vm.count( "log-stdout" ) == 0) && (vm.count( "foreground" ) == 0);
        int options = detach ? 0 : LOG_CONS | LOG_PERROR;
        options |= LOG_PID | LOG_NDELAY;
        openlog( "WFWFS", options, LOG_DAEMON );
        
        if( detach && daemon( 1, 0 ) ) {
            throw runtime_error( string("Failed to background process: ") + strerror( errno ) );
        }
        tmp_vm = vm;    // make a local copy so it can be restored to cmd-line options on reset
        while( true ) {
            try {
                Daemon daemon( tmp_vm );
                return daemon.run();
            }
            catch( Application::KillException ) {
                break;
            }
            catch( Application::ResetException ) {
                syslog( LOG_NOTICE, "Resetting WFWFS daemon" );
                tmp_vm = vm;
                continue;
            }
            catch( Application::RestartException ) {
                syslog( LOG_NOTICE, "Restarting WFWFS daemon" );
                continue;
            }
        }
    } catch( const exception &e ) {
        //syslog( LOG_EMERG, "Uncaught exception (fatal): %s", e.what() );
        cerr << "Uncaught exception (fatal): " << e.what() << endl;
    } catch( ... ) {
        //syslog( LOG_EMERG, "Uncaught exception (fatal)" );
        cerr << "Uncaught exception (fatal)"  << endl;
    }

    return EXIT_SUCCESS;

}

