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
#include "application.hpp"

#include "revision.hpp"
#include "version.hpp"

#include <iostream>
#include <vector>

#include <boost/filesystem.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace fs = boost::filesystem;
namespace bpt = boost::property_tree;
namespace bpo = boost::program_options;

using namespace wfwfs;
using namespace std;

string Application::executableName;
uint8_t Application::logLevel(0);

void Application::getOptions( bpo::options_description& options, const string& name ) {

    bpo::options_description general( "General Options" );
    general.add_options()
    ( "version,V", "Print version information and quit." )
    ( "copyright", "Print copyright information and quit." )
    ( "help,h", "Show command line options and quit." )
    ( "verbosity", bpo::value< int >(), "Specify verbosity level (0-8, 0 means no output)."
      " The environment variable WFWFS_VERBOSITY will be used as default if it exists." )
    ( "verbose,v", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "1" ), "" )
      ->composing(), "More output. (ignored if --verbosity is specified)" )
    ( "quiet,q", bpo::value<vector<string>>()->implicit_value( vector<string>( 1, "-1" ), "" )
      ->composing(), "Less output. (ignored if --verbosity is specified)" )

    ;

    bpo::options_description config( "Configuration" );
    config.add_options()
    ( "appname", bpo::value<string>()->default_value( name ),
      "Name used to identify this application." )
    ;

    options.add( general ).add( config );
}


pair<string, string> Application::customParser( const string& s ) { // custom parser to handle multiple -q/-v flags (e.g. -vvvv)

    if( s.find( "-v" ) == 0 || s.find( "--verbose" ) == 0 ) { //
        int count = std::count( s.begin(), s.end(), 'v' );
        while( count-- ) logLevel++;
    }
    else if( s.find( "-q" ) == 0 || s.find( "--quiet" ) == 0 ) { //
        int count = std::count( s.begin(), s.end(), 'q' );
        while( count-- ) logLevel--;
    }
    return make_pair( string(), string() );                 // no need to return anything, we handle the verbosity directly.
}


bpo::options_description& Application::parseCmdLine( int argc, char* argv[],  bpo::variables_map& vm,
                                                    bpo::options_description* programOptions,
                                                    bpo::positional_options_description *positionalOptions,
                                                    parserFunction custom_parser ) {
    static bpo::options_description all;
    
    try {
        if( programOptions ) {
            all.add( *programOptions );
        }

        fs::path tmpPath = fs::path(string(argv[0])).filename();
        Application::executableName = tmpPath.string();
        getOptions( all, Application::executableName );
        
        bpo::command_line_parser parser( argc, argv );
        parser.options( all );
        parser.extra_parser( customParser );

        if( positionalOptions ){
            parser.allow_unregistered().positional(*positionalOptions);
        }
        bpo::store( parser.run(), vm );

        if( custom_parser ) {
            custom_parser( all, vm );
        }

        // If e.g. --help was specified, just dump output and exit.
        checkGeneralOptions( all, vm );

        //bpo::store( bpo::parse_environment( all, Logger::environmentMap ), vm );
        
        //vm.notify();
        //bpo::variables_map::notify(vm);
    } catch( const exception& e ) {
        cout << "Failed to parse command-line. Reason: " << e.what() << endl;
        vm.insert( std::make_pair( "help", bpo::variable_value() ) );
    }
    
    return all;

}


void Application::checkGeneralOptions( bpo::options_description& desc, bpo::variables_map& vm ) {

    if( vm.count( "help" ) ) {
        cout << desc << endl;
        exit( 0 );
    }
    if( vm.count( "version" ) ) {
        if( logLevel ) {
            cout << "Version:  " << getLongVersionString() << endl;
            cout << "Commited: " << commitTime << endl;
            cout << "Compiled: " << buildTime << endl;
        } else {
            cout << getVersionString() << endl;
        }
        exit( 0 );
    }
    if( vm.count( "copyright" ) ) {
        cout << "Not implemented\n";
        exit( 0 );
    }
    if( vm.count( "tutorial" ) ) {
        cout << "Not implemented\n";
        exit( 0 );
    }
    if( vm.count( "sample" ) ) {
        cout << "Not implemented\n";
        exit( 0 );
    }

}


Application::Application( bpo::variables_map& vm, RunMode rm ) : runMode(rm), returnValue(0) {

    if( vm.count("settings") ) {
        settingsFile = vm["settings"].as<string>();
        if( fs::is_regular(settingsFile) ) {
            cout << "Loading file \"" <<  settingsFile << "\"" << endl;
            bpt::read_info( settingsFile, propTree );
        } else {
            cerr << "Failed to load file \"" << settingsFile << "\", starting with default settings." << endl;
            exit( 0 );
        }
    }
    applicationName = vm["appname"].as<string>();

}


Application::~Application( void ) {

}


string Application::getName( void ) const {
    return applicationName;
}


int Application::run( void ) {

    while( doWork() && !runMode ) {       // keep calling doWork() until it returns false, or shouldStop is raised

    }
    
    switch(runMode) {
        case RESET: throw ResetException();           // this will cause a complete reset (i.e. creation of a new Application instance)
        case RESTART: throw RestartException();           // this will cause a complete reset (i.e. creation of a new Application instance)
        //case EXIT:  throw KillException();            // exits the loop in main()
        default: ;
    }

    // No work left.
    return returnValue;
}


