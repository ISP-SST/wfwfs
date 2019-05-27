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

#include "arraystats.hpp"
#include "cell.hpp"
#include "filecam.hpp"
#include "fitswriter.hpp"
#include "pleoracam.hpp"
#include "recursepath.hpp"
#include "revision.hpp"
#include "stringutil.hpp"
#include "translators.hpp"

#include <atomic>
#include <chrono>
#include <cstring>
#include <functional>
#include <future>
#include <map>
#include <numeric>
#include <sys/resource.h> 
#include <queue> 

#include <boost/asio/time_traits.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <boost/property_tree/info_parser.hpp>

namespace bfs = boost::filesystem;
namespace bpo = boost::program_options;
namespace bpx = boost::posix_time;

using boost::algorithm::iequals;

using namespace wfwfs;
using namespace std;

mutex Daemon::globalMutex;
map<TcpConnection::Ptr, set<string>, PtrCompare<TcpConnection>> Daemon::subscriptions;


namespace {
    
    const set<string> auto_subscribe = { "exposure", "gain", "interval", "mode", "state" };
    
    bool draw_cells = false;
    
    atomic<int> nThreads(0);
    
    size_t histo_skip;
    string filename_prefix;
    string camera_id;
    vector<string> meta_cards;
    
    struct gui_meta {  // some volatile metadata that might change while daemon is running.
        gui_meta() : use(false) {};
        string observer;
        string target;
        string date_obs;
        int compress;
        bool use;
        string comments;
        size_t last_frame;      // just to keep track of framenumbers over consecutive bursts.
        operator std::string() const {
            return observer + ", " + target + ", " + date_obs + ", " + to_string(compress) + ", :" + comments;
        }
    } gui_info;

    bool show_timing(false);
    
    Array<uint8_t> cell_mask;
    
    struct CalibID {
        CalibID( int w, int h, int d, float e, float g ) : width(w), height(h), depth(d), exposure(e), gain(g) {}
        bool operator<(const CalibID& rhs) const {
            if( width != rhs.width ) return (width < rhs.width);
            if( height != rhs.height ) return (height < rhs.height);
            if( depth != rhs.depth ) return (depth < rhs.depth);
            if( fabs(gain-rhs.gain)/(gain+rhs.gain) > 0.1 ) return (gain < rhs.gain);
            if( fabs(exposure-rhs.exposure)/(exposure+rhs.exposure) > 0.1 ) return (exposure < rhs.exposure);
            return false;
        }
        // parameters to identify a unique camera configuration
        int width, height, depth;
        float exposure, gain;
    };
    struct CalibFile {
        CalibFile( const string& fn, const bpx::ptime& t ) : temp(0.0/0.0), rms(0.0), filename(canonical(bfs::path(fn))), time(t) {}
        bool operator<(const CalibFile& rhs) const {
            if( time != rhs.time ) return (time < rhs.time);
            return (filename < rhs.filename);
        }
        mutable float temp, rms;
        bfs::path filename;
        bpx::ptime time;
    };
    struct CalibSets {
        set<CalibFile> dd, ff;
    };
    map<CalibID,CalibSets> calib_archive;
    CalibID current_calib( 0, 0, 0, 0, 0 );
    
    CalibSets& getCalibSet( CalibID id ) {
        auto ret = calib_archive.emplace( id, CalibSets() );
        return ret.first->second;
    }
    
    CalibSets& getCalibSet( int w, int h, int d, float e, float g ) {
        return getCalibSet( CalibID(w,h,d,e,g) );
    }
    
    
    bool parse_path( const boost::filesystem::path& p ) {
        
        if( bfs::is_directory(p) ) return true;    // keep recursing
        
        if( bfs::is_regular(p) || bfs::is_symlink(p) ) {  // load and parse file
            try {
                string filename = p.filename().string();
                shared_ptr<Fits> hdr( new Fits( p.string() ) );
                const vector<string>& cards = hdr->primaryHDU.cards;
                int datamax = Fits::getValue<int>( cards, "DATAMAX" );
                float g = Fits::getValue<float>( cards, "DETGAIN" );
                float e = Fits::getValue<float>( cards, "TEXPOSUR" );
                size_t w = hdr->dimSize(0);
                size_t h = hdr->dimSize(1);
                bpx::ptime t = Fits::getValue<bpx::ptime>( cards, "DATE-OBS" );
                CalibSets& cs = getCalibSet( w, h, datamax, e, g );
                
                CalibFile cf( p.string(), t );
                
                if( filename.find("dark_") == 0 ) {
                    auto ret = cs.dd.emplace( cf );
                    if( !ret.second ) {
                        // do we want/need to modify anything if there is a collision ?
                    }
                } else if ( filename.find("flat_") == 0) {
                    auto ret = cs.ff.emplace( cf );
                    if( !ret.second ) { 
                        // do we want/need to modify anything if there is a collision ?
                    }
                }
            } catch( ... ) {}
        }
        
        return false;
    }
    
    
    void load_calibs( bfs::path p ) {
        calib_archive.clear();
        RecursePath rp( p, parse_path );
    }
    
    
    bool replace(std::string& str, const std::string& from, const std::string& to) {
        size_t start_pos = str.find(from);
        bool found = false;
        while(start_pos != std::string::npos) {
            str.replace(start_pos, from.length(), to);
            found = true;
            start_pos = str.find(from);
        }
        return found;
    }
    
    string cleanupCfg( string filename ) {

        std::ifstream in( filename, std::ios::in | std::ios::binary );
        if( !in ) {
            throw (errno);
        }

        string cfg;
        string line;
        while( in.good() ) {
            std::getline( in, line );
            replace( line, "=", " " );
            boost::trim(line);
            size_t found = line.find(" ");
            if( found && (found != string::npos) ) {
                string key = line.substr( 0, found );
                string value = line.substr( found, string::npos );
                boost::trim( value );
                if( value.find(" ") != string::npos ) value = "\"" + value + "\"";  // quote values with spaces
                line = key + " " + value;
            }
            cfg += line + "\n";
        }
        in.close();

        if( *(cfg.rbegin()) != '\n') cfg += "\n";       // possibly add a newline at end
        
        return cfg;
    }

    
    void add_meta( string key, string value ) {
        
        boost::trim(key);
        if( key.empty() ) {
            return;
        }
        
        string comment;
        string card;
        size_t pos = value.find('/');
        if( pos != string::npos ) {
            comment = value.substr( pos+1, string::npos );
            boost::trim(comment);
            value = value.substr(0,pos);
            boost::trim(value);
        }
        
        try {
            int val = boost::lexical_cast<int>(value);
            card = Fits::makeCard( key, val, comment );
        } catch( const boost::bad_lexical_cast& ) { }

        if( card.empty() ) {
            try {
                float val = boost::lexical_cast<float>(value);
                card = Fits::makeCard( key, val, comment );
            } catch( const boost::bad_lexical_cast& ) { }
        }
        
        if( card.empty() ) {
            card = Fits::makeCard( key, value, comment );
        }

        Fits::updateCard( meta_cards, card );
        
    }
    
    void add_default_meta( void ) {
        
        FitsWriter::clear_gmeta();
        FitsWriter::update_gmeta( "SOLARNET", 0.5 );
        FitsWriter::update_gmeta( "OBS_HDU", 1 );
        FitsWriter::update_gmeta( "EXTEND", true );
        FitsWriter::update_gmeta( "TELESCOP", "Swedish 1-meter Solar Telescope" );
        FitsWriter::update_gmeta( "INSTRUME", "WFWFS" );
        FitsWriter::update_gmeta( "OBSERVER", "SST Observer" );
        FitsWriter::update_gmeta( "ORIGIN", "Institute for Solar Physics" );
        FitsWriter::update_gmeta( "OBJECT", "Sun" );
        FitsWriter::update_gmeta( "DETECTOR", "cam0" );
        
        add_meta( "SOLARNET", "0.5" );
        add_meta( "OBS_HDU", "1" );
        add_meta( "TELESCOP", "Swedish 1-meter Solar Telescope" );
        add_meta( "INSTRUME", "WFWFS" );
        add_meta( "OBSERVER", "SST Observer" );
        add_meta( "ORIGIN", "Institute for Solar Physics" );
        add_meta( "OBJECT", "Sun" );
        add_meta( "DETECTOR", "cam0" );

/*        
        
        Fits::updateCard( hdr->primaryHDU.cards, Fits::makeCard( "CADENCE", "interval", "[s]" ) );
        Fits::updateCard( hdr->primaryHDU.cards, Fits::makeCard( "DETOFFS", 0, "[counts] or camera specific unit" ) );
        float temp = cam->get_temp();
        if( temp ) {
            Fits::updateCard( hdr->primaryHDU.cards, Fits::makeCard( "DETTEMP", temp,"[C] Current temperature of the detector" ) );
        }
        Fits::updateCard( hdr->primaryHDU.cards, Fits::makeCard( "FIELD", "name" ) );
        Fits::updateCard( hdr->primaryHDU.cards, Fits::makeCard( "DATAMIN", 0 ) );
*/
    }
    
    Daemon* currentD(nullptr);
    
}


Daemon& Daemon::get( void ) {
    
    lock_guard<mutex> lock( globalMutex );
    return *currentD;
    
}


void Daemon::set( Daemon* d ) {
    
    lock_guard<mutex> lock( globalMutex );
    currentD = d;
    
}


Daemon::Daemon( bpo::variables_map& vm ) : Application( vm, LOOP ), settings( vm ), has_light(false), timer( ioService ), ddPtr(nullptr), ggPtr(nullptr) {
    
    set( this );
    
    uint16_t nT = settings["threads"].as<uint16_t>();
    if( nT ) {
        // TODO
    }

    add_default_meta();
    
    init();
    
    uint16_t port = settings["port"].as<uint16_t>();
    if( port < 1024 ) {
        cerr << "Daemon:  using a port < 1024 requires root permissions, which this program should *not* have." << endl;
        stop();
        return;
    }

    try {
        server.reset( new TcpServer( ioService, port ) );
        server->setCallback( bind( &Daemon::connected, this, std::placeholders::_1 ) );
        server->accept();
    } catch ( const exception& e ) {
        server.reset();
        stop();
        throw;
    }
    

}


Daemon::~Daemon( void ) {
    stop();
}


void Daemon::reset( void ) {
    runMode = RESET;
    ioService.stop();
    pool.interrupt_all();
}


void Daemon::restart( void ) {
    runMode = RESTART;
    ioService.stop();
    pool.interrupt_all();
}


void Daemon::stop( void ) {
    runMode = EXIT;
    ioService.stop();
    pool.interrupt_all();
    if( cam ) cam->stop();
}


void Daemon::init(void) {

    if( !settings.count("config") ) {
        throw runtime_error( "No camera configuration supplied." );
    }

    string cfg = settings["config"].as<string>();
    if( !bfs::is_regular_file( cfg ) ) {
        throw runtime_error( string("No such file: ") + cfg );
    }

    stringstream filteredCfg;
    string tmpCfg = cleanupCfg( cfg );
    filteredCfg.write( tmpCfg.c_str(), tmpCfg.size() );

    boost::property_tree::ptree cfg_ptree;
    try {
        boost::property_tree::read_info( filteredCfg, cfg_ptree );
    } catch ( const exception& e ) {
        string msg = "Failed to parse config-file \"" + cfg + "\"\nboost reports: ";
        msg += e.what();
        throw runtime_error( msg.c_str() );
    } catch ( ... ) {
        string msg = "Failed to parse config-file \"" + cfg + "\".";
        throw runtime_error( msg.c_str() );
    }
    
    parsePropertyTree( cfg_ptree );
    
    FitsWriter::update_gmeta( "DETECTOR", cam->info.id );
    FitsWriter::update_gmeta( "DETMODEL", cam->info.model );
    FitsWriter::update_gmeta( "DETFIRM", cam->info.firmware_version );
    
    add_meta( "DETECTOR", cam->info.id );
    add_meta( "DETMODEL", cam->info.model );
    add_meta( "DETFIRM", cam->info.firmware_version );
    
    size_t frame_queue = cfg_ptree.get<size_t>( "frame_queue", 10 );
    fqueue.resize( cam->cfg.width, cam->cfg.height, frame_queue, cam->cfg.depth );
    
    updateCalibID();
    
    dd.resize( fqueue.width, fqueue.height );
    ff.resize( fqueue.width, fqueue.height );
    gg.resize( fqueue.width, fqueue.height );

    dd = 0.0;
    ff = 1.0;
    
    bfs::path tmp(outputdir+"/calib/");
    if( bfs::exists( tmp ) ) {
        try {
            load_calibs( tmp );
        } catch( ... ) {}
    } else {
        maybeCreateDir( tmp );
    }
    
    tmp = bfs::path(outputdir+"/calib/dark.fits");
    if( bfs::exists( tmp ) ) {
        try {
            readFile( tmp.string(), dd );
        } catch( ... ) {}
    }
    
    tmp = bfs::path(outputdir+"/calib/flat.fits");
    if( bfs::exists( tmp ) ) {
        try {
            readFile( tmp.string(), ff );
        } catch( ... ) {}
    }
    
    ggPtr = gg.get();
    ddPtr = dd.get();
    
    update_gain();
    cell_mask.resize( fqueue.width, fqueue.height );
    seeing.find_nominal_gridpoints( cam->get_size() );
    seeing.draw_cells( cell_mask );
    
    start_cam();
//    seeing.start_logs();
    
}


void Daemon::parsePropertyTree( boost::property_tree::ptree& cfg_ptree ) {
    
    outputdir = cfg_ptree.get<string>( "outputdir", "/data/" );
    
    histo_skip = cfg_ptree.get<size_t>( "histo_skip", 1 );

    for( auto& node : cfg_ptree ) {
        if( iequals( node.first, "META" ) ) {
            for( auto& meta : node.second ) {
                string value = meta.second.get_value<std::string>();
                string key = meta.first;
                add_meta( key, value );
            }
        } else if( iequals( node.first, "CAMERA" ) ) {
            string camType = node.second.get<string>( "camera_type", "" );
            if( boost::iequals( camType, "pleora" ) ) {
#ifdef WFWFS_WITH_PLEORA
                cam.reset( new PleoraCam( node.second ) );
#else
                throw runtime_error( "Support for Pleora camera not enabled/detected." );
#endif
            } else if( bfs::exists( camType ) ) {
                cam.reset( new FileCam( node.second ) );
            } else {
                throw runtime_error("camera_type not valid: \""+camType+"\"");
            }
        }
    }

    seeing.parsePropertyTree( cfg_ptree );
    cam->init();

    
}


void Daemon::queue_frame( const uint8_t* data, bpx::ptime ts ) {
    
    static double avg_interval(1.0);
    static bpx::ptime previous_frame( bpx::not_a_date_time );
    boost::this_thread::interruption_point();

    if( !previous_frame.is_not_a_date_time() ) {
        if( ts <= previous_frame ) {
            return;
        }
        bpx::time_duration elapsed(ts - previous_frame);
        double d = elapsed.total_microseconds()*1E-6;
        avg_interval *= 0.9;
        avg_interval += d * 0.1;
    }

    Frame& f = fqueue.getEmpty();
    f.timestamp = ts;
    if( show_timing ) {
        bpx::ptime start_time = bpx::microsec_clock::universal_time();
        cam->copyRawData( data, f.data, fqueue.depth );
        bpx::ptime read_time = bpx::microsec_clock::universal_time();
        copy_cell_data( f );
        bpx::ptime cell_time = bpx::microsec_clock::universal_time();
        fqueue.queue( f );
        bpx::ptime queue_time = bpx::microsec_clock::universal_time();
        if( has_light ) {
            seeing.process( avg_interval );
        }
        bpx::ptime proc_time = bpx::microsec_clock::universal_time();
        bpx::time_duration rt(read_time - start_time);
        bpx::time_duration ct(cell_time - read_time);
        bpx::time_duration qt(queue_time - cell_time);
        bpx::time_duration pt(proc_time - queue_time);
        bpx::time_duration tt(proc_time - start_time);
        string msg = boost::str( boost::format("\r  rt: %06d  ct: %06d  qt: %06d  pt: %06d  tot: %s  fps: %3.1lf (limit: %3.1f)")
                                 % rt.total_microseconds() % ct.total_microseconds() % qt.total_microseconds()
                                 % pt.total_microseconds() % to_simple_string(tt) % (1.0/avg_interval) % (1E6/tt.total_microseconds())
                               );
        cout << msg << flush;
    } else {
        cam->copyRawData( data, f.data, fqueue.depth );
        copy_cell_data( f );
        fqueue.queue( f );
        //seeing.process( ioService );
        if( has_light ) {
            ioService.post( bind(&Seeing::process, &seeing, avg_interval ) );
        }
    }
    
    if( !has_light ) {
        seeing.clear_frame_data();      // Discard frame-data that should not be processed.
    }
    
    previous_frame = ts;

}


void Daemon::get_frame( TcpConnection::Ptr connection, int x1, int y1, int x2, int y2, int scale, bool df, size_t fsel, bool do_histo ) {
    
    if( !cam ) {
        throw runtime_error( "get_frame: No camera running." );
    }
    
    x1 = max( x1, 0 );
    y1 = max( y1, 0 );
    scale = max( scale, 1 );

    x2 = min<int>( x2, fqueue.width/scale );
    y2 = min<int>( y2, fqueue.height/scale );

    int frameSize = (x2-x1)*(y2-y1);
    if( frameSize < 1 ) {
        throw runtime_error( "get_frame: illegal frame-size: ["+to_string(x2-x1)+","+to_string(y2-y1)+"] pixels." );
    }

    size_t imgSize = frameSize*(fqueue.depth <= 8 ? 1 : 2);
    
    string reply = boost::str( boost::format("OK image %lu %d %d %d %d %d") % imgSize % x1 % y1 % x2 % y2 % scale );
    
    Frame& f = fqueue.getFrame( 0 );
    
    size_t histSize(0);
    size_t histPixels(0);
    int histMask(0);
    if( !f.hist ) do_histo = false;
    
    if( do_histo ) {
        reply += " histogram";
        histSize = (1 << fqueue.depth);
        histMask = histSize-1;
        std::fill_n( f.hist, histSize, 0 );
    }
    
    size_t blockSize = imgSize + histSize * sizeof(uint32_t);
    
    unique_ptr<uint8_t[]> buf( new uint8_t[blockSize] );
    fill_n( buf.get(), blockSize, 0 );
    uint8_t* cmask = cell_mask.get();
    
    int maxval = (1 << fqueue.depth) - 1;
    
    if( fqueue.depth > 8 ) {
        uint16_t *in = reinterpret_cast<uint16_t*>( f.data );
        uint16_t *p = reinterpret_cast<uint16_t*>( buf.get() );
        for(int y = y1 * scale; y < y2 * scale; y += scale) {
            for(int x = x1 * scale; x < x2 * scale; x += scale) {
                size_t o = y * fqueue.width + x;
                if( !draw_cells || !cmask[o] ) {
                    if( df ) {
                        *p = dg_correct( in, o, ddPtr, ggPtr, maxval );
                    } else {
                        *p = in[o];
                    }
                }
                if( do_histo ) {
                    f.hist[ *p&histMask ]++;
                    histPixels++;
                }
                p++;
            }
        }
    } else {
        uint8_t *in = (uint8_t *)f.data;
        uint8_t *p = buf.get();
        for(int y = y1 * scale; y < y2 * scale; y += scale) {
            for(int x = x1 * scale; x < x2 * scale; x += scale) {
                size_t o = y * fqueue.width + x;
                if( !draw_cells || !cmask[o] ) {  
                    if( df ) {
                        *p = dg_correct( in, o, ddPtr, ggPtr, maxval );
                    } else {
                        *p = in[o];
                    }
                }
                if( do_histo ) {
                    f.hist[ *p&histMask ]++;
                    histPixels++;
                }
                p++;
            }
        }
    }

    if( histSize ) {
        uint32_t* histPtr = reinterpret_cast<uint32_t*>( buf.get()+imgSize );
        copy_n( f.hist, histSize, histPtr );
    }

    connection->writeline( reply );
    connection->syncWrite( buf.get(), blockSize );

}


string Daemon::make_filename( int scannum, long framenum, const string& state ) {
    
    bpx::ptime timestamp = bpx::second_clock::universal_time();
    if( gui_info.date_obs.size() < 11 || gui_info.date_obs[10] != 'T' ) {      
        gui_info.date_obs = to_iso_extended_string( timestamp );
    }
    string date = gui_info.date_obs.substr( 0, 10 );
    string time = gui_info.date_obs.substr( 11, 8 );

    // SST:    ${outputdir}/${date}/${instrument}/08:36:03/${detector}/sst_camXX_${scan}_${frame}_state.fits
    // WFWFS:  ${outputdir}/${date}/${time}/wfwfs_${scan}_${frame}.fits
    // or WFWFS:  ${outputdir}/${date}/wfwfs_${timestamp}.fits
    bfs::path filePath(outputdir);
    filePath /= date;
    filePath /= time;
    maybeCreateDir( filePath );

    string filename = filename_prefix;
    filename += boost::str( boost::format("_%05ld_%07ld") % scannum % framenum );
    //filename += "_" + to_iso_extended_string( timestamp );
    if( !state.empty() ) {
        filename += "_" + state;
    }
    filename += ".fits";
    filePath /= filename;
    
    return filePath.string();

}


string Daemon::make_filename( std::string tpl, int file_number, int frame_number ) {
    
    if( (frame_number >= 0) && contains( tpl, "%framenumber%", true ) ) {
        string tmp = boost::str( boost::format("%07ld") % frame_number );
        tpl = replace_n( tpl, "%framenumber%", tmp, 5 );
        tpl = replace_n( tpl, "%FRAMENUMBER%", tmp, 5 );
    }
    
    if( (file_number >= 0) && contains( tpl, "%filenumber%", true ) ) {
        string tmp = boost::str( boost::format("%05ld") % file_number );
        tpl = replace_n( tpl, "%filenumber%", tmp, 5 );
        tpl = replace_n( tpl, "%FILENUMBER%", tmp, 5 );
    }
    
    if( contains( tpl, "%date%", true ) || contains( tpl, "%time%", true ) ) {
        bpx::ptime timestamp = bpx::second_clock::universal_time();
        string timestamp_s = to_iso_extended_string( timestamp );
        string date = timestamp_s.substr( 0, 10 );
        string time = timestamp_s.substr( 11, 8 );
        tpl = replace_n( tpl, "%date%", date, 5 );
        tpl = replace_n( tpl, "%DATE%", date, 5 );
        tpl = replace_n( tpl, "%time%", time, 5 );
        tpl = replace_n( tpl, "%TIME%", time, 5 );
    }
    
    bfs::path filePath;
    if( isRelative(tpl) ) {
        filePath = bfs::path(outputdir);
        filePath /= tpl;
    } else {
        filePath = bfs::path(tpl);
    }
    
    return filePath.string();

}


void Daemon::save_fits( string filename_template, int nframes, uint32_t frames_per_file, bool compress,
                        size_t first_frame, string acc_filename ) {

    int nthreads = std::thread::hardware_concurrency();
 
    // Generate meta for some volatile info that might change between calls
    vector<string> cards;
    Fits::updateCard( cards, Fits::makeCard( "DATAMIN", 0 ) );
    Fits::updateCard( cards, Fits::makeCard( "DATAMAX", (1<<fqueue.depth)-1 ) );
    Fits::updateCard( cards, Fits::makeCard( "XPOSURE", cam->get_exposure(), "[s] Exposure time" ) );
    Fits::updateCard( cards, Fits::makeCard( "DETGAIN", cam->get_gain(), "[dB] or camera specific unit" ) );
    Fits::updateCard( cards, Fits::makeCard( "CADENCE", cam->get_interval(), "[s]" ) );
    Fits::updateCard( cards, Fits::makeCard( "DETOFFS", 0, "[counts] or camera specific unit" ) );
    float temp = cam->get_temp();
    if( isfinite(temp) ) {
        Fits::updateCard( cards, Fits::makeCard( "DETTEMP", temp,"[C] Current temperature of the detector" ) );
    }

    if( gui_info.use ) {
        if( !gui_info.observer.empty() ) FitsWriter::update_gmeta( "OBSERVER", gui_info.observer );
        if( !gui_info.target.empty() ) FitsWriter::update_gmeta( "OBJECT", gui_info.target );
        string comments = gui_info.comments;
        while( !comments.empty() ) {
            string comment = popword( comments, ";" );
            if( !comment.empty() ) {
                Fits::insertCard( cards, Fits::makeCard( "COMMENT", comment ) );
            }
        }
    }
    
    shared_ptr<FitsWriter> fw;
    fw.reset( new FitsWriter( fqueue, acc_filename, nthreads, compress ) );
    fw->save_meta( cards );
    
    int frame_count(0);
    int file_count(0);
    broadcast( "state", "OK state burst" );
    while( frame_count < nframes ) {
        string this_fn = make_filename( filename_template, file_count, frame_count );
        int nf = std::min<int>( nframes-frame_count, frames_per_file );
        fw->save( this_fn, first_frame, nf );
        file_count++;
        frame_count += frames_per_file;
    }
    broadcast( "state", "OK state ready" );

}


void Daemon::copy_cell_data( Frame& f  ) {
 
    if( fqueue.depth > 8 ) {
        seeing.copy_cell_data<uint16_t>(f, ddPtr, ggPtr, fqueue.width, (1<<fqueue.depth)-1);
    } else {
        seeing.copy_cell_data<uint8_t>(f, ddPtr, ggPtr, fqueue.width, (1<<fqueue.depth)-1);
    }

}


void Daemon::threadLoop( void ) {
    
    nThreads++;

    while( runMode == LOOP ) {
        try {
            boost::this_thread::interruption_point();
            ioService.run();
        } catch( const ThreadExit& e ) {
            break;
        } catch( const boost::thread_interrupted& ) {
            cout << "Daemon: Thread interrupted." << endl;
            break;
        } catch( exception& e ) {
            cerr << "Exception in thread: " << e.what() << endl;
        } catch( ... ) {
            cerr << "Unhandled exception in thread." << endl;
        }
    }

    nThreads--;
    
    // pool.remove_thread( tmap[boost::this_thread::get_id()] );
    

}


void Daemon::maintenance( void ) {


    timer.expires_from_now( boost::posix_time::seconds( 5 ) );
    
    static float last_intensity(0.05);
    
    float avg_intensity = seeing.get_avg_intensity();
    float ratio = avg_intensity/last_intensity;

    //if( fabs((avg_intensity-previous_intensity)/previous_intensity) > 0.3 ) { FIXME
    if( ratio < 0.2 ) {
        light( false );
    } else if( ratio > 5 ) {
        light( true );
    }

    last_intensity = avg_intensity;
    
    seeing.maintenance();
    
    timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
    
}


bool Daemon::doWork( void ) {

    try {
        
        // Add some threads for the async work.

        timer.expires_from_now( boost::posix_time::seconds( 5 ) );
        timer.async_wait( boost::bind( &Daemon::maintenance, this ) );
        
        setThreads( thread::hardware_concurrency() );

        while( !runMode ) {
            std::this_thread::sleep_for( std::chrono::milliseconds(100) );
        }
        
        // the io_service will keep running/blocking until stop is called, then wait for the threads to make a clean exit.
        pool.join_all();

    } catch( ... ) {
        throw;
    }

    return true;

}


void Daemon::setThreads( int nT ) {

    if( nT < 1 ) nT = 1;        // daemon will not respond to connections if it doesn't have threads
    int diff = nT - nThreads;
    if( diff < -nThreads ) diff = -nThreads;
    
    // FIXME  dynamically changing the number of sys-threads will cause threads to get stuck on pool.create_thread/remove_thread
    //diff = 50-nSysThreads;  // force 50 threads until bug is fixed
    
    if( diff > 0 ) {
        for( int i=0; i < diff; ++i ) {
            pool.create_thread( std::bind( &Daemon::threadLoop, this ) );
        }
    } else while( diff++ < 0 ) {
        ioService.post( [](){ throw ThreadExit(); } );
    }
    
}


void Daemon::connect( Host& host, TcpConnection::Ptr& conn ) {
    
//     if( host.connectName.empty() ) {
//         cerr << "Attempting to connect without a hostname." << endl;
//         return;
//     }
//     
//     if( !conn ) {
//         conn = TcpConnection::newPtr( ioService );
//     }
//     
//     try {
//         auto test WFWFS_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, will throw if not connected.
//         return;
//     } catch ( ... ) {
//         // if we get here the socket is disconnected, continue to reconnect below.
//     }
//     
//     if( conn->socket().is_open() ) {
//         conn->socket().close();
//     }
// 
// #ifdef DEBUG_
//     cout << "Attempting to connect to " << host.connectName << ":" << host.connectPort << endl;
// #endif
//     
//     try {
//         conn->connect( host.connectName, to_string(host.connectPort) );
//         if( conn->socket().is_open() ) {
// 
//             Command cmd;
// 
//             myInfo.info.peerType |= Host::TP_WORKER;
//             *conn << CMD_CONNECT;
//             *conn >> cmd;
//             if( cmd == CMD_AUTH ) {
//                 // implement
//             }
//             if( cmd == CMD_CFG ) {  // handshake requested
//                 *conn << myInfo.info;
//                 *conn >> host;
//                 *conn >> cmd;       // ok or err
//             }
//             if( cmd != CMD_OK ) {
//                 cerr << "Handshake with master failed  (server replied: " << cmd << ")" << endl;
//                 conn->socket().close();
//                 //myInfo.info.peerType &= ~Host::TP_WORKER;
//             }
//         }
// 
//     } catch ( const std::exception& e ) {
//         cerr << "Failed to connect to " << host.connectName << ":" << host.connectPort
//                 << "  reason:" << e.what() << endl;
//     } catch ( ... ) {
//         cerr << "Unhandled exception when connecting: " << __LINE__ << endl;
//     }
}


void Daemon::connected( TcpConnection::Ptr conn ) {
    
    if( !conn ) return;

    try {
        Host h;
        h.address = conn->getRemoteIP();
        h.port = conn->getRemotePort();
        addConnection( h, conn );
        conn->setCallback( bind( &Daemon::onMessage, this, std::placeholders::_1 ) );
        conn->setErrorCallback( bind( &Daemon::removeConnection, this, std::placeholders::_1 ) );
        conn->idle();
        return;
    }
    catch( const exception& e ) {
        cerr << "connected() Failed to process new connection. Reason: " << e.what() << endl;
    } catch( ... ) {
        cerr << "Daemon::connected() Unhandled exception." << endl;
    }
    conn->socket().close();
}


bool Daemon::processCmd( TcpConnection::Ptr conn, const string& cmd ) {
    
    if( cmd.empty() ) return true;
    
    bool disconnect( false );
    string line = cmd;
    string replyStr = "";
    string tag = "";

    string command = popword( line );

    string orig_line = line;       // backup since cam->parse_cmd() is destructive
    if( command == "quit" || command == "exit" ) {
        replyStr = "OK :Bye!";
        disconnect = true;
    } else if( command == "die" || command == "shutdown" ) {
        replyStr = "OK :Bye!";
        die();
        disconnect = true;
    } else if( command == "reset") {
        replyStr = "OK :Bye!";
        reset();
        disconnect = true;
    } else if( command == "restart" ) {
        replyStr = "OK :Bye!";
        restart();
        disconnect = true;
    } else if( command == "run" ) {
        replyStr = "OK run";
        run();
    } else if( command == "pause" ) {
        replyStr = "OK pause";
        pause();
    } else if( command == "segfault" ) {
        ioService.post( [](){
            std::this_thread::sleep_for( std::chrono::seconds(1) );  // delay for reply to be sent back
            *(int*)(8) = 0; }
        );
    } else if( command == "thumbnail" ) {
        thumbnail( conn );
    } else if( command == "test_threads" ) {
        int nT = pop<int>( line );
        int time = pop<int>( line );
        if( nT < 1 ) nT = 1;
        if( time < 1 ) time = 1;
        std::thread( std::bind(&Daemon::test_threads, this, nT, time) ).detach();
    } else if( command == "draw_cells" ) {
        draw_cells = !draw_cells;
        seeing.draw_cells( cell_mask );
        replyStr = "OK draw_cells "+to_string(draw_cells);
    } else if( command == "start" ) {
        string what = popword( line );
        if( what == "cam" ) {
            start_cam();
            replyStr = "OK start cam";
        } else if( what == "logs" ) {
            seeing.start_logs();
            replyStr = "OK start logs";
        } else if( what.empty() ) {
            seeing.start();
            replyStr = "OK start";
        } else {  // unrecognized
            if( !what.empty() ) replyStr = "start Huh?";
        }

    } else if( command == "stop" ) {
        string what = popword( line );
        if( what == "cam" ) {
            stop_cam();
            replyStr = "OK stop cam";
        } else if( what == "logs" ) {
            seeing.stop_logs();
            replyStr = "OK stop logs";
        } else if( what.empty() ) {
            seeing.stop();
            replyStr = "OK stop";
        } else {  // unrecognized
            if( !what.empty() ) replyStr = "stop Huh?";
        }

    } else if( command == "set" ) {
        string what = popword( line );
        tag = what;
        if( !(cam && cam->parse_cmd( command, what, line, replyStr )) ) {
            line = orig_line;       // cam->parse_cmd() is destructive
            what = popword( line );
            if( what == "histo_skip" ) {
                histo_skip = pop<size_t>( line );
                replyStr = boost::str( boost::format("OK histo_skip %lu") % histo_skip );
            } else if( what == "filename" ) {
                if( !line.empty() && line[0] == ':' ) {
                    line.erase( line.begin() );
                    filename_prefix = popword( line, "\t\n" );
                } else filename_prefix = popword( line );
                //replyStr = "OK filename " + filename_prefix;
            } else if( what == "cell" ) {
                int cid = pop<int>( line );
                int y = pop<int>( line );
                int x = pop<int>( line );
//                 lock_guard<mutex> lock(cellMutex);
//                 cells[cid].pos = PointI(y,x);
//                 make_sh_mask();
                replyStr = boost::str( boost::format("OK cell %d %d %d") % cid % x % y );
            } else if( what == "fits" ) {
                gui_info.observer = popword( line, "," );
                gui_info.target = popword( line, ", " );
                gui_info.date_obs = popword( line, ", " );
                gui_info.compress = pop<int>( line, ", " );
                gui_info.comments = popword(line);
                gui_info.last_frame = 0;
                gui_info.use = true;
                replyStr = "OK fits " + (string)gui_info;
            } else if( what == "threads" ) {
                int nT = pop<int>( line );
                setThreads( nT );
                usleep( 100000 );
                replyStr = boost::str( boost::format("OK threads %d") % nThreads );
            } else if( what == "image_scale" ) {
                float is = pop<float>( line );
                seeing.set_image_scale( is );
                replyStr = boost::str( boost::format("OK set image_scale %f") % is );
            } else if( what == "min_lock" ) {
                float ml = pop<float>( line );
                boost::trim( line );
                seeing.set_min_lock( ml, line );
                replyStr = boost::str( boost::format("OK set min_lock %f \"%s\"") % ml % line);
            } else if( what == "diam" ) {
                float d = pop<float>( line );
                seeing.set_diam( d );
                replyStr = boost::str( boost::format("OK set diam %f") % d );
            }  else if( what == "lambda" ) {
                float l = pop<float>( line );
                seeing.set_lambda( l );
                replyStr = boost::str( boost::format("OK set lambda %f") % l );
            } else if( what == "ravg" ) {
                float r = pop<float>( line );
                int id = pop<int>( line );
                seeing.set_ravg( r, id );
                if( id ) {
                    replyStr = boost::str( boost::format("OK set ravg %f at %d") % r % id );
                } else {
                    replyStr = boost::str( boost::format("OK set ravg %f") % r );
                }
            } else {  // unrecognized
                if( !what.empty() ) replyStr = "set Huh?";
            }
        }
                
        if( what == "exposure" || what == "gain" ) {    // settings might have changed
            updateCalibID();
        }

    } else if( command == "get" ) {
        string what = popword( line );
        tag = what;
        
        if( auto_subscribe.count(tag) ) subscribe( conn, tag );
        if( !(cam && cam->parse_cmd( command, what, line, replyStr )) ) {
            line = orig_line;       // cam->parse_cmd() is destructive
            what = popword( line );
            if( what == "shm" ) {
                replyStr = "OK shm meta_filename 0 ringbuffer_filename 0";  // FIXME: shared memory not implemented yet for this daemon
            } else if( what == "width" ) {
                replyStr = boost::str( boost::format("OK width %d") % fqueue.width );
            } else if( what == "height" ) {
                replyStr = boost::str( boost::format("OK height %d") % fqueue.height );
            } else if( what == "depth" ) {
                replyStr = boost::str( boost::format("OK depth %d") % fqueue.depth );
            } else if( what == "histo_skip" ) {
                replyStr = boost::str( boost::format("OK histo_skip %lu") % histo_skip );
            } else if( what == "filename" ) {
                replyStr = "OK filename " + filename_prefix;
            } else if( what == "mode" ) {
                replyStr = "OK mode master";
            } else if( what == "filename" ) {
                replyStr = "OK filename filename";
            } else if( what == "fits_camera" ) {
                replyStr = "OK fits_camera fits_camera";
            } else if( what == "state" ) {
                replyStr = "OK state ready";
            } else if( what == "calib" ) {
                replyStr = list_calib();
            } else if( what == "calibs" ) {
                replyStr = list_calibs();
            } else if(what == "dark") {
                bool do_shm WFWFS_UNUSED = (popword(line) == "shm");
                size_t blockSize = dd.nElements()*sizeof(float);
                conn->writeline( boost::str( boost::format("OK dark %d") % blockSize ) );
                conn->syncWrite( dd.getData().get(), blockSize );
            } else if(what == "flat") {
                bool do_shm WFWFS_UNUSED = (popword(line) == "shm");
                size_t blockSize = ff.nElements()*sizeof(float);
                conn->writeline( boost::str( boost::format("OK flat %d") % blockSize ) );
                conn->syncWrite( ff.getData().get(), blockSize );
            } else if(what == "cells") {
                size_t id = pop<int>( line );
                replyStr = "OK cells " + seeing.get_cells( id );
            } else if(what == "shifts") {
                replyStr = "OK shifts " + seeing.get_shifts();
            } else if(what == "ashifts") {
                size_t id = pop<int>( line );
                replyStr = "OK ashifts " + seeing.get_ashifts( id );
            } else if(what == "vars") {
                int duration = pop<int>( line );
                replyStr = "OK vars " + seeing.get_vars( duration );
            } else if(what == "meta") {
                replyStr = boost::str( boost::format("OK meta %d") % meta_cards.size() );
                int cnt(0);
                for( auto& c: meta_cards ) replyStr += "\n[" + to_string(cnt++) + "]: \""  + c + "\"";
            } else if( what == "threads" ) {
                replyStr = boost::str( boost::format("OK threads %d") % nThreads );
            } else {  // unrecognized
                if( !what.empty() ) replyStr = "get Huh?";
            }
        }
    } else if (command == "grab" || command == "stream") {
        int x1 = pop<int>( line );
        int y1 = pop<int>( line );
        int x2 = pop<int>( line );
        int y2 = pop<int>( line );
        int scale = pop<int>( line );
        bool df = false;
        bool do_histo = false;
        size_t fsel = 0;
        string option;
        while( !(option = popword(line)).empty() ) {
            if( option == "darkflat" ) {
                df = true;
            } else if( option == "select" ) {
                fsel = pop<int>( line );
            } else if( option == "histogram" ) {
                do_histo = true;
            }
        }

        get_frame( conn, x1, y1, x2, y2, scale, df, fsel, do_histo );

    } else if ( command == "oburst" ) {
        subscribe( conn, "state" );     // TBD: do we need to send back state ?
        int bcount = pop<int>( line );
        int scannum = line.empty() ? -1 : pop<int>( line );
        long framenum = line.empty() ? -1 : pop<int>( line );
        string statestr = popword(line);
        bool compress = false;
        string filename = make_filename( scannum, framenum, statestr );
        string acc_name = "";
        if( bcount < 0 ) {  // do accumulation
            bcount = std::abs(bcount);
            std::swap( filename, acc_name );
        }
        if( gui_info.use ) {
            compress = gui_info.compress;
        }
        ioService.post( bind(&Daemon::save_fits, this, filename, bcount, bcount, compress, 0, acc_name) );
    } else if ( command == "save" ) {
        int n_frames = pop<int>( line );
        int frames_per_file = line.empty() ? 20 : pop<int>( line );
        string filename = "%DATE%/data/%TIME%/wfwfs_%FRAMENUMBER%.fits";
        string acc_name = "";
        if( n_frames < 0 ) {  // do accumulation
            n_frames = std::abs(n_frames);
            std::swap( filename, acc_name );
        }
        save_fits(  make_filename(filename,-1,-1), n_frames, frames_per_file, true, 0,  make_filename(acc_name,-1,-1) );
        replyStr = "OK save";
    } else if(command == "dark") {
        //darkburst( pop<int>( line ) );
        darkburst( pop<int>( line ) );
        replyStr = "OK dark";
    } else if(command == "flat") {
        //flatburst( pop<int>( line ) );
        flatburst( pop<int>( line ) );
        replyStr = "OK flat";
    } else if( command == "light" ) {
        string state = popword( line );
        if( state == "on" ) {
            light( true );
        } else if( state == "off" ) {
            light( false );
        } else {
            light( !has_light );
        }
        replyStr = "OK light ";
        if( has_light ) replyStr += "on";
        else replyStr += "off";
    } else if(command == "load_calibs") {
        //flatburst( pop<int>( line ) );
        ioService.post( bind( &load_calibs, outputdir+"/calib/" ) );
    } else if( command == "load" ) {
        string what = popword( line );
        if( what == "calibs" ) {
            ioService.post( bind( &load_calibs, outputdir+"/calib/" ) );
        } else if( what == "dark" ) {
            int id = pop<int>( line );
            load_dark( id );
        } else if( what == "flat" ) {
            int id = pop<int>( line );
            load_flat( id );
        } 
    } else if( command == "meta" ) {
        string key = popword( line );
        add_meta( key, line );
    } else if( command == "timing" ) {
        show_timing = !show_timing;
    } else if( command == "adjust") {
        size_t id = pop<int>( line );
        replyStr = "OK adjust " + seeing.adjust_cells( id );
        seeing.draw_cells( cell_mask );
    } else if( command == "shift") {
        PointI s;
        s = pop<int>( line );
        boost::trim( line );
        if( !line.empty() ) {
            s.x = pop<int>( line );
            boost::trim( line ); 
        }
        size_t id = 0;
        if( !line.empty() ) {
            id = pop<int>( line );
        }
        replyStr = "OK shift " + seeing.shift_cells( s, id );
        seeing.draw_cells( cell_mask );
    } else if( command == "sub") {
        string tag = popword( line );
        if( !tag.empty() ) subscribe( conn, tag );
        replyStr = "OK sub " + tag;
    } else if( command == "unsub") {
        string tag = popword( line );
        unsubscribe( conn, tag );
        replyStr = "OK unsub " + tag;
    } else if( command == "broadcast") {
        string tag = popword( line );
        string msg = popword( line );
        broadcast( tag, msg );
        replyStr = "OK broadcast " + tag + " " + msg;
    } else {  // unrecognized
        if( !command.empty() ) replyStr = "Huh?";
    }

    if( !replyStr.empty() ) {
        conn->writeline( replyStr );
        if( !tag.empty() ) broadcast( tag, replyStr, conn );
    }

    if( disconnect ) {
        return false;
    }

    return true;
    
}


void Daemon::onMessage( TcpConnection::Ptr conn ) {

    try {
        //queue<string> lines;
        if( !conn ) {
            throw runtime_error("Broken connection.");
        }
        auto test WFWFS_UNUSED = conn->socket().remote_endpoint();  // check if endpoint exists, will throw if not connected.
        
        string line;
        size_t length;
        bool do_idle(false);
        while( conn && (length = conn->readline( line )) ) {
            do_idle = true;
            if( !conn->socket().is_open() || !processCmd( conn, line ) ) {
                do_idle = false;
            }
        }
        if( do_idle ) {
            conn->idle();
            return;
        }
    } catch( const exception& e ) {      // disconnected or unhandled error -> close socket and return.
        cout  << "onMessage: " << __LINE__ << "  error: " << e.what() << endl;
    } catch( ... ) {      // disconnected or unhandled error -> close socket and return.
        cout  << "onMessage: " << __LINE__ << "  unknown error!!!" << endl;
    }
    
    removeConnection(conn);

    
}


void Daemon::addConnection( const Host& remote_info, TcpConnection::Ptr& conn ) {
    
    lock_guard<mutex> lock( connMutex );
    connections[conn] = remote_info;
// cout << "nConn = " << connections.size() << "  conn: " << (void*)conn.get() << endl ;
}


void Daemon::removeConnection( TcpConnection::Ptr conn ) {

    lock_guard<mutex> lock( connMutex );
    auto connit = connections.find(conn);
    if( connit != connections.end() ) {
        connections.erase( connit );
    }
    unsubscribe( conn );
    conn->setErrorCallback(nullptr);
    conn->setCallback(nullptr);
    conn->socket().close();

}


void Daemon::subscribe( const TcpConnection::Ptr& conn, std::string tag ) {
    
    lock_guard<mutex> lock( globalMutex );
    auto& subs = subscriptions[conn];
    auto ret = subs.insert( tag );
    if( ret.second ) {
//        cout  << "subscribe: " << __LINE__ << "  conn: " << (void*)conn.get() << "  tag: \"" << tag << "\"   current:" <<printArray( subs, " subs" ) << endl;
    }
}


void Daemon::unsubscribe( const TcpConnection::Ptr& conn, std::string tag ) {
    
    lock_guard<mutex> lock( globalMutex );
    auto subsit = subscriptions.find(conn);
    if( subsit != subscriptions.end() ) {
        if( tag.empty() ) {
            subscriptions.erase( subsit );
        } else {
            subsit->second.erase(tag);
        }
    }
    auto& subs = subscriptions[conn];
    subs.erase( tag );
    
}


void Daemon::broadcast( std::string tag, std::string message, TcpConnection::Ptr skip ) {
    
    lock_guard<mutex> lock( globalMutex );
    for( auto it: subscriptions ) {
        if( it.first.get() == skip.get() ) continue;
        if( it.second.count( tag ) ) {
            it.first->writeline( message );
        }
    }
}


void Daemon::start_cam( void ) {
    
    if( !cam ) return;
    
    static const Camera::frame_cb_t cb = bind( &Daemon::queue_frame, this, placeholders::_1, placeholders::_2 );
    ioService.post( bind( &Camera::run, cam.get(), cb ) );

}


void Daemon::stop_cam( void ) {
    
    if( !cam ) return;
    
    cam->stop();
    
}


void Daemon::play( void ) {
    
    if( !cam ) return;
    
    // start logs & bursts

}


void Daemon::pause( void ) {
    
    // stop logs & bursts
    
}


void Daemon::light( bool state ) {
 
    if( has_light != state ) {
        if( state ) {
            seeing.zero_avgs();
            seeing.start_logs();
        } else {
            seeing.stop_logs();
        }
        has_light = state;
    }
}


string Daemon::list_calib( void ) {
    
    CalibSets& cs = getCalibSet( current_calib );
    string ret;
    if( cs.dd.size() ) ret += "Darks:\n";
    size_t i(0);
    for( auto& d: cs.dd ) {
        ret += "[" + to_string(i) + "]: " + d.filename.string() + "\n";
        i++;
    }
    
    if( cs.ff.size() ) ret += "Flats:\n";
    i = 0;
    for( auto& f: cs.ff ) {
        ret += "[" + to_string(i) + "]: " + f.filename.string() + "\n";
        i++;
    }

    if( ret.back() == '\n' ) ret.resize( ret.size()-1 );

    return ret;
    
}


string Daemon::list_calibs( void ) {

    string ret;
    for( auto& a: calib_archive ) {
        CalibSets& cs = a.second;
        string idString = to_string(a.first.width) + "," + to_string(a.first.height) + ","
        + to_string(a.first.depth) + "," + to_string(a.first.exposure) + "," + to_string(a.first.gain);
        if( cs.dd.size() ) ret += "Darks: (" +idString+ ")\n";
        size_t i(0);
        for( auto& d: cs.dd ) {
            ret += "[" + to_string(i) + "]: " + d.filename.string() + "\n";
            i++;
        }
        
        if( cs.ff.size() ) ret += "Flats: (" +idString+ ")\n";
        i = 0;
        for( auto& f: cs.ff ) {
            ret += "[" + to_string(i) + "]: " + f.filename.string() + "\n";
            i++;
        }
    }

    if( ret.back() == '\n' ) ret.resize( ret.size()-1 );

    return ret;
    
}


void Daemon::load_dark( int id ) {

    string calibPath = outputdir + "/calib/";
    CalibSets& cs = getCalibSet( current_calib );
    if( cs.dd.empty() ) return;
    
    vector<CalibFile> tmpDD;
    std::copy( cs.dd.begin(), cs.dd.end(), std::back_inserter(tmpDD) );
    if( id < 0 ) id = tmpDD.size()-1;       // id < 0 will load latest file
    if( (size_t)id < tmpDD.size() ) {
        bfs::path link_name(calibPath+"dark.fits");
        if( bfs::is_symlink(link_name) ) {
            bfs::remove(link_name);
        }
        
        bfs::create_symlink( bfs::relative(tmpDD[id].filename, calibPath), link_name );
        if( !bfs::exists(link_name) ) {
            readFile( tmpDD[id].filename.string(), dd );
            ddPtr = dd.get();
        }
        update_gain();
    }
    
//    broadcast( "dark", "OK dark 0");
//    broadcast( "flat", "OK flat 0");

}


void Daemon::load_flat( int id ) {

    string calibPath = outputdir + "/calib/";
    CalibSets& cs = getCalibSet( current_calib );
    if( cs.ff.empty() ) return;
    
    vector<CalibFile> tmpFF;
    std::copy( cs.ff.begin(), cs.ff.end(), std::back_inserter(tmpFF) );
    if( id < 0 ) id = tmpFF.size()-1;       // id < 0 will load latest file
    if( (size_t)id < tmpFF.size() ) {
        bfs::path link_name(calibPath+"flat.fits");
        if( bfs::is_symlink(link_name) ) {
            bfs::remove(link_name);
        }
        if( !bfs::exists(link_name) ) {
            bfs::create_symlink( bfs::relative(tmpFF[id].filename, calibPath), link_name );
        }
        readFile( tmpFF[id].filename.string(), ff );
        update_gain();
    }
    
//    broadcast( "flat", "OK flat 0");
    
}


void Daemon::updateCalibID( void ) {
    
    current_calib.depth = (1<<fqueue.depth)-1;
    current_calib.width = fqueue.width;
    current_calib.height = fqueue.height;
    current_calib.exposure = cam->get_exposure();
    current_calib.gain = cam->get_gain();
    
}


void Daemon::makeMeta( vector<string>& cards ) {
    
    cards.clear();
    cards.insert( cards.end(), meta_cards.begin(), meta_cards.end() );      // copy meta cards.

    // update/add some (maybe) modified values
    if( !gui_info.observer.empty() ) Fits::updateCard( cards, Fits::makeCard( "OBSERVER", gui_info.observer ) );
    if( !gui_info.target.empty() ) Fits::updateCard( cards, Fits::makeCard( "OBJECT", gui_info.target ) );
    Fits::updateCard( cards, Fits::makeCard( "DATAMIN", 0 ) );
    Fits::updateCard( cards, Fits::makeCard( "DATAMAX", (1<<fqueue.depth)-1 ) );
    Fits::updateCard( cards, Fits::makeCard( "XPOSURE", cam->get_exposure(), "[s] Exposure time" ) );
    Fits::updateCard( cards, Fits::makeCard( "DETGAIN", cam->get_gain(), "[dB] or camera specific unit" ) );
    Fits::updateCard( cards, Fits::makeCard( "CADENCE", cam->get_interval(), "[s]" ) );
    Fits::updateCard( cards, Fits::makeCard( "DETOFFS", 0, "[counts] or camera specific unit" ) );
    float temp = cam->get_temp();
    if( isfinite(temp) ) {
        Fits::updateCard( cards, Fits::makeCard( "DETTEMP", temp,"[C] Current temperature of the detector" ) );
    }

    string comments = gui_info.comments;
    while( !comments.empty() ) {
        string comment = popword( comments, ";" );
        if( !comment.empty() ) {
            Fits::insertCard( cards, Fits::makeCard( "COMMENT", comment ) );
        }
    }

//     Fits::updateCard( cards, Fits::makeCard( "FIELD", "name" ) );

}


int Daemon::accumulate( Array<uint32_t>& acc, size_t n ) {
    
    acc.zero();
    int cnt(0);
    
    if( n ) {
        broadcast( "state", "OK state waiting" );
        size_t id(0);
        while( n-- ) {
            Frame& f = fqueue.getFrame( id, true );
            if( id && (f.id != id) ) {
                cout << "Requested id = " << id << " but got id = " << f.id << endl;
            }
            id = f.id+1;
            fqueue.addFrame( f, acc.get() );
            cnt++;
        }
        broadcast( "state", "OK state ready" );
    }
    
    return cnt;    
    
}


void Daemon::do_accumulation( Array<float>& out, Array<uint32_t>& acc, string fn, size_t n ) {
    
    bfs::path fnPath(fn);
    bfs::path dirPath = fnPath.parent_path( );
    maybeCreateDir( dirPath );
    
    bpx::ptime beg_time = bpx::microsec_clock::universal_time();
    
    promise<int> prom;
    future<int> fut = prom.get_future();
    ioService.post( [&](){
        prom.set_value( accumulate( acc, n ) );
    } );

    
    string date = to_iso_extended_string( beg_time.date() );
    string timestamp = to_iso_extended_string( beg_time );
    size_t pos = timestamp.find_last_of('.');
    string timestamp_s = timestamp;
    if( pos != string::npos ) timestamp_s = timestamp.substr(0,pos);
    
    shared_ptr<Fits> hdr;
    if( !hdr ) {
        hdr.reset( new Fits() );
    }
    vector<string>& cards = hdr->primaryHDU.cards;
    makeMeta(cards);
    
    int datamax = (1<<fqueue.depth)-1;
    float exp_time = cam->get_exposure();
    float det_gain = cam->get_gain();
    int h = fqueue.height;
    int w = fqueue.width;
    
    
    Fits::updateCard( cards, Fits::makeCard( "DATE-OBS", timestamp_s ) );
    Fits::updateCard( cards, Fits::makeCard( "DATE", timestamp ) );
    Fits::updateCard( cards, Fits::makeCard( "DATE-BEG", timestamp, "Start time of summed observations" ) );
    
    int cnt = fut.get();     // sync-point, efter this the accumulation is completed.
    
    bpx::ptime end_time = bpx::microsec_clock::universal_time();
    timestamp = to_iso_extended_string(end_time);
    bpx::time_duration elapsed = (end_time - beg_time);
    end_time = beg_time + elapsed/2;
    Fits::updateCard( cards, Fits::makeCard( "DATE-AVG", to_iso_extended_string(end_time), "Average time of summed observations" ) );
    Fits::updateCard( cards, Fits::makeCard( "DATE-END", timestamp, "End time of summed observations" ) );
    
    out = acc;
    out /= cnt;

    Fits::updateCard( cards, Fits::makeCard( "FILENAME", fnPath.filename().string() ) );
    Fits::updateCard( cards, Fits::makeCard( "XPOSURE", cnt*exp_time, "[s] Total exposure time" ) );
    Fits::updateCard( cards, Fits::makeCard( "TEXPOSUR", exp_time, "[s] Single exposure time" ) );
    Fits::updateCard( cards, Fits::makeCard( "NSUMEXP", cnt, "Number of summed exposures" ) );
    
    Fits::write( fn, out, hdr );
            
    CalibFile cf( fn, beg_time );
    fn = fnPath.filename().string();
    if( fn.find("dark_") == 0 || fn.find("flat_") == 0 ) {
        CalibSets& cs = getCalibSet( w, h, datamax, exp_time, det_gain );
        if( fn.find("dark_") == 0 ) {
            auto ret = cs.dd.emplace( cf );
            if( !ret.second ) {
                // do we want/need to modify anything if there is a collision ?
            }
        } else if ( fn.find("flat_") == 0) {
            auto ret = cs.ff.emplace( cf );
            if( !ret.second ) { 
                // do we want/need to modify anything if there is a collision ?
            }
        }
        
        update_gain();
    }
}


void Daemon::darkburst( size_t n ) {
    
    string filename = "calib/%DATE%/darks/%TIME%/dark_%FRAMENUMBER%.fits";
    string acc_name = "calib/%DATE%/dark_%TIME%.fits";

    save_fits( make_filename(filename,-1,-1), n, 20, true, 0, make_filename(acc_name,-1,-1) );
    
    load_calibs( outputdir+"/calib/" );
    load_dark(-1);
    
}


void Daemon::flatburst( size_t n ) {
    
    string filename = "calib/%DATE%/flats/%TIME%/flat_%FRAMENUMBER%.fits";
    string acc_name = "calib/%DATE%/flat_%TIME%.fits";
    
    save_fits( make_filename(filename,-1,-1), n, 20, true, 0, make_filename(acc_name,-1,-1) );
    
    load_calibs( outputdir+"/calib/" );
    load_flat(-1);
    
}


void Daemon::update_gain( void ) {
    
    if( !ff.sameSizes(dd) ) {
        throw std::logic_error("Daemon::update_gain: FF & DD does not have the same dimensions." );
    }
    
    if( !ff.sameSizes(gg) ) {
        gg.resize( ff.dimensions(true) );
        ggPtr = gg.get();
    }

    size_t nElements = ff.nElements();
    std::transform( ff.get(), ff.get()+nElements, ddPtr, ggPtr, std::minus<float>() );
    
    ArrayStats stats;
    stats.getMinMaxMean( gg );
    transform( ggPtr, ggPtr+nElements, ggPtr, [&](float a) {
            if( a > 0.0 ) {
                double tmp = stats.max/a;
                if( tmp > 4 ) return 4.0;
                if( tmp < 0.1 ) return 0.1;
                return tmp;
            }
            return 0.0;
    } );
    
    stats.getStats(gg);
    gg /= stats.mean*0.99;   // NOTE: the 0.99 is just to approximately compensate for the intensity decrease from dd subtraction

    try {
        bfs::path tmpFN( outputdir );
        tmpFN /= "calib/";
        tmpFN /= "gaintable.fits";
        if( bfs::exists(tmpFN) && !bfs::remove(tmpFN) ) {
            cerr << boost::format( "Failed to remove existing file: %s" ) % tmpFN << endl;
            return;
        }

        Fits::write( tmpFN.string(), gg );
    } catch( const std::exception& ) {
        // TODO
    } catch ( ... ) {
        // ignore unrecognized exceptions.
    }
    
}


void Daemon::die(void) {

    stop();
    
}


void Daemon::die( TcpConnection::Ptr& conn, bool urgent ) {
    
    if( urgent ) {
        die();
        return;
    }

}


void Daemon::softExit( void ) {
/*    
    if( myInfo.info.peerType == Host::TP_MASTER ) {
        cerr << "Daemon::softExit, not implemented for master yet." << endl;
    } else {
        cout << "Slave will exit after the current job is completed." << endl;
        worker.exitWhenDone();
    }*/
}


void Daemon::reset( TcpConnection::Ptr& conn, bool urgent ) {
    
    if( urgent ) {  // received on a slave, just do it without replying
        reset();
        return;
    }

}


void Daemon::thumbnail( TcpConnection::Ptr conn ) {

    const size_t sz = 32*32;
    uint8_t thumb_buf[ sz ];
    uint8_t *out = thumb_buf;
    memset( thumb_buf, 0, sz );

    const size_t width = fqueue.width;
    const size_t height = fqueue.height;
    const size_t depth = fqueue.depth;
    
    size_t step, xoff, yoff;

    if( width < height ) {
        step = width / 32;
    } else {
        step = height / 32;
    }
    
    xoff = (width - step * 31) / 2;
    yoff = (height - step * 31) / 2;

    Frame& f = fqueue.getFrame( 0 );    // get latest frame
    if( depth > 8 ) {
        int shift = (depth - 8);
        for(int y=0; y<32; ++y) {
            uint16_t* f_row = reinterpret_cast<uint16_t*>(f.data) + (y*step+yoff)*width + xoff;
            for(int x=0; x<32; ++x ) {
                out[x] = f_row[x * step] >> shift;
            }
            out += 32;
        }
    } else {
        for( int y=0; y<32; ++y ) {
            uint8_t* f_row = f.data + (y*step+yoff)*width + xoff;
            for(int x=0; x < 32; ++x ) {
                out[x] = f_row[x * step];
            }
            out += 32;
        }
    }

    conn->writeline( "OK thumbnail" );
    conn->syncWrite( thumb_buf, sz );
    
}


int Daemon::test_threads( int nT, int t ) {
    
    int ret(0);
    if( nT>0 && t>0 ) {
        
        atomic<int> run_test(1);
        run_test = true;
        atomic<int> count(0);
        atomic<int> maxConcurrent(0);
        std::future<int> retFut = std::async( [&](){
            sleep(t);
            run_test.store(0);
            return 1;
        });
        
        for( int i(0); i<nT; ++i ) {
            ioService.post([&](){
                count++;
                maxConcurrent = std::max<int>(maxConcurrent,count);
                while( run_test.load() ) ;
                count--;
            });
        }
        ret = retFut.get();    //
        usleep(10000);
        
        cout << "test_threads:  max concurrent threads available: " << maxConcurrent << endl;
        
    }
    
    return ret;
    
}
