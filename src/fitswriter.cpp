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
#include "fitswriter.hpp"

#include "ricecompress.hpp"
#include "version.hpp"

#include <algorithm>
#include <atomic>
#include <functional>
#include <iostream>
#include <string.h>
#include <sys/time.h>
#include <unistd.h>
#include <endian.h>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using namespace std;
using namespace wfwfs;
namespace bfs = boost::filesystem;


// Initialize static members
int FitsWriter::activeCount(0);
int FitsWriter::totalCount(0);
list<shared_ptr<uint8_t>> FitsWriter::buffers;
vector<string> FitsWriter::globalMeta;
mutex FitsWriter::globalMtx;
mutex FitsWriter::bufMtx;
set<pair<boost::posix_time::ptime,boost::posix_time::ptime>> FitsWriter::saves;

namespace {

    mutex writeMtx;

    // Pad file to a specified boundary
    void pad_to( int fd, size_t boundary ) {
        size_t pos = lseek( fd, 0, SEEK_CUR) % boundary;
        if( pos ) {
            size_t pad = boundary - pos;
            unique_ptr<char[]> tmp( new char[pad] );
            memset( tmp.get(), 0, pad );
            size_t n = write( fd, tmp.get(), pad );
            if( n != pad ) {
                cerr << "FitsWriter: Failed to write padding: " << strerror(errno) << endl;
            }
        }
    }

    void write_meta( int fd, const vector<string>& meta ) {
        
        if( meta.empty() ) return;
        
        size_t metaSize = ((meta.size()-1)/36 + 1)*2880;
        std::unique_ptr<char[]> hdu( new char[ metaSize] );
        char* hPtr = hdu.get();
        memset( hPtr, ' ', metaSize );
        for( auto& c: meta ) {
            std::copy( c.begin(), c.end(), hPtr );
            hPtr += 80;
        }
        size_t n = write( fd, hdu.get(), metaSize );
        if( n != metaSize ) {	// write offset+size for rows in binary table
            cerr << "\nwrite_meta: Failed to write hdu: " << strerror(errno) << endl;
        }
    }

}


FitsWriter::FitsWriter( FrameQueue& FQ, const string& afn, int nT, bool compress )
    : running(true), do_fsync(false), do_compress(compress), index(0), nthreads(nT), npixels(FQ.width*FQ.height), fd(-1),
      do_acc(false), acc_filename(afn), first_ts(bpx::not_a_date_time), last_ts(bpx::not_a_date_time), fq(FQ) {

    bytes_per_pixel = (fq.depth-1)/8 + 1;
    frame_count = 0;

    if( !acc_filename.empty() ) {
        bfs::path fnPath( acc_filename );
        bfs::path dirPath = fnPath.parent_path();
        maybeCreateDir( dirPath );
        if( !dirPath.empty() && bfs::exists(dirPath) ) {
            do_acc = true;
            acc.resize( fq.height, fq.width );
            acc.zero();
        } else {
            // TODO throw
        }
    }

    lock_guard<mutex> lock( globalMtx );
    ++totalCount;
}


FitsWriter::~FitsWriter() {
    
    if( do_acc ) {
        write_acc();
    }
    
    if( do_write || do_acc ) {
        if( !first_ts.is_not_a_date_time() && !last_ts.is_not_a_date_time() ) {
            add_save( first_ts, last_ts );
        }
    }
    
    lock_guard<mutex> lock( globalMtx );
    --totalCount;
    
}


void FitsWriter::makeHdr( void ) {
    
    if( !hdr ) {
        hdr.reset( new Fits() );
    }
    
    vector<string>& cards = hdr->primaryHDU.cards;
    cards.clear();
    Fits::insertCard( cards, Fits::makeCard( "SIMPLE", true ) );
    Fits::insertCard( cards, Fits::makeCard( "BITPIX", (fq.depth > 8?16:8)) );
    if( do_compress ) {
        Fits::insertCard( cards, Fits::makeCard( "NAXIS", 0 ) );
    } else {
        Fits::insertCard( cards, Fits::makeCard( "NAXIS", (times.size()>1)?3:2 ) );
        Fits::insertCard( cards, Fits::makeCard( "NAXIS1", fq.width ) );
        Fits::insertCard( cards, Fits::makeCard( "NAXIS2", fq.height ) );
        if( times.size() > 1 ) {
            Fits::insertCard( cards, Fits::makeCard( "NAXIS3", times.size() ) );
        }
    }
    Fits::insertCard( cards, Fits::makeCard<string>( "EXTNAME", "Main" ) );
    Fits::insertCard( cards, Fits::makeCard<string>( "TAB_HDUS", "TABULATIONS;DATE-BEG" ) );
    
    string timestamp = to_iso_extended_string( first_ts );
    size_t pos = timestamp.find_last_of('.');
    string timestamp_s = timestamp;
    if( pos != string::npos ) timestamp_s = timestamp.substr(0,pos);
    Fits::insertCard( cards, Fits::makeCard<string>( "DATE", timestamp ) );
    Fits::insertCard( cards, Fits::makeCard<string>( "DATE-OBS", timestamp_s ) );
    bpx::time_duration elapsed = (last_ts - first_ts);
    first_ts += elapsed/2;
    Fits::insertCard( cards, Fits::makeCard( "DATE-AVG", to_iso_extended_string(first_ts), "Average time of observations" ) );
    Fits::insertCard( cards, Fits::makeCard( "DATE-END", to_iso_extended_string(last_ts), "End time of observations" ) );

    Fits::insertCard( cards, Fits::makeCard<string>( "WFWFSVER", getVersionString() ) );
    
    cards.insert( cards.end(), globalMeta.begin(), globalMeta.end() );      // copy global meta-data.
    cards.insert( cards.end(), extra_meta.begin(), extra_meta.end() );      // copy extra meta-data.

    Fits::removeCards( cards, "END" );                   // just in case it is not the last card, or if there are multiple.
    Fits::insertCard( cards, Fits::makeCard( "END" ) );

}


void FitsWriter::save( const string& filename, size_t& next_frame, int nF ) {

    nframes = nF;
    do_write = !filename.empty();
    do_compress = (do_compress && do_write && (bytes_per_pixel == 2) );
    nthreads = do_compress?nthreads:1;

    offsets.assign( 2*nframes, 0 );
    frames.clear();
    times.clear();

    if( do_write ) {        
        open_file( filename, nframes*fq.frameSize+50*2880 );
        bfs::path fnPath(filename);
        Fits::updateCard( extra_meta, Fits::makeCard<string>( "FILENAME", fnPath.filename().string() ) );

        makeHdr();
        hdrEnd = ((hdr->primaryHDU.cards.size()-1)/36 + 1)*2880;		// end of primary header, and possibly start of compressed header
        dataStart = hdrEnd;
        if ( do_compress ) {
            index = 0;
            dataStart += 2880+nframes*2*sizeof(int);
        }
        lseek( fd, dataStart, SEEK_SET );	// set file-pointer to where the frames should be saved.
    }
    
    running = true;
    for( int i=0; i<nthreads; ++i ) {
        threads.push_back( thread(bind( &FitsWriter::thread_run, this )) );
    }

    int cnt(0);
    while( cnt < nframes ) {
        LockedFrame lf( fq.getFrame( next_frame, true ) );
        if( next_frame && (lf.frame.id != next_frame) ) {
            cout << "Requested id = " << next_frame << " but got id = " << lf.frame.id << endl;
        }
        next_frame = lf.frame.id+1;
        push( lf.frame );
        cnt++;
    }
    
    wait();

    if( do_write ) {        
        close_file();
    }
    
}


void FitsWriter::wait( void ) {

    {
        unique_lock<mutex> lock(queueMtx);
        while( !frames.empty() ) cond.wait_for(lock,std::chrono::milliseconds(1));
        running = false;
    }
    cond.notify_all();
    for( auto& t: threads ) t.join();

    threads.clear();
    
}


void FitsWriter::thread_run(void) {

    shared_ptr<uint8_t> f;
    unique_ptr<uint8_t[]> cData;
    if( do_compress ) {
        cData.reset( new uint8_t[ fq.frameSize ] );		// compressed storage of the same size as a raw frame
    }

    while( running ) {
        if( (f = pop()) ) {
            unique_lock<mutex> lock(queueMtx);
            int ind = index;
            index += 2;
            if( do_acc ) {
                fq.addFrame( f.get(), acc.get() );
                frame_count++;
            }
            lock.unlock();
            if( do_write ) {
                if( do_compress ) {
                    int frameOffset(0);
                    int16_t* tmpPtr = reinterpret_cast<int16_t*>(f.get());
                    // TODO: set blocksize to width (i.e. process one row at a time) when/if the v. 4.0 FITS standard allows it.
                    int cSize = rice_comp16( tmpPtr, npixels, cData.get(), fq.frameSize, 32 );
                    if( cSize > 0 ) {
                        lock_guard<mutex> wlock(writeMtx);
                        frameOffset = lseek( fd, 0, SEEK_CUR );
                        int count = write( fd, cData.get(), cSize );
                        if( count != cSize ) {
                            fprintf( stderr, "FitsWriter: write failed for frame #%d, count=%d  cSize=%d: %s\n",
                                     (ind/2), count, cSize, strerror(errno) );
                        }
                    } else {
                        fprintf( stderr, "FitsWriter: compression failed for frame #%d.\n", (ind/2) );
                        cSize = 0;
                    }
                    lock.lock();
                    offsets[ind] = cSize;
                    offsets[ind+1] = frameOffset;
                } else {
                    lock.unlock();
                    if( bytes_per_pixel == 2 ) {
                        int16_t* ptr = reinterpret_cast<int16_t*>( f.get() );
                        for( int i=0; i<npixels; ++i ) ptr[i] = htobe16(ptr[i]);
                    } else if( bytes_per_pixel == 4 ) {
                        int32_t* ptr = reinterpret_cast<int32_t*>( f.get() );
                        for( int i=0; i<npixels; ++i ) ptr[i] = htobe32(ptr[i]);
                    }
                    {
                        lock_guard<mutex> wlock(writeMtx);
                        size_t count = write( fd, f.get(), fq.frameSize );
                        if( count != fq.frameSize ) {
                            fprintf( stderr, "FitsWriter: write failed for frame #%d, count=%zu  frameSize=%zu: %s\n",
                                     (ind/2), count, fq.frameSize, strerror(errno) );
                        }
                    }
                }
            }
            return_buf( f );
        } else {
            unique_lock<mutex> lock(queueMtx);
            cond.wait(lock);
        }
    }

}


void FitsWriter::push( const Frame& f ) {
    if( f.data ) {
        shared_ptr<uint8_t> buf = get_buf( fq.frameSize );
        memcpy( buf.get(), f.data, fq.frameSize );
        {
            lock_guard<mutex> lock(queueMtx);
            frames.push_back( buf );
            times.push_back( f.timestamp );
            if( first_ts.is_not_a_date_time() || (first_ts > f.timestamp) ) {
                first_ts = f.timestamp;
            }
            if( last_ts.is_not_a_date_time() || (last_ts < f.timestamp) ) {
                last_ts = f.timestamp;
            }
        }
        cond.notify_all();
    }
}


void FitsWriter::push( void* data, bpx::ptime ts ) {
    if( data ) {
        shared_ptr<uint8_t> buf = get_buf( fq.frameSize );
        memcpy( buf.get(), data, fq.frameSize );
        {
            lock_guard<mutex> lock(queueMtx);
            frames.push_back( buf );
            times.push_back( ts );
        }
        cond.notify_all();
    }
}


shared_ptr<uint8_t> FitsWriter::pop() {
    lock_guard<mutex> lock(queueMtx);
    if( frames.empty() ) return nullptr;
    shared_ptr<uint8_t> ret = frames.front();
    frames.pop_front();
    return ret;
}


void FitsWriter::open_file( const string& filename, size_t sz ) {
    
    bfs::path fnPath(filename);
    bfs::path dirPath = fnPath.parent_path();
    maybeCreateDir( dirPath );

    fd = open( filename.c_str(), O_WRONLY|O_CREAT, 0664 ); // |O_EXCL
    if( fd < 0 ) {
        //connection->write(format("ERROR :Could not open output file '%s': %s", filename.c_str(), strerror(errno)));
        return;
    }
    
    // Allocate enough space for the whole file
    if( sz ) {
        posix_fallocate( fd, 0, sz );
    }
    
}


void FitsWriter::close_file( void ) {

    makeHdr();      // call makeHdr again to store the correct sizes.
    
    // Pad to a multiple of 2880 bytes, because FITS requires that
    pad_to( fd, 2880 );
    write_exptime_table();
    
    // Pad to a multiple of 2880 bytes, because FITS requires that
    pad_to( fd, 2880 );
    
    size_t pos = lseek( fd, 0, SEEK_CUR );		// save EOF position
    lseek( fd, 0, SEEK_SET );			        // move back to primary HDU
    
    write_meta( fd, hdr->primaryHDU.cards );

    if( do_compress ) {
        
        maxRowSize = 0;
        pcount = 0;
        
        int dataStart = INT32_MAX;
        for( size_t i=0; i<offsets.size(); i+=2 ) {
            if( offsets[i] > maxRowSize ) maxRowSize = offsets[i];
            if( offsets[i+1] && (offsets[i+1] < dataStart) ) dataStart = offsets[i+1];
            pcount += offsets[i];
        }
        
        for( size_t i=0; i<offsets.size(); i+=2 ) {
            offsets[i] = htobe32(offsets[i]);
            if( offsets[i+1] ) offsets[i+1] = htobe32(offsets[i+1]-dataStart);
        }
        
        
        const vector<string> comp_meta = {
            Fits::makeCard( "XTENSION", "BINTABLE", "binary table extension" ),
            Fits::makeCard( "BITPIX", 8, "8-bit data" ),
            Fits::makeCard( "NAXIS", 2, "2-dimensional binary table" ),
            Fits::makeCard( "NAXIS1", 8, "width of table in bytes" ),
            Fits::makeCard( "NAXIS2", nframes, "number of rows in table" ),
            Fits::makeCard( "PCOUNT", pcount, "size of special data area" ),
            Fits::makeCard( "GCOUNT", 1, "one data group (required)" ),
            Fits::makeCard( "TFIELDS", 1, "number of fields in each row" ),
            Fits::makeCard( "TTYPE1", "COMPRESSED_DATA", "label for field #1" ),
            Fits::makeCard( "TFORM1", "1PB("+to_string(maxRowSize)+")", "data format of field: variable length array" ),
            // start of mandatory keywords for tiled image compression (10.1.1 in v. 4.0 of the FITS standard)
            Fits::makeCard( "ZIMAGE", true, "extension contains compressed image" ),
            Fits::makeCard( "ZTILE1", fq.width, "size of tiles" ),
            Fits::makeCard( "ZTILE2", fq.height, "size of tiles" ),
            Fits::makeCard( "ZTILE3", 1, "size of tiles" ),
            Fits::makeCard( "ZCMPTYPE", "RICE_1", "compression algorithm used" ),
            Fits::makeCard( "ZNAME1", "BLOCKSIZE" ),
            Fits::makeCard( "ZVAL1", 32 ),
            Fits::makeCard( "ZNAME2", "BYTEPIX" ),
            Fits::makeCard( "ZVAL2", 2 ),
            Fits::makeCard( "ZSIMPLE", true ),
            Fits::makeCard( "ZBITPIX", 16, "bitpix of original image" ),
            Fits::makeCard( "ZNAXIS", 3, "naxis of original image" ),
            Fits::makeCard( "ZNAXIS1", fq.width, "naxis1 of original image" ),
            Fits::makeCard( "ZNAXIS2", fq.height, "naxis2 of original image" ),
            Fits::makeCard( "ZNAXIS3", nframes, "naxis3 of original image" ),
            Fits::makeCard( "END" )
        };

        lseek( fd, hdrEnd, SEEK_SET );		// move to end of primary header
        write_meta( fd, comp_meta );

        size_t sz = offsets.size()*sizeof(int);
        size_t n = write( fd, offsets.data(), sz );
        if( n != sz ) {	// write offset+size for rows in binary table
            cerr << "FitsWriter: Failed to write row information: " << strerror(errno) << endl;
        }
    }

    if( ftruncate( fd, pos ) ) {		// truncate file
        cerr << "FitsWriter: Failed to truncate file at position " << pos << ": " << strerror(errno) << endl;
    }
    
    if( do_fsync ) {	// Do an fsync if requested, and close the file
        fsync(fd);
    }
    
    close(fd);
    fd = -1;
    
}


void FitsWriter::write_exptime_table( void ) {

    vector<string> tmeta;
    tmeta.reserve(16);
    Fits::insertCard( tmeta, Fits::makeCard<string>( "XTENSION", "TABLE   ") );
    Fits::insertCard( tmeta, Fits::makeCard( "BITPIX", 8 ) );
    Fits::insertCard( tmeta, Fits::makeCard( "NAXIS", 2 ) );
    Fits::insertCard( tmeta, Fits::makeCard( "NAXIS1", 26 ) );
    Fits::insertCard( tmeta, Fits::makeCard( "NAXIS2", (int)times.size() ) );
    Fits::insertCard( tmeta, Fits::makeCard( "PCOUNT", 0) );
    Fits::insertCard( tmeta, Fits::makeCard( "GCOUNT", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "TFIELDS", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "TTYPE1", "DATE-BEG") );
    Fits::insertCard( tmeta, Fits::makeCard( "TBCOL1", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "TFORM1", "A26") );
    Fits::insertCard( tmeta, Fits::makeCard( "TUNIT1", "time") );
    Fits::insertCard( tmeta, Fits::makeCard( "EXTNAME", "TABULATIONS") );
    Fits::insertCard( tmeta, Fits::makeCard( "SOLARNET", 0.5) );
    Fits::insertCard( tmeta, Fits::makeCard( "OBS_HDU", 1) );
    Fits::insertCard( tmeta, Fits::makeCard( "END" ) );
    
    write_meta( fd, tmeta );
    pad_to( fd, 2880 );

    for( auto& t: times ) {
        string ts = to_iso_extended_string( t );
        ts.resize( 26, ' ' );
        if( write( fd, ts.data(), 26 ) != 26 ) {
            cerr << "write_exptime_table: Failed to write timestamp: " << strerror(errno) << endl;
        }
    }

}


void FitsWriter::write_acc( void ) {
    
    if( !frame_count ) {
        return;
    }
    
    bfs::path fnPath( acc_filename );

    Array<float> out;
    out = acc;
    if( frame_count ) {
        out *= 1.0/frame_count;
    }
    
    do_compress = false;
    times.clear();
    
    makeHdr();
    
    vector<string>& cards = hdr->primaryHDU.cards;
    Fits::removeCards( cards, "SIMPLE" );
    Fits::removeCards( cards, "BITPIX" );
    Fits::removeCards( cards, "NAXES" );
    Fits::removeCards( cards, "NAXIS1" );
    Fits::removeCards( cards, "NAXIS2" );
    Fits::removeCards( cards, "TAB_HDUS" );
    
    float exp_time = Fits::getValue<float>( cards, "XPOSURE" );
    if( exp_time != 0.0 ) {
        Fits::updateCard( cards, Fits::makeCard( "XPOSURE", frame_count*exp_time, "[s] Total exposure time" ) );
        Fits::updateCard( cards, Fits::makeCard( "TEXPOSUR", exp_time, "[s] Single exposure time" ) );
    }
    Fits::updateCard( cards, Fits::makeCard( "NSUMEXP", frame_count, "Number of summed exposures" ) );
    
    string timestamp = to_iso_extended_string( first_ts );
    size_t pos = timestamp.find_last_of('.');
    string timestamp_s = timestamp;
    if( pos != string::npos ) timestamp_s = timestamp.substr(0,pos);
    Fits::updateCard( cards, Fits::makeCard<string>( "DATE-OBS", timestamp_s ) );
    Fits::updateCard( cards, Fits::makeCard<string>( "DATE", timestamp ) );
    Fits::updateCard( cards, Fits::makeCard( "DATE-BEG", timestamp, "Start time of summed observations" ) );
    bpx::time_duration elapsed = (last_ts - first_ts);
    first_ts += elapsed/2;
    Fits::updateCard( cards, Fits::makeCard( "DATE-AVG", to_iso_extended_string(first_ts), "Average time of summed observations" ) );
    Fits::updateCard( cards, Fits::makeCard( "DATE-END", to_iso_extended_string(last_ts), "End time of summed observations" ) );
    Fits::updateCard( cards, Fits::makeCard<string>( "FILENAME", fnPath.filename().string() ) );

    try {
        if( bfs::exists(fnPath) && !bfs::remove(fnPath) ) {
            cerr << boost::format( "Failed to remove existing file: %s" ) % fnPath << endl;
            return;
        }
        Fits::write( fnPath.string(), out, hdr );
    } catch( const std::exception& ) {
        // TODO
    } catch ( ... ) {
        // ignore unrecognized exceptions.
    }


}


std::shared_ptr<uint8_t> FitsWriter::get_buf( size_t N ) {
    shared_ptr<uint8_t> buf;
    {
        lock_guard<mutex> lock(bufMtx);
        if( !buffers.empty() ) {
            buf = buffers.front();
            buffers.pop_front();
        }
    }
    if(!buf) {
        buf.reset( new uint8_t[N], [](uint8_t*& p) { delete[] p; p=nullptr; });
    }
    return buf;
}


void FitsWriter::return_buf( std::shared_ptr<uint8_t>& buf ) {
    lock_guard<mutex> block(bufMtx);
    buffers.push_back(buf);
    buf.reset();
}


void FitsWriter::clear_bufs( void ) {
    lock_guard<mutex> block(bufMtx);
    buffers.clear();
}


string FitsWriter::get_saves( void ) {
    
    static const bpx::ptime epoch_time( boost::gregorian::date(1970,1,1) ); 
    string ret;
    
    lock_guard<mutex> lock( globalMtx );
    for( auto& s: saves ) {
        bpx::time_duration tmp = s.first - epoch_time;
        size_t nSecs = tmp.total_seconds();
        size_t nMicros = s.first.time_of_day().total_microseconds() - s.first.time_of_day().total_seconds()*1000000L;
        string line = boost::str( boost::format("%ld.%06ld") % nSecs % nMicros );
        tmp = s.second - epoch_time;
        nSecs = tmp.total_seconds();
        nMicros = s.second.time_of_day().total_microseconds() - s.second.time_of_day().total_seconds()*1000000L;
        line += boost::str( boost::format("    %ld.%06ld") % nSecs % nMicros );
        ret += line + "\n";
    }

    return ret;
}


void FitsWriter::add_save( boost::posix_time::ptime from, boost::posix_time::ptime to ) {
    
    lock_guard<mutex> lock( globalMtx );
    saves.emplace( make_pair(from, to) );
    
}
