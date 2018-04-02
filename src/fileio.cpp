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
#include "fileio.hpp"

#include "filefits.hpp"
#include "stringutil.hpp"

#include <iostream>
#include <future>
#include <map>
#include <mutex>

using namespace wfwfs;
using namespace std;

ErrorHandling wfwfs::errorHandling = EH_PRINT;

namespace {

    mutex fileMutex;




    template <typename T>
    map<string, shared_ptr<T>>& getFileCache( void ) {
        static map<string, shared_ptr<T>> cache;
        return cache;
    }

    template <typename T>
    shared_ptr<T> getFile( const string& fn ) {
        static map<string, shared_ptr<T>> cache;
        {
            unique_lock<mutex> lock( fileMutex );
            auto found = cache.find( fn );
            if( found != cache.end() ) {
                return found->second;
            }
        }


    }

}


Format wfwfs::readFmt( const string& filename ) {

    ifstream strm( filename, ifstream::binary );
    if( strm ) {
        uint32_t magic;
        strm.read( reinterpret_cast<char*>( &magic ), sizeof(uint32_t) );
        if( strm.good() && (strm.gcount()==sizeof(uint32_t)) ) {
            switch( magic ) {
//                 case Ana::MAGIC_ANA: ;
//                 case Ana::MAGIC_ANAR: return FMT_ANA;
//                 case FileMomfbd::MAGIC_MOMFBD8: ;	// Fall through
//                 case FileMomfbd::MAGIC_MOMFBD10: ;	// Fall through
//                 case FileMomfbd::MAGIC_MOMFBD11: return FMT_MOMFBD;
#ifdef WFWFS_WITH_FITS
                case Fits::MAGIC_FITS: return FMT_FITS;
#endif
                //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
                default: std::ios_base::failure("readFmt needs to be implemented for this file-type: \"" +filename+"\""); 
            }
        } else throw std::ios_base::failure("Failed to read file: "+filename);
    } else throw std::ios_base::failure("Failed to open file: "+filename);

    return FMT_NONE;

}


Format wfwfs::guessFmt( const string& filename ) {

    size_t pos = filename.find_last_of('.');
    if( pos != string::npos && pos < filename.length() ) {
        string ext = filename.substr(pos+1);
        std::transform(ext.begin(), ext.end(),ext.begin(), ::toupper);
        if( ext == "F0" || ext == "FZ" ) {
            return FMT_ANA;
        } else if( ext == "FITS" ) {
            return FMT_FITS;
        } else if( ext == "MOMFBD" ) {
            return FMT_MOMFBD;
        }
    }
    return FMT_NONE;

}





template <typename T>
void wfwfs::getOrRead( const string& fn, shared_ptr<T>& data ) {

    //static auto& cache = getFileCache<T>();

    cout << "getOrRead(" << fn << ")     "  << endl;

   // future<shared_ptr<T>> async( getFile<T>, fn, 100 );

    /*
        {
            unique_lock<mutex> lock( fileMutex );
            auto found = cache.find( fn );
            if( found != cache.end() ) {
                data = found->second;
                return;
            }
        }*/




}

// template <typename T>
// void wfwfs::getOrRead2( const string& fn, shared_ptr<wfwfs::Image<T>>& im ) {
// //void wfwfs::getOrRead(const string& fn, wfwfs::Image<T>::Ptr& im) {
//     cout << "getOrRead2(" << fn << ")" << endl;
// }

//template void wfwfs::getOrRead( const string&, typename wfwfs::Image<int16_t>::Ptr& );
// template void wfwfs::getOrRead<int32_t>(const string&, typename wfwfs::Image<int32_t>::Ptr&);
// template void wfwfs::getOrRead<float>(const string&, typename wfwfs::Image<float>::Ptr&);
// template void wfwfs::getOrRead<double>(const string&, typename wfwfs::Image<double>::Ptr&);

// template void wfwfs::getOrRead2( const string&, shared_ptr<wfwfs::Image<int16_t>>& );
// template void wfwfs::getOrRead2( const string&, shared_ptr<wfwfs::Image<int32_t>>& );
// template void wfwfs::getOrRead2( const string&, shared_ptr<wfwfs::Image<float>>& );
// template void wfwfs::getOrRead2( const string&, shared_ptr<wfwfs::Image<double>>& );
// 

shared_ptr<wfwfs::FileMeta> wfwfs::getMeta(const string& fn, bool size_only) {

    Format fmt = readFmt(fn);
    shared_ptr<wfwfs::FileMeta> meta;
    
    switch(fmt) {
//         case FMT_ANA: {
//             meta.reset( new Ana(fn) );
//             break;
//         }
#ifdef WFWFS_WITH_FITS
        case FMT_FITS: {
            meta.reset( new Fits(fn) );
            break;
        }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
        default: {
            string msg = "getMeta(arr) needs to be implemented for this file-type: " + to_string(fmt)
                       + "   \"" + fn + "\"";
            throw runtime_error(msg);
        }
    }
    
    return meta;
    
}


void wfwfs::readFile( const string& filename, char* data, shared_ptr<FileMeta>& meta ) {

    try {
        Format fmt = readFmt(filename);
        switch(fmt) {
//             case FMT_ANA: {
//                 if( !meta ) {
//                     meta.reset( new Ana() );
//                 }
//                 shared_ptr<Ana> hdr = static_pointer_cast<Ana>(meta);
//                 if( hdr ) {
//                     Ana::read( filename, data, hdr );
//                 } else {
//                     string msg = "readFile(string,char*,shared_ptr<FileMeta>&) failed to cast meta-pointer into Ana type.";
//                     throw runtime_error(msg);
//                 }
//                 break;
//             }
#ifdef WFWFS_WITH_FITS
            case FMT_FITS: {
                if( !meta ) {
                    meta.reset( new Fits() );
                }
                shared_ptr<Fits> hdr = static_pointer_cast<Fits>(meta);
                if( hdr ) {
                    hdr->read( filename );
                    Fits::read( hdr, data );
                } else {
                    string msg = "readFile(string,char*,shared_ptr<FileMeta>&) failed to cast meta-pointer into Fits type.";
                    throw runtime_error(msg);
                }
                break;
            }
#endif
        //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "readFile(string,char*,shared_ptr<FileMeta>&) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }


}

template <typename T>
void wfwfs::readFile( const string& filename, wfwfs::Array<T>& data ) {
    
    try {
        Format fmt = readFmt(filename);
        switch(fmt) {
//             case FMT_ANA: {
//                 shared_ptr<Ana> hdr(new Ana());
//                 Ana::read(filename,data,hdr);
//                 break;
//             }
#ifdef WFWFS_WITH_FITS
            case FMT_FITS: {
                shared_ptr<Fits> hdr(new Fits());
                Fits::read(filename,data,hdr);
                break;
            }
#endif
            //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "readFile(str,Array<T>) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }

}
template void wfwfs::readFile( const string& filename, wfwfs::Array<uint8_t>& data );
template void wfwfs::readFile( const string& filename, wfwfs::Array<int16_t>& data );
template void wfwfs::readFile( const string& filename, wfwfs::Array<int32_t>& data );
template void wfwfs::readFile( const string& filename, wfwfs::Array<int64_t>& data );
template void wfwfs::readFile( const string& filename, wfwfs::Array<float>& data );
template void wfwfs::readFile( const string& filename, wfwfs::Array<double>& data );


// template <typename T>
// void wfwfs::readFile( const string& filename, wfwfs::Image<T>& image, bool metaOnly ) {
//     
//     try {
//         Format fmt = readFmt(filename);
//         switch(fmt) {
//             case FMT_ANA: Ana::read(filename, image, metaOnly); break;
// #ifdef WFWFS_WITH_FITS
//             case FMT_FITS: Fits::read(filename, image, metaOnly); break;
// #endif
//             //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
//             default: {
//                 string msg = "readFile(string,Image<T>,bool) needs to be implemented for this file-type: " + to_string(fmt)
//                            + "   \"" + filename + "\"";
//                 throw runtime_error(msg);
//             }
//         }
//     } catch ( std::exception& e ) {
//         if( errorHandling == EH_PRINT ) {
//             cerr << "FileIO Error: " << e.what() << endl;
//         } else throw;
//     }
// 
// }
// template void wfwfs::readFile( const string& filename, wfwfs::Image<uint8_t>& data, bool );
// template void wfwfs::readFile( const string& filename, wfwfs::Image<int16_t>& data, bool );
// template void wfwfs::readFile( const string& filename, wfwfs::Image<int32_t>& data, bool );
// template void wfwfs::readFile( const string& filename, wfwfs::Image<int64_t>& data, bool );
// template void wfwfs::readFile( const string& filename, wfwfs::Image<float>& data, bool );
// template void wfwfs::readFile( const string& filename, wfwfs::Image<double>& data, bool );


template <typename T>
void wfwfs::writeFile( const string& filename, wfwfs::Array<T>& data ) {
    
    try {
        Format fmt = guessFmt(filename);
        switch(fmt) {
//             case FMT_ANA: Ana::write(filename, data); break;
#ifdef WFWFS_WITH_FITS
            case FMT_FITS: Fits::write(filename, data); break;
#endif
            //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
            default: {
                string msg = "writeFile(string,Array<T>) needs to be implemented for this file-type: " + to_string(fmt)
                           + "   \"" + filename + "\"";
                throw runtime_error(msg);
            }
        }
    } catch ( std::exception& e ) {
        if( errorHandling == EH_PRINT ) {
            cerr << "FileIO Error: " << e.what() << endl;
        } else throw;
    }


}
template void wfwfs::writeFile( const string& filename, wfwfs::Array<uint8_t>& data );
template void wfwfs::writeFile( const string& filename, wfwfs::Array<int16_t>& data );
template void wfwfs::writeFile( const string& filename, wfwfs::Array<int32_t>& data );
template void wfwfs::writeFile( const string& filename, wfwfs::Array<int64_t>& data );
template void wfwfs::writeFile( const string& filename, wfwfs::Array<float>& data );
template void wfwfs::writeFile( const string& filename, wfwfs::Array<double>& data );


// template <typename T>
// void wfwfs::writeFile( const string& filename, wfwfs::Image<T>& image ) {
//     
//     try {
//         Format fmt = guessFmt(filename);
//         switch(fmt) {
//             case FMT_ANA: Ana::write(filename, image); break;
// #ifdef WFWFS_WITH_FITS
//             case FMT_FITS: Fits::write(filename, image); break;
// #endif
//             //case Ncdf::MAGIC_NCDF: return FMT_NCDF;
//             default: {
//                 string msg = "writeFile(string,Image<T>) needs to be implemented for this file-type: " + to_string(fmt)
//                            + "   \"" + filename + "\"";
//                 throw runtime_error(msg);
//             }
//         }
//     } catch ( std::exception& e ) {
//         if( errorHandling == EH_PRINT ) {
//             cerr << "FileIO Error: " << e.what() << endl;
//         } else throw;
//     }
// 
// 
// }
// template void wfwfs::writeFile( const string& filename, wfwfs::Image<uint8_t>& data );
// template void wfwfs::writeFile( const string& filename, wfwfs::Image<int16_t>& data );
// template void wfwfs::writeFile( const string& filename, wfwfs::Image<int32_t>& data );
// template void wfwfs::writeFile( const string& filename, wfwfs::Image<int64_t>& data );
// template void wfwfs::writeFile( const string& filename, wfwfs::Image<float>& data );
// template void wfwfs::writeFile( const string& filename, wfwfs::Image<double>& data );


/*void wfwfs::loadFiles( const vector<string>& filenames, char* out, size_t frameSize, uint8_t nThreads,
                             double* averages, double* times, string progressMsg ) {
    
    size_t nImages = filenames.size();
    
    atomic<size_t> imgIndex(0);
    
    vector<thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myIndex;
                shared_ptr<FileMeta> myMeta;
                while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                    char* myPtr = out + myIndex*frameSize;
                    try {
                        readFile( filenames[myIndex], myPtr, myMeta );
                        if( times ) times[ myIndex ] = myMeta->getAverageTime().time_of_day().total_nanoseconds()*1E-9;
                    } catch (const exception& e ) {
                        if( !progressMsg.empty() ) printProgress( "\nloadFiles: " + string(e.what()) + " at file #" + to_string(myIndex) + "\n", -1);
                        memset( myPtr, 0, frameSize );  // zero image and continue.
                    }
                    if( !progressMsg.empty() ) printProgress( progressMsg, (myIndex*100.0/(nImages-1)));
                }
            }));
    }
    for (auto& th : threads) th.join();

}*/


void wfwfs::loadFiles( const vector<string>& filenames, char* out, size_t frameSize, uint8_t nThreads,
                             function<void(char*,size_t,shared_ptr<FileMeta>&)> postLoad) {
    
    size_t nImages = filenames.size();
    
    atomic<size_t> imgIndex(0);

    mutex mtx;
    list<string> msgs;
    vector<thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myIndex;
                shared_ptr<FileMeta> myMeta;
                while( (myIndex=imgIndex.fetch_add(1)) < nImages ) {
                    char* myPtr = out + myIndex*frameSize;
                    string fn;
                    {
                        lock_guard<mutex> lock(mtx);
                        fn = filenames[myIndex];
                    }
                    try {
                        readFile( fn, myPtr, myMeta );
                        postLoad( myPtr, myIndex, myMeta );
                    } catch (const exception& e ) {
                        string msg = "file #" + to_string(myIndex) + "\"" + fn + "\": " + e.what();  
                        memset( myPtr, 0, frameSize );  // zero image and continue.
                        lock_guard<mutex> lock(mtx);
                        msgs.push_back( msg );
                    }
                }
            }));
    }
    for (auto& th : threads) th.join();
    if( msgs.size() ) {
        string msg = "FileIO Error: loadFiles(vector<string>,char*,size_t,uint8_t,callback_t)";
        for( auto& m: msgs ) msg += "\n\t  " + m;
        if( errorHandling == EH_PRINT ) {
            cerr << msg << endl;
        } else throw runtime_error( msg );
    }
}


void wfwfs::sumFiles( const std::vector<std::string>& filenames, double* out, size_t frameSize, uint8_t nThreads,
                            preSumCallback preSum ) {
    
/*    
    size_t nImages = filenames.size();
    mutex mtx;

    atomic<size_t> imgIndex(0);
    atomic<size_t> threadIndex(0);
    auto sumFunc = [=,&mtx]( double* a ) {
        std::unique_lock<mutex> lock(mtx);
        for( size_t b=0; b<nPixels; ++b ) summedData[b] += a[b];
    };
    //proc.append( std::bind( sumFunc, std::placeholders::_1 ) );
    
    std::vector<std::thread> threads;
    for( uint8_t i=0; i<nThreads; ++i ) {
        threads.push_back( std::thread(
            [&](){
                size_t myImgIndex;
                size_t myThreadIndex = threadIndex.fetch_add(1);
                double* mySumPtr = sumPtr+myThreadIndex*nPixels;
                double* myTmpPtr = tmpPtr+myThreadIndex*nPixels;
                while( (myImgIndex=imgIndex.fetch_add(1)) < nImages ) {
                    initFuncs[myImgIndex](myTmpPtr,myThreadIndex);
                    for( size_t i=0; i<nPixels; ++i ) mySumPtr[i] += myTmpPtr[i];
                    if( kw.verbose > 1 ) printProgress( statusString, (myImgIndex*100.0/(nImages-1)));
                }
                sumFunc(mySumPtr);
            }));
    }

    for (auto& th : threads) th.join();
*/
}
