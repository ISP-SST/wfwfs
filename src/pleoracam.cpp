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
#include "pleoracam.hpp"

#ifdef WFWFS_WITH_PLEORA

#include "datautil.hpp"

#include <emmintrin.h>
#include <tmmintrin.h>

#include <PvSystem.h>
#include <PvDeviceGEV.h>
#include <PvStreamGEV.h>

#include <iostream>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/thread/thread.hpp>

using namespace wfwfs;
using namespace std;

namespace bpx = boost::posix_time;


namespace {
    
    
    __attribute__ ((target ("default")))
    static void baslerToImg10( uint8_t *data, size_t width, size_t height ) {
        // Basler format: AA ba BB (AA, BB = high nibbles, a, b = low nibble)
        size_t n = width * height / 2;
        uint8_t *p = (data + (width * height * 3) / 2) - 3;
        uint16_t *q = (uint16_t *)(data + (width * height * 2)) - 2;
        for(size_t i = 0; i < n; i++) {
            q[0] = (p[0] << 2) | (p[1] & 0x3);
            q[1] = (p[2] << 2) | ((p[1] >> 4) & 0x3);
            p -= 3;
            q -= 2;
        }
    }


    __attribute__ ((target ("ssse3")))
    static void baslerToImg10( uint8_t *data, size_t width, size_t height ) {
        // Basler format: AA ba BB (AA, BB = high nibbles, a, b = low nibble)
        size_t n = width * height / 8; // We're going to process 12 bytes = 8 pixels per step

        // We have to start at the end
        uint8_t *p = data + (width * height * 3) / 2;
        __m128i *q = (__m128i *)(data + width * height * 2);

        // Masks to pick the right bytes. -1 means copy in a zero.
        const __m128i sh = _mm_setr_epi8(-1,  0, -1,  2, -1,  3, -1,  5, -1,  6, -1,  8, -1,  9, -1, 11);
        const __m128i sl = _mm_setr_epi8( 1, -1, -1, -1,  4, -1, -1, -1,  7, -1, -1, -1, 10, -1, -1, -1);

        for(size_t i = 0; i < n; i++) {
            p -= 12;
            q -= 1;
            __m128i a = _mm_loadu_si128((__m128i *)p); // load 12 bytes (= 8x 12-bit)                   AA ba BB CC dc DD .. ..
            __m128i b = _mm_slli_epi16(a, 6);          // make a shifted copy to align the low nibble   .. a. .. .. c. .. .. ..
            __m128i c = _mm_slli_epi16(a, 2);          // make a shifted copy to align the low nibble   .. b. .. .. d. .. .. ..
            __m128i h = _mm_shuffle_epi8(a, sh);       // put the high nibbles in the right place       00 AA 00 BB 00 CC 00 DD
            __m128i l = _mm_shuffle_epi8(b, sl);       // first set of low nibbles                      a. 00 00 00 c. 00 00 00
            __m128i k = _mm_shuffle_epi8(c, sl);       // second set of low nibbles                     00 00 b. 00 00 00 d. 00
            __m128i x = _mm_or_si128(h, _mm_or_si128(l, k)); // or together                             a. AA b. BB c. CC d. DD
            *q        = _mm_srli_epi16(x, 6);          // shift right every 16 bit word by 6 bits       Aa 0A Bb 0B Cc 0C Dd 0D
        }
    }



    __attribute__ ((target ("default")))
    static void baslerToImg12( uint8_t *data, size_t width, size_t height ) {
        // Basler format: AA ba BB (AA, BB = high nibbles, a, b = low nibble)
        size_t n = width * height / 2;
        uint8_t *p = (data + (width * height * 3) / 2) - 3;
        uint16_t *q = (uint16_t *)(data + (width * height * 2)) - 2;
        for(size_t i = 0; i < n; i++) {
            q[0] = (p[0] << 4) | (p[1] & 0xf);
            q[1] = (p[2] << 4) | (p[1] >> 4);
            p -= 3;
            q -= 2;
        }
    }


    __attribute__ ((target ("ssse3")))
    static void baslerToImg12( uint8_t *data, size_t width, size_t height ) {
        // Basler format: AA ba BB (AA, BB = high nibbles, a, b = low nibble)
        size_t n = width * height / 8; // We're going to process 12 bytes = 8 pixels per step

        // We have to start at the end
        uint8_t *p = data + (width * height * 3) / 2;
        __m128i *q = (__m128i *)(data + width * height * 2);

        // Masks to pick the right bytes. -1 means copy in a zero.
        const __m128i sh = _mm_setr_epi8(-1,  0, -1,  2, -1,  3, -1,  5, -1,  6, -1,  8, -1,  9, -1, 11);
        const __m128i sl = _mm_setr_epi8( 1, -1, -1, -1,  4, -1, -1, -1,  7, -1, -1, -1, 10, -1, -1, -1);
        const __m128i sk = _mm_setr_epi8(-1, -1,  1, -1, -1, -1,  4, -1, -1, -1,  7, -1, -1, -1, 10, -1);

        for(size_t i = 0; i < n; i++) {
            p -= 12;
            q -= 1;
            __m128i a = _mm_loadu_si128((__m128i *)p); // load 12 bytes (= 8x 12-bit)                   AA ba BB CC dc DD .. ..
            __m128i b = _mm_slli_epi16(a, 4);          // make a shifted copy to align the low nibble   .. a. .. .. c. .. .. ..
            __m128i h = _mm_shuffle_epi8(a, sh);       // put the high nibbles in the right place       00 AA 00 BB 00 CC 00 DD
            __m128i l = _mm_shuffle_epi8(b, sl);       // first set of low nibbles                      a. 00 00 00 c. 00 00 00
            __m128i k = _mm_shuffle_epi8(a, sk);       // second set of low nibbles                     00 00 b. 00 00 00 d. 00
            __m128i x = _mm_or_si128(h, _mm_or_si128(l, k)); // or together                             a. AA b. BB c. CC d. DD
            *q        = _mm_srli_epi16(x, 4);          // shift right every 16 bit word by 4 bits       Aa 0A Bb 0B Cc 0C Dd 0D
        }
    }


}


PleoraCam::PleoraCam( boost::property_tree::ptree& pt ) : Camera(pt),
                                                          lDevice(nullptr),
                                                          lDeviceParams(nullptr),
                                                          lStream(nullptr) {

    string path = pt.get<string>( "camera_type", "" );

    info.model = "Pleora camera";
    info.firmware_version = "TODO";
    info.serial = "12345";

 
}


PleoraCam::~PleoraCam( void ) {

}


void PleoraCam::init( void ) {
    
    PvSystem lSystem;
    PvResult lResult;
    
    //lSystem.SetDetectionTimeout(10);
    if( !maybeLog( lSystem.Find(), "Find" ) ) {
        abort();
    }

    // Detect, select device.
    vector<const PvDeviceInfo *> lDIVector;
    for ( uint32_t i(0); i<lSystem.GetInterfaceCount(); i++ ) {
        const PvInterface *lInterface = lSystem.GetInterface( i );
        if ( lInterface != NULL ) {
            for ( uint32_t j = 0; j < lInterface->GetDeviceCount(); j++ ) {
                const PvDeviceInfo *lDI = lInterface->GetDeviceInfo( j );
                if ( lDI != NULL ) {
                    lDIVector.push_back( lDI );
                }
            }
        }
    }

    if ( lDIVector.empty() ) {
        fprintf( stderr, "pleora: Failed to detect device(s).\n" );
        abort();
    }

    const PvDeviceInfo *lSelectedDI = nullptr;
    for( auto& di: lDIVector ) {
        PvString lConnectionID = di->GetConnectionID();
        lDevice = PvDevice::CreateAndConnect( lConnectionID, &lResult );
        if ( !lResult.IsOK() ) {
            PvDevice::Free( lDevice );
            lDevice = nullptr;
        } else {
            lSelectedDI = di;
            break;
        }
    }
    
    if ( lSelectedDI == nullptr ) {
        fprintf( stderr, "pleora: Failed to connect to any device.\n" );
        abort();
    }
    
    info.model = lSelectedDI->GetVendorName().GetAscii();
    info.model += ":";
    info.model += lSelectedDI->GetModelName().GetAscii();
    info.serial = lSelectedDI->GetSerialNumber().GetAscii();
    info.firmware_version = lSelectedDI->GetVersion().GetAscii();
    
    lDeviceParams = lDevice->GetParameters();
//     for( uint32_t i(0); i<lDeviceParams->GetCount(); ++i ) {
//         PvGenParameter* par = lDeviceParams->Get(i);
//         if( !par) continue;
//     cout << par->GetName().GetAscii() << " - " << par->ToString().GetAscii() << endl;
//     }
    
    int64_t tmp64, valMin, valMax;

    maybeLog( lDeviceParams->GetIntegerRange( "Width", valMin, valMax ), "getting Width-range" );
    cfg.width = restrictValue( cfg.width, valMin, valMax );
    maybeLog( lDeviceParams->SetIntegerValue( "Width", cfg.width ), "setting Width" );
    if( maybeLog( lDeviceParams->GetIntegerValue( "Width", tmp64 ), "getting Width" ) ) {
        cfg.width = tmp64;
    }
    int64_t maxWidth = valMax;
    
    maybeLog( lDeviceParams->GetIntegerRange( "Height", valMin, valMax ), "getting Height-range" );
    cfg.height = restrictValue( cfg.height, valMin, valMax );
    maybeLog( lDeviceParams->SetIntegerValue( "Height", cfg.height ), "setting Height" );
    if( maybeLog( lDeviceParams->GetIntegerValue( "Height", tmp64 ), "getting Height" ) ) {
        cfg.height = tmp64;
    }
    int64_t maxHeight = valMax;

//     tmp64 = config.getint( "offset_x", 0 );
//     maybeLog( lDeviceParams->GetIntegerRange( "OffsetX", valMin, valMax ), "getting OffsetX-range" );
//     tmp64 = limit( tmp64, valMin, valMax );
//     maybeLog( lDeviceParams->SetIntegerValue( "OffsetX", tmp64 ), "setting OffsetX" );
// 
//     tmp64 = config.getint( "offset_y", 0 );
//     maybeLog( lDeviceParams->GetIntegerRange( "OffsetY", valMin, valMax ), "getting OffsetY-range" );
//     tmp64 = limit( tmp64, valMin, valMax );
//     maybeLog( lDeviceParams->SetIntegerValue( "OffsetY", tmp64 ), "setting OffsetY" );
    
    switch( cfg.depth ) {
        case 8: maybeLog( lDeviceParams->SetEnumValue( "PixelFormat", "Mono8" ), "setting PixelFormat:Mono8" ); break;
        case 10: maybeLog( lDeviceParams->SetEnumValue( "PixelFormat", "Mono10Packed" ), "setting PixelFormat:Mono10Packed" ); break;
        case 12: maybeLog( lDeviceParams->SetEnumValue( "PixelFormat", "Mono12Packed" ), "setting PixelFormat:Mono12Packed" ); break;
        default: maybeLog( lDeviceParams->SetEnumValue( "PixelFormat", "Mono8" ), "setting PixelFormat:Mono8" );
        fprintf( stderr, "pleora: unrecognized bit-depth (%d). Using Mono8!\n", cfg.depth );
        cfg.depth = 8;
    }
    
    //syslog( LOG_INFO, "Found %zux%zu GigE camera, running at %zux%zu", maxWidth, maxHeight, cfg.width, cfg.height );
cout << boost::format("Found %dx%d GigE camera, running at %dx%d") % maxWidth % maxHeight % cfg.width % cfg.height << endl;
// cout << boost::str( boost::format("Found %zux%zu GigE camera, running at %zux%zu")
//         % maxWidth % maxHeight % cfg.width % cfg.height ) << endl;
    
    PvString tmpS;
    lDeviceParams->GetEnumValue( "PixelFormat", tmpS );

    lStream = PvStream::CreateAndOpen( lSelectedDI->GetConnectionID(), &lResult );
    maybeLog( lResult, "creating stream object" );
    if( lStream == nullptr ) {
        fprintf( stderr, "pleora: Failed to create stream object.\n" );
        abort();
    }
    
   // If this is a GigE Vision device, configure GigE Vision specific streaming parameters
    PvDeviceGEV* lDeviceGEV = dynamic_cast<PvDeviceGEV *>( lDevice );
    if ( lDeviceGEV != NULL ) {
        PvStreamGEV *lStreamGEV = static_cast<PvStreamGEV *>( lStream );
        lDeviceGEV->NegotiatePacketSize();
        lDeviceGEV->SetStreamDestination( lStreamGEV->GetLocalIPAddress(), lStreamGEV->GetLocalPort() );
    }

    CreateStreamBuffers();
    
    //syslog( LOG_INFO, "PixelFormat: %s  nFrames: %ld", tmpS.GetAscii(), nframes );
cout << boost::format("PixelFormat: %s  nFrames: %ld") % tmpS.GetAscii() % cfg.n_buf << endl;
    
    for( auto& b: lBuffers ) {
        lStream->QueueBuffer( b.get() );
    }


    maybeLog( lDevice->StreamEnable(), "enabling stream" );
    maybeLog( lDeviceParams->ExecuteCommand( "AcquisitionStart" ), "starting acquisition" );
   
    //camera_get_temp();
    fprintf( stderr, "Camera temp: %f degrees.\n", state.temp );
    

}


void PleoraCam::cleanup( void ) {
    
    running = false;
    
    maybeLog( lDeviceParams->ExecuteCommand( "AcquisitionStop" ), "stopping acquisition" );
    maybeLog( lDevice->StreamDisable(), "disabling stream" );
    maybeLog( lStream->AbortQueuedBuffers(), "aborting queued buffers" );
    while ( lStream->GetQueuedBufferCount() > 0 ) {
        PvBuffer *lBuffer = NULL;
        PvResult lOperationResult;
        lStream->RetrieveBuffer( &lBuffer, &lOperationResult );
    }
    
    if( lStream ) {
        maybeLog( lStream->Close(), "closing stream" );
        PvStream::Free( lStream );
        lStream = nullptr;
    }
    
    if( lDevice ) {
        maybeLog( lDevice->Disconnect(), "disconnecting device" );
        PvDevice::Free( lDevice );
        lDevice = nullptr;
    }
    
    lBuffers.clear();
    
 
}


void PleoraCam::run( frame_cb_t cb ) {
    
    running = true;
    if( cb ) qcb = cb;

    update_temp();
    bpx::ptime last_temp_update = bpx::microsec_clock::universal_time();
    const bpx::time_duration temp_interval( 0, 0, 1, 0 );

    while( running ) {
        boost::this_thread::interruption_point();
        PvBuffer *lBuffer = NULL;
        PvResult lOperationResult;
        uint8_t* data = nullptr;
        // Retrieve next buffer
        if( maybeLog( lStream->RetrieveBuffer( &lBuffer, &lOperationResult, 10000 ), "retreive buffer" ) ) {
            if( running && maybeLog( lOperationResult, "RetrieveBuffer" ) ) {
                PvPayloadType lType = lBuffer->GetPayloadType();
                if( lType == PvPayloadTypeImage ) {
                    PvImage *lImage = lBuffer->GetImage();
                    data = lImage->GetDataPointer();
                }
            } 
        } 
        bpx::ptime now = bpx::microsec_clock::universal_time();
        if( now > last_temp_update+temp_interval ) {
            update_temp();
            last_temp_update = now;
        }
        if( !lBuffer || !data ) {
            // failure, try again after 5 ms
            usleep(5000);
        } else {
            if( cfg.depth == 10 ) {
                baslerToImg10( data, cfg.width, cfg.height );
            } else if( cfg.depth == 12 ) {
                baslerToImg12( data, cfg.width, cfg.height );
            }
            if( qcb ) {     // if we have a callback defined, execute it
                qcb( data, now );
            }
        }
        if( lBuffer ) lStream->QueueBuffer( lBuffer );
    }

}


void PleoraCam::stop( void ) {
    
    running = false;
    
}


double PleoraCam::get_exposure( void ) const {
    double exposure = 0;
    if( lDevice && lDeviceParams ) {
        lock_guard<mutex> h(mtx);
        double micro_exp;
        if( maybeLog( lDeviceParams->GetFloatValue( "ExposureTime", micro_exp ), "get_exposure" ) ) {
            exposure = micro_exp * 1e-6;
        }
    }
    return exposure;
}


void PleoraCam::set_exposure( double value ) {
    if( lDevice && lDeviceParams ) {
        lock_guard<mutex> h(mtx);
        double micro_exp = value * 1000000.0;
        if( maybeLog( lDeviceParams->SetFloatValue( "ExposureTime", micro_exp ), "set_exposure" ) ) {
            state.exposure = value;
        }
    }
}


void PleoraCam::set_interval( double value ) {

    if( lDevice && lDeviceParams ) {
        lock_guard<mutex> h(mtx);
        if( isfinite(value) && value > 0 ) {
            if( state.interval == 0.0 ) lDeviceParams->SetBooleanValue( "AcquisitionFrameRateEnable", true );
            double maxFPS;
            lDeviceParams->GetFloatValue( "AcquisitionFrameRateMax", maxFPS );
            double fps = std::min( 1.0 / value, maxFPS );
            if( maybeLog( lDeviceParams->SetFloatValue( "AcquisitionFrameRate", fps ), "set_interval" ) ) {
                state.interval = 1.0 / fps;
            }
        } else {
            lDeviceParams->SetBooleanValue( "AcquisitionFrameRateEnable", false );
            state.interval = 0.0;
        }
    }
}

void PleoraCam::set_gain(double value) {
	lock_guard<mutex> h(mtx);
// 	tPvUint32 raw = value;
// 	PvAttrUint32Set(handle, "GainValue", raw);
// 	PvAttrUint32Get(handle, "GainValue", &raw);
// 	gain = raw;
}


void PleoraCam::copyRawData( const uint8_t* in, uint8_t* out, uint8_t outDepth ) const {
    
    int shift = 0;
    size_t nEl = cfg.width*cfg.height;
    size_t totalSize = nEl;
    if( cfg.depth > 8 ) totalSize *= 2;
    if( cfg.depth == outDepth ) {
        memcpy( out, in, totalSize );
    } else if( outDepth > cfg.depth ) {
        shift = (outDepth - cfg.depth);
        if( cfg.depth > 8 ) {
            const uint16_t* in16 = reinterpret_cast<const uint16_t*>(in);
            std::transform( in16, in16+totalSize, reinterpret_cast<uint16_t*>(out), [shift]( uint16_t a ) { return a<<shift; } );
        } else {
            if( outDepth > 8 ) {
                std::transform( in, in+totalSize, reinterpret_cast<uint16_t*>(out), [shift]( uint8_t a ) { return static_cast<uint16_t>(a)<<shift; } );
            } else {
                std::transform( in, in+totalSize, out, [shift]( uint8_t a ) { return a<<shift; } );
            }
        }
    } else {
        shift = (cfg.depth - outDepth);
        if( cfg.depth > 8 ) {
            const uint16_t* in16 = reinterpret_cast<const uint16_t*>(in);
            if( outDepth > 8 ) {
                std::transform( in16, in16+totalSize, reinterpret_cast<uint16_t*>(out), [shift]( uint16_t a ) { return a>>shift; } );
            } else {
                std::transform( in16, in16+totalSize, out, [shift]( uint16_t a ) { return a>>shift; } );
            }
        } else {
            std::transform( in, in+totalSize, out, [shift]( uint8_t a ) { return a>>shift; } );
        }
    }

}


void PleoraCam::update_temp( void ) {

    float temp = 0.0/0.0;
    if( lDevice && lDeviceParams ) {
        lock_guard<mutex> h(mtx);
        double tmpD;
        if( maybeLog( lDeviceParams->GetFloatValue( "DeviceTemperature", tmpD ), "get_temp" ) ) {
            temp = tmpD;
        }
    }
    
    if( isfinite( temp ) ) {
        state.temp = temp;
    }

}


bool PleoraCam::maybeLog( const PvResult& res, string descr ) const {
    if( !res.IsOK() ) {
        const string msg = "pleora: " + descr + string(" failed. Reason: ") + res.GetCodeString().GetAscii();
        //syslog( LOG_NOTICE, "%s", msg.c_str() );
        cerr << "Pleora:notice: " << msg << endl;
        //lock_guard<mutex> lock(mtx);
        msgList.push_back( msg );
        return false;
    }
    return true;
}


void PleoraCam::CreateStreamBuffers( void ) {
    uint32_t lSize = cfg.width*cfg.height*((cfg.depth==8)?1:2);
    //lock_guard<mutex> lock(mtx);
    cfg.n_buf = std::min<uint64_t>( lStream->GetQueuedBufferMaximum(), cfg.n_buf );
    for( int i(0); i < cfg.n_buf; ++i ) {
        shared_ptr<PvBuffer> lBuffer(new PvBuffer);
        lBuffer->Alloc( static_cast<uint32_t>( lSize ) );
        lBuffers.push_back( lBuffer );
    }
}




#endif // WFWFS_WITH_PLEORA
