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
#include "model.hpp"

#include "arraystats.hpp"
#include "utils.hpp"

using namespace wfwfs;
using namespace std;


Model& Model::get(void) {
    static Model model;
    return model;
}


void Model::generate( PointI frameSize, float scale, PointF opticalAxis ) {
    get().generate_( frameSize, scale, opticalAxis );
}


template <typename T>
void Model::make_mask( const Array<T>& ff, bool save ) {
    get().make_mask_( ff, save );
}
template void Model::make_mask( const Array<float>& ff, bool );


void Model::generate_( PointI frameSize, float scale, PointF opticalAxis ) {
    
    lock_guard<mutex> lock(mtx);
    subfields.clear();
    const int N = 2*frameSize.max()/scale;      // Just something that makes sure we cover the frame

    double max_radial_distance(0);              // Limit the subfields to the "rings" that fit on the middle horizontal row.
    int i=0;
    while( 2*i*base1.x*scale < frameSize.x ) {
        max_radial_distance = i++*base1.x*scale*1.05;
    }
    
    for( int i=-N; i<=N; ++i ) {
        for( int j=-N; j<=N; ++j ) {
            PointF pos = base1*i + base2*j;
            pos *= scale;
            if( sqrt(pos.norm()) < max_radial_distance ) {
                pos += opticalAxis;
                PointI nc(i,j);
                subfields.emplace( nc, SF_Info(nc,pos) );
            }
        }
    }

}


template <typename T>
void Model::make_mask_( const Array<T>& ff, bool save ) {
    
    size_t sizeY = ff.dimSize(0);
    size_t sizeX = ff.dimSize(1);
    size_t nElements = ff.nElements();
    
    if( save ) mask.clear();
    
    if( !nElements ) {
        return;
    }
    
    Array<T> tmpFF;
    ff.copy(tmpFF);
    T* ffPtr = tmpFF.get();
    
    Array<uint8_t> tmpMask( sizeY, sizeX );
    uint8_t* tmpMaskPtr = tmpMask.get();
    tmpMask.zero();

    Array<uint8_t> tmpMask2( sizeY, sizeX );
    uint8_t* tmpMaskPtr2 = tmpMask2.get();
    tmpMask2.zero();
    
    ArrayStats stats;
    stats.getMedian( ffPtr, nElements );
    
    if( stats.min == stats.max ) {
        return;
    }
    
    // a hard threshold to suppress noise and weak illumination areas.
    double thres = 0.75*stats.mean;
    std::transform( ffPtr, ffPtr+nElements, ffPtr, [&]( const T& a ) {
        if( a < thres ) return T(0);
        return a;
    });

    // apply an adaptive threshold to detect prominent features.
    double maxval = 1.0;
    int blockSize = 13;
    double offset = 2.0;
    threshold( ffPtr, tmpMaskPtr, sizeY, sizeX, cv::ADAPTIVE_THRESH_MEAN_C, cv::THRESH_BINARY_INV, maxval, blockSize, offset );

    // Invert mask.
    // At this point the mask is a black image with white frames around the features, inverting it allows us to use maskConnected below.
    std::transform( tmpMaskPtr, tmpMaskPtr+nElements, tmpMaskPtr, [&]( const uint8_t& a ) {
        return !a;
    });

    PointD avg_size(0);
    size_t nFittedSubfields(0);
    const double minArea(10000);
    vector<cv::Point2f> refPoints, targetPoints;

    for( auto& s: subfields ) {
        PointF pos = s.second.nominal_position;
        PointI posI = pos + PointF(0.5);  // round coordinates to nearest integer.
        tmpMask2.zero();
        maskConnected2( tmpMaskPtr, tmpMaskPtr2, sizeY, sizeX, posI.y, posI.x, uint8_t(1) );
        wfwfs::morph_smooth( tmpMaskPtr2, sizeY, sizeX, 7 , false );
        PointI first(0),last(0);
        bounding_rect( tmpMaskPtr2, sizeY, sizeX, first, last );
        PointF sz = (last-first)+1;
        if( sz.area() > minArea ) {
            refPoints.push_back( cv::Point2f( pos.x, pos.y ) );
            s.second.real_position = first + sz/2;
            targetPoints.push_back( cv::Point2f( s.second.real_position.x, s.second.real_position.y ) );
            avg_size += sz;
            nFittedSubfields++;
        }
    }

    if( nFittedSubfields > 3 ) {
        avg_size /= nFittedSubfields;
    } else {    // Not enough subfields detected/fitted, return.
        return;
    }

    subfield_size = avg_size-10;        // Add some margin around to avoid the edges of the subfields (which has large gradients) 
    
    // Find affine transformation that maps nominal to measured coordinates.
    cv::Mat H = findHomography( refPoints, targetPoints, CV_LMEDS );    // CV_RANSAC
    
    targetPoints.clear();
    refPoints.clear();
    for( auto& s: subfields ) {     // clear & re-add so that also subfields not detected above gets mapped
        refPoints.push_back( cv::Point2f(s.second.nominal_position.x,s.second.nominal_position.y) );
    }
    
    try {
        // Map nominal coordinates according to new/fitted model
        perspectiveTransform( refPoints, targetPoints, H );
    } catch( ... ) {
        targetPoints.clear();
    }

    if( save && (subfields.size() == targetPoints.size()) ) {
        if( !tmpFF.sameSizes(mask) ) {
            mask.resize( tmpFF.dimensions(true) );
        }
        mask.zero();
        uint8_t* maskPtr = mask.get();
        size_t cnt(0);
        for( auto& s: subfields ) {
            PointF& pos = s.second.real_position;
            pos.x = targetPoints[cnt].x;
            pos.y = targetPoints[cnt++].y;
            PointI offset = pos - subfield_size/2;
            if( maskPtr ) {
                for( int yy(0); yy<=subfield_size.y; ++yy ) {
                    for( int xx(0); xx<=subfield_size.x; ++xx ) {
                        maskPtr[ (offset.y+yy)*sizeX + offset.x+xx ] = 1;
                    }
                }
            }
        }
    }
    
}
template void Model::make_mask_( const Array<float>& ff, bool );


template <typename T>
void Model::draw_subfields( Array<T>& img, const T value ) const {

    lock_guard<mutex> lock(mtx);
    PointI hs = subfield_size/2;
    PointI hs2 = hs*PointI(1,-1);
    for( const auto& s: subfields ) {
        PointI pos = s.second.nominal_position;
        if( s.second.real_position.max() > 0 ) {
            pos = s.second.real_position;
        }
        PointI c1 = pos + hs;
        PointI c2 = pos + hs2;
        PointI c3 = pos - hs;
        PointI c4 = pos - hs2;
        draw_line( img, c1, c2, value, 2 );
        draw_line( img, c1, c4, value, 2 );
        draw_line( img, c3, c2, value, 2 );
        draw_line( img, c3, c4, value, 2 );
    }
}
template void Model::draw_subfields( Array<uint8_t>& img, const uint8_t ) const;
