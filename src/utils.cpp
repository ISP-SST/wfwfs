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
#include "utils.hpp"

#include "arraystats.hpp"

#include <functional>
#include <map>
#include <set>
#include <math.h>
#include <queue>

#include <gsl/gsl_multifit.h>

#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>


using namespace wfwfs;
using namespace std;

namespace {

    const int maxDistance = 64;
    const int maxDistance2 = maxDistance*maxDistance;
    const double deltasqr = 4;
    const double beta = 2;
    
    double distMap[2*maxDistance2+1];
    
    const double* getDistanceMap (void) {
        memset(distMap,0,(2*maxDistance2+1)*sizeof(double));
        for(int i=0; i<=2*maxDistance2; ++i) distMap[i] = pow (i + deltasqr, -beta);
        return distMap;
    }

    template <typename T>
    void maskConnected( const T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int startY, unsigned int startX, T threshold) {
        if (startY >= sizeY || startX >= sizeX) return;
        if ( !mask[startY][startX] && (data[startY][startX] > threshold)) {       // new point inserted, check neighbours
            mask[startY][startX] = 1;
            maskConnected(data, mask, sizeY, sizeX, startY + 1, startX, threshold);
            maskConnected(data, mask, sizeY, sizeX, startY - 1, startX, threshold);
            maskConnected(data, mask, sizeY, sizeX, startY, startX + 1, threshold);
            maskConnected(data, mask, sizeY, sizeX, startY, startX - 1, threshold);
        }
    }

}


template <typename T>
void wfwfs::maskConnected2( const T* data, uint8_t* mask, size_t sizeY, size_t sizeX, unsigned int startY, unsigned int startX, T threshold ) {
    
    if ( startY >= sizeY || startX >= sizeX ) return;
    
    Array<uint8_t> tmpMask( sizeY, sizeX );
    uint8_t* tmpPtr = tmpMask.get();
    
    size_t nEl = sizeX*sizeY;
    std::transform( data, data+nEl, tmpPtr, [=](const T a) {
        if( a < threshold ) return 0;
        return 1;
    });
    
    std::queue<size_t> offsets;
    offsets.push( startY*sizeX+startX );
    tmpPtr[ startY*sizeX+startX ] = 2;     // always set starting point.
    while( !offsets.empty() ) {
        size_t offset = offsets.front();
        offsets.pop();
        if( tmpPtr[ offset ] != 2 ) continue;
        tmpPtr[ offset ] = 2;
        if( (offset+1 < nEl) && tmpPtr[ offset+1 ] == 1 ) { offsets.push( offset+1 ); tmpPtr[ offset+1 ] = 2; }
        if( (offset > 0) && tmpPtr[ offset-1 ] == 1 ) { offsets.push( offset-1 ); tmpPtr[ offset-1 ] = 2; }
        if( (offset+sizeX < nEl) && tmpPtr[ offset+sizeX ] == 1 ) { offsets.push( offset+sizeX ); tmpPtr[ offset+sizeX ] = 2; }
        if( (offset >= sizeX) && tmpPtr[ offset-sizeX ] == 1 ) { offsets.push( offset-sizeX ); tmpPtr[ offset-sizeX ] = 2; }
    }
    
    std::transform( tmpPtr, tmpPtr+nEl, mask, [](const uint8_t a) {
        if( a == 2 ) return 1;
        return 0;
    });
    
}
template void wfwfs::maskConnected2( const float*, uint8_t*, size_t, size_t, unsigned int, unsigned int, float);
template void wfwfs::maskConnected2( const uint8_t*, uint8_t*, size_t, size_t, unsigned int, unsigned int, uint8_t);


template <typename T>
void wfwfs::apodize( T* data, size_t n, const T& target ) {

    double mean = 0.5 * static_cast<double>(*data + target);
    double amplitude = 0.5 * static_cast<double>(*data - target);
    double step = M_PI / (n-1);

    for( size_t i = 1; i < n; ++i ) {
        *(data + i) = mean + amplitude * cos( i * step );
    }

}
template void wfwfs::apodize( int16_t*, size_t, const int16_t&);
template void wfwfs::apodize( int32_t*, size_t, const int32_t&);
template void wfwfs::apodize( float*, size_t, const float&);
template void wfwfs::apodize( double*, size_t, const double&);


template <typename T>
vector<double> wfwfs::fitPlane( const T* inPtr, size_t sizeY, size_t sizeX, size_t strideY ) {
    
    if( strideY == 0 ) strideY = sizeX;     //  default is to assume a dense datablock, i.e. y-stride = x-size.
    
    int n = sizeY * sizeX;
    int nParams = 3;                        //  fit a plane as:   z = a*x + b*y + c

    gsl_vector *data = gsl_vector_alloc( n );
    gsl_vector *coeff = gsl_vector_alloc( nParams );
    gsl_matrix *X = gsl_matrix_alloc( n, nParams );
    gsl_matrix *covar = gsl_matrix_alloc( nParams, nParams );

    for( size_t i = 0; i < sizeY; ++i ) {
        for( size_t j = 0; j < sizeX; ++j ) {
            int inOffset = i*strideY + j;
            int xOffset = i*sizeX + j;
            //cout << "i = " << i  << "   j = " << j  << "   offset = " << xOffset << endl;
            gsl_matrix_set( X, xOffset, 0, j );
            gsl_matrix_set( X, xOffset, 1, i );
            gsl_matrix_set( X, xOffset, 2, 1 );
            gsl_vector_set( data, xOffset, static_cast<double>(inPtr[inOffset]) );
        }
    }

    double chisq;
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, nParams);
    gsl_multifit_linear (X, data, coeff, covar, &chisq, work);
    gsl_multifit_linear_free (work);
    gsl_vector_free (data);
    gsl_matrix_free (X);
    gsl_matrix_free (covar);

    vector<double> ret( coeff->data, coeff->data+nParams );
    gsl_vector_free (coeff);

    return ret;
}
template vector<double> wfwfs::fitPlane ( const uint8_t* , size_t, size_t, size_t );
template vector<double> wfwfs::fitPlane ( const uint16_t* , size_t, size_t, size_t );
template vector<double> wfwfs::fitPlane ( const int16_t* , size_t, size_t, size_t );
template vector<double> wfwfs::fitPlane ( const int32_t* , size_t, size_t, size_t );
template vector<double> wfwfs::fitPlane ( const float* , size_t, size_t, size_t );
template vector<double> wfwfs::fitPlane ( const double* , size_t, size_t, size_t );


        
template <typename T>
wfwfs::Array<T> wfwfs::fitPlane (const wfwfs::Array<T>& in, bool subtract_mean, double* coeffs) {

    int ySize = in.dimSize(0);
    int xSize = in.dimSize(1);

    vector<double> cc = fitPlane( in.ptr(), ySize, xSize );

    double a = cc[0];
    double b = cc[1];
    double c = cc[2];

    if (subtract_mean) {
        c = 0;                                          // ignore c-coefficient (mean), to just fit the tilts, not the offset.
    }

    if( coeffs ) {
        coeffs[0] = a;
        coeffs[1] = b;
        coeffs[2] = c;
    }

    wfwfs::Array<T> ret(ySize, xSize);
    T* retPtr = ret.ptr();
    for (int i = 0; i < ySize; ++i) {
        double y = static_cast<double>(i)+0.5;
        for (int j = 0; j < xSize; ++j) {
            double x = static_cast<double>(j)+0.5;
            int offset = i * xSize + j;
            retPtr[offset] = a*x + b*y + c;
        }
    }
    
    return ret;
}
template wfwfs::Array<float> wfwfs::fitPlane (const wfwfs::Array<float>&, bool, double*);
template wfwfs::Array<double> wfwfs::fitPlane (const wfwfs::Array<double>&, bool, double*);


template <typename T>
void wfwfs::connectedRegion( const T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int y_start, unsigned int x_start, T threshold) {
    memset(*mask, 0, sizeY * sizeX);
    maskConnected( data, mask, sizeY, sizeX, y_start, x_start, threshold );
}
template void wfwfs::connectedRegion(const float**, uint8_t**, size_t, size_t, uint, uint, float);
template void wfwfs::connectedRegion(const double**, uint8_t**, size_t, size_t, uint, uint, double);


template <typename T>
void wfwfs::smooth(T** data, size_t sizeY, size_t sizeX, size_t nY, size_t nX) {
    if (nY == 0 || nX == 0) return;
    T** tmp = newArray<T> (sizeY, sizeX);
    size_t n = sizeY*sizeX;
    std::fill_n( *tmp, n, T(0) );
    for (unsigned int y = 0; y < sizeY; ++y) {
        int yl = std::max<int>(y-nY, 0);
        int yh = std::min(y+nY, sizeY);
        for (unsigned int x = 0; x < sizeX; ++x) {
            int xl = std::max<int>(x-nX, 0);
            int xh = std::min(x+nX, sizeX);
            int cnt = 0;
            for (int yy=yl; yy < yh; ++yy) {
                for (int xx=xl; xx < xh; ++xx) {
                    tmp[y][x] += data[yy][xx];
                    cnt++;
                }
            }
            tmp[y][x] /= cnt;
        }
    }
    std::copy_n( *tmp, n, *data );
    delArray (tmp);
}
template void wfwfs::smooth (float**, size_t, size_t, size_t, size_t);
template void wfwfs::smooth (double**, size_t, size_t, size_t, size_t);
template void wfwfs::smooth (complex_t**, size_t, size_t, size_t, size_t);


template <typename T>
void wfwfs::ScharmerFilter (T** data, double** q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff) {

    static const double hi = 1.0;
    static const double lo = 0.2;
    double noiseFactor = noise_power*sizeY*sizeX;
    
    double** pt = newArray<double> (sizeY, sizeX);
    double** rr = newArray<double> (sizeY, sizeX);
    
    std::transform (*data, *data + sizeY * sizeX, *q2_inv, *pt, [](const T& a, const double& b ){ return norm(a)*b; } );
    
    smooth (pt, sizeY, sizeX, 1, 1);
    
    std::transform (*pt, *pt + sizeY * sizeX, *pt, [noiseFactor](const double& a){ return noiseFactor/a;} );
    
    for (unsigned int y=0; y<sizeY; ++y) {
        rr[y][0] = pt[y][0];
        for (unsigned int x=1; x<sizeX; ++x) {
            rr[y][x] = pt[y][sizeX-x];
        }
    }
    
    double* tmp = newArray<double> (sizeY);
    for(unsigned int x=0; x<sizeX; ++x){                        // rr==filter
        rr[0][x] = max((1.0-0.5*(rr[0][x]+pt[0][x])),0.0);
        for(unsigned int y=1; y<sizeY;++y) tmp[y]=rr[y][x];       // temporary storage
        for(unsigned int y=1; y<sizeY;++y) rr[y][x]=max((1.0-0.5*(tmp[sizeY-y]+pt[y][x])), 0.0);
    }
    delArray (tmp);
    
    for(unsigned int y=0; y<sizeY; ++y) {
        for(unsigned int x=0; x<sizeX; ++x){
            if(rr[y][x]<lo) rr[y][x] = 0.0;
            if(rr[y][x]>hi) rr[y][x] = 1.0;
        }
    }
    
    unsigned int yHalf = sizeY/2;
    unsigned int xHalf = sizeX/2;
    rr[yHalf][xHalf] = 1.0;                   // ; DC gain = 1
    
    connectedRegion(rr, sizeY, sizeX, yHalf, xHalf, 0.0);
    smooth(rr, sizeY, sizeX, 4, 4);
    
    frequency_cutoff *= frequency_cutoff;
    for(unsigned int y=0; y<sizeY; ++y) {
        for(unsigned int x=0; x<sizeX; ++x) {
            if( (sqr(y-yHalf)+sqr(x-xHalf)) > frequency_cutoff){
                data[y][x] = 0.0;
            } else{
                data[y][x] *= rr[y][x];
            }
        }
    }
    delArray (pt);
    delArray (rr);
    
}
template void wfwfs::ScharmerFilter (float**, double**, size_t, size_t, double, double);
template void wfwfs::ScharmerFilter (double**, double**, size_t, size_t, double, double);
template void wfwfs::ScharmerFilter (complex_t**, double**, size_t, size_t, double, double);


double wfwfs::inv_dist_wght (float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    int xl = std::max (0L, static_cast<int64_t> (posX - maxDistance));
    int xh = std::min (sizeX, posX + maxDistance + 1);
    int yl = std::max (0L, static_cast<int64_t> (posY - maxDistance));
    int yh = std::min (sizeY, posY + maxDistance + 1);


    double weight = 0.0, res = 0.0;
    for (int y = yl; y < yh; ++y)
        for (int x = xl; x < xh; ++x)
            if (a[y][x]) {
                double c = pow ( (double) sqr (x - posX) + (double) sqr (y - posY) + deltasqr, -beta);
                res += c * a[y][x];
                weight += c;
            }
    return res / weight;
}


template <typename T>
double wfwfs::inverseDistanceWeight (T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    static const double* const inverseDistanceSquared = getDistanceMap();     // will only be initialized once.

    // TODO: verify this function, results look weird
    int64_t beginX = std::max (0L, static_cast<int64_t> (posX - maxDistance));
    int64_t endX = std::min (sizeX, posX + maxDistance+1);
    int64_t beginY = std::max (0L, static_cast<int64_t> (posY - maxDistance));
    int64_t endY = std::min (sizeY, posY + maxDistance+1);

    double normalization = 0.0, weightedSum = 0.0;
    for (int y=beginY; y < endY; ++y) {
        int y2 = (y-posY)*(y-posY);
        for (int x = beginX; x < endX; ++x) {
            int x2 = (x-posX)*(x-posX);
            //if( x2+y2 > maxDistance2 ) break;
            if( array[y][x] ) {
                double tmp = inverseDistanceSquared[y2+x2];
                weightedSum += tmp * array[y][x];
                normalization += tmp;
            }
        }
    }
    if( normalization ) {
        return weightedSum / normalization;
    }
    return 0.0;
}
template double wfwfs::inverseDistanceWeight (unsigned char**, size_t, size_t, size_t, size_t);
template double wfwfs::inverseDistanceWeight (short**, size_t, size_t, size_t, size_t);
template double wfwfs::inverseDistanceWeight (int**, size_t, size_t, size_t, size_t);
template double wfwfs::inverseDistanceWeight (float**, size_t, size_t, size_t, size_t);
template double wfwfs::inverseDistanceWeight (double**, size_t, size_t, size_t, size_t);

template <typename T>
double wfwfs::horizontalInterpolation (T** array, size_t sizeY, size_t sizeX, size_t posY, size_t posX) {

    T* ptr = array[posY];

    //map the 5 pixel surrounding as bits in a byte
    int val = 0;
    if ( (posX > 1)) val |= ( (ptr[posX - 2] > 0) << 4);
    if ( (posX > 0))   val |= ( (ptr[posX - 1] > 0) << 3);
    if ( (posX+1 < sizeX))   val |= ( (ptr[posX + 1] > 0) << 1);
    if ( (posX+2 < sizeX)) val |= ( (ptr[posX + 2] > 0));
    //now select based on the number
    switch (val) {
        case (10) :     // = 0 1 x 1 0
        case (11) :     // = 0 1 x 1 1
        case (26) :     // = 1 1 x 1 0
        case (27) :     // = 1 1 x 1 1
            return (ptr[posX-1] + ptr[posX+1]) / 2;
        case (18) :     // = 1 0 x 1 0
        case (19) :     // = 1 0 x 1 1
            return (ptr[posX-2] + 2 * ptr[posX+1]) / 3;
        case (9) :      // = 0 1 x 0 1
        case (25) :     // = 1 1 x 0 1
            return (2 * ptr[posX-1] + ptr[posX+2]) / 3;
        default:
            return inverseDistanceWeight<T> (array, sizeY, sizeX, posY, posX);
    }

}
template double wfwfs::horizontalInterpolation (float**, size_t, size_t, size_t, size_t);
template double wfwfs::horizontalInterpolation (double**, size_t, size_t, size_t, size_t);


template <typename T>
void wfwfs::apodizeInPlace( T** data, size_t nRows, size_t nCols, size_t rowBlend, size_t colBlend, size_t rowMargin, size_t colMargin ) {

    if( rowBlend+rowMargin > nCols/2 ) {
        rowBlend = nCols/2-rowMargin;
    }
    if( colBlend+colMargin > nRows/2 ) {
        colBlend = nRows/2-colMargin;
    }
    
    size_t sz = std::max( rowBlend+rowMargin, colBlend+colMargin ) + 2;
    double* tmp = new double[ sz ];
    
    memset( tmp, 0, sz*sizeof(T) );
    wfwfs::apodize( tmp+rowMargin, rowBlend+2, 1.0 );
    for( size_t c=0; c<nCols; ++c ) {
        for( size_t r=0; r<rowBlend+rowMargin; ++r ) {
            data[r][c] *= tmp[r+1];
            data[nRows-r-1][c] *= tmp[r+1];
        }
    }
    memset( tmp, 0, sz*sizeof(T) );
    wfwfs::apodize( tmp+colMargin, colBlend+2, 1.0 );
    for( size_t rr=0; rr<nRows; ++rr ) {
        for( size_t c=0; c<colBlend+colMargin; ++c ) {
            data[rr][c] *= tmp[c+1];
            data[rr][nCols-c-1] *= tmp[c+1];
        }
    }
    delete[] tmp;

}
template void wfwfs::apodizeInPlace(int16_t**, size_t, size_t, size_t, size_t, size_t, size_t);
template void wfwfs::apodizeInPlace(int32_t**, size_t, size_t, size_t, size_t, size_t, size_t);
template void wfwfs::apodizeInPlace(float**, size_t, size_t, size_t, size_t, size_t, size_t);
template void wfwfs::apodizeInPlace(double**, size_t, size_t, size_t, size_t, size_t, size_t);
template void wfwfs::apodizeInPlace(complex_t**, size_t, size_t, size_t, size_t, size_t, size_t);


void wfwfs::make_mask( cv::InputArray image, cv::InputOutputArray mask, double thres, int smooth, bool filterLarger, bool invert) {
    
    using namespace cv;
    
    Mat src = image.getMat();
    Mat dst = mask.getMat();
    
    CV_Assert( !src.empty() );
    CV_Assert( !dst.empty() );

    Mat graySrc;
    src.convertTo( graySrc, CV_32FC1 );
    
    threshold( graySrc, graySrc, thres, 1, THRESH_BINARY ); // THRESH_BINARY,THRESH_TOZERO
    
    double minValue, maxValue;
    minMaxLoc( graySrc, &minValue, &maxValue );
    if( maxValue == minValue ) {        // uniform input, no need to filter.
        if( invert ) {
            graySrc = 1 - graySrc;
        }
        graySrc.assignTo( dst, dst.type() );
        return;
    }

    if( smooth < 2 && !invert ) {         // if no filtering was requested, return 
        graySrc.assignTo( dst, dst.type() );
        return;
    }    
    
    Mat filteredMask( graySrc.rows, graySrc.cols, CV_8UC1 );
    graySrc.convertTo( filteredMask, CV_8UC1, 1 );
    
    if( smooth > 1 ) {
        Mat originalMask = filteredMask.clone();
        // shape:  enum { MORPH_RECT, MORPH_CROSS, MORPH_ELLIPSE };
        // operations: enum { MORPH_ERODE, MORPH_DILATE, MORPH_OPEN, MORPH_CLOSE, MORPH_GRADIENT, MORPH_TOPHAT, MORPH_BLACKHAT };
        cv::Point anchor( smooth/2, smooth/2 );
        Mat element = getStructuringElement( MORPH_RECT, Size( smooth, smooth ), anchor );
        morphologyEx( filteredMask, filteredMask, MORPH_CLOSE, element, anchor, 1, BORDER_REFLECT_101 );
        if( !(smooth & 1) ) anchor -= cv::Point( 1, 1 );    // move anchor-point for the reverse transform, so that asymmetry-shifts cancel
        morphologyEx( filteredMask, filteredMask, MORPH_OPEN, element, anchor, 1, BORDER_REFLECT_101 ); 
        bitwise_xor( filteredMask, originalMask, filteredMask );
        if( filterLarger ) {
            bitwise_or( 1-filteredMask, originalMask, filteredMask );
        } else {
            bitwise_xor( filteredMask, originalMask, filteredMask );
        }
    }
    
    if( invert ) {
        filteredMask = 1 - filteredMask;
    }
    
    filteredMask.assignTo( dst, dst.type() );

}


void wfwfs::morph_smooth( cv::InputOutputArray data, int smooth, bool filterLarger ) {
    
    using namespace cv;
    
    Mat dataMat = data.getMat();
    CV_Assert( !dataMat.empty() );
    
    Mat tmpImg( dataMat.rows, dataMat.cols, CV_8UC1 );
    dataMat.convertTo( tmpImg, CV_8UC1, 1 );
    
    if( smooth > 1 ) {
        Mat originalMask = tmpImg.clone();
//         if( invert ) {
//             cv::bitwise_not( tmpImg, tmpImg );  
//         }
        // shape:  enum { MORPH_RECT, MORPH_CROSS, MORPH_ELLIPSE };
        // operations: enum { MORPH_ERODE, MORPH_DILATE, MORPH_OPEN, MORPH_CLOSE, MORPH_GRADIENT, MORPH_TOPHAT, MORPH_BLACKHAT };
        cv::Point anchor( smooth/2, smooth/2 );
        Mat element = getStructuringElement( MORPH_RECT, Size( smooth, smooth ), anchor );
        morphologyEx( tmpImg, tmpImg, MORPH_CLOSE, element, anchor, 1, BORDER_REFLECT_101 );
        if( !(smooth & 1) ) anchor -= cv::Point( 1, 1 );    // move anchor-point for the reverse transform, so that asymmetry-shifts cancel
        morphologyEx( tmpImg, tmpImg, MORPH_OPEN, element, anchor, 1, BORDER_REFLECT_101 ); 
        bitwise_xor( tmpImg, originalMask, tmpImg );
        if( filterLarger ) {
            bitwise_or( 1-tmpImg, originalMask, tmpImg );
        } else {
            bitwise_xor( tmpImg, originalMask, tmpImg );
        }
    }
    
    tmpImg.assignTo( dataMat, dataMat.type() );

}


namespace {

    PointF lineIntersection( const cv::Vec4i &line1, const cv::Vec4i &line2 ) {
        
        PointF intersection(NAN);
        
        double A1 = (line1[2] - line1[0]);
        double B1 = (line1[3] - line1[1]);
        double C1 = (A1*line1[0] + B1*line1[1]);

        double A2 = (line2[2] - line2[0]);
        double B2 = (line2[3] - line2[1]);
        double C2 = (A2*line2[0] + B2*line2[1]);

        double det = (A1*B2)-(A2*B1);
        if( fabs(det) > 1E-5 ) {
            intersection.x = static_cast<float>(((C1*B2)-(C2*B1)) / (det));
            intersection.y = static_cast<float>(((C2*A1)-(C1*A2)) / (det));
        }

        return intersection;
    }

    
}


std::vector<PointF> wfwfs::detect_corners( cv::InputArray data, double rho, double theta, int threshold, double minLineLength, double maxLineGap, int smooth ) {
    
    std::vector<PointF> ret;
    
    using namespace cv;
    
    Mat dataMat = data.getMat();
    CV_Assert( !dataMat.empty() );
    
    Mat tmpImg;
    Canny( dataMat, tmpImg, 1, 1, 3 );
    
    if( smooth > 1 ) {
        Mat originalMask = tmpImg.clone();
        cv::Point anchor( smooth/2, smooth/2 );
        Mat element = getStructuringElement( MORPH_RECT, Size( smooth, smooth ), anchor );
        morphologyEx( tmpImg, tmpImg, MORPH_CLOSE, element, anchor, 1, BORDER_REFLECT_101 );
        if( !(smooth & 1) ) anchor -= cv::Point( 1, 1 );    // move anchor-point for the reverse transform, so that asymmetry-shifts cancel
        morphologyEx( tmpImg, tmpImg, MORPH_OPEN, element, anchor, 1, BORDER_REFLECT_101 ); 
    }
    
    vector<Vec4i> lines;  
    HoughLinesP( tmpImg, lines, rho, theta, threshold, minLineLength, maxLineGap );
    for( size_t i(0); i<lines.size(); ++i ) {
        for( size_t j(i+1); j<lines.size(); ++j ) {
            wfwfs::PointF intersection = lineIntersection( lines[i], lines[j] );
            if( intersection.isfinite() && (intersection.min() > 0) && (intersection.x < dataMat.rows) && (intersection.y < dataMat.cols) ) {
                ret.push_back( intersection );
            }
        }
    }
    
    return ret;
    
}


void wfwfs::bounding_rect( cv::InputOutputArray data, PointI& first, PointI& last ) {
    
    using namespace cv;
    
    Mat dataMat = data.getMat();
    CV_Assert( !dataMat.empty() );
    
    Mat tmpImg;
    Canny( dataMat, tmpImg, 0, 1, 5 );
    
    vector< vector<cv::Point> > contours;
    findContours( tmpImg, contours, RETR_EXTERNAL, CHAIN_APPROX_NONE );     // RETR_LIST, CHAIN_APPROX_SIMPLE
    
    if( contours.empty() ) return;
    
    if( contours.size() == 1 ) {
        auto c = contours.front();
        Rect r = boundingRect( c );
        if( (r.height > 0) && (r.width > 0) ) {
            first = PointI( r.y, r.x );
            last = first + PointI( r.height-1, r.width-1 );
        }
    }

}


cv::Mat wfwfs::getFloatMat( cv::InputOutputArray data ) {

    using namespace cv;
    
    Mat ret;
    try {
        Mat dataMat = data.getMat();
        CV_Assert( !dataMat.empty() );
        dataMat.convertTo( ret, CV_32FC1 );
        double minValue, maxValue;
        cv::minMaxLoc( ret, &minValue, &maxValue);
        ret = (ret - minValue) / (maxValue - minValue);
    } catch( cv::Exception& e ) {
        std::cerr << "OpenCV error: " << e.msg << std::endl;
    }

    return ret;

}

template <typename T>
double wfwfs::flat2gain( const T* in, T* out, size_t sizeY, size_t sizeX, size_t stride, const uint8_t* mask2, int smooth ) {
    
    if( smooth == 0 ) smooth = 7;
    
    float bad = 1.0;
    float mx = 4.0;
    float mn = 0.1;
    
    ArrayStats stats;
    stats.getMedian( in, sizeY, sizeX, stride, mask2 );
    
    size_t nEl = sizeY*sizeX;
    Array<T> tmp( 2*sizeY, sizeX );
    T* tmpPtr1 = tmp.get();
    T* tmpPtr2 = tmpPtr1+nEl;

    Array<uint8_t> tmpMask( sizeY, sizeX );
    uint8_t* maskPtr = tmpMask.get();
    
    for( size_t yy(0); yy<sizeY; ++yy ) {
        for( size_t xx(0); xx<sizeX; ++xx ) {
            size_t offset1 = yy*stride+xx;
            size_t offset2 = yy*sizeX+xx;
            if( in[offset1] < (stats.median*1E-5) ) tmpPtr1[offset2] = stats.median*1E-5;
            tmpPtr1[offset2] = stats.median/in[offset1];
            if( !isfinite(tmpPtr1[offset2]) ) tmpPtr1[offset2] = 0;
        }
    }
    

    cv::Mat tmpMat1( sizeY, sizeX, cvType<T>(), tmpPtr1 );
    cv::Mat tmpMat2( sizeY, sizeX, cvType<T>(), tmpPtr2 );
    cv::Mat maskMat( sizeY, sizeX, cvType<uint8_t>(), maskPtr );
    double sigma = smooth/( 4 * sqrt(log(2.0)) );
    cv::GaussianBlur( tmpMat1, tmpMat2, cv::Size(), sigma, sigma, cv::BORDER_REFLECT_101 );
    
    for( size_t i(0); i<nEl; ++i ) {
        maskPtr[i] = (tmpPtr1[i] >= mn) && (tmpPtr1[i] <= mx) &&
            ((tmpPtr1[i]-tmpPtr2[i]) < bad) && isfinite(tmpPtr1[i]); 
    }

    cv::Point anchor(-1,-1); //( 3, 3 );
    cv::Mat element = cv::getStructuringElement( cv::MORPH_RECT, cv::Size( 5, 5 ), anchor );
    //cv::morphologyEx( maskMat, maskMat, cv::MORPH_CLOSE, element, anchor, 1, cv::BORDER_REFLECT_101 );
    //anchor -= cv::Point( 1, 1 );
    cv::morphologyEx( maskMat, maskMat, cv::MORPH_OPEN, element, anchor, 1, cv::BORDER_REFLECT_101 );
    for( size_t yy(0); yy<sizeY; ++yy ) {
        for( size_t xx(0); xx<sizeX; ++xx ) {
            size_t offset1 = yy*stride+xx;
            size_t offset2 = yy*sizeX+xx;
            if( maskPtr[offset2] && isfinite(tmpPtr1[offset2]) ) {
                out[offset1] = stats.median/in[offset1];
                if( !isfinite(out[offset1]) ) out[offset1] = 0;
            } else {
                out[offset1] = 0;
            }
        }
    }
     
    return stats.median;
}
template double wfwfs::flat2gain( const float*, float*, size_t, size_t, size_t, const uint8_t*, int );
template double wfwfs::flat2gain( const double*, double*, size_t, size_t, size_t, const uint8_t*, int );


