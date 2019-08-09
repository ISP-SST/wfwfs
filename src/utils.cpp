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

#include <functional>
#include <map>
#include <set>
#include <math.h>

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
    void maskConnected(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int yy, unsigned int xx, unsigned int y_start, unsigned int x_start, T threshold) {
        if (yy >= sizeY || xx >= sizeX) return;
        if ( !mask[yy][xx] && ((yy == y_start && xx == x_start) || data[yy][xx] > threshold)) {       // new point inserted, check neighbours
            mask[yy][xx] = 1;
            maskConnected(data, mask, sizeY, sizeX, yy + 1, xx, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy - 1, xx, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy, xx + 1, y_start, x_start, threshold);
            maskConnected(data, mask, sizeY, sizeX, yy, xx - 1, y_start, x_start, threshold);
        }
    }

}

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
void wfwfs::connectedRegion(T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int y_start, unsigned int x_start, T threshold) {
    memset(*mask, 0, sizeY * sizeX);
    maskConnected(data, mask, sizeY, sizeX, y_start, x_start, y_start, x_start, threshold);
}
template void wfwfs::connectedRegion(float**, uint8_t**, size_t, size_t, uint, uint, float);
template void wfwfs::connectedRegion(double**, uint8_t**, size_t, size_t, uint, uint, double);


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


