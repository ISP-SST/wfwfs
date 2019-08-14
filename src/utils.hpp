#ifndef WFWFS_UTIL_HPP
#define WFWFS_UTIL_HPP
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

#include "array.hpp"
#include "point.hpp"

#include <thread>

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/calib3d/calib3d.hpp>

namespace wfwfs {

        template <typename T> void apodize( T* data, size_t n, const T& target );


        template <typename T>
        std::vector<double> fitPlane( const T* in, size_t sizeY, size_t sizeX, size_t strideY=0 );

        template <typename T>
        wfwfs::Array<T> fitPlane( const wfwfs::Array<T>& in, bool subtract_mean=false, double* coeffs=nullptr );


        template <typename T>
        T median( std::vector<T> in ) {
            std::nth_element( in.begin(), in.begin() + in.size() / 2, in.end() );
            return *( in.begin() + in.size() / 2 );
        }


        template <typename T>
        T median( const wfwfs::Array<T>& in ) {
            wfwfs::Array<T> tmp = in.copy();  // nth_element is destructive, so make a deep copy.
            std::nth_element( tmp.begin(), tmp.mid(), tmp.end() );
            return *tmp.mid();
        }

        template <typename T>
        void maskConnected2( const T* data, uint8_t* mask, size_t sizeY, size_t sizeX, unsigned int startY, unsigned int startX, T threshold=T(0));

        template <typename T>
        void connectedRegion(const T** data, uint8_t** mask, size_t sizeY, size_t sizeX, unsigned int y, unsigned int x, T threshold=T(0));
        template <typename T>
        void connectedRegion(T** data, size_t sizeY, size_t sizeX, unsigned int y, unsigned int x, T threshold=T(0)) {
            uint8_t** mask = wfwfs::newArray<uint8_t>(sizeY,sizeX);
            connectedRegion(const_cast<const T**>(data), mask, sizeY, sizeX, y, x, threshold);
            std::transform(*data, *data+sizeY*sizeX, *mask, *data, wfwfs::multiply<T,uint8_t>());
            wfwfs::delArray(mask);
        }
        template <typename T>
        void connectedRegion(T* data, size_t sizeY, size_t sizeX, unsigned int y, unsigned int x, T threshold=T(0)) {
            T** ptr = wfwfs::makePointers(data, sizeY, sizeX);
            connectedRegion(ptr, sizeY, sizeX, y, x, threshold);
            wfwfs::delPointers(ptr);
        }
       
        template <typename T>
        void smooth(T** data, size_t sizeY, size_t sizeX, size_t nY, size_t nX);
        template <typename T>
        void smooth(T** data, size_t sizeY, size_t sizeX, size_t n) { smooth(data, sizeY, sizeX, n, n); }
        template <typename T>
        void smooth(T* data, size_t sizeY, size_t sizeX, size_t nY, size_t nX) {
            T** ptr = wfwfs::makePointers(data, sizeY, sizeX);
            smooth(ptr, sizeY, sizeX, nY, nX);
            wfwfs::delPointers(ptr);
        }
        template <typename T>
        void smooth(T* data, size_t sizeY, size_t sizeX, size_t n) { smooth(data, sizeY, sizeX, n, n); }
        
        template <typename T>
        void ScharmerFilter (T** data, double** q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff);
        template <typename T>
        void ScharmerFilter(T* data, double* q2_inv, size_t sizeY, size_t sizeX, double noise_power, double frequency_cutoff) {
            T** ptr = wfwfs::makePointers(data, sizeY, sizeX);
            double** qPtr = wfwfs::makePointers(q2_inv, sizeY, sizeX);
            ScharmerFilter(ptr, qPtr, sizeY, sizeX, noise_power, frequency_cutoff);
            wfwfs::delPointers(ptr);
            wfwfs::delPointers(qPtr);
        }
        
        template <typename T>
        double inverseDistanceWeight( T**, size_t sizeY, size_t sizeX, size_t posY, size_t posX );

        double inv_dist_wght( float **a, size_t sizeY, size_t sizeX, size_t posY, size_t posX );

        template <typename T>
        double horizontalInterpolation( T**, size_t sizeY, size_t sizeX, size_t posY, size_t posX );


        template <typename T, typename U>
        double chisq( const wfwfs::Array<T>& a, const wfwfs::Array<U>& b ) {
            if( !a.sameSizes( b ) ) {
                throw std::logic_error( "Array dimensions does not match." );
            }
            size_t nElements = a.nElements();
            if( nElements == 0 ) return 0.0;
            double chisq = 0;
            typename wfwfs::Array<U>::const_iterator bit = b.begin();
            for( auto &avalue : a ) {
                double tmp = avalue - *bit++;
                chisq += tmp * tmp;
            }
            return chisq / static_cast<double>( nElements );
        }


        template <typename T, typename U, typename V>
        double chisq( const wfwfs::Array<T>& a, const wfwfs::Array<U>& b, const wfwfs::Array<V>& weight ) {
            if( !a.sameSizes( b ) || !a.sameSizes( weight ) ) {
                throw std::logic_error( "Array dimensions does not match." );
            }
            double tmp, chisq = 0;
            typename wfwfs::Array<U>::const_iterator bit = b.begin();
            typename wfwfs::Array<V>::const_iterator wit = weight.begin()--;
            size_t count(0);
            for( auto &avalue : a ) {
                if( *++wit ) {
                    tmp = ( avalue - *bit++ ) * ( *wit );
                    chisq += tmp * tmp;
                    ++count;
                }
            }
            if( count == 0 ) return 0.0;
            return chisq / static_cast<double>( count );
        }

        template <typename T>
        void clipImage( wfwfs::Array<T>& img, const std::vector<int16_t> clip, bool symmetricClip=false )  {
            size_t nDims = img.nDimensions();
            if( nDims < 2 ) return;
            std::vector<size_t> finalClip;
            size_t nImgs = 1;
            if( nDims > 2 ) {
                for( size_t i=0; i<nDims-2; ++i ) {
                    size_t sz = img.dimSize(i);
                    nImgs *= sz;
                    finalClip.push_back(0);
                    finalClip.push_back(sz-1);
                }
            }
            std::vector<int16_t> tmpClip = clip;
            if( tmpClip.size() == 2 ) {        // apply same values to both dimensions.
                tmpClip.insert( tmpClip.end(), clip.begin(), clip.end() );
            }
            
            if( tmpClip.size() == 4 ) {
                bool flipX = false, flipY = false;
                // we have the y (row/slow) dimension first, momfbd cfg-files (and thus alignClip) has x first.
                if ( tmpClip[0] > tmpClip[1] ) {
                    std::swap( tmpClip[0], tmpClip[1] );
                    flipX = true;
                }
                if ( tmpClip[2] > tmpClip[3] ) {
                    std::swap( tmpClip[2], tmpClip[3] );
                    flipY = true;
                }
                for( auto & index : tmpClip )
                    --index;       // NOTE: momfbd cfg files uses 1-based indexes, internally we start with 0.
                size_t sy = tmpClip[3] - tmpClip[2] + 1;
                size_t sx = tmpClip[1] - tmpClip[0] + 1;
                if( symmetricClip ) {
                    const std::vector<size_t>& dims = img.dimensions();
                    int skewY = (dims[0] - sy) / 2  - tmpClip[2];
                    int skewX = (dims[1] - sx) / 2  - tmpClip[0];
                    tmpClip[0] += skewX;
                    tmpClip[1] += skewX;
                    tmpClip[2] += skewY;
                    tmpClip[3] += skewY;
                }
                
                finalClip.push_back(tmpClip[2]);
                finalClip.push_back(tmpClip[3]);
                finalClip.push_back(tmpClip[0]);
                finalClip.push_back(tmpClip[1]);
                img.setLimits( finalClip );
                img.trim();

                if( flipX || flipY ) {
                    std::shared_ptr<T*> arrayPtr = img.reshape(nImgs*sy, sx);
                    T** imgPtr = arrayPtr.get();
                    for( size_t i=0; i<nImgs; ++i) {
                        if (flipX) wfwfs::reverseX(imgPtr, sy, sx);
                        if (flipY) wfwfs::reverseY(imgPtr, sy, sx);
                        imgPtr += sy;
                    }
                }
            }
        }
        
        
        template <typename T, typename Predicate>
        void fillPixels( wfwfs::Array<T>& array, T fillValue, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            for( auto & value : array ) {
                if( predicate( value ) ) value = fillValue;
            }
        }


        template <typename T, typename Predicate, typename MaskType=uint8_t>
        void fillPixels( T** array, size_t sy, size_t sx, std::function<double( size_t, size_t )> filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ), MaskType** mask=nullptr  ) {
            std::map<size_t, T> tmp;
            T* ptr = *array;
            size_t offset = 0;
            for( size_t y = 0; y < sy; ++y ) {
                for( size_t x = 0; x < sx; ++x ) {
                    if( (!mask || !mask[y][x]) && predicate( array[y][x] ) ) tmp.insert( std::pair<size_t, T>( offset, filler( y, x ) ) );
                    ++offset;
                }
            }
            size_t cnt = 0;
            for( auto &it : tmp ) {
                ptr[it.first] = it.second;
                ++cnt;
            }
        }


        template <typename T, typename FillFunction, typename Predicate>
        void fillPixels( wfwfs::Array<T>& array, FillFunction* filler, Predicate predicate = std::bind2nd( std::less_equal<T>(), 0 ) ) {
            for( auto it = array.begin(); it != array.end(); ++it ) {
                if( predicate( *it ) ) *it = filler( it );
            }
        }

        template <typename T, typename U=uint8_t>
        void fillPixels( T** image, size_t sy, size_t sx, U** mask=nullptr ) {
            
            size_t offset(0);
            const size_t nPixels = sy*sx;

            std::map<size_t, T> values;
            size_t o;
            while( (o = offset++) < nPixels ) {
                if( !mask || mask[0][o%nPixels] ) {
                    size_t y = o/sy;
                    size_t x = o%sy;
                    values.insert( std::pair<size_t, T>( o, inverseDistanceWeight( image, sy, sx, y, x ) ) );
                    //values.insert( std::pair<size_t, T>( o, horizontalInterpolation( image, sy, sx, y, x ) ) );
                }
            }
            for( auto &it: values ) {
                image[0][it.first] = it.second;
            }

        }

        
        template <typename T, typename U=uint8_t>
        void fillPixels( T*** image, size_t nImages, size_t sy, size_t sx, U** mask=nullptr, unsigned int nThreads=std::thread::hardware_concurrency() ) {
            
            std::atomic<size_t> offset(0);
            const size_t nPixels = sy*sx;
            
            auto nextbad = [&](void) {
                size_t o;
                while( (o = offset.fetch_add(1)) < nPixels ) {
                    if( !mask || mask[0][o%nPixels]) return o;
                }
                return nPixels;
            };

            std::vector<std::thread> threads;
            for( unsigned int t=0; t<nThreads; ++t ) {
                threads.push_back( std::thread(
                    [&](){
                        size_t myOffset;
                        std::map<size_t, T> values;
                        while( (myOffset=nextbad()) < nPixels ) {
                            size_t y = myOffset/sy;
                            size_t x = myOffset%sy;
                            values.insert( std::pair<size_t, T>( myOffset, inverseDistanceWeight( image, sy, sx, y, x ) ) );
                        }

                        for( auto &it: values ) {
                            image[0][it.first] = it.second;
                        }
                    }));
            }
            for( auto& th : threads ) th.join();

        }
                
        template <typename T>
        void apodizeInPlace( T** data, size_t nRows, size_t nCols, size_t rowBlend, size_t colBlend, size_t rowMargin=0, size_t colMargin=0 );
        template <typename T>
        void apodizeInPlace( wfwfs::Array<T>& array, size_t blendRegion, size_t margin=0 ) {
            size_t nRows = array.dimSize(0);
            size_t nCols = array.dimSize(1);
            std::shared_ptr<T*> tmp = array.reshape( nRows, nCols );
            apodizeInPlace( tmp.get(), nRows, nCols, blendRegion, blendRegion, margin, margin );
        }
        template <typename T>
        void apodizeInPlace( wfwfs::Array<T>& array, size_t rowBlend, size_t colBlend, size_t rowMargin, size_t colMargin ) {
            size_t nRows = array.dimSize(0);
            size_t nCols = array.dimSize(1);
            std::shared_ptr<T*> tmp = array.reshape( nRows, nCols );
            apodizeInPlace( tmp.get(), nRows, nCols, rowBlend, colBlend, rowMargin, colMargin );
        }

        template <typename T>
        wfwfs::Array<T> apodize( const wfwfs::Array<T>& in, size_t blendRegion, size_t margin=0 ) {
            if( !blendRegion && !margin ) return in;       // nothing to do
            wfwfs::Array<T> array;
            in.copy( array );
            apodizeInPlace( array, blendRegion, margin );
            return std::move(array);
        }
        template <typename T>
        wfwfs::Array<T> apodize( const wfwfs::Array<T>& in, size_t rowBlend, size_t colBlend, size_t rowMargin, size_t colMargin ) {
            if( !rowBlend && !colBlend && !rowMargin && !colMargin ) return in;       // nothing to do
            wfwfs::Array<T> array;
            in.copy( array );
            apodizeInPlace( array, rowBlend, colBlend, rowMargin, colMargin );
            return std::move(array);
        }

        template <typename T>
        void img_trim( T**& img, size_t& imgRows, size_t& imgCols, float threshold=1E-6 ) {
            std::vector<T> rowSums( imgCols, 0 );
            std::vector<T> colSums( imgRows, 0 );
            for ( size_t c=0; c < imgCols; ++c ) {
                for ( size_t r=0; r < imgRows; ++r ) {
                    T val = img[r][c];
                    rowSums[c] += val;
                    colSums[r] += val;
                }
            }

            size_t firstCol(0), lastCol(imgCols-1), firstRow(0), lastRow(imgRows-1);
            while ( firstCol < lastCol && rowSums[firstCol] < threshold ) ++firstCol;
            while ( lastCol && rowSums[lastCol] < threshold  ) --lastCol;
            while ( firstRow < lastRow && colSums[firstRow] < threshold  ) ++firstRow;
            while ( lastRow && colSums[lastRow] < threshold  ) --lastRow;

            size_t imgRows2 = lastRow-firstRow+1;
            size_t imgCols2 = lastCol-firstCol+1;

            if( imgRows2 != imgRows || imgCols2 != imgCols ) {
                T** tmp = wfwfs::newArray<T>( imgRows2, imgCols2 );
                for( size_t r=0; r<imgRows2; ++r ) {
                    memcpy( tmp[r], img[r+firstRow]+firstCol, imgCols2*sizeof(T) );
                }
                wfwfs::delArray( img );
                img = tmp;
                imgRows = imgRows2;
                imgCols = imgCols2;
            }

        }
        

        template <typename T>
        void draw_line( Array<T>& img, PointI pos1, PointI pos2, const T value=T(), size_t width=1 ) {

            if( width == 0 ) return;
            
            PointD diff = pos2 - pos1;
            size_t length = sqrt(diff.norm());
            double nrm = 1.0/length;
            size_t stride = img.dimSize(1);
            T* ptr = img.get();
            PointI delta(1,0);
            if( fabs(diff.x) < fabs(diff.y) ) {
                delta = PointI(0,1);
            }
            for( size_t w(0); w<width; ++w ) {
                PointI pos = pos1;
                if( w ) {
                    pos += delta*((w+1)/2)*pow(-1.0,w+1);
                    
                }
                for( size_t i(0); i<=length; ++i ) {
                    PointI tmp = pos + diff*nrm*i;
                    size_t offset = static_cast<size_t>(tmp.y*stride + tmp.x);
                    if( offset>img.nElements() ) continue;
                    ptr[offset] = value;
                }
            }
        }

        void make_mask( cv::InputArray image, cv::InputOutputArray mask, double thres=0, int smooth=5, bool filterLarger=false, bool invert=false );
        void morph_smooth( cv::InputOutputArray mask, int smooth=5, bool filterLarger=false );
        std::vector<PointF> detect_corners( cv::InputArray image, double rho=1, double theta=CV_PI/180, int threshold=1, double minLineLength=100, double maxLineGap=10, int smooth=0 );
        void bounding_rect( cv::InputOutputArray image, PointI& first, PointI& last );
        cv::Mat getFloatMat( cv::InputOutputArray data );

        template <typename T> inline int cvType(void) { return CV_8UC1; }
        template<> inline int cvType<int16_t>(void) { return CV_16SC1; }
        template<> inline int cvType<int32_t>(void) { return CV_32SC1; }
        template<> inline int cvType<float>(void) { return CV_32FC1; }
        template<> inline int cvType<double>(void) { return CV_64FC1; }


        template <typename T, typename U>
        void make_mask( const T* input, U* mask, size_t ySize, size_t xSize, double thres=0, int smooth=5, bool filterLarger=false, bool invert=false ) {
            cv::Mat inMat( ySize, xSize, cvType<T>(), const_cast<T*>(input) );
            cv::Mat maskMat( ySize, xSize, cvType<U>(), mask );
            make_mask( inMat, maskMat, thres, smooth, filterLarger, invert );
        }
        
        template <typename T>
        void morph_smooth( T* data, size_t ySize, size_t xSize, int smooth=5, bool filterLarger=false ) {
            cv::Mat dataMat( ySize, xSize, cvType<T>(), data );
            morph_smooth( dataMat, smooth, filterLarger );
        }
        

        template <typename T>
        std::vector<PointF> detect_corners( T* data, size_t ySize, size_t xSize,
                                            double rho=1, double theta=CV_PI/180, int threshold=1, double minLineLength=100, double maxLineGap=10, int smooth=0 ) {
            cv::Mat dataMat( ySize, xSize, cvType<T>(), data );
            return detect_corners( dataMat, rho, theta,  threshold, minLineLength, maxLineGap, smooth );
        }
        
        template <typename T>
        void bounding_rect( T* data, size_t ySize, size_t xSize, PointI& first, PointI& last ) {
            cv::Mat dataMat( ySize, xSize, cvType<T>(), data );
            bounding_rect( dataMat, first, last );
        }
        
        template <typename T>
        cv::Mat getFloatMat( T* data, size_t ySize, size_t xSize ) {
            cv::Mat dataMat( ySize, xSize, cvType<T>(), data );
            return getFloatMat( dataMat );
        }

        template <typename T, typename U>
        double threshold( T* data, U* out, size_t ySize, size_t xSize, double thresh, double maxval=1, int type=0 ) {
            double ret = 0.0;
            cv::Mat fMat = getFloatMat( data, ySize, xSize );
            cv::Mat outMat( ySize, xSize, cvType<U>(), out );
            cv::Mat bMat( fMat.size(), CV_8UC1 );
            cv::Mat tmpMat( fMat.size(), CV_8UC1 );
            fMat.convertTo( bMat, CV_8UC1, 255 );
            ret = threshold( bMat, tmpMat, thresh, maxval, type );
            // void adaptiveThreshold(InputArray src, OutputArray dst, double maxValue, int adaptiveMethod, int thresholdType, int blockSize, double C)
            // double threshold(InputArray src, OutputArray dst, double thresh, double maxval, int type)
            // cv::threshold( ySize, xSize, cvType<T>(), data );
            tmpMat.convertTo( outMat, cvType<U>() );
            return ret;

        }

        template <typename T, typename U>
        void threshold( T* data, U* out, size_t ySize, size_t xSize, int adaptiveMethod=cv::ADAPTIVE_THRESH_GAUSSIAN_C,
                        int thresholdType=cv::THRESH_BINARY, double maxval=1, int blockSize=5, double offset=0 ) {
            cv::Mat fMat = getFloatMat( data, ySize, xSize );
            cv::Mat outMat( ySize, xSize, cvType<U>(), out );
            cv::Mat bMat( fMat.size(), CV_8UC1 );
            cv::Mat tmpMat( fMat.size(), CV_8UC1 );
            fMat.convertTo( bMat, CV_8UC1, 255 );
            //ret = threshold( bMat, tmpMat, thresh, maxval, type );
            adaptiveThreshold( bMat, tmpMat, maxval, adaptiveMethod, thresholdType, blockSize, offset );
            // double threshold(InputArray src, OutputArray dst, double thresh, double maxval, int type)
            // cv::threshold( ySize, xSize, cvType<T>(), data );
            tmpMat.convertTo( outMat, cvType<U>() );
        }


        template <typename T>
        double flat2gain( const T* in, T* out, size_t sizeY, size_t sizeX, size_t stride, const uint8_t* mask=nullptr, int smooth=0 );
        template <typename T>
        double flat2gain( const T* in, T* out, size_t sizeY, size_t sizeX, const uint8_t* mask=nullptr, int smooth=0 ) {
            return flat2gain( in, out, sizeY, sizeX, sizeX, mask, smooth );
        }
        template <typename T> double flat2gain( const T* in, T* out, size_t count, const uint8_t* mask=nullptr, int smooth=0 ) {
            return flat2gain( in, out, 1, count, count, mask, smooth );
        }
        template <typename T> double flat2gain( const wfwfs::Array<T>& in, wfwfs::Array<T>& out, int smooth=0) {
            if( in.dense() ) return flat2gain( in.ptr(), out.ptr(), in.nElements() );
            else {
                wfwfs::Array<T> tmp;
                in.copy(tmp);
                return flat2gain( tmp.get(), out.ptr(), in.nElements(), nullptr, smooth );
            }
        }

        
}   // wfwfs


#endif  // WFWFS_UTIL_HPP
