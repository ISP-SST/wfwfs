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
#include "seeing.hpp"

#include <boost/algorithm/string.hpp>


using namespace wfwfs;
using namespace std;

namespace bpx = boost::posix_time;


using boost::algorithm::iequals;

float Seeing::pixelsize(8.0E-6);
float Seeing::arcsecs_per_pixel(0.21);
float Seeing::meters_per_pixel(0.000463519);
float Seeing::diam(0.135);
float Seeing::diam_px(291.25);
float Seeing::lambda(500e-9);
double Seeing::radians_per_pixel(0);
double Seeing::dimm_K(0);

        
Seeing::Seeing( void ) : dsID(1) {

    precalculate();
    
}


void Seeing::parsePropertyTree( boost::property_tree::ptree& cfg_ptree ) {
    
    string output_base = cfg_ptree.get<string>( "outputdir", "/data/" );

    
    for( auto& node : cfg_ptree ) {
        if( iequals( node.first, "DIMM" ) ) {
            DimmSet ds( dsID++ );
            ds.parsePropertyTree( node.second );
            dimm_sets.push_back(std::move(ds));
        }
    }

    for( auto& node : cfg_ptree ) { // do the logs separately to ensure all DimmSets are already loaded/created.
        if( iequals( node.first, "LOG" ) ) {
            string outputdir = node.second.get<string>( "outputdir", output_base+"/logs/" );
            string name = node.second.get<string>( "name", "log" );
            int interval = node.second.get<int>( "interval", 1 );
            string columns = node.second.get<string>( "columns", "" );
            if( columns.empty() ) {
                bool first(true);
                for( auto& ds: dimm_sets ) {
                    if( !first ) columns += " ";
                    columns += ds.get_name();
                }
            }
            SeeingLog sl( name, interval );
            sl.setDir( outputdir );
            sl.addColumns( columns );
            logs.push_back( std::move(sl) );
        }
    }

    arcsecs_per_pixel = cfg_ptree.get<float>( "arcsecs_per_pixel", arcsecs_per_pixel );
    diam = cfg_ptree.get<float>( "subaperture_diameter", diam );
    diam_px = cfg_ptree.get<float>( "subaperture_pixels", diam_px );
    lambda = cfg_ptree.get<float>( "wavelength", lambda );
    pixelsize = cfg_ptree.get<float>( "pixelsize", pixelsize );

    precalculate();

    
}


template <typename T>
void Seeing::copy_cell_data( Frame& f, float* dd, float* gg, size_t stride, int maxval ) {

    auto buf = get_buffer();
    T* in = reinterpret_cast<T*>( f.data );
    auto out = buf.get();
    size_t ds_offset(0);

    memset( out, 0, get_data_size() );
        
    for( auto& ds: dimm_sets ) {
        ds.add_frame_data( f.timestamp, buf, ds_offset );
        
        
        auto& cells = ds.get_cells();
        size_t cs = ds.get_cell_size();
        size_t cs2 = cs*cs;
        size_t imgStride = cs2;
        size_t c_offset(0);
        for( size_t n(0); n<cells.size(); ++n ) {
            const Cell& c = cells[n];
            auto outPtr = out + ds_offset + c_offset;
            size_t this_cell_size = cs;
            if( n == ds.get_ref_cell_id() ) {
                this_cell_size = ds.get_ref_cell_size();
                outPtr += ds.get_max_shift()*(cs+1);
            }
            PointI pos = c.pos;
            pos -= this_cell_size/2;
            for( size_t y(0); y < this_cell_size; ++y ) {    // for each cell-row
                size_t offset = (pos.y+y)*stride + pos.x;
                if( dd && gg ) {
                    dg_correct( in+offset, this_cell_size, outPtr, dd+offset, gg+offset, maxval );
                } else {
                    std::copy_n( in+offset, this_cell_size, outPtr );
                }
                outPtr += cs;
            }
            c_offset += imgStride;
        }
        ds_offset += c_offset;
    }

}
template void Seeing::copy_cell_data<uint8_t>( Frame&, float*, float*, size_t, int );
template void Seeing::copy_cell_data<uint16_t>( Frame&, float*, float*, size_t, int );


void Seeing::process( boost::asio::io_service& ios ) {
    
    uint64_t tmp[256*256];
    
    for( auto& ds: dimm_sets ) {
//        ios.post( [&](){
//            uint64_t tmp[256*256];
            ds.measure_shifts( tmp );
            ds.calculate_dimms();
//            ds.calculate_r0();
//        });
    }
    
}


/*
 * These calculations are based on 'The ESO differential image motion monitor'
 * by M. Sarazin and F. Roddier (1989). From its eqs. 13 and 14 we get that:
 * 
 *   r0_l = (2.0*(lambda^2)*(0.179*D^(-1/3)-0.0968*d^(-1/3))/var_l)^(3/5)
 *   r0_t = (2.0*(lambda^2)*(0.179*D^(-1/3)-0.1450*d^(-1/3))/var_t)^(3/5)
 * 
 * where the wavelength and distances are in meters, and the variance in angular
 * units.
 * 
 * By defining
 *   A =  2.0*(lambda^2)*(0.179*D^(-1/3))
 * We get
 *   r0_l = (A*(1-0.54078*(d/D)^(-1/3))/var_l)^(3/5)
 *   r0_t = (A*(1-0.81006*(d/D)^(-1/3))/var_t)^(3/5)
 * We can go a bit further and also define
 *   B = A^(3/5) = ((lambda^2)*(0.358*D^(-1/3)))(3/5)
 *   r0_l = B*((1-0.54078*(d/D)^(-1/3))/var_l)^(3/5)
 *   r0_t = B*((1-0.81006*(d/D)^(-1/3))/var_t)^(3/5)
 * This B remains fixed, and can thus be pre-calculated and only the varying
 * part needs to be computed for each evaluation.
 * For convenience, we also include the factor for converting the variance
 * from pixels to angular units:
 *   K = B / radians_per_pixel^(6/5) = ((lambda/radians_per_pixel)^2*(0.358*D^(-1/3)))(3/5)
 */


void Seeing::precalculate( void ) {
    
    radians_per_pixel = arcsecs_per_pixel / 3600 * M_PI / 180;
    meters_per_pixel = diam / diam_px;
    dimm_K = pow( 0.358*sqr(lambda/radians_per_pixel)*pow(diam, -0.3333), 0.6 );

}


// input is in pixels.
PointD Seeing::apply_dimm_equations( PointD var, float separation ) {
    
    double tmpD = pow( diam_px/separation, 0.3333 );
    PointD ret;

    // longitudinal part
    ret.x = pow( (1 - 0.54078*tmpD)/var.x, 0.6 );
    // transverse part
    ret.y = pow( (1 - 0.81006*tmpD)/var.y, 0.6 );
    
    ret *= dimm_K;
    
    return ret;
    
}


// input is in pixels.
PointD Seeing::simple_dimm_equations( PointD var, float separation ) {

    PointD ret;
    
    var *= radians_per_pixel*radians_per_pixel;     // convert variance from pixel-units to radians
    separation *= diam/diam_px;                     // convert sub-aperture separation from pixels to meters
    
    // longitudinal part
    ret.x = pow(2.0 * sqr(lambda) * (0.179 * pow(diam,-1./3) - 0.0968 * pow(separation,-1./3)) / var.x, 3./5);
    // transverse part
    ret.y = pow(2.0 * sqr(lambda) * (0.179 * pow(diam,-1./3) - 0.1450 * pow(separation,-1./3)) / var.y, 3./5);
    
    return ret;
    
}


void Seeing::start_logs( void ) {
    
    for( auto& l: logs ) {
        l.run( dimm_sets );
    }
    
}


void Seeing::stop_logs( void ) {
    
    for( auto& l: logs ) {
        l.stop();
    }
    
}


string Seeing::adjust_cells( size_t id ) {

    bool changed(false);
    for( auto& ds: dimm_sets ) {
        if( !id || (id == ds.get_id()) ) {
           changed |= ds.adjust_cells();
        }
    }
    if( changed ) {
        return get_cells(id);
    }
    return "";
}


string Seeing::shift_cells( PointI s, size_t id ) {

    for( auto& ds: dimm_sets ) {
        if( !id || (id == ds.get_id()) ) {
           ds.shift_cells(s);
        }
    }

    return get_cells(id);

}


string Seeing::get_cells( size_t id ) const {
    string ret;
    bool first(true);
    for( auto& ds: dimm_sets ) {
        if( !id || (id == ds.get_id()) ) {
            if( !first ) ret += "\n";
            ret += printArray( ds.get_cells(), ds.get_name() );
            first = false;
        }
    }
    return ret;
}


string Seeing::get_shifts( size_t id ) const {
    string ret;
    bool first(true);
    for( auto& ds: dimm_sets ) {
        if( !id || (id == ds.get_id()) ) {
            if( !first ) ret += "\n";
            ret += printArray( ds.get_avg_shifts(), ds.get_name() );
            first = false;
        }
    }
    return ret;
}


template <typename T>
void Seeing::draw_cell( Array<T>& img, PointI pos, uint16_t cell_size, bool mark ) const {

    size_t stride = img.dimSize(1);
    pos -= cell_size/2;
    int offset = pos.y*stride + pos.x;
    if( offset < 0 || (offset+cell_size*stride)>img.nElements() ) return;   // don't draw outside img
    T* ptr = img.get() + offset;
    for( size_t y=0; y<cell_size; ++y ) {
        if( y > 1U && y+2 < cell_size ) {
            ptr[0] = ptr[cell_size-1] = 1;
            ptr[1] = ptr[cell_size-2] = 1;
            if( mark ) {
                ptr[y] = ptr[cell_size-2-y] = 1;
            }
        } else {
            for( size_t x=0; x<cell_size; ++x ) {
                ptr[x] = 1;
            }
        }
        ptr += stride;
    }
    
}


template <typename T>
void Seeing::draw_cells( Array<T>& img ) const {

    img.zero();
    for( const auto& ds: dimm_sets ) {
        size_t ref_id = ds.get_ref_cell_id();
        auto& cells = ds.get_cells();
        for( size_t i(0); i<cells.size(); ++i ) {
            bool is_ref = (i == ref_id);
            int sz = is_ref?ds.get_ref_cell_size():ds.get_cell_size();
            draw_cell( img, cells[i].pos, sz, is_ref );
        }
    }

}
template void Seeing::draw_cells( Array<uint8_t>& ) const;


size_t Seeing::get_data_size( void ) const {
    
    size_t sz(0);
    for( auto& ds: dimm_sets ) {
        sz += ds.get_data_size();
    }
    return sz;
    
}


DimmSet::dimm_data_t Seeing::get_buffer( void ) {

     for( auto& b: buffers ) {
        if( b.use_count() == 1 ) {
            return b;
        }

    }
    size_t sz = get_data_size();
    if( sz == 0 ) {
        throw std::runtime_error("Seeing::get_buffer_size() returned 0!!");
    }
    DimmSet::dimm_data_t new_buffer( new DimmSet::data_t[sz], []( DimmSet::data_t*& p){ delete[] p; p=nullptr; } );
    buffers.push_back( new_buffer );
    return new_buffer;
    
}


void Seeing::clear_buffers( void ) {
    
    std::lock_guard<std::mutex> lock(mtx);
    buffers.clear();
    
}
