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

#include "utils.hpp"

#include <numeric>

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

float Seeing::int_weight(1.0);


Seeing::Seeing( void ) : subtract_plane(false), dsID(1), avg_intensity(0.0) {

    precalculate();
    
}


void Seeing::find_nominal_gridpoints( PointI detector_size ) {
    
     for( auto& ds: dimm_sets ) {
            ds.find_nominal_gridpoints( diam_px, detector_size );
    }
   
}


void Seeing::parsePropertyTree( boost::property_tree::ptree& cfg_ptree ) {
    
    string output_base = cfg_ptree.get<string>( "outputdir", "/data/" );
    arcsecs_per_pixel = cfg_ptree.get<float>( "arcsecs_per_pixel", arcsecs_per_pixel );
    diam = cfg_ptree.get<float>( "subaperture_diameter", diam );
    diam_px = cfg_ptree.get<float>( "subaperture_pixels", diam_px );
    lambda = cfg_ptree.get<float>( "wavelength", lambda );
    pixelsize = cfg_ptree.get<float>( "pixelsize", pixelsize );
    subtract_plane = cfg_ptree.count("subtract_plane");
    
    precalculate();

    for( auto& node : cfg_ptree ) {
        if( iequals( node.first, "DIMM" ) ) {
            DimmSet ds( dsID++ );
            ds.parsePropertyTree( node.second );
            dimm_sets.push_back(std::move(ds));
        } else if( iequals( node.first, "LOG" ) ) {
            SeeingLog sl;
            sl.parsePropertyTree( node.second, output_base );
            logs.push_back( std::move(sl) );
        } else if( iequals( node.first, "SAVE" ) ) {
            AutoSave sv;
            sv.parsePropertyTree( node.second, output_base );
            saves.push_back( std::move( sv ) );
        }
    }
    
    string default_columns;
    for( auto& ds: dimm_sets ) {    // add default logging (all defined DIMM sets)
        default_columns += ds.get_name() + " ";
    }
    boost::trim( default_columns );
    if( !default_columns.empty() ) {
        for( auto& l: logs ) {
            if( !l.hasColumns() ) {
                l.addColumns( default_columns );
            }
        }
    }

    
}


template <typename T>
void Seeing::copy_cell_data( Frame& f, float* dd, float* gg, size_t stride, int maxval ) {

    auto buf = get_buffer();
    T* in = reinterpret_cast<T*>( f.data );
    auto out = buf.get();
    size_t ds_offset(0);

    memset( out, 0, get_data_size() );
    
    double total_sum(0.0);
    size_t nTotalSummed(0);
        
    for( auto& ds: dimm_sets ) {
        ds.add_frame_data( f.timestamp, buf, ds_offset );
        
        
        auto& cells = ds.get_cells();
        size_t cs = ds.get_cell_size();
        size_t cs2 = cs*cs;
        size_t imgStride = cs2;
        size_t c_offset(0);
        for( size_t n(0); n<cells.size(); ++n ) {
            const Cell& c = cells[n];
            double cell_sum(0.0);
            size_t n_cell_summed(0);
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
                cell_sum = std::accumulate( outPtr, outPtr+this_cell_size, cell_sum );
                n_cell_summed += this_cell_size;
                outPtr += cs;
            }
            total_sum += cell_sum;
            nTotalSummed += n_cell_summed;
            if( n_cell_summed ) cell_sum /= n_cell_summed;
            if( subtract_plane ) {  // subtract plane
                outPtr = out + ds_offset + c_offset;
                if( n == ds.get_ref_cell_id() ) {
                    outPtr += ds.get_max_shift()*(cs+1);
                }
                vector<double> slopes = fitPlane( outPtr, this_cell_size, this_cell_size, cs );
                for( int y(0); y < (int)this_cell_size; ++y ) {    // for each cell-row
                    double delta_y = y*slopes[1];
                    for( int x(0); x < (int)this_cell_size; ++x ) {
                        outPtr[x] -= x*slopes[0] + delta_y + slopes[2];     // subtract plane & avg
                    }
                    outPtr += cs;
                }
            } else {    // subtract avg
                outPtr = out + ds_offset + c_offset;
                if( n == ds.get_ref_cell_id() ) {
                    outPtr += ds.get_max_shift()*(cs+1);
                }
                for( size_t y(0); y < this_cell_size; ++y ) {    // for each cell-row
                    for( size_t x(0); x < this_cell_size; ++x ) {
                        outPtr[x] -= cell_sum;
                    }
                    outPtr += cs;
                }
            }
            c_offset += imgStride;
        }
        ds_offset += c_offset;
    }
    
    if( nTotalSummed ) {
        total_sum /= nTotalSummed;
    }
    
    const float mix = 0.01;
    avg_intensity *= (1.0-mix);
    avg_intensity += mix*total_sum*int_weight;
    
}
template void Seeing::copy_cell_data<uint8_t>( Frame&, float*, float*, size_t, int );
template void Seeing::copy_cell_data<uint16_t>( Frame&, float*, float*, size_t, int );


void Seeing::process( double avg_interval ) {
    
    uint64_t tmp[256*256];
    
    for( auto& ds: dimm_sets ) {
//        ios.post( [&](){
//            uint64_t tmp[256*256];
            ds.measure_shifts( tmp, avg_interval );
            ds.calculate_dims();
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


void Seeing::start_dimms( void ) {
    
    for( auto& d: dimm_sets ) {
        d.start();
    }
    
}


void Seeing::stop_dimms( void ) {
    
    for( auto& d: dimm_sets ) {
        d.stop();
    }
    
}


void Seeing::start_saves( void ) {
    
    for( auto& s: saves ) {
        s.start( dimm_sets );
    }
    
}


void Seeing::stop_saves( void ) {
    
    for( auto& s: saves ) {
        s.stop();
    }
    
}


void Seeing::start_logs( void ) {
    
    for( auto& l: logs ) {
        l.start( dimm_sets );
    }
    
}


void Seeing::stop_logs( void ) {
    
    for( auto& l: logs ) {
        l.stop();
    }
    
}


void Seeing::start( void ) {
    
    timestamp = bpx::second_clock::universal_time();

    start_dimms();
    start_logs();
    start_saves();
    
}


void Seeing::stop( void ) {
    
    stop_logs();
    stop_saves();
    stop_dimms();

}


void Seeing::maintenance( void ) {
    
    bpx::ptime now = bpx::second_clock::universal_time();  
    bool new_day = now.date() != timestamp.date();
    
    if( new_day ) {
        for( auto& ds: dimm_sets ) {
            ds.reset();
        }
        for( auto& l: logs ) {
            l.reopen();
        }
        timestamp = now;
    }
    
    
}


void Seeing::zero_avgs( void ) {
    
    for( auto& ds: dimm_sets ) {
        ds.zero_avgs();
    }
    
}


void Seeing::set_ravg( float r, int ds_id  ) {
    
    for( auto& ds: dimm_sets ) {
        if( (ds_id < 0) || (ds_id == ds.get_id()) ) {
            ds.set_ravg( r );
        }
    }
    
}


void Seeing::set_min_lock( float ml, std::string tag ) {
    
    if( ml < 0.0 || ml > 1.0 ) return;
    
    for( auto& ds: dimm_sets ) {
        if( tag.empty() || (tag == ds.get_name()) ) {
           ds.set_min_lock( ml );
        }
    }
    
    for( auto& l: logs ) {
        if( tag.empty() || (tag == l.get_name()) ) {
           l.set_min_lock( ml );
        }
    }
    
}


string Seeing::adjust_cells( int ds_id ) {

    bool changed(false);
    for( auto& ds: dimm_sets ) {
        if( (ds_id < 0) || (ds_id == ds.get_id()) ) {
           changed |= ds.adjust_cells();
        }
    }
    if( changed ) {
        return get_cells(ds_id);
    }
    return "";
}


string Seeing::shift_cell( PointI s, int cell_id, int ds_id ) {

    for( auto& ds: dimm_sets ) {
        if( (ds_id < 0) || (ds_id == ds.get_id()) ) {
           ds.shift_cell(cell_id,s);
        }
    }

    return get_cells(ds_id);

}


string Seeing::shift_cells( PointI s, int ds_id ) {

    for( auto& ds: dimm_sets ) {
        if( (ds_id < 0) || (ds_id == ds.get_id()) ) {
           ds.shift_cells(s);
        }
    }

    return get_cells(ds_id);

}


string Seeing::get_cells( int ds_id ) const {
    
    string ret;
    
    for( auto& ds: dimm_sets ) {
        if( (ds_id < 0) || (ds_id == ds.get_id()) ) {
            ret += "\n";
            ret += printArray( ds.get_cells(), ds.get_name() );
        }
    }
    
    return ret;
}


string Seeing::get_shifts( void ) const {
    
    string ret;
    for( auto& ds: dimm_sets ) {
        ret += "\n";
        ret += ds.get_name() + "=[";
        for( auto& dm: ds.get_differential_motion() ) {
            ret += (string)(dm.first) + ":" + (string)dm.second.getShift( ) + ",";
        }
        ret += "]";
    }
    
    return ret;
    
}


string Seeing::get_ashifts( int ds_id ) const {
    
    string ret;
    for( auto& ds: dimm_sets ) {
        if( (ds_id < 0) || (ds_id == ds.get_id()) ) {
            ret += "\n";
            ret += ds.get_name() + "=[";
            for( auto& as: ds.get_avg_shifts() ) {
                ret += to_string(as.first) + ":" + (string)as.second + ",";
            }
            ret += "]";
        }
    }
    
    return ret;
    
}


string Seeing::get_vars( int duration ) const {
    
    string ret;
    duration = std::min<int>( duration, 1 );
    
    bpx::ptime end = bpx::microsec_clock::universal_time();
    bpx::ptime begin = end - bpx::time_duration( 0, 0, duration, 0 );
    
    for( auto& ds: dimm_sets ) {
        ret += "\n";
        ret += ds.get_name() + "=[";
        for( auto& dm: ds.get_differential_motion() ) {
            ret += (string)(dm.first) + ":" + (string)dm.second.getVariance( begin, end ) + ",";
        }
        ret += "]";
    }
    
    return ret;
}


string Seeing::get_locks( void ) const {
    
    string ret;
    for( auto& ds: dimm_sets ) {
        auto& locks = ds.get_locks();
        ret += "\n";
        ret += ds.get_name() + "=[";
        for( auto& l: locks ) {
            ret += to_string(l.first) + ":" + to_string(l.second)+ ",";
        }
        ret += "]";
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

namespace {
    
    const PointF base1( 0, 1 );                   // Define base-vectors for the microlense array .
    const PointF base2( 0.5*sqrt(3), 0.5 );       // We assume a subfield is located approximately at the center of the CCD, hexagonal symmetry,
                                                  // and that the main row is aligned horizontally.

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

    for( int i=-5; i<5; ++i ) {
        for( int j=-5; j<5; ++j ) {
            PointF pos_px = base1*i + base2*j;
            pos_px *= diam_px;
            if( sqrt(pos_px.norm()) < 1040 ) {
                pos_px += PointI(1040,1040);
                draw_cell( img, pos_px, 2, false );
            }
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

void Seeing::clear_frame_data( void ) {
    
    for( auto& ds: dimm_sets ) {
        ds.clear_frame_data();
    }

}


void Seeing::clear_buffers( void ) {
    
    std::lock_guard<std::mutex> lock(mtx);
    buffers.clear();
    
}
