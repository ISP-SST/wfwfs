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
#include "dimmset.hpp"

#include "sads.hpp"
#include "seeing.hpp"
#include "translators.hpp"

//#define PRINT_DEBUG

using namespace wfwfs;
using namespace std;

namespace bpx = boost::posix_time;

namespace {
    
    const PointF base1( 0, 1 );                   // Define base-vectors for the microlense array .
    const PointF base2( 0.5*sqrt(3), 0.5 );       // We assume a subfield is located approximately at the center of the CCD, hexagonal symmetry,
                                                  // and that the main row is aligned horizontally.

}


DimmSet::DimmSet(size_t i) : id(i), ref_cell_size(64), cell_size(76), subcell_size(0),
                             max_shift(6), ref_cell(-1), cadence(5), duration(5),
                             running_average(20.0), min_lock(0.0),
                             last_dimm(bpx::not_a_date_time),
                             last_r0(bpx::not_a_date_time) {

    frame_data.reserve(1000);
    
}


DimmSet::DimmSet( const DimmSet& rhs ) : id(rhs.id),
                                    name(rhs.name),
                                    cells(rhs.cells),
                                    subcells(rhs.subcells),
                                    avg_shifts(rhs.avg_shifts),
                                    ref_cell_size(rhs.ref_cell_size),
                                    cell_size(rhs.cell_size),
                                    subcell_size(rhs.subcell_size),
                                    max_shift(rhs.max_shift),
                                    ref_cell(rhs.ref_cell),
                                    cadence(rhs.cadence),
                                    duration(rhs.duration),
                                    running_average(rhs.running_average),
                                    min_lock(rhs.min_lock),
                                    last_dimm(rhs.last_dimm),
                                    last_r0(rhs.last_r0),
                                    differential_motion(rhs.differential_motion),
                                    frame_data(rhs.frame_data) {


}


DimmSet::DimmSet( DimmSet&& rhs ) : id(rhs.id),
                                    name(std::move(rhs.name)),
                                    cells(std::move(rhs.cells)),
                                    subcells(std::move(rhs.subcells)),
                                    avg_shifts(std::move(rhs.avg_shifts)),
                                    ref_cell_size(rhs.ref_cell_size),
                                    cell_size(rhs.cell_size),
                                    subcell_size(rhs.subcell_size),
                                    max_shift(rhs.max_shift),
                                    ref_cell(rhs.ref_cell),
                                    cadence(rhs.cadence),
                                    duration(rhs.duration),
                                    running_average(rhs.running_average),
                                    min_lock(rhs.min_lock),
                                    last_dimm(rhs.last_dimm),
                                    last_r0(rhs.last_r0),
                                    differential_motion(std::move(rhs.differential_motion)),
                                    frame_data(std::move(rhs.frame_data)) {


}


void DimmSet::find_nominal_gridpoints( float scale, PointI detector_size ) {
    
    detector_size /= 2;

    // Find the nearest grid-point
    for( auto& c: cells ) {
        float min_dist = std::numeric_limits<float>::max();
        for( int i=-5; i<5; ++i ) {
            for( int j=-5; j<5; ++j ) {
                PointF pos_px = base1*i + base2*j;
                pos_px *= scale;
                pos_px += detector_size;
                pos_px -= c.pos;
                float this_dist = sqrt( pos_px.norm() );
                if( this_dist < min_dist ) {
                    min_dist = this_dist;
                    c.grid_pos = PointI(j,i);
                }
            }
        }
    }
    
    // pre-calculate the separations of all subfield combinations
    size_t nCells = cells.size();
    for( size_t i(0); i<nCells; ++i ) {
        if ( i == ref_cell ) continue;
        for( size_t j(i+1); j<nCells; ++j ) {
            if ( j == ref_cell ) continue;
            PointI pair(i,j);
            pair_info& pi = differential_motion[ pair ];
            pi.diff = cells[j].grid_pos - cells[i].grid_pos;
            pi.diff = base1*pi.diff.x + base2*pi.diff.y;
            pi.diff *= scale;
            pi.separation = sqrt( pi.diff.norm() );                 // Save nominal separation
            //pi.diff = cells[j].pos - cells[i].pos;                  // Store true difference vector
        }
    }


}


void DimmSet::parsePropertyTree( boost::property_tree::ptree& cfg_ptree ) {
    
    name = cfg_ptree.get<string>( "name", "seeing_"+to_string(id) );
    
    ref_cell_size = cfg_ptree.get<size_t>( "cell_size", ref_cell_size );
    max_shift = cfg_ptree.get<uint8_t>( "max_shift", max_shift );
    cell_size = ref_cell_size + 2*max_shift;
    
    cadence = cfg_ptree.get<uint16_t>( "cadence", cadence );
    duration = cfg_ptree.get<uint16_t>( "duration", duration );
    running_average = cfg_ptree.get<float>( "running_average", running_average );
    min_lock = cfg_ptree.get<float>( "min_lock", 0.0 );

    ref_cell = cfg_ptree.get<size_t>( "ref_cell", ref_cell );
    
    PointI cell_offset(0,0);
    vector<int32_t> co = cfg_ptree.get< vector<int32_t> >( "cell_offset", vector<int32_t>() );
    if( co.size() == 2 ) {
        cell_offset = PointI( co[0], co[1] );
    }
    
    vector<uint32_t> cpv = cfg_ptree.get< vector<uint32_t> >( "cell_pos", vector<uint32_t>() );
    int cell_id(0);
    if( !cpv.empty() ) {
        if( cpv.size()%2 == 0 ) {
            for( size_t i(0); i<cpv.size()-1; i += 2) {
                cells.push_back( Cell(PointI( cpv[i], cpv[i+1] )+cell_offset, cell_id++) );
            }
        } else {
            cout << "Error: cell_pos must contain an even number of values (it's row/column pairs)" << endl;
        }
    } else {
        throw runtime_error("Error: cell_pos must be defined for each dimm-set!!");
    }
    
    try {
        subcell_size = cfg_ptree.get<size_t>( "subcell_size", subcell_size );
        PointI subcell_offset(0,0);
        subcell_offset += (ref_cell_size-subcell_size)/2;
        vector<int32_t> scpv = cfg_ptree.get< vector<int32_t> >( "subcell_pos", vector<int32_t>() );
        if( !scpv.empty() ) {
            if( scpv.size()%2 == 0 ) {
                for( size_t i(0); i<scpv.size()-1; i += 2) {
                    PointI subcell_pos( scpv[i], scpv[i+1] );
                    subcell_pos += subcell_offset;
                    subcells.push_back( Cell( subcell_pos, cell_id++ ) );
                }
            } else {
                cout << "Error: subcell_pos must contain an even number of values (it's row/column pairs)" << endl;
            }
        }
    } catch( ... ) {
        cout << "Error: failed to parse subcell configuration." << endl;
    }

    if( ref_cell > cells.size() ) ref_cell = 0;    // use first cell as reference if none was specified


}



Cell& DimmSet::get_ref_cell( void ) {

    if( ref_cell < cells.size() ) {
        return cells[ref_cell];
    } else if( !cells.empty() ) {   // if no ref_cell was supplied, use the first cell
        return cells[0];
    }
    
    throw runtime_error("DimmSet: No cells supplied.");
    
}


bool DimmSet::adjust_cells( void ) {
    
    bool changed(false);
    for( auto& c: cells ) {
        if( avg_shifts.count(c.ID) && (avg_shifts[c.ID].max_abs() > 0.5) ) {
            PointI tmpPos = c.pos;
            PointF tmpAF = avg_shifts[c.ID];
            tmpAF.round();
            tmpPos += tmpAF;
            if( tmpPos != c.pos ) {
                c.pos = tmpPos;
                avg_shifts[c.ID] = 0;
                changed = true;
            }
        }
    }
    return changed;
}


void DimmSet::shift_cells( PointI s ) {

    for( auto& c: cells ) {
        c.pos += s;
    }

}


void DimmSet::check( void ) {

    // TODO
    
}


size_t DimmSet::get_data_size( void ) const {
    
    return cell_size*cell_size*cells.size();
    
}


void DimmSet::set_ravg( float r ) {
    
    running_average = r;
    
}


void DimmSet::add_frame_data( bpx::ptime ts, dimm_data_t buf, size_t offset ) {
    
    lock_guard<mutex> lock( mtx );
    frame_data.push_back( {ts, buf, offset} );
    
    if( last_r0.is_not_a_date_time() ) {    // force wait until we have enough data
        last_r0 = ts + bpx::time_duration( 0, 0, duration-cadence, 0 );
    }
    
}


void DimmSet::measure_shifts( uint64_t* tmp, double avg_interval ) {

    size_t cs2 = cell_size*cell_size;
    size_t refOffset = ref_cell*cs2 + max_shift*(cell_size+1);
    size_t N = 2*max_shift+1;
    
    float mix = 1.0 - avg_interval/running_average;

    vector<frame_data_t> fd;
    {
        unique_lock<mutex> lock( mtx );
        std::swap( frame_data, fd );
    }

    for( auto& data: fd ) {
        const data_t* dataPtr = data.buf.get() + data.offset;
        const data_t* refPtr = dataPtr + refOffset;
        map<int,PointF> shifts;
        for( auto& c: cells ) {
            int cellID = c.ID;
            if( cellID == ref_cell ) continue;
            const data_t* cellPtr = dataPtr+ cellID*cs2;
            avg_shifts[cellID] *= mix;
            PointF cell_shift(0,0);
            bool min_found = sads( refPtr, ref_cell_size, cell_size, cellPtr, cell_size, cell_size, cell_shift, tmp );
            if( cell_shift.max_abs() > max_shift ) {
                min_found = false;
            }
            if( min_found ) {
                avg_shifts[cellID] += cell_shift*(1.0-mix);
            } else {
                // if the minumum was not well-determined (at the edge of search-area), set to NAN
                cell_shift = NAN;
            }
            shifts[cellID] = cell_shift;
            if( !subcells.empty() ) {
                static const size_t subcell_search = 3;                         // TBD: hardcoded range for subfield search?
                static const int subcell_shift = (max_shift - subcell_search); 
                if( subcell_shift <= 0 ) {
                    continue;
                }
                size_t scSize = subcell_size + 2*subcell_search;
                PointI cell_offset;
                if( min_found ) {                                               // restrict search to the vicinity of the larger FOV
                    cell_offset.x = std::round( cell_shift.x );
                    cell_offset.y = std::round( cell_shift.y );
                }
                for( auto& sc: subcells ) {
                    int subcellID = c.ID + 100*sc.ID;                           // Create a unique ID for the cell/subcell combination.
                    avg_shifts[subcellID] *= mix;
                    size_t subcellOffset = sc.pos.y*cell_size + sc.pos.x;
                    const data_t* subRefPtr = refPtr + subcellOffset;
                    subcellOffset += (cell_offset.y+subcell_shift)*cell_size + cell_offset.x+subcell_shift;
                    const data_t* subcellPtr = cellPtr + subcellOffset;
                    if( subcellOffset + scSize*(cell_size+1) > cs2 ) {          // This shift will be too large, i.e. the array will go out of bounds!
                        continue;
                    }
                    PointF subcell_delta;
                    bool ok = sads( subRefPtr, subcell_size, cell_size, subcellPtr, scSize, cell_size, subcell_delta, tmp );
                    if( ok ) {
                        if( min_found ) {
                            subcell_delta -= cell_offset;
                        }
                        avg_shifts[subcellID] += subcell_delta*(1.0-mix);
                    } else {
                        // if the minumum was not well-determined (at the edge of search-area), set to NAN
                        subcell_delta = NAN;
                    }
                    shifts[subcellID] = subcell_delta;
                }
            }
        }
        unique_lock<mutex> lock( mtx );
        cell_shifts[ data.timestamp ] = std::move(shifts);
    }
    
}


void DimmSet::calculate_dim( boost::posix_time::ptime ts, const map<int,PointF>& shifts, const PointI& ij ) {

    pair_info& pi = differential_motion[ ij ];
    // The DIMM method is only valid for distances larger than approximately 2x the subaperture diameter.
    if( pi.separation < 2*Seeing::diam_px ) return;
    if( false && (pi.diff.min_abs() > 0.1*pi.diff.max_abs()) ) {    // NOTE: this is just to mimic the AO code, i.e. only allow the "almost" horizontal/vertical pairs
        return;
    }
    const PointD avg_ji = avg_shifts[ij.x] - avg_shifts[ij.y];
    try {
        PointD differential_shift = avg_ji + shifts.at(ij.y) - shifts.at(ij.x);
        auto& data = pi.data[ ts ];
        if( differential_shift.isfinite() && (differential_shift.min_abs() > 0.0) ) {
            differential_shift = differential_shift.projectOnto( pi.diff, true );          // Convert to longitudinal/transverse components w.r.t. the pair-separaion (and preserve norm),
                                                                                           //   the longitudinal part is stored in "x", transverse in "y"

            data.squared_shift = differential_shift*differential_shift;
            data.ok = true;
        } else {
            data.ok = false;
        }
    } catch ( ... ) { }
    
}


void DimmSet::calculate_dims( void ) {

    lock_guard<mutex> lock( mtx );
    size_t nCells = cells.size();

    for( auto& it: cell_shifts ) {
        if( last_dimm.is_not_a_date_time() || it.first > last_dimm ) {
            last_dimm = it.first;
            const map<int,PointF>& shifts = it.second;
            for( size_t i(0); i<nCells; ++i ) {
                if( i == ref_cell ) continue;
                for( size_t j(i+1); j<nCells; ++j ) {
                    if( j == ref_cell ) continue;
                    PointI cell_pair(i,j);
                    if( subcells.empty() ) {
                        calculate_dim( last_dimm, shifts, cell_pair );
                    } else {
                        for( auto& sc: subcells ) {
                            PointI subcell_pair = cell_pair;
                            subcell_pair += 100*sc.ID;
                            if( differential_motion.count(subcell_pair) == 0 ) {
                                differential_motion[ subcell_pair ] = differential_motion[ cell_pair ];
                            }
                            calculate_dim( last_dimm, shifts, subcell_pair );
                        }
                    }
                }
            }
        }
    }
    
    cell_shifts.clear();

}


void DimmSet::calculate_r0( void ) {

    lock_guard<mutex> lock( mtx );

    bpx::ptime next = last_r0 + bpx::time_duration( 0, 0, cadence, 0 );
    
    if( last_r0.is_not_a_date_time() || last_dimm >= next ) {
        PointF this_r0 = calculate_r0( last_dimm, duration );
        last_r0 = last_dimm;
        r0[ last_dimm ] = this_r0;
        cout << to_time_t(last_r0) << " " << ((this_r0.x+this_r0.y)/2) << " " << this_r0.x << " " << this_r0.y << endl;
        
    }
}


PointD DimmSet::calculate_r0( bpx::ptime& end, uint16_t dur, float ml ) {

    if( end.is_not_a_date_time() ) {
        end = bpx::microsec_clock::universal_time();
    }
    bpx::ptime begin = end - bpx::time_duration( 0, 0, dur, 0 );
    
#ifdef PRINT_DEBUG
    cout << "calculate_r0:    " << __LINE__ << "  dur = " << dur << "\nTime: " << to_simple_string(begin) << "  -> " << to_simple_string(end) << endl;
    vector<PointI> pair_ids;
    vector<size_t> pair_n;
    vector<size_t> pair_good;
    vector<PointF> pair_var;
    vector<PointF> pair_r0;
    vector<float> pair_sep;
#endif
    
    if( ml < 0.0 ) {    // min_lock not provided in call, use the local value.
        ml = min_lock;
    }
    
    PointD r0_sum(0,0);
    size_t nPairs(0);
    
    const PointD variance_cutoff = (2*max_shift)*(2*max_shift)/12.0;       // Take a uniform distribution over the search-space as limiting value.

    for( auto& dm: differential_motion ) {
        auto& pi = dm.second;
        // The DIMM method is only valid for distances larger than 2x the subaperture diameter.
        if( pi.separation < 2*Seeing::diam_px ) continue;
        
        PointD variance(0,0);
        size_t nTotalData(0);
        size_t nGoodData(0);
        auto it = pi.data.rbegin();
        while( (it != pi.data.rend()) && (it->first > end) ) it++;
        while( (it != pi.data.rend()) && (it->first > begin) ) {
            nTotalData++;
            if( it->second.ok ) {
                variance += it->second.squared_shift;
                nGoodData++;
            }
            it++;
        }

        if( !nTotalData || !nGoodData || (variance.min() <= 0.0) ) continue;        // No data available.
        if( static_cast<double>(nGoodData)/nTotalData < ml ) continue;              // Less data than required.

        variance /= nGoodData;

        variance = variance.min( variance_cutoff );                           // Restrict variance to < variance_cutoff
        
        PointD r0 = Seeing::apply_dimm_equations( variance, pi.separation );

#ifdef PRINT_DEBUG
        pair_ids.push_back( dm.first );
        pair_n.push_back(nTotalData);
        pair_good.push_back(nGoodData);
        pair_r0.push_back(r0);
        pair_var.push_back(variance);
        pair_sep.push_back(pi.separation);
#endif
        if( !r0.isfinite() ) continue;

        r0_sum += r0;
        nPairs++;
        
    }

    if( nPairs ) {
        r0_sum /= nPairs;
    }
    
#ifdef PRINT_DEBUG
    cout << printArray( pair_ids, "  pairs", 5 ) << endl;
    cout << printArray( pair_n, "  count", 5 ) << printArray( pair_good, "  good", 5 ) << endl;
    cout << printArray( pair_var, "    var", 5 ) << endl;
    cout << printArray( pair_sep, "    sep", 5 ) << endl;
    cout << printArray( pair_r0,  "     r0", 5 ) << endl;
    cout << "     R0=" << r0_sum << " => " << ((r0_sum.x+r0_sum.y)/2) << endl;
#endif
    
    return r0_sum;
    
}


PointD DimmSet::pair_info::getShift(void) const {
    
    PointD ret;
    
    if( !data.empty() ) {
        ret = data.rbegin()->second.squared_shift;
    }
    
    return ret;
    
}


PointD DimmSet::pair_info::getVariance( bpx::ptime begin, bpx::ptime end ) const {

    PointD variance(0,0);
    size_t nData(0);

    auto it = data.rbegin();
    while( (it != data.rend()) && (it->first > end) ) it++;
    while( (it != data.rend()) && (it->first > begin) ) {
        if( it->second.ok ) {
            variance += it->second.squared_shift;
            nData++;
        }
        it++;
    }

    if( nData ) {
        variance /= nData;
    }

    return variance;
    
}


