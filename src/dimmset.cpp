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
                                    avg_shifts(rhs.avg_shifts),
                                    ref_cell_size(rhs.ref_cell_size),
                                    cell_size(rhs.cell_size),
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
                                    avg_shifts(std::move(rhs.avg_shifts)),
                                    ref_cell_size(rhs.ref_cell_size),
                                    cell_size(rhs.cell_size),
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
    if( !cpv.empty() ) {
        if( cpv.size()%2 == 0 ) {
            for( size_t i(0); i<cpv.size()-1; i += 2) {
                cells.push_back( PointI( cpv[i], cpv[i+1] )+cell_offset );
            }
        } else {
            cout << "Error: cell_pos must contain an even number of values (it's row/column pairs)" << endl;
        }
    } else {
        throw runtime_error("Error: cell_pos must be defined for each dimm-set!!");
    }
    
    if( ref_cell > cells.size() ) ref_cell = 0;    // use first cell as reference if none was specified

    avg_shifts.resize( cells.size(), 0.0 );

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
    size_t ns = avg_shifts.size();
    if( ns == cells.size() ) {
        for( size_t i(0); i<ns; ++i ) {
            avg_shifts[i].round();
            PointI tmp = cells[i].pos;
            tmp += avg_shifts[i];
            if( tmp != cells[i].pos ) {
                cells[i].pos = tmp;
                avg_shifts[i] = 0;
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


void DimmSet::measure_shifts( uint64_t* tmp ) {

    size_t nCells = cells.size();
    size_t cs2 = cell_size*cell_size;
    size_t refOffset = ref_cell*cs2;
    
    float framerate = 10;   // FIXME
    float mix = powf(0.05, 1.0 / framerate / running_average);
    
    vector<frame_data_t> fd;
    {
        unique_lock<mutex> lock( mtx );
        std::swap( frame_data, fd );
    }
    for( auto& data: fd ) {
        const data_t* dataPtr = data.buf.get() + data.offset;
        const data_t* refPtr = dataPtr + refOffset;
        vector<PointF> shifts( nCells, 0.0 );
        for( size_t n(0); n<nCells; ++n ) {
            if( n != ref_cell ) {
                const data_t* cPtr = dataPtr+n*cs2;
                bool ok = sads( refPtr, ref_cell_size, cPtr, cell_size, shifts[n], tmp );
                if( !ok ) {
                    // if the minumum was not well-determined (at the edge of search-area), set to NAN to exclude from average.
                    shifts[n] = NAN;
                } else {
                    // update running averages
                    avg_shifts[n] *= mix;
                    avg_shifts[n] += shifts[n]*(1.0-mix);
                }
            }
        }
        unique_lock<mutex> lock( mtx );
        cell_shifts[ data.timestamp ] = std::move(shifts);
    }
    
}


/*PointF DimmSet::get_avg_shift( size_t cell_id, bpx::ptime end ) const {
    
    bpx::ptime begin = end - bpx::time_duration( 0, 0, running_average, 0 );
    PointF ret(0,0);
    size_t count(0);

    for( auto it: cell_shifts ) {
        if( (it.first > begin) && (it.first < end) &&
            (cell_id < it.second.size()) && it.second[cell_id].isfinite()) {
            ret += it.second[cell_id];
            count++;
        }
    }
    if( count ) {
        ret /= count;
    }

    return ret;
    
}


vector<PointF> DimmSet::get_avg_shifts( bpx::ptime end ) const {
    
    bpx::ptime begin = end - bpx::time_duration( 0, 0, running_average, 0 );
    size_t nCells = cells.size();
    vector<PointF> ret( nCells, {0.0, 0.0} );
    size_t count(0);

    for( auto it: cell_shifts ) {
        if( (it.first > begin) && (it.first < end) && (nCells == it.second.size()) ) {
            transform( it.second.begin(), it.second.end(), ret.begin(), ret.begin(),
                [&count]( const PointF& a, const PointF& b ) {
                    if( a.isfinite() ) {
                        count++;
                        return a+b;
                    }
                    return b;
                }
            );
        }
    }
    if( count ) {
        for( auto& i: ret ) i /= count;
    }

    return ret;
    
}*/


void DimmSet::calculate_dimms( void ) {

    lock_guard<mutex> lock( mtx );
    size_t nCells = cells.size();

    for( auto& it: cell_shifts ) {
        if( last_dimm.is_not_a_date_time() || it.first > last_dimm ) {

            last_dimm = it.first;
            vector<PointF> shifts = get_avg_shifts();
            if( shifts.size() < cells.size() ) shifts.resize( cells.size(), 0.0 );

            for( size_t i(0); i<nCells; ++i ) {
                if ( i == ref_cell ) continue;
                for( size_t j(i+1); j<nCells; ++j ) {
                    if ( j == ref_cell ) continue;
                    PointD diff = (cells[j].pos+shifts[j]) - (cells[i].pos+shifts[i]);
                    if( diff.min_abs() > 0.1*diff.max_abs() ) {    // FIXME: this is just to mimic the AO code, i.e. only allow the "almost" horizontal/vertical pairs
                        //continue;
                    }
                    PointI pair(i,j);
                    pair_info& pi = differential_motion[ pair ];
                    auto& data = pi.data[ it.first ];
                    
                    data.separation = sqrt( diff.norm() );
                    PointD differential_shift = (it.second[i] - shifts[i])
                                              - (it.second[j] - shifts[j]);
                    if( !differential_shift.isfinite() ) {
                        data.ok = false;
                        continue;
                    }
                    if( true ) {    // 
                        differential_shift = differential_shift.projectOnto( diff, true );          // convert to local longitudinal/transverse components (and preserve norm)
                                                                                                    // The longitudinal part is stored in "x", transverse in "y"
                    } else {                                                
                        if( diff.x < diff.y ) {                                                     // Use the old way (hardcoded longitudinal/transverse directions)
                            std::swap( differential_shift.x, differential_shift.y );
                        }
                    }
                    PointD differential_shift_change = differential_shift - pi.previous_shift;      // Estimate noise by summing change in motion squared
                    data.shift = differential_shift;
                    data.squared_shift = differential_shift*differential_shift;
                    data.variance = differential_shift*differential_shift;
                    data.rms = differential_shift_change*differential_shift_change;
 

                    pi.previous_shift = differential_shift;                                         // ... and save for next call.
                    
                }
            }

        }
    }

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


PointD DimmSet::calculate_r0( bpx::ptime& end, uint16_t dur ) {

    if( end.is_not_a_date_time() ) {
        end = get_last_frametime();
    }
    bpx::ptime begin = end - bpx::time_duration( 0, 0, dur, 0 );

#ifdef PRINT_DEBUG
    cout << "calculate_r0:    " << __LINE__ << "  Time: " << to_simple_string(begin) << "  -> " << to_simple_string(end) << endl;
    vector<PointI> pair_ids;
    vector<size_t> pair_n;
    vector<PointF> pair_r0;
    vector<float> pair_dist;
#endif
    
    PointD r0_sum(0,0);
    size_t nPairs(0);
    size_t nTotalData(0);
    
    for( auto& dm: differential_motion ) {
        auto& pi = dm.second;
        PointD variance(0,0);
        float separation_sum(0);
        size_t nPairData(0);
        for( auto& d: pi.data ) {
            if( d.second.ok && (d.first > begin) && (d.first <= end) ) {
                variance += d.second.squared_shift;
                separation_sum += d.second.separation;
                nPairData++;
            }
        }
        if( nPairData < 10 ) continue;
        
        separation_sum /= nPairData;
        variance /= nPairData;
        
        // The DIMM method is only valid for distances larger than 2x the subaperture diameter.
        if( separation_sum < 2*Seeing::diam_px ) continue;

        PointD r0 = Seeing::apply_dimm_equations( variance, separation_sum );
        
#ifdef PRINT_DEBUG
        pair_ids.push_back( dm.first );
        pair_n.push_back(nPairData);
        pair_r0.push_back(r0);
        pair_dist.push_back(separation_sum);
#endif

        r0_sum += r0;
        nPairs++;
        nTotalData += nPairData;
        
    }

    if( nPairs ) {
        r0_sum /= nPairs;
    }
    
#ifdef PRINT_DEBUG
    cout << printArray( pair_ids, "  pairs", 5 ) << endl;
    cout << printArray( pair_n, "  count", 5 ) << endl;
    cout << printArray( pair_dist, "   dist", 5 ) << endl;
    cout << printArray( pair_r0, "     r0", 5 ) << endl;
    cout << "     R0=" << r0_sum << " => " << ((r0_sum.x+r0_sum.y)/2) << endl;
#endif
    
    return r0_sum;
    
}


boost::posix_time::ptime DimmSet::get_last_frametime( void ) {
    
    if( cell_shifts.empty() ) {
        return bpx::not_a_date_time;
    }
    return cell_shifts.rbegin()->first;
    
}
