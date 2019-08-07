#ifndef WFWFS_DIMMSET_HPP
#define WFWFS_DIMMSET_HPP
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

#include "cell.hpp"

#include <atomic>
#include <string>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>


namespace wfwfs {

    class DimmSet {
    public:
        struct pair_info {              // container for statistics/motion a given pair of cells
            PointF diff;
            float separation;
            struct snapshot_t {
                snapshot_t() : ok(false) {};
                PointF squared_shift;           
                bool ok;
            };
            PointD getShift(void) const;
            PointD getVariance( boost::posix_time::ptime begin, boost::posix_time::ptime end ) const;
            std::map<boost::posix_time::ptime, snapshot_t> data;
        };
        
        typedef int32_t data_t;
        typedef std::shared_ptr<data_t> dimm_data_t;

        DimmSet( int i );
        DimmSet( const DimmSet& );
        DimmSet( DimmSet&& );
        
        void find_nominal_gridpoints( float scale, PointI detector_size );
        
        void parsePropertyTree( boost::property_tree::ptree& cfg_ptree );
        
        int get_id( void ) const { return id; }
        const std::string& get_name( void ) const { return name; }
        
        Cell& get_ref_cell( void );
        std::vector<Cell>& get_cells( void ) { return cells; };
        const std::vector<Cell>& get_cells( void ) const { return cells; };
        bool adjust_cells(void);
        void shift_cell( int cell_id, PointI shift );
        void shift_cells(PointI);
        void check( void );
        
        uint16_t get_cell_size( void ) const { return cell_size; }
        size_t get_ref_cell_id( void ) const { return ref_cell; }
        uint16_t get_ref_cell_size( void ) const { return ref_cell_size; }
        uint8_t get_max_shift( void ) const { return max_shift; }
        size_t get_data_size( void ) const;
        
        void clear_frame_data( void );
        void zero_avgs( void );
        void set_ravg( float );
        void set_min_lock( float ml ){ min_lock=ml; };
        
        void add_frame_data( boost::posix_time::ptime, dimm_data_t, size_t );
        
        void measure_shifts( uint64_t*, double );
        PointF get_avg_shift( int i ) const { if( avg_shifts.count(i) ) return avg_shifts.at(i); return PointF(0.0); }
        const std::map<int,PointF>& get_avg_shifts( void ) const { return avg_shifts; }
        const std::map<int,PointF> get_shifts( void ) const { if( !cell_shifts.empty() ) return cell_shifts.rbegin()->second; return std::map<int,PointF>(); }
        const std::map<PointI,pair_info>& get_differential_motion( void ) const { return differential_motion; }
        void calculate_dim( boost::posix_time::ptime, const std::map<int,PointF>&, const PointI& );
        void calculate_dims( void );
        void calculate_r0( void );
        PointD calculate_r0( boost::posix_time::ptime& );
        PointD get_r0( boost::posix_time::ptime, uint16_t span ) const;
        const std::map<int,float>& get_locks( void ) const { return locks; };
        std::chrono::high_resolution_clock::time_point get_next_time( void );
        
        void start( void );
        void stop( void );
        void reset( void );

    private:
        
        const int id;
        std::string name;
        mutable std::mutex mtx;
        std::thread trd;
        std::atomic<bool> running;
        
        std::vector<Cell> cells;
        std::vector<Cell> subcells;
        std::map<int,PointF> avg_shifts;
        std::map<int,float> locks;
        uint16_t ref_cell_size;
        uint16_t cell_size;
        uint16_t subcell_size;
        uint8_t max_shift;
        int ref_cell;
        
        uint16_t interval;                          // (s) How often r0 should be calculated
        uint16_t duration;                          // (s) How long timespan to accumulate statistics in the r0 calculation
        float running_average;                      // (s) Timescale for the moving average of the (slowly varying) image motion for the cells
        float min_lock;                             // (fraction) Minimum lock-ratio for including a pair in r0-calculation

        boost::posix_time::ptime last_dimm;
        boost::posix_time::ptime last_r0;           // When the latest r0 calculation was done

        std::map<PointI,pair_info> differential_motion;
        
        std::map<boost::posix_time::ptime, std::map<int,PointF> > cell_shifts;
        
        std::map<boost::posix_time::ptime, PointF > r0;
        
        struct frame_data_t {
            boost::posix_time::ptime timestamp;
            dimm_data_t buf;
            size_t offset;
        };
        std::vector<frame_data_t> frame_data;

        
    };
    

}   // wfwfs


#endif // WFWFS_DIMMSET_HPP
