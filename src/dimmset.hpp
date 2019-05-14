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

#include <string>
#include <memory>
#include <mutex>
#include <vector>

#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>


namespace wfwfs {

    class DimmSet {
        
    public:
        typedef uint16_t data_t;
        typedef std::shared_ptr<data_t> dimm_data_t;

        DimmSet( size_t i );
        DimmSet( const DimmSet& );
        DimmSet( DimmSet&& );
        
        void parsePropertyTree( boost::property_tree::ptree& cfg_ptree );
        
        size_t get_id( void ) const { return id; }
        const std::string& get_name( void ) const { return name; }
        
        Cell& get_ref_cell( void );
        const std::vector<Cell>& get_cells( void ) const { return cells; };
        bool adjust_cells(void);
        void shift_cells(PointI);
        
        uint16_t get_cell_size( void ) const { return cell_size; }
        size_t get_ref_cell_id( void ) const { return ref_cell; }
        uint16_t get_ref_cell_size( void ) const { return ref_cell_size; }
        uint8_t get_max_shift( void ) const { return max_shift; }
        size_t get_data_size( void ) const;
        
        void set_ravg( float );
        void set_min_lock( float ml ){ min_lock=ml; };
        
        void add_frame_data( boost::posix_time::ptime, dimm_data_t, size_t );
        
        void measure_shifts( uint64_t* );
        PointF get_avg_shift( size_t i ) const { if( i<avg_shifts.size() ) return avg_shifts[i]; return PointF(0.0); }
        std::vector<PointF> get_avg_shifts( void ) const { return avg_shifts; }
        void calculate_dimms( void );
        void calculate_r0( void );
        PointD calculate_r0( boost::posix_time::ptime&, uint16_t );
        boost::posix_time::ptime get_last_frametime( void );
        
    private:
        const size_t id;
        std::string name;
        std::mutex mtx;
        
        std::vector<Cell> cells;
        std::vector<PointF> avg_shifts;
        uint16_t ref_cell_size;
        uint16_t cell_size;
        uint8_t max_shift;
        size_t ref_cell;
        
        uint16_t cadence;                           // (s) How often r0 should be calculated
        uint16_t duration;                          // (s) How long timespan to accumulate statistics in the r0 calculation
        float running_average;                      // (s) Timescale for the moving average of the (slowly varying) image motion for the cells
        float min_lock;                             // (fraction) Minimum lock-ratio for including a pair in r0-calculation
        boost::posix_time::ptime last_dimm;
        boost::posix_time::ptime last_r0;           // When the latest r0 calculation was done

        struct pair_info {              // container for statistics/motion a given pair of cells
            PointF previous_shift;      // remember the relative shift of the pair
            struct snapshot_t {
                snapshot_t() : separation(0), shift(0,0), variance(0,0), rms(0,0), ok(true) {};
                float separation;
                PointF shift;           
                PointF squared_shift;           
                PointF variance;
                PointF rms;
                bool ok;
            };
            std::map<boost::posix_time::ptime, snapshot_t> data;
        };
        std::map<PointI,pair_info> differential_motion;
        
        std::map<boost::posix_time::ptime, std::vector<PointF> > cell_shifts;
        
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
