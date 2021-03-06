#ifndef WFWFS_SEEING_HPP
#define WFWFS_SEEING_HPP
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
#include "dimmset.hpp"
#include "frame.hpp"

#include "autosave.hpp"
#include "seeinglog.hpp"

#include <mutex>

#include <boost/asio.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/property_tree/ptree.hpp>


namespace wfwfs {
    
    class Seeing  {
        
    public:

        Seeing( void );
        
        void find_nominal_gridpoints( PointI detector_size );
        void parsePropertyTree( boost::property_tree::ptree& );
        
        template <typename T>
        void copy_cell_data( Frame&, float*, float*, size_t, int );
        
        void process( double );
        
        void precalculate( void );
        static PointD apply_dimm_equations( PointD var, float separation );
        
        void start_dimms( void );
        void stop_dimms( void );
        void start_saves( void );
        void stop_saves( void );
        void start_logs( void );
        void stop_logs( void );
        void start( void );
        void stop( void );
        void maintenance( void );
        void zero_avgs( void );
        
        void set_image_scale( float s ) { arcsecs_per_pixel = s; precalculate(); }
        void set_diam( float d ) { diam = d; precalculate(); }
        void set_lambda( float l ) { lambda = l; precalculate(); }
        void set_ravg( float r, int ds_id=-1 );
        void set_min_lock( float ml, std::string tag="" );
        
        std::string adjust_cells( int ds_id=-1 );
        std::string shift_cell( PointI, int cell_id, int ds_id=-1 );
        std::string shift_cells( PointI, int ds_id=-1 );
        std::string get_cells( int ds_id=-1 ) const;
        std::string get_shifts( void ) const;
        std::string get_ashifts( int ds_id=-1 ) const;
        std::string get_vars( int duration=0 ) const;
        std::string get_locks( void ) const;
        
        float get_avg_intensity( void ) const { return avg_intensity; }
        
        template <typename T>
        void draw_cell( Array<T>& img, PointI c, uint16_t cell_size, bool mark=false ) const;
        template <typename T>
        void draw_cells( Array<T>& img ) const;

        size_t get_data_size( void ) const;
        DimmSet::dimm_data_t get_buffer( void );
        void clear_frame_data( void );
        void clear_buffers( void );


        // Global properties of the telescope + camera
        static float pixelsize;             // sst/WFWFS = 8 \mu
        static float diam;                  // Sub-aperture diameter    (m)                 sst = 0.135 (theoretical), at CCD: 2.33 mm  (measured on 2019-04-23: 0.13245)
        static float diam_px;               // Sub-aperture diameter    (pixels)            sst = 291.25
        static float arcsecs_per_pixel;     // Image-scale              (arcsec/pixel)      sst = 0.21
        static float meters_per_pixel;      // Image scale              (m/pixel)           sst = 0.000463519
        static float lambda;                // Wavelength               (m)                 default = 500e-9
        static double radians_per_pixel;    // Pre-computed factor for converting units from pixels to radians
        static double dimm_K;               // pre-computed factor used in the DIMM equations
        
        static float int_weight;
        bool subtract_plane;
        
    private:

        std::vector<DimmSet> dimm_sets;
        std::vector<DimmSet::dimm_data_t> buffers;
        size_t dsID;
        float avg_intensity;                // Store current intensity, as measured within the DIMM cells (stored as a fraction of max-intensity)

        boost::posix_time::ptime timestamp;
        std::vector<AutoSave> saves;
        std::vector<SeeingLog> logs;
        
        std::mutex mtx;
        
    };

}   // wfwfs


#endif // WFWFS_SEEING_HPP
