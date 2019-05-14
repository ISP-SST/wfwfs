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

#include "seeinglog.hpp"

#include <mutex>

#include <boost/asio.hpp>
#include <boost/property_tree/ptree.hpp>


namespace wfwfs {
    
    class Seeing  {
        
    public:

        Seeing( void );
        
        void parsePropertyTree( boost::property_tree::ptree& );
        
        template <typename T>
        void copy_cell_data( Frame&, float*, float*, size_t, int );
        
        void process( boost::asio::io_service& );
        
        void precalculate( void );
        static PointD apply_dimm_equations( PointD var, float separation );
        static PointD simple_dimm_equations( PointD var, float separation );
        
        void start_logs( void );
        void stop_logs( void );
        
        void set_image_scale( float s ) { arcsecs_per_pixel = s; precalculate(); }
        void set_diam( float d ) { diam = d; precalculate(); }
        void set_lambda( float l ) { lambda = l; precalculate(); }
        
        std::string adjust_cells( size_t id=0 );
        std::string shift_cells( PointI, size_t id=0 );
        std::string get_cells( size_t id=0 ) const;
        std::string get_shifts( size_t id=0 ) const;
        
        template <typename T>
        void draw_cell( Array<T>& img, PointI c, uint16_t cell_size, bool mark=false ) const;
        template <typename T>
        void draw_cells( Array<T>& img ) const;

        size_t get_data_size( void ) const;
        DimmSet::dimm_data_t get_buffer( void );
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
        
    private:

        std::vector<DimmSet> dimm_sets;
        std::vector<DimmSet::dimm_data_t> buffers;
        size_t dsID;

        std::vector<SeeingLog> logs;
        
        std::mutex mtx;
        
    };

}   // wfwfs


#endif // WFWFS_SEEING_HPP
