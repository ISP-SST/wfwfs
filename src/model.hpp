#ifndef WFWFS_MODEL_HPP
#define WFWFS_MODEL_HPP
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

#include <mutex>



namespace wfwfs {
        
        // Container-class for the theoretical/nominal subfield pattern.
        class Model {

        public:
            
            struct SF_Info {
                SF_Info() : enabled(true), median_gain(0) {};
                SF_Info( PointI mc) : model_coordinates(mc), enabled(true), median_gain(0) {};
                SF_Info( PointI mc, PointF nc ) : model_coordinates(mc), nominal_position(nc), enabled(true), median_gain(0) {};
                const PointF& get_position( void ) const { if( real_position.min()>0 ) return real_position; return nominal_position; };
                PointI model_coordinates;
                PointF nominal_position;
                PointF real_position;
                bool enabled;
                mutable double median_gain;
            };
            
            static Model& get(void);
            static void generate( PointI frameSize, float scale, PointF opticalAxis );
            template <typename T>
            static void make_mask( const Array<T>&, bool save=true );
            static const Array<uint8_t>& get_mask( void ){ return get().mask; };
            static const std::map<PointI,SF_Info>& get_subfields( void ) { return get().subfields; };
            static bool has_mask(void) { return (get().mask.nElements() > 0); };
            static size_t size(void) { return get().subfields.size(); };
            static PointF& get_subfield_size(void) { return get().subfield_size; };
            static void setBase( PointF b1, PointF b2 ) { get().base1 = b1; get().base2 = b2; };
            
            template <typename T>
            void draw_subfields( Array<T>& img, const T value=T() ) const;
            
        private:
            
            Model() : base1(0,1), base2(0.5*sqrt(3),0.5), subfield_size(150,180) {};
            
            void generate_( PointI frameSize, float scale, PointF opticalAxis );
            template <typename T>
            void make_mask_( const Array<T>&, bool save=true );

            PointF base1;       // Define base-vectors for the subfield grid. For a rectangular grid, use base1(1,0) & base2(0,1).
            PointF base2;       // For a hexagonal layout, as in the WFWFS, use base1(0,1) & base2(0.5*sqrt(3),0.5)
            PointF subfield_size;

            mutable std::mutex mtx;
            std::map<PointI,SF_Info> subfields;
            Array<uint8_t> mask;

        };      // Model
        

}               // namespace wfwfs


#endif // WFWFS_MODEL_HPP
