/*
 * mptrap.hpp
 *
 *  Created on: Jan, 2016
 *      Author: Erwin Aertbelien
 *
*  Implementation of trapezoidal motion profile including derivatives.
*
* Copyright 2016 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
*
* Licensed under the EUPL, Version 1.1 only (the "Licence");
* You may not use this work except in compliance with the Licence.
* You may obtain a copy of the Licence at:
*
* http://ec.europa.eu/idabc/eupl 
*
* Unless required by applicable law or agreed to in writing, software 
* distributed under the Licence is distributed on an "AS IS" basis,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the Licence for the specific language governing permissions and 
* limitations under the Licence.
*/

#ifndef KDL_MPTRAP_HPP
#define KDL_MPTRAP_HPP

#include "expressiontree_expressions.hpp"
#include "expressiontree_mimo.hpp"
#include <stdexcept>

namespace KDL {

    /**
     * MPTrap implements a trapezoidal motionprofile
     * ( implementation class, not be used externally )
     */

    struct MPTrap {
        double amax;   ///< maximum acceleration
        double vmax;   ///< maximum velocity
        double spos;   ///< start position
        double epos;   ///< end position
        double duration;
        double t1;
        double t2;
        double s;
        double dT;
        double d_t1_d_dpos;
        double d_t2_d_dpos;
        double d_duration_d_dpos;

        MPTrap();

        void setPlan(double spos, double epos, double vmax, double amax);

        double planMinDuration();
        double adaptDuration(double new_duration);
        /**
         * derivatives only to be used when this is the critical motion profile
         */
        void compute_derivs();
        double pos(double time);
        double d_pos_d_time(double time);
        double d_pos_d_spos(double time); 
        double d_pos_d_epos(double time);
        double d_pos_d_t(double time);
    };


} // end of namespace KDL
#endif
