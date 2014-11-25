/*
 * expressiontree_motprof.hpp
 *
 *  Created on: Oct 15, 2013
 *      Author: Erwin Aertbelien
 *
* expressiongraph library
* 
* Copyright 2014 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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

#ifndef KDL_EXPRESSIONTREE_MOTPROF_HPP
#define KDL_EXPRESSIONTREE_MOTPROF_HPP

#include <kdl/expressiontree_expressions.hpp>

namespace KDL {

/**
 * \param [in] startpos starting value of position (or path parameter s)
 * \param [in] endpos ending value of position (or path parameter s)
 * \param [in] maxvel maximum velocity (absolute value)) 
 * \param [in] maxacc maximum acceleration (absolute value)) 
 *
 * \returns duration of the motion profile.
 */
double trapezoidal_mp_duration(double startpos,double endpos,double maxvel,double maxacc);

/**
 * create an expression graph representing a trapezoidal motion profile
 * \param [in] argument input, typically representing time.
 * \param [in] startpos starting value of position (or path parameter s)
 * \param [in] endpos ending value of position (or path parameter s)
 * \param [in] maxvel maximum velocity (absolute value)) 
 * \param [in] maxacc maximum acceleration (absolute value)) 
 *
 * \returns expression graph representing this motion profile. 
 */
Expression<double>::Ptr trapezoidal_mp(
        Expression<double>::Ptr argument, 
        double startpos, double endpos, double maxvel, double maxacc
);

/**
 * create an expression graph representing an half trapezoidal motion profile
 * \param [in] argument input, typically representing time.
 * \param [in] startpos starting value of position (or path parameter s)
 * \param [in] endpos ending value of position (or path parameter s)
 * \param [in] maxvel maximum velocity (absolute value)) 
 * \param [in] maxacc maximum acceleration (absolute value)) 
 * \param [in] starting true if starting with null velocity, false otherwise
 * 
 * \returns expression graph representing this motion profile. 
 */
Expression<double>::Ptr half_trapezoidal_mp(
        Expression<double>::Ptr argument, 
        double startpos, double endpos, double maxvel, double maxacc,
        bool starting
);
 
/**
 * create an expression graph representing a trapezoidal motion profile
 * \param [in] argument input, typically representing time.
 * \param [in] startpos starting value of position (or path parameter s)
 * \param [in] endpos ending value of position (or path parameter s)
 * \param [in] maxvel maximum velocity (absolute value)) 
 * \param [in] maxacc maximum acceleration (absolute value)) 
 *
 * \returns expression graph representing this motion profile. 
 */
Expression<double>::Ptr trapezoidal_mp_fixed_duration(
        Expression<double>::Ptr argument, 
        double startpos, double endpos, double maxvel, double maxacc,
        double newduration
);
 

/**
 * Expression graph node representing a trapezoidal ( or halftrapezoidal or rectangular) motion
 * profile.  Typically initialized with factory functions such as trapezoidal_mp.
 */
class Trapezoidal_Double:
	public UnaryExpression<double, double>
{
public:
	typedef UnaryExpression<double,double> UnExpr;

    // For "running" a motion profile :
    double a1,a2,a3; // coef. from ^0 -> ^2 of first part
    double b1,b2,b3; // of 2nd part
    double c1,c2,c3; // of 3th part
    double duration;
    double t1,t2;

    // specification of the motion profile :
    double maxvel;
    double maxacc;
    double startpos;
    double endpos;
    double time;

public:
    
	Trapezoidal_Double() {}

	Trapezoidal_Double(
        double _a1,double _a2,double _a3,
        double _b1,double _b2,double _b3,
        double _c1,double _c2,double _c3,
        double _duration, double _t1, double _t2,
		const  Expression<double>::Ptr& arg );

	virtual double value();

	virtual double derivative(int i);

    virtual  Expression<double>::Ptr derivativeExpression(int i);

    virtual  Expression<double>::Ptr clone();
};

} // end of namespace KDL
#endif
