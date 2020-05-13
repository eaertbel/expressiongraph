/**
 * expressiontree_quaternion.hpp
 * expressiongraph library
 * 
 * Copyright 2019 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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
 
     **********************************************************************************
     *   QUATERNION OPERATORS/FUNCTIONS AND THEIR NOTATION:
     *
     *   arguments of the operators follow the following notation:
     *      q : quaternion (w,vec) with real part w and a vector part vec.
     *      u : pure quaternion (i.e. real part w==0), represented by Vector or vector
     *      qn : unit quaternion (i.e. modulo (norm) == 1 ), to indicate that only
     *           unit quaternions are allowed as argument.
     *      s  : scalar (double value)
     *      v  : vector
     *      vn : normalized vector
     *      R  : a rotation matrix
     *
     *   q1 + q1        : quaternion sum
     *   q1 - q1        : quaternion difference 
     *   -q             : unary negation, i.e. (-w, -vec)
     *   q*s , s*q      : product scalar-quaternion
     *   q*u , u*q      : product vector-quaternion,  vector u is interpreted as pure 
     *                    quaternion (0, u) 
     *   q1*q2          : quaternion product (w,vec)
     *                    with real part   w   : q1.w*q2.w - dot(q1.vec,q2.vec)  
     *                         vector part vec : q1.w*q2.vec + q2.w*q1.vec + cross(q1.vec, q2.vec)
     *   dot(q1,q2)     : quaternion dot product, i.e. q1.w*q2.w + dot(q1.vec,q2.vec)
     *   q/s            : division of a quaternion by a scalar (w/s, vec/s)
     *   conj(q)        : quaternion conjugate, i.e. (w,-vec)
     *   w(q)           : real part of quaternion == (q+conj(q))/2
     *   vec(q)         : vector part of quaternion == (q-conj(q))/2
     *   squared_norm(q) : conj(q)*q
     *   norm(q)        : sqrt( conj(q)*q )
     *   inv(q)         : quaternion inverse, i.e. inv(q)*q == 1,  inv(q) == conj(q)/squared_norm(q)
     *   normalized(q)  : returns a unit quaternion, q/norm(q), 
     *                    for norm(q)==0, returns quaternion with w=1 and vec=0,0,0.
     *   normalized(v)  : returns a unit vector, v/norm(v), for norm(v)==0, returns vec=0,0,0 )
     *   apply(qn,v)    : computes efficiently the rotation of vector v by unit quaternion qn.
     *                    assumes that the unit quaternion qn has norm==1.
     *   exp(q)         : quaternion exponential of q, resulting in a not necessarily normalized
     *                    quaternion.  The exponential of (w,vec=0,0,0) is equal 
     *                    to exp(w)*Quaternion(1,(0,0,0))
     *   exp(u)         : quaternion exponential of a pure quaternion represented by a vector u
     *                    resulting in a normalized quaternion.  The exponential of u=0,0,0
     *                    is equal to a quaternion with w=1 vec=(0,0,0).
     *   log(q)         : quaternion logarithm of q. -inf if norm(q)==0 (similar to log(0) ).
     *                    if norm(vec(q)) very small, resuts in log(norm(q), Vector(0,0,0)).
     *   logUnit(q)     : quaternion logarighm of a unit quaternion, resulting in a pure quaternion
     *                    represented by a vector.
     *                    if norm(vec(q)) very small, resuts in Vector(0,0,0).
     *   pow(q,s)       : quaternion q to the power s.  s is a scalar. equal to exp(s*log(q)).
     *   powUnit(qn,s)  : unit quaternion q to the power s.  Returns a normalized quaternion.               
     *                    equal to exp( s* logUnit(q) )
     *   toQuat(v)      : interprets the vector as a normal vector indicating the rotation axis,
     *                    multiplied by a full (NOT half!) angle (i.e. "axis-angle" representation) and returns
     *                    the equivalent quaternion.
     *   toQuat(v, s): returns a quaternion corresponding to a rotation around v (non necessarily
     *                     normalized) with the angle given by s.
     *   axisAngle(qn)  : returns the axis-angle representation for a given unit quaternion, i.e.
     *                    the rotation axis and the full angle (NOT half!) around that rotation axis.
     *   slerp(q1,q2,s) : spherical linear interpolation between q1 (corresponding to s==0) and q2 
     *                    (corresponding to s==1).
     *   slerpUnit(qn1,qn2,s) : spherical linear interpolation between UNIT quaternion q1 (corresponding to s==0) and 
     *                      UNIT quaternion q2 (corresponding to s==1).
     *   diff(q1,q2)     : q2-q1
     *   diffUnit(q1,q2) : difference over the sphere, i.e. rotational velocity  that you need to go from q1 to q2 in 
     *                     one unit of time, i.e. angle-axis to go from q1 to q2
     *   toQuat(R)      : returns the quaternion corresponding to the given rotation matrix
     *   toRot(qn)       : returns the rotation matrix corresponding to the given unit quaternion.
     *   quaternion(s,v),
     *   quaternion(v,s) : returns quaternion with real part s and vector part v
     *
     ****************************************************************************************************
**/



#ifndef KDL_EXPRESSIONTREE_QUATERNION_HPP
#define KDL_EXPRESSIONTREE_QUATERNION_HPP
#include <expressiongraph/quat.hpp>
#include <expressiongraph/quat_io.hpp>
#include <expressiongraph/expressiontree_expressions.hpp>
#include <expressiongraph/expressiontree_double.hpp>
#include <expressiongraph/expressiontree_vector.hpp>

/**
 *
 * Quaternion class for both value and derivative
 *
 */

namespace KDL {
   

    /**********************************************************************************
     *   OPERATORS
     **********************************************************************************/



Expression<Quaternion>::Ptr conj( Expression<Quaternion>::Ptr a);
Expression<double>::Ptr w( Expression<Quaternion>::Ptr a);
Expression<Vector>::Ptr vec( Expression<Quaternion>::Ptr a);
Expression<double>::Ptr squared_norm( Expression<Quaternion>::Ptr a);
Expression<double>::Ptr norm( Expression<Quaternion>::Ptr a);
Expression<Quaternion>::Ptr operator+( Expression<Quaternion>::Ptr a1, Expression<Quaternion>::Ptr a2 );
Expression<Quaternion>::Ptr operator-( Expression<Quaternion>::Ptr a1, Expression<Quaternion>::Ptr a2 );
Expression<Quaternion>::Ptr operator-( Expression<Quaternion>::Ptr a);
Expression<Quaternion>::Ptr operator* ( Expression<Quaternion>::Ptr a1, Expression<double>::Ptr a2 );
Expression<Quaternion>::Ptr operator* ( Expression<double>::Ptr a2 , Expression<Quaternion>::Ptr a1);
Expression<Quaternion>::Ptr operator* ( Expression<Quaternion>::Ptr a1, Expression<Quaternion>::Ptr a2 );
Expression<Quaternion>::Ptr operator* ( Expression<Vector>::Ptr a1, Expression<Quaternion>::Ptr a2 );
Expression<Quaternion>::Ptr operator* ( Expression<Quaternion>::Ptr a1, Expression<Vector>::Ptr a2 );
Expression<Quaternion>::Ptr operator/ ( Expression<Quaternion>::Ptr a1, Expression<double>::Ptr a2 );
Expression<Quaternion>::Ptr inv( Expression<Quaternion>::Ptr q);
Expression<double>::Ptr dot ( Expression<Quaternion>::Ptr a1, Expression<Quaternion>::Ptr a2 );
Expression<Quaternion>::Ptr normalized ( Expression<Quaternion>::Ptr a1);
Expression<Vector>::Ptr normalized ( Expression<Vector>::Ptr a1);
Expression<Vector>::Ptr apply ( Expression<Quaternion>::Ptr q, Expression<Vector>::Ptr v);
Expression<Quaternion>::Ptr exp ( Expression<Quaternion>::Ptr a1);
Expression<Quaternion>::Ptr exp ( Expression<Vector>::Ptr a1);
Expression<Quaternion>::Ptr log ( Expression<Quaternion>::Ptr a1);
Expression<Vector>::Ptr logUnit ( Expression<Quaternion>::Ptr a1);
Expression<Quaternion>::Ptr pow ( Expression<Quaternion>::Ptr q, Expression<double>::Ptr s);
Expression<Quaternion>::Ptr powUnit ( Expression<Quaternion>::Ptr q, Expression<double>::Ptr s);
Expression<Vector>::Ptr diffUnit( Expression<Quaternion>::Ptr q1, Expression<Quaternion>::Ptr q2);

/*
    inline Expression<Quaternion>::Ptr diff( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2) {
        return (q2-q1);
    }
*/


Expression<Quaternion>::Ptr toQuat( Expression<Vector>::Ptr angle_and_axis) ;
Expression<Quaternion>::Ptr toQuat( Expression<Vector>::Ptr axis, Expression<double>::Ptr angle);
Expression<Vector>::Ptr axisAngle( Expression<Quaternion>::Ptr q);
Expression<Quaternion>::Ptr slerp( Expression<Quaternion>::Ptr q1, 
        Expression<Quaternion>::Ptr q2, 
        Expression<double>::Ptr s);
Expression<Quaternion>::Ptr slerpUnit( 
        Expression<Quaternion>::Ptr q1, 
        Expression<Quaternion>::Ptr q2, 
        Expression<double>::Ptr s);
Expression<Quaternion>::Ptr toQuat( Expression<Rotation>::Ptr a);
Expression<Rotation>::Ptr toRot( Expression<Quaternion>::Ptr a);
Expression<Quaternion>::Ptr quaternion ( 
        Expression<double>::Ptr a1, 
        Expression<Vector>::Ptr a2 );
Expression<Quaternion>::Ptr quaternion ( 
        Expression<Vector>::Ptr a1, 
        Expression<double>::Ptr a2 );


}; // namespace KDL

#endif



