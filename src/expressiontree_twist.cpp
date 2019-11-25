/*
 * expressiontree_twist.cpp
 *
 *  Created on: Sept., 2012
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

#include <kdl/expressiontree_vector.hpp>
#include <kdl/expressiontree_rotation.hpp>
#include <kdl/expressiontree_twist.hpp>
#include <kdl/frames.hpp>

namespace KDL {

Expression<Twist>::Ptr Twist_VectorVector::derivativeExpression(int i)
{
    return twist( argument1->derivativeExpression(i), argument2->derivativeExpression(i) );
}


Expression<Twist>::Ptr Negate_Twist::derivativeExpression(int i) {
    return -argument->derivativeExpression(i);
}
Expression<Vector>::Ptr Velocity_Twist::derivativeExpression(int i) {
    return transvel(argument->derivativeExpression(i));
}
Expression<Vector>::Ptr RotVelocity_Twist::derivativeExpression(int i) {
    return rotvel(argument->derivativeExpression(i));
}

Expression<Twist>::Ptr Addition_TwistTwist::derivativeExpression(int i) {
    return argument1->derivativeExpression(i) + argument2->derivativeExpression(i);
}

Expression<Twist>::Ptr Subtraction_TwistTwist::derivativeExpression(int i) {
    return argument1->derivativeExpression(i) - argument2->derivativeExpression(i);
}

Expression<Twist>::Ptr Composition_RotationTwist::derivativeExpression(int i) {
    Expression<Rotation>::Ptr a = cached<Rotation>( argument1 ); 
    Expression<Twist>::Ptr    b = cached<Twist>( argument2 ); 
    Expression<Vector>::Ptr  da = cached<Vector>(argument1->derivativeExpression(i));
    Expression<Twist>::Ptr   db = cached<Twist>(argument2->derivativeExpression(i));
    return twist( a*transvel(db) + da*(a*transvel(b)),
        			a*rotvel(db) + da*(a*rotvel(b))
	);
}
Expression<Twist>::Ptr Multiplication_TwistDouble::derivativeExpression(int i) {
   return argument1 * argument2->derivativeExpression(i) +
          argument1->derivativeExpression(i) * argument2;
}

Expression<Twist>::Ptr RefPoint_TwistVector::derivativeExpression(int i) {
    Expression<Twist>::Ptr da = cached<Twist>( argument1->derivativeExpression(i) );
    return twist( transvel(da) + rotvel(da)*argument2 + 
                  rotvel(argument1)*argument2->derivativeExpression(i),
                  rotvel(da) );
}

} // end of namespace KDL

