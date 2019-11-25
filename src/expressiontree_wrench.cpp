/*
 * expressiontree_wrench.cpp
 *
 *  Created on: Sept. 2012 
 *      Author: Erwin Aertbelien - Wouter Bancken
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

#include <kdl/expressiontree_wrench.hpp>
#include <kdl/expressiontree_vector.hpp>
#include <kdl/expressiontree_rotation.hpp>
#include <kdl/frames.hpp>

namespace KDL {
Expression<Wrench>::Ptr Wrench_VectorVector::derivativeExpression(int i)
{
    return wrench( argument1->derivativeExpression(i), argument2->derivativeExpression(i) );
}
Expression<Vector>::Ptr Force_Wrench::derivativeExpression(int i) {
    return force(argument->derivativeExpression(i));
}
Expression<Vector>::Ptr Torque_Wrench::derivativeExpression(int i) {
    return torque(argument->derivativeExpression(i));
}



Expression<Wrench>::Ptr Negate_Wrench::derivativeExpression(int i) {
    return -argument->derivativeExpression(i);
}
Expression<Wrench>::Ptr Addition_WrenchWrench::derivativeExpression(int i) {
    return argument1->derivativeExpression(i) + argument2->derivativeExpression(i);
}
Expression<Wrench>::Ptr Subtraction_WrenchWrench::derivativeExpression(int i) {
    return argument1->derivativeExpression(i) - argument2->derivativeExpression(i);
}

Expression<Wrench>::Ptr Composition_RotationWrench::derivativeExpression(int i) {
    Expression<Rotation>::Ptr a = cached<Rotation>(argument1);    
    Expression<Vector>::Ptr da  = cached<Vector>(argument1->derivativeExpression(i));    
    Expression<Wrench>::Ptr b   = cached<Wrench>(argument2);    
    Expression<Wrench>::Ptr db  = cached<Wrench>(argument2->derivativeExpression(i));    
    return wrench( a*force(db) + da*(a*force(b)),
                   a*torque(db) + da*(a*torque(b)) );
}

Expression<Wrench>::Ptr  Multiplication_WrenchDouble::derivativeExpression(int i) {
    return argument1*argument2->derivativeExpression(i) +
           argument1->derivativeExpression(i)*argument2;
}

Expression<Wrench>::Ptr RefPoint_WrenchVector::derivativeExpression(int i) {
    Expression<Wrench>::Ptr a   = cached<Wrench>(argument1);    
    Expression<Wrench>::Ptr da  = cached<Wrench>(argument1->derivativeExpression(i));    
    Expression<Vector>::Ptr b   = cached<Vector>(argument2);    
    Expression<Vector>::Ptr db  = cached<Vector>(argument2->derivativeExpression(i));    
        return wrench( force(da),
                       torque(da) + force(da)*b + force(a)*db );
}



} // end of namespace KDL

