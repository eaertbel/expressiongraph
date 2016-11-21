/*
 * expressiontree_rotation.cpp
 *
 *  Created on: Sept., 2012
 *      Author: Erwin Aertbelien
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
#include <kdl/frames.hpp>
#include "util.hpp"

namespace KDL {

Expression<Vector>::Ptr RotVec_Double::derivativeExpression(int i) {
        Expression<Vector>::Ptr arg1 = cached<Vector>( argument1 );
        Expression<double>::Ptr   arg2 = cached<double>( argument2 );
        int nr = getDep2<Vector,double>(i,argument1, argument2);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } 
        if (nr==2) {
            return argument1*argument2->derivativeExpression(i);
        }
        if (nr==3) {
            return argument1->derivativeExpression(i)*argument2;
        } else {
            return argument1*argument2->derivativeExpression(i) + argument1->derivativeExpression(i)*argument2;
        }
}



Expression<Vector>::Ptr Rot_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return Constant( axis ) * argument->derivativeExpression(i);
        }
}

Expression<Vector>::Ptr RotX_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return KDL::vector( argument->derivativeExpression(i), Constant(0.0), Constant(0.0) );
        }
}

Expression<Vector>::Ptr RotY_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return KDL::vector( Constant(0.0), argument->derivativeExpression(i), Constant(0.0) );
        }
}

Expression<Vector>::Ptr RotZ_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return KDL::vector( Constant(0.0), Constant(0.0), argument->derivativeExpression(i)  );
        }
}

Expression<Vector>::Ptr Inverse_Rotation::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return -( inv( argument ) * argument->derivativeExpression(i) );
        }
}

Expression<Vector>::Ptr Composition_RotationRotation::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument1,argument2);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } if (nr==2) {
            return argument1*argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i);
        } else {
            return argument1*argument2->derivativeExpression(i) + argument1->derivativeExpression(i);
        }
}

Expression<Vector>::Ptr Composition_RotationVector::derivativeExpression(int i) {
        Expression<Rotation>::Ptr arg1 = cached<Rotation>( argument1 );
        Expression<Vector>::Ptr arg2 = cached<Vector>( argument2 );
        int nr = getDep2<Rotation,Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Vector::Zero());
        } if (nr==2) {
           return arg1*argument2->derivativeExpression(i);
        } if (nr==3) {
           return argument1->derivativeExpression(i) * (arg1 * arg2);
        } else {
           return argument1->derivativeExpression(i) * (arg1 * arg2) + arg1*argument2->derivativeExpression(i);
        }
}

Expression<Vector>::Ptr UnitX_Rotation::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return argument->derivativeExpression(i) * unit_x( argument );
        }
}
Expression<Vector>::Ptr UnitY_Rotation::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return argument->derivativeExpression(i) * unit_y( argument );
        }
}
Expression<Vector>::Ptr UnitZ_Rotation::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return argument->derivativeExpression(i) * unit_z( argument );
        }
}

Expression<Vector>::Ptr Get_Rotation_Vector::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return argument->derivativeExpression(i);
        }
}

Expression<Vector>::Ptr Get_RPY_Rotation::derivativeExpression(int i) {
        int nr = getDep<Rotation>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            Expression<Vector>::Ptr omega = argument->derivativeExpression(i);
            Expression<Vector>::Ptr rpy   = getRPY( argument);
            Expression<double>::Ptr ca    = cos( coord_z(rpy) );
            Expression<double>::Ptr sa    = sin( coord_z(rpy) );
            Expression<double>::Ptr cb    = cos( coord_y(rpy) );
            Expression<double>::Ptr sb    = sin( coord_y(rpy) );
            return vector(
                    ca/cb*coord_x(omega) +    sa/cb*coord_y(omega),
                      -sa*coord_x(omega) +       ca*coord_y(omega),
                 ca*sb/cb*coord_x(omega) + sa*sb/cb*coord_y(omega) + coord_z(omega)
            );
        }
}





} // end of namespace KDL
