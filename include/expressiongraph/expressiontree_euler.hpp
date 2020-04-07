/*
 * expressiontree_euler.cpp
 *
 *  Created on: Aug 7, 2012
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

#ifndef KDL_EXPRESSIONTREE_EULER_HPP
#define KDL_EXPRESSIONTREE_EULER_HPP

#include "expressiontree_expressions.hpp"
#include <kdl/expressiontree_rotation.hpp>



namespace KDL {

class EulerZYX_double_double_double:
    public TernaryExpression<KDL::Rotation, double, double, double>
{
    double ca,sa,cb,sb,cc,sc;
public:
    typedef TernaryExpression<KDL::Rotation, double,double,double> TernExpr;
public:
    EulerZYX() {}
    EulerZYX( const typename Expression<double>::Ptr& arg1,
              const typename Expression<double>::Ptr& arg2,
              const typename Expression<double>::Ptr& arg3):
                TernExpr("eulerzyx",arg1,arg2,arg3)
                {}

    virtual KDL::Rotation value() {
        double arg;
        arg = argument1->value();
        ca = cos(arg);sa=sin(arg); 
        arg = argument2->value();
        cb = cos(arg);sb=sin(arg); 
        arg = argument3->value();
        cc = cos(arg);sc=sin(arg); 
        return KDL::Rotation( ca*cb, ca*sb*sc-cc*sa, sa*sc+ca*cc*sb,
                              cb*sa, ca*cc+sa*sb*sc, cc*sa*sb-ca*sc,
                              -sb  ,  cb*sc,         cb*cc );
    }

    virtual KDL::Vector derivative(int i) {
        double da = argument1->derivative(i);
        double db = argument2->derivative(i);
        double dc = argument3->derivative(i);
        return KDL::Vector( -sa*db + ca*cb*dc ,
                             ca*db + cb*sa*dc ,
                             da - sb*dc );
        // vector(0,0,dc) 
        //   + Rot_Z(a)*vector(0,db,0) 
        //   + Rot_Z(a)*Rot_Y(b)*vector(dc,0,0)
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual typename TernExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new EulerZYX(argument1->clone(),argument2->clone(),argument3->clone())
        );
        return expr;
    }
};



} // namespace KDL
#endif
