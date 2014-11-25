/*
 * expressiontree_frame.cpp
 *
 *  Created on: September 
 *      Author: Wouter Bancken - Erwin Aertbelien
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

#include <kdl/expressiontree_double.hpp>
#include <kdl/expressiontree_vector.hpp>
#include <kdl/expressiontree_rotation.hpp>
#include <kdl/expressiontree_frame.hpp>
#include <kdl/expressiontree_twist.hpp>
#include <kdl/frames.hpp>
#include "util.hpp"

namespace KDL {
/*
Expression<Twist>::Ptr ChangeCoordinateFrame_FrameRotation::derivativeExpression(int i) {
}
*/
Expression<Twist>::Ptr Frame_RotationVector::derivativeExpression(int i) {
        int nr = getDep2<Rotation,Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant<Twist>(Twist::Zero());
        } if (nr==2) {
            return twist(
                argument2->derivativeExpression(i),
                Constant(Vector::Zero()) 
            );  
        } if (nr==3) {
            return twist(
                Constant(Vector::Zero()),
                argument1->derivativeExpression(i)
            );  
        } else {
            return twist(
                argument2->derivativeExpression(i),
                argument1->derivativeExpression(i)
            );  
        }
}


Expression<Twist>::Ptr Inverse_Frame::derivativeExpression(int i)
{
        int nr = getDep<Frame>(i,argument);
        if (nr==1) {
            return Constant(Twist::Zero());
        } else {
            Expression<Twist>::Ptr da = cached<Twist>( argument->derivativeExpression(i) ); 
            Expression<Frame>::Ptr val = cached<Frame>( argument  );
            Expression<Vector>::Ptr   rotvelda = cached<Vector>(rotvel(da));
            return inv(rotation(val))*twist( rotvelda*origin(val)-transvel(da), -rotvelda);
        }
}


Expression<Twist>::Ptr Composition_FrameFrame::derivativeExpression(int i)
{
        Expression<Frame>::Ptr a  = cached<Frame>(argument1 );
        Expression<Frame>::Ptr b  = cached<Frame>(argument2 );
        Expression<Twist>::Ptr da = cached<Twist>(argument1->derivativeExpression(i));
        Expression<Twist>::Ptr db = cached<Twist>(argument2->derivativeExpression(i));
        Expression<Vector>::Ptr rotvelda     = cached<Vector>(rotvel(da));
        Expression<Rotation>::Ptr rotationa  = cached<Rotation>(rotation(a));
        int nr = getDep<Frame>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Twist::Zero());
        } if (nr==2) {
            return twist(  rotationa*transvel(db) ,
                           rotationa*rotvel(db)
                   );
        } if (nr==3) {
            return twist(  rotvelda*(rotationa*origin(b)) + transvel(da),
                           rotvelda 
                   );
        } else {
            return twist(  rotvelda*(rotationa*origin(b)) + rotationa*transvel(db) + transvel(da),
                           rotvelda + rotationa*rotvel(db)
                   );
        }
}

Expression<Vector>::Ptr Composition_FrameVector::derivativeExpression(int i) 
{
        Expression<Rotation>::Ptr aRot = cached<Rotation>( rotation( argument1 ) );
        Expression<Twist>::Ptr da = cached<Twist>( argument1->derivativeExpression(i) );

        int nr = getDep2<Frame,Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Vector::Zero());
        } if (nr==2) {
            return aRot*argument2->derivativeExpression(i);
        } if (nr==3) {
            return rotvel(da)*( aRot* argument2 ) + transvel(da);
        } else {
            return rotvel(da)*( aRot* argument2 ) +
                   aRot*argument2->derivativeExpression(i) +
                   transvel(da);
        }
}


Expression<Vector>::Ptr Origin_Frame::derivativeExpression(int i)
{
        int nr = getDep<Frame>(i,argument);
        if (nr==1) {
            return Constant(Vector::Zero());
        } else {
            return transvel( argument->derivativeExpression(i) );
        }
}


Expression<Vector>::Ptr Rotation_Frame::derivativeExpression(int i)
{
        int nr = getDep<Frame>(i,argument);
        if (nr==1) {
            return Constant(Vector::Zero());
        } else {
            return rotvel( argument->derivativeExpression(i) );
        }
}






}; // namespace KDL
