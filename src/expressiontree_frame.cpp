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

#include <expressiongraph/expressiontree_double.hpp>
#include <expressiongraph/expressiontree_vector.hpp>
#include <expressiongraph/expressiontree_rotation.hpp>
#include <expressiongraph/expressiontree_frame.hpp>
#include <expressiongraph/expressiontree_twist.hpp>
#include <kdl/frames.hpp>
#include "util.hpp"

namespace KDL {


/*
class ChangeCoordinateFrame_FrameRotation:
    public BinaryExpression<Frame,Frame,Rotation>
{
public:

    ChangeCoordinateFrame_FrameRotation(  Expression<Frame>::Ptr a1, 
                           Expression<Rotation>::Ptr a2):
        BinaryExpression<Frame,Rotation,Vector>("change coordinate frame",a1,a2) {}

    virtual Frame value() {
        
        return argument1->value().ChangeCoordinateFrame(argument2->value());
    } 

    virtual Twist derivative(int i) {
        xxx
        return Twist(
            argument2->derivative(i),
            argument1->derivative(i)
        );
    } 

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual typename Expression<Frame>::Ptr clone() {
        typename Expression<Frame>::Ptr expr(
            new ChangeCoordinateFrame_FrameRotation(argument1->clone(), argument2->clone())
        );
        return expr;
    } 
};

inline Expression<Frame>::Ptr change_coordinate_frame( Expression<Frame>::Ptr a1, 
                                Expression<Rotation>::Ptr a2) {
    Expression<Frame>::Ptr expr(
        new ChangeCoordinateFrame_FrameRotation( a1, a2 )
    );
    return expr;
}
*/
class Frame_RotationVector:
    public BinaryExpression<Frame,Rotation,Vector>
{
public:

    Frame_RotationVector(){}
    Frame_RotationVector(  Expression<Rotation>::Ptr a1, 
                           Expression<Vector>::Ptr a2):
        BinaryExpression<Frame,Rotation,Vector>("frame",a1,a2) {}

    virtual Frame value() {
        return Frame(argument1->value(),argument2->value());
    } 

    virtual Twist derivative(int i) {
        return Twist(
            argument2->derivative(i),
            argument1->derivative(i)
        );
    } 

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual Expression<Frame>::Ptr clone() {
         Expression<Frame>::Ptr expr(
            new Frame_RotationVector(argument1->clone(), argument2->clone())
        );
        return expr;
    } 
};

Expression<Frame>::Ptr frame( Expression<Rotation>::Ptr a1, 
                                Expression<Vector>::Ptr a2) {
    Expression<Frame>::Ptr expr(
        new Frame_RotationVector( a1, a2 )
    );
    return expr;
}

Expression<Frame>::Ptr frame( Expression<Rotation>::Ptr a1 ) {
    Expression<Frame>::Ptr expr(
        new Frame_RotationVector( a1, Constant(Vector::Zero() ))
    );
    return expr;
}

Expression<Frame>::Ptr frame(  Expression<Vector>::Ptr a2) {
    Expression<Frame>::Ptr expr(
        new Frame_RotationVector( Constant(Rotation::Identity()), a2 )
    );
    return expr;
}


//Inverse Frame
class Inverse_Frame:
    public UnaryExpression<KDL::Frame, KDL::Frame>
{
public:
    typedef UnaryExpression<KDL::Frame, KDL::Frame> UnExpr;
    KDL::Frame val;
public:
    Inverse_Frame() {}
    Inverse_Frame(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("inv",arg)
                {}

    virtual KDL::Frame value() {
        val = argument->value();
        return val.Inverse();
    }

    virtual KDL::Twist derivative(int i) {
    	KDL::Twist da = argument->derivative(i);
        return KDL::Twist(val.M.Inverse(da.rot*val.p-da.vel), -val.M.Inverse(da.rot));
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Frame>::Ptr expr(
            new Inverse_Frame( argument->clone())
        );
        return expr;
    }
};

Expression<KDL::Frame>::Ptr inv ( Expression<KDL::Frame>::Ptr a) {
    Expression<KDL::Frame>::Ptr expr(
        new Inverse_Frame(a)
    );
    return expr;
}

//Composition Frame Frame
class Composition_FrameFrame:
	public BinaryExpression<KDL::Frame, KDL::Frame, KDL::Frame>
{
public:
	typedef BinaryExpression<KDL::Frame,KDL::Frame,KDL::Frame> BinExpr;
	KDL::Frame arg1value;
	KDL::Frame arg2value;
public:
	Composition_FrameFrame() {}
	Composition_FrameFrame(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("compose",arg1,arg2)
				{}

	virtual KDL::Frame value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Twist derivative(int i){
		KDL::Twist da = argument1->derivative(i);
		KDL::Twist db = argument2->derivative(i);
		return KDL::Twist(
					da.rot*(arg1value.M*arg2value.p)+arg1value.M * db.vel + da.vel,
					da.rot + arg1value.M * db.rot
		);
	}

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Frame>::Ptr expr(
            new Composition_FrameFrame( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Frame>::Ptr operator* ( Expression<KDL::Frame>::Ptr a1, Expression<KDL::Frame>::Ptr a2 ) {
	Expression<KDL::Frame>::Ptr expr(new Composition_FrameFrame( a1, a2 ));
	return expr;
}

//Composition Frame Vector
class Composition_FrameVector:
	public BinaryExpression<KDL::Vector, KDL::Frame, KDL::Vector>
{
public:
	typedef BinaryExpression<KDL::Vector,KDL::Frame,KDL::Vector> BinExpr;
	KDL::Frame arg1value;
	KDL::Vector arg2value;
public:
	Composition_FrameVector() {}
	Composition_FrameVector(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("transform",arg1,arg2)
				{}

	virtual KDL::Vector value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Vector derivative(int i){
		KDL::Twist da = argument1->derivative(i);
		return da.rot*(arg1value.M*arg2value)+arg1value.M*argument2->derivative(i)+da.vel;
	}

    virtual Expression<Vector>::Ptr derivativeExpression(int i);
    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Composition_FrameVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Vector>::Ptr operator* ( Expression<KDL::Frame>::Ptr a1, Expression<KDL::Vector>::Ptr a2 ) {
	Expression<KDL::Vector>::Ptr expr(new Composition_FrameVector( a1, a2 ));
	return expr;
}

//Origin Frame
class Origin_Frame:
    public UnaryExpression<KDL::Vector, KDL::Frame>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Frame> UnExpr;
    Origin_Frame() {}
    Origin_Frame(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("origin",arg)
                {}

    virtual KDL::Vector value() {
    	return argument->value().p;
    }

    virtual KDL::Vector derivative(int i) {
    	return argument->derivative(i).vel;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Origin_Frame( argument->clone())
        );
        return expr;
    }
};

Expression<KDL::Vector>::Ptr origin ( Expression<KDL::Frame>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Origin_Frame( a )
    );
    return expr;
}

//Rotation Frame
class Rotation_Frame:
    public UnaryExpression<KDL::Rotation, KDL::Frame>
{
public:
    typedef UnaryExpression<KDL::Rotation, KDL::Frame> UnExpr;
    Rotation_Frame() {}
    Rotation_Frame(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("rotation",arg)
                {}

    virtual KDL::Rotation value() {
    	return argument->value().M;
    }

    virtual KDL::Vector derivative(int i) {
    	return argument->derivative(i).rot;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new Rotation_Frame( argument->clone())
        );
        return expr;
    }
};

Expression<KDL::Rotation>::Ptr rotation( Expression<KDL::Frame>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new Rotation_Frame(a)
    );
    return expr;
}

// just because the analog exists for rotation matrices, we define it
// also for rotation matrices:
Expression<Vector>::Ptr apply ( 
        Expression<Frame>::Ptr F, Expression<Vector>::Ptr v) {
    return F*v; 
}









































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
