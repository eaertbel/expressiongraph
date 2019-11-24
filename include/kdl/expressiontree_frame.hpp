/*
 * expressiontree_frame.hpp
 *
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

#ifndef KDL_EXPRESSIONTREE_FRAME_HPP
#define KDL_EXPRESSIONTREE_FRAME_HPP

#include <kdl/expressiontree_expressions.hpp>
#include <kdl/frames.hpp>
/*
 * \todo : complete documentation
 * Frame related expression tree classes.
 * Operations have the same semantics as the operations in frames.hpp,
 * (naming is similar, except for member functions that are replaced by functions) 
 *
 * Frame unary operations:
 *     - inv
 * Frame binary operations:
 *     - "*"
 *     -
 * Corresponding Frame unary Expression classes:
 *     -  Inverse_Frame
 * Corresponding Frame binary Expression classes:
 *     - Composition_FrameFrame
 *     - Composition_FrameVector
 *
 * 
 *      - Inverse(Jacobian<Frame>)        			returns Jacobian<Frame>
 *      - Jacobian<Frame> * Jacobian<Frame>     	returns Jacobian<Frame>
 *      - Jacobian<Frame> * Jacobian<Vector>    	returns Jacobian<Vector>
 *      - Origin( Jacobian<Frame> )       			returns Jacobian<Vector>
 *      - RotMat( Jacobian<Frame> )      			returns Jacobian<Rotation>
 *      change_coordinate_frame ...
 * Constructor:
 *    Frame_RotationVector -> frame(Expression<Rotation>::Ptr, Expression<Vector>::Ptr )
 *                         -> frame( Expression<Vector>::Ptr )
 *                         -> frame( Expression<Rotation>::Ptr )
 *                     
 */

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

inline Expression<Frame>::Ptr frame( Expression<Rotation>::Ptr a1, 
                                Expression<Vector>::Ptr a2) {
    Expression<Frame>::Ptr expr(
        new Frame_RotationVector( a1, a2 )
    );
    return expr;
}

inline Expression<Frame>::Ptr frame( Expression<Rotation>::Ptr a1 ) {
    Expression<Frame>::Ptr expr(
        new Frame_RotationVector( a1, Constant(Vector::Zero() ))
    );
    return expr;
}

inline Expression<Frame>::Ptr frame(  Expression<Vector>::Ptr a2) {
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

inline Expression<KDL::Frame>::Ptr inv ( Expression<KDL::Frame>::Ptr a) {
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

inline Expression<KDL::Frame>::Ptr operator* ( Expression<KDL::Frame>::Ptr a1, Expression<KDL::Frame>::Ptr a2 ) {
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

inline Expression<KDL::Vector>::Ptr operator* ( Expression<KDL::Frame>::Ptr a1, Expression<KDL::Vector>::Ptr a2 ) {
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

inline Expression<KDL::Vector>::Ptr origin ( Expression<KDL::Frame>::Ptr a) {
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

inline Expression<KDL::Rotation>::Ptr rotation( Expression<KDL::Frame>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new Rotation_Frame(a)
    );
    return expr;
}

// just because the analog exists for rotation matrices, we define it
// also for rotation matrices:
inline Expression<Vector>::Ptr apply ( 
        Expression<Frame>::Ptr F, Expression<Vector>::Ptr v) {
    return F*v; 
}

}; // namespace KDL
#endif
