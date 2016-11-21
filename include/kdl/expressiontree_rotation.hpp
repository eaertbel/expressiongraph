/*
 * expressiontree_rotation.cpp
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

#ifndef KDL_EXPRESSIONTREE_ROTATION_HPP
#define KDL_EXPRESSIONTREE_ROTATION_HPP

#include <kdl/expressiontree_expressions.hpp>
#include <kdl/expressiontree_vector.hpp>
#include <kdl/expressiontree_double.hpp>
#include <kdl/frames.hpp>


namespace KDL {
/*
 * Rotation
 *      - RotX(Expression<double>::Ptr)          					returns Expression<Rotation>::Ptr
 *      - RotY(Expression<double>::Ptr)          					returns Expression<Rotation>::Ptr
 *      - RotZ(Expression<double>::Ptr)          					returns Expression<Rotation>::Ptr
 *      - Inverse(Expression<Rotation>::Ptr)     					returns Expression<Rotation>::Ptr
 *      - Expression<Rotation>::Ptr*Expression<Rotation::Ptr> 		returns Expression<Rotation>::Ptr
 *      - Expression<Rotation>::Ptr*Expression<Vector::Ptr>			returns Expression<Vector>::Ptr
 *      - UnitX( Expression<Rotation>::Ptr)      					returns Expression<Vector>::Ptr
 *      - UnitY( Expression<Rotation>::Ptr)      					returns Expression<Vector>::Ptr
 *      - UnitZ( Expression<Rotation>::Ptr)      					returns Expression<Vector>::Ptr
 *      - construct_rotation_from_vectors                           returns Expression<Rotation>::Ptr
 */
// Rot Double
class Rot_Double:
    public UnaryExpression<KDL::Rotation, double>
{
public:
    typedef UnaryExpression<KDL::Rotation, double> UnExpr;
    KDL::Vector axis;
public:
    Rot_Double() {}
    Rot_Double( const KDL::Vector& _axis,
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("rot",arg),
                axis(_axis/_axis.Norm())
                {}

    virtual KDL::Rotation value() {
        return KDL::Rotation::Rot2(axis,argument->value());
    }

    virtual KDL::Vector derivative(int i) {
        return axis*argument->derivative(i);
    }

    virtual void print(std::ostream& os) const {
        os << "rot( Vector("<<axis[0]<<","<<axis[1] <<"," <<axis[2] << "),";
        argument->print(os);
        os << ")";  
    }
    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new Rot_Double(axis, argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr rot(const KDL::Vector& axis, Expression<double>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new Rot_Double(axis, a )
    );
    return expr;
}

/**
 *
 * Rot Double
 *   argument1: rotation vector (should be normalized!)
 *   argument2: rotation angle
 * result:
 *   Rotation matrix corresponding to rotating with angle around the given rotation vector.
 */

class RotVec_Double:
    public BinaryExpression<KDL::Rotation, KDL::Vector, double>
{
public:
    typedef BinaryExpression<KDL::Rotation, KDL::Vector, double> BinExpr;
    KDL::Vector axis_value;
    double      angle_value;
public:
    RotVec_Double() {}
    RotVec_Double( 
                const  BinExpr::Argument1Expr::Ptr& arg1,
                const  BinExpr::Argument2Expr::Ptr& arg2
                ):
        BinExpr("rotVec",arg1,arg2)
        {}

    virtual KDL::Rotation value() {
        axis_value = argument1->value();
        angle_value= argument2->value();
        return KDL::Rotation::Rot2(axis_value,angle_value);
    }

    virtual KDL::Vector derivative(int i) {
        return axis_value*argument2->derivative(i) + argument1->derivative(i)*angle_value;
    }
    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new RotVec_Double(argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr rotVec(Expression<Vector>::Ptr a, Expression<double>::Ptr b) {
    Expression<KDL::Rotation>::Ptr expr(
        new RotVec_Double(a, b )
    );
    return expr;
}



// RotX Double
class RotX_Double:
    public UnaryExpression<KDL::Rotation, double>
{
public:
    typedef UnaryExpression<KDL::Rotation, double> UnExpr;
public:
    RotX_Double() {}
    RotX_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("rot_x",arg)
                {}
    virtual KDL::Rotation value() {
        return KDL::Rotation::RotX(argument->value());
    }

    virtual KDL::Vector derivative(int i) {
        return KDL::Vector(argument->derivative(i),0,0);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new RotX_Double( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr rot_x( Expression<double>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new RotX_Double( a )
    );
    return expr;
}

// RotY Double
class RotY_Double:
    public UnaryExpression<KDL::Rotation, double>
{
public:
    typedef UnaryExpression<KDL::Rotation, double> UnExpr;
public:
    RotY_Double() {}
    RotY_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("rot_y",arg)
                {}
    virtual KDL::Rotation value() {
        return KDL::Rotation::RotY(argument->value());
    }

    virtual KDL::Vector derivative(int i) {
        return KDL::Vector(0,argument->derivative(i),0);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new RotY_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr rot_y( Expression<double>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new RotY_Double( a )
    );
    return expr;
}

// RotZ Double
class RotZ_Double:
    public UnaryExpression<KDL::Rotation, double>
{
public:
    typedef UnaryExpression<KDL::Rotation, double> UnExpr;
public:
    RotZ_Double() {}
    RotZ_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("rot_z",arg)
                {}
    virtual KDL::Rotation value() {
        return KDL::Rotation::RotZ(argument->value());
    }

    virtual KDL::Vector derivative(int i) {
        return KDL::Vector(0,0,argument->derivative(i));
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new RotZ_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr rot_z( Expression<double>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new RotZ_Double(a)
    );
    return expr;
}

//Inverse Rotation
class Inverse_Rotation:
    public UnaryExpression<KDL::Rotation, KDL::Rotation>
{
public:
    typedef UnaryExpression<KDL::Rotation, KDL::Rotation> UnExpr;
    KDL::Rotation val;
public:
    Inverse_Rotation() {}
    Inverse_Rotation(
            const  UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("inv",arg)
            {}

    virtual KDL::Rotation value() {
        val = argument->value();
        return val.Inverse();
    }

    virtual KDL::Vector derivative(int i) {
        return val.Inverse(-(argument->derivative(i)));
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new Inverse_Rotation( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr inv (Expression<KDL::Rotation>::Ptr a) {
    Expression<KDL::Rotation>::Ptr expr(
        new Inverse_Rotation(a)
    );
    return expr;
}

//Composition Rotation Rotation
class Composition_RotationRotation:
    public BinaryExpression<KDL::Rotation, KDL::Rotation, KDL::Rotation>
{
public:
    typedef BinaryExpression<KDL::Rotation, KDL::Rotation, KDL::Rotation> BinExpr;
    KDL::Rotation arg1value;
public:
    Composition_RotationRotation() {}
    Composition_RotationRotation(
                    const  BinExpr::Argument1Expr::Ptr& arg1,
                    const  BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("compose",arg1,arg2)
                    {}

    virtual KDL::Rotation value() {
        arg1value = argument1->value();
        return arg1value*argument2->value();
    }

    virtual KDL::Vector derivative(int i) {
        return arg1value*argument2->derivative(i) + argument1->derivative(i);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Rotation>::Ptr expr(
            new Composition_RotationRotation( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr operator* ( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Rotation>::Ptr a2 ) {
    Expression<KDL::Rotation>::Ptr expr(new Composition_RotationRotation( a1, a2 ));
    return expr;
}

//Composition Rotation Vector
class Composition_RotationVector:
	public BinaryExpression<KDL::Vector, KDL::Rotation, KDL::Vector>
{
public:
	typedef BinaryExpression<KDL::Vector,KDL::Rotation,KDL::Vector> BinExpr;
	KDL::Rotation arg1value;
	KDL::Vector arg2value;
public:
	Composition_RotationVector() {}
	Composition_RotationVector(
			const   BinExpr::Argument1Expr::Ptr& arg1,
			const   BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("transform",arg1,arg2)
				{}

	virtual KDL::Vector value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Vector derivative(int i){
		return argument1->derivative(i) * (arg1value * arg2value) + arg1value*argument2->derivative(i);
	}

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Composition_RotationVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr operator* ( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Vector>::Ptr a2 ) {
	Expression<KDL::Vector>::Ptr expr(new Composition_RotationVector(a1, a2));
	return expr;
}

//UnitX Rotation
class UnitX_Rotation:
    public UnaryExpression<KDL::Vector, KDL::Rotation>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Rotation> UnExpr;
    KDL::Vector val;
public:
    UnitX_Rotation() {}
    UnitX_Rotation(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("unit_x",arg)
                {}

    virtual KDL::Vector value() {
    	val = argument->value().UnitX();
    	return val;
    }

    virtual KDL::Vector derivative(int i) {
    	return argument->derivative(i) * val;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new UnitX_Rotation( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr unit_x ( Expression<KDL::Rotation>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new UnitX_Rotation( a )
    );
    return expr;
}

//UnitY Rotation
class UnitY_Rotation:
    public UnaryExpression<KDL::Vector, KDL::Rotation>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Rotation> UnExpr;
    KDL::Vector val;
public:
    UnitY_Rotation() {}
    UnitY_Rotation(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("unit_y",arg)
                {}

    virtual KDL::Vector value() {
    	val = argument->value().UnitY();
    	return val;
    }

    virtual KDL::Vector derivative(int i) {
    	return argument->derivative(i) * val;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new UnitY_Rotation( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr unit_y ( Expression<KDL::Rotation>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new UnitY_Rotation( a )
    );
    return expr;
}

//UnitZ Rotation
class UnitZ_Rotation:
    public UnaryExpression<KDL::Vector, KDL::Rotation>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Rotation> UnExpr;
    KDL::Vector val;
public:
    UnitZ_Rotation() {}
    UnitZ_Rotation(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("unit_z",arg)
                {}

    virtual KDL::Vector value() {
    	val = argument->value().UnitZ();
    	return val;
    }

    virtual KDL::Vector derivative(int i) {
    	return argument->derivative(i) * val;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new UnitZ_Rotation( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr unit_z ( Expression<KDL::Rotation>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new UnitZ_Rotation( a )
    );
    return expr;
}

class Construct_Rotation:
    public TernaryExpression<KDL::Rotation, KDL::Vector, KDL::Vector, KDL::Vector>
{
public:
    typedef TernaryExpression<KDL::Rotation, KDL::Vector, KDL::Vector, KDL::Vector> TExpr;
public:
    Construct_Rotation() {}  
    Construct_Rotation( const TExpr::Argument1Expr::Ptr& arg1,
                        const TExpr::Argument2Expr::Ptr& arg2,
                        const TExpr::Argument3Expr::Ptr& arg3) : TExpr("construct_rotation",arg1,arg2,arg3) {}  
    virtual KDL::Rotation value() {
        return KDL::Rotation(argument1->value(), argument2->value(), argument3->value());
    }
    virtual KDL::Vector derivative(int i) {
        KDL::Rotation R(argument1->value(), argument2->value(), argument3->value());
        KDL::Rotation Rd(argument1->derivative(i), argument2->derivative(i), argument3->derivative(i));
        KDL::Rotation omegax = Rd*R.Inverse();
        KDL::Vector tmp(  (omegax(2,1)-omegax(1,2))/2.0,  ( omegax(0,2)-omegax(2,0))/2.0, (omegax(1,0)-omegax(0,1))/2.0 );  
        return tmp;
    } 
    virtual Expression<Vector>::Ptr derivativeExpression(int i) {
        assert( 0 /*not yet implemented */ );
    }
    virtual TExpr::Ptr clone() {
        TExpr::Ptr expr( new Construct_Rotation( argument1->clone(), argument2->clone(), argument3->clone()));
        return expr;
    }
};

inline Expression<KDL::Rotation>::Ptr construct_rotation_from_vectors( Expression<KDL::Vector>::Ptr a, Expression<KDL::Vector>::Ptr b, Expression<KDL::Vector>::Ptr c) {
    Expression<KDL::Rotation>::Ptr expr(
        new Construct_Rotation(a,b,c) 
    );
    return expr;
}

class Get_Rotation_Vector:
    public UnaryExpression<KDL::Vector, KDL::Rotation>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Rotation> UnExpr;
    KDL::Vector val;
public:
    Get_Rotation_Vector() {}
    Get_Rotation_Vector(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("getRotVec",arg)
                {}

    virtual KDL::Vector value() {
    	val = argument->value().GetRot();
    	return val;
    }

    virtual KDL::Vector derivative(int i) {
    	return argument->derivative(i);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Get_Rotation_Vector( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr getRotVec( Expression<KDL::Rotation>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Get_Rotation_Vector( a )
    );
    return expr;
}

class Get_RPY_Rotation:
    public UnaryExpression<KDL::Vector, KDL::Rotation>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Rotation> UnExpr;
    KDL::Vector val;
    double m00,m01,m02,m10,m11,m12,m20,m21,m22;
public:
    Get_RPY_Rotation() {}
    Get_RPY_Rotation(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("getRPY",arg)
                {}

    virtual KDL::Vector value() {
    	argument->value().GetRPY(val[0],val[1],val[2]);
        double ca = cos(val[2]);
        double sa = sin(val[2]);
        double cb = cos(val[1]);
        double sb = sin(val[1]);
        m00 = ca/cb;   m01 = sa/cb;    m02 = 0.0;
        m10 = -sa;     m11 = ca;       m12 = 0.0;
        m20 = ca*sb/cb;m21 = sa*sb/cb; m22 = 1.0; 
    	return val;
    }

    virtual KDL::Vector derivative(int i) {
        Vector omega = argument->derivative(i);
        Vector result;
        result[0] = m00*omega[0] + m01*omega[1] + m02*omega[2];
        result[1] = m10*omega[0] + m11*omega[1] + m12*omega[2];
        result[2] = m20*omega[0] + m21*omega[1] + m22*omega[2];
    	return result;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Get_RPY_Rotation( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr getRPY( Expression<KDL::Rotation>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Get_RPY_Rotation( a )
    );
    return expr;
}



} // end of namespace KDL
#endif
