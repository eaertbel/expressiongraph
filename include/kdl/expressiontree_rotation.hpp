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
#include <kdl/frames.hpp>


namespace KDL {
/*
 * Rotation
 *      - RotX(Jacobian<double>)          					returns Jacobian<Rotation>
 *      - RotY(Jacobian<double>)          					returns Jacobian<Rotation>
 *      - RotZ(Jacobian<double>)          					returns Jacobian<Rotation>
 *      - Inverse(Jacobian<Rotation>)     					returns Jacobian<Rotation>
 *      - Jacobian<Rotation>*Jacobian<Rotation> 			returns Jacobian<Rotation>
 *      - Jacobian<Rotation>*Jacobian<Vector>				returns Jacobian<Vector>
 *      - UnitX( Jacobian<Rotation>)      					returns Jacobian<Vector>
 *      - UnitY( Jacobian<Rotation>)      					returns Jacobian<Vector>
 *      - UnitZ( Jacobian<Rotation>)      					returns Jacobian<Vector>
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

} // end of namespace KDL
#endif
