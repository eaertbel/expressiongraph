/*
 * expressiontree_twist.hpp
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

#ifndef KDL_EXPRESSIONTREE_TWIST_HPP
#define KDL_EXPRESSIONTREE_TWIST_HPP

#include <kdl/expressiontree_expressions.hpp>
#include <kdl/frames.hpp>

namespace KDL {
/*
 * Twist operations:
 *     -Twist                	returns Twist
 *     Twist+Twist           	returns Twist
 *     Twist-Twist           	returns Twist
 *     Rotation*Twist        	returns Twist
 *     Twist*double          	returns Twist
 *     double*Twist          	returns Twist
 *     RefPoint(Twist,Vector) 	returns Twist
 *     Stiffness*Twist			returns Wrench
 *     vel() returns a Vector corresponding to translational velocity
 *     rot() returns a Vector corresponding to rotational velocity
 */


//Negate Twist
class Twist_VectorVector:
    public BinaryExpression<Twist, Vector, Vector>
{
public:
    Twist_VectorVector() {}
    Twist_VectorVector( const Expression<Vector>::Ptr& arg1,
                  const Expression<Vector>::Ptr& arg2  ):
                BinaryExpression<Twist,Vector,Vector>("twist",arg1,arg2)
                {}

    virtual KDL::Twist value() {
        return Twist(argument1->value(),argument2->value());
    }

    virtual KDL::Twist derivative(int i) {
        return Twist(argument1->derivative(i),argument2->derivative(i));
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual Expression<Twist>::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Twist_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Twist>::Ptr twist( Expression<KDL::Vector>::Ptr a, Expression<KDL::Vector>::Ptr b) {
    Expression<KDL::Twist>::Ptr expr(
        new Twist_VectorVector( a,b )
    );
    return expr;
}



//Negate Twist
class Negate_Twist:
    public UnaryExpression<KDL::Twist, KDL::Twist>
{
public:
    typedef UnaryExpression<KDL::Twist, KDL::Twist> UnExpr;
public:
    Negate_Twist() {}
    Negate_Twist(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("negation",arg)
                {}

    virtual KDL::Twist value() {
    	return -argument->value();
    }

    virtual KDL::Twist derivative(int i) {
    	return -argument->derivative(i);
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Negate_Twist( argument->clone())
        );
        return expr;
    }
};
inline Expression<KDL::Twist>::Ptr operator-( Expression<KDL::Twist>::Ptr a) {
    Expression<KDL::Twist>::Ptr expr(
        new Negate_Twist( a )
    );
    return expr;
}




class Velocity_Twist:
    public UnaryExpression<Vector, Twist>
{
public:
    typedef UnaryExpression<Vector, Twist> UnExpr;
public:
    Velocity_Twist() {}
    Velocity_Twist(const Expression<Twist>::Ptr& arg):
                UnExpr("transvel",arg)
                {}

    virtual Vector value() {
    	return argument->value().vel;
    }

    virtual Vector derivative(int i) {
    	return argument->derivative(i).vel;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<Vector>::Ptr expr(
            new Velocity_Twist( argument->clone())
        );
        return expr;
    }
};

inline Expression<Vector>::Ptr transvel( Expression<Twist>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Velocity_Twist( a )
    );
    return expr;
}

class RotVelocity_Twist:
    public UnaryExpression<Vector, Twist>
{
public:
    typedef UnaryExpression<Vector, Twist> UnExpr;
    RotVelocity_Twist() {}
    RotVelocity_Twist(
                  Expression<Twist>::Ptr arg):
                UnExpr("rotvel",arg)
                {}

    virtual Vector value() {
    	return argument->value().rot;
    }

    virtual Vector derivative(int i) {
    	return argument->derivative(i).rot;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<Vector>::Ptr expr(
            new RotVelocity_Twist( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr rotvel( Expression<KDL::Twist>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new RotVelocity_Twist( a )
    );
    return expr;
}


//Addition Twist Twist
class Addition_TwistTwist:
    public BinaryExpression<KDL::Twist, KDL::Twist, KDL::Twist>
{
public:
    typedef BinaryExpression<KDL::Twist, KDL::Twist, KDL::Twist> BinExpr;
    Addition_TwistTwist() {}
    Addition_TwistTwist(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("addition",arg1,arg2)
                    {}

    virtual KDL::Twist value() {
        return argument1->value() + argument2->value();
    }

    virtual KDL::Twist derivative(int i) {
        return argument1->derivative(i) + argument2->derivative(i);
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Addition_TwistTwist( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Twist>::Ptr operator+( Expression<KDL::Twist>::Ptr a1, Expression<KDL::Twist>::Ptr a2) {
    Expression<KDL::Twist>::Ptr expr(
        new Addition_TwistTwist( a1, a2 )
    );
    return expr;
}

//Subtraction Twist Twist
class Subtraction_TwistTwist:
    public BinaryExpression<KDL::Twist, KDL::Twist, KDL::Twist>
{
public:
    typedef BinaryExpression<KDL::Twist, KDL::Twist, KDL::Twist> BinExpr;
    Subtraction_TwistTwist() {}
    Subtraction_TwistTwist(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("subtraction",arg1,arg2)
                    {}

    virtual KDL::Twist value() {
        return argument1->value() - argument2->value();
    }

    virtual KDL::Twist derivative(int i) {
        return argument1->derivative(i) - argument2->derivative(i);
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Subtraction_TwistTwist( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Twist>::Ptr operator-( Expression<KDL::Twist>::Ptr a1, Expression<KDL::Twist>::Ptr a2) {
    Expression<KDL::Twist>::Ptr expr(
        new Subtraction_TwistTwist( a1, a2 )
    );
    return expr;
}

//Composition Rotation Twist
class Composition_RotationTwist:
    public BinaryExpression<KDL::Twist, KDL::Rotation, KDL::Twist>
{
public:
    typedef BinaryExpression<KDL::Twist, KDL::Rotation, KDL::Twist> BinExpr;
    KDL::Rotation arg1value;
    KDL::Twist arg2value;
    Composition_RotationTwist(){}
    Composition_RotationTwist(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("transform",arg1,arg2)
                    {}

    virtual KDL::Twist value() {
        arg1value = argument1->value();
        arg2value = argument2->value();
        return arg1value * arg2value;
    }

    virtual KDL::Twist derivative(int i) {
    	KDL::Vector da = argument1->derivative(i);
    	KDL::Twist db = argument2->derivative(i);
        return KDL::Twist(
        			arg1value*db.vel + da*(arg1value*arg2value.vel),
        			arg1value * db.rot + da * (arg1value * arg2value.rot)
        		);
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Composition_RotationTwist( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Twist>::Ptr operator*( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Twist>::Ptr a2) {
    Expression<KDL::Twist>::Ptr expr(
        new Composition_RotationTwist( a1, a2 )
    );
    return expr;
}

//Composition Twist Double
class Multiplication_TwistDouble:
    public BinaryExpression<KDL::Twist, KDL::Twist, double>
{
public:
    typedef BinaryExpression<KDL::Twist, KDL::Twist, double> BinExpr;
    KDL::Twist arg1value;
    double arg2value;
    Multiplication_TwistDouble() {}
    Multiplication_TwistDouble(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("multiplication",arg1,arg2)
                    {}

    virtual KDL::Twist value() {
        arg1value = argument1->value();
        arg2value = argument2->value();
        return arg1value * arg2value;
    }

    virtual KDL::Twist derivative(int i) {
    	return arg1value*argument2->derivative(i) + argument1->derivative(i) * arg2value;
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Multiplication_TwistDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Twist>::Ptr operator*( Expression<KDL::Twist>::Ptr a1, Expression<double>::Ptr a2) {
    Expression<KDL::Twist>::Ptr expr(
        new Multiplication_TwistDouble( a1, a2 )
    );
    return expr;
}
inline Expression<KDL::Twist>::Ptr operator*( Expression<double>::Ptr a1, Expression<KDL::Twist>::Ptr a2) {
    Expression<KDL::Twist>::Ptr expr(
        new Multiplication_TwistDouble( a2, a1 )
    );
    return expr;
}

//RefPoint Twist Vector
class RefPoint_TwistVector:
    public BinaryExpression<KDL::Twist, KDL::Twist, KDL::Vector>
{
public:
    typedef BinaryExpression<KDL::Twist, KDL::Twist, KDL::Vector> BinExpr;
    KDL::Twist arg1value;
    KDL::Vector arg2value;
    RefPoint_TwistVector() {}
    RefPoint_TwistVector(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("ref_point",arg1,arg2)
                    {}

    virtual KDL::Twist value() {
    	arg1value = argument1->value();
    	arg2value = argument2->value();
        return arg1value.RefPoint(arg2value);
    }

    virtual KDL::Twist derivative(int i) {
    	KDL::Twist da = argument1->derivative(i);
    	return KDL::Twist(
    			da.vel + da.rot * arg2value + arg1value.rot * argument2->derivative(i),
    			da.rot
    	);
    }

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new RefPoint_TwistVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Twist>::Ptr ref_point ( Expression<KDL::Twist>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<KDL::Twist>::Ptr expr(
        new RefPoint_TwistVector( a1, a2 )
    );
    return expr;
}

/**
 * /todo Finish composition with stiffness
 *
//Composition Stiffness, Twist
class Composition_StiffnessTwist:
	public BinaryExpression<KDL::Wrench, KDL::Stiffness, KDL::Twist>
{
public:
	typedef BinaryExpression<KDL::Wrench,KDL::Stiffness,KDL::Twist> BinExpr;
private:
	KDL::Stiffness arg1value;
	KDL::Twist arg2value;
public:
	Composition_StiffnessTwist(
			const typename BinExpr::Argument1Expr::Ptr& arg1,
			const typename BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("composition",arg1,arg2)
				{}

	virtual KDL::Wrench value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Wrench derivative(int i){
		return arg1value * argument2->derivative(i) + argument1->derivative(i)*arg2value;
	}

    virtual typename BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Composition_StiffnessTwist( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Wrench>::Ptr operator* ( Expression<KDL::Stiffness>::Ptr a1, Expression<KDL::Twist>::Ptr a2 ) {
	Expression<KDL::Wrench>::Ptr expr(new Composition_StiffnessTwist( a1, a2 ));
	return expr;
}
*/

} // end of namespace KDL

#endif
