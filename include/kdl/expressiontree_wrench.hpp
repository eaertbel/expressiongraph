/*
 * expressiontree_wrench.hpp
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

#ifndef KDL_EXPRESSIONTREE_WRENCH_HPP
#define KDL_EXPRESSIONTREE_WRENCH_HPP

#include <kdl/expressiontree_expressions.hpp>
#include <kdl/frames.hpp>

namespace KDL {

/*
 * Wrench operations:
 *     -Wrench                		returns Wrench
 *     Wrench+Wrench           		returns Wrench
 *     Wrench-Wrench           		returns Wrench
 *     Rotation*Wrench        		returns Wrench
 *     Wrench*double          		returns Wrench
 *     double*Wrench          		returns Wrench
 *     RefPoint(Wrench,Vector) 		returns Wrench
 *     Inverse(Stiffness, Wrench)	returns Twist
 */

class Wrench_VectorVector:
    public BinaryExpression<Wrench, Vector, Vector>
{
public:
    Wrench_VectorVector(){}
    Wrench_VectorVector( const Expression<Vector>::Ptr& arg1,
                  const Expression<Vector>::Ptr& arg2  ):
                BinaryExpression<Wrench,Vector,Vector>("wrench",arg1,arg2)
                {}

    virtual Wrench value() {
        return Wrench(argument1->value(),argument2->value());
    }

    virtual Wrench derivative(int i) {
        return Wrench(argument1->derivative(i),argument2->derivative(i));
    }

    virtual Expression<Wrench>::Ptr derivativeExpression(int i);

    virtual Expression<Wrench>::Ptr clone() {
        Expression<Wrench>::Ptr expr(
            new Wrench_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<Wrench>::Ptr wrench( Expression<KDL::Vector>::Ptr a, Expression<KDL::Vector>::Ptr b) {
    Expression<KDL::Wrench>::Ptr expr(
        new Wrench_VectorVector( a,b )
    );
    return expr;
}

class Force_Wrench:
    public UnaryExpression<Vector, Wrench>
{
public:
    typedef UnaryExpression<Vector, Wrench> UnExpr;
public:
    Force_Wrench() {}
    Force_Wrench(const Expression<Wrench>::Ptr& arg):
                UnExpr("force",arg)
                {}

    virtual Vector value() {
    	return argument->value().force;
    }

    virtual Vector derivative(int i) {
    	return argument->derivative(i).force;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<Vector>::Ptr expr(
            new Force_Wrench( argument->clone())
        );
        return expr;
    }
};

inline Expression<Vector>::Ptr force( Expression<Wrench>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Force_Wrench( a )
    );
    return expr;
}

class Torque_Wrench:
    public UnaryExpression<Vector, Wrench>
{
public:
    typedef UnaryExpression<Vector, Wrench> UnExpr;
public:
    Torque_Wrench(){}
    Torque_Wrench(
                  Expression<Wrench>::Ptr arg):
                UnExpr("torque",arg)
                {}

    virtual Vector value() {
    	return argument->value().torque;
    }

    virtual Vector derivative(int i) {
    	return argument->derivative(i).torque;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<Vector>::Ptr expr(
            new Torque_Wrench( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Vector>::Ptr torque( Expression<KDL::Wrench>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Torque_Wrench( a )
    );
    return expr;
}


//Negate Wrench
class Negate_Wrench:
    public UnaryExpression<KDL::Wrench, KDL::Wrench>
{
public:
    typedef UnaryExpression<KDL::Wrench, KDL::Wrench> UnExpr;
public:
    Negate_Wrench(){}
    Negate_Wrench(
                const   UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("negation",arg)
                {}

    virtual KDL::Wrench value() {
    	return -argument->value();
    }

    virtual KDL::Wrench derivative(int i) {
    	return -argument->derivative(i);
    }
    virtual Expression<Wrench>::Ptr derivativeExpression(int i);

    virtual   UnExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Negate_Wrench( argument->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Wrench>::Ptr operator-( Expression<KDL::Wrench>::Ptr a) {
    Expression<KDL::Wrench>::Ptr expr(
        new Negate_Wrench( a )
    );
    return expr;
}

//Addition Wrench Wrench
class Addition_WrenchWrench:
    public BinaryExpression<KDL::Wrench, KDL::Wrench, KDL::Wrench>
{
public:
    typedef BinaryExpression<KDL::Wrench, KDL::Wrench, KDL::Wrench> BinExpr;
public:
    Addition_WrenchWrench(){}
    Addition_WrenchWrench(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("addition",arg1,arg2)
                    {}

    virtual KDL::Wrench value() {
        return argument1->value() + argument2->value();
    }

    virtual KDL::Wrench derivative(int i) {
        return argument1->derivative(i) + argument2->derivative(i);
    }
    virtual Expression<Wrench>::Ptr derivativeExpression(int i);


    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Addition_WrenchWrench( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Wrench>::Ptr operator+( Expression<KDL::Wrench>::Ptr a1, Expression<KDL::Wrench>::Ptr a2) {
    Expression<KDL::Wrench>::Ptr expr(
        new Addition_WrenchWrench( a1, a2 )
    );
    return expr;
}

//Subtraction Wrench Wrench
class Subtraction_WrenchWrench:
    public BinaryExpression<KDL::Wrench, KDL::Wrench, KDL::Wrench>
{
public:
    typedef BinaryExpression<KDL::Wrench, KDL::Wrench, KDL::Wrench> BinExpr;
public:
    Subtraction_WrenchWrench(){}
    Subtraction_WrenchWrench(
                    const   BinExpr::Argument1Expr::Ptr& arg1,
                    const   BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("subtraction",arg1,arg2)
                    {}
    virtual KDL::Wrench value() {
        return argument1->value() - argument2->value();
    }

    virtual KDL::Wrench derivative(int i) {
        return argument1->derivative(i) - argument2->derivative(i);
    }
    virtual Expression<Wrench>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Subtraction_WrenchWrench( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Wrench>::Ptr operator-( Expression<KDL::Wrench>::Ptr a1, Expression<KDL::Wrench>::Ptr a2) {
    Expression<KDL::Wrench>::Ptr expr(
        new Subtraction_WrenchWrench( a1, a2 )
    );
    return expr;
}

//Composition Rotation Wrench
class Composition_RotationWrench:
	public BinaryExpression<KDL::Wrench, KDL::Rotation, KDL::Wrench>
{
public:
	typedef BinaryExpression<KDL::Wrench,KDL::Rotation,KDL::Wrench> BinExpr;
	KDL::Rotation arg1value;
	KDL::Wrench arg2value;
public:
	Composition_RotationWrench() {}
	Composition_RotationWrench(
			const   BinExpr::Argument1Expr::Ptr& arg1,
			const   BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("transform",arg1,arg2)
				{}

	virtual KDL::Wrench value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Wrench derivative(int i){
		KDL::Vector da = argument1->derivative(i);
		KDL::Wrench db = argument2->derivative(i);
		return KDL::Wrench(
				arg1value * db.force + da*(arg1value*arg2value.force),
				arg1value * db.torque + da * (arg1value*arg2value.torque)
		);
	}
    virtual Expression<Wrench>::Ptr derivativeExpression(int i);


    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Composition_RotationWrench( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Wrench>::Ptr operator* ( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Wrench>::Ptr a2 ) {
	Expression<KDL::Wrench>::Ptr expr(new Composition_RotationWrench( a1, a2 ));
	return expr;
}

//Composition Wrench Double
class Multiplication_WrenchDouble:
	public BinaryExpression<KDL::Wrench, KDL::Wrench, double>
{
public:
	typedef BinaryExpression<KDL::Wrench,KDL::Wrench,double> BinExpr;
	KDL::Wrench arg1value;
	double arg2value;
public:
	Multiplication_WrenchDouble() {}
	Multiplication_WrenchDouble(
			const   BinExpr::Argument1Expr::Ptr& arg1,
			const   BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("multiplication",arg1,arg2)
				{}

	virtual KDL::Wrench value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Wrench derivative(int i){
		return arg1value*argument2->derivative(i) + argument1->derivative(i)*arg2value;
	}
    virtual Expression<Wrench>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Multiplication_WrenchDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Wrench>::Ptr operator* ( Expression<KDL::Wrench>::Ptr a1, Expression<double>::Ptr a2 ) {
	Expression<KDL::Wrench>::Ptr expr(new Multiplication_WrenchDouble( a1, a2 ));
	return expr;
}
inline Expression<KDL::Wrench>::Ptr operator* ( Expression<double>::Ptr a1, Expression<KDL::Wrench>::Ptr a2 ) {
	Expression<KDL::Wrench>::Ptr expr(new Multiplication_WrenchDouble( a2, a1 ));
	return expr;
}

//RefPoint Wrench Vector
class RefPoint_WrenchVector:
	public BinaryExpression<KDL::Wrench, KDL::Wrench, KDL::Vector>
{
public:
	typedef BinaryExpression<KDL::Wrench,KDL::Wrench,KDL::Vector> BinExpr;
	KDL::Wrench arg1value;
	KDL::Vector arg2value;
public:
	RefPoint_WrenchVector(){}
	RefPoint_WrenchVector(
			const   BinExpr::Argument1Expr::Ptr& arg1,
			const   BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("ref_point",arg1,arg2)
				{}

	virtual KDL::Wrench value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value.RefPoint(arg2value);
	}

	virtual KDL::Wrench derivative(int i){
		KDL::Wrench da = argument1->derivative(i);
		return KDL::Wrench(	da.force,
						da.torque + da.force*arg2value + arg1value.force*argument2->derivative(i));
	}
    virtual Expression<Wrench>::Ptr derivativeExpression(int i);

    virtual   BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new RefPoint_WrenchVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<KDL::Wrench>::Ptr ref_point ( Expression<KDL::Wrench>::Ptr a1, Expression<KDL::Vector>::Ptr a2 ) {
	Expression<KDL::Wrench>::Ptr expr(new RefPoint_WrenchVector( a1, a2 ));
	return expr;
}

/*  **********************************************
//Inverse Stiffness Wrench
class Inverse_StiffnessWrench:
	public BinaryExpression<KDL::Twist, KDL::Stiffness, KDL::Wrench>
{
public:
	typedef BinaryExpression<KDL::Twist,KDL::Stiffness,KDL::Wrench> BinExpr;
private:
	KDL::Stiffness arg1value;
	KDL::Wrench arg2value;
public:
	Inverse_StiffnessWrench(
			const typename BinExpr::Argument1Expr::Ptr& arg1,
			const typename BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("inverse",arg1,arg2)
				{}

	virtual KDL::Twist value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value.Inverse(arg2value);
	}

	virtual KDL::Twist derivative(int i){ //TODO: Check implementation! The old implementation asserts false in the VV case.
		return arg1value.Inverse(argument2->derivative(i));
	}
    virtual Expression<Twist>::Ptr derivativeExpression(int i) {
        assert( 0 && "derivativeExpression for InverseStiffnessWrench NOT IMPLEMENTED (YET)");
    }


    virtual typename BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Inverse_StiffnessWrench( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Twist>::Ptr inv ( Expression<KDL::Stiffness>::Ptr a1, Expression<KDL::Wrench>::Ptr a2 ) {
	Expression<KDL::Twist>::Ptr expr(new Inverse_StiffnessWrench( a1, a2 ));
	return expr;
}
**************************************************** */

} // end of namespace KDL

#endif
