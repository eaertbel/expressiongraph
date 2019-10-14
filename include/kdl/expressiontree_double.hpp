/*
 * expressiontree_doubles.hpp
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

#ifndef KDL_EXPRESSIONTREE_DOUBLE_HPP
#define KDL_EXPRESSIONTREE_DOUBLE_HPP

#include <kdl/expressiontree_expressions.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>


namespace KDL {

/*
 * ExpressionTree classes and operation that involve "double"
 *
 * Operations:
 * 		Binary: + - * / atan2
 * 		Unary:  - sin cos tan asin acos exp log sqrt atan abs fmod
 * Corresponding Expression objects:
 *      Binary:
 *      -  Addition_DoubleDouble 
 *      -  Subtraction_DoubleDouble 
 *      -  Multiplication_DoubleDouble 
 *      -  Division_DoubleDouble
 *      -  Atan2_DoubleDouble
 *      Unary:
 *      -  Negate_Double
 *      -  Sin_Double
 *      -  Cos_Double
 *      -  Tan_Double
 *      -  Asin_Double
 *      -  Acos_Double
 *      -  Exp_Double
 *      -  Exp_Double
 *      -  Log_Double
 *      -  Atan_Double
 *      -  Abs_Double
 *      -  Fmod_Double
 */

// +
class Addition_DoubleDouble:
	public BinaryExpression<double, double, double>
{
public:
	typedef BinaryExpression<double,double,double> BinExpr;
public:
    
	Addition_DoubleDouble() {}

	Addition_DoubleDouble(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("add",arg1,arg2)
				{}

	virtual double value() {
		return argument1->value() + argument2->value();
	}

	virtual double derivative(int i){
		return (argument1->derivative(i) + argument2->derivative(i));
	}

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Addition_DoubleDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr operator+ ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
    if (isConstantZero(a1)) {
        return checkConstant<double>(a2);
    } 
    if (isConstantZero(a2)) {
        return checkConstant<double>(a1);
    } 
	Expression<double>::Ptr expr(new Addition_DoubleDouble( a1, a2 ));
    return expr;
}

// -
class Subtraction_DoubleDouble:
	public BinaryExpression<double, double, double>
{
public:
	typedef BinaryExpression<double,double,double> BinExpr;
public:
	Subtraction_DoubleDouble() {}

	Subtraction_DoubleDouble(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("sub",arg1,arg2)
				{}

	virtual double value() {
		return argument1->value() - argument2->value();
	}

	virtual double derivative(int i){
		return (argument1->derivative(i) - argument2->derivative(i));
	}

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Subtraction_DoubleDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

// *
class Multiplication_DoubleDouble:
	public BinaryExpression<double, double, double>
{
public:
	typedef BinaryExpression<double,double,double> BinExpr;

	double arg1value;
	double arg2value;

	Multiplication_DoubleDouble(){}
	Multiplication_DoubleDouble(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("mul",arg1,arg2)
				{}

	virtual double value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual double derivative(int i){
		return (arg1value * argument2->derivative(i) + argument1->derivative(i) * arg2value);
	}

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Multiplication_DoubleDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr operator* ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
    if (isConstantZero(a1)) {
        return Constant<double>(0);
    } 
    if (isConstantOne(a1)) {
        return a2;
    }
    if (isConstantZero(a2)) {
        return Constant<double>(0); 
    } 
    if (isConstantOne(a2)) {
        return a1;
    }
	Expression<double>::Ptr expr(new Multiplication_DoubleDouble( a1, a2 ));
	return expr;
}

// /
class Division_DoubleDouble:
	public BinaryExpression<double, double, double>
{
public:
	typedef BinaryExpression<double,double,double> BinExpr;

	double arg1value;
	double arg2value;
	
    Division_DoubleDouble() {}

    Division_DoubleDouble(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("div",arg1,arg2)
				{}

	virtual double value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value / arg2value;
	}

	virtual double derivative(int i){
		return ((arg2value*argument1->derivative(i) - arg1value * argument2->derivative(i))/arg2value/arg2value);
	}

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Division_DoubleDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr operator/ ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
    if (isConstantZero(a1)) {
        return Constant<double>(0);
    } 
	Expression<double>::Ptr expr(new Division_DoubleDouble( a1, a2 ));
	return expr;
}

// atan2
class Atan2_DoubleDouble:
	public BinaryExpression<double, double, double>
{
public:
	typedef BinaryExpression<double,double,double> BinExpr;
	double arg1value;
	double arg2value;

public:
	Atan2_DoubleDouble() {}
	Atan2_DoubleDouble(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("atan2",arg1,arg2)
				{}

	virtual double value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return atan2(arg1value,arg2value);
	}

	virtual double derivative(int i){
		double numerator = -arg1value*argument2->derivative(i)+arg2value*argument1->derivative(i);
		double denominator = arg1value * arg1value + arg2value * arg2value;
		return numerator/denominator;
	}

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Atan2_DoubleDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr atan2 ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
	Expression<double>::Ptr expr(new Atan2_DoubleDouble( a1, a2 ));
	return expr;
}

// -
class Negate_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
public:
    Negate_Double() {}

    Negate_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("negate",arg)
                {}

    virtual double value() {
        return -argument->value();
    }

    virtual double derivative(int i) {
        return -(argument->derivative(i));
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Negate_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr operator-( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Negate_Double(a)
    );
    return expr;
}

inline Expression<double>::Ptr operator- ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
    if (isConstantZero(a1)) {
        return checkConstant<double>(-a2);
    } 
    if (isConstantZero(a2)) {
        return checkConstant<double>(a1);
    } 
	Expression<double>::Ptr expr(new Subtraction_DoubleDouble( a1, a2 ));
	return expr;
}

// sin
class Sin_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Sin_Double() {}

    Sin_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("sin",arg)
                {}

    virtual double value() {
        val = argument->value();
        return sin(val);
    }

    virtual double derivative(int i) {
        return cos(val) * argument->derivative(i);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Sin_Double( argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr sin( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Sin_Double( a )
    );
    return expr;
}

// cos
class Cos_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Cos_Double() {}
    Cos_Double(
                const UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("cos",arg)
                {}

    virtual double value() {
        val = argument->value();
        return cos(val);
    }

    virtual double derivative(int i) {
         return -sin(val) * argument->derivative(i);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Cos_Double( argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr cos( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Cos_Double( a )
    );
    return expr;
}

//tan
class Tan_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;

    double val;

public:
    Tan_Double() {}
    Tan_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("tan",arg)
                {}

    virtual double value() {
        val = argument->value();
        return tan(val);
    }

    virtual double derivative(int i) {
    	double c = cos(val);
        return argument->derivative(i)/c/c;
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);


    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Tan_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr tan (Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Tan_Double(a)
    );
    return expr;
}

// asin
class Asin_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Asin_Double() {}
    Asin_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("asin",arg)
                {}

    virtual double value() {
        val = argument->value();
        return asin(val);
    }

    virtual double derivative(int i) {
        return (argument->derivative(i))/sqrt(1- val * val);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Asin_Double( argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr asin( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Asin_Double( a )
    );
    return expr;
}

// acos
class Acos_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;

public:
    Acos_Double() {}
    Acos_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("acos",arg)
                {}

    virtual double value() {
        val = argument->value();
        return acos(val);
    }

    virtual double derivative(int i) {
    	return -(argument->derivative(i))/sqrt(1- val*val);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Acos_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr acos( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Acos_Double(a)
    );
    return expr;
}

// exp
class Exp_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Exp_Double() {}
    Exp_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("exp",arg)
                {}

    virtual double value() {
        val = exp(argument->value());
        return val;
    }

    virtual double derivative(int i) {
    	return val * argument->derivative(i);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Exp_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr exp( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Exp_Double(a)
    );
    return expr;
}

// log
class Log_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Log_Double() {}
    Log_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("log",arg)
                {}

    virtual double value() {
        val = argument->value();
        return log(val);
    }

    virtual double derivative(int i) {
    	return argument->derivative(i)/val;
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Log_Double( argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr log( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Log_Double( a )
    );
    return expr;
}

//sqrt
class Sqrt_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Sqrt_Double() {}
    Sqrt_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("sqrt",arg)
                {}
    virtual double value() {
        val = sqrt(argument->value());
        return val;
    }

    virtual double derivative(int i) {
        return 0.5/val*argument->derivative(i);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Sqrt_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr sqrt( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Sqrt_Double( a )
    );
    return expr;
}

// atan
class Atan_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Atan_Double() {}
    Atan_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("atan",arg)
                {}
    virtual double value() {
        val = argument->value();
        return atan(val);
    }

    virtual double derivative(int i) {
        return argument->derivative(i)/(1+val*val);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Atan_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr atan( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Atan_Double( a )
    );
    return expr;
}

// abs
class Abs_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Abs_Double() {}
    Abs_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("abs",arg)
                {}
    virtual double value() {
        val = argument->value();
        return fabs(val);
    }

    virtual double derivative(int i) {
        return (KDL::sign(val)*argument->derivative(i));
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Abs_Double( argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr abs( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Abs_Double( a )
    );
    return expr;
}

// fmod
class Fmod_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double denominator;
public:
    Fmod_Double() {}
    Fmod_Double(const UnExpr::ArgumentExpr::Ptr& arg,
                double _denominator):
                UnExpr("fmod",arg),
                denominator(_denominator)
                {}
    virtual double value() {
        return fmod(argument->value(), denominator);
    }

    virtual double derivative(int i) {
        // Note: This is a simplified version of the derivative which
        // ignores that fmod(x,b)=infinity for x=n*b with n being an integer.
        return argument->derivative(i);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Fmod_Double(argument->clone(), denominator)
        );
        return expr;
    }
};

inline Expression<double>::Ptr fmod( Expression<double>::Ptr a, double b) {
    Expression<double>::Ptr expr(
        new Fmod_Double( a, b )
    );
    return expr;
}

//sqrt
class Sqr_Double:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
public:
    Sqr_Double() {}
    Sqr_Double(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("sqr",arg)
                {}
    virtual double value() {
        val = argument->value();
        return val*val;
    }

    virtual double derivative(int i) {
        return 2*val*argument->derivative(i);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Sqr_Double(argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr sqr(Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Sqr_Double( a )
    );
    return expr;
}


template <typename R>
class Conditional_double:
    public TernaryExpression<R,double,R,R>
{
public:
    bool condition;

    Conditional_double() {}

    Conditional_double(  typename Expression<double>::Ptr a1, 
                         typename Expression<R>::Ptr a2,
                         typename Expression<R>::Ptr a3):
        TernaryExpression<R,double,R,R>("conditional",a1,a2,a3) {
        }

    virtual R value() {
        condition = (this->argument1->value() >= 0);
        if (condition) {
            return this->argument2->value(); 
        } else {
            return this->argument3->value(); 
        }
    } 

    virtual typename AutoDiffTrait<R>::DerivType derivative(int i) {
        if (condition) {
            return this->argument2->derivative(i);
        } else {
            return this->argument3->derivative(i);
        }
    } 

    virtual typename Expression<typename AutoDiffTrait<R>::DerivType >::Ptr derivativeExpression(int i) {
        typename Expression<typename AutoDiffTrait<R>::DerivType>::Ptr expr( 
            new Conditional_double<typename AutoDiffTrait<R>::DerivType>( 
                this->argument1, 
                this->argument2->derivativeExpression(i), 
                this->argument3->derivativeExpression(i) )
        );
        return expr;
    }

    virtual typename Expression<R>::Ptr clone() {
        typename Expression<R>::Ptr expr(
            new Conditional_double(this->argument1->clone(), this->argument2->clone(), this->argument3->clone())
        );
        return expr;
    } 
};


/**
 * returns a2 when a1 >= 0 otherwise returns a3.
 *
 * The user is responsible for the proper use of the function w.r.t. the derivatives.
 * - When the value is discontinuous when argument==0, then the derivative is not defined, only the left and right limit
 * of the derivative to the value.
 * - For a continuous derivative, the derivative of the 2nd argument and the 3th arguments should be the same when the
 * first argument == 0.
 */
template <typename R>
inline typename Expression<R>::Ptr conditional( typename Expression<double>::Ptr a1, 
                                       typename Expression<R>::Ptr a2,
                                       typename Expression<R>::Ptr a3) {
    if (!a1 || !a2 || !a3) {
        throw std::out_of_range("conditional: null pointer is given as one of the arguments");
    }
    std::set<int> vset;
    a1->getDependencies(vset);
    if (vset.empty()) {
        double value = a1->value();
        if (value >= 0) {
            return a2;
        } else {
            return a3;
        }
    } else {
        typename Expression<R>::Ptr expr(
            new Conditional_double<R>( a1, a2, a3 )
        );
        return expr;
    }
}





template <typename R>
class NearZero_double:
    public TernaryExpression<R,double,R,R>
{
public:
    bool condition;
    double tolerance;

    NearZero_double() {}

    NearZero_double(  typename Expression<double>::Ptr a1, 
                         typename Expression<R>::Ptr a2,
                         typename Expression<R>::Ptr a3,
                         double _tolerance):
        TernaryExpression<R,double,R,R>("near_zero",a1,a2,a3),
        tolerance(_tolerance) {
        }

    virtual R value() {
        double val = this->argument1->value();
        condition = (-tolerance<=val) && (val<=tolerance);
        if (condition) {
            return this->argument2->value(); 
        } else {
            return this->argument3->value(); 
        }
    } 

    virtual typename AutoDiffTrait<R>::DerivType derivative(int i) {
        if (condition) {
            return this->argument2->derivative(i);
        } else {
            return this->argument3->derivative(i);
        }
    } 

    virtual typename Expression<typename AutoDiffTrait<R>::DerivType>::Ptr derivativeExpression(int i) {
        typename Expression<typename AutoDiffTrait<R>::DerivType>::Ptr expr( 
            new NearZero_double<typename AutoDiffTrait<R>::DerivType>( 
                this->argument1, 
                this->argument2->derivativeExpression(i), 
                this->argument3->derivativeExpression(i),
                tolerance )
        );
        return expr;
    }

    virtual typename Expression<R>::Ptr clone() {
        typename Expression<R>::Ptr expr(
            new NearZero_double(this->argument1->clone(), this->argument2->clone(), this->argument3->clone(), tolerance)
        );
        return expr;
    } 
};


/**
 * checks whether an expression tree is near zero.
 * @param a1 expression tree in double 
 * @param tolerance 
 * @param a2 is returned when the value of a1 is near zero:  -tolerance <= a1->value() <= tolerance
 * @param a3 is returned when the value of a1 is not near zero.
 *
 * The user is responsible for the proper use of the function w.r.t. the derivatives.
 * - When the value is discontinuous when the condition is true then the derivative is not defined, only the left and right limit
 * of the derivative to the value.
 * - For a continuous derivative, the derivative of the 2nd argument and the 3th arguments should be the same when the
 * the condition is true. 
 */
template <typename R>
inline typename Expression<R>::Ptr near_zero( typename Expression<double>::Ptr a1, 
                                              double tolerance,
                                       typename Expression<R>::Ptr a2,
                                       typename Expression<R>::Ptr a3) {
    typename Expression<R>::Ptr expr(
        new NearZero_double<R>( a1, a2, a3,tolerance )
    );
    return expr;
}



/**
 * return the largest of the two expressions
 */
inline Expression<double>::Ptr maximum( Expression<double>::Ptr a1, Expression<double>::Ptr a2) {
    return conditional<double>(a2-a1, a2, a1);
}

/**
 * return the smallest of the two expressions
 */
inline Expression<double>::Ptr minimum( Expression<double>::Ptr a1, Expression<double>::Ptr a2) {
    return conditional<double>(a1-a2, a2, a1);
}

/**
 * returns the expression a, saturated between lower and upper values.
 */
inline Expression<double>::Ptr saturate( Expression<double>::Ptr a, double lower, double upper) {
    return minimum( Constant<double>(upper), maximum(Constant<double>(lower), a) );
}


template <class R>
class BlockWave:
    public UnaryExpression<R,double>
{
public:
    double period;
    R level1;
    R level2;


    BlockWave() {}

    BlockWave(  typename Expression<double>::Ptr a1, double _period, const R& _level1, const R& _level2):
        UnaryExpression<R,double>("BlockWave",a1),period(_period),level1(_level1),level2(_level2) {
        }

    virtual R value() {
        double t = this->argument->value();
        double s = t/period;
        s -= floor(s);
        if (s<=0.5) {
            return level1; 
        } else {
            return level2; 
        }
    } 

    virtual typename AutoDiffTrait<R>::DerivType derivative(int i) {
        return AutoDiffTrait<R>::zeroDerivative();
    } 

    virtual typename Expression<typename AutoDiffTrait<R>::DerivType>::Ptr derivativeExpression(int i) {
        return Constant(AutoDiffTrait<R>::zeroDerivative());
    }

    virtual typename Expression<R>::Ptr clone() {
        typename Expression<R>::Ptr expr(
            new BlockWave<R>(this->argument->clone(), period, level1, level2 )
        );
        return expr;
    } 
};

class BlockWave_double:
    public UnaryExpression<double,double>
{
public:
    double period;
    double level1;
    double level2;


    BlockWave_double() {}

    BlockWave_double(   Expression<double>::Ptr a1, double _period, double _level1, double _level2):
        UnaryExpression<double,double>("BlockWave",a1),period(_period),level1(_level1),level2(_level2) {
        }

    virtual double value() {
        double t = this->argument->value();
        double s = t/period;
        s -= floor(s);
        if (s<=0.5) {
            return level1; 
        } else {
            return level2; 
        }
    } 

    virtual double derivative(int i) {
        return 0.0;
    } 

    virtual  Expression<double>::Ptr derivativeExpression(int i) {
        return Constant(0.0);
    }

    virtual  Expression<double>::Ptr clone() {
         Expression<double>::Ptr expr(
            new BlockWave_double(this->argument->clone(), period, level1, level2 )
        );
        return expr;
    } 
};

/**
 * An expression tree node for a blockwave signal.
 *
 * @code
 *
 *      -----   ----- 
 *  ----|   |---|   |---
 * 
 * @endcode
 *   
 * It has the following parameters:
 * - period : time between two rising flanks.
 * - level1 : value of the level of the signal at time 0
 * - level2 : value of the  level of the signal staring at time 0 + period/2
 *
 */
template <class R>
inline typename Expression<R>::Ptr blockwave( typename Expression<double>::Ptr a1,double period, const R& level1, const R& level2) {
    typename Expression<R>::Ptr expr(
        new BlockWave<R>( a1, period, level1, level2 )
    );
    return expr;
}


extern boost::mt19937 rng; 

class NormalDistributedNoise_double:
    public UnaryExpression<double,double>
{

    boost::normal_distribution<> nd;

    boost::variate_generator<boost::mt19937&, 
                           boost::normal_distribution<> > noise;

public:
    double stddev;

    NormalDistributedNoise_double():stddev(1.0),nd(0.0,1.0), noise(rng, nd) {}

    NormalDistributedNoise_double( Expression<double>::Ptr a1, double _stddev):
        UnaryExpression<double,double>("normal_distributed_noise",a1), 
        stddev(_stddev), 
        nd(0.0,_stddev),
        noise(rng,nd) {}

    virtual double value() {
        double time = this->argument->value();
        if (time==0) {
            return 0.0;
        } else {
            return noise();
        }
    } 

    virtual double derivative(int i) {
        return 0.0;
    }
    virtual  Expression<double>::Ptr derivativeExpression(int i) {
        return Constant(0.0);
    }

    virtual  Expression<double>::Ptr clone() {
         Expression<double>::Ptr expr(
            new NormalDistributedNoise_double(this->argument->clone(), stddev )
        );
        return expr;
    } 
};
/**
 * An expression tree node for normal-distributed random noise
 * (with zero derivative)
 *   
 * @param time  It should have time as an argument, just to ensure updates to the computed value at
 * each time interval.
 *
 * @param stddev The standard deviation of the generated random noise 
 *
 * if time == 0, then the noise is also zero ( such that the convergence criteria for the initialization procedure
 * still works).
 */
inline typename Expression<double>::Ptr normaldistributednoise( typename Expression<double>::Ptr a1,double stddev) {
    typename Expression<double>::Ptr expr(
        new NormalDistributedNoise_double( a1, stddev )
    );
    return expr;
}




} // end of namespace KDL
#endif
