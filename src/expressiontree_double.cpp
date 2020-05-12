/*
 * expressiontree_doubles.cpp
 *
 *  Created on: Sept., 2012
 *      Author: Erwin Aertbelien - Wouter Bancken
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
#include <iostream>
#include "util.hpp"

namespace KDL {

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

Expression<double>::Ptr operator+ ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
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

Expression<double>::Ptr operator* ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
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

Expression<double>::Ptr operator/ ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
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

Expression<double>::Ptr atan2 ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
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

Expression<double>::Ptr operator-( Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Negate_Double(a)
    );
    return expr;
}

Expression<double>::Ptr operator- ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 ) {
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

Expression<double>::Ptr sin( Expression<double>::Ptr a) {
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

Expression<double>::Ptr cos( Expression<double>::Ptr a) {
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

Expression<double>::Ptr tan (Expression<double>::Ptr a) {
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

Expression<double>::Ptr asin( Expression<double>::Ptr a) {
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

Expression<double>::Ptr acos( Expression<double>::Ptr a) {
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

Expression<double>::Ptr exp( Expression<double>::Ptr a) {
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

Expression<double>::Ptr log( Expression<double>::Ptr a) {
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

Expression<double>::Ptr sqrt( Expression<double>::Ptr a) {
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

Expression<double>::Ptr atan( Expression<double>::Ptr a) {
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

Expression<double>::Ptr abs( Expression<double>::Ptr a) {
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
        return ::fmod(argument->value(), denominator);
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

Expression<double>::Ptr fmod( Expression<double>::Ptr a, double b) {
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

Expression<double>::Ptr sqr(Expression<double>::Ptr a) {
    Expression<double>::Ptr expr(
        new Sqr_Double( a )
    );
    return expr;
}


boost::mt19937 rng; 

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
typename Expression<double>::Ptr normaldistributednoise( typename Expression<double>::Ptr a1,double stddev) {
    typename Expression<double>::Ptr expr(
        new NormalDistributedNoise_double( a1, stddev )
    );
    return expr;
}



    Expression<double>::Ptr Addition_DoubleDouble::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument1,argument2);
        if (nr==1) {
                return Constant<double>(0);
        } else if (nr==2) {
                return argument2->derivativeExpression(i);
        } else if (nr==3) {
                return argument1->derivativeExpression(i);
        } else  {
                return argument1->derivativeExpression(i) + argument2->derivativeExpression(i);
        }
    }

    Expression<double>::Ptr Subtraction_DoubleDouble::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument1,argument2);
        switch (nr) {
            case 1:
                return Constant<double>(0);
            case 2:
                return -argument2->derivativeExpression(i);
            case 3:
                return argument1->derivativeExpression(i);
            default:
                return argument1->derivativeExpression(i) - argument2->derivativeExpression(i);
        }

    }


    Expression<double>::Ptr Multiplication_DoubleDouble::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument1,argument2);
        switch (nr) {
            case 1:
                return Constant<double>(0.0);
            case 2:
		        return argument1 * argument2->derivativeExpression(i);
            case 3:
		        return argument1->derivativeExpression(i) * argument2;
            default:
		        return (argument1 * argument2->derivativeExpression(i) + argument1->derivativeExpression(i) * argument2);
        }
    }

    Expression<double>::Ptr Division_DoubleDouble::derivativeExpression(int i) {
        Expression<double>::Ptr arg2 = cached<double>(argument2);

        int nr = getDep<double>(i,argument1,argument2);
        switch (nr) {
            case 1:
                return Constant<double>(0.0);
            case 2:
		        return (- argument2->derivativeExpression(i) * argument1) / arg2 /arg2;
            case 3:
		        return (arg2 * argument1->derivativeExpression(i) ) / arg2 /arg2;
            default:
		        return (arg2 * argument1->derivativeExpression(i) - argument2->derivativeExpression(i) * argument1) / arg2 /arg2;
        }

    }
    Expression<double>::Ptr Atan2_DoubleDouble::derivativeExpression(int i) {
        Expression<double>::Ptr arg1 = cached<double>( argument1 );
        Expression<double>::Ptr arg2 = cached<double>( argument2 );
        int nr = getDep<double>(i,argument1,argument2);
        if (nr==1) {
            return Constant<double>(0.0);
        } if (nr==2) {
            return (-arg1*argument2->derivativeExpression(i)  ) /
               (arg1*arg1 + arg2*arg2 );
        } if (nr==3) {
            return ( arg2*argument1->derivativeExpression(i) ) /
               (arg1*arg1 + arg2*arg2 );
        } else {
            return (-arg1*argument2->derivativeExpression(i) + arg2*argument1->derivativeExpression(i) ) /
               (arg1*arg1 + arg2*arg2 );
        }
    }


    Expression<double>::Ptr  Negate_Double::derivativeExpression(int i) {
        return -argument->derivativeExpression(i);
    }

    Expression<double>::Ptr Sin_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return cos(argument)*argument->derivativeExpression(i);
        }
    }

    Expression<double>::Ptr Cos_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return -sin(argument)*argument->derivativeExpression(i);
        }


    }

    Expression<double>::Ptr Tan_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<double>::Ptr c = cached<double>( cos(argument));
            return argument->derivativeExpression(i)/(c*c);
        }


    }

    Expression<double>::Ptr Asin_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<double>::Ptr val = cached<double>( argument);
            return argument->derivativeExpression(i) / sqrt( Constant(1.0) - val*val );
        }

    }

    Expression<double>::Ptr Acos_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<double>::Ptr val = cached<double>( argument);
            return -argument->derivativeExpression(i) / sqrt( Constant(1.0) - val*val );
        }
    }

    Expression<double>::Ptr Exp_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<double>::Ptr val = cached<double>( exp(argument) );
            return val * argument->derivativeExpression(i);
        }
    }

    Expression<double>::Ptr  Log_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return argument->derivativeExpression(i) / argument;
        }


    }

    Expression<double>::Ptr Sqrt_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return Constant(0.5)/sqrt(argument)*argument->derivativeExpression(i);
        }
    }

    Expression<double>::Ptr Sqr_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return Constant(2.0)*argument*argument->derivativeExpression(i);
        }
    }

    Expression<double>::Ptr Atan_Double::derivativeExpression(int i) {
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<double>::Ptr val = cached<double>( argument );
            return argument->derivativeExpression(i) / (Constant(1.0)+val*val);
        }
    }

    Expression<double>::Ptr Abs_Double::derivativeExpression(int i) {
        return conditional<double>(argument,
                 argument->derivativeExpression(i), 
                 -argument->derivativeExpression(i) 
               );
     } 

    Expression<double>::Ptr Fmod_Double::derivativeExpression(int i) {
        // Note: This is a simplified version of the derivative which
        // ignores that fmod(x,b)=infinity for x=n*b with n being an integer.
        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return argument->derivativeExpression(i);
        }
    }

template class Conditional_double<double>;
template class Conditional_double<Vector>;
template class Conditional_double<Rotation>;
template class Conditional_double<Frame>;
template class Conditional_double<Twist>;
template class Conditional_double<Wrench>;
template class Conditional_double<Quaternion>;

template class NearZero_double<double>;
template class NearZero_double<Vector>;
template class NearZero_double<Frame>;
template class NearZero_double<Rotation>;
template class NearZero_double<Twist>;
template class NearZero_double<Wrench>;
template class NearZero_double<Quaternion>;


} // end of namespace KDL
