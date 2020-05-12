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

#include "expressiontree_expressions.hpp"
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
Expression<double>::Ptr operator+ ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 );

Expression<double>::Ptr operator* ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 );

Expression<double>::Ptr operator/ ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 );

Expression<double>::Ptr atan2 ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 );

Expression<double>::Ptr operator-( Expression<double>::Ptr a);

Expression<double>::Ptr operator- ( Expression<double>::Ptr a1, Expression<double>::Ptr a2 );

Expression<double>::Ptr sin( Expression<double>::Ptr a);

Expression<double>::Ptr cos( Expression<double>::Ptr a);

Expression<double>::Ptr tan (Expression<double>::Ptr a);

Expression<double>::Ptr asin( Expression<double>::Ptr a);

Expression<double>::Ptr acos( Expression<double>::Ptr a);

Expression<double>::Ptr exp( Expression<double>::Ptr a);

Expression<double>::Ptr log( Expression<double>::Ptr a);

Expression<double>::Ptr sqrt( Expression<double>::Ptr a);

Expression<double>::Ptr atan( Expression<double>::Ptr a);

Expression<double>::Ptr abs( Expression<double>::Ptr a);

Expression<double>::Ptr fmod( Expression<double>::Ptr a, double b);

Expression<double>::Ptr sqr(Expression<double>::Ptr a);

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
        throw NullPointerException();
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
typename Expression<double>::Ptr normaldistributednoise( typename Expression<double>::Ptr a1,double stddev);

extern template class Conditional_double<double>;
extern template class Conditional_double<Vector>;
extern template class Conditional_double<Rotation>;
extern template class Conditional_double<Frame>;
extern template class Conditional_double<Twist>;
extern template class Conditional_double<Wrench>;
extern template class Conditional_double<Quaternion>;

extern template class NearZero_double<double>;
extern template class NearZero_double<Vector>;
extern template class NearZero_double<Frame>;
extern template class NearZero_double<Rotation>;
extern template class NearZero_double<Twist>;
extern template class NearZero_double<Wrench>;
extern template class NearZero_double<Quaternion>;



} // end of namespace KDL
#endif
