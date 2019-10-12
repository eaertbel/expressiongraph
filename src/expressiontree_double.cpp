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

#include <kdl/expressiontree_double.hpp>
#include <iostream>
#include "util.hpp"

namespace KDL {


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
        assert( 0 && "Not yet implemented, need an ExpressionTree sign function"); 
        //return val * argument->derivativeExpression(i);
        return Expression<double>::Ptr();
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

boost::mt19937 rng; 

} // end of namespace KDL
