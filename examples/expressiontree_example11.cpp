/*
 * expressiontree_example11.cpp
 *
 *  Created on: Sept., 2012
 *      Author: Erwin Aertbelien 
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

#include <expressiongraph/expressiontree.hpp>

/*
 * Example that demonstrates the evaluation of a simple expression in doubles for
 * both value() and derivative()
 * and the use of derivativeExpression() to compute a derivative to an arbitrary order.
 */
int main(int argc, char* argv[]) {
	using namespace KDL;

    // note:
    //  input(0) represent variable 0, it gives back the input values given by
    //           the setInputValue.. type methods.  It always has type
    //           Expression<double>::Ptr.
    //  Constant(xx) represent a constant value (i.e. with derivatives zero).
    //           It has the type of xx.  If you specify Constant(1), you'll
    //           have the type Expression<int>::Ptr, which is probably not what
    //           you want. (There are no implicit type conversions)
    Expression<double>::Ptr result   = Constant(2.0)*input(0)*sin(Constant(2*M_PI)*input(0));
    Expression<double>::Ptr resultd  = result->derivativeExpression(0);
    Expression<double>::Ptr resultdd = resultd->derivativeExpression(0);
    
    // The columns printer are:
    //
    // time /  derivative using result / derivative using resultd /  2nd derivative using resultd / 3th derivative using resultdd
    //
    for (double t=0.0;t<10.0;t+=0.01) {
        // be sure to call the methods in this order:
        // it is important to first call setInputValue(), then value() and then derivative()
        // (since value() can cache computations to be used by derivative(), for efficiency reasons )
        result->setInputValue(t);
        resultd->setInputValue(t);
        double val = result->value();
        double deriv = result->derivative(0); // normally you would not use this but resultd->value()
        double vald  = resultd->value();
        double deriv2 = resultd->derivative(0);
        double deriv3 = resultdd->derivative(0);
        std::cout << t << "\t" << val << "\t" << deriv << "\t" << vald << "\t" << deriv2 << "\t" << deriv3 <<std::endl;
    }
    // you can also use the display() function :
    // you will have to specify the template parameter explicitly, because
    // C++ cannot deduce the type with display() :
    //      display<double>(std::cout, result);
    //      std::cout << std::endl;
	return 0;
}
