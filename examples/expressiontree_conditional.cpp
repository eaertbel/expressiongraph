/*
 * expressiontree_conditional.cpp
 *
 *  Created on: Dec., 2012
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

#include <kdl/expressiontree.hpp>
#include <fstream>

/*
 * Example that demonstrates the use of the "conditional" function.
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;

    Expression<double>::Ptr f =
        KDL::conditional<double>(
            input(0)-Constant(5.0),
            (input(0)-Constant(5.0))*(input(0)-Constant(5.0)),
            Constant(0.0)
        );
    Expression<double>::Ptr fd = f->derivativeExpression(0);
    for (double t=0.0;t<=10;t+=0.01) {
       f->setInputValue(t);
       fd->setInputValue(t);
       cout << t << "\t" 
            << f->value() << "\t"
            << f->derivative(0) << "\t"
            << fd->value() << "\t"
            << fd->derivative(0) << "\n";
    }
    return 0;
}
