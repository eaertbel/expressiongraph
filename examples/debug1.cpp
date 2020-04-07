/**
 * \file expressiontree_example6.cpp
 * \brief example of using ExpressionTree in combination with KDL::Chain.
 *
 * \Author: 2014, Erwin Aertbelien 
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
#include <fstream>
#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>

/**
 * numerically comparing to ways to compute a Jacobian of a Chain.
 *  (combining expressiontree_example5.cpp and expressiontree_example6.cpp)
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;

    Expression<double>::Ptr a = input(1);
    Expression<double>::Ptr b = input(2);

    Expression<double>::Ptr r = acos( cos( a ) )-a;

    cout << "\nExpression 2 in prefix notation : \n";
    cout << "\n\n";
    
    a->setInputValue(1,0.1);
    b->setInputValue(1,0.2);

    cout << "value r " << r->value() << endl; 
    cout << "derivative towards 1  " << r->derivative(1) << endl; 
    cout << "derivative towards 2  " << r->derivative(2) << endl; 
	return 0;
}
