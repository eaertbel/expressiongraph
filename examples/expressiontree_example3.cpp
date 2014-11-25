/*
 * expressiontree_example3.cpp
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

#include <kdl/expressiontree.hpp>

/**
 * An example of a <double> expression in multiple variables.
 */
int main(int argc, char* argv[]) {
    using namespace KDL;
    using namespace std;
    Expression<double>::Ptr expr = sin(input(0)) + cos(input(1));
    std::vector<double> inp(2);
    inp[0] = 1;
    inp[1] = 2;
    expr->setInputValues(inp);
    cout << "Expression Tree in prefix notation : \n";
    expr->debug_printtree();
    cout << "\n\nValue                " << expr->value()       << endl;
    cout << "Derivative towards 0 " << expr->derivative(0) << endl;
    cout << "Derivative towards 1 " << expr->derivative(1) << endl;
    Expression<double>::Ptr exprd0 = expr->derivativeExpression(0);
    cout << "derivative towards 0 " << endl;
    exprd0->print(cout);
    Expression<double>::Ptr exprd1 = expr->derivativeExpression(1);
    cout << endl << "derivative towards 1 " << endl;
    exprd1->print(cout);
    cout << endl;
	return 0;
}
