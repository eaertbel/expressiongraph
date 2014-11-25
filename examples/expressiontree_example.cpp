/*
 * expressiontree_example.cpp
 *
 *  Created on: Aug 7, 2012
 *      Author: Wouter Bancken
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

/*
 * Example on how to use VariableType
 */
int main(int argc, char* argv[]) {
	using namespace KDL;

    VariableType<Rotation>::Ptr v = Variable<Rotation>(0,2);
    v->setValue( Rotation::EulerZYX(0.1, 0.2, 0.3) ); 
    v->setJacobian(0,  KDL::Vector(0,0,1) );  // be carefull: "0" refers to the 0-th column of the Jacobian, not to the 0-th variable.
    v->setJacobian(1,  KDL::Vector(1,0,0) );
    std::cout << std::endl;
    std::cout << v->ndx.size() << std::endl;
	Expression<Rotation>::Ptr result;
	result = Constant(Rotation::EulerZYX(0.01, 0.02, 0.03))*v; 

	std::cout << "expression tree:\n";
	result->debug_printtree();
	std::cout << std::endl;

	std::cout << "\n\nvalue \n";
	std::cout << result->value();
	std::cout << "\n\nderivative 0 \n";
	std::cout << result->derivative(0);
	std::cout << "\n\nderivative 1 \n";
	std::cout << result->derivative(1);
	std::cout << "\n\nderivative 2 \n";
	std::cout << result->derivative(2);
    std::cout << " \n\nnumber of derivatives ";
    std::cout << result->number_of_derivatives();
	std::cout << "\n\n";

	/* second example: */
    VariableType<Vector>::Ptr vec = Variable<Vector>(1,1);
    vec->setValue( KDL::Vector(4,5,6) );
    vec->setJacobian(0, KDL::Vector(0,1,0) );
	Expression<double>::Ptr result2 = norm(vec);

	std::cout << "expression tree:\n";
	result2->debug_printtree();
	std::cout << std::endl;

	std::cout << "\n\nvalue \n";
	std::cout << result2->value();
	std::cout << "\n\nderivative 0 \n";
	std::cout << result2->derivative(0);
	std::cout << "\n\nderivative 1 \n";
	std::cout << result2->derivative(1);
	std::cout << "\n\n";
	return 0;
}
