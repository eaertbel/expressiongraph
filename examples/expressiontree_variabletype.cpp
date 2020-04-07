/*
 * expressiontree_variabletype.cpp
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

#include <expressiongraph/expressiontree.hpp>

/*
 * Tutorial on the use of VariableType
 * It shows:
 *  - How does it behave w.r.t. cloning
 *  - How you change the value of expressions.
 *  - It is important to always depend on some variable, even when the VariableType
 *    really expresses a constant value. Otherwise you cannot indicate a change of the
 *    VariableType
 */
int main(int argc, char* argv[]) {
	using namespace KDL;

    //---------------------------------------------------------------------
    //1. The use of VariableType:

    // create a VariableType dependend on two variables starting from index 0 
    // (i.e. input(0) and input(1) variables
    VariableType<Rotation>::Ptr v = Variable<Rotation>(0,2);
    v->setValue( Rotation::EulerZYX(0.0, 0.0, 0.0) ); 
    v->setJacobian(0,  KDL::Vector(0,0,1) );  // be carefull: "0" refers to the 0-th column of the Jacobian, not to the 0-th variable.
    v->setJacobian(1,  KDL::Vector(1,0,0) );
    std::cout << std::endl;
    // display the amount of dependend variables
    std::cout << v->ndx.size() << std::endl;
    // list the dependencies (can be done for any expression):
    std::cout << "Dependencies ";
    std::set<int> s;
    v->getDependencies(s);
    for (std::set<int>::iterator it = s.begin();it!=s.end();++it) {
        std::cout << *it << "   ";
    }
    std::cout << std::endl;
    // use it in an expression (where there is a cached node):
	Expression<Rotation>::Ptr result;
	result = cached<Rotation>(Constant(Rotation::EulerZYX(0.0, -M_PI/2, 0.0))*v); 
    std::cout << "The result before change the variable " << std::endl;
    std::cout << result->value() << std::endl;
    // change the value of the VariableType
    v->setValue( Rotation::EulerZYX(0.3,0.3,0.2));
    std::cout << "VariableType value : " << std::endl;
    std::cout << v->value() << std::endl;
    std::cout << "The next result is wrong :" << std::endl;
    std::cout << result->value() << std::endl;  // WRONG RESULT !  The expression is not warned about the change in VariableType!
                                                // and the cached variable is therefore not updated.
    std::cout << "The next result is correct :" << std::endl;
    result->setInputValue(0,0.0);
    std::cout << result->value() << std::endl;  // by touching one of the dependant variables, the cached variables are updated and now
                                  // the result is correct.
    std::cout << "expression tree:\n";
	result->print(std::cout);
	std::cout << std::endl;
    // In conclusion, the dependant variables of a VariableType serve two purposes:
    // a) indicated dependencies and when to update values
    // b) indicate variables for which there is/can be non-zero derivatives


    //---------------------------------------------------------------------

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

    // cloning the expression:
    Expression<Rotation>::Ptr cloned_result = result->clone();
    std::cout << "Value of result " << std::endl;
    std::cout << result->value() << std::endl;
    std::cout << "Value of cloned result "<< std::endl; 
    std::cout << cloned_result->value() << std::endl;
    std::cout << "Change the original VariableType " << std::endl;
    v->setValue( Rotation::EulerZYX(0.5, 1.0, -0.2) ); 
    result->setInputValue(0,0);
    cloned_result->setInputValue(0,0);
    std::cout << "Value of result after clone and change of variable value" << std::endl;
    std::cout << result->value() << std::endl;
    std::cout << "Value of cloned result after clone and change of variable value"<< std::endl; 
    std::cout << "(the value is now decoupled because the clone) " << std::endl;
    std::cout << cloned_result->value() << std::endl;
    std::cout << "But a clone of a VariableType can be updated with the state of the " << std::endl;
    std::cout << "original where it is cloned from" << std::endl;
    cloned_result->update_variabletype_from_original();
    std::cout << cloned_result->value() << std::endl; 
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
