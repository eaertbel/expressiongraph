/*
* expressiontree_trait.cpp
*
*  Created on: Sept., 2020
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

int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;
    Expression<double>::Ptr e1 = cached<double>( sin(input(1))*input(2));
    cout << "e1->getResultTYpe() " << (int)e1->getResultType() << endl;

    Expression<Vector>::Ptr e2 = KDL::vector(input(1), input(2), input(3));
    cout << "e2->getResultTYpe() " << (int)e2->getResultType() << endl;

    ExpressionBase::Ptr base_e1 = e1;
    ExpressionBase::Ptr base_e2 = e2;
    cout << "base_e1->getResultType() " << (int) base_e1->getResultType() << endl;
    cout << "base_e2->getResultType() " << (int) base_e2->getResultType() << endl;
	return 0;
}
