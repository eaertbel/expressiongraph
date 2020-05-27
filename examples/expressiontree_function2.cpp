
/*
* expressiontree_function2.cpp
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
#include <expressiongraph/expressiontree_function.hpp>
#include <expressiongraph/expressiontree_integral.hpp>


int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;

    // Definition of the function:
    auto myfunc = boost::make_shared<FunctionDefinition>("myfunc") ;
    auto arg1 = make_FunctionParameter<double>(myfunc,"arg1");
    myfunc->setBodyExpression<double>( Constant<double>(2)*arg1 + Constant<double>(1.0) ); 
    
    //double a = 0.04;
    //auto ea = Constant<double>(a);
    //myfunc->setBodyExpression<double>( Constant<double>(2)*sin(ea*arg1)/ea ); 
    //myfunc->setBodyExpression<double>( Constant<double>(2)); 

    FunctionEvaluation<double>::Ptr funceval =  make_FunctionEvaluation<double>(myfunc);
    funceval->addTypeCheckedArgument<double>("arg1",input(1));
    funceval->setInputValue(1,1.0);
    cout << "value    : " << funceval->value() << endl;

    auto I = make_Integral( myfunc, Constant<double>(0),Constant<double>(10), 1E-5, 5, 10);
    cout << "Integral node " << I->value() << endl;
	return 0;
}
