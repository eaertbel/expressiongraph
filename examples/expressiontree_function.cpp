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
#include <expressiongraph/expressiontree_function.hpp>









int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;


    // Definition of the function:
    auto myfunc = boost::make_shared<FunctionDefinition>() ;
    myfunc->setName("myfunc");
    myfunc->addParam("arg1",ExpressionType::expression_double);
    myfunc->addParam("arg2",ExpressionType::expression_double);
    Expression<double>::Ptr arg1 = functionParameter<double>(myfunc,"arg1");
    Expression<double>::Ptr arg2 = functionParameter<double>(myfunc,"arg2");
 
    myfunc->set_body_expression( Constant<double>(2.0)*arg1+Constant<double>(3.0)*arg2+input(4)  ); 
    
    myfunc->get_body_expression<double>()->print(std::cout);
    std::cout << std::endl;
    // Using the function: 
    auto funceval = new FunctionEvaluation<double>(myfunc);
    funceval->addTypeCheckedArgument("arg1", input(2));
    funceval->addTypeCheckedArgument("arg2", input(3));
    funceval->setInputValue(2,1.0);
    funceval->setInputValue(3,2.0);
    funceval->setInputValue(4,2.0);
    cout << "funceval->setInputValue(2,1.0);" << endl;
    cout << "funceval->setInputValue(3,2.0);" << endl;
    cout << "funceval->setInputValue(4,2.0);" << endl;
 
    cout << "Constant<double>(2.0)*arg1+Constant<double>(3.0)*arg2+input(4)" << endl;
    cout << "arg 1 : input(2) " << endl;
    cout << "arg 2 : input(3) " << endl;
    cout << funceval->value() << endl;
    cout << "gradient 1  : " << funceval->derivative(1) << endl;
    cout << "gradient 2  : " << funceval->derivative(2) << endl;
    cout << "gradient 3  : " << funceval->derivative(3) << endl;
    cout << "gradient 4  : " << funceval->derivative(4) << endl;

    auto funceval2 = new FunctionEvaluation<double>(myfunc);
    funceval2->addTypeCheckedArgument("arg1", input(2)*input(2));
    funceval2->addTypeCheckedArgument("arg2", input(3)*input(3));
    funceval2->setInputValue(2,1.0);
    funceval2->setInputValue(3,2.0);
    funceval2->setInputValue(4,2.0);
    cout << "Constant<double>(2.0)*arg1+Constant<double>(3.0)*arg2+input(4)" << endl;
    cout << "arg 1 : input(2)*input(2)" << endl;
    cout << "arg 2 : input(3)*input(3)" << endl;
    cout << funceval2->value() << endl;
    cout << "gradient 1  : " << funceval2->derivative(1) << endl;
    cout << "gradient 2  : " << funceval2->derivative(2) << endl;
    cout << "gradient 3  : " << funceval2->derivative(3) << endl;
    cout << "gradient 4  : " << funceval2->derivative(4) << endl;


    cout << funceval->value() << endl;
    cout << funceval2->value() << endl;
    cout << "gradient 2  : " << funceval->derivative(2) << endl;
    cout << "gradient 2  : " << funceval2->derivative(2) << endl;
    cout << "gradient 2  : " << funceval->derivative(2) << endl;

    //--------------------------------------------------------------
    {
        auto myfunc = boost::make_shared<FunctionDefinition>() ;
        myfunc->setName("myfunc");
        myfunc->addParam("arg1",ExpressionType::expression_vector);
        myfunc->addParam("arg2",ExpressionType::expression_double);
        Expression<Vector>::Ptr arg1 = functionParameter<Vector>(myfunc,"arg1");
        Expression<double>::Ptr arg2 = functionParameter<double>(myfunc,"arg2");
     
        myfunc->set_body_expression( norm(arg1)*arg2 ); 
        
        myfunc->get_body_expression<double>()->print(std::cout);
        std::cout << std::endl;
        // Using the function: 
        auto funceval = new FunctionEvaluation<double>(myfunc);
        funceval->addTypeCheckedArgument("arg1", KDL::vector(input(1),input(2),input(3)));
        funceval->addTypeCheckedArgument("arg2", input(4));
        funceval->setInputValue(1,3.0);
        funceval->setInputValue(2,1.0);
        funceval->setInputValue(3,2.0);
        funceval->setInputValue(4,2.0);
        cout << "funceval->setInputValue(1,3.0);" << endl;
        cout << "funceval->setInputValue(2,1.0);" << endl;
        cout << "funceval->setInputValue(3,2.0);" << endl;
        cout << "funceval->setInputValue(4,2.0);" << endl;
        cout << "myfunc->set_body_expression( norm(arg1)*arg2 ); " << endl;
        cout << "funceval->addTypeCheckedArgument(arg1, vector(input(1),input(2),input(3)));" << endl;
        cout << "funceval->addTypeCheckedArgument(arg2, input(4));" << endl;
        cout << funceval->value() << endl;
        cout << "gradient 1  : " << funceval->derivative(1) << endl;
        cout << "gradient 2  : " << funceval->derivative(2) << endl;
        cout << "gradient 3  : " << funceval->derivative(3) << endl;
        cout << "gradient 4  : " << funceval->derivative(4) << endl;
    } 
	//--------------------------------------------------------------
    {
        cout << "-----------------------------------------------------------------" << endl;
        cout << "Definition without dummy body : " << endl;
        auto factorial = boost::make_shared<FunctionDefinition>() ;
        factorial->setName("factorial");
        factorial->addParam("N",ExpressionType::expression_double);
        auto argN = functionParameter<double>(factorial,"N");
        // dummy body expression during the construction of the recursive call is necessary: 
        factorial->set_body_expression( Constant<double>(1.0) );
 
        cout << "creating evaluator node : " << endl;
        auto fact_eval = boost::make_shared<FunctionEvaluation<double>>(factorial);
        fact_eval->addTypeCheckedArgument("N", argN-Constant<double>(1.0));
        cout << "setting body : " << endl;
        factorial->set_body_expression( KDL::conditional<double>(argN, argN*fact_eval, Constant<double>(1.0)) ); 
        cout << endl; 

        cout << "parameter N: ";
        argN->print(cout);
        cout << endl;
        cout << "fact_eval: ";
        fact_eval->print(cout); 
        cout << endl;
        cout << "setInputValue()" << endl; 
        fact_eval ->setInputValue(1,3.0); 
        cout << "Body Expression ";
        factorial->get_body_expression<double>()->print(cout);
        cout << endl;
        // Using the function: 
        cout << "---------------------" << endl;
        cout << "Evaluation of the function : " << endl;
        auto fact_eval2 = boost::make_shared<FunctionEvaluation<double>>(factorial);
        fact_eval2->addTypeCheckedArgument("N",input(2));
        fact_eval2->setInputValue(2,3.0);
        fact_eval2->print(cout); cout << endl;
        cout << "Value : " << fact_eval2->value()  << endl;
        cout << "gradient 1  : " << fact_eval2->derivative(1) << endl;
        cout << "gradient 2  : " << fact_eval2->derivative(2) << endl;
    } 
	
	return 0;
}
