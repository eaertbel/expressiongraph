/*
* expressiontree_function1.cpp
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

void tst1() {
        using namespace KDL;
	    using namespace std;


        // Definition of the function:
        auto myfunc = boost::make_shared<FunctionDefinition>("myfunc") ;
        auto arg1 = make_FunctionParameter<double>(myfunc,"arg1");
        auto arg2 = make_FunctionParameter<double>(myfunc,"arg2");
        myfunc->setBodyExpression<double>( Constant<double>(2.0)*arg1+cached<double>(Constant<double>(3.0)*arg2+input(4))  ); 
      
        // constructing call object: 
        //Expression<double>::Ptr funceval (new FunctionEvaluation<double>(myfunc,{input(2), input(3)}));
        FunctionEvaluation<double>::Ptr funceval =  make_FunctionEvaluation<double>(myfunc);
        funceval->addTypeCheckedArgument<double>("arg1",input(2));
        funceval->addTypeCheckedArgument<double>("arg2",input(3));
        Expression<double>::Ptr result = funceval;
        //,{input(2), input(3)}));
        // Using the expression: 
        result->setInputValue(2,1.0);
        result->setInputValue(3,2.0);
        result->setInputValue(4,2.0);
        cout << "value       : " << result->value() << endl;
        cout << "gradient 1  : " << result->derivative(1) << endl;
        cout << "gradient 2  : " << result->derivative(2) << endl;
        cout << "gradient 3  : " << result->derivative(3) << endl;
        cout << "gradient 4  : " << result->derivative(4) << endl;
        // constructing another call object: 
        cout << "Evaluate 2 myfunc " << endl;
        auto funceval2 = make_FunctionEvaluation<double>(myfunc,{input(2)*input(2), input(3)*input(3)});
        
        // Using the expression: 
        funceval2->setInputValue(2,1.0);
        funceval2->setInputValue(3,2.0);
        funceval2->setInputValue(4,2.0);
        cout << funceval2->value() << endl;
        cout << "gradient 1  : " << funceval2->derivative(1) << endl;
        cout << "gradient 2  : " << funceval2->derivative(2) << endl;
        cout << "gradient 3  : " << funceval2->derivative(3) << endl;
        cout << "gradient 4  : " << funceval2->derivative(4) << endl;
    }

// ==================== Testing Rotation:
void tst2() { 
         using namespace KDL;
	    using namespace std;


        cout << "define function " << endl;
        cout << "  FunctionDefinition(...)" << endl;
        auto myfunc = make_FunctionDefinition("myfunc") ;
        cout << "  make_FunctionParameter(...) " << endl;
        auto arg1 = make_FunctionParameter<Rotation>(myfunc,"R");
        cout << "  setBodyExpression(...) " << endl;
        myfunc->setBodyExpression<Rotation>(  arg1*arg1  ); 
        
        cout << "evaluate function " << endl; 
        Expression<Rotation>::Ptr argR = rot_x(input(1));
        auto funceval = make_FunctionEvaluation<Rotation>(myfunc, {argR}); 
        funceval->setInputValue(1,0.1);
        cout << "value " << funceval->value() << endl;
        cout << "gradient 1 " << funceval->derivative(1) << endl;
}
void tst3() {
        using namespace KDL;
	    using namespace std;

        cout << "================= DEFINITION =================" << std::endl; 
        // Definition of the function:
        auto myfunc = boost::make_shared<FunctionDefinition>("myfunc") ;
        auto arg1 = make_FunctionParameter<double>(myfunc,"arg1");
        myfunc->setBodyExpression<double>( arg1*cached<double>("cache_arg1",arg1)  ); 
     
        cout << "================= EVALUATION NODE =================" << std::endl; 
        // constructing call object: 
        //Expression<double>::Ptr funceval (new FunctionEvaluation<double>(myfunc,{input(2), input(3)}));
        FunctionEvaluation<double>::Ptr funceval =  make_FunctionEvaluation<double>(myfunc);
        funceval->addTypeCheckedArgument<double>("arg1",input(2));
        Expression<double>::Ptr result = funceval;
        //,{input(2), input(3)}));
        // Using the expression: 
        cout << "================= EVALUATION =================" << std::endl; 
        result->setInputValue(2,2.0);
        cout << "value       : " << result->value() << endl;
        cout << "gradient 2  : " << result->derivative(2) << endl;
    }



int main(int argc, char* argv[]) {
    tst3();
    return 0;
}

