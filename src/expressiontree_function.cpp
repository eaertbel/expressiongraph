/*
* expressiongraph library
* 
* Copyright 2020 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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

#include <expressiongraph/expressiontree_function.hpp>
namespace KDL {


    template <typename T>
    void add_if_type( Arguments& args, ExpressionType vartype) {
        if (AutoDiffTrait<T>::expr_type==vartype) {
            args.push_back( Constant<T>(AutoDiffTrait<T>::zeroValue())  );
        }
    }


    /*-----------------------------------------------------------------*/
    /*                 FunctionDefinition                              */
    /*-----------------------------------------------------------------*/

    FunctionDefinition::FunctionDefinition() {
        // push the default arguments, such that the parameters
        // always point to something that doesn't cause a segmentation 
        // error.
        argstack.push_back( &default_arguments );
    }


    void FunctionDefinition::setName(const std::string& name) {
        this->name = name;
    }

    std::string FunctionDefinition::getName() const {
        return name;
    }

    void FunctionDefinition::addParam(const std::string& name, ExpressionType vartype) {
        argnames.insert( argnames.end(), name);
        argtypes.insert( argtypes.end(), vartype);
        name2ndx[name] = argnames.size()-1;
        add_if_type<double>(default_arguments, vartype); 
        add_if_type<Vector>(default_arguments, vartype); 
        add_if_type<Rotation>(default_arguments, vartype); 
        add_if_type<Frame>(default_arguments, vartype); 
        add_if_type<Wrench>(default_arguments, vartype); 
        add_if_type<Twist>(default_arguments, vartype); 
        add_if_type<Quaternion>(default_arguments, vartype); 
    }

    int FunctionDefinition::getParamNdx(const std::string& name) const {
        auto p = name2ndx.find(name);
        if (p==name2ndx.end()) {
            return -1;
        } else {
            return p->second;
        }
    }

    ExpressionType FunctionDefinition::getParamType(int i) const {
        return argtypes[i];
    }

    const std::string& FunctionDefinition::getParamName(int i) const {
        return argnames[i];
    }

    int FunctionDefinition::getNrOfParam() const {
        return argtypes.size();
    }

    void FunctionDefinition::set_body_expression( ExpressionBase::Ptr _body_expr) {
        body_expr = _body_expr;
    }


    ExpressionType FunctionDefinition::getResultType() const {
        if (body_expr) {
            return body_expr->getResultType();
        } else {
            return ExpressionType::unknown;
        }
    }
 
    void FunctionDefinition::pushArgStack(Arguments* st) {
        argstack.push_back(st);
    }

    Arguments* FunctionDefinition::topArgStack() {
        if (argstack.size()!=0) {
            return argstack.back();
        } else {
            return nullptr;
        }
    }

    void FunctionDefinition::popArgStack() {
        argstack.pop_back();
    }

    //template Expression<double>::Ptr functionParameter<double>(FunctionDefinition::Ptr, const std::string&);

 
    /*-----------------------------------------------------------------*/
    /*                 FunctionEvaluation                              */
    /*-----------------------------------------------------------------*/
/*    
    template<typename T>
    FunctionEvaluation<T>::FunctionEvaluation(
                FunctionDefinition::Ptr definition
            ): N_aryExpression<T>::N_aryExpression(definition->getName(), definition->argtypes.size())
    {
        this->definition(definition);
        body_expr = definition->get_body_expression<T>(); 
        if (body_expr==nullptr) {
            throw WrongTypeException();
        }
    }

    template<typename T> template<typename ArgType>
    int FunctionEvaluation<T>::addTypeCheckedArgument(const std::string& name, boost::shared_ptr<KDL::Expression<ArgType>> arg) {
        // find argument name:
        auto p = definition->name2ndx.find(name);
        if (p==definition->name2ndx.end()) {
            return -1;  
        }
        addTypeCheckedArgument(p->second, cached<ArgType>(arg));
    }

    template<typename T>
    int FunctionEvaluation<T>::addTypeCheckedArgument(int idx, ExpressionBase::Ptr arg) {
        // check range of idx:
        if ((idx<0)||(idx>=definition->getNrOfParam())) return -2;
        // type check
        if (arg->getResultType()!=definition->getParamType(idx)) return -3;
        // call superclass:
        return N_aryExpression<T>::setArgument(idx,arg);
    }
    
    template<typename T>
    typename FunctionEvaluation<T>::ResultType FunctionEvaluation<T>::value() {
        definition->pushArgStack(&N_aryExpression<T>::arguments);     
        body_expr->update_variabletype_from_original(); // invalidate all cache
        ResultType result = body_expr->value();
        definition->popArgStack();
        return result;
    }

    template<typename T>
    typename FunctionEvaluation<T>::DerivType FunctionEvaluation<T>::derivative(int i) {
        definition->pushArgStack(&N_aryExpression<T>::arguments);
        body_expr->update_variabletype_from_original(); // invalidate all cache
        ResultType val = body_expr->value();  // need because another function call can get in between
        DerivType  der = body_expr->derivative(i);
        definition->popArgStack();
        return der;
    }

    template<typename T>
    typename FunctionEvaluation<T>::DerivExprType FunctionEvaluation<T>::derivativeExpression(int i) {
    }


    template<typename T>
    typename Expression<T>::Ptr FunctionEvaluation<T>::clone() {    
        // we do not clone the function definition itself.
        typename Expression<T>::Ptr expr = boost::make_shared< FunctionEvaluation<T> >(definition);
        for (int i=0;i<this->nrOfArguments();++i) {
            expr->addTypeCheckedArgument(i, this->getArgument(i)->clone() ); 
        }
        return expr;
    }

    template<typename T>
    FunctionEvaluation<T>::~FunctionEvaluation() {
    }
    */
} // namespace KDL

