#ifndef EXPRESSIONGRAPH_FUNCTION_51566a22_0d16_4ad0_b673_fa9a6d3cecec
#define EXPRESSIONGRAPH_FUNCTION_51566a22_0d16_4ad0_b673_fa9a6d3cecec

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

#include "expressiontree_expressions.hpp"
//#include "expressiontree_n_ary.hpp"
#include <string> 
#include <boost/variant.hpp>
#include <unordered_map>
#include <map>
#include <vector>

namespace KDL {

    typedef std::vector<ExpressionBase::Ptr> Arguments;

    class FunctionDefinition {
            std::string name;
            ExpressionBase::Ptr         body_expr;   // pointer to an expression for the body of the function
            std::vector<std::string>   argnames;    // argument names
            std::vector<ExpressionType>  argtypes;    // argument variable types
            std::vector< Arguments*>   argstack;    // stack of pointers to the current value of the arguments. 
            std::unordered_map<std::string, int> name2ndx;
            Arguments default_arguments;
        public:
            typedef boost::shared_ptr<FunctionDefinition > Ptr;

            FunctionDefinition(const std::string& name);
            virtual void setName(const std::string& name);
            virtual std::string getName() const;

            virtual void addParam(const std::string& name, ExpressionType vartype); 
            virtual int getParamNdx(const std::string& name) const;
            virtual std::string getParamName(int i) const;
            virtual ExpressionType getParamType(int i) const;
            virtual int getNrOfParam() const;
           
            virtual void pushArgStack(Arguments* st);
            virtual Arguments* topArgStack();
            virtual void popArgStack(); 
     
            template <typename T>
            void setBodyExpression( typename Expression<T>::Ptr _body_expr) {
                body_expr = _body_expr;
            }
 

            virtual ExpressionType getResultType() const;            

            template <typename T>
            typename Expression<T>::Ptr getBodyExpression() {
                if (AutoDiffTrait<T>::expr_type != body_expr->getResultType()) {
                    return nullptr;
                } else {
                    return boost::dynamic_pointer_cast< Expression<T> >( body_expr);
                }
            }
            virtual ~FunctionDefinition(){}
    };

    inline typename FunctionDefinition::Ptr make_FunctionDefinition(const std::string& name ) {
        return boost::make_shared< FunctionDefinition>(name);
    } 
    
    std::ostream& operator<< (std::ostream& os, const FunctionDefinition& arg); 

    template<typename T>
    class FunctionParameter :
        public FunctionType<T> 
    {
            FunctionDefinition::Ptr definition;
            int                     index;
        public:
           typedef T ResultType;
           typedef typename AutoDiffTrait<T>::DerivType DerivType;
           /**
            * i: index of the parameters, assumed to be always valid (checked and handled by
            * the function that creates this object).
            */
           FunctionParameter( FunctionDefinition::Ptr  _definition, int i):
             FunctionType<T>(_definition->getParamName(i)),
             definition(_definition),
             index(i)  
            {
                // type check:
                if ( AutoDiffTrait<T>::expr_type!=definition->getParamType(i) ) {
                    throw ArgumentWrongTypeException();
                }
           }
           virtual ResultType value() {
                #ifdef EG_LOG_CACHE
                    std::cerr << "FunctionParameter " << this->name << ":value() entered" << std::endl;
                #endif
                typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                EG_ASSERT(e!=nullptr);
                ResultType retval = e->value();
                #ifdef EG_LOG_CACHE
                    std::cerr << "FunctionParameter " << this->name << ":value() =" << retval << std::endl;
                    std::cerr << "FunctionParameter " << this->name << ":value() returned" << std::endl;
                #endif
                return retval;
           }

           virtual DerivType derivative(int i) {
                #ifdef EG_LOG_CACHE
                    std::cerr << "FunctionParameter " << this->name << ":derivative("<<i<<") entered" << std::endl;
                #endif
                typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                EG_ASSERT(e!=nullptr);
                DerivType retval = e->derivative(i);
                #ifdef EG_LOG_CACHE
                    std::cerr << "FunctionParameter " << this->name << ":value() =" << retval << std::endl;
                    std::cerr << "FunctionParameter " << this->name << ":derivative("<<i<<") returned" << std::endl;
                #endif
                return retval;
          }
            virtual bool isConstant() const {
                return false;
            }
 
            virtual void getDependencies(std::set<int>& varset) {
                typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                EG_ASSERT(e!=nullptr);
                e->getDependencies(varset);
                /*if (varset.size()==0) {
                    varset.insert(0); // force a dependency on var. 0 (time) to avoid optimizing the expression away during construction
                                      // of function evaluation object
                }*/
            }

            virtual void getScalarDependencies(std::set<int>& varset) {
                typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                EG_ASSERT(e!=nullptr);
                e->getScalarDependencies(varset);
            }

            virtual void getRotDependencies(std::set<int>& varset) {
                typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                EG_ASSERT(e!=nullptr);
                e->getRotDependencies(varset);
            }

           virtual void setInputValues(const std::vector<double>& values) {}
           virtual void setInputValue(int variable_number, double val) {}
           virtual void setInputValue(int variable_number, const Rotation& val){}
           virtual void update_variabletype_from_original(){}

           virtual int number_of_derivatives(){
                typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                EG_ASSERT( e  );
                return e->number_of_derivatives();
            } 
            virtual void resize_nr_of_derivatives() {
               typename Expression<T>::Ptr e = boost::dynamic_pointer_cast< Expression<T> > (   (*definition->topArgStack())[index] );
                e->resize_nr_of_derivatives();
                #ifdef EG_LOG_CACHE
                    std::cerr << "FunctionParameter " << this->name << " resize_nr_of_derivatives to " <<  e->number_of_derivatives() << " of ";
                    e->print(std::cerr); std::cerr << e << std::endl; 
                #endif
            }
  
           virtual typename Expression<T>::Ptr clone() {
                return  boost::make_shared< FunctionParameter<T> >(this->definition, this->index);
           }   

           ~FunctionParameter() {} 
    };
   
    /**
     * returns nullpointer if parameter name does not exist for this definition.
     *
    template <typename T>
    typename Expression<T>::Ptr functionParameter( 
        typename FunctionDefinition::Ptr definition, 
        const std::string& name);*/
 
    /**
     * returns nullpointer if parameter name does not exist for this definition.
     *
    template<typename T>
    typename Expression<T>::Ptr functionParameter( 
        FunctionDefinition::Ptr definition, 
        const std::string& name) 
    {
        int ndx=definition->getParamNdx(name);
        if (ndx<0) return nullptr;
        return boost::make_shared< FunctionParameter<T> >( definition, ndx );
    }
    */

    template<typename T>
    class FunctionEvaluation: public Expression<T> {

        public:

            using Ptr            = boost::shared_ptr< FunctionEvaluation<T> >;
            using ResultType     = T;
            using DerivType      = typename AutoDiffTrait<T>::DerivType;
            using DerivExprType  = typename Expression< DerivType >::Ptr;

        private:
            int argcount; 
            bool invalidated;    
            ResultType                         cached_value;
            std::map<int, DerivType>           cached_derivatives;
            Arguments                          arguments;
            FunctionDefinition::Ptr            definition;
            boost::shared_ptr< Expression<T> > body_expr;   // typed copy of definition->body_expr


            // assuming valid ndx:
            void addArg(int ndx, boost::shared_ptr<ExpressionBase> arg ) {
                // type check
                if (arg->getResultType()!=definition->getParamType(ndx)) {
                    throw ArgumentWrongTypeException();
                }
                if (!arguments[ndx])  argcount++;   // if not already added, increase the argcount
                arguments[ndx] = arg;               // if already specified, overwrite
                if (argcount==definition->getNrOfParam()) {
                    // number of derivatives and dependencies may have changed:
                    std::set<int> varset;
                    definition->pushArgStack( &arguments );
                    body_expr->resize_nr_of_derivatives();
                    body_expr->getDependencies(varset);
                    definition->popArgStack();
                    cached_derivatives.clear();
                        #ifdef EG_LOG_CACHE
                        std::cout << "cached_derivative indices : ";
                        #endif
                    for (auto el: varset) {
                            #ifdef EG_LOG_CACHE
                            std::cout << el << "\t";
                            #endif
                        cached_derivatives[el] = AutoDiffTrait<DerivType>::zeroValue(); 
                    }
                        #ifdef EG_LOG_CACHE
                        std::cout << std::endl;
                        #endif
                }
                invalidated=true;
            }
            void ensure_computed() {
                if (invalidated) {
                    if (argcount!=definition->getNrOfParam()) {
                        throw WrongNumberOfArgumentsException(); 
                    }
                    definition->pushArgStack( &arguments );     
                    body_expr->update_variabletype_from_original(); // invalidate all cache
                    cached_value = body_expr->value();
                        #ifdef EG_LOG_CACHE
                        std::cout << "cached_derivative indices during ensure_compute : ";
                        #endif
                    for (auto &pair: cached_derivatives) {
                             #ifdef EG_LOG_CACHE
                             std::cout << pair.first << "\t";
                             #endif
                         pair.second = body_expr->derivative(pair.first);
                             #ifdef EG_LOG_CACHE
                             std::cout << "cached_derivative[" << pair.first << "]=" << pair.second << "\t\t";
                             #endif
                    }
                        #ifdef EG_LOG_CACHE
                        std::cout << std::endl;
                        #endif
                    definition->popArgStack();
                    invalidated=false;    
                } else {
                }
            } 

        public:
            
            explicit FunctionEvaluation( 
                FunctionDefinition::Ptr _definition
            ): Expression<T>(_definition->getName()),
               definition( _definition ),
               arguments(  _definition->getNrOfParam() ),
               invalidated(true)
            {
                body_expr = _definition->getBodyExpression<T>(); 
                if (body_expr==nullptr) {
                    throw BodyWrongTypeException();
                }   
                argcount = 0;
            }

            FunctionEvaluation(FunctionDefinition::Ptr _definition, std::initializer_list<ExpressionBase::Ptr> list):
                 Expression<T>(_definition->getName()),
                 definition( _definition ),
                 arguments(  _definition->getNrOfParam() ),
                 invalidated(true)
            {
                body_expr = _definition->getBodyExpression<T>(); 
                if (body_expr==nullptr) {
                    throw BodyWrongTypeException();
                }
                argcount = 0;
                int idx=0;
                if (list.size() != definition->getNrOfParam()) {
                    throw WrongNumberOfArgumentsException(); 
                }
                for (auto it = list.begin(); it != list.end(); ++it) {
                   addArg(idx, *it);
                   idx++;
                }
            }
            template<typename ArgType>
            void addTypeCheckedArgument(
                const std::string& name, 
                boost::shared_ptr<Expression<ArgType>> arg
            ) {
                int ndx = definition->getParamNdx(name);
                if (ndx==-1) {
                    throw ArgumentNameException();
                }
                arg = cached<ArgType>("argument "+name,arg);
                addArg(ndx,arg);
             }

            ResultType value() {
                ensure_computed();
                return cached_value; 
            }

            DerivType derivative(int i) {
                auto p = cached_derivatives.find(i);
                if (p!=cached_derivatives.end()) {
                        #ifdef EG_LOG_CACHE
                        std::cout << "derivative("<<i<<") ensure_compute used" << std::endl;
                        #endif
                    ensure_computed();
                    return cached_derivatives[i];
                } else {
                        #ifdef EG_LOG_CACHE
                        std::cout << "derivative("<<i<<") ensure_compute NOT used" << std::endl;
                        #endif
                    return AutoDiffTrait<DerivType>::zeroValue();
                }
            }

            virtual DerivExprType derivativeExpression(int i) {
                throw NotImplementedException();
            }

            virtual typename Expression<T>::Ptr clone() {
                // we do not clone the function definition itself.
                auto expr = boost::make_shared< FunctionEvaluation<T> >(definition);
                for (int i=0;i<arguments.size();++i) {
                    expr->addArg(i, arguments[i]->cloneBase() ); 
                }
                return expr;
            }

            virtual void setInputValues(const std::vector<double>& values) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->setInputValues(values); 
                }
                body_expr->setInputValues(values);
                invalidated=true;
            }

            virtual void setInputValue(int variable_number, double val) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->setInputValue(variable_number,val); 
                }
                body_expr->setInputValue(variable_number,val);
                invalidated=true;
            }

            virtual void setInputValue(int variable_number, const Rotation& val) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->setInputValue(variable_number,val); 
                }
                body_expr->setInputValue(variable_number,val);
                invalidated=true;
            }

            virtual int number_of_derivatives() {
                int n=0;
                for (unsigned int i=0;i<arguments.size();++i) {
                    n = std::max( n, arguments[i]->number_of_derivatives() ); 
                }
                n = std::max( n, body_expr->number_of_derivatives() ); 
                return n;
            } 

            virtual void resize_nr_of_derivatives() {
                #ifdef EG_LOG_CACHE
                    std::cerr << "FunctionEvaluation " << this->name << " resize_nr_of_derivatives" << std::endl; 
                #endif
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->resize_nr_of_derivatives(); 
                }
                body_expr->resize_nr_of_derivatives();
                invalidated=true;
            }

            virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
                typename Expression<Frame>::Ptr a;
                for (unsigned int i=0;i<arguments.size();++i) {
                    a = arguments[i]->subExpression_Frame(name);
                    if (a) {
                        return a;
                    }
                }
                return a;
            }

            virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
                typename Expression<Rotation>::Ptr a;
                for (unsigned int i=0;i<arguments.size();++i) {
                    a = arguments[i]->subExpression_Rotation(name);
                    if (a) {
                        return a;
                    }
                }
                return a;
            }

            virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
                typename Expression<Vector>::Ptr a;
                for (unsigned int i=0;i<arguments.size();++i) {
                    a = arguments[i]->subExpression_Vector(name);
                    if (a) {
                        return a;
                    }
                }
                return a;
            }

            virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
                typename Expression<Wrench>::Ptr a;
                for (unsigned int i=0;i<arguments.size();++i) {
                    a = arguments[i]->subExpression_Wrench(name);
                    if (a) {
                        return a;
                    }
                }
                return a;
            }

            virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
                typename Expression<Twist>::Ptr a;
                for (unsigned int i=0;i<arguments.size();++i) {
                    a = arguments[i]->subExpression_Twist(name);
                    if (a) {
                        return a;
                    }
                }
                return a;

            }

            virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
                typename Expression<double>::Ptr a;
                for (unsigned int i=0;i<arguments.size();++i) {
                    a = arguments[i]->subExpression_Double(name);
                    if (a) {
                        return a;
                    }
                }
                return a;
            }

            // the optimizer only knows the arguments, the function itself is always recomputed.
            virtual void addToOptimizer(ExpressionOptimizer& opt) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->addToOptimizer(opt);
                }
            }
            virtual bool isConstant() const {
                return false;
            }
 

            virtual void getDependencies(std::set<int>& varset) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->getDependencies(varset);
                }
                body_expr->getDependencies(varset);
            }

            virtual void getScalarDependencies(std::set<int>& varset) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->getScalarDependencies(varset);
                }
                body_expr->getScalarDependencies(varset);
            }

            virtual void getRotDependencies(std::set<int>& varset) {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->getRotDependencies(varset);
                }
                body_expr->getRotDependencies(varset);
            }

            virtual void debug_printtree() {
                std::cout << Expression<ResultType>::name << "(";
                for (unsigned int i=0;i<arguments.size();++i) {
                    if (i!=0) std::cout << ",";
                    arguments[i]->debug_printtree();
                }
                std::cout << ")";
            }
            virtual void print(std::ostream& os) const {
                os << "" << Expression<ResultType>::name << "(";
                for (unsigned int i=0;i<arguments.size();++i) {
                    if (i!=0) os << ",";
                    arguments[i]->print(os);
                }
                os << ")";
            }

            virtual void update_variabletype_from_original() {
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->update_variabletype_from_original();
                }
                body_expr->update_variabletype_from_original(); 
                invalidated=true;
            }

            virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
                thisnode=(size_t)this;
                of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
                   << COLOR_OPERATION << ",color=black]\n";
                std::vector<pnumber> argnode( arguments.size());
                for (unsigned int i=0;i<arguments.size();++i) {
                    arguments[i]->write_dotfile_update(of,argnode[i]);
                }
                for (unsigned int i=0;i<arguments.size();++i) {
                    of << "S"<<thisnode<<" -> " << "S"<<argnode[i] << "\n"; //rev
                }
            }
            virtual void write_dotfile_init() {
                for (unsigned int i =0;i<arguments.size();++i) {
                    arguments[i]->write_dotfile_init();
                }
            } 

            ~FunctionEvaluation() {
            }
    };
 
    template<typename T>
    typename FunctionParameter<T>::Ptr make_FunctionParameter(FunctionDefinition::Ptr definition, const std::string& name ) {
        definition->addParam( name, AutoDiffTrait<T>::expr_type);
        int ndx=definition->getParamNdx(name);
        return boost::make_shared< FunctionParameter<T> >( definition, ndx );
    } 
  

    template<typename T>
    typename FunctionEvaluation<double>::Ptr make_FunctionEvaluation( 
                FunctionDefinition::Ptr _definition
            ) {
       FunctionEvaluation<double>::Ptr e = boost::make_shared<FunctionEvaluation<double>>(_definition);
       return e;
    }
      

    template<typename T>
    typename FunctionEvaluation<T>::Ptr make_FunctionEvaluation( 
                FunctionDefinition::Ptr _definition,
                std::initializer_list<ExpressionBase::Ptr> list
            ) {
       return boost::make_shared<FunctionEvaluation<T>>(_definition, list);
    }
      

} // namespace KDL
#endif
