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


#include <expressiongraph/expressiontree_expressions.hpp>
#include <expressiongraph/expressiontree_function.hpp>
#include <expressiongraph/expressiontree_integral.hpp>
#include <expressiongraph/integrator.hpp>

#define EG_LOG_INTEGRALNODE
namespace KDL {

    inline bool contains(const std::set<int>& s, const std::vector<int>& vec) {
        for (auto el : vec) {
            if (s.count(el)>0) return true;
        }
        return false;
    }

    class IntegralNode:public FunctionType<double> {
            FunctionDefinition::Ptr integrand;
            Expression<double>::Ptr lower; 
            Expression<double>::Ptr upper; 
            double epsilon; 
            int minRecDepth; 
            int maxRecDepth; 
            bool invalidated;    
            bool first_time;
            double                cached_value;
            std::map<int, double> cached_derivatives;
            std::vector<int>      index_map; // for result
            double*               result;    // contains result for integration.
            boost::shared_ptr<IntegratorAdaptiveSimpson> integrator;
            FunctionEvaluation<double>::Ptr func_eval;
            VariableType<double>::Ptr var;
            IntegrandType func;

            std::set<int> set_lower;
            std::set<int> set_upper;
            std::set<int> set_integrand;

            void ensure_initialized();
            void ensure_computed();

        public:

            IntegralNode( 
                FunctionDefinition::Ptr _integrand, 
                Expression<double>::Ptr _lower, 
                Expression<double>::Ptr _upper, 
                double _epsilon, 
                int _minRecDepth, 
                int _maxRecDepth );
            virtual double value();
            virtual double derivative(int i);
            virtual Expression<double>::Ptr derivativeExpression(int i);
            virtual Expression<double>::Ptr clone();
            virtual void setInputValues(const std::vector<double>& values);
            virtual void setInputValue(int variable_number, double val);
            virtual void setInputValue(int variable_number, const Rotation& val);
            virtual int number_of_derivatives();
            virtual void resize_nr_of_derivatives();
            virtual void update_variabletype_from_original();
            virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name);
            virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name);
            virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name);
            virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name);
            virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name);
            virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name);
            virtual void addToOptimizer(ExpressionOptimizer& opt);
            virtual bool isConstant() const;
            virtual void getDependencies(std::set<int>& varset);
            virtual void getScalarDependencies(std::set<int>& varset);
            virtual void getRotDependencies(std::set<int>& varset);
            virtual void debug_printtree();
            virtual void print(std::ostream& os) const;
            virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode);
            virtual void write_dotfile_init();
            virtual ~IntegralNode();
    };

    void IntegralNode::ensure_initialized() {
        if (first_time) {
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << "initialize() entered" << std::endl;
            #endif
            // getting dependency information 
            lower->getDependencies(set_lower);
            upper->getDependencies(set_upper);
            integrand->getBodyExpression<double>()->getDependencies(set_integrand); 
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << "dependency of lower limit : ";
            #endif
            for (auto &&e : set_lower) {     
                cached_derivatives[e] = 0.0;
                #ifdef EG_LOG_INTEGRALNODE
                std::cout << e << " " <<std::endl;
                #endif
            }
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << std::endl;
            std::cout << "dependency of upper limit : ";
            #endif
            for (auto &&e : set_upper) {
                cached_derivatives[e] = 0.0;
                #ifdef EG_LOG_INTEGRALNODE
                std::cout << e << " " <<std::endl;
                #endif
            }
            int ndx=1;
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << std::endl;
            std::cout << "dependency of integrand : " << std::endl;
            #endif
            index_map.resize( set_integrand.size()+1 );
            for (auto &&e : set_integrand) {
                #ifdef EG_LOG_INTEGRALNODE
                std::cout << "         index_map[ndx] with ndx=" << ndx << " and index_map[ndx]= "<<e << std::endl;
                #endif
                cached_derivatives[e] = 0.0;
                index_map[ndx] = e;
                ndx = ndx + 1;
            }
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << std::endl;
            std::cout << "Is integrand constant() " << integrand->getBodyExpression<double>()->isConstant() << std::endl;
            #endif
            cached_value=0.0;
            invalidated=true;
            integrator = boost::make_shared< IntegratorAdaptiveSimpson >( 
                epsilon, 
                minRecDepth, 
                maxRecDepth, 
                set_integrand.size()+1
            );

            var = Variable<double>({});
            func_eval.reset( new FunctionEvaluation<double>(integrand,{var}));
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << "Function evaluation node "<< std::endl;
            std::cout << "    " << func_eval << std::endl;
            std::cout << "    isConstant() " << func_eval->isConstant() << std::endl;
            std::set<int> fset;
            func_eval->getDependencies(fset); 
            std::cout << "    dependencies : "; 
            for (auto && e: fset) {
                std::cout << e<<" ";
            }
            std::cout << std::endl;
            #endif
            int vecsize = set_integrand.size()+1;
            result = new double[vecsize];
            func =  [=](double arg,double* result)->void {
                var->setValue(arg);         
                func_eval->update_variabletype_from_original();
                result[0] = func_eval->value();
                #ifdef EG_LOG_INTEGRALNODE
                    std::cout << "evaluate at " << arg << " resulting in " << result[0];
                #endif
                int i=1;
                for (auto&& ndx : set_integrand) {
                    result[i] = func_eval->derivative(ndx);
                    #ifdef EG_LOG_INTEGRALNODE
                    std::cout << "\t("<<ndx<<")" << result[i];
                    #endif
                    i = i + 1;
                }
                #ifdef EG_LOG_INTEGRALNODE
                std::cout << std::endl;
                #endif
                return;
            };
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << "initialize() exited" << std::endl;
            #endif
            first_time=false;
        }
    }
    void IntegralNode::ensure_computed() {
        if (invalidated) {
            ensure_initialized();
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << "recompute value and derivatives :";
            #endif
            for (auto el: cached_derivatives) {
                el.second = 0.0;
            }
            // Integral(f,a,b)->deriv(i) = Integral( f->deriv(i),a,b) + f( b)*b->deriv(i) - f(a)*a->deriv(i) 
            double lowerval = lower->value();
            double upperval = upper->value();
            int errcode =  integrator->integrate(func,lowerval,upperval,result);
            cached_value = result[0];
            for (auto &e: cached_derivatives) {
                e.second = 0.0;
                #ifdef EG_LOG_INTEGRALNODE
                std::cout << e.first;
                #endif
            }
            #ifdef EG_LOG_INTEGRALNODE
            std::cout << std::endl;
            std::cout << "returned result : ";
            for (int i=0;i<set_integrand.size()+1;++i) {
                std::cout << result[i] <<"\t";
            }
            std::cout << std::endl;
            #endif
            for (int i=1;i<index_map.size()+1;++i) {
                #ifdef EG_LOG_INTEGRALNODE
                std::cout << "         index_map[i] with i=" << i << " and index_map[i]= "<< index_map[i] << std::endl;
                #endif
                cached_derivatives[index_map[i]] = result[i];
            } 
            for (auto i: set_lower) {
                cached_derivatives[i] += -integrator->get_fa_value()*lower->derivative(i);
            }
            for (auto i: set_upper) {
                cached_derivatives[i] +=  integrator->get_fb_value()*upper->derivative(i);
            }
            invalidated=false;
        }
    }

    IntegralNode::IntegralNode( 
        FunctionDefinition::Ptr _integrand, 
        Expression<double>::Ptr _lower, 
        Expression<double>::Ptr _upper, 
        double _epsilon, 
        int _minRecDepth, 
        int _maxRecDepth ):
            integrand(_integrand),
            lower(_lower),
            upper(_upper),
            epsilon(_epsilon),
            minRecDepth(_minRecDepth),
            maxRecDepth(_maxRecDepth),
            first_time(true),
            invalidated(true)
     {
        #ifdef EG_LOG_INTEGRALNODE
        std::cout << "create Integral Node" << std::endl;
        std::cout << "    Integrand : ";
        integrand->getBodyExpression<double>()->print(std::cout);
        std::cout << std::endl;
        std::cout << "    Lower limit : ";
        lower->print(std::cout);
        std::cout << std::endl;
        std::cout << "    Upper limit : ";
        upper->print(std::cout);
        std::cout << std::endl;
        #endif 
        if (integrand->getNrOfParam()!=1) {
            throw WrongNumberOfArgumentsException();
        }
    }

    double IntegralNode::value() {
        ensure_computed();
        return cached_value; 
    }

    double IntegralNode::derivative(int i) {
        ensure_initialized();
        auto p = cached_derivatives.find(i);
        if (p!=cached_derivatives.end()) {
            ensure_computed();
            #ifdef EG_LOG_INTEGRALNODE
                std::cout << "derivative("<<i<<") cached value returned " <<  cached_derivatives[i] << std::endl;
            #endif
            return cached_derivatives[i];
        } else {
            #ifdef EG_LOG_INTEGRALNODE
                std::cout << "derivative("<<i<<") not a dependency zero returned" << std::endl;
            #endif
            return 0.0;
        }
    }

    Expression<double>::Ptr IntegralNode::derivativeExpression(int i) {
        throw NotImplementedException();
    }

    Expression<double>::Ptr IntegralNode::clone() {
        return boost::make_shared< IntegralNode >(integrand, lower,upper,epsilon,minRecDepth,maxRecDepth);
    }

    void IntegralNode::setInputValues(const std::vector<double>& values) {
        throw NotImplementedException();
        /*if (contains(set_integrand,values)) {
            func_eval->setInputValues(values);
            invalidated=true;
        }
        if (contains(set_lower,values)) {
            lower->setInputValues(values);
            invalidated=true;
        }
        if (contains(set_upper,values)) {
            upper->setInputValues(values);
            invalidated=true;
        }*/
    }

    void IntegralNode::setInputValue(int variable_number, double val) {
        ensure_initialized();
        if (set_integrand.count(variable_number)>0) {
            func_eval->setInputValue(variable_number,val);
            invalidated=true;
        }
        if (set_lower.count(variable_number)>0) {
            lower->setInputValue(variable_number,val);
            invalidated=true;
        }
        if (set_upper.count(variable_number)>0) {
            upper->setInputValue(variable_number,val);
            invalidated=true;
        }
    }

    void IntegralNode::setInputValue(int variable_number, const Rotation& val) {
        ensure_initialized();
        if ((set_integrand.count(variable_number)>0) ||
            (set_integrand.count(variable_number+1)>0) ||
            (set_integrand.count(variable_number+2)>0)
           ) {
            func_eval->setInputValue(variable_number,val);
            invalidated=true;
        }
        if ((set_lower.count(variable_number)>0) ||
            (set_lower.count(variable_number+1)>0) ||
            (set_lower.count(variable_number+2)>0)
           ) {
            lower->setInputValue(variable_number,val);
            invalidated=true;
        }
        if ((set_upper.count(variable_number)>0)   ||
            (set_upper.count(variable_number+1)>0) ||
            (set_upper.count(variable_number+2)>0)
           ) {
            upper->setInputValue(variable_number,val);
            invalidated=true;
        }
        func_eval->setInputValue(variable_number,val);
        lower->setInputValue(variable_number,val);
        upper->setInputValue(variable_number,val);
    }

    int IntegralNode::number_of_derivatives() {
        ensure_initialized();
        return std::max( { 
                func_eval->number_of_derivatives(),
                lower->number_of_derivatives(),
                upper->number_of_derivatives()
        });
    } 

    void IntegralNode::resize_nr_of_derivatives() {
        ensure_initialized();
        func_eval->resize_nr_of_derivatives();
        lower->number_of_derivatives();
        upper->number_of_derivatives();
    }
    void IntegralNode::update_variabletype_from_original() {
        ensure_initialized();
        func_eval->update_variabletype_from_original();
        lower->update_variabletype_from_original();
        upper->update_variabletype_from_original();
        invalidated=true;
    }

    bool IntegralNode::isConstant() const {
        return false;
    }
    void IntegralNode::getDependencies(std::set<int>& varset) {
        ensure_initialized();
        lower->getDependencies(varset);
        upper->getDependencies(varset);
        integrand->getBodyExpression<double>()->getDependencies(varset); 
    }

    void IntegralNode::getScalarDependencies(std::set<int>& varset) {
        ensure_initialized();
        lower->getDependencies(varset);
        upper->getDependencies(varset);
        integrand->getBodyExpression<double>()->getDependencies(varset); 
    }

    void IntegralNode::getRotDependencies(std::set<int>& varset) {
        ensure_initialized();
        lower->getDependencies(varset);
        upper->getDependencies(varset);
        integrand->getBodyExpression<double>()->getDependencies(varset); 
    }

    typename Expression<Frame>::Ptr IntegralNode::subExpression_Frame(const std::string& name) {
        throw NotImplementedException();
//        typename Expression<Frame>::Ptr a;
//        for (unsigned int i=0;i<arguments.size();++i) {
//            a = arguments[i]->subExpression_Frame(name);
//            if (a) {
//                return a;
//            }
//        }
//        return a;
    }

    typename Expression<Rotation>::Ptr IntegralNode::subExpression_Rotation(const std::string& name) {
        throw NotImplementedException();
//        typename Expression<Rotation>::Ptr a;
//        for (unsigned int i=0;i<arguments.size();++i) {
//            a = arguments[i]->subExpression_Rotation(name);
//            if (a) {
//                return a;
//            }
//        }
//        return a;
    }

    typename Expression<Vector>::Ptr IntegralNode::subExpression_Vector(const std::string& name) {
        throw NotImplementedException();
//        typename Expression<Vector>::Ptr a;
//        for (unsigned int i=0;i<arguments.size();++i) {
//            a = arguments[i]->subExpression_Vector(name);
//            if (a) {
//                return a;
//            }
//        }
//        return a;
    }

    typename Expression<Wrench>::Ptr IntegralNode::subExpression_Wrench(const std::string& name) {
        throw NotImplementedException();
//        typename Expression<Wrench>::Ptr a;
//        for (unsigned int i=0;i<arguments.size();++i) {
//            a = arguments[i]->subExpression_Wrench(name);
//            if (a) {
//                return a;
//            }
//        }
//        return a;
    }

    typename Expression<Twist>::Ptr IntegralNode::subExpression_Twist(const std::string& name) {
        throw NotImplementedException();
//        typename Expression<Twist>::Ptr a;
//        for (unsigned int i=0;i<arguments.size();++i) {
//            a = arguments[i]->subExpression_Twist(name);
//            if (a) {
//                return a;
//            }
//        }
//        return a;

    }

    typename Expression<double>::Ptr IntegralNode::subExpression_Double(const std::string& name) {
        throw NotImplementedException();
 //       typename Expression<double>::Ptr a;
 //       for (unsigned int i=0;i<arguments.size();++i) {
 //           a = arguments[i]->subExpression_Double(name);
 //           if (a) {
 //               return a;
 //           }
 //       }
 //       return a;
    }

    // the optimizer only knows the arguments, the function itself is always recomputed.
    void IntegralNode::addToOptimizer(ExpressionOptimizer& opt) {
        throw NotImplementedException();
//        for (unsigned int i=0;i<arguments.size();++i) {
//            arguments[i]->addToOptimizer(opt);
//        }
    }


    void IntegralNode::debug_printtree() {
        throw NotImplementedException();
        //std::cout << Expression<ResultType>::name << "(";
        //for (unsigned int i=0;i<arguments.size();++i) {
        //    if (i!=0) std::cout << ",";
        //    arguments[i]->debug_printtree();
       // }
        //std::cout << ")";
    }
    void IntegralNode::print(std::ostream& os) const {
        os << "Integral(lower=";lower->print(os);
        os << ", upper=";upper->print(os);
        os << ", integrand = "; func_eval->print(os);
        os << ")"; 
        //os << "" << Expression<ResultType>::name << "(";
        //for (unsigned int i=0;i<arguments.size();++i) {
        //    if (i!=0) os << ",";
        //    arguments[i]->print(os);
        //}
        //os << ")";
    }

    void IntegralNode::write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        throw NotImplementedException();
        //thisnode=(size_t)this;
        //of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
        //   << COLOR_OPERATION << ",color=black]\n";
        //std::vector<pnumber> argnode( arguments.size());
        //for (unsigned int i=0;i<arguments.size();++i) {
        //    arguments[i]->write_dotfile_update(of,argnode[i]);
       // }
        //for (unsigned int i=0;i<arguments.size();++i) {
         //   of << "S"<<thisnode<<" -> " << "S"<<argnode[i] << "\n"; //rev
        //}
    }
    void IntegralNode::write_dotfile_init() {
        throw NotImplementedException();
        //for (unsigned int i =0;i<arguments.size();++i) {
        //    arguments[i]->write_dotfile_init();
        //}
    } 
    IntegralNode::~IntegralNode(){
        delete [] result;
    }



    /**
     *
     * Integral(f,a,b)->deriv(i) = Integral( f->deriv(i),a,b) + f( b)*b->deriv(i) - f(a)*a->deriv(i) 
     *
     * bijv.   integral(0,s, A*sin(t) ) = A(1 - cos(s)) = A*sin(s)*s->deriv(i) + (1-cos(s))*A->deriv(i)
     * deriv s: 
     *       Integral(0,s, A->deriv(i)*sin(t)) + A*sin(s)*s->deriv(i) 
     *
     * Recognizing common patterns?
     *    Explicitly maintain the dependencies of the integrand, lower limit and upper limit.
     * TERMINOLOGY:
     *    \int_a^b { f(s) \mathrm{d} s}  is called the integral
     *    integrand: f(s)
     *    variable of integration: s
     *    a and b : the limits of the integral: lower limit and upper limit.
     *
     * Integral routine:
     *   - first argument and only argument is the variable of integration
     *   - a vector containing everything that needs to be integrated:
     *        - from the value of the integrad
     *        - from the derivaives of the integrand
     *
     * Integral(f,a,b)->deriv(i) = Integral( f->deriv(i),a,b) + f( b)*b->deriv(i) - f(a)*a->deriv(i)  
     */
    Expression<double>::Ptr make_Integral( 
                FunctionDefinition::Ptr _integrand, 
                Expression<double>::Ptr _lower, 
                Expression<double>::Ptr _upper, 
                double _epsilon, 
                int _minRecDepth, 
                int _maxRecDepth ) {
        return boost::make_shared<IntegralNode>(
            _integrand,_lower,_upper,_epsilon,_minRecDepth,_maxRecDepth
        );
    }



} // namespace


