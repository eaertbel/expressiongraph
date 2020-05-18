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
            double maxstepsize; 
            int maxRecDepth; 
            bool invalidated;    
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

            void ensure_computed();

        public:

            IntegralNode( 
                FunctionDefinition::Ptr _integrand, 
                Expression<double>::Ptr _lower, 
                Expression<double>::Ptr _upper, 
                double _epsilon, 
                double _maxstepsize, 
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
/*
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
             ~FunctionEvaluation();
*/ 
            virtual ~IntegralNode();
    };

    void IntegralNode::ensure_computed() {
        if (invalidated) {
            for (auto el: cached_derivatives) {
                el.second = 0.0;
            }
            // Integral(f,a,b)->deriv(i) = Integral( f->deriv(i),a,b) + f( b)*b->deriv(i) - f(a)*a->deriv(i) 
            double lowerval = lower->value();
            double upperval = upper->value();
            int errcode =  integrator->integrate(func,lowerval,upperval,result);
            cached_value = result[0];
            for (int i=0;i<index_map.size();++i) {
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
        double _maxstepsize, 
        int _maxRecDepth ):
            integrand(_integrand),
            lower(_lower),
            upper(_upper),
            epsilon(_epsilon),
            maxstepsize(_maxstepsize),
            maxRecDepth(_maxRecDepth)
     {
        if (integrand->getNrOfParam()!=1) {
            throw WrongNumberOfArgumentsException();
        }
        // getting dependency information and establish 
        lower->getDependencies(set_lower);
        upper->getDependencies(set_upper);
        integrand->getBodyExpression<double>()->getDependencies(set_integrand); 
        for (auto &&e : set_lower) {     
            cached_derivatives[e] = 0.0;
        }
        for (auto &&e : set_upper) {
            cached_derivatives[e] = 0.0;
        }
        int ndx=1;
        index_map.resize( set_integrand.size() );
        for (auto &&e : set_integrand) {
            cached_derivatives[e] = 0.0;
            index_map[ndx] = e;
            ndx = ndx + 1;
        }
        cached_value=0.0;
        invalidated=true;
        integrator = boost::make_shared< IntegratorAdaptiveSimpson >( 
            epsilon, 
            maxstepsize, 
            maxRecDepth, 
            set_integrand.size()+1
        );

        var = Variable<double>({});
        func_eval.reset( new FunctionEvaluation<double>(integrand,{var}));
        int vecsize = set_integrand.size()+1;
        result = new double[vecsize];
        func =  [=](double arg,double* result) {
            var->setValue(arg);         
            result[0] = func_eval->value();
            int i=1;
            for (auto&& ndx : set_integrand) {
                result[i] = func_eval->derivative(ndx);
                i = i + 1;
            }
            return;
        };
    }

    double IntegralNode::value() {
        ensure_computed();
        return cached_value; 
    }

    double IntegralNode::derivative(int i) {
        auto p = cached_derivatives.find(i);
        if (p!=cached_derivatives.end()) {
            ensure_computed();
            return cached_derivatives[i];
        } else {
            return 0.0;
        }
    }

    Expression<double>::Ptr IntegralNode::derivativeExpression(int i) {
    }

    Expression<double>::Ptr IntegralNode::clone() {
        return boost::make_shared< IntegralNode >(integrand, lower,upper,epsilon,maxstepsize,maxRecDepth);
    }

    void IntegralNode::setInputValues(const std::vector<double>& values) {
        assert(0 && "NOT SUPPORTED");
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
        return std::max( { 
                func_eval->number_of_derivatives(),
                lower->number_of_derivatives(),
                upper->number_of_derivatives()
        });
    } 

    void IntegralNode::resize_nr_of_derivatives() {
        func_eval->resize_nr_of_derivatives();
        lower->number_of_derivatives();
        upper->number_of_derivatives();
    }
    void IntegralNode::update_variabletype_from_original() {
        func_eval->update_variabletype_from_original();
        lower->update_variabletype_from_original();
        upper->update_variabletype_from_original();
        invalidated=true;
    }

/*
    typename Expression<Frame>::Ptr IntegralNode::subExpression_Frame(const std::string& name) {
        typename Expression<Frame>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Frame(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    typename Expression<Rotation>::Ptr IntegralNode::subExpression_Rotation(const std::string& name) {
        typename Expression<Rotation>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Rotation(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    typename Expression<Vector>::Ptr IntegralNode::subExpression_Vector(const std::string& name) {
        typename Expression<Vector>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Vector(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    typename Expression<Wrench>::Ptr IntegralNode::subExpression_Wrench(const std::string& name) {
        typename Expression<Wrench>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Wrench(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    typename Expression<Twist>::Ptr IntegralNode::subExpression_Twist(const std::string& name) {
        typename Expression<Twist>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Twist(name);
            if (a) {
                return a;
            }
        }
        return a;

    }

    typename Expression<double>::Ptr IntegralNode::subExpression_Double(const std::string& name) {
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
    void IntegralNode::addToOptimizer(ExpressionOptimizer& opt) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->addToOptimizer(opt);
        }
    }


    void IntegralNode::getDependencies(std::set<int>& varset) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->getDependencies(varset);
        }
        body_expr->getDependencies(varset);
    }

    void IntegralNode::getScalarDependencies(std::set<int>& varset) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->getScalarDependencies(varset);
        }
        body_expr->getScalarDependencies(varset);
    }

    void IntegralNode::getRotDependencies(std::set<int>& varset) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->getRotDependencies(varset);
        }
        body_expr->getRotDependencies(varset);
    }

    void IntegralNode::debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        for (unsigned int i=0;i<arguments.size();++i) {
            if (i!=0) std::cout << ",";
            arguments[i]->debug_printtree();
        }
        std::cout << ")";
    }
    void IntegralNode::print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        for (unsigned int i=0;i<arguments.size();++i) {
            if (i!=0) os << ",";
            arguments[i]->print(os);
        }
        os << ")";
    }

    void IntegralNode::write_dotfile_update(std::ostream& of, pnumber& thisnode) {
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
    void IntegralNode::write_dotfile_init() {
        for (unsigned int i =0;i<arguments.size();++i) {
            arguments[i]->write_dotfile_init();
        }
    } 

    IntegralNode::~FunctionEvaluation() {
    }
*/ 
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
                double _maxstepsize, 
                int _maxRecDepth ) {
        return boost::make_shared<IntegralNode>(
            _integrand,_lower,_upper,_epsilon,_maxstepsize,_maxRecDepth
        );
    }



} // namespace


