#ifndef EXPRESSIONGRAPH_TF_ITASC_OUTPUT
#define EXPRESSIONGRAPH_TF_ITASC_OUTPUT

/*
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

#include <kdl/expressiontree_expressions.hpp>

namespace KDL { 

/**
 * Implements a MIMO (multiple input, multiple output) expressiongraph node that is always cached.
 * This node deals correctly with the expression optimizer.
 * \caveat If you override one of the methods, be sure to call the MIMO methods also
 * \caveat the overriding class is responsible for filling inputDouble/inputFrame/inputTwist with the correct input expression graphs.
 */
class MIMO: public CachedExpression {
    std::vector<boost::weak_ptr<MIMO> >  queue_of_clones;
    bool                                          dot_already_written;
public:

    typedef boost::shared_ptr<MIMO> Ptr;

    const std::string                             name;
    std::vector< Expression<double>::Ptr >        inputDouble;
    std::vector< Expression<Frame>::Ptr >         inputFrame;
    std::vector< Expression<Twist>::Ptr >         inputTwist;
    bool                                          cached;

    MIMO();

    MIMO( const std::string& name);
    
    virtual void setInputValues(const std::vector<double>& values);
    virtual void setInputValue(int variable_number, double val);
    virtual void setInputValue(int variable_number, const Rotation& val);
    
    virtual int number_of_derivatives();

    virtual Expression<Frame>::Ptr subExpression_Frame(const std::string& name);
    virtual  Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name);
    virtual  Expression<Vector>::Ptr subExpression_Vector(const std::string& name);
    virtual  Expression<Twist>::Ptr subExpression_Twist(const std::string& name);
    virtual  Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name);
    virtual  Expression<double>::Ptr subExpression_Double(const std::string& name);


    virtual void addToOptimizer(ExpressionOptimizer& opt);

    virtual void getDependencies(std::set<int>& varset);
    virtual void getScalarDependencies(std::set<int>& varset);
    virtual void getRotDependencies(std::set<int>& varset);
    
    virtual void invalidate_cache();

    virtual void debug_printtree();

    virtual void print(std::ostream& os) const;

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode);
    virtual void write_dotfile_init();

    // to be overwritten by somebody who derives from MIMO
    virtual MIMO::Ptr clone() = 0;


    // only to be called by MIMO_Output, typically not overridden
    virtual MIMO::Ptr getClone(int count);

    virtual ~MIMO();
};


/**
 * \caveat if you override one of the methods, be sure to call the MIMO_Output in the overridden method.
 */
template<class ResultType>
class MIMO_Output: public Expression<ResultType> {
    int                nr_of_clones;
public:
    MIMO::Ptr          mimo;

    MIMO_Output() {}
    MIMO_Output(const std::string& name, MIMO::Ptr _mimo):
                    Expression<double>(name),
                    nr_of_clones(-1),
                    mimo(_mimo)
                {}
 
    virtual void setInputValues(const std::vector<double>& values) {
        mimo->setInputValues(values);
    }

    virtual void setInputValue(int variable_number, double val) {
        mimo->setInputValue(variable_number,val);
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        mimo->setInputValue(variable_number,val);
    }

    virtual int number_of_derivatives() {
        return mimo->number_of_derivatives();
    } 

    virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        return mimo->subExpression_Frame(name);
    }

    virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        return mimo->subExpression_Rotation(name);
    }

    virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        return mimo->subExpression_Vector(name);
    }

    virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        return mimo->subExpression_Twist(name);
    }
    virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        return mimo->subExpression_Wrench(name);
    }
    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        return mimo->subExpression_Double(name);
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        mimo->addToOptimizer(opt);
    }

    virtual void getDependencies(std::set<int>& varset) {
        mimo->getDependencies(varset);
    }

    virtual void getScalarDependencies(std::set<int>& varset) {
        mimo->getScalarDependencies(varset);
    }

    virtual void getRotDependencies(std::set<int>& varset) {
        mimo->getRotDependencies(varset);
    }


    virtual void debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        mimo->debug_printtree();
        std::cout << ")";
    }
    // typically for this type of objects, derivativeExpression is not allowed,
    // (but you can always override it)
    virtual typename Expression<ResultType>::Ptr derivativeExpression(int i) {
        assert( 0 /* derivativeExpression of MIMO_Output not allowed */);
        return typename Expression<ResultType>::Ptr();
    }

    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        mimo->print(os);
        os << ")";
    }

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
           << COLOR_OPERATION << ",color=black]\n";
        pnumber argnode1;
        mimo->write_dotfile_update(of,argnode1);
        of << "S"<<thisnode<<" -> " << "S"<<argnode1<< "\n"; //rev
    }
    virtual void write_dotfile_init() {
        mimo->write_dotfile_init();
    } 
    virtual MIMO::Ptr getMIMOClone() {
        nr_of_clones++;
        return mimo->getClone(nr_of_clones);
    }
};

} // namespace KDL
#endif
