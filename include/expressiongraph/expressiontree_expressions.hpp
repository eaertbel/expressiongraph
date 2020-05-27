/**
 * @file expressiontree_expressions.hpp
 * @brief contains basic infrastructure for expression trees.
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

#ifndef KDL_EXPRESSIONTREE_EXPRESSIONS_HPP
#define KDL_EXPRESSIONTREE_EXPRESSIONS_HPP

#include <kdl/frames.hpp>
#include <kdl/frames_io.hpp>
//#include <kdl/framevel.hpp>
#include <kdl/stiffness.hpp>
#include "quat.hpp"
#include "quat_io.hpp"
#include <boost/smart_ptr.hpp>
#include <kdl/utilities/utility.h>
#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <ostream>
#include <set>
#include <list>
#include <cmath>
#include "expressiontree_traits.hpp"
#include "expressiontree_exceptions.hpp"
// colorscheme:
#define COLOR_OPERATION "\"#5CCCCC\""
#define COLOR_CACHED    "\"#FF7400\""
#define COLOR_LEAF      "\"#FFB273\""

#include <stdint.h>
typedef uintptr_t  pnumber;

namespace KDL {


class ExpressionOptimizer;


/**
 * start a dotfile 
 *
 * Use it as follows:
 *  - write_dotfile_start(of)
 *  - e1->write_dotfile_init()
 *  - e2->write_dotfile_init()
 *  - e1->write_dotfile_update(of)
 *  - e2->write_dotfile_update(of)
 *  - write_dotfile_end(of)
 *
 * You can use e->write_dotfile(of) if you only have to write one expression to the dot-file.
 *
 * \param [in] of stream to output to. 
 * 
 */
inline void write_dotfile_start(std::ostream& of) {
     of << "digraph expressiontree { \n"
       << "rankdir=BT\n\n"; // rev
}

/**
 * end a dotfile
 * \param [in] of stream to output to. 
 */
inline void write_dotfile_end(std::ostream& of) {
    of << "}";
}

template <typename T>
class Expression;


/**
 * Definition of all methods of Expression<T> whose interface does not depend on T.
 */
class ExpressionBase {
public:
   typedef boost::shared_ptr< ExpressionBase > Ptr;


   /**
     * Fills in the input values for this expression. 
     * This method call is passed through to all underlying nodes of the expression tree.
     * A Node decides for itself whether it can use any of the input values.  The default
     * action it to do nothing and pass it to the descendants.
     */
    virtual void setInputValues(const std::vector<double>& values)=0;

    /**
     * Fills in the input values for this expression. 
     * This method call is passed through to all underlying nodes of the expression tree.
     * A Node decides for itself whether it can use any of the input values.  The default
     * action it to do nothing and pass it to the descendants.
     *
     * The effect of calling this method is the same as repeatedly calling 
     *    setInputValue(ndx[i], values[i])
     * for all i = 0.. size-1
     *
     * @param ndx a vector of indices, 
     * @param values a vector of values, should be the same size as ndx.
     */
    virtual void setInputValues(const std::vector<int>& ndx,const std::vector<double>& values) {
        EG_ASSERT(ndx.size()==values.size());
        for (size_t i=0;i<ndx.size();++i) {
            setInputValue(ndx[i],values[i]);
        }
    }
 
   /**
     * Fills in the input values for this expression. 
     * This method call is passed through to all underlying nodes of the expression tree.
     * A Node decides for itself whether it can use any of the input values.  The default
     * action it to do nothing and pass it to the descendants.
     *
     * The effect of calling this method is the same as repeatedly calling 
     *    setInputValue(ndx[i], values[i])
     * for all i = 0.. size-1
     *
     * @param ndx a vector of indices, 
     * @param values a vector of values, should be the same size as ndx.
     */
    virtual void setInputValues(const std::vector<int>& ndx,const Eigen::VectorXd& values) {
        EG_ASSERT(ndx.size()==(size_t)values.rows());
        for (size_t i=0;i<ndx.size();++i) {
            setInputValue(ndx[i],values[i]);
        }
    }

   /**
     * Identical to the other setInputValues(ndx,values) method, but with
     * an std::vector<Rotation> for the values.
     * @param ndx a vector of indices, of course, these should be at least 3 indices apart 
     * @param values a vector of values, should be the same size as ndx.
     */
    virtual void setInputValues(const std::vector<int>& ndx,const std::vector<Rotation>& values) {
        EG_ASSERT(ndx.size()==values.size());
        for (size_t i=0;i<ndx.size();++i) {
            setInputValue(ndx[i],values[i]);
        }
    }

    /**
     * Fills in the input values for this expression. 
     * This method call is passed through to all underlying nodes of the expression tree.
     * A Node decides for itself whether it can use any of the input values.  The default
     * action it to do nothing and pass it to the descendants.
     * This routines allows to set the variable_number-th variable with the geven value "val"
     */
    virtual void setInputValue(int variable_number, double val) = 0;

    /**
     * Fills in the input values for this expression.  This method assigns a Rotation variable.
     * This method call is passed through to all underlying nodes of the expression tree.
     * A Node decides for itself whether it can use any of the input values.  The default
     * action it to do nothing and pass it to the descendants.
     * This routines allows to set the variable_number-th variable with the geven value "val"
     */
    virtual void setInputValue(int variable_number, const Rotation& val) = 0;


    /**
     * Fills in the input values for this expression. 
     * This method call is passed through to all underlying nodes of the expression tree.
     * A Node decides for itself whether it can use any of the input values.  The default
     * action it to do nothing and pass it to the descendants.
     * This routines allows to set the variable_number-th variable with the geven value "val"
     */
    virtual void setInputValue(double val) {
        setInputValue(0,val);
    } 

   
    /**
     * add the expression to the given ExpressionOptimizer.  ExpressionOptimizer can this
     * later use to optimize setInputValue() and cache-update strategies.
     */
    virtual void addToOptimizer(ExpressionOptimizer& opt) = 0;


    /**
     * determines wether an expression is constant, e.g. can be used to determine 
     * whether an expression can be folded by evaluating this expression and store a
     * cached representation.
     * This is not exactly the same then having an empty dependency set.
     * E.g. for when using make_constant(), or using variable with no Jacobian specified.
     * Also used in the implementation of expressiongraph functions (where we have to be
     * careful with our dependencies and folding of expressions).
     *
     * SEMANTICS: expression is certified to be constant.
     * If it is false, we are not sure whether it is constant.
     */
    virtual bool isConstant() const = 0;
    /**
     * INVARIANT CONDITION:
     *    The derivative towards an index that is not in the varset returned by this function
     *    is garanteed to be zero.
     *
     * adds the dependencies on (all types) input variables numbers to the given set of variables.
     * For rotational input variables, three dependencies are returned.
     * For VariableType nodes, the during construction specified variable numbers are returned.
     * \warning getDependencies is not the union of getScalarDependencies and getRotDependencies, since
     *          also VariableType nodes influence getDependencies, and these are of neither InputType or InputRotationType
     */
    virtual void getDependencies(std::set<int>& varset) = 0;

    /**
     * adds the dependencies on scalar input variables numbers to the given set of variables.
     * rotational input variables are ignored
     */
    virtual void getScalarDependencies(std::set<int>& varset) = 0;

    /**
     * adds the dependencies on rotational input variables numbers to the given set of variables.
     * Scalar input variables are ignored.
     */
    virtual void getRotDependencies(std::set<int>& varset) = 0;

    /**
     * For all VariableType objects in the expression, update their value and Jacobian
     * from the original VariableType (when they  cloned).
     */
    virtual void update_variabletype_from_original() = 0;
    /**
     * get a named subexpression.
     * Cached nodes can be given a name.  A reference to these nodes
     * can be obtained using this name.  There are no checks for duplicate names, the behavior
     * is undefined with duplicate names.
     * Returns a null pointer if no subexpression with the name is found. 
     * Returns also a null pointer if the type of the cached node does not correspond to the
     * type asked.
     * 
     * example usage:  a = subExpression_Frame("blablabla"); if (a) { action_when_subexpr_exists }
     */
    virtual boost::shared_ptr<Expression<Frame> > subExpression_Frame(const std::string& name) = 0; 
    virtual boost::shared_ptr<Expression<Rotation> > subExpression_Rotation(const std::string& name) = 0; 
    virtual boost::shared_ptr<Expression<Vector> > subExpression_Vector(const std::string& name) = 0; 
    virtual boost::shared_ptr<Expression<Twist> > subExpression_Twist(const std::string& name) = 0; 
    virtual boost::shared_ptr<Expression<Wrench> > subExpression_Wrench(const std::string& name) = 0; 
    virtual boost::shared_ptr<Expression<double> > subExpression_Double(const std::string& name) = 0; 
 
    virtual void debug_printtree()=0;

    virtual void print(std::ostream& os) const = 0;

    virtual void write_dotfile(std::ostream& of) = 0;

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) = 0;

    virtual void write_dotfile_init() = 0;

    virtual int number_of_derivatives() = 0;

    /**
     * resize datastructures (CachedType) when the number of derivatives has changed.
     * (Only increases the number of derivatives)
     * Necessary because of binding of the function arguments when calling an
     * expressiongraph function can cause the number of derivatives to change.
     */
    virtual void resize_nr_of_derivatives() = 0;
    /**
     * if the expression is a scalar variable (i.e. InputType object) this
     * method will return the variable index number.
     * Otherwise it returns -1.
     */
    virtual int isScalarVariable() = 0;
    
    virtual ExpressionType getResultType() = 0;

    virtual boost::shared_ptr< ExpressionBase> cloneBase() = 0; 

    virtual ~ExpressionBase() {}
};



/** The following classes define the structure of the expression tree
 * an expression can have a type (e.g. double, KDL::Vector, etc....)
 * it can encapsulate a unary function ( e.g. sin(x) )
 * or a binary function ( e.g. dot(x,y) for the dotproduct between vectors)
 *
 * Two of the special types are :
 *   ConstantType encapsulates a constant value, i.e. a value whose derivatives are zero.
 *   InputType
 *
 * An expression can have a (values of) derivatives towards an arbitrary number of variables.
 *   When the derivative is not given, it is zero.
 *
 *   IMPORTANT:  value() always has to be called before derivative(int i) on an expression.
 *     (see class definition).  
 *     A Typical pattern to call the methods is
 *      - expr->setInputValues(q) 
 *      - a = expr->value()
 *      - ad1 = expr->derivative(0)
 *      - ad2 = expr->derivative(1)
 *
 *
 * There are two special types of expression graph nodes that respond to setInputValue(...):
 *  - InputType
 *  - InputRotationType
 * All other nodes just pass calls to setInputValue() and setInputValues() along or ignore them.
 *
 * Example. Input variables can be scalar values or variables of the type Rotation.
 * A rotation variable is equivalent to 3 scalar values. e.g. you can
 * have the following variables at the following index values:
 *   - 0 : Rotation
 *   - 3 : Scalar
 *   - 4 : Scalar
 * In this way, you can request the derivatives for indices 0 up to 4.  
 * Index 0 corresponds to the x component of the angular velicty.
 * Index 1 corresponds to the y component of the angular velicty.
 * Index 2 corresponds to the z component of the angular velicty.
 * For the positions, however, you specify a Rotation matrix.
 *
 * You have to use the correct type (Scalar/Rotation).  If setInputValue(ndx,scalar) is
 * called on an InputRotation object, it is ignored. (and vice-versa).
 */

template< typename ResultType >
class Expression: public ExpressionBase {
public:
    std::string name;

    typedef typename boost::shared_ptr< Expression<ResultType> > Ptr;
    typedef typename AutoDiffTrait<ResultType>::DerivType DerivType;

    Expression() {}

    Expression(const std::string& _name) : name(_name) {};
    /**
     * returns the value of the expression tree.
     */
    virtual ResultType value()      = 0;


    /**
     * returns the derivative of the expression tree towards variable
     * i. It can be assumed that value() is always called before
     * derivative().  For efficiency reasons the call to value()
     * should cache results needed for the call to derivative(i).
     */
    virtual DerivType derivative(int i) = 0;


    /**
     * returns towards how many variables the derivative is computed.
     * This function is only introduced for efficiency reasons
     */
    virtual int number_of_derivatives() = 0;

    /**
     * resizes the data structure of some of the nodes to accomodate a change in
     * the number of derivatives.  Necessary for expressiongraph functions where
     * the binding of arguments can change the number of derivatives.  This changes
     * the CachedType nodes and possibly others.
     */
    virtual void resize_nr_of_derivatives() = 0;

    virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
        throw NotImplementedException();
    }

   
    /**
     * makes a deep copy of the expression.
     */
    virtual boost::shared_ptr<Expression<ResultType> > clone() = 0;
    virtual boost::shared_ptr< ExpressionBase> cloneBase() {
        return clone();
    } 


    virtual void debug_printtree() {
        std::cout << name;
    }
    virtual void print(std::ostream& os) const {
        os << name;
    } 

    /**
     * \param [in] of stream to output to. 
     */
    virtual void write_dotfile(std::ostream& of) {
        of << "digraph expressiontree { \n"
           << "rankdir=BT\n\n"; // rev
        pnumber argnode;
        write_dotfile_init();
        write_dotfile_update(of,argnode);  
        of << "}";
    }

    /** 
     * write this node to the .dot file (has to be initialized)
     * \param [in] thisnode: returns the value of this node,
     */
    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << name << "\",shape=box,style=filled,fillcolor=" << COLOR_OPERATION << ",color=black]\n";
    }

    /**
     * write this node to the .dot file (has to be initialized)
     */
    virtual void write_dotfile_update(std::ostream& of) {
        pnumber argnode;
        write_dotfile_update(of,argnode);
    }
 
    /**
     * hook to be used if a node wants to be initialized before write_dotfile is called.
     * (can be overridden by nodes inheriting from this base class).
     */
    virtual void write_dotfile_init() {
    }

    virtual int isScalarVariable() {
        return -1;
    }

    virtual ExpressionType getResultType() {
        return AutoDiffTrait<ResultType>::expr_type;
    }

    virtual ~Expression() {
    }
};

template<typename T>
typename Expression<T>::Ptr checkConstant( typename Expression<T>::Ptr a );
 

template< typename ResultType, typename T>
class UnaryExpression: public Expression<ResultType> {
public:
    typedef T ArgumentType;
    typedef Expression<ArgumentType> ArgumentExpr;
    typename ArgumentExpr::Ptr argument;

    UnaryExpression() {}

    UnaryExpression( const std::string& name,const typename ArgumentExpr::Ptr& ptr):
            Expression<ResultType>(name)
    {
        argument = checkConstant<T>(ptr);
    }

    virtual void setInputValues(const std::vector<double>& values) {
        argument->setInputValues(values);
    }
    virtual void setInputValue(int variable_number, double val) {
        argument->setInputValue(variable_number,val);
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        argument->setInputValue(variable_number,val);
    }

    virtual int number_of_derivatives() {
        return argument->number_of_derivatives();
    } 
    virtual void resize_nr_of_derivatives() {
        argument->resize_nr_of_derivatives();
    }

    virtual Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        return argument->subExpression_Frame(name);
    }
    virtual Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        return argument->subExpression_Rotation(name);
    }
    virtual Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        return argument->subExpression_Vector(name);
    }
    virtual Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        return argument->subExpression_Twist(name);
    }
    virtual Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        return argument->subExpression_Wrench(name);
    }
    virtual Expression<double>::Ptr subExpression_Double(const std::string& name) {
        return argument->subExpression_Double(name);
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        argument->addToOptimizer(opt);
    }
    virtual bool isConstant() const {
        return argument->isConstant();
    }
    virtual void getDependencies(std::set<int>& varset) {
        argument->getDependencies(varset);
    }
    virtual void getScalarDependencies(std::set<int>& varset) {
        argument->getScalarDependencies(varset);
    }
    virtual void getRotDependencies(std::set<int>& varset) {
        argument->getRotDependencies(varset);
    }

    virtual void update_variabletype_from_original() {
        argument->update_variabletype_from_original();
    }
 
    virtual void debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        argument->debug_printtree();
        std::cout << ")";
    }

    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        argument->print(os);
        os << ")";
    }
    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
           << COLOR_OPERATION << ",color=black]\n";
        pnumber argnode;
        argument->write_dotfile_update(of,argnode);
        of << "S"<<thisnode<<" -> " << "S"<<argnode<< "\n"; // rev.
    }
    virtual void write_dotfile_init() {
        argument->write_dotfile_init();
    } 
};


template< typename ResultType, typename T1,typename  T2>
class BinaryExpression: public Expression<ResultType> {
public:
    typedef Expression<T1>                      Argument1Expr;
    typedef Expression<T2>                      Argument2Expr;

    typename Argument1Expr::Ptr argument1;
    typename Argument2Expr::Ptr argument2;

    BinaryExpression() {}

    BinaryExpression( const std::string& name,
                      const typename Argument1Expr::Ptr& arg1ptr,
                      const typename Argument2Expr::Ptr& arg2ptr):
            Expression<ResultType>(name),
            argument1(checkConstant<T1>(arg1ptr)),
            argument2(checkConstant<T2>(arg2ptr)) {
    }

    virtual void setInputValues(const std::vector<double>& values) {
        argument1->setInputValues(values);
        argument2->setInputValues(values);
    }

    virtual void setInputValue(int variable_number, double val) {
        argument1->setInputValue(variable_number,val);
        argument2->setInputValue(variable_number,val);
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        argument1->setInputValue(variable_number,val);
        argument2->setInputValue(variable_number,val);
    }

    virtual int number_of_derivatives() {
        int n1 = argument1->number_of_derivatives();
        int n2 = argument2->number_of_derivatives();
        return n1 > n2 ? n1 : n2;
    } 
    virtual void resize_nr_of_derivatives() {
        argument1->resize_nr_of_derivatives();
        argument2->resize_nr_of_derivatives();
    }


    virtual Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        typename Expression<Frame>::Ptr a;
        a = argument1->subExpression_Frame(name);
        if (a) {
            return a;
        }
        return argument2->subExpression_Frame(name);
    }
    virtual Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        typename Expression<Rotation>::Ptr a;
        a = argument1->subExpression_Rotation(name);
        if (a) {
            return a;
        }
        return argument2->subExpression_Rotation(name);
    }
    virtual Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        typename Expression<Vector>::Ptr a;
        a = argument1->subExpression_Vector(name);
        if (a) {
            return a;
        }
        return argument2->subExpression_Vector(name);
    }
    virtual  Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        typename Expression<Twist>::Ptr a;
        a = argument1->subExpression_Twist(name);
        if (a) {
            return a;
        }
        return argument2->subExpression_Twist(name);
    }
    virtual Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        typename Expression<Wrench>::Ptr a;
        a = argument1->subExpression_Wrench(name);
        if (a) {
            return a;
        }
        return argument2->subExpression_Wrench(name);
    }
    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        typename Expression<double>::Ptr a;
        a = argument1->subExpression_Double(name);
        if (a) {
            return a;
        }
        return argument2->subExpression_Double(name);
    }


    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        argument1->addToOptimizer(opt);
        argument2->addToOptimizer(opt);
    }
    virtual bool isConstant() const {
        return (argument1->isConstant()) && (argument2->isConstant());
    }
 
    virtual void getDependencies(std::set<int>& varset) {
        argument1->getDependencies(varset);
        argument2->getDependencies(varset);
    }

    virtual void getScalarDependencies(std::set<int>& varset) {
        argument1->getScalarDependencies(varset);
        argument2->getScalarDependencies(varset);
    }

    virtual void getRotDependencies(std::set<int>& varset) {
        argument1->getRotDependencies(varset);
        argument2->getRotDependencies(varset);
    }

    virtual void update_variabletype_from_original() {
        argument1->update_variabletype_from_original();
        argument2->update_variabletype_from_original();
    }


    virtual void debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        argument1->debug_printtree();
        std::cout << ",";
        argument2->debug_printtree();
        std::cout << ")";
    }

    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        argument1->print(os);
        os << ",";
        argument2->print(os);
        os << ")";
    }

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
           << COLOR_OPERATION << ",color=black]\n";
        pnumber argnode1,argnode2;
        argument1->write_dotfile_update(of,argnode1);
        argument2->write_dotfile_update(of,argnode2);
        of << "S"<<thisnode<<" -> " << "S"<<argnode1<< "\n"; //rev
        of << "S"<<thisnode<<" -> " << "S"<<argnode2<< "\n"; //rev
    }
    virtual void write_dotfile_init() {
        argument1->write_dotfile_init();
        argument2->write_dotfile_init();
    } 
};

template< typename ResultType, typename T1,typename  T2,typename T3>
class TernaryExpression: public Expression<ResultType> {
public:
    typedef Expression<T1>                      Argument1Expr;
    typedef Expression<T2>                      Argument2Expr;
    typedef Expression<T3>                      Argument3Expr;

    typename Argument1Expr::Ptr argument1;
    typename Argument2Expr::Ptr argument2;
    typename Argument3Expr::Ptr argument3;



    TernaryExpression() {}

    TernaryExpression( const std::string& name,
                      const typename Argument1Expr::Ptr& arg1ptr,
                      const typename Argument2Expr::Ptr& arg2ptr,
                      const typename Argument3Expr::Ptr& arg3ptr
                      ):
        Expression<ResultType>(name),
        argument1(checkConstant<T1>(arg1ptr)),
        argument2(checkConstant<T2>(arg2ptr)),
        argument3(checkConstant<T3>(arg3ptr)) {
    }

    virtual void setInputValues(const std::vector<double>& values) {
        argument1->setInputValues(values);
        argument2->setInputValues(values);
        argument3->setInputValues(values);
    }

    virtual void setInputValue(int variable_number, double val) {
        argument1->setInputValue(variable_number,val);
        argument2->setInputValue(variable_number,val);
        argument3->setInputValue(variable_number,val);
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        argument1->setInputValue(variable_number,val);
        argument2->setInputValue(variable_number,val);
        argument3->setInputValue(variable_number,val);
    }

    virtual int number_of_derivatives() {
        int n1 = argument1->number_of_derivatives();
        int n2 = argument2->number_of_derivatives();
        int n3 = argument3->number_of_derivatives();
        return std::max(n1,std::max(n2,n3));
 
    } 
    virtual void resize_nr_of_derivatives() {
        argument1->resize_nr_of_derivatives();
        argument2->resize_nr_of_derivatives();
        argument3->resize_nr_of_derivatives();
    }


    virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        typename Expression<Frame>::Ptr a;
        a = argument1->subExpression_Frame(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Frame(name);
        if (a) {
            return a;
        }
        return argument3->subExpression_Frame(name);
    }

    virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        typename Expression<Rotation>::Ptr a;
        a = argument1->subExpression_Rotation(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Rotation(name);
        if (a) {
            return a;
        }
        return argument3->subExpression_Rotation(name);
    }

    virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        typename Expression<Vector>::Ptr a;
        a = argument1->subExpression_Vector(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Vector(name);
        if (a) {
            return a;
        }
        return argument3->subExpression_Vector(name);
    }

    virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        typename Expression<Wrench>::Ptr a;
        a = argument1->subExpression_Wrench(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Wrench(name);
        if (a) {
            return a;
        }
        return argument3->subExpression_Wrench(name);
    }

    virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        typename Expression<Twist>::Ptr a;
        a = argument1->subExpression_Twist(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Twist(name);
        if (a) {
            return a;
        }
        return argument3->subExpression_Twist(name);
    }

    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        typename Expression<double>::Ptr a;
        a = argument1->subExpression_Double(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Double(name);
        if (a) {
            return a;
        }
        return argument3->subExpression_Double(name);
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        argument1->addToOptimizer(opt);
        argument2->addToOptimizer(opt);
        argument3->addToOptimizer(opt);
    }

    virtual bool isConstant() const {
      return
        argument1->isConstant() &&
        argument2->isConstant() &&
        argument3->isConstant();
    }

    virtual void getDependencies(std::set<int>& varset) {
        argument1->getDependencies(varset);
        argument2->getDependencies(varset);
        argument3->getDependencies(varset);
    }

    virtual void getScalarDependencies(std::set<int>& varset) {
        argument1->getScalarDependencies(varset);
        argument2->getScalarDependencies(varset);
        argument3->getScalarDependencies(varset);
    }

    virtual void getRotDependencies(std::set<int>& varset) {
        argument1->getRotDependencies(varset);
        argument2->getRotDependencies(varset);
        argument3->getRotDependencies(varset);
    }

    virtual void update_variabletype_from_original() {
        argument1->update_variabletype_from_original();
        argument2->update_variabletype_from_original();
        argument3->update_variabletype_from_original();
    }

    virtual void debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        argument1->debug_printtree();
        std::cout << ",";
        argument2->debug_printtree();
        std::cout << ",";
        argument3->debug_printtree();
        std::cout << ")";
    }
    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        argument1->print(os);
        os << ",";
        argument2->print(os);
        os << ",";
        argument3->print(os);
        os << ")";
    }

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
           << COLOR_OPERATION << ",color=black]\n";
        pnumber argnode1,argnode2,argnode3;
        argument1->write_dotfile_update(of,argnode1);
        argument2->write_dotfile_update(of,argnode2);
        argument3->write_dotfile_update(of,argnode3);
        of << "S"<<thisnode<<" -> " << "S"<<argnode1<< "\n"; //rev
        of << "S"<<thisnode<<" -> " << "S"<<argnode2<< "\n"; //rev
        of << "S"<<thisnode<<" -> " << "S"<<argnode3<< "\n"; //rev
    }
    virtual void write_dotfile_init() {
        argument1->write_dotfile_init();
        argument2->write_dotfile_init();
        argument3->write_dotfile_init();
    } 
};

template< typename ResultType, typename T1,typename  T2,typename T3, typename T4>
class QuaternaryExpression: public Expression<ResultType> {
public:
    typedef Expression<T1>                      Argument1Expr;
    typedef Expression<T2>                      Argument2Expr;
    typedef Expression<T3>                      Argument3Expr;
    typedef Expression<T4>                      Argument4Expr;

    typename Argument1Expr::Ptr argument1;
    typename Argument2Expr::Ptr argument2;
    typename Argument3Expr::Ptr argument3;
    typename Argument4Expr::Ptr argument4;



    QuaternaryExpression() {}

    QuaternaryExpression( const std::string& name,
                      const typename Argument1Expr::Ptr& arg1ptr,
                      const typename Argument2Expr::Ptr& arg2ptr,
                      const typename Argument3Expr::Ptr& arg3ptr,
                      const typename Argument4Expr::Ptr& arg4ptr
                      ):
        Expression<ResultType>(name),
        argument1(checkConstant<T1>(arg1ptr)),
        argument2(checkConstant<T2>(arg2ptr)),
        argument3(checkConstant<T3>(arg3ptr)),
        argument4(checkConstant<T4>(arg4ptr)) {
    }

    virtual void setInputValues(const std::vector<double>& values) {
        argument1->setInputValues(values);
        argument2->setInputValues(values);
        argument3->setInputValues(values);
        argument4->setInputValues(values);
    }

    virtual void setInputValue(int variable_number, double val) {
        argument1->setInputValue(variable_number,val);
        argument2->setInputValue(variable_number,val);
        argument3->setInputValue(variable_number,val);
        argument4->setInputValue(variable_number,val);
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        argument1->setInputValue(variable_number,val);
        argument2->setInputValue(variable_number,val);
        argument3->setInputValue(variable_number,val);
        argument4->setInputValue(variable_number,val);
    }

    virtual int number_of_derivatives() {
        int n1 = argument1->number_of_derivatives();
        int n2 = argument2->number_of_derivatives();
        int n3 = argument3->number_of_derivatives();
        int n4 = argument3->number_of_derivatives();
        return std::max(std::max(n1,n2),std::max(n3,n4));
 
    } 
    virtual void resize_nr_of_derivatives() {
        argument1->resize_nr_of_derivatives();
        argument2->resize_nr_of_derivatives();
        argument3->resize_nr_of_derivatives();
        argument4->resize_nr_of_derivatives();
    }


    virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        typename Expression<Frame>::Ptr a;
        a = argument1->subExpression_Frame(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Frame(name);
        if (a) {
            return a;
        }
        a = argument3->subExpression_Frame(name);
        if (a) {
            return a;
        }
        return argument4->subExpression_Frame(name);
    }

    virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        typename Expression<Rotation>::Ptr a;
        a = argument1->subExpression_Rotation(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Rotation(name);
        if (a) {
            return a;
        }
        a = argument3->subExpression_Rotation(name);
        if (a) {
            return a;
        }
        return argument4->subExpression_Rotation(name);
    }

    virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        typename Expression<Vector>::Ptr a;
        a = argument1->subExpression_Vector(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Vector(name);
        if (a) {
            return a;
        }
        a = argument3->subExpression_Vector(name);
        if (a) {
            return a;
        }
        return argument4->subExpression_Vector(name);
    }

    virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        typename Expression<Wrench>::Ptr a;
        a = argument1->subExpression_Wrench(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Wrench(name);
        if (a) {
            return a;
        }
        a = argument3->subExpression_Wrench(name);
        if (a) {
            return a;
        }
        return argument4->subExpression_Wrench(name);
    }

    virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        typename Expression<Twist>::Ptr a;
        a = argument1->subExpression_Twist(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Twist(name);
        if (a) {
            return a;
        }
        a = argument3->subExpression_Twist(name);
        if (a) {
            return a;
        }
        return argument4->subExpression_Twist(name);
    }

    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        typename Expression<double>::Ptr a;
        a = argument1->subExpression_Double(name);
        if (a) {
            return a;
        }
        a = argument2->subExpression_Double(name);
        if (a) {
            return a;
        }
        a = argument3->subExpression_Double(name);
        if (a) {
            return a;
        }
        return argument4->subExpression_Double(name);
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        argument1->addToOptimizer(opt);
        argument2->addToOptimizer(opt);
        argument3->addToOptimizer(opt);
        argument4->addToOptimizer(opt);
    }

    virtual bool isConstant() const {
      return
        argument1->isConstant() &&
        argument2->isConstant() &&
        argument3->isConstant() &&
        argument4->isConstant();
    }

    virtual void getDependencies(std::set<int>& varset) {
        argument1->getDependencies(varset);
        argument2->getDependencies(varset);
        argument3->getDependencies(varset);
        argument4->getDependencies(varset);
    }

    virtual void getScalarDependencies(std::set<int>& varset) {
        argument1->getScalarDependencies(varset);
        argument2->getScalarDependencies(varset);
        argument3->getScalarDependencies(varset);
        argument4->getScalarDependencies(varset);
    }

    virtual void getRotDependencies(std::set<int>& varset) {
        argument1->getRotDependencies(varset);
        argument2->getRotDependencies(varset);
        argument3->getRotDependencies(varset);
        argument4->getRotDependencies(varset);
    }

    virtual void update_variabletype_from_original() {
        argument1->update_variabletype_from_original();
        argument2->update_variabletype_from_original();
        argument3->update_variabletype_from_original();
        argument4->update_variabletype_from_original();
    }

    virtual void debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        argument1->debug_printtree();
        std::cout << ",";
        argument2->debug_printtree();
        std::cout << ",";
        argument3->debug_printtree();
        std::cout << ",";
        argument4->debug_printtree();
        std::cout << ")";
    }
    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        argument1->print(os);
        os << ",";
        argument2->print(os);
        os << ",";
        argument3->print(os);
        os << ",";
        argument4->print(os);
        os << ")";
    }

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
           << COLOR_OPERATION << ",color=black]\n";
        pnumber argnode1,argnode2,argnode3,argnode4;
        argument1->write_dotfile_update(of,argnode1);
        argument2->write_dotfile_update(of,argnode2);
        argument3->write_dotfile_update(of,argnode3);
        argument4->write_dotfile_update(of,argnode4);
        of << "S"<<thisnode<<" -> " << "S"<<argnode1<< "\n"; //rev
        of << "S"<<thisnode<<" -> " << "S"<<argnode2<< "\n"; //rev
        of << "S"<<thisnode<<" -> " << "S"<<argnode3<< "\n"; //rev
        of << "S"<<thisnode<<" -> " << "S"<<argnode4<< "\n"; //rev
    }
    virtual void write_dotfile_init() {
        argument1->write_dotfile_init();
        argument2->write_dotfile_init();
        argument3->write_dotfile_init();
        argument4->write_dotfile_init();
    } 
};


template <typename _ResultType>
class FunctionType: public Expression<_ResultType> {
public:
    bool dot_already_written;
    typedef _ResultType ResultType;
    typedef typename AutoDiffTrait<_ResultType>::DerivType DerivType;

    FunctionType():dot_already_written(false) {}
    FunctionType(const std::string& name):
        Expression<_ResultType>(name),dot_already_written(false) {}

    virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        return typename Expression<Frame>::Ptr();
    }

    virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        return typename Expression<Rotation>::Ptr();
    }

    virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        return typename Expression<Vector>::Ptr();
    }

    virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        return typename Expression<Twist>::Ptr();
    }

    virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        return typename Expression<Wrench>::Ptr();
    }

    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        return typename Expression<double>::Ptr();
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
    }

    virtual void getDependencies(std::set<int>& varset) {
    }
    virtual void getScalarDependencies(std::set<int>& varset) {
    }
    virtual void getRotDependencies(std::set<int>& varset) {
    }

    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name;
    }

    virtual void write_dotfile_init() {
        dot_already_written = false;
    }



    void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        if (!dot_already_written) {
            dot_already_written = true;
            thisnode=(size_t)this;
            of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
               << COLOR_LEAF << ",color=black]\n";
        } else {
            thisnode=(size_t)this;
        }
    }
};

namespace detail {
    inline void print(std::ostream& os, double val) {
        os << "constant(" << val << ")"; 
    }

    inline void print(std::ostream& os, const KDL::Vector& v) {
        os << "constant(Vector(" << v.x() << "," << v.y() << "," << v.z()  << "))"; 
    }

    inline void print(std::ostream& os, const KDL::Rotation& R) {
        os << "constant(Rotation(" << R(0,0) << "," << R(0,1) << "," << R(0,2)  << "," 
                                                     << R(1,0) << "," << R(1,1) << "," << R(1,2)  << "," 
                                                     << R(2,0) << "," << R(2,1) << "," << R(2,2)  << "))" ;
    }

    inline void print(std::ostream& os, const Quaternion&  val) {
            os << "constant( Quaternion(";
            KDL::operator<<(os, val);
            os << "))"; 
    }

    inline void print(std::ostream& os, const KDL::Frame& F) {
        os <<  "constant(Frame(Rotation(" << F.M(0,0) << "," << F.M(0,1) << "," << F.M(0,2)  << "," 
                                            << F.M(1,0) << "," << F.M(1,1) << "," << F.M(1,2)  << "," 
                                            << F.M(2,0) << "," << F.M(2,1) << "," << F.M(2,2)  << "),";
        os << "Vector("<<F.p(0) << "," << F.p(1) << "," << F.p(2) << ")))";
    }

    inline void print(std::ostream& os, const KDL::Twist& t) {
        os <<  "constant(Twist(Vector(" << t.vel(0) << "," << t.vel(1) << "," << t.vel(2)  << ")," 
                                   << "Vector(" << t.rot(0) << "," << t.rot(1) << "," << t.rot(2)  << ")))";
    }

    inline void print(std::ostream& os, const KDL::Wrench& t) {
        os <<  "constant(Wrench(Vector(" << t.force(0) << "," << t.force(1) << "," << t.force(2)  << ")," 
                                   << "Vector(" << t.torque(0) << "," << t.torque(1) << "," << t.torque(2)  << ")))";
    }

    template <int n, int m>
    inline void print(std::ostream& os, const Eigen::Matrix<double,n,m>&  t) {
        os <<  "constant(Matrix["<<n<<","<<m<<"](" << t << "))";
    }

}// namespace detail

template<typename ResultType>
inline typename Expression<ResultType>::Ptr Constant( const ResultType& _val );
 
template <typename ResultType>
class ConstantType: public FunctionType<ResultType> {
public:
    ResultType val;
    typedef typename AutoDiffTrait<ResultType>::DerivType DerivType;
    ConstantType() {}
    ConstantType(const ResultType& _val):
        FunctionType<ResultType>("constant"),val(_val) {
    }

    virtual void setInputValues(const std::vector<double>& values) {
    }

    virtual void setInputValue(int variable_number, double val) {
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
    }

    virtual ResultType value() {
        return val;
    }

    virtual DerivType derivative(int i) {
        return AutoDiffTrait<ResultType>::zeroDerivative();
    }

    virtual int number_of_derivatives() {
        return 0;
    }   
    virtual void resize_nr_of_derivatives() {
    }

    virtual bool isConstant() const {
        return true;
    }

    virtual void print(std::ostream& os) const {
        detail::print(os,val);
    }
    virtual void update_variabletype_from_original() {}

    virtual typename Expression<ResultType>::Ptr clone() {
        typename Expression<ResultType>::Ptr expr(
            new ConstantType( val )
        );
        return expr;
    }
    virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
        return Constant( AutoDiffTrait<DerivType>::zeroDerivative() );
    }

};

/** utility function to create ConstantType */
template<typename ResultType>
inline typename Expression<ResultType>::Ptr Constant( const ResultType& _val ) {
   typename Expression<ResultType>::Ptr cnst(
        new ConstantType<ResultType>( _val )
   );
   return cnst;
}



class InputType : public FunctionType<double> {
public:
    char name_buffer[32];
    int    variable_number;
    double val;
    InputType() {}
    /**
     * defaultvalue specifies the initial value of the variable "value".  The value
     * when it is not filled in using setInputValues() method ( because it is not called, or
     * because variable_number does not fall in the range of values specified by setInputValues()
     *
     */
    InputType(int _variable_number, double _defaultvalue):
        FunctionType<double>("input"),
        variable_number(_variable_number),
        val(_defaultvalue) {
            EG_ASSERT_MSG( variable_number >= 0,"Variable number for InputType should be >=0");
            sprintf(name_buffer,"input(%d)",variable_number);
            name = name_buffer;
    }

    virtual void setInputValues(const std::vector<double>& values) {
        if (variable_number < (int)values.size()) {
            val = values[variable_number];
        }
    }

    virtual void setInputValue(int _variable_number, double _val) {
        if (variable_number == _variable_number) {
            val = _val;
        }
    }
 
    virtual void setInputValue(int variable_number, const Rotation& val) {
    }

    virtual double value() {
        return val;
    } 

    virtual void addToOptimizer(ExpressionOptimizer& opt);
    virtual bool isConstant() const {
        return false;
    }
    virtual void getDependencies(std::set<int>& varset) {
        varset.insert( variable_number);
    }
    virtual void getScalarDependencies(std::set<int>& varset) {
        varset.insert( variable_number);
    }
    virtual void getRotDependencies(std::set<int>& varset) {
    }

    virtual void update_variabletype_from_original() {
    }

    virtual double derivative(int i) {
        if (variable_number == i) {
            return 1.0;
        } else {
            return 0.0;
        }
    } 
 
    virtual Expression<DerivType>::Ptr derivativeExpression(int i) {
        if (variable_number== i) {
            return Constant(1.0);
        } else {
            return Constant(0.0);
        }
    }

    virtual int number_of_derivatives() {
        return variable_number+1;
    };
    virtual void resize_nr_of_derivatives() {}



    virtual int isScalarVariable() {
        return variable_number;
    }
    /**
     * \warn  Default value for the cloned object will be the value of the original InputType object.
     */
    virtual Expression<ResultType>::Ptr clone() {
        Expression<ResultType>::Ptr expr(
            new InputType( variable_number, val)
        );
        return expr;
    }

}; 

/** utility functions to create InputType 
 * creates a variable for which the value returns inputvariable with variablenumber
 * and whose derivative with number derivative_number is equal to 1.
 */
inline Expression<double>::Ptr input(int variable_number, double default_value  ) {
   Expression<double>::Ptr var(
        new InputType( variable_number,default_value )
   );
   return var;
}

inline Expression<double>::Ptr input(int variable_number ) {
   Expression<double>::Ptr var(
        new InputType( variable_number, 0.0 )
   );
   return var;
}


class InputRotationType : public FunctionType<Rotation> {
public:
    char name_buffer[32];
    int    variable_number;
    Rotation  val;
    InputRotationType() {}
    /**
     * defaultvalue specifies the initial value of the variable "value".  The value
     * when it is not filled in using setInputValues() method ( because it is not called, or
     * because variable_number does not fall in the range of values specified by setInputValues()
     *
     */
    InputRotationType(int _variable_number, const Rotation& _defaultvalue):
        FunctionType<Rotation>("input"),
        variable_number(_variable_number),
        val(_defaultvalue) {
            EG_ASSERT_MSG( variable_number >= 0,"Variable number for input type should be >=0");
            sprintf(name_buffer,"input(%d)",variable_number);
            name = name_buffer;
    }

    virtual void setInputValues(const std::vector<double>& values) {
    }

    virtual void setInputValue(int _variable_number, double _val) {
    }
 
    virtual void setInputValue(int _variable_number, const Rotation& _val) {
        if (variable_number == _variable_number) {
            val = _val;
        }
    }

    virtual Rotation value() {
        return val;
    } 

    virtual void addToOptimizer(ExpressionOptimizer& opt);
    
    virtual bool isConstant() const {
        return false;
    }

    virtual void getDependencies(std::set<int>& varset) {
        varset.insert( variable_number);
        varset.insert( variable_number+1);
        varset.insert( variable_number+2);
    }

    virtual void getScalarDependencies(std::set<int>& varset) {
    }
    virtual void getRotDependencies(std::set<int>& varset) {
        varset.insert( variable_number);
    }

    virtual void update_variabletype_from_original() {
    }



    virtual Vector derivative(int i) {
        if (variable_number == i) {
            return Vector(1,0,0);
        } if (variable_number+1 == i) {
            return Vector(0,1,0);
        } if (variable_number+2 == i) {
            return Vector(0,0,1);
        } else {
            return Vector(0,0,0);
        }
    } 
 
    virtual Expression<DerivType>::Ptr derivativeExpression(int i) {
        if (variable_number == i) {
            return Constant(Vector(1,0,0));
        } if (variable_number+1 == i) {
            return Constant(Vector(0,1,0));
        } if (variable_number+2 == i) {
            return Constant(Vector(0,0,1));
        } else {
            return Constant(Vector(0,0,0));
        }

    }

    virtual int number_of_derivatives() {
        return variable_number+3;
    };
    virtual void resize_nr_of_derivatives() {}
 
    /**
     * \warn  Default value for the cloned object will be the value of the original InputType object.
     */
    virtual Expression<ResultType>::Ptr clone() {
        Expression<ResultType>::Ptr expr(
            new InputRotationType( variable_number, val)
        );
        return expr;
    }

}; 

/**
 * creates a variable for which the value returns a Rotation, with the given variable number
 * ( and variable number + 1 and variable number + 2).  The derivative is a rotational velocity (type Vector).
 * The derivative towards variable number corresponds to Vector(1,0,0) 
 * The derivative towards variable number+1 corresponds to Vector(0,1,0) 
 * The derivative towards variable number+2 corresponds to Vector(0,0,1) 
 * \param [in] variable_number the variable number to start with
 * \param [in] default_value   a default value for the position value (i.e. a Rotation Matrix).
 * \return an expression graph of type rotation. 
 */

inline Expression<Rotation>::Ptr inputRot(int variable_number, const Rotation& default_value  ) {
   Expression<Rotation>::Ptr var(
        new InputRotationType( variable_number,default_value )
   );
   return var;
}
/**
 * creates a variable for which the value returns a Rotation, with the given variable number
 * ( and variable number + 1 and variable number + 2).  The derivative is a rotational velocity (type Vector).
 * The derivative towards variable number corresponds to Vector(1,0,0) 
 * The derivative towards variable number+1 corresponds to Vector(0,1,0) 
 * The derivative towards variable number+2 corresponds to Vector(0,0,1) 
 * \param [in] variable_number the variable number to start with
 * \param [in] default_value   a default value for the position value (i.e. a Rotation Matrix).
 * \return an expression graph of type rotation. 
 */
inline Expression<Rotation>::Ptr inputRot(int variable_number ) {
   Expression<Rotation>::Ptr var(
        new InputRotationType( variable_number, Rotation::Identity() )
   );
   return var;
}

class CachedExpression {
    public:
        virtual void getDependencies(std::set<int>& varset)=0;
        virtual void invalidate_cache() = 0;
        virtual void addToOptimizer(ExpressionOptimizer& opt);
};


/**
 * caches the first max_number_of_var results for derivative(i) and value()
 * such that needless computation is avoided when expressions are reused.
 * derivatives of variables with variable number > max_number_of_var are not cached.
 */
template <typename ResultType>
class CachedType: public Expression<ResultType>, public CachedExpression {
public:
    typedef typename AutoDiffTrait<ResultType>::DerivType DerivType;
    typename Expression<ResultType>::Ptr argument;
    ResultType val;
    std::vector<DerivType> deriv;
    bool dot_already_written;
    std::vector<bool> cached_deriv;
    bool cached_value;
    std::string cached_name;

    CachedType() {}
    /**
     * caches the first the results for derivative(i) and value()
     * (to avoid unnecessary computations)
     */
    CachedType(typename Expression<ResultType>::Ptr _argument, const std::string& _name):
        Expression<ResultType>("cached"),
        argument(checkConstant<ResultType>(_argument)),
        deriv(_argument->number_of_derivatives()), 
        cached_deriv(_argument->number_of_derivatives()),
        cached_value(false),
        cached_name(_name) {
        //#ifdef EG_LOG_CACHE
        //    std::cerr << "Constructed CachedType("<<this->cached_name<<")"<<std::endl;
        //#endif
    }

    virtual ResultType value() {
       if (cached_value) {
            #ifdef EG_CHECK_CACHE
                EG_ASSERT_MSG( val == argument->value(),"CHECK ON CachedType FAILED" );
            #endif
            #ifdef EG_LOG_CACHE
                std::cerr<<"CachedType " << this->cached_name << " : cached value used="<< val <<std::endl;
            #endif
       
            return val;
        } else {
            val          = argument->value();
            #ifdef EG_LOG_CACHE
                std::cerr<<"CachedType " << this->cached_name << " : val updated with argument->value()="<< val <<std::endl;
            #endif
            cached_value = true;
            return val;
        }
    }

    virtual void invalidate_cache() {
        #ifdef EG_LOG_CACHE
            std::cerr<<"CachedType " << this->cached_name << " invalidate_cache() has been called"<<std::endl;
        #endif
        //std::cout << "invalidate cache of " << cached_name << std::endl;
        cached_value = false;
        fill_n(cached_deriv.begin(), deriv.size(), false);
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        CachedExpression::addToOptimizer(opt);
        argument->addToOptimizer(opt);
    }
    virtual bool isConstant() const {
        return argument->isConstant();
    }
    virtual void getDependencies(std::set<int>& varset) {
        argument->getDependencies(varset);
    }
    virtual void getScalarDependencies(std::set<int>& varset) {
        argument->getScalarDependencies(varset);
    }
    virtual void getRotDependencies(std::set<int>& varset) {
        argument->getRotDependencies(varset);
    }

    virtual void update_variabletype_from_original() {
        invalidate_cache();
        argument->update_variabletype_from_original();
    }


    virtual DerivType derivative(int i) {
        EG_ASSERT(i>=0);
        if (i < (int)deriv.size() ) {
           if (cached_deriv[i]) {
                #ifdef EG_CHECK_CACHE
                    EG_ASSERT_MSG( deriv[i] == argument->derivative(i),"CHECK ON CachedType FAILED");
                #endif
                #ifdef EG_LOG_CACHE
                    std::cerr<<"CachedType " << this->cached_name << " : derivative("<<i<<") call with cached_value="<<cached_value
                            <<" and deriv.size="<<(int)deriv.size()<<" and val="<<val << "(CACHED)"<<std::endl;
                #endif
                return deriv[i];
            } else {
                deriv[i] = argument->derivative(i);
                cached_deriv[i] = true; 
                #ifdef EG_LOG_CACHE
                    std::cerr<<"CachedType " << this->cached_name << " : derivative("<<i<<") call with cached_value="<<cached_value
                            <<" and deriv.size="<<(int)deriv.size()<<" and val="<<val << "(RECOMPUTED)"<<std::endl;
                #endif
                return deriv[i];
            }
        } else {
            #ifdef EG_CHECK_CACHE
                EG_ASSERT_MSG( AutoDiffTrait<ResultType>::zeroDerivative() == argument->derivative(i), "CHECK ON CachedType FAILED" );
            #endif
            #ifdef EG_LOG_CACHE
                std::cerr<<"CachedType " << this->cached_name << " : derivative("<<i<<") returned zero because > deriv.size="<<(int)deriv.size()<< std::endl;
            #endif
            return AutoDiffTrait<ResultType>::zeroDerivative();
        }
    }

    virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
        // or should it be cached(...)
        if (this->name.size()==0) {
            typename Expression<DerivType>::Ptr retval( new CachedType<DerivType>(argument->derivativeExpression(i),""));
            return retval;
        } else {
            typename Expression<DerivType>::Ptr retval( new CachedType<DerivType>(argument->derivativeExpression(i),std::string(this->name)+"(deriv)" ));
            return retval;
        }
    }

    virtual int number_of_derivatives() {
        return deriv.size();
    }

    virtual void resize_nr_of_derivatives() {
       argument->resize_nr_of_derivatives();
        int new_size = std::max( (int)deriv.size(), argument->number_of_derivatives() );
        cached_deriv.resize(new_size);
        deriv.resize(new_size);
        #ifdef EG_LOG_CACHE
            std::cerr<<"CachedType " << this->cached_name << " : resize_nr_of_derivatives() to " << new_size << std::endl;
        #endif
    }

    virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        if (cached_name == name) { 
            //std::cout << "matched"<< std::endl;
            return type_comparison<typename Expression<ResultType>::Ptr, typename Expression<Frame>::Ptr >::return_if_equal(argument);
        } else {
            return argument->subExpression_Frame(name);
        }
    }

    virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        if (cached_name == name) { 
            //std::cout << "matched"<< std::endl;
            return type_comparison<typename Expression<ResultType>::Ptr, typename Expression<Rotation>::Ptr >::return_if_equal(argument);
        } else {
            return argument->subExpression_Rotation(name);
        }
 
    }

    virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        if (cached_name == name) { 
            //std::cout << "matched"<< std::endl;
            return type_comparison<typename Expression<ResultType>::Ptr, typename Expression<Vector>::Ptr >::return_if_equal(argument);
        } else {
            return argument->subExpression_Vector(name);
        }
 
    }

    virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        if (cached_name == name) { 
            //std::cout << "matched"<< std::endl;
            return type_comparison<typename Expression<ResultType>::Ptr, typename Expression<Twist>::Ptr >::return_if_equal(argument);
        } else {
            return argument->subExpression_Twist(name);
        }
 
    }
    
    virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        if (cached_name == name) { 
            //std::cout << "matched"<< std::endl;
            return type_comparison<typename Expression<ResultType>::Ptr, typename Expression<Wrench>::Ptr >::return_if_equal(argument);
        } else {
            return argument->subExpression_Wrench(name);
        }
 
    }

    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        if (cached_name == name) { 
            //std::cout << "matched"<< std::endl;
            return type_comparison<typename Expression<ResultType>::Ptr, typename Expression<double>::Ptr >::return_if_equal(argument);
        } else {
            return argument->subExpression_Double(name);
        }
    }


    virtual typename Expression<ResultType>::Ptr clone() {
        typename Expression<ResultType>::Ptr expr(
            new CachedType<ResultType>( argument->clone(),this->cached_name)
        );
        return expr;
    }

    virtual void setInputValues(const std::vector<double>& values) {
        cached_value=false;
        fill_n(cached_deriv.begin(), deriv.size(), false);
        argument->setInputValues(values);
    } 

    virtual void setInputValue(int variable_number, double val) {
        cached_value=false;
        if (variable_number < (int)deriv.size()) {
            fill_n(cached_deriv.begin(), deriv.size(), false);
        }
        argument->setInputValue(variable_number,val);
    } 
    virtual void setInputValue(int variable_number, const Rotation& val) {
        cached_value=false;
        if (variable_number < (int)deriv.size()) {
            fill_n(cached_deriv.begin(), deriv.size(), false);
        }
        argument->setInputValue(variable_number,val);
    }

    virtual void debug_printtree() {
        std::cout << "cached(";
        argument->debug_printtree();
        std::cout << ")";
    }

    virtual void print(std::ostream& os) const {
        os << "cached(";
        argument->print(os);
        os << ")";
    }
    virtual void write_dotfile_init() {
        dot_already_written = false;
        argument->write_dotfile_init();
    }

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        if (!dot_already_written) {
            dot_already_written=true;
            thisnode=(size_t)this;
            of << "S"<<thisnode<<"[label=\"cached("
               << cached_name
               << ")\",shape=box,style=filled,fillcolor="
               << COLOR_CACHED << ",color=black]\n";
            pnumber argnode;
            argument->write_dotfile_update(of,argnode);
            of << "S"<<thisnode<<" -> " << "S"<<argnode<< "\n"; //rev
        } else {
            thisnode=(size_t)this;
        }
    }
    ~CachedType() {
        //#ifdef EG_LOG_CACHE
        //    std::cerr << "destructor CachedType("<<this->cached_name<<")"<<std::endl;
        //#endif
    }
};
/*
// template specialization that overrides the previous declaration when ResultType and T match:
template <>
inline typename Expression<Frame>::Ptr CachedType<Frame>::subExpression_Frame(const std::string& name) {
        if ( this->name == name ) {
            typename Expression<Frame>::Ptr result(this);
            return result;
        } else {
            return typename Expression<Frame>::Ptr();
        }
}
*/


/** utility function to create VariableType */
template<typename ResultType>
inline typename Expression<ResultType>::Ptr cached( typename Expression<ResultType>::Ptr argument ) {
   typename Expression<ResultType>::Ptr cach(
        new CachedType<ResultType>( argument,"" )
   );
   return cach;
}

/** utility function to create VariableType */
template<typename ResultType>
inline typename Expression<ResultType>::Ptr cached( const std::string& name,typename Expression<ResultType>::Ptr argument ) {
   typename Expression<ResultType>::Ptr cach(
        new CachedType<ResultType>( argument,name )
   );
   return cach;
}

/**
 * computing the numerical derivative for an expression tree
 * ( not using the derivative() function )
 * \param expr [in] expression to compute the numerical derivative.
 * \param towards_var [in]  variable number of the variable towards the derivative will be taken.
 * \param value [in]  value for the towards_var-th variable. 
 * \param h [in] interval over which to take the numerical derivative. 
 */
template<typename ResultType>
typename AutoDiffTrait<ResultType>::DerivType
inline numerical_derivative( typename Expression<ResultType>::Ptr expr, int towards_var, double value, double h=1E-7) {
    ResultType a,b;
    double val;
    val = value - h;
    expr->setInputValue(towards_var,val);
    a = expr->value();


    val = value + h;
    expr->setInputValue(towards_var,val);
    b = expr->value();

    return KDL::diff(a,b,1.0)/(2.0*h);
}


/**
 * output to an ostream:
 */
template<typename ResultType>
inline std::ostream& display( std::ostream& os, typename Expression<ResultType>::Ptr expr ) {
   os << "Value : ";
   os << expr->value() << "\n";
   for (int i=0;i<expr->number_of_derivatives();++i) {
       os << "Derivative towards var " << i << " : : " << expr->derivative(i) << "\n";
   }
   os << "\n";
   return os;
}

template <typename R>
class MakeConstantType: public UnaryExpression<R,R> {
public:
    typedef UnaryExpression<R,R> UnExpr;
    typedef typename AutoDiffTrait<R>::DerivType DerivType;

    MakeConstantType() {}

    MakeConstantType(typename UnExpr::ArgumentExpr::Ptr arg ):
        UnaryExpression<R,R>("make_constant",arg) {}

    virtual R value() {
        return this->argument->value();
    }
    virtual DerivType derivative(int i) {
        return AutoDiffTrait<R>::zeroDerivative();
    }
    virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
        return Constant(  AutoDiffTrait<R>::zeroDerivative());
    }
    virtual typename Expression<R>::Ptr clone() {
        typename Expression<R>::Ptr expr(
            new MakeConstantType( this->argument->clone() )
        );
        return expr;
    }
};

template<typename R>
inline typename Expression<R>::Ptr make_constant( const typename Expression<R>::Ptr& arg ) {
   typename Expression<R>::Ptr cnst(
        new MakeConstantType<R>( arg )
   );
   return cnst;
}

/**
 * This expressiongraph node represents a constant value corresponding to 
 * the initial value of an expression. 
 *
 * Caveat:
 *  - The initial value (time==0) should have occured before, otherwise its value is not cached.
 *    (this is satisfied when time is always increasing and starts from zero)
 *  - The resulting expression has derivatives zero, i.e. it is a constant.
 */
template<typename T>
class InitialValueType: public BinaryExpression<T,double,T> {
        T initial_value;
    public:
        typedef BinaryExpression<T,double,T> BinExpr;
        typedef typename AutoDiffTrait<T>::DerivType DerivType; 

        InitialValueType() {}

        InitialValueType(Expression<double>::Ptr time_var, typename Expression<T>::Ptr arg):
            BinExpr("initial_value", time_var, arg) {}

        InitialValueType(Expression<double>::Ptr time_var, typename Expression<T>::Ptr arg, const T& _initial_value):
            BinExpr("initial_value", time_var, arg), initial_value(_initial_value) {}

        virtual T value() {
            double time=this->argument1->value();
            if (time==0) {
                initial_value = this->argument2->value();
            }
            return initial_value; 
        } 
        
        virtual DerivType derivative(int i) {
            return AutoDiffTrait<T>::zeroDerivative();
        }
        virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
            return Constant(  AutoDiffTrait<T>::zeroDerivative());
        }
        virtual typename Expression<T>::Ptr clone() {
            typename Expression<T>::Ptr expr(
                new InitialValueType(this->argument1->clone(), this->argument2->clone(), initial_value)
            );
            return expr;
        }
};

template<typename R>
inline typename Expression<R>::Ptr initial_value( const typename Expression<double>::Ptr time, const typename Expression<R>::Ptr& arg ) {
   typename Expression<R>::Ptr e(
        new InitialValueType<R>( time,arg )
   );
   return e;
}

/**
 * - This class allows to optimize setInputValue() routines and optimizes
 * the use of the cached nodes.
 * - prepare should always be called before the addInput/addCached_xxx routines.
 * - addInput/addCached_xxx should be called before setInputValues(..) 
 * - setInputValue called on this object replaces all calls to the individual
 *   expression trees.
 *
 * Typical usage:
 * @code
 *   Expression<Frame> e1 = ....
 *   Expression<Vector> e2 = ....
 *    ...
 *   ExpressionOptimizer opt;
 *   opt.prepare(variable_list);
 *   e1->addToOptimizer(opt);
 *   e2->addToOptimizer(opt);
 *   ...
 *   opt.setInputValues(values);   
 * @endcode
 * 
 * the last call to opt.setInputValues(..) replaces calling e1.setInputValues(..) and e2.setInputValues(..)
 * It has the following effects:
 * - the expression tree does not need to be traversed to find input(..) nodes.
 * - calling e2.setInputValues(..) does not invalidate the cached expressions in e1.setInputValues(..)
 * - cached nodes are put in dirty state without traversing the expression tree.
 * - cached nodes that only depend on variables not in the variable_list will not be put in dirty state.  The optimizer
 *   supposes that variables not in the variable_list are constant.
 *
 * @warning: ExpressionOptimizer supposes that the references to the expression trees remain valid as long as the expression optimizer
 *           is used.  Due to technical reasons (an object does not know its own reference count), it does not increase the reference
 *           count in the smart pointer that is normally used to handle expression trees.
 *
 * The ExpressionOptimizer works as follows:
 *  - all Input objects that are relevant inside the expressions are registered.
 *  - all relevant Cached objects are also registered.
 *  - during setInputValue, all registered input values are set to the appropriate value
 *    and all registered Cached objects are invalidated.
 */
class ExpressionOptimizer {
    typedef std::list<InputType*>  ListInput;
    typedef std::vector<ListInput> Inputs;
    typedef std::set<int>          InputSet;
    typedef std::list<InputRotationType*> ListRotInput;
    typedef std::vector<ListRotInput>     RotInputs;
    typedef std::set<int>                 RotInputSet;


    Inputs                                inputs;          ///< for each input variable number, a list of InputType objects with this variable number
    std::vector<int>                      inputvarnr;      ///< input variable numbers

    RotInputs                             rotinputs;       ///< for each input variable number, a list of InputRotationType objects with this variable number
    std::vector<int>                      rotinputvarnr;   ///< input variable numbers corresponding to Rotation variables.

    InputSet                              inputset;        ///< std::set that contains all involved variable numbers (for scalar + Rotation)

    std::set<CachedExpression*>             cached;        ///< set of cached of objects that will have to be invalidated
    std::vector<CachedExpression*>          v_cached;      ///< list of cached of objects that will have to be invalidated

public:
        /**
         * configure the optimizer for input of these variable numbers.
         * @param inputvarnr a list of (integer) variable numbers (corresponding to scalar variables)
         */
        void prepare(const std::vector<int>& inputvarnr);

        /**
         * configure the optimizer for input of these variable numbers.
         * @param inputvarnr a list of (integer) variable numbers (corresponding to scalar variables)
         * @param rotinputvarnr a list of (integer) variable numbers (corresponding to Rotation variables)
         * \warning each Rotation variable takes 3 variable numbers.
         */
        void prepare(const std::vector<int>& inputvarnr, const std::vector<int>& rotinputvarnr);


        /**
         * registers an input object to be refreshed.
         */
        void addInput(InputType* obj);

        /**
         * registers an input object to be refreshed.
         */
        void addInput(InputRotationType* obj);

        /// registers a Cached object
        void addCached( CachedExpression* obj); 

        /**
         * set the input values for all of the involved expressions, the values given correspond to the
         * inputvarnr given with the prepare method.
         */
        void setInputValues(const std::vector<double>& values);
        /**
         * set the input values for all of the involved expressions (both scalar and Rotation), the values given correspond to the
         * inputvarnr and rotinputvarnr given with the prepare method.
         */
        void setInputValues(const std::vector<double>& values, const std::vector<Rotation>& rotvalues);

        /**
         * set the input values for all of the involved expressions, the values given correspond to the
         * inputvarnr given with the prepare method.
         */
        void setInputValues(const Eigen::VectorXd& values);
        
        /**
         * set the input values for all of the involved expressions (both scalar and Rotation), the values given correspond to the
         * inputvarnr and rotinputvarnr given with the prepare method.
         */
        void setInputValues(const Eigen::VectorXd& values, const std::vector<Rotation>& rotvalues);
};

inline void InputType::addToOptimizer(ExpressionOptimizer& opt) {
        //std::cout << "calling addinput " << std::endl;
        opt.addInput(this);
}

inline void InputRotationType::addToOptimizer(ExpressionOptimizer& opt) {
        //std::cout << "calling addinput " << std::endl;
        opt.addInput(this);
}

inline void CachedExpression::addToOptimizer(ExpressionOptimizer& opt) {
    opt.addCached(this);
}

template<typename T>
inline typename Expression<T>::Ptr checkConstant( typename Expression<T>::Ptr a ) {
        if (!a) {
            throw NullPointerException();
        }
        if (a->isConstant()) {
            return Constant( a->value() );
        } else {
            return a;
        }
        /*
        std::set<int> vset;
        a->getDependencies(vset);
        if (vset.empty()) {
            return Constant( a->value() );
        } else {
            return a;
        }*/
}

/*
inline ExpressionBase::Ptr checkConstant(typename ExpressionBase::Ptr a) {
        if (!a) {
            throw NullPointerException();
        }
        std::set<int> vset;
        a->getDependencies(vset);
        if (vset.empty()) {
            switch (a->getResultType()) {
                case ExpressionType::expression_double: {
                    Expression<double>::Ptr p = boost::dynamic_pointer_cast<Expression<double> >(a);
                    return Constant<double>( p->value() );
                    break;
                }
                case ExpressionType::expression_vector: {
                    Expression<Vector>::Ptr p = boost::dynamic_pointer_cast<Expression<Vector> >(a);
                    return Constant<Vector>( p->value() );
                    break;
                }
                case ExpressionType::expression_rotation: {
                    Expression<Rotation>::Ptr p = boost::dynamic_pointer_cast<Expression<Rotation> >(a);
                    return Constant<Rotation>( p->value() );
                    break;
                }
                case ExpressionType::expression_frame: {
                    Expression<Frame>::Ptr p = boost::dynamic_pointer_cast<Expression<Frame> >(a);
                    return Constant<Frame>( p->value() );
                    break;
                }
                case ExpressionType::expression_twist: {
                    Expression<Twist>::Ptr p = boost::dynamic_pointer_cast<Expression<Twist> >(a);
                    return Constant<Twist>( p->value() );
                    break;
                }
                case ExpressionType::expression_wrench: {
                    Expression<Wrench>::Ptr p = boost::dynamic_pointer_cast<Expression<Wrench> >(a);
                    return Constant<Wrench>( p->value() );
                    break;
                }
                case ExpressionType::expression_quaternion: {
                    Expression<Quaternion>::Ptr p = boost::dynamic_pointer_cast<Expression<Quaternion> >(a);
                    return Constant<Quaternion>( p->value() );
                    break;
                }
                default:
            }
        } else {
            return a;
        }
}
*/

template<typename T>
inline bool isConstant( typename Expression<T>::Ptr a) {
    if (!a) {
        throw NullPointerException();
    }
    return a->isConstant();
    /*
    std::set<int> vset;
    a->getDependencies(vset);
    return vset.empty();*/
}

inline bool isConstantZero( Expression<double>::Ptr a) {
    if (!a) {
        throw NullPointerException();
    }
    return a->isConstant() && (a->value()==0);
    /*std::set<int> vset;
    a->getDependencies(vset);
    return vset.empty() && (a->value()==0);*/
}

inline bool isConstantOne( Expression<double>::Ptr a) {
    if (!a) {
        throw NullPointerException();
    }
    return a->isConstant() && (1-EG_EPS_EQUAL <= a->value()) && (a->value() <= 1+EG_EPS_EQUAL);
    /*std::set<int> vset;
    a->getDependencies(vset);
    double eps = 1E-16;
    return vset.empty() && (1-eps <= a->value()) && (a->value() <= 1+eps);*/
}


/**
 * Explicite template instantiations:
 */
extern template class Expression<double>;
extern template class Expression<Vector>;
extern template class Expression<Frame>;
extern template class Expression<Rotation>;
extern template class Expression<Twist>;
extern template class Expression<Wrench>;
extern template class Expression<Quaternion>;

extern template class UnaryExpression<double,double>;
extern template class UnaryExpression<Vector,Vector>;
extern template class UnaryExpression<Frame,Frame>;
extern template class UnaryExpression<Rotation,Rotation>;
extern template class UnaryExpression<Twist,Twist>;
extern template class UnaryExpression<Wrench,Wrench>;

extern template class UnaryExpression<Quaternion,Quaternion>;
extern template class UnaryExpression<double,Quaternion>;
extern template class UnaryExpression<Vector,Quaternion>;
extern template class UnaryExpression<Quaternion,Vector>;
extern template class UnaryExpression<Quaternion,Rotation>;
extern template class UnaryExpression<Rotation,Quaternion>;

extern template class BinaryExpression<double,double,double>;
extern template class BinaryExpression<Vector,Vector,Vector>;
extern template class BinaryExpression<Rotation,Rotation,Rotation>;
extern template class BinaryExpression<Frame,Frame,Frame>;
extern template class BinaryExpression<Twist,Twist,Twist>;
extern template class BinaryExpression<Wrench,Wrench,Wrench>;
extern template class BinaryExpression<Quaternion,Quaternion,Quaternion>;

extern template class BinaryExpression<Quaternion,Quaternion,Vector>;
extern template class BinaryExpression<Quaternion,Vector,Quaternion>;
extern template class BinaryExpression<double,Quaternion,Quaternion>;
extern template class BinaryExpression<Vector,Quaternion,Vector>;
extern template class BinaryExpression<Quaternion,double,Vector>;


extern template class TernaryExpression<double,double,double,double>;
extern template class TernaryExpression<Vector,double,double,double>;
extern template class TernaryExpression<Rotation,double,double,double>;

extern template class QuaternaryExpression<double,double,double,double,double>;

extern template class FunctionType<double>;
extern template class FunctionType<Vector>;
extern template class FunctionType<Rotation>;
extern template class FunctionType<Frame>;
extern template class FunctionType<Twist>;
extern template class FunctionType<Wrench>;
extern template class FunctionType<Quaternion>;

extern template class ConstantType<double>;
extern template class ConstantType<Vector>;
extern template class ConstantType<Rotation>;
extern template class ConstantType<Frame>;
extern template class ConstantType<Twist>;
extern template class ConstantType<Wrench>;
extern template class ConstantType<Quaternion>;


extern template class CachedType<double>;
extern template class CachedType<Vector>;
extern template class CachedType<Rotation>;
extern template class CachedType<Frame>;
extern template class CachedType<Twist>;
extern template class CachedType<Wrench>;
extern template class CachedType<Quaternion>;


extern template class MakeConstantType<double>;
extern template class MakeConstantType<Vector>;
extern template class MakeConstantType<Rotation>;
extern template class MakeConstantType<Frame>;
extern template class MakeConstantType<Twist>;
extern template class MakeConstantType<Wrench>;
extern template class MakeConstantType<Quaternion>;



} // end of namespace KDL
#endif
