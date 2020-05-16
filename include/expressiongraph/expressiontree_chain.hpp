/*
 * expressiontree_chain.hpp
 *
 *  Created on: Sept, 2012
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

#ifndef KDL_EXPRESSIONTREE_CHAIN_HPP
#define KDL_EXPRESSIONTREE_CHAIN_HPP

#include <expressiongraph/expressiontree_expressions.hpp>
#include <kdl/jntarray.hpp>
#include <kdl/chain.hpp>

namespace KDL {

/**
 * provides a mapping of a KDL::Chain to an ExpressionTree.
 */
class Expression_Chain:
	public FunctionType<Frame>
{
    Chain chain;
    std::vector<int> jointndx_to_segmentndx;
    std::vector<double> jval;
    std::vector<Frame> T_base_jointroot;
    std::vector<Frame> T_base_jointtip;
    Frame              T_base_head;
    std::vector<Twist> jacobian;
    std::vector<bool>  cached_deriv;
    bool cached;
    int _number_of_derivatives;
    int index_of_first_joint;
    void initialize();
public:
    /**
     * - provide a chain, This chain will be copied and used within this Expression chain object.
     * - provide a mapping of the joints of the kinematic chain to the Expression tree variable number.
     * - also provides the default value of a joint, it is not filled in using jointndx_to_varndx, this
     *   allows to keep certain joints constant of a kinematic chain.
     * - This class caches its results until a new call to setInputValue(s) (no need to use CachedType for
     *   this functionality).
     * - The jointndx_to_varndx mechanism is parallel to the mechanism already built-in into chain, but more
     *   powerfull, and is necessary to combine different chains in an expression.  Additional constructors 
     *   are/will be provided for ease of use.
     * - The constructors are not real-time.
     * - This classes caches its result, no need to add a CachedType node.  Can efficiently be called
     *   multiple times.
     *
     * See documentation of other constructor.
     * uses a jointndx_to_varndx mapping that maps the joints of a chain
     * in a sequence starting at varndx_of_first_joint.
     */   
    Expression_Chain( const Chain& _chain, int varndx_of_first_joint );

    virtual void setInputValues(const std::vector<double>& values);

    virtual void setInputValue(int var, double value);
    virtual void setInputValue(int var, const Rotation& value) {}

	virtual Frame value();

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        assert( 0 /* for now chain cannot be used together with optimization */ );
    }

    virtual void getDependencies(std::set<int>& varset) {
        for (size_t i=0;i<chain.getNrOfJoints();++i) {
            varset.insert(index_of_first_joint+i);
        }
    }
    virtual void getScalarDependencies(std::set<int>& varset) {
        getDependencies(varset);
    }
    virtual void getRotDependencies(std::set<int>& varset) {
    }

    virtual void update_variabletype_from_original() {}

    /**
     * This class caches the computation of the derivative
     */
	virtual Twist derivative(int var_ndx);

    /**
     * compute 2nd derivative, first deriving towards first_var, then towards second_var
     * An additional method specific for Expression_Chain.
     * The second derivative is not cached.
     */
    virtual Twist derivative(int first_var,int second_var);

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    /**
     * @param column: column to gives back.
     * @param jval_dot: time derivative of the joint variables
     *                  not specified joint variable velocities are
     *                  zero.
     */
    virtual Twist derivative_dot(int column, const std::vector<double>& jval_dot);
   
    virtual int number_of_derivatives() { 
        return _number_of_derivatives;
    }

    virtual void resize_nr_of_derivatives() {
    }

    virtual  Expression_Chain::Ptr clone();

    //virtual void write_dotfile_helper(std::ostream& of, size_t& thisnode, size_t& counter);


};

class Expression_Chain_Derivative:
    public FunctionType<Twist>
{
    boost::shared_ptr<Expression_Chain> argument;
    int var_ndx;
protected:
    Expression_Chain_Derivative(boost::shared_ptr<Expression_Chain> arg, int i);
public:

    virtual void setInputValues(const std::vector<double>& values);
    virtual void setInputValue(int var, double value);
    virtual void setInputValue(int var, const Rotation& value) {}
	virtual Twist value();

	virtual Twist derivative(int var_ndx);

    virtual Expression<Twist>::Ptr derivativeExpression(int i);

    virtual int number_of_derivatives() { 
        return argument->number_of_derivatives();
    }
    virtual void resize_nr_of_derivatives() {}
    virtual void update_variabletype_from_original() {}

    virtual Expression<Twist>::Ptr clone();

    friend class Expression_Chain;
};

Expression<Frame>::Ptr kinematic_chain(const Chain& chain, int index_of_first_joint );





#if 0
class Expression_Manipulability:
    public Expression<double>
{
	Expression_Manipulability(Expression_Chain::Ptr chainexpr);

    Expression_Chain( const Chain& _chain, int varndx_of_first_joint );

    virtual void setInputValues(const std::vector<double>& values);

	virtual Frame value();

	virtual Twist derivative(int var_ndx);

    virtual Twist derivative(int first_var,int second_var);

    virtual Twist derivative_dot(int column, const std::vector<double>& jval_dot);

    virtual int number_of_derivatives() { 
        return _number_of_derivatives;
    }

    virtual typename Expression_Chain::Ptr clone();


};
#endif

}; // namespace KDL
#endif
