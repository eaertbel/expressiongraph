#ifndef EXPRESSIONGRAPH_REFVARIABLETYPE_HPP
#define EXPRESSIONGRAPH_REFVARIABLETYPE_HPP

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
#include <Eigen/Dense>
#include <map>
#include <vector>

namespace KDL {

/**
 * This class is a variable expressiongraph node. 
 * 
 * It stores a pointer to a value and a derivative of arbitrary geometric type,
 * and it returns this value and derivative when it is used as an expression graph node. 
 *
 * This class is useful to get external information inside an expression graph.
 *
 * \warning logically it depends on the (during construction) specified variable numbers. 
 *          It assumes it remains unchanged, if setInputValue is not called on the 
 *          specified variable numbers.  It will return a zero derivative for other variable numbers,
 *          beside the variable numbers specified during construction.  
 *          When asked for derivatives, it returns the corresponding
 *          columns of its Jacobian data structure.  The data passed with setInputValue() calls is not used.
 *          It is however necessary to call setInputValue, in order to achieve appropriate behavior of CachedType and
 *          the ExpressionOptimizer.
 * \warning getScalarDependencies and getRotDependencies return both an empty set, since they are used
 *          to determine setInputValue(..) calls.  getDependencies returns the during construction specified variable numbers,
 *          since getDependencies is used for dependencies and derivatives.
 * \warning The semantics of this node are:  a value with a given first order derivative and with all higher order derivatives equal to zero.
 *          The behavior of derivativeExpression(..) is according to these semantics.
 */
template <typename ResultType>
class VariableType: public FunctionType<ResultType> {
public:
    typedef typename AutoDiffTrait<ResultType>::DerivType       DerivType;
    typedef ResultType                                          ValueType;
    typedef typename std::vector<DerivType>                     JacobianType;
    
public:
    typedef typename boost::shared_ptr<VariableType> Ptr;

    /// smart pointer to the current value:
    ValueType val;
    /// smart pointer to the current Jacobian:
    JacobianType deriv;
    /// ndx[i] is the variable number that corresponds to the ith column of the Jacobian.
    /// This is the same format as in setInputValues(ndx,vals)
    std::vector<int> ndx;
    /// maps a variable_number ( in ndx) to a column of the Jacobian deriv. 
    std::map<int, int> ndx_map;

    int maxderiv;

    VariableType() {} 
   
    VariableType(const std::vector<int>& _ndx):
        FunctionType<ResultType>("variable"),deriv(_ndx.size()),ndx(_ndx) {
            maxderiv = 0;
            for (size_t i=0;i<ndx.size();++i) {
                ndx_map[ ndx[i] ] = i;
                maxderiv = std::max( maxderiv , ndx[i] );
            }
            maxderiv += 1;
    }

    /**
     * sets the value of this Variable to _val.
     */
    virtual void setValue(const ValueType& _val) {
        val = _val;
    }
    /**
     * sets the i-th component of the Jacobian of this Variable to _d.
     * This corresponds to the derivative towards the variable ndx[i].
     */
    virtual void setJacobian(int i, const DerivType& _d) {
        deriv[i] = _d;
    }  

    virtual void setInputValue(int variable_number, double val) {
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
    }


    virtual void setInputValues(const std::vector<double>& values) {
    }

    virtual ResultType value() {
        return val;
    }

    virtual void getDependencies(std::set<int>& varset) {
        for (size_t i=0;i<ndx.size();++i) {
            varset.insert( ndx[i]);
        }
    }
    virtual void getScalarDependencies(std::set<int>& varset) {
    }
    
    virtual void getRotDependencies(std::set<int>& varset) {
    }

    virtual DerivType derivative(int i) {
        std::map<int,int>::iterator it = ndx_map.find(i);     
        if (it!=ndx_map.end()) {
            return deriv[it->second];
        } else {
            return AutoDiffTrait<DerivType>::zeroDerivative();
        }
    }

    virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
        std::map<int,int>::iterator it = ndx_map.find(i);     
        if (it!=ndx_map.end()) {
            return Constant( deriv[it->second] );
        } else {
            return Constant(AutoDiffTrait<DerivType>::zeroDerivative());
        } 
    }

    virtual int number_of_derivatives() {
        return maxderiv;
    }

    /**
     * A clone does not clone the value and derivative that is this node refers to.
     * Cloning also the value and derivative would make no sense, because it would point to a value that nobody can change. 
     */
    virtual typename Expression<ResultType>::Ptr clone() {
        typename Expression<ResultType>::Ptr expr(
            new VariableType<ResultType>( ndx)
        );
        return expr;
    }

    virtual ~VariableType() {}
};


/**
 * typically when you create a variable, you'll store a separate copy of the VariableType pointer and
 * cast it to Expression<T>::Ptr for further use in the expression.
 */
template <typename T>
inline typename VariableType<T>::Ptr Variable( const std::vector<int>& ndx) 
{
        typename KDL::VariableType<T>::Ptr tmp(
            new VariableType<T>(  ndx )
        );
        return tmp;
}


/**
 * typically when you create a variable, you'll store a separate copy of the VariableType pointer and
 * cast it to Expression<T>::Ptr for further use in the expression.
 */
template <typename T>
inline typename VariableType<T>::Ptr Variable( int startndx, int nrofderiv)
{
        std::vector<int> ndx(nrofderiv);
        for (int i=0;i<nrofderiv;++i) {
            ndx[i] = startndx+i;
        }
        return Variable<T>(ndx);
}



template<typename ValueType>
class Callback {
    public:
        typedef typename AutoDiffTrait<ValueType>::DerivType       DerivType;
        typedef typename boost::shared_ptr<Callback> Ptr;
        virtual void compute(ValueType& value, std::vector<DerivType>& jacobian) = 0;
        virtual Callback::Ptr clone()        = 0;
        virtual ~Callback() {}
};

/**
 * This class is a variable expressiongraph node that uses a callback function
 * to get its values and jacobians.
 * 
 * It stores a pointer to a value and a derivative of arbitrary geometric type,
 * and it returns this value and derivative when it is used as an expression graph node. 
 *
 * This class is useful to get external information inside an expression graph.
 *
 * 
 *
 * \warning It is necessary to call setInputValue such that the node is invalidated and a new computation is requested
 *          from the callback function, in order to achieve appropriate behavior of CachedType and
 *          the ExpressionOptimizer.  A call for any input variable number invalidates the cache.
 * \warning The semantics of this node are:  a value with a given first order derivative and with all higher order derivatives equal to zero.
 *          The behavior of derivativeExpression(..) is according to these semantics.
 *
template <typename ResultType>
class CallbackNode: public FunctionType<ResultType>, public CachedExpression {
public:
    typedef typename AutoDiffTrait<ResultType>::DerivType       DerivType;
    typedef ResultType                                          ValueType;
    typedef typename std::vector<DerivType>                     JacobianType;
    
public:
    typedef typename boost::shared_ptr<CallbackNode> Ptr;

    /// smart pointer to the current value:
    ValueType val;
    /// smart pointer to the current Jacobian:
    JacobianType deriv;
    /// ndx[i] is the variable number that corresponds to the ith column of the Jacobian.
    /// This is the same format as in setInputValues(ndx,vals)
    std::vector<int> ndx;
    /// maps a variable_number ( in ndx) to a column of the Jacobian deriv. 
    std::map<int, int> ndx_map;

    int maxderiv;
    typename Callback<ResultType>::Ptr cb;
    bool cached_value;

    CallbackNode() {} 
  
    CallbackNode(typename Callback<ResultType>::Ptr _cb,std::vector<int> _ndx):
        FunctionType<ResultType>("callbacknode"),deriv(_ndx.size()),ndx(_ndx),cb(_cb) {
            maxderiv = 0;
            for (size_t i=0;i<ndx.size();++i) {
                ndx_map[ ndx[i] ] = i;
                maxderiv = std::max( maxderiv , ndx[i] );
            }
            maxderiv += 1;
            invalidate_cache();
    }

   virtual void invalidate_cache() {
        cached_value = false;
    }


    virtual void setInputValue(int variable_number, double val) {
        invalidate_cache();
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        invalidate_cache();
    }

    virtual void setInputValues(const std::vector<double>& values) {
        invalidate_cache();
    }

    virtual ResultType value() {
        if (!cached_value) {
            cb->compute(val,deriv);
            cached_value = true;
        }
        return val;
    }

    virtual void getDependencies(std::set<int>& varset) {
        for (size_t i=0;i<ndx.size();++i) {
            varset.insert( ndx[i]);
        }
    }

    virtual DerivType derivative(int i) {
        if (!cached_value) {
            cb->compute(val,deriv);
            cached_value = true;
        }
        std::map<int,int>::iterator it = ndx_map.find(i);     
        if (it!=ndx_map.end()) {
            return deriv[it->second];
        } else {
            return AutoDiffTrait<DerivType>::zeroDerivative();
        }
    }

    virtual typename Expression<DerivType>::Ptr derivativeExpression(int i) {
        if (!cached_value) {
            cb->compute(val,deriv);
            cached_value = true;
        }
        std::map<int,int>::iterator it = ndx_map.find(i);     
        if (it!=ndx_map.end()) {
            return Constant( deriv[it->second] );
        } else {
            return Constant(AutoDiffTrait<DerivType>::zeroDerivative());
        } 
    }

    virtual int number_of_derivatives() {
        return maxderiv;
    }

    **
     * A clone does not clone the value and derivative that is this node refers to.
     * Cloning also the value and derivative would make no sense, because it would point to a value that nobody can change. 
     *
    virtual typename Expression<ResultType>::Ptr clone() {
        typename Expression<ResultType>::Ptr expr(
            new CallbackNode<ResultType>(cb->clone(), ndx)
        );
        return expr;
    }
};


**
 * typically when you create a variable, you'll store a separate copy of the CallbackNode pointer and
 * cast it to Expression<T>::Ptr for further use in the expression.
 *
template <typename T>
inline typename CallbackNode<T>::Ptr create_callback_node( std::vector<int> ndx) 
{
        typename KDL::CallbackNode<T>::Ptr tmp(
            new CallbackNode<T>(  ndx )
        );
        return tmp;
}


**
 * typically when you create a variable, you'll store a separate copy of the CallbackNode pointer and
 * cast it to Expression<T>::Ptr for further use in the expression.
 *
template <typename T>
inline typename CallbackNode<T>::Ptr create_callback_node( int startndx, int nrofderiv)
{
        std::vector<int> ndx(nrofderiv);
        for (int i=0;i<nrofderiv;++i) {
            ndx[i] = startndx+i;
        }
        return CallbackNode<T>(ndx);
}
*/

}; // namespace KDL
#endif

