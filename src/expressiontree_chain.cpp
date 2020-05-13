/*
 * expressiontree_chain.cpp
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

#include <expressiongraph/expressiontree_chain.hpp>
#include <iostream>
namespace KDL {

Expression<Frame>::Ptr kinematic_chain(const Chain& chain, int index_of_first_joint ) {
	Expression<Frame>::Ptr expr( new Expression_Chain( chain, index_of_first_joint ) );
    return expr;
}



Expression_Chain::Expression_Chain( const Chain& _chain, int _index_of_first_joint ):
    FunctionType("kinematic_chain"),
    chain(_chain),
    jointndx_to_segmentndx( _chain.getNrOfJoints()),
    jval( _chain.getNrOfJoints() ),
    T_base_jointroot( _chain.getNrOfJoints()),
    T_base_jointtip(  _chain.getNrOfJoints()),
    jacobian( _chain.getNrOfJoints() ),
    cached_deriv( _chain.getNrOfJoints() ),
    cached(false),
    index_of_first_joint(_index_of_first_joint)
{
    using namespace std;
    fill(jval.begin(),jval.end(),0.0);
    fill(cached_deriv.begin(),cached_deriv.end(),false);
    _number_of_derivatives=index_of_first_joint+_chain.getNrOfJoints();
}

void Expression_Chain::setInputValues(const std::vector<double>& values) {
    int index_of_last_joint = index_of_first_joint+chain.getNrOfJoints()-1;
    cached = false;
    int N = min( values.size()-1 , index_of_last_joint );
    for (int i=index_of_first_joint;i<=N;++i) {
        jval[i-index_of_first_joint] = values[i];
    }
    // redo all derivatives, not only the ones who's value is changed:
    fill(cached_deriv.begin(),cached_deriv.end(),false);
}


void Expression_Chain::setInputValue(int var, double value) {
    if (var < index_of_first_joint) return;
    if (var >= index_of_first_joint + (int)chain.getNrOfJoints() ) return;
    cached = false;
    jval[var-index_of_first_joint] = value;
    // redo all derivatives, not only the ones who's value is changed:
    fill(cached_deriv.begin(),cached_deriv.end(),false);

}
Frame Expression_Chain::value() {
    if (cached) {
        return T_base_head;
    }
    unsigned int jointndx=0;
    T_base_head = Frame::Identity(); // frame w.r.t. base of head
    for (unsigned int i=0;i<chain.getNrOfSegments();i++) {
        Segment segment = chain.getSegment(i);
        if (segment.getJoint().getType()!=Joint::None) {
            T_base_jointroot[jointndx] = T_base_head;
            T_base_head                = T_base_head * segment.pose(jval[jointndx]);
            T_base_jointtip[jointndx]  = T_base_head;
            jointndx_to_segmentndx[jointndx] = i;
            jointndx++;
        } else {
            T_base_head = T_base_head * segment.pose(0.0);
        }
    }
    cached = true;
    return T_base_head;
}

Twist Expression_Chain::derivative(int var_ndx) {
    int jointndx = var_ndx - index_of_first_joint;
    if (jointndx < 0) {
        return Twist::Zero();
    }
    if (jointndx >= (int)chain.getNrOfJoints() ) {
        return Twist::Zero();
    }
    if (cached_deriv[jointndx]) {
        return jacobian[jointndx];
    }
    int segmentndx = jointndx_to_segmentndx[ jointndx ];
    Segment segment=chain.getSegment(segmentndx);
    jacobian[jointndx]     = ( T_base_jointroot[jointndx].M * segment.twist(jval[jointndx],1.0) ).RefPoint( T_base_head.p - T_base_jointtip[jointndx].p);
    cached_deriv[jointndx] = true;
    return jacobian[jointndx];
}
/**
 *
    inline KDL::Twist jacobian_derivative(int i,int j) {
        if (j<=i) {
            KDL::Twist t;
            KDL::Vector omega_j (KDL::Vector(jac(3,j),jac(4,j),jac(5,j)));
            KDL::Vector vel_i(KDL::Vector(jac(0,i),jac(1,i),jac(2,i)));
            KDL::Vector omega_i(KDL::Vector(jac(3,i),jac(4,i),jac(5,i)));
            t.vel =  omega_j * vel_i;
            t.rot =  omega_j * omega_i;
            return t;
        } else {
            KDL::Twist t;
            KDL::Vector vel_j(KDL::Vector(jac(0,j),jac(1,j),jac(2,j)));
            KDL::Vector omega_i(KDL::Vector(jac(3,i),jac(4,i),jac(5,i)));
            t.vel =  -vel_j * omega_i;
            t.rot =  KDL::Vector::Zero();
            return t;
        }
    }
*/

Twist Expression_Chain::derivative(int first_var,int second_var) {
    int jointndx1 = first_var - index_of_first_joint;
    if (jointndx1 < 0) {
        return Twist::Zero();
    }
    if (jointndx1 >= (int)chain.getNrOfJoints() ) {
        return Twist::Zero();
    }
    int jointndx2 = second_var - index_of_first_joint;
    if (jointndx2 < 0) {
        return Twist::Zero();
    }
    if (jointndx2 >= (int)chain.getNrOfJoints() ) {
        return Twist::Zero();
    }
    if (!cached_deriv[jointndx1]) derivative(first_var); 
    if (!cached_deriv[jointndx2]) derivative(second_var); 
    if (jointndx2<=jointndx1) {
        KDL::Twist t;
        KDL::Vector omega_j (jacobian[jointndx2].rot);
        KDL::Vector vel_i(jacobian[jointndx1].vel);
        KDL::Vector omega_i(jacobian[jointndx1].rot);
        t.vel =  omega_j * vel_i;
        t.rot =  omega_j * omega_i;
        return t;
    } else {
        KDL::Twist t;
        KDL::Vector vel_j(jacobian[jointndx2].vel);
        KDL::Vector omega_i(jacobian[jointndx1].rot);
        t.vel =  -vel_j * omega_i;
        t.rot =  KDL::Vector::Zero();
        return t;
    }

}

Twist Expression_Chain::derivative_dot(int column, const std::vector<double>& jval_dot) {
    int jointndx = column - index_of_first_joint;
    if (jointndx < 0) {
        return Twist::Zero();
    }
    if (jointndx >= (int)chain.getNrOfJoints() ) {
        return Twist::Zero();
    }
    if (!cached_deriv[jointndx]) derivative(column); 
    Twist t; 
    for (int i=0;i<(int)jval_dot.size();++i) {
        t += jval_dot[i]*derivative(column,index_of_first_joint+i);
    } 
    return t;
}

Expression<Frame>::Ptr Expression_Chain::clone() {
    Expression<Frame>::Ptr expr(
        new Expression_Chain( chain, index_of_first_joint )
    );
    return expr;
}

Expression<Twist>::Ptr Expression_Chain::derivativeExpression(int i) 
{
    boost::shared_ptr<Expression_Chain> chain( this );
    Expression<Twist>::Ptr expr( new Expression_Chain_Derivative(chain,i));
    return expr;
}
/*
void Expression_Chain::write_dotfile_helper(std::ostream& of, size_t& thisnode,size_t& counter) {
        thisnode=counter++;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor=\"#FF7400\",color=black]\n";
    }
*/


Expression_Chain_Derivative::Expression_Chain_Derivative(boost::shared_ptr<Expression_Chain> arg, int i):
    FunctionType<Twist>("chain_derivative"),
    argument(arg),
    var_ndx(i) {
}

void Expression_Chain_Derivative::setInputValues(const std::vector<double>& values) {
    argument->setInputValues(values);
}

void Expression_Chain_Derivative::setInputValue(int var, double value) {
    argument->setInputValue(var,value);
}

Twist Expression_Chain_Derivative::value() {
    return argument->derivative(var_ndx);
}

Twist Expression_Chain_Derivative::derivative(int i) {
    return argument->derivative(var_ndx,i);
}

Expression<Twist>::Ptr Expression_Chain_Derivative::derivativeExpression(int i) {
    throw NotImplementedException();
}

Expression<Twist>::Ptr Expression_Chain_Derivative::clone() {
    throw NotImplementedException();
    //"Implementation of Expression_Chain_Derivative::clone() is not correct and preliminary");
    //Expression<Twist>::Ptr expr( new Expression_Chain_Derivative( argument, var_ndx) ); 
    //return expr;  
} 
}; // end of namespace KDL
