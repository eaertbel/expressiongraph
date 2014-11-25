/**
 * \file expressiontree_example6.cpp
 * \brief example of using ExpressionTree in combination with KDL::Chain.
 *
 * \Author: Sept. 2012, Erwin Aertbelien 
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

#include <kdl/expressiontree.hpp>
#include <fstream>
#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>

KDL::Twist numerical_derivative_jacobian( boost::shared_ptr<KDL::Expression_Chain>& expr, std::vector<double> values, int first,int second, double dt=1E-7) {
    KDL::Twist a,b;
    double val = values[second];

    values[second] = val-dt;
    expr->setInputValues(values);
    expr->value(); // always call value() before derivative() 
    a = expr->derivative(first); 

    values[second] = val+dt;
    expr->setInputValues(values);
    expr->value();
    b = expr->derivative(first); 

    values[second] = val;
    return diff(a,b,2*dt);
}

KDL::Twist numerical_timederivative_jacobian( boost::shared_ptr<KDL::Expression_Chain>& expr, std::vector<double> values, std::vector<double> values_dot,int column, double dt=1E-7) {
    std::vector<double> val(values.size());
    KDL::Twist a,b;
    for (int i=0;i<val.size();++i) {
        val[i] = values[i] - dt*values_dot[i];
    }

    expr->setInputValues(val);
    expr->value(); // always call value() before derivative() 
    a = expr->derivative(column); 


    for (int i=0;i<val.size();++i) {
        val[i] = values[i] + dt*values_dot[i];
    }
    expr->setInputValues(val);
    expr->value();
    b = expr->derivative(column); 

    return diff(a,b,2*dt);
}

/**
 * computing derivatives of a Jacobian: dJ/dq
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;
    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;

    /**
     * The same chain is defined twice. Note that the definitions of chain
     * are different:  
     *    chain: first a joint rotation and then the frame transformation.
     *    frame in the expression tree: first translation and then rotation4.
     */
    Chain chain;
    chain.addSegment(Segment("Segment 0", Joint("Joint 0", Joint::RotZ),Frame(Vector(0,0,L1))));
    chain.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotZ),Frame(Vector(0,0,L2))));
    chain.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotZ),Frame(Vector(0,0,L3))));
    chain.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotX),Frame(Vector(0,0,L4))));
    chain.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),Frame(Vector(0,0,0))));

    // you cannot use the kinematic_chain helper function because you need the specific type
    // of Expression_chain instead of the base class Expression<>
    boost::shared_ptr<Expression_Chain> kinchain( new Expression_Chain(chain, 0 ));   
    std::vector<double> joints(kinchain->number_of_derivatives()); 
    // specify input values:  
    int N=10;
    joints[0] = 0*M_PI*0.8/N;
    joints[1] = 1*M_PI*0.8/N;
    joints[2] = 2*M_PI*0.8/N;
    joints[3] = 3*M_PI*0.8/N;
    joints[4] = 4*M_PI*0.8/N;
    joints[5] = 5*M_PI*0.8/N;
    joints[6] = 6*M_PI*0.8/N;
    kinchain->setInputValues( joints );
    cout << "Value : \n";
    Frame tf = kinchain->value();
    cout << tf << "\n";
    cout << "Columns of the Jacobian \n";
    for (int i=0;i<kinchain->number_of_derivatives();++i) {
        cout << "Derivative " << i << " : " << kinchain->derivative(i) << endl;
    }

    cout << "\n\n\nComputing the partial derivative of the Jacobian\n";
    for (int i=0;i<kinchain->number_of_derivatives();++i) {
        for (int j=0;j<kinchain->number_of_derivatives();++j) {
            cout << "Derivative towards  " << i << " and then towards " << j << " : " << kinchain->derivative(i,j) << endl;
            cout << "Numerically         " << i << " towards          " << j << " : " 
                 << numerical_derivative_jacobian(kinchain,joints,i,j,1E-7)  << endl;
        }
        std::cout << endl;
    }

   
    cout << "\n\n\nComputing the time derivative of the Jacobian\n";
    std::vector<double> joints_dot(kinchain->number_of_derivatives());
    joints_dot[0] = 0*M_PI*0.1;
    joints_dot[1] = 1*M_PI*0.1;
    joints_dot[2] = 2*M_PI*0.1;
    joints_dot[3] = 3*M_PI*0.1;
    joints_dot[4] = 4*M_PI*0.1;
    joints_dot[5] = 5*M_PI*0.1;
    joints_dot[6] = 6*M_PI*0.1;
    kinchain->value();  // for all ExpressionTree objects, always call value(), even if you do not need it.
    for (int i=0;i<kinchain->number_of_derivatives();++i) {
        cout << "J_dot column " << i << " : " << kinchain->derivative_dot(i,joints_dot) << "\n";
        cout << "numerical    " << i << " : " << numerical_timederivative_jacobian(kinchain,joints,joints_dot,i) << "\n";
    } 
	return 0;
}
