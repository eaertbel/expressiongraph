/*
 * expressiontree_example11.cpp
 *
 *  Created on: Sept., 2012
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

#include <kdl/expressiontree.hpp>

/*
 * Example that demonstrates the evaluation of a simple expression in doubles for
 * both value() and derivative()
 * and the use of derivativeExpression() to compute a derivative to an arbitrary order.
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
    using namespace std;
    cout << " -----------------------------------------------\n";
    cout << "Specifying a kinematic chain in two different ways\n";
    cout << " -----------------------------------------------\n";

    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;
    Chain chain;
    chain.addSegment(Segment("Segment 0", Joint("Joint 0", Joint::RotZ),Frame(Vector(0,0,L1))));
    chain.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotZ),Frame(Vector(0,0,L2))));
    chain.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotZ),Frame(Vector(0,0,L3))));
    chain.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotX),Frame(Vector(0,0,L4))));
    chain.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),Frame(Vector(0,0,0))));

    Expression<Frame>::Ptr kinchain1 = kinematic_chain( chain, 0 );   
    Expression<Frame>::Ptr kinchain2 = cached<Frame>(
        frame( rot_z(input(0)),  Constant(Vector(0,0,L1))) *
        frame( rot_x(input(1)) ) *
        frame( rot_z(input(2) ) ) *
        frame( rot_x(input(3) ), Constant(Vector(0,0,L2))) *
        frame( rot_z(input(4) ), Constant(Vector(0,0,L3))) *
        frame( rot_x(input(5) ) ) *
        frame( rot_z(input(6) ), Constant(Vector(0,0,L4)) )
    );

    cout << " -----------------------------------------------\n";
    cout << "Specifying joint values\n";
    cout << " -----------------------------------------------\n";
    std::vector<double> joints(7);
    int N=10;
    joints[0] = 0*M_PI*0.8/N;
    joints[1] = 1*M_PI*0.8/N;
    joints[2] = 2*M_PI*0.8/N;
    joints[3] = 3*M_PI*0.8/N;
    joints[4] = 4*M_PI*0.8/N;
    joints[5] = 5*M_PI*0.8/N;
    joints[6] = 6*M_PI*0.8/N;
 

    cout << " -----------------------------------------------\n";
    cout << "evaluating:\n";
    cout << " -----------------------------------------------\n";
    kinchain1->setInputValues(joints);
    kinchain2->setInputValues(joints);
    display<Frame>(cout,kinchain2);
    for (int i=0;i<7;++i) {
        Expression<Twist>::Ptr deriv = kinchain2->derivativeExpression(i);
        cout << "-"<< endl;
        kinchain2->debug_printtree();
        deriv->setInputValues(joints);
        cout << "Derivative towards " << i << endl;
        display<Twist>(cout,deriv);
    }
	return 0;
}
