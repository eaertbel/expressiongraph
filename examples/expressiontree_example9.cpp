/**
 * \file expressiontree_example9.cpp
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

#include <expressiongraph/expressiontree.hpp>
#include <fstream>
#include <boost/timer.hpp>
#include <iostream>
#include <iomanip>

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

    // two robots,  1 m apart
    // note that chain is copied, so you can reuse it freely.
    Expression<Frame>::Ptr F_w_robot1 = cached<Frame>( Constant(Frame(Vector(0,0,0)))   * kinematic_chain(chain, 0 ));   
    Expression<Frame>::Ptr F_w_robot2 = cached<Frame>( Constant(Frame(Vector(1.0,0,0))) * kinematic_chain(chain, 7 ));   


    // ** constraint 1: distance between 2 robot ee is L,
    double L=1.0;
    Expression<double>::Ptr constraint1 = norm( origin(F_w_robot2) - origin(F_w_robot1) ) - Constant(L);

    // ** constraint 2&3: robot1 ee is on a cylinder, ee xy-plane lies on surface of cylinder.
    Frame F_w_cylinder( Rotation::Identity(), Vector(0.1,-0.95,0.0)); // center point at the base of the cylinder. 
    double  radius  = 0.1;                                      // radius of cylinder
    Expression<Frame>::Ptr tmpF_cylinder_robot1 = cached<Frame>( Constant(F_w_cylinder.Inverse()) * F_w_robot1 );
    // ee point on the cylinder:
    Expression<double>::Ptr constraint2         =  coord_x(origin(tmpF_cylinder_robot1))*coord_x(origin(tmpF_cylinder_robot1))  +
                                                   coord_y(origin(tmpF_cylinder_robot1))*coord_y(origin(tmpF_cylinder_robot1)) 
                                                   - Constant(radius*radius);
    // ee on surface of cylinder 
    Expression<double>::Ptr constraint3         = dot( unit_z( rotation(F_w_robot1) ), Constant(F_w_cylinder.M.UnitZ()));

    std::vector<double> joints(constraint1->number_of_derivatives()); 
    // specify input values:  
    joints[0] = 0;
    joints[1] = 20*deg2rad;
    joints[2] = 10*deg2rad;
    joints[3] = 30*deg2rad;
    joints[4] = 0;
    joints[5] = 30*deg2rad;
    joints[6] = 0;
    joints[7] = 0;
    joints[8] = 20*deg2rad;
    joints[9] = 0;
    joints[10] = 30*deg2rad;
    joints[11] = 0;
    joints[12] = 30*deg2rad;
    joints[13] = 0;

    cout << "--- Distance constraint ---\n";
    constraint1->setInputValues( joints );
    display<double>(cout,constraint1);
    cout << "\n";
    cout << "--- Cylinder constraint 2 ---\n";
    constraint2->setInputValues( joints );
    display<double>(cout,constraint2);
    cout << "--- Cylinder constraint 3 ---\n";
    constraint3->setInputValues( joints );
    display<double>(cout,constraint3);

    cout << "--- dependencies ---- \n";
    std::set<int> varset;
    constraint3->getDependencies(varset);     
    for (std::set<int>::iterator it=varset.begin(); it!= varset.end(); ++it) {
        cout << *it << "\t";
    }
    cout << "\n";
    std::set<int> varset2;
    constraint1->getDependencies(varset2);     
    for (std::set<int>::iterator it=varset2.begin(); it!= varset2.end(); ++it) {
        cout << *it << "\t";
    }
    cout << "\n";
	
	return 0;
}
