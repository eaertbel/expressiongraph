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

    Chain chain;
    chain.addSegment(Segment("Segment 0", Joint("Joint 0", Joint::RotZ),Frame(Vector(0,0,L1))));
    chain.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotZ),Frame(Vector(0,0,L2))));
    chain.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotZ),Frame(Vector(0,0,L3))));
    chain.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotX),Frame(Vector(0,0,L4))));
    chain.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),Frame(Vector(0,0,0))));

    Expression<Frame>::Ptr F_w_robot1 = kinematic_chain(chain, 0 );   
    Frame F_w_plane( Rotation::Identity(), Vector(0,0,0.5) );

    cout << "--------------------------------------" << endl;
    cout << " A constraint for the distance to a plane " << endl;
    cout << " along a beam along the z-axis at the end effector " << endl;
    cout << "--------------------------------------" << endl;

    cout << "------ Constraint, directly computed " << endl;
    // signed distance to plane divided by cos(angle), 
    // problematic if beam // to plane. 
    Expression<double>::Ptr distance = dot( Constant(F_w_plane.M.UnitZ()),  Constant( F_w_plane.p) - origin(F_w_robot1) ) /
                                  dot( Constant(F_w_plane.M.UnitZ()),  unit_z(rotation(F_w_robot1) ) ); 
    display<double>(cout,distance);

    cout << "------ Constraint, iTasc style " << endl;
    // follows iTasc paper, figure 6
    // introduces 6 chi_f joints (7->12).
    Expression<Frame>::Ptr o1 = Constant(F_w_plane); 
    Expression<Frame>::Ptr o2 = F_w_robot1;
    Expression<Frame>::Ptr f1 = o1 * frame( KDL::vector(input(7),input(8),Constant(0.0))); 
    Expression<Frame>::Ptr f2 = o2 * frame( KDL::vector(Constant(0.0),Constant(0.0),input(12)) );

    Expression<Frame>::Ptr loop_closure = f1*frame(rot_x(input(9)))*frame(rot_y(input(10)))*frame(rot_z(input(11)))*inv(f2); 

    Expression<double>::Ptr distance_2 = input(12);
    cout << "loop closure " << endl;
    display<Frame>(cout,loop_closure);
    cout << "distance " << endl;
    display<double>(cout,distance_2);

	return 0;
}
