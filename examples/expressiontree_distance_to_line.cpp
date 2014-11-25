/*
 * expressiontree_example.cpp
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
#include <fstream>

/*
 * Example that demonstrates the evaluation of an expression in KDL Types ( Frames, Rotation, Vector) for
 * both value() and derivative(). 
 *
 * This example also shows that ExpressionTree can be used as an alternative to compute the Jacobian for
 * a kinematic chain. 
 * @TODO: provide adaptor for Chain
 * @TODO: name collision between std::vector and KDL::vector (expressiontrees) can be confusing 
 *        for unexperienced users.
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


    Expression<Frame>::Ptr kinchain = kinematic_chain(chain, 0 );   
    Frame F_w_plane( Rotation::Identity(), Vector(0,0,0.5) );


    // specify input values:  
    std::vector<double> joints(7);
    joints[0] = 0;
    joints[1] = M_PI/6;
    joints[2] = 0;
    joints[3] = 0.3;
    joints[4] = 0;
    joints[5] = 0.5;
    joints[6] = 0;
    kinchain->setInputValues( joints );

    Frame line( Rotation::Identity(), Vector(0.2,0.2,0.0)); // z axis of frame corresponds to line.

    Expression<Vector>::Ptr tmp = cached<Vector>( Constant(line.Inverse()) * origin(kinchain) );
    Expression<double>::Ptr distance_to_line = sqrt( coord_x(tmp)*coord_x(tmp) + coord_z(tmp)*coord_z(tmp) );
    cout << "value " << distance_to_line->value() << "\n";
    for (int i=0;i<7;++i) {
        cout << "jac towards joint " << i << " : " << distance_to_line->derivative(i) << endl;
    }
    
    cout << " visualize with 'dot -Tpdf expressiontree_distance_to_line.dot > expressiontree_distance_to_line.pdf'" << endl;
    std::ofstream of("expressiontree_distance_to_line.dot");
    distance_to_line->write_dotfile(of);
	return 0;
}
