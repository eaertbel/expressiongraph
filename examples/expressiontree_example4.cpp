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

#include <expressiongraph/expressiontree.hpp>
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
    Expression<Frame>::Ptr kinchain = cached<Frame>(
        frame( rot_z(input(0)),  Constant(Vector(0,0,L1))) *
        frame( rot_x(-input(1)) ) *
        frame( rot_z(input(2) ) ) *
        frame( rot_x(input(3) ), Constant(Vector(0,0,L2))) *
        frame( rot_z(input(4) ), Constant(Vector(0,0,L3))) *
        frame( rot_x(-input(5) ) ) *
        frame( rot_z(input(6) ), Constant(Vector(0,0,L4)) )
    );

    std::cout << "Expression in prefix notation : \n";
    kinchain->debug_printtree();
    cout << "\n\n";

    // specify input values:  
    std::vector<double> joints(7);
    joints[0] = 0;
    joints[1] = M_PI/6;
    joints[2] = 0;
    joints[3] = 0;
    joints[4] = 0;
    joints[5] = 0;
    joints[6] = 0;
    kinchain->setInputValues( joints );
    
    cout << "--- ex. 1 : Pose of the end effector of the kin. chain ---" << endl;
    // remember: always call setInputValues() before value() 
    cout << "Position of end effector : " << kinchain->value() << endl;
    // remember: always call value() before derivative()
    for (int i=0;i<7;++i) {
        cout << "jac towards joint " << i << kinchain->derivative(i) << endl;
    }

    cout << "--- ex. 2 : Distance to point [0.2,0.2,0.8] ---" << endl;
    Vector refpoint(0.2,0.2,0.8);

    Expression<double>::Ptr distance = norm( origin(kinchain) - Constant(refpoint) );

    cout << "value " << distance->value() << "\n";
    for (int i=0;i<7;++i) {
        cout << "jac towards joint " << i << " : " << distance->derivative(i) << endl;
    }


    cout << "--- ex. 3 : Distance to Plane given by x-y axis and point in plane [0,0,0.8] ---" << endl;
    Frame plane( Rotation::Identity(), Vector(0,0,0.8) );  // x-y axis of frame corresponds to plane.

    Expression<double>::Ptr distance_to_plane =
        coord_z(Constant(plane.Inverse()) * origin(kinchain) );

    cout << "value " << distance_to_plane->value() << "\n";
    for (int i=0;i<7;++i) {
        cout << "jac towards joint " << i << " : " << distance_to_plane->derivative(i) << endl;
    }
    cout << "--- ex. 4 : Distance to a vertical line through [0.2,0.2,0] ---" << endl;
    Frame line( Rotation::Identity(), Vector(0.2,0.2,0.0)); // z axis of frame corresponds to line.

    Expression<Vector>::Ptr tmp = cached<Vector>( Constant(line.Inverse()) * origin(kinchain) );
    Expression<double>::Ptr distance_to_line = sqrt( coord_x(tmp)*coord_x(tmp) + coord_z(tmp)*coord_z(tmp) );
    cout << "value " << distance_to_line->value() << "\n";
    for (int i=0;i<7;++i) {
        cout << "jac towards joint " << i << " : " << distance_to_line->derivative(i) << endl;
    }
    
    cout << "--- ex. 5 : write an expressiontree to a graphviz .dot file " << endl;
    cout << " visualize with 'dot -Tpdf expressiontree_example4.dot > t.pdf'" << endl;
    std::ofstream of("expressiontree_example4.dot");
    distance_to_line->write_dotfile(of);
	return 0;
}
