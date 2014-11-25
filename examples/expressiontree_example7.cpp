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

/**
 * numerically comparing to ways to compute a Jacobian of a Chain.
 *  (combining expressiontree_example5.cpp and expressiontree_example6.cpp)
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

  

    std::cout << "\nExpression 1 in prefix notation : \n";
    kinchain1->debug_printtree();
    std::cout << "\nExpression 2 in prefix notation : \n";
    kinchain2->debug_printtree();
    cout << "\n\n";
     
    Frame line( Rotation::Identity(), Vector(0.2,0.2,0.0)); // z axis of frame corresponds to line.

    Expression<Vector>::Ptr tmp1 = cached<Vector>( Constant(line.Inverse()) * origin(kinchain1) );
    Expression<Vector>::Ptr tmp2 = cached<Vector>( Constant(line.Inverse()) * origin(kinchain2) );
    Expression<double>::Ptr distance_to_line1 = sqrt( coord_x(tmp1)*coord_x(tmp1) + coord_z(tmp1)*coord_z(tmp1) );
    Expression<double>::Ptr distance_to_line2 = sqrt( coord_x(tmp2)*coord_x(tmp2) + coord_z(tmp2)*coord_z(tmp2) );
    std::cout << "numer of derivatives 1 " << distance_to_line1->number_of_derivatives() << std::endl;
    std::cout << "numer of derivatives 2 " << distance_to_line2->number_of_derivatives() << std::endl;
    boost::timer timer;
    int N = 10;
    std::vector<double> joints(distance_to_line1->number_of_derivatives());
        // specify input values:  
        joints[0] = 0*M_PI*0.8/N;
        joints[1] = 1*M_PI*0.8/N;
        joints[2] = 2*M_PI*0.8/N;
        joints[3] = 3*M_PI*0.8/N;
        joints[4] = 4*M_PI*0.8/N;
        joints[5] = 5*M_PI*0.8/N;
        joints[6] = 6*M_PI*0.8/N;
        distance_to_line1->setInputValues( joints );
        distance_to_line2->setInputValues( joints );
        cout << "Value : \n";\
        cout << "chain     \t,\t pure expr \n";
        cout << distance_to_line1->value() << "\t,\t" << distance_to_line2->value() << "\n";
        cout << "Derivatives : \n";
        cout << "chain     \t,\t num.deriv.\t,\tpure expr \t,\tnum. deriv \n";
        for (int i=0;i<joints.size();++i) {
            cout << setw(10) << distance_to_line1->derivative(i) << "\t,\t" 
                 << setw(10) << numerical_derivative<double>(distance_to_line1,i,joints[i],1E-10)<< "\t,\t"
                 << setw(10) << distance_to_line2->derivative(i) <<"\t,\t" 
                 << setw(10) << numerical_derivative<double>(distance_to_line2,i,joints[i],1E-10)<< "\n";
        }
	return 0;
}
