/*
 * expressiontree_example3.cpp
 *
 *  Created on: Sept., 2012
 *      Author: Erwin Aertbelien 
 *
*/
#include <expressiongraph/expressiontree.hpp>

/**
 * An example of a <double> expression in multiple variables.
 */
int main(int argc, char* argv[]) {
    using namespace KDL;
    using namespace std;
    
    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;

    // identical to the kinematic chain of taskfunctions3, but using
    // expression trees instead of kinematic chain,
    // and it exposes the intermediate elbow.
    Expression<Frame>::Ptr elbow = cached<Frame>(
        frame( rot_z(input(0)),  Constant(Vector(0,0,L1))) *
        frame( rot_x(input(1)) ) *
        frame( rot_z(input(2) ) ) *
        frame( rot_x(input(3) ), Constant(Vector(0,0,L2))));
    Expression<Frame>::Ptr wrist = cached<Frame>(
        elbow *
        frame( rot_z(input(4) ), Constant(Vector(0,0,L3))));
    Expression<Frame>::Ptr kinchain = cached<Frame>(
        wrist *
        frame( rot_x(input(5) ) ) *
        frame( rot_z(input(6) ), Constant(Vector(0,0,L4)) )
    );

    std::vector<double> joints(7);
    joints[0] = 0*M_PI*0.08;
    joints[1] = 30*M_PI/180.0;
    joints[2] = 0*M_PI/180.0;
    joints[3] = 120*M_PI/180.0;
    joints[4] = 0*M_PI/180.0;
    joints[5] = 30*M_PI/180.0;
    joints[6] = 0*M_PI/180.0;

    kinchain->setInputValues(joints);
    cout << "pose and Jacobian at the given joint position " << endl;
    display<Frame>(cout, kinchain );
    // visualize as a tree:
    ofstream of("tutorial2.dot");
    kinchain->write_dotfile(of);
    of.close();

    
	return 0;
}
