/*
 * expressiontree_example3.cpp
 *
 *  Created on: Sept., 2012
 *      Author: Erwin Aertbelien 
 *
***************************************************************************/

#include <kdl/expressiontree.hpp>

/**
 * This tutorial is not finished.  Correctness is not verified !
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

    
    // supose there is a camera mounted at the mounting plate of the end effector
    Expression<Frame>::Ptr camera = kinchain*Constant(Frame::Identity());
    // with its main viewing axis along Z of the mounting plate:
    // the X-axis of the mounting plate corresponds to a horizontal line in the image.

    // horizon should be horizontal: 
    Expression<double>::Ptr e_1 = dot(unit_x(rotation(camera)), Constant(Vector(0,0,1)));
   
    // pointing towards a given point p, distance point p to center-axis camera:
    Expression<Vector>::Ptr p = Constant(Vector(0.9,0.9,0)); 
    Expression<Vector>::Ptr d = cached<Vector>( unit_z(rotation(camera))); 
    Expression<double>::Ptr e_2 = norm(d*(d*(p-origin(camera))) ); 

   
    // camera image:
    Expression<Vector>::Ptr p_cam = cached<Vector>(inv(camera)*p);
    Expression<double>::Ptr x_n   = coord_x(p_cam)/coord_z(p_cam); 
    Expression<double>::Ptr y_n   = coord_y(p_cam)/coord_y(p_cam); 

    // camera parameters:
    double W  = 640;
    double H  = 480;
    double alpha_x = 20/180.0*M_PI; // openingsangle for x
    // camera model (http://en.wikipedia.org/wiki/Camera_resectioning):
    Expression<double>::Ptr fc1 =   Constant( W / tan(alpha_x) );
    Expression<double>::Ptr fc2 =   fc1; 
    Expression<double>::Ptr alpha = Constant( 0.0 ); 
    Expression<double>::Ptr u0 =    Constant( W/2 ); 
    Expression<double>::Ptr v0 =    Constant( H/2 ); 

    Expression<double>::Ptr u = fc1*x_n + alpha*fc1*y_n + u0;
    Expression<double>::Ptr v =           fc2*y_n + v0;

    cout << "point in image, x-coordinate " << u->value() << endl; 
    cout << "point in image, y-coordinate " << v->value() << endl; 

    cout << "pose and Jacobian at the given joint position " << endl;
    display<Frame>(cout, kinchain );
    // visualize as a tree:
    ofstream of("tutorial3.dot");
    kinchain->write_dotfile(of);
    of.close();

    
	return 0;
}
