/**
 * \file expressiontree_example6.cpp
 * \brief example of using ExpressionTree in combination with KDL::Chain.
 *
 * \Author: Sept. 2012, Erwin Aertbelien 
***************************************************************************/

#include <expressiongraph/expressiontree.hpp>
#include <fstream>
#include <boost/timer.hpp>

/*
 * An example that investigates the computation time of expression tree's.
 * takes a 30 dof kin chain and some additional expressions and evaluates
 * the computation time.
 * 
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

    // base at [0 0 1], tool at [0 0 0.3] wrt mounting plate: 
    // chain is copied when using kinematic_chain(...)
    Expression<Frame>::Ptr kinchain1 = cached<Frame>(
         Constant(Frame(Vector(0,0,1))) * kinematic_chain( chain, 0 ) * Constant(Frame(Vector(0,0,0.3)))
    );
    // base at [0 0 0], tool at [0 0 0.3] wrt mounting plate: 
    Expression<Frame>::Ptr kinchain2 = cached<Frame>(
         Constant(Frame(Vector(0,0,0))) * kinematic_chain( chain, 7 ) * Constant(Frame(Vector(0,0,0.3)))
    );
 
   
    // z-axes of the end effector of the 2 robots are perpendicular to each other: 
    Expression<double>::Ptr perpendicular = dot(unit_z( rotation( kinchain1 )), unit_z(rotation( kinchain2 ) ));
    Expression<double>::Ptr angle = acos(dot(unit_z( rotation( kinchain1 )), unit_z(rotation( kinchain2 ) )));

    // write to a .dot file for visualization:
    std::ofstream of("expressiontree_perpendicular.dot");
    std::ofstream of2("expressiontree_angle.dot");
    perpendicular->write_dotfile(of);
    angle->write_dotfile(of2);

	
	return 0;
}
