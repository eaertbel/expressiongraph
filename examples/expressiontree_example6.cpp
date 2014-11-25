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

include <kdl/expressiontree.hpp>
#include <fstream>
#include <boost/timer.hpp>

/*
 * An example that investigates the computation time of expression tree's.
 * takes a 30 dof kin chain and some additional expressions and evaluates
 * the computation time.
 * 
 * computation for this kin. chain using expression trees is 10us
 * on a 2011 Dell Latitude Laptop.
 * @TODO : cloning on cached works unexpectedly and will result into multiple instances
 *         of the same.
 * @TODO : cached is not represented well in the graphviz plots.  One instance is represented
 *         as multiple instances.
 * @TODO: we are dealing with DAG not with tree's !
 * @TODO: handle/forbid copy constructors and assignment operators for all ExpressionTree classes.
 * @TODO: adapt license to 
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;
    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;
    Chain chain;
    chain.addSegment(Segment("Segment 0", Joint("Joint 0", Joint::RotZ),Frame(Vector(0.0,0.0,L1))));
    chain.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotX),Frame(Vector(0.0,0.0,0.0))));
    chain.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotZ),Frame(Vector(0.0,0.0,0.0))));
    chain.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),Frame(Vector(0.0,0.0,L2))));
    chain.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotZ),Frame(Vector(0.0,0.0,L3))));
    chain.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotX),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 7", Joint("Joint 7", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 8", Joint("Joint 8", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 9", Joint("Joint 9", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 10", Joint("Joint 10", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 11", Joint("Joint 11", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 12", Joint("Joint 12", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 13", Joint("Joint 13", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 14", Joint("Joint 14", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 15", Joint("Joint 15", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 16", Joint("Joint 16", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 17", Joint("Joint 17", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 18", Joint("Joint 18", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 19", Joint("Joint 19", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 20", Joint("Joint 20", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 21", Joint("Joint 21", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 22", Joint("Joint 22", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 23", Joint("Joint 23", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 24", Joint("Joint 24", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 25", Joint("Joint 25", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 26", Joint("Joint 26", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 27", Joint("Joint 27", Joint::RotY),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 28", Joint("Joint 28", Joint::RotZ),Frame(Vector(0.0,0.0,L4))));
    chain.addSegment(Segment("Segment 29", Joint("Joint 29", Joint::RotY),Frame(Vector(0.0,0.0,L4))));

    Expression<Frame>::Ptr kinchain = Constant( Frame(Vector(0,0,0.1)))*kinematic_chain( chain, 0 );   
    std::cout << "Expression in prefix notation : \n";
    kinchain->debug_printtree();
    cout << "\n\n";
     
    Frame line( Rotation::Identity(), Vector(0.2,0.2,0.0)); // z axis of frame corresponds to line.

    Expression<Vector>::Ptr tmp = cached<Vector>( Constant(line.Inverse()) * origin(kinchain) );
    Expression<double>::Ptr distance_to_line = coord_x(tmp)*coord_x(tmp) + coord_z(tmp)*coord_z(tmp);
    std::cout << "numer of derivatives " << distance_to_line->number_of_derivatives() << std::endl;
    boost::timer timer;
    int N = 100000;
    std::vector<double> joints(distance_to_line->number_of_derivatives());
    for (int n=0;n<N;++n) {
        // specify input values:  
        joints[0] = 0;
        joints[1] = n*M_PI*0.8/N;
        joints[2] = 0;
        joints[3] = 0;
        joints[4] = 0;
        joints[5] = 0;
        joints[6] = 0;
        for (int i=7;i<30;++i) {
            joints[i]=0.03;
        }
        distance_to_line->setInputValues( joints );
        distance_to_line->value();
        for (int i=0;i<joints.size();++i) {
            distance_to_line->derivative(i);
        }
    } 
    cout << "time per evaluation in microseconds : "
              << timer.elapsed()*1000000.0/N << " us " << endl;
    std::ofstream of("expressiontree_example6.dot");
    distance_to_line->write_dotfile(of);
	
	return 0;
}
