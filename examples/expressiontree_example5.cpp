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
 * @TODO: kinfam_io.hpp does not work for input.
 * @TODO: split up into .hpp and .cpp for all classes where this is possible.
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
        frame( rot_y(input(6) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(7) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(8) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(9) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(10) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(11) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(12) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(13) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(14) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(15) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(16) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(17) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(18) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(19) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(20) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(21) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(22) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(23) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(24) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(25) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(26) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(27) ), Constant(Vector(0,0,L4)) )*
        frame( rot_z(input(28) ), Constant(Vector(0,0,L4)) )*
        frame( rot_y(input(29) ), Constant(Vector(0,0,L4)) )
    );

    std::cout << "Expression in prefix notation : \n";
    kinchain->debug_printtree();
    cout << "\n\n";
     
    Frame line( Rotation::Identity(), Vector(0.2,0.2,0.0)); // z axis of frame corresponds to line.

    Expression<Vector>::Ptr tmp = cached<Vector>( Constant(line.Inverse()) * origin(kinchain) );
    Expression<double>::Ptr distance_to_line = sqrt( coord_x(tmp)*coord_x(tmp) + coord_z(tmp)*coord_z(tmp) );
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
    std::ofstream of("expressiontree_example5.dot");
    distance_to_line->write_dotfile(of);
	
	return 0;
}
