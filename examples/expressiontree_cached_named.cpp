/**
 * \file expressiontree_named_cached.cpp
 * \brief example of Expression<> classes, and more specifically cached and named nodes.
 *
 * The cached node remembers the evaluation of a node, such that it can be reused later on.
 * The cached value remains available until setInputValue(s) is called.
 * 
 * One can give a name to a cached node. In an expression containing this cached node,
 * one can recover the subexpression using a subexpression<Type>(complete_expr, name) call.
 *
 * \Author: Dec 2012, Erwin Aertbelien 
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

int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;

    // the following example can be useful to understand the source code:
    cout << type_comparison<Frame,double>::return_if_equal(Frame::Identity()) << endl;
    cout << type_comparison<double,double>::return_if_equal(1.0) << endl;

    // and now for the real example:
    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;
    Expression<Frame>::Ptr robot = cached<Frame>(
        "ee",
        cached<Frame>("wrist",
            cached<Frame>("elbow",
                cached<Frame>("shoulder",
                    frame( rot_z(input(0)),  Constant(Vector(0,0,L1))) *
                    frame( rot_x(input(1)) ) *
                    frame( rot_z(input(2) ) )
                ) *
                frame( rot_x(input(3) ), Constant(Vector(0,0,L2)))
            ) *
            frame( rot_z(input(4) ), Constant(Vector(0,0,L3)))
        ) *
        frame( rot_x(input(5) ) ) *
        frame( rot_z(input(6) ), Constant(Vector(0,0,L4)) )
    );

    // write to a .dot file for visualization:
    std::ofstream of("expressiontree_named_cached1.dot");
    robot->write_dotfile(of);
    Expression<Frame>::Ptr elbow = robot->subExpression_Frame("elbow");
    if (elbow) {
        std::ofstream of("expressiontree_named_cached2.dot");
        elbow->write_dotfile(of);
    } else {
        cerr << "did not find a subexpression with the name 'elbow' " << endl;
    }
    Expression<Frame>::Ptr wrist = robot->subExpression_Frame("wrist");
    if (wrist) {
        std::ofstream of("expressiontree_named_cached3.dot");
        wrist->write_dotfile(of);
    } else {
        cerr << "did not find a subexpression with the name 'wrist' " << endl;
    }

	return 0;
}
