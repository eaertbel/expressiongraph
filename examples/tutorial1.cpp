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

    // define expression:
    Expression<double>::Ptr expr = Constant(2.0)*cos(input(0)) + Constant(3.0)*sin(input(1));
    std::vector<double> inp(2);
    inp[0] = 1;
    inp[1] = 2;

    // compute value and derivatives:
    expr->setInputValues(inp);
    cout << "Expression Tree in prefix notation : \n";
    expr->print(std::cout);
    cout << "\n\nValue                " << expr->value()       << endl;
    cout << "Derivative towards var. 0 " << expr->derivative(0) << endl;
    cout << "Derivative towards var. 1 " << expr->derivative(1) << endl;

    // manually compute the derivative:
    double deriv1 = -2*sin(inp[0]);
    double deriv2 = 3*cos(inp[1]);
    cout << "Symbolically computed derivative towards var. 0 " << deriv1 << endl;
    cout << "Symbolically computed derivative towards var. 1 " << deriv2 << endl;

    // visualize as a tree:
    std::ofstream of("tutorial1.dot");
    expr->write_dotfile(of);
    of.close();
    
	return 0;
}
