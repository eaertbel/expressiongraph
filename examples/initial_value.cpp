#include <kdl/expressiontree.hpp>

/*
 * Example on how to use VariableType
 */
int main(int argc, char* argv[]) {
	using namespace KDL;

    Expression<double>::Ptr time = input(0);
    
    Expression<double>::Ptr x = cached<double>( cos(time) );

    Expression<double>::Ptr y = x-initial_value<double>(time,x);

    for (double t=0;t < 2.0; t+= 0.1) {
        x->setInputValue(0,t);
        y->setInputValue(0,t);
        double val = y->value();
        double valx = x->value();
        double dval = y->derivative(0);
        std::cout << t << "\t" << val << "\t" << valx << "\t" << dval << std::endl; 
    }
    
	return 0;
}
