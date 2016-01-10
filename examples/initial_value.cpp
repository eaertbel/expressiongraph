#include <kdl/expressiontree.hpp>

/*
 * Example on how to use the initial_value function.
 * The initial_value functions takes time and another argument x.
 * When time <= 0, the argument x is cached.
 * When time >= 0, the cached argument is used instead of the current value of x.
 *
 * So, to have the correct semantics (i.e. an initial value of an expression), time has
 * to monotoneously increase and start from zero, in your program.
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
