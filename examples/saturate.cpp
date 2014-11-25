#include <kdl/expressiontree.hpp>
#include <fstream>
#include <boost/timer.hpp>

int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;

    Expression<double>::Ptr v =  saturate(input(1),0.5,1.5);

    Expression<double>::Ptr dv = v->derivativeExpression(1);    
    for (double t=0.0; t < 2.0; t+=0.05 ) {
        // remember: the order of these method calls matter:
        v->setInputValue(1,t);
        double val = v->value();
        double dval = v->derivative(1);
        // also evaluate the derivativeExpression, to compare it with the previous:
        dv->setInputValue(1,t);
        double dval2= dv->value();
        double ddval= dv->derivative(1);
        cout << t << "\t" << val << "\t" << dval << "\t" << dval2 << "\t" << ddval <<"\n"; 
    }
	return 0;
}
