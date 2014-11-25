
#include <kdl/expressiontree_motprof.hpp>


int main(int argc,char* argv[]) {
    using namespace KDL;
    using namespace std;

    Expression<double>::Ptr time = input(1);
  
    double startpos = 1.0;
    double endpos   = -4.0;
    double maxvel   = 2.0;
    double maxacc   = 3.0; 
    Expression<double>::Ptr s = trapezoidal_mp( time, startpos,endpos, maxvel, maxacc );   
    double duration           = trapezoidal_mp_duration( startpos,endpos, maxvel, maxacc );
    cerr << "duration = " << duration << endl;
    for (double t=0.0;t<4.0;t+=0.01) {
        s->setInputValue(1,t);
        double sval = s->value();
        double sder = s->derivative(1);
        cout << t << "\t" << sval << "\t" << sder << "\n";
    }
    
    return 0;
}
