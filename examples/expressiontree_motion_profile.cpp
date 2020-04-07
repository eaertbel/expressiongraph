/**
 * OBSOLETE:
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
**/

#include <expressiongraph/expressiontree_motionprofiles.hpp>

int main(int argc, char* argv[]) {
    using namespace KDL;
    using namespace std;

    MotionProfileTrapezoidal::Ptr mp = create_motionprofile_trapezoidal();
    Expression<double>::Ptr time    = input(1);
    mp->setProgressExpression(time);

    Expression<double>::Ptr start1  = input(2);
    Expression<double>::Ptr end1    = input(3);
    Expression<double>::Ptr maxvel1 = input(4);
    Expression<double>::Ptr maxacc1 = input(5);
    mp->addOutput(start1,end1,maxvel1,maxacc1);

    Expression<double>::Ptr start2  = input(6);
    Expression<double>::Ptr end2    = input(7);
    Expression<double>::Ptr maxvel2 = input(8);
    Expression<double>::Ptr maxacc2 = input(9);
    mp->addOutput(start2,end2,maxvel2,maxacc2);

    Expression<double>::Ptr output1 = get_output_profile( mp, 0);
    Expression<double>::Ptr output2 = get_output_profile( mp, 1);
    Expression<double>::Ptr duration = get_duration( mp);

    // setInputValue for one of the outputs is enough...
    output1->setInputValue(2,0.0); 
    output1->setInputValue(3,4.0); 
    output1->setInputValue(4,1.0); 
    output1->setInputValue(5,0.25); 

    output1->setInputValue(6,3.0); 
    output1->setInputValue(7,0.0); 
    output1->setInputValue(8,0.8); 
    output1->setInputValue(9,0.8); 
    double dt = 0.01;
    for (double time = 0; time < duration->value(); time+=dt) {
        output1->setInputValue(1,time);
        cout << time << "\t" << output1->value() << "\t" << output1->derivative(1) << "\t" 
                             << output2->value() << "\t" << output2->derivative(1) << endl; 
    } 
    cerr << mp->mp[0].s << endl;
    cerr << mp->mp[1].s << endl;



    return 0;
}

