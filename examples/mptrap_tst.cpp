//
// small program to verify functionality of MPTrap
// 
// Typically this class is not used directly, but the expressiongraph node
// MotionProfileTrapezoidal is used.
//
#include <kdl/mptrap.hpp>

using namespace KDL;
using namespace std;

int main() {
    MPTrap mp;
    double vmax = 1.0;
    double amax = 0.5;
    mp.setPlan(0.0, 4.0, vmax, amax);
    double duration = mp.planMinDuration();
    mp.compute_derivs();
    //mp.adaptDuration(8);
    /*for (double time=0.;time < mp.duration; time+=0.01) {
        cout << time << "\t" << mp.pos(time) << "\t" << mp.d_pos_d_time(time) << "\n";    
    }*/

    cout << " ==== derivatives towards end position :" << endl;
    double dx = 1E-8;
    MPTrap mp2;
    mp2.setPlan(0.0, 4.0+dx, vmax, amax);
    mp2.planMinDuration();
    cout << "at time=1 :    " << mp.pos(1) << endl;
    cout << "at time=1+dt : " << mp2.pos(1) << endl;
    cout << "num. der. :    " << (mp2.pos(1)-mp.pos(1))/dx << endl;
    cout << "autodiff. :    " << mp.d_pos_d_epos(1) << endl;
    
    cout << "at time=3 :    " << mp.pos(3) << endl;
    cout << "at time=3+dt : " << mp2.pos(3) << endl;
    cout << "num. der. :    " << (mp2.pos(3)-mp.pos(3))/dx << endl;
    cout << "autodiff. :    " << mp.d_pos_d_epos(3) << endl;

    cout << "num. der. duration " << (mp2.duration - mp.duration)/dx << endl;
    cout << "autodiff duration " << mp.d_duration_d_dpos << endl;
    cout << "at time=5 :    " << mp.pos(5) << endl;
    cout << "at time=5+dt : " << mp2.pos(5) << endl;
    cout << "num. der. :    " << (mp2.pos(5)-mp.pos(5))/dx << endl;
    cout << "autodiff. :    " << mp.d_pos_d_epos(5) << endl;
    cout << " ==== derivatives towards start position :" << endl;
    dx = 1E-8;
    mp2.setPlan(0.0+dx, 4.0, vmax, amax);
    mp2.planMinDuration();
    cout << "at time=1 :    " << mp.pos(1) << endl;
    cout << "at time=1+dt : " << mp2.pos(1) << endl;
    cout << "num. der. :    " << (mp2.pos(1)-mp.pos(1))/dx << endl;
    cout << "autodiff. :    " << mp.d_pos_d_spos(1) << endl;
    
    cout << "at time=3 :    " << mp.pos(3) << endl;
    cout << "at time=3+dt : " << mp2.pos(3) << endl;
    cout << "num. der. :    " << (mp2.pos(3)-mp.pos(3))/dx << endl;
    cout << "autodiff. :    " << mp.d_pos_d_spos(3) << endl;

    cout << "num. der. duration " << (mp2.duration - mp.duration)/dx << endl;
    cout << "autodiff duration " << -mp.d_duration_d_dpos << endl;
    cout << "at time=5 :    " << mp.pos(5) << endl;
    cout << "at time=5+dt : " << mp2.pos(5) << endl;
    cout << "num. der. :    " << (mp2.pos(5)-mp.pos(5))/dx << endl;
    cout << "autodiff. :    " << mp.d_pos_d_spos(5) << endl;
 
    return 0;
}
