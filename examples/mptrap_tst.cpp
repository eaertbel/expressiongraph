// small program to verify functionality of MPTrap
#include <kdl/mptrap.hpp>

using namespace KDL;
using namespace std;

int main() {
    MPTrap mp;
    double vmax = 1.0;
    double amax = 0.5;
    mp.setPlan(0.0, 2.0, vmax, amax);
    double duration = mp.planMinDuration();
    cout << duration << endl;
    return 0;
}
