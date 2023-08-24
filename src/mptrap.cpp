#include <kdl/mptrap.hpp>
#include <limits>
#include <math.h>
namespace KDL {

MPTrap::MPTrap():
    amax(1),vmax(1),spos(0),epos(0),duration(0),t1(0),t2(0),
    s(1),dT(0), d_t1_d_dpos(0),d_t2_d_dpos(0), d_duration_d_dpos(0)
{}

void MPTrap::setPlan(double _spos, double _epos, double _vmax, double _amax) {
    spos = _spos;
    epos = _epos;
    vmax = _vmax;
    amax = _amax;
}

double MPTrap::planMinDuration() {
    double dpos  = epos - spos;
    s            = dpos >= 0 ? 1  : -1;
    if (fabs(dpos) < 1E-7) {
        dpos = s*1E-7;
    }
    t1           = vmax/amax;
    double dx    = s*amax*t1*t1/2.0;
    dT           = (dpos -2*dx)*s/vmax;
    if (dT>0) {
        duration = 2*t1 + dT;
        t2       = duration - t1;
    } else {
        t1       = sqrt( dpos*s/amax);
        t2       = t1;
        duration = 2*t1;
    }
    return duration;
}

void MPTrap::adaptDuration(double new_duration) {
    double f = duration/new_duration;
    t1       = t1 / f;
    t2       = t2 / f;
    duration = new_duration;
    vmax     = vmax * f;
    amax     = amax * f*f;
}

/**
 * derivatives only to be used when this is the critical motion profile
 */
void MPTrap::compute_derivs() {
    if (dT>0) {
        d_duration_d_dpos  = s/vmax;
        d_t1_d_dpos        = 0.;
        d_t2_d_dpos        = d_duration_d_dpos;
    } else {
        d_t1_d_dpos       = -0.5/t1*s/amax;
        d_t2_d_dpos       = d_t1_d_dpos;
        d_duration_d_dpos = 2*d_t1_d_dpos;
    }

}
double MPTrap::pos(double time) {
    if (time <= 0) {
        return spos;
    } else if (time < t1) {
        return spos + s*amax*time*time/2.;
    } else if (time < t2) {
        return spos + s*amax*t1*t1/2. + s*vmax*(time-t1);
    } else if (time < duration) {
        return epos - s*amax*(duration-time)*(duration-time)/2.;
    } else {
        return epos;
    }
}

double MPTrap::d_pos_d_time(double time) {
    if (time <= 0.) {
        return 0.;
    } else if (time < t1) {
        return s*amax*time;
    } else if (time < t2) {
        return s*vmax;;
    } else if (time < duration) {
        return s*amax*(duration-time);
    } else {
        return 0.;
    }
}

double MPTrap::d_pos_d_spos(double time) {
    if (time <= 0) {
        return 1.;
    } else if (time < t1) {
        return 1.;
    } else if (time < t2) {
        return 1. - s*amax*t1*d_t1_d_dpos + s*vmax*d_t1_d_dpos;
    } else if (time < duration) {
        return s*amax*(duration-time)*d_duration_d_dpos;
    } else {
        return 0;
    }

}

double MPTrap::d_pos_d_epos(double time) {
    if (time <= 0) {
        return 0;
    } else if (time < t1) {
        return 0;
    } else if (time < t2) {
        return s*amax*t1*d_t1_d_dpos - s*vmax*d_t1_d_dpos;
    } else if (time < duration) {
        return 1. - s*amax*(duration-time)*d_duration_d_dpos;
    } else {
        return 1;
    }
}


} // end of namespace KDL
