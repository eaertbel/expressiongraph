#include <kdl/expressiontree_motionprofiles.hpp>
#include <kdl/mptrap.hpp>

namespace KDL {

const double mp_eps = 1E-14;
        
        double amax;   ///< maximum acceleration
        double vmax;   ///< maximum velocity
        double spos;   ///< start position
        double epos;   ///< end position
        double duration;
        double t1;
        double t2;
        double s;
        double dT;
        double d_t1_d_dpos;
        double d_t2_d_dpos;
        double d_duration_d_dpos;

        MPTrap::MPTrap():
            amax(1),vmax(1),spos(0),epos(0),duration(0),t1(0),t2(0),
            s(1),dT(0), d_t1_d_dpos,d_t2_d_dpos, d_duration_d_dpos
        {}

        void MPTrap::setPlan(double _spos, double _epos, double _vmax, double _amax) {
            spos = _spos;
            epos = _epos;
            vmax = _vmax;
            amax = _amax;
        }

        double MPTrap::planMinDuration() {
            double dpos  = epos - spos;
            double s     = dpos >= 0 ? 1  : -1;
            t1           = vmax/amax;
            dx           = s*amax*t1*t1/2.0;
            dT           = (dpos -2*dx)*s/amax;
            if (dT>0) {
                duration = 2*t1 + dT;
                t2       = duration - t1;
            } else {
                t1       = sqrt( dpos*s/amax);
                t2       = t1;
                duration = 2*t1;
            }
        }

        double MPTrap::adaptDuration(double new_duration) {
            double f = duration/new_duration;
            t1       = t1 / f;
            t2       = t2 / f;
            duration = new_duration;
            vmax     = vmax * f;
            amax     = amax * f;
        }

        /**
         * derivatives only to be used when this is the critical motion profile
         */
        void MPTrap::compute_derivs() {
            if (dT>0) {
                d_duration_d_dpos  = s/amax;
                d_t1_d_dpos        = 0.;
                d_t2_d_dpos        = d_duration_d_dpos;
            } else {
                d_t1_d_dpos       = -0.5/t1*s/amax;
                d_t2_d_dpos       = d_t1_d_dpos;
                d_duration_d_dpos = 2*d_t1_d_dpos;
            }

        }
        double MPTrap::pos(double time) {
            if (time < 0) {
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
            if (time < 0.) {
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
            if (time < 0) {
                return 1.;
            } else if (time < t1) {
                return 1.;
            } else if (time < t2) {
                return 1. - s*amax*t1*d_t1_d_dpos. + s*vmax*d_t1_d_dpos;
            } else if (time < duration) {
                return s*amax*(duration-time)*d_duration_d_dpos;
            } else {
                return 0;
            }

        }

        double MPTrap::d_pos_d_epos(double time) {
            if (time < 0) {
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



/**
 * To reduce the complexity:
 *  plan in function of
 *      dpos = epos-spos
 *      amax
 *      assume vmax==1
 *  => in function of two variables, the rest can be
 *     done be scaling and shifting.
 *  
 *  planner node:
 *     input: vmax, amax, deltapos
 *     output: t1,t2,duration
 *
 * computing the input:
 *     amax_normalized = amax / vmax / vmax;
 */

MotionProfileTrapezoidal::MotionProfileTrapezoidal() {
    // reasonable default values:
    inputDouble.push_back( Constant<double>(1.0));
    inputDouble.push_back( Constant<double>(1E10));
    inputDouble.push_back( input(1));
}

void MotionProfileTrapezoidal::setMaximumVelocity( const Expression<double>::Ptr& mvel) {
    inputDouble[ idx_maxvel ] = mvel;
}
Expression<double>::Ptr MotionProfileTrapezoidal::getMaximumVelocity() {
    return inputDouble[ idx_maxvel ];
}

void MotionProfileTrapezoidal::setMaximumAcceleration( const Expression<double>::Ptr& macc) {
    inputDouble[ idx_maxacc] = macc;
}
Expression<double>::Ptr MotionProfileTrapezoidal::getMaximumAcceleration() {
    return inputDouble[ idx_maxacc];
}

void MotionProfileTrapezoidal::setProgressExpression(const Expression<double>::Ptr& s) {
    inputDouble[ idx_progrvar] = s;
}
Expression<double>::Ptr MotionProfileTrapezoidal::getProgressExpression() {
    return inputDouble[ idx_progrvar];
}

void MotionProfileTrapezoidal::addOutput( const Expression<double>::Ptr& startv, const Expression<double>::Ptr& endv) {
    inputDouble.push_back( startv );
    inputDouble.push_back( endv );
}


void MotionProfileTrapezoidal::compute() {

}

MIMO::Ptr MotionProfileTrapezoidal::clone() {
    MotionProfileTrapezoidal::Ptr tmp =
        boost::make_shared< MotionProfileTrapezoidal > ();
    tmp->setMaximumVelocity( getMaximumVelocity()->clone());
    tmp->setMaximumAcceleration( getMaximumAcceleration()->clone());
    tmp->setProgressExpression(  getProgressExpression()->clone());
    return tmp;
}  
int MotionProfileTrapezoidal::nrOfOutputs() {
    return (inputDouble.size()-3)/2;
}
Expression<double>::Ptr MotionProfileTrapezoidal::getStartValue(int idx) {
    if ( (0 <= idx)&&( idx < nrOfOutputs() ) ) {
        return inputDouble[idx*2+3];
    } else {
        throw std::out_of_range("MotionProfileTrapezoidal::getStartValue argument out of range");
    }
}
Expression<double>::Ptr MotionProfileTrapezoidal::getEndValue(int idx) {
    if ( (0<=idx) && (idx < nrOfOutputs() ) ) {
        return inputDouble[idx*2+4];
    } else {
        throw std::out_of_range("MotionProfileTrapezoidal::getStartValue argument out of range");
    }
}


MotionProfileTrapezoidalOutput::MotionProfileTrapezoidalOutput(MIMO::Ptr m, int _outputnr):
    MIMO_Output<double>("MotionProfileTrapezoidal", m), outputnr(_outputnr) {
    }

double MotionProfileTrapezoidalOutput::value() {
    MotionProfileTrapezoidal::Ptr p = 
        boost::static_pointer_cast<MotionProfileTrapezoidal>(mimo);
    p->compute();
    return 0; 
}
double MotionProfileTrapezoidalOutput::derivative(int i) {
    MotionProfileTrapezoidal::Ptr p = 
        boost::static_pointer_cast<MotionProfileTrapezoidal>(mimo);
    p->compute();
    return 0;
}
MIMO_Output<double>::Ptr MotionProfileTrapezoidalOutput::clone() {
    MotionProfileTrapezoidalOutput::Ptr tmp(
            new MotionProfileTrapezoidalOutput( getMIMOClone(), outputnr));
    return tmp;
}

MotionProfileTrapezoidalOutput::Ptr 
create_motionprofile_trapezoidal_output( MotionProfileTrapezoidal::Ptr m, int i) {
    MotionProfileTrapezoidalOutput::Ptr tmp( new MotionProfileTrapezoidalOutput(m, i));
    return tmp;
}









}; // end namespace KDL
