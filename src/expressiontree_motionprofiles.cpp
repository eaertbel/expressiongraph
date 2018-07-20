#include <kdl/expressiontree_motionprofiles.hpp>
#include <kdl/mptrap.hpp>
#include <limits>

namespace KDL {

const double mp_eps = 1E-14;
#if 0
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
            double s     = dpos >= 0 ? 1  : -1;
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

        double MPTrap::adaptDuration(double new_duration) {
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

#endif

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

static const int idx_progrvar=0;

static const int grp_offset  =1; ///< index where the groups start 
static const int idx_maxvel  =0; ///< index of max. velocity (rel. to group)
static const int idx_maxacc  =1; ///< index of max. acceleration (rel. to group)
static const int idx_startval=2; ///< index of starting value (rel. to group)
static const int idx_endval  =3; ///< index of end value (rel. to group)
static const int grp_size    =4; ///< 4 values to store for each output


MotionProfileTrapezoidal::MotionProfileTrapezoidal():
   MIMO("MotionProfileTrapezoidal") {
    // reasonable default values:
    inputDouble.push_back( input(1));               // idx_progrvar
    critical_output = -1;
}

void MotionProfileTrapezoidal::setProgressExpression(const Expression<double>::Ptr& s) {
    inputDouble[ idx_progrvar] = s;
}
Expression<double>::Ptr MotionProfileTrapezoidal::getProgressExpression() {
    return inputDouble[ idx_progrvar];
}

void MotionProfileTrapezoidal::addOutput( const Expression<double>::Ptr& startv, const Expression<double>::Ptr& endv, const Expression<double>::Ptr& maxvel, const Expression<double>::Ptr& maxacc) {
    inputDouble.push_back( maxvel );
    inputDouble.push_back( maxacc );
    inputDouble.push_back( startv );
    inputDouble.push_back( endv );
    mp.push_back( MPTrap() );
}
int MotionProfileTrapezoidal::nrOfOutputs() {
    return (inputDouble.size()-grp_offset)/grp_size;
}

void MotionProfileTrapezoidal::compute() {
    if (cached) return;
    if (mp.size()==0) return; // nothing to do.
    int idx=grp_offset;
    double duration=-1;
    critical_output=-1;
    for (int i=0;i<mp.size();++i) {
        mp[i].setPlan( 
               inputDouble[idx+idx_startval]->value(),
               inputDouble[idx+idx_endval]->value(),
               inputDouble[idx+idx_maxvel]->value(),
               inputDouble[idx+idx_maxacc]->value());
        double duration_i = mp[i].planMinDuration();
        if (duration_i > duration) {
            critical_output = i;
            duration        = duration_i;
        }
        idx += grp_size;
    }
    // POST: duration_i  >= 0 & mp.size()!=0  ==>  0 <= criticial_output < mp.size()
    // scale the noncritical motion profiles to the length of the longest:
    for (int i=0;i<mp.size();++i) {
        if (i!=critical_output) {
            mp[i].adaptDuration(duration);
        }
        mp[i].compute_derivs();   // or you can use the scaled version of the critical derives (exactly the same)
    }   
    progrvar_value = inputDouble[idx_progrvar]->value();
    cached=true; 
}

MIMO::Ptr MotionProfileTrapezoidal::clone() {
    MotionProfileTrapezoidal::Ptr tmp =
        boost::make_shared< MotionProfileTrapezoidal > ();
    tmp->setProgressExpression( 
            getProgressExpression()->clone() 
    );
    for (int i=0;i<nrOfOutputs();++i) {
        tmp->addOutput( 
                getStartValue(i)->clone(), 
                getEndValue(i)->clone(), 
                getMaxVelocity(i)->clone(), 
                getMaxAcceleration(i)->clone() );
    }
    return tmp;
}  

Expression<double>::Ptr MotionProfileTrapezoidal::getStartValue(int idx) {
    if ( (0 <= idx)&&( idx < nrOfOutputs() ) ) {
        return inputDouble[idx*grp_size+grp_offset+idx_startval];
    } else {
        throw std::out_of_range("MotionProfileTrapezoidal::getStartValue argument out of range");
    }
}

Expression<double>::Ptr MotionProfileTrapezoidal::getEndValue(int idx) {
    if ( (0<=idx) && (idx < nrOfOutputs() ) ) {
        return inputDouble[idx*grp_size+grp_offset+idx_endval];
    } else {
        throw std::out_of_range("MotionProfileTrapezoidal::getEndValue argument out of range");
    }
}

Expression<double>::Ptr MotionProfileTrapezoidal::getMaxVelocity(int idx) {
    if ( (0<=idx) && (idx < nrOfOutputs() ) ) {
        return inputDouble[idx*grp_size+grp_offset+idx_maxvel];
    } else {
        throw std::out_of_range("MotionProfileTrapezoidal::getMaxVelocity argument out of range");
    }
}

Expression<double>::Ptr MotionProfileTrapezoidal::getMaxAcceleration(int idx) {
    if ( (0<=idx) && (idx < nrOfOutputs() ) ) {
        return inputDouble[idx*grp_size+grp_offset+idx_maxacc];
    } else {
        throw std::out_of_range("MotionProfileTrapezoidal::getMaxAcceleration argument out of range");
    }
}




MotionProfileTrapezoidalOutput::MotionProfileTrapezoidalOutput(MIMO::Ptr m, int _outputnr):
    MIMO_Output<double>(_outputnr==-1?"Duration":"Output", m), outputnr(_outputnr) {
        MotionProfileTrapezoidal::Ptr p = 
            boost::static_pointer_cast<MotionProfileTrapezoidal>(m);
        if ( (outputnr<-1) || (outputnr>= p->nrOfOutputs() ) ) {
            throw std::out_of_range("MotionProfileTrapezoidal::non existing output requested");
        }
        idx_base= grp_offset + grp_size*outputnr;
    }

double MotionProfileTrapezoidalOutput::value() {
    MotionProfileTrapezoidal::Ptr p = 
        boost::static_pointer_cast<MotionProfileTrapezoidal>(mimo);
    p->compute();
    if (outputnr!=-1) {
        // output profile: 
        return p->mp[outputnr].pos( p->progrvar_value) ; 
    } else {
        // duration
        return p->mp[p->critical_output].duration;
    }
}
double MotionProfileTrapezoidalOutput::derivative(int i) {
    MotionProfileTrapezoidal::Ptr p = 
        boost::static_pointer_cast<MotionProfileTrapezoidal>(mimo);
    p->compute();
    if (outputnr!=-1) {
    // value is already called on all inputDouble's in compute
        double t= p->progrvar_value;
        return  
            p->mp[outputnr].d_pos_d_spos(t)  *  p->inputDouble[idx_base+idx_startval]->derivative(i) +
            p->mp[outputnr].d_pos_d_epos(t)  *  p->inputDouble[idx_base+idx_endval]->derivative(i) +
            p->mp[outputnr].d_pos_d_time(t)  *  p->inputDouble[idx_progrvar]->derivative(i);
    } else {
        return
            -p->mp[outputnr].d_duration_d_dpos  *  p->inputDouble[idx_base+idx_startval]->derivative(i) +
            p->mp[outputnr].d_duration_d_dpos  *  p->inputDouble[idx_base+idx_endval]->derivative(i);
    }
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
