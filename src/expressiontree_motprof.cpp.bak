/*
 * \file expressiontree_motprof.cpp
 *
 * \author Erwin Aertbelien
 *
* expressiongraph library
* 
* Copyright 2014 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
*
* Licensed under the EUPL, Version 1.1 only (the "Licence");
* You may not use this work except in compliance with the Licence.
* You may obtain a copy of the Licence at:
*
* http://ec.europa.eu/idabc/eupl 
*
* Unless required by applicable law or agreed to in writing, software 
* distributed under the Licence is distributed on an "AS IS" basis,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the Licence for the specific language governing permissions and 
* limitations under the Licence.
*/

#include <kdl/expressiontree_motprof.hpp>


namespace KDL {



double trapezoidal_mp_duration(double startpos,double endpos,double maxvel,double maxacc) {
    double s       = sign(endpos-startpos);
    double t1      = maxvel/maxacc;
    double deltax1 = s*maxacc*sqr(t1)/2.0;
    double deltaT  = (endpos-startpos-2.0*deltax1) / (s*maxvel);
    if (deltaT > 0.0) { // plan a complete profile :
        return 2*t1+deltaT;
    } else { // plan an incomplete profile :
        t1         = sqrt((endpos-startpos)/s/maxacc);
        return 2.0*t1;
    }
}

Expression<double>::Ptr trapezoidal_mp(
        Expression<double>::Ptr argument, 
        double startpos, double endpos, double maxvel, double maxacc
) {
    double s       = sign(endpos-startpos);
    double t1      = maxvel/maxacc;
    double deltax1 = s*maxacc*sqr(t1)/2.0;
    double deltaT  = (endpos-startpos-2.0*deltax1) / (s*maxvel);
    double duration,t2;
    if (deltaT > 0.0) { // plan a complete profile :
        duration   = 2*t1+deltaT;
        t2         = duration - t1;
    } else { // plan an incomplete profile :
        t1         = sqrt((endpos-startpos)/s/maxacc);
        t2         = t1;
        duration   = 2*t1;
    }
    double a3 = s*maxacc/2.0;
    double a2 = 0;
    double a1 = startpos;

    double b3 = 0;
    double b2 = a2+2*a3*t1 - 2.0*b3*t1;
    double b1 = a1+t1*(a2+a3*t1) - t1*(b2+t1*b3);

    double c3 = -s*maxacc/2.0;
    double c2 = b2+2*b3*t2 - 2.0*c3*t2;
    double c1 = b1+t2*(b2+b3*t2) - t2*(c2+t2*c3);
    Expression<double>::Ptr tmp( new Trapezoidal_Double(a1,a2,a3,b1,b2,b3,c1,c2,c3,duration,t1,t2,argument));
    return tmp;
}


Expression<double>::Ptr half_trapezoidal_mp(
        Expression<double>::Ptr argument, 
        double startpos, double endpos, double maxvel, double maxacc,
        bool starting
)
{
  Expression<double>::Ptr tmp;
  double t1,t2;
  double s = sign(endpos-startpos);
  double duration =  s*(endpos-startpos)/maxvel+maxvel/maxacc/2.0;
  if(starting) {
    t1 = 0;
    t2 = maxvel/maxacc;
    tmp = Expression<double>::Ptr( new 
         Trapezoidal_Double(startpos,
                            0,
                            0,
                            startpos,
                            -s*maxacc*t1,
                            s*maxacc/2,
                            endpos - s*maxvel*duration,
                            s*maxvel,
                            0,
                            duration,t1,t2,argument));
  } else {
    t1 = duration-maxvel/maxacc;
    t2 = duration;
    tmp = Expression<double>::Ptr( new 
         Trapezoidal_Double(startpos,
                            s*maxvel,
                            0,
                            endpos - s*maxacc*t2*t2/2,
                            s*maxacc*t2,
                            -s*maxacc/2,
                            endpos,
                            0,
                            0,
                            duration,t1,t2,argument));
  }
  
  
  return tmp;
}

Expression<double>::Ptr trapezoidal_mp_fixed_duration(
        Expression<double>::Ptr argument, 
        double startpos, double endpos, double maxvel, double maxacc,
        double newduration
) {
    double s       = sign(endpos-startpos);
    double t1      = maxvel/maxacc;
    double deltax1 = s*maxacc*sqr(t1)/2.0;
    double deltaT  = (endpos-startpos-2.0*deltax1) / (s*maxvel);
    double duration,t2;
    if (deltaT > 0.0) { // plan a complete profile :
        duration   = 2*t1+deltaT;
        t2         = duration - t1;
    } else { // plan an incomplete profile :
        t1         = sqrt((endpos-startpos)/s/maxacc);
        t2         = t1;
        duration   = 2*t1;
    }
    double a3 = s*maxacc/2.0;
    double a2 = 0;
    double a1 = startpos;

    double b3 = 0;
    double b2 = a2+2*a3*t1 - 2.0*b3*t1;
    double b1 = a1+t1*(a2+a3*t1) - t1*(b2+t1*b3);

    double c3 = -s*maxacc/2.0;
    double c2 = b2+2*b3*t2 - 2.0*c3*t2;
    double c1 = b1+t2*(b2+b3*t2) - t2*(c2+t2*c3);

    double factor = duration / newduration;
    if (factor > 1) {
        return Expression<double>::Ptr();
    }  
    Expression<double>::Ptr tmp( 
        new Trapezoidal_Double(
            a1,a2*factor,a3*factor*factor,
            b1,b2*factor,b3*factor*factor,
            c1,c2*factor,c3*factor*factor,
            newduration,t1/factor,t2/factor,argument));
    return tmp;
}







Trapezoidal_Double::Trapezoidal_Double(
    double _a1,double _a2,double _a3,
    double _b1,double _b2,double _b3,
    double _c1,double _c2,double _c3,
    double _duration, double _t1, double _t2,
    const  Expression<double>::Ptr& arg):
            UnExpr("trapezoidal",arg),
            a1(_a1),a2(_a2),a3(_a3),
            b1(_b1),b2(_b2),b3(_b3),
            c1(_c1),c2(_c2),c3(_c3),
            duration(_duration),t1(_t1),t2(_t2) {
        startpos = a1;
        endpos   = c1+duration*(c2+c3*duration);
    }


double Trapezoidal_Double::value() {
    time = argument->value();
    if (time < 0) {
        return startpos;
    } else if (time<t1) {
        return a1+time*(a2+a3*time);
    } else if (time<t2) {
        return b1+time*(b2+b3*time);
    } else if (time<=duration) {
        return c1+time*(c2+c3*time);
    } else {
        return endpos;
    }
}

double Trapezoidal_Double::derivative(int i){
    double dtimedt = argument->derivative(i);
    if (time < 0) {
        return 0;
    } else if (time<t1) {
        return (a2+2*a3*time)*dtimedt;
    } else if (time<t2) {
        return (b2+2*b3*time)*dtimedt;
    } else if (time<=duration) {
        return (c2+2*c3*time)*dtimedt;
    } else {
        return 0;
    }
}

Expression<double>::Ptr Trapezoidal_Double::derivativeExpression(int i){
    Expression<double>::Ptr expr(
        new Trapezoidal_Double( a2,2.0*a3,0.0,
                                b2,2.0*b3,0.0,
                                c2,2.0*c3,0.0,
                                duration,t1,t2, argument)
    );
    return expr;
}


 Expression<double>::Ptr Trapezoidal_Double::clone() {
    Expression<double>::Ptr expr(
        new Trapezoidal_Double( a1,a2,a3,b1,b2,b3,c1,c2,c3,duration,t1,t2, argument->clone())
    );
    return expr;
}


} // end of namespace KDL




