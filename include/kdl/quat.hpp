/**
 * quat.hpp
 * expressiongraph library
 *
 * Basic encapsulation of a quaternion
 *
 * Copyright 2019 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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

#ifndef KDL_QUAT_HPP
#define KDL_QUAT_HPP

#include <kdl/frames.hpp>

#define QUAT_EPS 1E-14

namespace KDL {

    class Quaternion {
        public:
            Quaternion(const Vector& _vec): w(0.0),vec(_vec) {}
            Quaternion(double _w): w(_w),vec(Vector::Zero()) {}
            Quaternion(double _w,const Vector& _vec): w(_w),vec(_vec) {}
            Quaternion(double _w, double _x, double _y, double _z) : w(_w), vec(Vector(_x,_y,_z)) {}
            Quaternion() {}
            double w;
            Vector vec;

            inline static Quaternion Zero() {
                return Quaternion(0,Vector(0,0,0));
            }
            inline static Quaternion Identity() {
                return Quaternion(1,Vector(0,0,0));
            }
            void normalize();
    };
    inline Quaternion normalized( Quaternion d) {
        d.normalize();
        return d;
    }

    inline Vector normalized( Vector v) {
        v.Normalize();
        return v;
    }
    
    inline double squared_norm( const Quaternion& q) {
        return q.w*q.w +  dot(q.vec,q.vec);
    }
 
    inline double squared_norm( const Vector& q) {
        return dot(q,q);
    }

    inline double norm( const Quaternion& q) {
        return sqrt( squared_norm(q) );
    }

    inline void Quaternion::normalize() {
        double n = norm(*this);
        if (n < QUAT_EPS) {
            *this = Quaternion(1.,0.,0.,0.);
        } else {
            w/=n;
            vec = vec/n;
        }
    }


    inline Quaternion conj( const Quaternion& q) {
        return Quaternion(q.w, -q.vec);
    }
    inline Quaternion operator+(const Quaternion& q1, const Quaternion& q2) {
            return Quaternion(
                    q1.w + q2.w,
                    q1.vec + q2.vec
            );
    }

    inline Quaternion operator-(const Quaternion& q1, const Quaternion& q2) {
            return Quaternion(
                    q1.w - q2.w,
                    q1.vec - q2.vec
            );
    }

    inline Quaternion operator-( const Quaternion& q) {
        return Quaternion( -q.w, -q.vec);
    }


    inline Quaternion operator*(const Quaternion& q1, const Quaternion& q2) {
            return Quaternion(
                    q1.w*q2.w - dot( q1.vec, q2.vec),
                    q1.w*q2.vec + q2.w*q1.vec + q1.vec * q2.vec 
                );
    }

    inline Quaternion operator*(const Quaternion& q1, double s) {
            return Quaternion(
                    q1.w*s,
                    q1.vec*s 
            );
    }

    inline Quaternion operator*(double s, const Quaternion& q1) {
            return Quaternion(
                    q1.w*s,
                    q1.vec*s 
            );
    }

    inline Quaternion operator*(const Quaternion& q1, const Vector& p2) {
        return Quaternion(
            -dot(q1.vec, p2),
            q1.w*p2 + q1.vec*p2
        );
    }
    
    inline Quaternion operator*(const Vector& p1,const Quaternion& q2) {
        return Quaternion(
            -dot(p1,q2.vec),
            q2.w*p1 + p1*q2.vec
        );
    }


    inline Quaternion operator/(const Quaternion& q1, double s) {
            return Quaternion(
                    q1.w/s,
                    q1.vec/s 
            );
    }

    inline Quaternion inv( const Quaternion& q) {
        return conj(q)/squared_norm(q);
    }
 

    inline double dot(const Quaternion& q1, const Quaternion& q2) {
        return q1.w*q2.w + dot(q1.vec,q2.vec);
    }

    /**
     * uses the UNIT quaternion q to rotate a point given by vector p 
     */
    inline Vector apply(const Quaternion& q, const Vector& p) {
        Vector a = q.vec * p + q.w * p;
        Vector result = q.vec*a;
        result = result+result +p;
        return result;
    }

    /**
     * extracts the rotation axis and the angle from a UNIT quaternion.
     */
    inline void axis_angle(const Quaternion& q, Vector& axis, double& angle) {
        double n = q.vec.Norm(); 
        if (n > 0) {
            axis = q.vec/n;
            angle = 2*atan2(n, q.w);
        } else {
            angle = 0;
            axis  = Vector(1,0,0);
        }
    }

    /**
     * exponent of a Quaternion
     */
    inline Quaternion exp( const Quaternion& q ) {
        double n = q.vec.Norm();
        double c;  // cos
        double sc; // sinc
        if (n!=0) {
            c = cos(n);
            sc=sin(n)/n;
        } else {
            c =1;
            sc=1;
        }
        double e=exp(q.w);
        return Quaternion(e*c, e*q.vec*sc);
    }

    /**
     * exponent of a Quaternion with real part (w)  == 0, resulting into
     * a unit Quaternion.
     */
    inline Quaternion exp( const Vector& v) {
        double n = v.Norm();
        double c;  // cos 
        double sc; // sinc
        if (n!=0) {
            sc = sin(n)/n;
            c  = cos(n);
        } else {
            sc=1;
            c =1;
        }
        return Quaternion(c, v*sc);
    }

    /**
     * transform the angel*axis representation (also called Rodriguez parameters) to
     * Quaternion representation
     */
    inline Quaternion toQuat( const Vector& angle_and_axis ) {
        return exp(angle_and_axis/2.0);
    }

    /**
     * Constructs Quaternion from a rotation axis (not necessarily normalized) and an angle.
     * If norm(axis) is extremely small (QUAT_EPS), the axis is assumed to be (1,0,0)
     */
    inline Quaternion toQuat(const Vector& axis, double angle) {
        Vector n_axis = axis;
        double n2 = dot(axis,axis);
        if (n2 < (QUAT_EPS*QUAT_EPS) ) {
            n_axis = Vector(1,0,0);
        } else {
            n_axis = n_axis / sqrt( n2 );
        }
        double halfangle = angle/2.0;
        return Quaternion( cos(halfangle), sin(halfangle)*n_axis );
    }
 
    // q2 = exp( q )
    // q2 = exp(w)*( c + s/n*v)
    // ||q|| = exp(w) 
    // w = log ||q||
    // angle = atan2( n, q.w) 
    
    inline Quaternion log( const Quaternion& q) {
        double n = norm(q);
        double nv = q.vec.Norm();
        double eps = 1E-15; 
        Vector qvec_n;
        if (nv > eps ) {
            qvec_n = q.vec/nv;
            double angle = atan2(nv,q.w); 
            return Quaternion( log(n), qvec_n*angle);
        } else {        
            return Quaternion( log(n), Vector(0,0,0));
        }
    }

    /**
     * logarithm of a unit quaternion.
     * returns the vector part of the quaternion, the real part is equal to zero.
     */
    inline Vector logUnit( const Quaternion& q) {
        double n = norm(q);
        double nv = q.vec.Norm();
        double eps = 1E-15; 
        Vector qvec_n;
        if (nv > eps ) {
           qvec_n = q.vec/nv;
        } else {        
           qvec_n = Vector(1,0,0);
        }
        double angle = atan2(nv,q.w); 
        return qvec_n*angle;
    }


    /**
     * returns rotational axis multiplied by angle
     */
    inline Vector axis( const Quaternion& q) {
        return logUnit(q)*2.0;
    }
 


    inline Quaternion pow( const Quaternion& q, double t) {
        return exp( t * log(q) );
    }

    inline Quaternion powUnit( const Quaternion& q, double t) {
        return exp( t * logUnit(q) );
    }

    /**
     * Spherical interpolation between q0 (s==0) and q1 (s==1)
     */
    inline Quaternion slerp( const Quaternion& q0, const Quaternion& q1, double s ) {
        double c = dot(q0, q1);
        if (c < 0.0) {
            // half angle is over PI/2, flip one of the quaternions to get shortest path
            // from q0 to q1 
            return q0 * pow( -inv(q0)*q1, s);
        } else {
            return q0 * pow( inv(q0)*q1, s);
        }
    }

    /**
     * The rotational velocity vector to go from q1 to q2 in one unit of time.
     */ 
    inline Vector diffUnit( const Quaternion& q1, const Quaternion& q2 ) {
        double c = dot(q1, q2);
        if (c < 0.0) {
            return logUnit(-q2*conj(q1) )*2.0;
        } else {
            return logUnit(q2*conj(q1) )*2.0;
        }
    }
    inline Vector diffUnit( const Quaternion& q1, const Quaternion& q2, double dt ) {
        double c = dot(q1, q2);
        if (c < 0.0) {
            return logUnit(-q2*conj(q1) )*2.0/dt;
        } else {
            return logUnit(q2*conj(q1) )*2.0/dt;
        }
    }
/*
    inline Quaternion diff( const Quaternion& q1, const Quaternion& q2 ) {
        return (q2-q1);
    }
    inline Quaternion diff( const Quaternion& q1, const Quaternion& q2, double dt ) {
        return (q2-q1)/dt;
    }
*/





    inline Rotation toRot(Quaternion q) {
        Rotation R;
        q.normalize(); 
        const double tx  = 2.0*q.vec[0];
        const double ty  = 2.0*q.vec[1];
        const double tz  = 2.0*q.vec[2];
        const double twx = tx*q.w;
        const double twy = ty*q.w;
        const double twz = tz*q.w;
        const double txx = tx*q.vec[0];
        const double txy = ty*q.vec[0];
        const double txz = tz*q.vec[0];
        const double tyy = ty*q.vec[1];
        const double tyz = tz*q.vec[1];
        const double tzz = tz*q.vec[2];

        R(0,0) = 1.0-(tyy+tzz);
        R(0,1) = txy-twz;
        R(0,2) = txz+twy;
        R(1,0) = txy+twz;
        R(1,1) = 1.0-(txx+tzz);
        R(1,2) = tyz-twx;
        R(2,0) = txz-twy;
        R(2,1) = tyz+twx;
        R(2,2) = 1.0-(txx+tyy);

        return R;
    }

    inline Quaternion toQuat(const Rotation& R) {
        // This algorithm comes from  "Quaternion Calculus and Fast Animation",
        // Ken Shoemake, 1987 SIGGRAPH course notes
        Quaternion q;
        double t = R(0,0)+R(1,1)+R(2,2);
        if (t > 0) {
            t         = sqrt(t + 1.0);
            q.w       = 0.5*t;
            t         = 0.5/t;
            q.vec[0]  = (R(2,1) - R(1,2)) * t;
            q.vec[1]  = (R(0,2) - R(2,0)) * t;
            q.vec[2]  = (R(1,0) - R(0,1)) * t;
        } else {
            int i = 0;
            if (R(1,1) > R(0,0))
                i = 1;
            if (R(2,2) > R(i,i))
                i = 2;
            int j    = (i+1)%3;
            int k    = (j+1)%3;
            t        = sqrt( R(i,i)-R(j,j)-R(k,k) + 1.0 );
            q.vec[i] = 0.5* t;
            t        = 0.5/t;
            q.w      = ( R(k,j)-R(j,k) ) * t;
            q.vec[j] = ( R(j,i)+R(i,j) ) * t;
            q.vec[k] = ( R(k,i)+R(i,k) ) * t;
        }
        return q;
    }

} // namespace KDL
#endif 
