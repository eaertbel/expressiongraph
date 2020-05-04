/**
 * expressiontree_quaternion.hpp
 * expressiongraph library
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
 
     **********************************************************************************
     *   QUATERNION OPERATORS/FUNCTIONS AND THEIR NOTATION:
     *
     *   arguments of the operators follow the following notation:
     *      q : quaternion (w,vec) with real part w and a vector part vec.
     *      u : pure quaternion (i.e. real part w==0), represented by Vector or vector
     *      qn : unit quaternion (i.e. modulo (norm) == 1 ), to indicate that only
     *           unit quaternions are allowed as argument.
     *      s  : scalar (double value)
     *      v  : vector
     *      vn : normalized vector
     *      R  : a rotation matrix
     *
     *   q1 + q1        : quaternion sum
     *   q1 - q1        : quaternion difference 
     *   -q             : unary negation, i.e. (-w, -vec)
     *   q*s , s*q      : product scalar-quaternion
     *   q*u , u*q      : product vector-quaternion,  vector u is interpreted as pure 
     *                    quaternion (0, u) 
     *   q1*q2          : quaternion product (w,vec)
     *                    with real part   w   : q1.w*q2.w - dot(q1.vec,q2.vec)  
     *                         vector part vec : q1.w*q2.vec + q2.w*q1.vec + cross(q1.vec, q2.vec)
     *   dot(q1,q2)     : quaternion dot product, i.e. q1.w*q2.w + dot(q1.vec,q2.vec)
     *   q/s            : division of a quaternion by a scalar (w/s, vec/s)
     *   conj(q)        : quaternion conjugate, i.e. (w,-vec)
     *   w(q)           : real part of quaternion == (q+conj(q))/2
     *   vec(q)         : vector part of quaternion == (q-conj(q))/2
     *   squared_norm(q) : conj(q)*q
     *   norm(q)        : sqrt( conj(q)*q )
     *   inv(q)         : quaternion inverse, i.e. inv(q)*q == 1,  inv(q) == conj(q)/squared_norm(q)
     *   normalized(q)  : returns a unit quaternion, q/norm(q), 
     *                    for norm(q)==0, returns quaternion with w=1 and vec=0,0,0.
     *   normalized(v)  : returns a unit vector, v/norm(v), for norm(v)==0, returns vec=0,0,0 )
     *   apply(qn,v)    : computes efficiently the rotation of vector v by unit quaternion qn.
     *                    assumes that the unit quaternion qn has norm==1.
     *   exp(q)         : quaternion exponential of q, resulting in a not necessarily normalized
     *                    quaternion.  The exponential of (w,vec=0,0,0) is equal 
     *                    to exp(w)*Quaternion(1,(0,0,0))
     *   exp(u)         : quaternion exponential of a pure quaternion represented by a vector u
     *                    resulting in a normalized quaternion.  The exponential of u=0,0,0
     *                    is equal to a quaternion with w=1 vec=(0,0,0).
     *   log(q)         : quaternion logarithm of q. -inf if norm(q)==0 (similar to log(0) ).
     *                    if norm(vec(q)) very small, resuts in log(norm(q), Vector(0,0,0)).
     *   logUnit(q)     : quaternion logarighm of a unit quaternion, resulting in a pure quaternion
     *                    represented by a vector.
     *                    if norm(vec(q)) very small, resuts in Vector(0,0,0).
     *   pow(q,s)       : quaternion q to the power s.  s is a scalar. equal to exp(s*log(q)).
     *   powUnit(qn,s)  : unit quaternion q to the power s.  Returns a normalized quaternion.               
     *                    equal to exp( s* logUnit(q) )
     *   toQuat(v)      : interprets the vector as a normal vector indicating the rotation axis,
     *                    multiplied by a full (NOT half!) angle (i.e. "axis-angle" representation) and returns
     *                    the equivalent quaternion.
     *   toQuat(v, s): returns a quaternion corresponding to a rotation around v (non necessarily
     *                     normalized) with the angle given by s.
     *   axisAngle(qn)  : returns the axis-angle representation for a given unit quaternion, i.e.
     *                    the rotation axis and the full angle (NOT half!) around that rotation axis.
     *   slerp(q1,q2,s) : spherical linear interpolation between q1 (corresponding to s==0) and q2 
     *                    (corresponding to s==1).
     *   slerpUnit(qn1,qn2,s) : spherical linear interpolation between UNIT quaternion q1 (corresponding to s==0) and 
     *                      UNIT quaternion q2 (corresponding to s==1).
     *   diff(q1,q2)     : q2-q1
     *   diffUnit(q1,q2) : difference over the sphere, i.e. rotational velocity  that you need to go from q1 to q2 in 
     *                     one unit of time, i.e. angle-axis to go from q1 to q2
     *   toQuat(R)      : returns the quaternion corresponding to the given rotation matrix
     *   toRot(qn)       : returns the rotation matrix corresponding to the given unit quaternion.
     *   quaternion(s,v),
     *   quaternion(v,s) : returns quaternion with real part s and vector part v
     *
     ****************************************************************************************************
**/



#ifndef KDL_EXPRESSIONTREE_QUATERNION_HPP
#define KDL_EXPRESSIONTREE_QUATERNION_HPP
#include <expressiongraph/quat.hpp>
#include <expressiongraph/quat_io.hpp>
#include <expressiongraph/expressiontree_expressions.hpp>
#include <expressiongraph/expressiontree_double.hpp>
#include <expressiongraph/expressiontree_vector.hpp>

/**
 *
 * Quaternion class for both value and derivative
 *
 */

namespace KDL {

    class Quaternion_conj:
        public UnaryExpression<Quaternion,Quaternion>
    {   
    public:
        typedef UnaryExpression<Quaternion,Quaternion> UnExpr;
 
        Quaternion_conj() {}

        Quaternion_conj(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("conj",arg) {}
        
        virtual Quaternion value() {
            Quaternion a = argument->value();   
            return Quaternion(a.w, -a.vec);
        }

        virtual Quaternion derivative(int i) {
            Quaternion da = argument->derivative(i);
            return Quaternion(da.w, -da.vec);
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
            return boost::make_shared<Quaternion_conj>(argument->clone());
        } 
    };

    class Quaternion_w:
        public UnaryExpression<double,Quaternion>
    {   
    public:
        typedef UnaryExpression<double,Quaternion> UnExpr;
 
        Quaternion_w() {}

        Quaternion_w(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("quaternion_w",arg) {}
        
        virtual double value() {
            return argument->value().w;
        }

        virtual double derivative(int i) {
            return argument->derivative(i).w;
        }

        virtual Expression<double>::Ptr derivativeExpression(int i);

        virtual Expression<double>::Ptr clone() {
            return boost::make_shared<Quaternion_w>(argument->clone());
        } 
    };

    class Quaternion_vec:
        public UnaryExpression<Vector,Quaternion>
    {   
    public:
        typedef UnaryExpression<Vector,Quaternion> UnExpr;
 
        Quaternion_vec() {}

        Quaternion_vec(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("quaternion_vec",arg) {}
        
        virtual Vector value() {
            return argument->value().vec;
        }

        virtual Vector derivative(int i) {
            return argument->derivative(i).vec;
        }

        virtual Expression<Vector>::Ptr derivativeExpression(int i);

        virtual Expression<Vector>::Ptr clone() {
            return boost::make_shared<Quaternion_vec>(argument->clone());
        } 
    };

    class Quaternion_squared_norm:
        public UnaryExpression<double,Quaternion>
    {   
    public:
        Quaternion q;
        typedef UnaryExpression<double,Quaternion> UnExpr;
 
        Quaternion_squared_norm() {}

        Quaternion_squared_norm(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("squared_norm",arg) {}
        
        virtual double value() {
            q = argument->value();
            return q.w*q.w + dot(q.vec,q.vec);
        }

        virtual double derivative(int i) {
            Quaternion dq;
            dq = argument->derivative(i);
            return 2*q.w*dq.w + 2*dot(q.vec, dq.vec); 
        }

        virtual Expression<double>::Ptr derivativeExpression(int i);

        virtual Expression<double>::Ptr clone() {
            return boost::make_shared<Quaternion_squared_norm>(argument->clone());
        } 
    };

    class Quaternion_norm:
        public UnaryExpression<double,Quaternion>
    {   
    public:
        Quaternion q;
        double     n;
        typedef UnaryExpression<double,Quaternion> UnExpr;
 
        Quaternion_norm() {}

        Quaternion_norm(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("norm",arg) {}
        
        virtual double value() {
            q = argument->value();
            n = sqrt(q.w*q.w + dot(q.vec,q.vec));  
            return n;
        }

        virtual double derivative(int i) {
            Quaternion dq = argument->derivative(i);
            return (q.w*dq.w + dot(q.vec, dq.vec) ) / n ; 
        }

        virtual Expression<double>::Ptr derivativeExpression(int i);

        virtual Expression<double>::Ptr clone() {
            return boost::make_shared<Quaternion_norm>(argument->clone());
        } 
    };

    class Sum_QuaternionQuaternion:
        public BinaryExpression<Quaternion,Quaternion, Quaternion>
    {
        public:
            Quaternion arg1value;
            Quaternion arg2value;
            typedef BinaryExpression<Quaternion,Quaternion, Quaternion> BinExpr;
            Sum_QuaternionQuaternion(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("plus",arg1,arg2)
				{}
           	virtual Quaternion value() {
		        return argument1->value() + argument2->value();
	        }

	        virtual Quaternion derivative(int i){
                return argument1->derivative(i) + argument2->derivative(i);
            } 

            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Sum_QuaternionQuaternion>( argument1->clone(), argument2->clone());
            }
    };

    class Minus_QuaternionQuaternion:
        public BinaryExpression<Quaternion,Quaternion, Quaternion>
    {
        public:
            Quaternion arg1value;
            Quaternion arg2value;
            typedef BinaryExpression<Quaternion,Quaternion, Quaternion> BinExpr;
            Minus_QuaternionQuaternion(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("minus",arg1,arg2)
				{}
           	virtual Quaternion value() {
		        return argument1->value() - argument2->value();
	        }

	        virtual Quaternion derivative(int i){
                return argument1->derivative(i) - argument2->derivative(i);
            } 

            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Minus_QuaternionQuaternion>( argument1->clone(), argument2->clone());
            }
    };

    class Quaternion_Minus:
        public UnaryExpression<Quaternion,Quaternion>
    {   
    public:
        typedef UnaryExpression<Quaternion,Quaternion> UnExpr;
 
        Quaternion_Minus() {}

        Quaternion_Minus(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("minus",arg) {}
        
        virtual Quaternion value() {
            Quaternion a = argument->value();   
            return Quaternion(-a.w, -a.vec);
        }

        virtual Quaternion derivative(int i) {
            Quaternion da = argument->derivative(i);
            return Quaternion(-da.w, -da.vec);
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
            return boost::make_shared<Quaternion_Minus>(argument->clone());
        } 
    };

    class Mult_QuaternionQuaternion:
        public BinaryExpression<Quaternion,Quaternion, Quaternion>
    {
        public:
            Quaternion arg1value;
            Quaternion arg2value;
            typedef BinaryExpression<Quaternion,Quaternion, Quaternion> BinExpr;
            Mult_QuaternionQuaternion(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("mult",arg1,arg2)
				{}
           	virtual Quaternion value() {
		        arg1value = argument1->value();
		        arg2value = argument2->value();
                return arg1value * arg2value;
	        }

	        virtual Quaternion derivative(int i){
		        Quaternion da1 = argument1->derivative(i);
		        Quaternion da2 = argument2->derivative(i);
                return da1*arg2value + arg1value*da2;
            } 
            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Mult_QuaternionQuaternion>( argument1->clone(), argument2->clone());
            }
    };

    class Mult_QuaternionVector:
        public BinaryExpression<Quaternion,Quaternion, Vector>
    {
        public:
            Quaternion arg1value;
            Vector arg2value;
            typedef BinaryExpression<Quaternion,Quaternion, Vector> BinExpr;
            Mult_QuaternionVector(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("mult",arg1,arg2)
				{}
           	virtual Quaternion value() {
		        arg1value = argument1->value();
		        arg2value = argument2->value();
                return arg1value * arg2value;
	        }

	        virtual Quaternion derivative(int i){
		        Quaternion da1 = argument1->derivative(i);
		        Vector da2 = argument2->derivative(i);
                return da1*arg2value + arg1value*da2;
            } 
            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Mult_QuaternionVector>( argument1->clone(), argument2->clone());
            }
    };

    class Mult_VectorQuaternion:
        public BinaryExpression<Quaternion,Vector,Quaternion>
    {
        public:
            Vector arg1value;
            Quaternion arg2value;
            typedef BinaryExpression<Quaternion,Vector,Quaternion> BinExpr;
            Mult_VectorQuaternion(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("mult",arg1,arg2)
				{}
           	virtual Quaternion value() {
		        arg1value = argument1->value();
		        arg2value = argument2->value();
                return arg1value * arg2value;
	        }

	        virtual Quaternion derivative(int i){
		        Vector da1 = argument1->derivative(i);
		        Quaternion da2 = argument2->derivative(i);
                return da1*arg2value + arg1value*da2;
            } 
            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Mult_VectorQuaternion>( argument1->clone(), argument2->clone());
            }
    };


    class Mult_QuaternionScalar:
        public BinaryExpression<Quaternion,Quaternion, double>
    {
        public:
            Quaternion arg1value;
            double         arg2value;
            typedef BinaryExpression<Quaternion,Quaternion, double> BinExpr;

            Mult_QuaternionScalar(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("mult",arg1,arg2)
				{}

           	virtual Quaternion value() {
		        arg1value = argument1->value();
		        arg2value = argument2->value();
                return arg1value * arg2value;
	        }

	        virtual Quaternion derivative(int i){
		        Quaternion da1 = argument1->derivative(i);
		        double          da2 = argument2->derivative(i);
                return da1*arg2value + arg1value*da2;
            } 
            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Mult_QuaternionScalar>( argument1->clone(), argument2->clone());
            }
    };

    class Div_QuaternionScalar:
        public BinaryExpression<Quaternion,Quaternion, double>
    {
        public:
            Quaternion arg1value;
            double         arg2value;
            typedef BinaryExpression<Quaternion,Quaternion, double> BinExpr;

            Div_QuaternionScalar(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("division",arg1,arg2)
				{}

           	virtual Quaternion value() {
		        arg1value = argument1->value();
		        arg2value = argument2->value();
                return arg1value / arg2value;
	        }

	        virtual Quaternion derivative(int i){
		        Quaternion da1 = argument1->derivative(i);
		        double          da2 = argument2->derivative(i);
                return da1/arg2value - arg1value*da2/arg2value/arg2value;
            } 
            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Div_QuaternionScalar>( argument1->clone(), argument2->clone());
            }
    };

    class Dot_QuaternionQuaternion:
        public BinaryExpression<double,Quaternion, Quaternion>
    {
        public:
            Quaternion     a1;
            Quaternion     a2;
            typedef BinaryExpression<double, Quaternion,Quaternion> BinExpr;

            Dot_QuaternionQuaternion(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("dot",arg1,arg2)
				{}

           	virtual double value() {
		        a1 = argument1->value();
		        a2 = argument2->value();
                return dot(a1,a2);
	        }

	        virtual double derivative(int i){
		        Quaternion da1 = argument1->derivative(i);
		        Quaternion da2 = argument2->derivative(i);
                return a1.w * da2.w + da1.w * a2.w +
                       dot( a1.vec, da2.vec) + dot(da1.vec, a2.vec);
            } 
            virtual Expression<double>::Ptr derivativeExpression(int i);
            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Dot_QuaternionQuaternion>( argument1->clone(), argument2->clone());
            }
    };

    class Normalize_Quaternion:
        public UnaryExpression<Quaternion,Quaternion>
    {   
    public:
        typedef UnaryExpression<Quaternion,Quaternion> UnExpr;
        Quaternion result;
        double n;
        bool   iszero;
        Quaternion q;
 
        Normalize_Quaternion() {}

        Normalize_Quaternion(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("normalize",arg) {}
 
        virtual Quaternion value() {
            q   = argument->value();
            double sqn = q.w*q.w + dot(q.vec,q.vec);
            iszero = sqn <= (QUAT_EPS*QUAT_EPS);
            if (iszero) {
                n = 0.0;
                return Quaternion(1.0,0.0,0.0,0.0);
            } else {
                n = sqrt(sqn); 
                result = Quaternion( q.w/n, q.vec/n);
                return result;
            } 
        }

        virtual Quaternion derivative(int i) {
            // some effort to compute it efficiently, since we are going
            // to use this function a lot...:
            if (iszero) {
                // breaks down at this point, we just do a best effort:
                return Quaternion(0,0,0,0);
            }
            Quaternion Dq = argument->derivative(i);
            double Dn     = (q.w*Dq.w + dot(q.vec, Dq.vec))/n;
            return (Dq - result*Dn)/n;
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
            return boost::make_shared<Normalize_Quaternion>(argument->clone());
        } 
    };

    class Normalize_Vector:
        public UnaryExpression<Vector,Vector>
    {   
    public:
        typedef UnaryExpression<Vector,Vector> UnExpr;
        Vector result;
        Vector q;
        double n;
        bool   iszero;
 
        Normalize_Vector() {}

        Normalize_Vector(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("normalize",arg) {}
 
        virtual Vector value() {
            q   = argument->value();
            double sqn = dot(q,q);
            iszero = sqn <= (QUAT_EPS*QUAT_EPS);
            if (iszero) {
                n = 0.0;
                return Vector(0.0,0.0,0.0);
            } else {
                n = sqrt(sqn); 
                result = q/n;
                return result;
            } 
        }

        virtual Vector derivative(int i) {
            // some effort to compute it efficiently, since we are going
            // to use this function a lot...:
            if (iszero) {
                // breaks down at this point, we just do a best effort:
                return Vector(0,0,0);
            }
            // n*2 = 2*q*Dq
            // 2*n*dn = 2*q*Dq
            // dn = q*dq/n
            // q/n = r  -> r*n = q
            // Dq = Dr*n + r*Dn
            // Dr = (Dq - r*Dn)/n
            Vector Dq = argument->derivative(i);
            double Dn = dot(q, Dq)/n;
            return (Dq - result*Dn)/n;
        }

        virtual Expression<Vector>::Ptr derivativeExpression(int i);

        virtual Expression<Vector>::Ptr clone() {
            return boost::make_shared<Normalize_Vector>(argument->clone());
        } 
    };

    class Apply_QuaternionVector:
        public BinaryExpression<Vector,Quaternion, Vector>
    {
        public:
            Quaternion     q;
            Vector         p;
            Vector         a,b;
            typedef BinaryExpression<Vector, Quaternion,Vector> BinExpr;

            Apply_QuaternionVector(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("apply",arg1,arg2)
				{}

           	virtual Vector value() {
		        q = argument1->value();
		        p = argument2->value();
                a = q.vec*p + q.w*p;
                b = q.vec*a;
                return b+b+p;
	        }

	        virtual Vector derivative(int i){
		        Quaternion dq = argument1->derivative(i);
		        Vector dp     = argument2->derivative(i);
                Vector     da = q.vec*dp+dq.vec*p + dq.w*p + q.w*dp;
                Vector     db = q.vec*da + dq.vec*a;
                return db+db+dp;
            } 
            virtual Expression<Vector>::Ptr derivativeExpression(int i);
            virtual Expression<Vector>::Ptr clone() {
                return boost::make_shared<Apply_QuaternionVector>( argument1->clone(), argument2->clone());
            }
    };

    class Quaternion_exp:
        public UnaryExpression<Quaternion,Quaternion>
    {   
    public:
        typedef UnaryExpression<Quaternion, Quaternion> UnExpr;
        Quaternion q;       
        Quaternion expq_unit;  
        double n;     
        double c;     
        double s;     
        double sc;     
        double e;     
        double e_sc;
 
        Quaternion_exp() {}

        Quaternion_exp(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("exp",arg) {}
        
        virtual Quaternion value() {
            q = argument->value();
            n = q.vec.Norm();
            if (n!=0) {
                c = cos(n);
                s = sin(n);
                sc= s/n;
            } else {
                c = 1;
                s = 0;
                sc= 1;
            }
            e   = exp(q.w);
            e_sc= e*sc;
            return Quaternion(e*c, e_sc*q.vec); 
        }

        virtual Quaternion derivative(int i) {
            Quaternion d_q  = argument->derivative(i);
            if (n==0) {
                return Quaternion(e*d_q.w, 0.0, 0.0, 0.0);
            } else {
                double d_n  = dot(q.vec,d_q.vec)/n;
                double d_sc = ( c - sc)/n*d_n;
                double d_e  = e * d_q.w;
                double d_c  = -s*d_n;   
                double d_e_sc = e*d_sc + d_e*sc;    
                return Quaternion( d_e*c + e*d_c,  e_sc*d_q.vec + d_e_sc*q.vec ); 
            }
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
           return boost::make_shared<Quaternion_exp>(argument->clone());
        } 
    };

    class Quaternion_expv:
        public UnaryExpression<Quaternion,Vector>
    {   
    public:
        typedef UnaryExpression<Quaternion,Vector> UnExpr;
        Vector q;       
        Quaternion expq;       
        Quaternion expq_unit;  
        double n;     
        double c;     
        double s;     
        double sc;     
 

 
        Quaternion_expv() {}

        Quaternion_expv(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("exp",arg) {}
        
        virtual Quaternion value() {
            q = argument->value();
            n = q.Norm();
            if (n!=0) {
                c = cos(n);
                s = sin(n);
                sc= s/n;
            } else {
                c = 1;
                s = 0;
                sc= 1;
            }
            return Quaternion(c, sc*q); 
        }

        virtual Quaternion derivative(int i) {
            Vector d_q = argument->derivative(i);
            if (n==0) {
                return Quaternion(1.0, 0.0, 0.0, 0.0);
            } else {
                double d_n  = dot(q,d_q)/n;
                double d_sc = ( c - sc)/n*d_n;
                double d_c  = -s*d_n;   
                return Quaternion( d_c,  sc*d_q + d_sc*q ); 
            }
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
           return boost::make_shared<Quaternion_expv>(argument->clone());
        } 
    };


    class Quaternion_log:
        public UnaryExpression<Quaternion,Quaternion>
    {   
    public:
        typedef UnaryExpression<Quaternion, Quaternion> UnExpr;
        Quaternion q;       
        Vector qvec_n;
        double m; // modulo
        double n; // norm(q.vec)
        double angle;
        
 
        Quaternion_log() {}

        Quaternion_log(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("log",arg) {}
        
        virtual Quaternion value() {
            // m==0 remains an error
            q = argument->value();
            m = norm(q);
            n = q.vec.Norm();
            if (n > QUAT_EPS ) {
                qvec_n = q.vec/n;
                angle = atan2(n,q.w); 
                return Quaternion( log(m), qvec_n*angle);
            } else {        
                return Quaternion( log(m), Vector(0,0,0));
            }
        }

        virtual Quaternion derivative(int i) {
            // m==0 remains an error
            Quaternion d_q = argument->derivative(i);
            double d_m      = dot(q,d_q)/m;
            if (n > QUAT_EPS) {
                double d_n      = dot(q.vec, d_q.vec)/n;
                double d_angle  = ( d_n * q.w - n * d_q.w ) / ( n*n + q.w*q.w );
                Vector d_qvec_n = d_q.vec/n - q.vec/n/n*d_n;
                return Quaternion( d_m/m , d_qvec_n*angle + qvec_n*d_angle ); 
            } else {
                return Quaternion( d_m/m , Vector(0,0,0) ); 
            }
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
           return boost::make_shared<Quaternion_log>(argument->clone());
        } 
    };

    class Quaternion_logUnit:
        public UnaryExpression<Vector,Quaternion>
    {   
    public:
        typedef UnaryExpression<Vector, Quaternion> UnExpr;
        Quaternion q;       
        Vector qvec_n;
        double n; // norm(q.vec)
        double angle;
      
        Quaternion_logUnit() {}

        Quaternion_logUnit(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("logUnit",arg) {}
        
        virtual Vector value() {
            q = argument->value();
            n = q.vec.Norm();
            if (n > QUAT_EPS ) {
                qvec_n = q.vec/n;
                angle = atan2(n,q.w); 
                return qvec_n*angle;
            } else {        
                return Vector(0,0,0);
            }
        }

        virtual Vector derivative(int i) {
            Quaternion d_q = argument->derivative(i);
            if (n>QUAT_EPS) {
                double d_n      = dot(q.vec, d_q.vec)/n;
                double d_angle  = ( d_n * q.w - n * d_q.w );
                Vector d_qvec_n = d_q.vec/n - q.vec/n/n*d_n;
                return  d_qvec_n*angle + qvec_n*d_angle ; 
            } else {
                return Vector(0,0,0);
            }
        }

        virtual Expression<Vector>::Ptr derivativeExpression(int i);

        virtual Expression<Vector>::Ptr clone() {
           return boost::make_shared<Quaternion_logUnit>(argument->clone());
        } 
    };



    class QuaternionFromRotation:
        public UnaryExpression<Quaternion,Rotation>
    {   
    public:
        typedef UnaryExpression<Quaternion, Rotation> UnExpr;
        Quaternion q;       
 
        QuaternionFromRotation() {}

        QuaternionFromRotation(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("quaternionFromRotation",arg) {}
        
        virtual Quaternion value() {
            q = toQuat(argument->value());
            return q;
        }

        virtual Quaternion derivative(int i) {
            return Quaternion(0.0, argument->derivative(i)/2.0)*q;
        }

        virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

        virtual Expression<Quaternion>::Ptr clone() {
           return boost::make_shared<QuaternionFromRotation>(argument->clone());
        } 
    };

    class RotationFromQuaternion:
        public UnaryExpression<Rotation,Quaternion>
    {   
    public:
        typedef UnaryExpression<Rotation, Quaternion> UnExpr;
        Quaternion q;       
 
        RotationFromQuaternion() {}

        RotationFromQuaternion(const UnExpr::ArgumentExpr::Ptr& arg):
            UnExpr("rotationFromQuaternion",arg) {}
        
        virtual Rotation value() {
            //q = argument->value(); 
            // store 2*inverse(q):
            //double N = q.w*q.w + dot(q.vec,q.vec);
            //q.w   = 2.0*q.w/N;
            //q.vec = -2.0*q.vec/N;
            q = argument->value();
            return toRot( q );
        }

        virtual Vector derivative(int i) {
            return (argument->derivative(i)*conj(q)).vec*2.0;
        }

        virtual Expression<Vector>::Ptr derivativeExpression(int i);

        virtual Expression<Rotation>::Ptr clone() {
            return boost::make_shared<RotationFromQuaternion>(argument->clone());
        } 
    };

    class Quaternion_Scalar_Vector:
        public BinaryExpression<Quaternion,double, Vector>
    {
        public:
            typedef BinaryExpression<Quaternion,double, Vector> BinExpr;
            Quaternion_Scalar_Vector(
            			const  BinExpr::Argument1Expr::Ptr& arg1,
			            const  BinExpr::Argument2Expr::Ptr& arg2
                ): BinExpr("quaternion",arg1,arg2)
				{
                }
           	virtual Quaternion value() {
                return Quaternion(argument1->value(), argument2->value());
	        }

	        virtual Quaternion derivative(int i){
                return Quaternion(argument1->derivative(i), argument2->derivative(i));
            } 

            virtual Expression<Quaternion>::Ptr derivativeExpression(int i);

            virtual BinExpr::Ptr clone() {
                return boost::make_shared<Quaternion_Scalar_Vector>( argument1->clone(), argument2->clone());
            }
    };


    

    /**********************************************************************************
     *   OPERATORS
     **********************************************************************************/



    inline Expression<Quaternion>::Ptr conj( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_conj>(a);
    }

    inline Expression<double>::Ptr w( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_w>(a);
    }

    inline Expression<Vector>::Ptr vec( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_vec>(a);
    }

    inline Expression<double>::Ptr squared_norm( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_squared_norm>(a);
    }

    inline Expression<double>::Ptr norm( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_norm>(a);
    }

    inline Expression<Quaternion>::Ptr operator+( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Sum_QuaternionQuaternion>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr operator-( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Minus_QuaternionQuaternion>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr operator-( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_Minus>(a);
    }
 
    inline Expression<Quaternion>::Ptr operator* ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<double>::Ptr a2 ) {
        return boost::make_shared<Mult_QuaternionScalar>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr operator* ( 
            Expression<double>::Ptr a2 ,
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Mult_QuaternionScalar>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr operator* ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Mult_QuaternionQuaternion>( a1, a2 );
    }
    inline Expression<Quaternion>::Ptr operator* ( 
            Expression<Vector>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Mult_VectorQuaternion>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr operator* ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Vector>::Ptr a2 ) {
        return boost::make_shared<Mult_QuaternionVector>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr operator/ ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<double>::Ptr a2 ) {
        return boost::make_shared<Div_QuaternionScalar>( a1, a2 );
    }
    
    inline Expression<Quaternion>::Ptr inv( Expression<Quaternion>::Ptr q) {
        return conj(q) / squared_norm(q);
    }

    inline Expression<double>::Ptr dot ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Dot_QuaternionQuaternion>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr normalized ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Normalize_Quaternion>( a1);
    }

    inline Expression<Vector>::Ptr normalized ( 
            Expression<Vector>::Ptr a1) {
        return boost::make_shared<Normalize_Vector>( a1);
    }

    inline Expression<Vector>::Ptr apply ( 
            Expression<Quaternion>::Ptr q, Expression<Vector>::Ptr v) {
        return boost::make_shared<Apply_QuaternionVector>(q,v);
    }

    inline Expression<Quaternion>::Ptr exp ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Quaternion_exp>( a1);
    }

    inline Expression<Quaternion>::Ptr exp ( 
            Expression<Vector>::Ptr a1) {
        return boost::make_shared<Quaternion_expv>( a1);
    }


    inline Expression<Quaternion>::Ptr log ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Quaternion_log>( a1);
    }

    inline Expression<Vector>::Ptr logUnit ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Quaternion_logUnit>( a1);
    }

    inline Expression<Quaternion>::Ptr pow ( 
            Expression<Quaternion>::Ptr q, Expression<double>::Ptr s) {
        return exp( s*log(q) );
    }

    inline Expression<Quaternion>::Ptr powUnit ( 
            Expression<Quaternion>::Ptr q, Expression<double>::Ptr s) {
        return exp( s*logUnit(q) );
    }

    inline Expression<Vector>::Ptr diffUnit( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2) {
        return conditional<Vector>( dot(q1,q2), 
                    logUnit(q2*conj(q1) )*Constant<double>(2.0),
                    logUnit(-q2*conj(q1) )*Constant<double>(2.0)
               );
    }
/*
    inline Expression<Quaternion>::Ptr diff( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2) {
        return (q2-q1);
    }
*/


    inline Expression<Quaternion>::Ptr toQuat( 
            Expression<Vector>::Ptr angle_and_axis) 
    {
        return exp(angle_and_axis*Constant<double>(1.0/2.0));
    }

    inline Expression<Quaternion>::Ptr toQuat( 
            Expression<Vector>::Ptr axis, Expression<double>::Ptr angle) 
    {
        return toQuat( normalized(axis)*angle );
        //return exp( normalized(axis)*(angle/Constant<double>(2.0)) );
    }

    inline Expression<Vector>::Ptr axisAngle( Expression<Quaternion>::Ptr q) {
        return logUnit(q)*Constant<double>(2.0);
    }

    inline Expression<Quaternion>::Ptr slerp( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2, 
            Expression<double>::Ptr s) {
        return conditional<Quaternion>( dot(q1,q2), 
                    q1*pow( inv(q1)*q2 , s),
                    q1*pow( -inv(q1)*q2 , s)
               );
    }

    inline Expression<Quaternion>::Ptr slerpUnit( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2, 
            Expression<double>::Ptr s) {
        return conditional<Quaternion>( dot(q1,q2), 
                    q1*powUnit( conj(q1)*q2 , s),
                    q1*powUnit( -conj(q1)*q2 , s)
               );
    }

    inline Expression<Quaternion>::Ptr toQuat( Expression<Rotation>::Ptr a) {
        return boost::make_shared<QuaternionFromRotation>( a );
    }

    inline Expression<Rotation>::Ptr toRot( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<RotationFromQuaternion>( a );
    }

    inline Expression<Quaternion>::Ptr quaternion ( 
            Expression<double>::Ptr a1, 
            Expression<Vector>::Ptr a2 ) {
        return boost::make_shared<Quaternion_Scalar_Vector>( a1, a2 );
    }

    inline Expression<Quaternion>::Ptr quaternion ( 
            Expression<Vector>::Ptr a1, 
            Expression<double>::Ptr a2 ) {
        return boost::make_shared<Quaternion_Scalar_Vector>( a2, a1 );
    }






    

}; // namespace KDL

#endif



