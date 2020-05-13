/*
 * expressiontree_quaternion.cpp
 *
 *  Created on: September 
 *      Author: Erwin Aertbelien
 *
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
*/

#include <expressiongraph/expressiontree_quaternion.hpp>
#include "util.hpp"

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


 
Expression<Quaternion>::Ptr Quaternion_conj::derivativeExpression(int i) {
        //int nr = getDep<Quaternion>(i,argument);
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<Quaternion>(Quaternion(0,0,0,0));
        } else {
            return conj( argument->derivativeExpression(i) );
        }
}

Expression<Quaternion>::Ptr Quaternion_Scalar_Vector::derivativeExpression(int i) {
        int nr = getDep2(i,argument1,argument2);
        if (nr==1) {
            return Constant<Quaternion>(Quaternion::Zero());
        } if (nr==2) {
            return quaternion(Constant(0.0), argument2->derivativeExpression(i)); 
        } if (nr==3) {
            return quaternion(Constant(Vector::Zero()), argument1->derivativeExpression(i)); 
        } else {
            return quaternion( argument1->derivativeExpression(i),
                               argument2->derivativeExpression(i) );
        }

}

Expression<double>::Ptr Quaternion_w::derivativeExpression(int i) {
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant(0.0);
        } else {
            return w(argument->derivativeExpression(i)); 
        }
}

Expression<Vector>::Ptr Quaternion_vec::derivativeExpression(int i) {
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return vec(argument->derivativeExpression(i));
        }
}

Expression<Quaternion>::Ptr Minus_QuaternionQuaternion::derivativeExpression(int i){
        int nr = getDep<Quaternion>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Quaternion::Zero());
        } if (nr==2) {
            return  -argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i);
        } else {
            return argument1->derivativeExpression(i) - argument2->derivativeExpression(i);
        }
}

Expression<Quaternion>::Ptr Quaternion_Minus::derivativeExpression(int i){
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<Quaternion>(Quaternion(0,0,0,0));
        } else {
            return -argument->derivativeExpression(i);
        }
}

Expression<Quaternion>::Ptr Div_QuaternionScalar::derivativeExpression(int i){
        int nr = getDep2(i,argument1,argument2);
        Expression<double>::Ptr a2 = cached<double>(argument2);
        if (nr==1) {
            return Constant<Quaternion>(Quaternion::Zero());
        } if (nr==2) {
            return  - argument1*a2->derivativeExpression(i)/(a2*a2);
        } if (nr==3) {
            return argument1->derivativeExpression(i)/ a2;
        } else {
            return  argument1->derivativeExpression(i)/a2
                    - argument1*a2->derivativeExpression(i)/(a2*a2);
        }
}

Expression<double>::Ptr Dot_QuaternionQuaternion::derivativeExpression(int i){
        int nr = getDep2(i,argument1,argument2);
        if (nr==1) {
            return Constant<double>(0.0);
        } if (nr==2) {
            return  dot(argument1,argument2->derivativeExpression(i));
        } if (nr==3) {
            return  dot(argument2,argument1->derivativeExpression(i));
        } else {
            return  dot(argument1, argument2->derivativeExpression(i))
                    + dot(argument1->derivativeExpression(i),argument2) ;
        }
}

Expression<Quaternion>::Ptr Normalize_Quaternion::derivativeExpression(int i){
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<Quaternion>(Quaternion::Zero());
        } else {
            Expression<double>::Ptr n = cached<double>(norm(argument)+Constant(QUAT_EPS));
            Expression<Quaternion>::Ptr result = cached<Quaternion>(argument/n);
            Expression<double>::Ptr dn = dot( argument, argument->derivativeExpression(i) )/n;
            return (argument->derivativeExpression(i)-result*dn)/n;
        }
}

Expression<Vector>::Ptr Normalize_Vector::derivativeExpression(int i){
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            Expression<double>::Ptr n = cached<double>(Constant(1.0)/(norm(argument)+Constant(QUAT_EPS)));
            Expression<Vector>::Ptr v = cached<Vector>(argument);
            return (argument*n)->derivativeExpression(i);
        }
}

Expression<Vector>::Ptr Apply_QuaternionVector::derivativeExpression(int i){
        int nr = getDep2(i,argument1,argument2);
        Expression<Quaternion>::Ptr q = cached<Quaternion>(argument1);
        if (nr==1) {
            return Constant(Vector::Zero());
        } else {
            return vec((cached<Quaternion>(q* argument2)*conj(q))->derivativeExpression(i));
        }
}

Expression<Quaternion>::Ptr Quaternion_exp::derivativeExpression(int i){
    int nr = getDep(i,argument);
    if (nr==1) {
        return Constant(Quaternion::Zero());
    } else {
        Expression<Quaternion>::Ptr q = cached<Quaternion>( 
            argument
        );
        Expression<Quaternion>::Ptr d_q = cached<Quaternion>( 
            argument->derivativeExpression(i)
        );
        Expression<double>::Ptr n  = cached<double>( 
            norm(vec(q)) +Constant(QUAT_EPS) 
        );
        Expression<double>::Ptr d_n  = cached<double>( 
            dot( vec(q), vec(d_q))/n
        );
        Expression<double>::Ptr c    = cached<double>( cos(n)  );
        Expression<double>::Ptr s    = cached<double>( sin(n)  );
        Expression<double>::Ptr d_c  = cached<double>( -s*d_n  );
        Expression<double>::Ptr sc   = cached<double>( s/n     );
        Expression<double>::Ptr d_sc = cached<double>( 
            (c-sc)/n*d_n
        );
        Expression<double>::Ptr e    = cached<double>( exp(w(q))   );
        Expression<double>::Ptr d_e  = cached<double>( e*w(d_q)    );
        Expression<double>::Ptr e_sc = cached<double>( e*sc        );
        Expression<double>::Ptr d_e_sc  = cached<double>( 
            e*d_sc +d_e*sc
        );
        return cached<Quaternion>( quaternion(
                d_e*c + e*d_c,
                e_sc*vec(d_q) + d_e_sc*vec(q)
                ));
    }
}

Expression<Quaternion>::Ptr Quaternion_expv::derivativeExpression(int i){
    int nr = getDep(i,argument);
    if (nr==1) {
        return Constant(Quaternion::Zero());
    } else {
        Expression<Vector>::Ptr q = cached<Vector>( argument);
        Expression<Vector>::Ptr d_q = cached<Vector>( 
            argument->derivativeExpression(i)
        );
        Expression<double>::Ptr n  = cached<double>( 
            norm(q) +Constant(QUAT_EPS) 
        );
        Expression<double>::Ptr d_n  = cached<double>( 
            dot( q, d_q)/n
        );
        Expression<double>::Ptr c    = cached<double>( cos(n)  );
        Expression<double>::Ptr s    = cached<double>( sin(n)  );
        Expression<double>::Ptr d_c  = cached<double>( -s*d_n  );
        Expression<double>::Ptr sc   = cached<double>( s/n     );
        Expression<double>::Ptr d_sc = cached<double>( 
            (c-sc)/n*d_n
        );
        return cached<Quaternion>( quaternion(
                d_c,
                sc*d_q + d_sc*q
                ));
    }
}

Expression<Quaternion>::Ptr Quaternion_log::derivativeExpression(int i) {
    int nr = getDep(i,argument);
    if (nr==1) {
        return Constant(Quaternion::Zero());
    } else {
        Expression<Quaternion>::Ptr q = cached<Quaternion>( argument);
        Expression<Vector>::Ptr qv = cached<Vector>( vec(q) );
        Expression<double>::Ptr qw = cached<double>( w(q) );
        Expression<Quaternion>::Ptr d_q = cached<Quaternion>( 
            argument->derivativeExpression(i)
        );
        Expression<Vector>::Ptr d_qv = cached<Vector>( 
            vec(d_q)
        );
        Expression<double>::Ptr m  = cached<double>( 
            norm(q) +Constant(QUAT_EPS) 
        );
        Expression<double>::Ptr d_m  = cached<double>( 
            dot(q,d_q)/m
        );
        Expression<double>::Ptr n  = cached<double>( 
            norm(qv) +Constant(QUAT_EPS) 
        );
        Expression<double>::Ptr d_n  = cached<double>( 
            dot(qv, d_qv)/n
        );
        Expression<double>::Ptr n_inv = cached<double>( 
            Constant(1.0)/n 
        ); 
        Expression<Vector>::Ptr qvec_n  = cached<Vector>( 
            qv*n_inv 
        );
        Expression<Vector>::Ptr d_qvec_n  = cached<Vector>( 
            d_qv*n_inv - qv*cached<double>(d_n/(n*n))
        );
        Expression<double>::Ptr angle = cached<double>( 
            atan2(n,qw )
        );
        Expression<double>::Ptr d_angle = cached<double>( 
            ( d_n * qw - n * w(d_q) ) / ( n*n + qw*qw )
        );
        return cached<Quaternion>( quaternion(
            d_m/m,
            d_qvec_n*angle + qvec_n*d_angle 
        ));
    }
}

Expression<Vector>::Ptr Quaternion_logUnit::derivativeExpression(int i){
    int nr = getDep(i,argument);
    if (nr==1) {
        return Constant(Vector::Zero());
    } else {
        Expression<Quaternion>::Ptr q = cached<Quaternion>( argument);
        Expression<Vector>::Ptr qv = cached<Vector>( vec(q) );
        Expression<double>::Ptr qw = cached<double>( w(q) );
        Expression<Quaternion>::Ptr d_q = cached<Quaternion>( 
            argument->derivativeExpression(i)
        );
        Expression<Vector>::Ptr d_qv = cached<Vector>( 
            vec(d_q)
        );
        Expression<double>::Ptr n  = cached<double>( 
            norm(qv) +Constant(QUAT_EPS) 
        );
        Expression<double>::Ptr d_n  = cached<double>( 
            dot(qv, d_qv)/n
        );
        Expression<double>::Ptr n_inv = cached<double>( 
            Constant(1.0)/n 
        ); 
        Expression<Vector>::Ptr qvec_n  = cached<Vector>( 
            qv*n_inv 
        );
        Expression<Vector>::Ptr d_qvec_n  = cached<Vector>( 
            d_qv*n_inv - qv*cached<double>(d_n/(n*n))
        );
        Expression<double>::Ptr angle = cached<double>( 
            atan2(n,qw )
        );
        Expression<double>::Ptr d_angle = cached<double>( 
            ( d_n * qw - n * w(d_q) ) 
        );
        return cached<Vector>( 
            d_qvec_n*angle + qvec_n*d_angle 
        );
    }
}

Expression<Quaternion>::Ptr 
Sum_QuaternionQuaternion::derivativeExpression(int i) {
        int nr = getDep<Quaternion>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Quaternion::Zero());
        } if (nr==2) {
            return  argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i);
        } else {
            return argument1->derivativeExpression(i) + argument2->derivativeExpression(i);
        }
}


Expression<Quaternion>::Ptr 
Mult_QuaternionQuaternion::derivativeExpression(int i) {
        int nr = getDep<Quaternion>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Quaternion::Zero());
        } if (nr==2) {
            return argument1 * argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i) * argument2;
        } else {
            return argument1 * argument2->derivativeExpression(i) + argument1->derivativeExpression(i)*argument2;
        }
}

Expression<Quaternion>::Ptr 
Mult_QuaternionVector::derivativeExpression(int i) {
        int nr = getDep2<Quaternion,Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant(Quaternion::Zero());
        } if (nr==2) {
            return argument1 * argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i) * argument2;
        } else {
            return argument1 * argument2->derivativeExpression(i) + argument1->derivativeExpression(i)*argument2;
        }
}

Expression<Quaternion>::Ptr 
Mult_VectorQuaternion::derivativeExpression(int i) {
    int nr = getDep2<Vector,Quaternion>(i,argument1,argument2);
    if (nr==1) {
        return Constant(Quaternion::Zero());
    } if (nr==2) {
        return argument1 * argument2->derivativeExpression(i);
    } if (nr==3) {
        return argument1->derivativeExpression(i) * argument2;
    } else {
        return argument1 * argument2->derivativeExpression(i) + argument1->derivativeExpression(i)*argument2;
    }
}


Expression<Quaternion>::Ptr 
Mult_QuaternionScalar::derivativeExpression(int i) {
    int nr = getDep2<Quaternion,double>(i,argument1,argument2);
    if (nr==1) {
        return Constant(Quaternion::Zero());
    } if (nr==2) {
        return argument1 * argument2->derivativeExpression(i);
    } if (nr==3) {
        return argument1->derivativeExpression(i) * argument2;
    } else {
        return argument1 * argument2->derivativeExpression(i) + argument1->derivativeExpression(i)*argument2;
    }
}

Expression<Quaternion>::Ptr QuaternionFromRotation::derivativeExpression(int i) {
    int nr = getDep(i,argument);
    if (nr==1) {
        return Constant<Quaternion>(Quaternion::Zero());
    } else {
        return cached<Quaternion>(
            quaternion( Constant(0.0), argument->derivativeExpression(i)*Constant(0.5))*toQuat(argument)
        );
    }
}


Expression<Vector>::Ptr 
RotationFromQuaternion::derivativeExpression(int i) {
    int nr = getDep(i,argument);
    if (nr==1) {
        return Constant<Vector>(Vector::Zero());
    } else {
        return cached<Vector>(
                vec( argument->derivativeExpression(i) *conj(argument) )*Constant(2.0)
              );
    }
}

Expression<double>::Ptr 
Quaternion_squared_norm::derivativeExpression(int i) {
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return cached<double>( Constant(2.0)*dot(argument, argument->derivativeExpression(i)) );
        }
}

Expression<double>::Ptr 
Quaternion_norm::derivativeExpression(int i) {
        int nr = getDep(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            return cached<double>( dot(argument, argument->derivativeExpression(i)) / norm(argument) );
        }

}


//====================================== operators =================================================
    Expression<Quaternion>::Ptr conj( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_conj>(a);
    }

    Expression<double>::Ptr w( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_w>(a);
    }

    Expression<Vector>::Ptr vec( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_vec>(a);
    }

    Expression<double>::Ptr squared_norm( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_squared_norm>(a);
    }

    Expression<double>::Ptr norm( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_norm>(a);
    }

    Expression<Quaternion>::Ptr operator+( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Sum_QuaternionQuaternion>( a1, a2 );
    }

    Expression<Quaternion>::Ptr operator-( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Minus_QuaternionQuaternion>( a1, a2 );
    }

    Expression<Quaternion>::Ptr operator-( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<Quaternion_Minus>(a);
    }
 
    Expression<Quaternion>::Ptr operator* ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<double>::Ptr a2 ) {
        return boost::make_shared<Mult_QuaternionScalar>( a1, a2 );
    }

    Expression<Quaternion>::Ptr operator* ( 
            Expression<double>::Ptr a2 ,
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Mult_QuaternionScalar>( a1, a2 );
    }

    Expression<Quaternion>::Ptr operator* ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Mult_QuaternionQuaternion>( a1, a2 );
    }
    Expression<Quaternion>::Ptr operator* ( 
            Expression<Vector>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Mult_VectorQuaternion>( a1, a2 );
    }

    Expression<Quaternion>::Ptr operator* ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Vector>::Ptr a2 ) {
        return boost::make_shared<Mult_QuaternionVector>( a1, a2 );
    }

    Expression<Quaternion>::Ptr operator/ ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<double>::Ptr a2 ) {
        return boost::make_shared<Div_QuaternionScalar>( a1, a2 );
    }
    
    Expression<Quaternion>::Ptr inv( Expression<Quaternion>::Ptr q) {
        return conj(q) / squared_norm(q);
    }

    Expression<double>::Ptr dot ( 
            Expression<Quaternion>::Ptr a1, 
            Expression<Quaternion>::Ptr a2 ) {
        return boost::make_shared<Dot_QuaternionQuaternion>( a1, a2 );
    }

    Expression<Quaternion>::Ptr normalized ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Normalize_Quaternion>( a1);
    }

    Expression<Vector>::Ptr normalized ( 
            Expression<Vector>::Ptr a1) {
        return boost::make_shared<Normalize_Vector>( a1);
    }

    Expression<Vector>::Ptr apply ( 
            Expression<Quaternion>::Ptr q, Expression<Vector>::Ptr v) {
        return boost::make_shared<Apply_QuaternionVector>(q,v);
    }

    Expression<Quaternion>::Ptr exp ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Quaternion_exp>( a1);
    }

    Expression<Quaternion>::Ptr exp ( 
            Expression<Vector>::Ptr a1) {
        return boost::make_shared<Quaternion_expv>( a1);
    }


    Expression<Quaternion>::Ptr log ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Quaternion_log>( a1);
    }

    Expression<Vector>::Ptr logUnit ( 
            Expression<Quaternion>::Ptr a1) {
        return boost::make_shared<Quaternion_logUnit>( a1);
    }

    Expression<Quaternion>::Ptr pow ( 
            Expression<Quaternion>::Ptr q, Expression<double>::Ptr s) {
        return exp( s*log(q) );
    }

    Expression<Quaternion>::Ptr powUnit ( 
            Expression<Quaternion>::Ptr q, Expression<double>::Ptr s) {
        return exp( s*logUnit(q) );
    }

    Expression<Vector>::Ptr diffUnit( 
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


    Expression<Quaternion>::Ptr toQuat( 
            Expression<Vector>::Ptr angle_and_axis) 
    {
        return exp(angle_and_axis*Constant<double>(1.0/2.0));
    }

    Expression<Quaternion>::Ptr toQuat( 
            Expression<Vector>::Ptr axis, Expression<double>::Ptr angle) 
    {
        return toQuat( normalized(axis)*angle );
        //return exp( normalized(axis)*(angle/Constant<double>(2.0)) );
    }

    Expression<Vector>::Ptr axisAngle( Expression<Quaternion>::Ptr q) {
        return logUnit(q)*Constant<double>(2.0);
    }

    Expression<Quaternion>::Ptr slerp( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2, 
            Expression<double>::Ptr s) {
        return conditional<Quaternion>( dot(q1,q2), 
                    q1*pow( inv(q1)*q2 , s),
                    q1*pow( -inv(q1)*q2 , s)
               );
    }

    Expression<Quaternion>::Ptr slerpUnit( 
            Expression<Quaternion>::Ptr q1, 
            Expression<Quaternion>::Ptr q2, 
            Expression<double>::Ptr s) {
        return conditional<Quaternion>( dot(q1,q2), 
                    q1*powUnit( conj(q1)*q2 , s),
                    q1*powUnit( -conj(q1)*q2 , s)
               );
    }

    Expression<Quaternion>::Ptr toQuat( Expression<Rotation>::Ptr a) {
        return boost::make_shared<QuaternionFromRotation>( a );
    }

    Expression<Rotation>::Ptr toRot( Expression<Quaternion>::Ptr a) {
        return boost::make_shared<RotationFromQuaternion>( a );
    }

    Expression<Quaternion>::Ptr quaternion ( 
            Expression<double>::Ptr a1, 
            Expression<Vector>::Ptr a2 ) {
        return boost::make_shared<Quaternion_Scalar_Vector>( a1, a2 );
    }

    Expression<Quaternion>::Ptr quaternion ( 
            Expression<Vector>::Ptr a1, 
            Expression<double>::Ptr a2 ) {
        return boost::make_shared<Quaternion_Scalar_Vector>( a2, a1 );
    }





} // namespace KDL


