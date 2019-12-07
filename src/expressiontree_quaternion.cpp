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

#include <kdl/expressiontree_quaternion.hpp>
#include "util.hpp"

namespace KDL {

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


} // namespace KDL


