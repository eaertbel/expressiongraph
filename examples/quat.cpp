//#include <Eigen/Geometry>
#include <iostream>
#include <kdl/frames.hpp>
#include <kdl/quat_io.hpp>
#include <kdl/frames_io.hpp>
#include <kdl/conversions.hpp>
#include <kdl/quat.hpp>
#include <kdl/quat_io.hpp>
#include <kdl/expressiontree_quaternion.hpp>

using namespace std;
//using namespace Eigen;
using namespace KDL;

void test_operations(Quaternion q1, Quaternion q2) {
        cout << "================================================" << endl;
        cout << "Basic quaternion operations "<< endl;
        cout << "================================================" << endl;
     
        cout << "q1 " << q1 << endl;
        cout << "q2 " << q2 << endl;
        q1.normalize();
        q2.normalize();
        cout << "q1 normalized " << q1 << endl;
        cout << "q2 normalized " << q2 << endl;

        Rotation R1 = toRot(q1);
        cout << "R1 (from q1) \n" << R1 << endl;
        cout << "check orthormal \n" << R1*R1.Inverse() << endl;
        cout << "toQuat(R1) " << toQuat(R1) << endl;
        Rotation R2 = toRot(q2);
        cout << "R2 (from q2) \n" << R1 << endl;
        cout << "check orthonormal\n" << R2*R2.Inverse() << endl;

        cout << "R1*R2 \n" << (R1*R2) << endl;
        cout << "toRot(q1*q2)\n" << toRot(q1*q2) << endl;
  
        Vector a(-1,2,-3);
        cout<< "a " << a << " and norm of a " << a.Norm() << endl;
        cout << "R1*a " << (R1*a) << " and norm of R1*a " << (R1*a).Norm() << endl;
        cout << "conj(q) * a * q " << q1*Quaternion(0,a)*conj(q1) << endl;
        cout << "conj(q) * a * q " << q1*Quaternion(0,a)*inv(q1) << endl;
        cout << "appy(a,a) " << apply(q1,a) << endl;

        Vector axis,axis1,axis2;
        double angle,angle1,angle2;
        axis_angle(q1,axis,angle);
        cout << "axis 1 " << axis << " and angle 1 " << angle << endl;
        axis_angle(q2,axis,angle);
        cout << "axis 2 " << axis << " and angle 1 " << angle << endl;
        axis_angle( q1*q2, axis, angle );
        cout << "axis 1+2 " << axis << " and angle 1+2 " << angle << endl;
        cout << "back to quaternion " << toQuat(axis,angle) << endl;  

        angle = acos( dot(q1,q2) )*2.0;
        cout << "angle between (not so numerically stable ) : " << angle << endl;
        axis_angle( inv(q1)*q2, axis1, angle1 );
        cout << "alternative comp angle between  (numerically stable) : " << angle1 << endl;
        cout << "angle1 - angle " << angle1-angle << endl;

        Quaternion eq = exp(q1);
        cout << "exp(q1) " << eq << endl;
        cout << " log(exp(q1)) " << log(eq) << endl;
        cout << " q1 - log(exp(q1)) " << q1 - log(eq) << endl;
} 
 
void test_diff(double dt) {
    cout << "================================================" << endl;
    cout << "omega->dq, integrate q, compute R, differentiate R ->omega with dt="<< dt << endl;
    cout << "================================================" << endl;
    Quaternion q1(0.5,0.5,0.3,0.2);
    q1.normalize();
    cout << "q1 (normalized) " << q1 << endl;
    Vector omega1( 1.0,2.0,3.0 );
    cout << "omega " << omega1 << endl;
    // compute dq/dt = 0.5* w_bar * q
    Quaternion dq   = Quaternion(0, 0.5*omega1[0],0.5*omega1[1],0.5*omega1[2]) * q1;
    cout << "dq (0.5 * w_bar * q) = " << dq << endl;
    Quaternion q2     = q1 + dt*dq;
    cout << "q2 "<< q2 << endl;
    cout << "norm(q2) " << norm(q2) << endl;
    KDL::Rotation R1  = toRot(q1);
    KDL::Rotation R2  = toRot(q2);
    KDL::Vector omega = diff(R1,R2)/dt;   
    cout << "omega computed via rot. matrices " << omega << endl;

    //w = 2*dq*conj(q)
    Quaternion omegaq = 2.0* dq * conj(q1); 
    //Quaterniond omegaq( Quaterniond(2*dq.w(),2*dq.x(),2*dq.y(),2*dq.z()) * q1.conjugate() );
    cout << "omega computed from dq/dt (2*dq*conj(q1) ) " << omegaq << endl; 
    cout << "dot(q,dq) " << dot(q1,dq) << endl; 

    dq = toQuat(omega1*dt);
    cout << "Exact integration : dq=toQuat(omega1*dt) = " << dq << endl;
    cout << "Exact integration : alternative dq = exp(omega1*dt) = " << exp(Quaternion(0,omega1*dt/2)) << endl;
    q2 = dq*q1;
    cout << "q2 = q1*dq " <<  q2 << endl;
    R1  = toRot(q1);
    R2  = toRot(q2);
    omega = diff( R1, R2 )/dt;
    cout << "omega computed via rot. matrices forr exact diff  " << omega  << endl;
    cout << "diff using log " << log(q2*inv(q1) )/dt*2.0 << endl;
    cout << "diff using logUnit " << logUnit(q2*inv(q1) )/dt*2.0 << endl;
    cout << "diff using diff " << diffUnit(q1,q2)/dt << endl;
} 


void test_exp_log_2(const Quaternion& q) {

    cout << "Quaternion q = " << q << endl;
    Quaternion eq = exp(q);
    cout << "exp(q) " << eq << endl; 
    cout << "log(exp(q)) " << log(eq) << endl;
    cout << "log(exp(q))-q " << log(eq)-q << endl;
    cout << "norm of vec component/2/pi " << (log(eq)-q).vec.Norm()/2/PI << endl;
    Quaternion lq = log(q);
    cout << "log(q) " << lq << endl; 
    cout << "exp(log(q)) " << exp(lq) << endl; 
    cout << "exp(log(q))-q " << exp(lq)-q << endl; 
    cout << "norm of vec component/2/pi " << (exp(lq)-q).vec.Norm()/2/PI << endl;
    cout << endl;

    cout << "Normalized: " << endl;
    Quaternion q2 = q;
    q2.normalize();

    cout << "Vector part from q " << q.vec << endl;
    eq = exp(q2.vec);
    cout << "exp(q.vec) " << eq << endl; 
    cout << "log(exp(q2)) " << logUnit(eq) << endl;
    cout << "log(exp(q2))-q2 " << logUnit(eq)-q2.vec << endl;
    cout << "norm of vec component/2/pi " << (logUnit(eq)-q2.vec).Norm()/2/PI << endl;
    Vector lq2 = logUnit(q2);
    cout << "log(q2) " << lq2 << endl; 
    cout << "exp(log(q2)) " << exp(lq2) << endl; 
    cout << "exp(log(q2))-q2 " << exp(lq2)-q2 << endl; 
    cout << "norm of vec component/2/pi " << (exp(lq2)-q2).vec.Norm()/2/PI << endl;
    cout << endl;
}

void test_exp_log() {
    cout << "================================================" << endl;
    cout << "Testing log and exp " << endl;
    cout << "================================================" << endl;
    cout << "as expected, sometimes the difference is proportional to 2*PI " << endl;
    cout << "because the quaternions can only reprsent angles from 0 to 2*pi " << endl;
    test_exp_log_2( Quaternion(1,2,3,4) ); 
    test_exp_log_2( Quaternion(1,0,0,0) ); 
    test_exp_log_2( Quaternion(0,1,2,3) ); 
    test_exp_log_2( Quaternion(0,10,20,30) ); 
}


void test_diff(Quaternion q1, Quaternion q2, Quaternion q3) {
    cout << "================================================" << endl;
    cout << "Testing diff " << endl;
    cout << "================================================" << endl;
    Vector axis; double angle;
    axis_angle(q1,axis, angle);
    cout << "q1 : " << q1 << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl;
    axis_angle(q2,axis, angle);
    cout << "q2 : " << q2 << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl;
    axis_angle(q3,axis, angle);
    cout << "q3 : " << q3 << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl;
    axis = diffUnit(q1,q2);
    cout << "diff(q1,q2) : "  << axis << " norm [DEG] " << axis.Norm()/PI*180.0 << endl;
    axis = diffUnit(q1,-q2);
    cout << "diff(q1,-q2) : " << axis << " norm [DEG] " << axis.Norm()/PI*180.0 << endl;
    axis = diffUnit(q1,q3);
    cout << "diff(q1,q3) : "  << axis << " norm [DEG] " << axis.Norm()/PI*180.0 << endl;
    axis = diffUnit(q1,-q3);
    cout << "diff(q1,-q3) : " << axis << " norm [DEG] " << axis.Norm()/PI*180.0 << endl;
} 

void test_slerp(Quaternion q1, Quaternion q2) {
    cout << "================================================" << endl;
    cout << "Testing slerp: always along the smallest angle ! " << endl;
    cout << "================================================" << endl;
    Vector axis; double angle;
    axis_angle(q1,axis, angle);
    cout << "q1 : " << q1 << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl;
    axis_angle(q2,axis, angle);
    cout << "q2 : " << q2 << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl;
    for (double s=0.0; s <= 1.0; s+=0.1) {
        Quaternion q = slerp(q1, q2, s);
        axis_angle(q,axis, angle);
        cout << "q[s="<<s<<"] : " << q << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl;
    }
}


void test_toQuat(Rotation R) {
    cout << "================================================" << endl;
    cout << "Testing toQuat/toRot " << endl;
    cout << "================================================" << endl;
    cout << "Rotation matrix " << endl;
    cout << R << endl;  
    Quaternion q = toQuat(R);
    Vector axis; double angle;
    axis_angle(q,axis,angle);
    cout << "q : " << q << " ( angle (DEG) : " << angle/PI*180.0 << " and axis : " << axis << ")" << endl; 
    Rotation R2 = toRot(q);
    cout << "back to rotation matrix " << endl;
    cout << R2 << endl;
    cout << "inv(R)*R2 "<< endl << R.Inverse()*R2 << endl;
}


void test_exp_prop() {
    cout << "================================================" << endl;
    cout << "Testing properties of exponentials " << endl;
    cout << "================================================" << endl;
    Quaternion q1(1,2,3,4); 
    q1 = normalized(q1);
    Vector omega(1,2,3);
    //Quaternion dq = 0.5*Quaternion(omega)*q1;
    Quaternion dq(1,0,0,0);
    double dt = 1E-7;
    Quaternion deriv = (exp(q1+dq*dt) - exp(q1-dq*dt))/(2*dt);
    cout << "exp value         " << exp(q1) << endl;
    cout << "exp deriv numeric " << deriv << endl;

    Expression<Quaternion>::Ptr qe = exp( Constant(q1) + Constant(dq)*input(1) );
    qe->setInputValue(1,0.0);
    cout << "exp value using expr : " << qe->value() << endl; 
    cout << "exp deriv using expr : " << qe->derivative(1) << endl;

    deriv = (log(q1+dq*dt) - log(q1-dq*dt))/(2*dt);
    cout << "log value         " << exp(q1) << endl;
    cout << "log deriv numeric " << deriv << endl;
    qe = log( Constant(q1) + Constant(dq)*input(1) );
    qe->setInputValue(1,0.0);
    cout << "log value using expr : " << qe->value() << endl; 
    cout << "log deriv using expr : " << qe->derivative(1) << endl;

}

int main() {
    Rotation R1 = Rotation::EulerZYX(0.1,0,0);
    Rotation R2 = Rotation::EulerZYX(0.1+0.001,0,0);
    Rotation R3 = Rotation::EulerZYX(0.1,0,0);
    cout << R1 << endl;
    cout << R2 << endl;
    cout << R3 << endl;
    cout << diff(R1,R2,1.0) << endl;
    cout << diff(R1,R3,1.0) << endl;
    {
        Vector v1(1,2,3);
        Vector v2(3,2,3);
        Quaternion q1(0.5,v1[0],v1[1],v1[2]);
        Quaternion q2(0.5,v2[0],v2[1],v2[2]);
        test_operations(q1,q2);
    }
    {
        Quaternion q1 = toQuat(Vector(1,0,0), 0.9*PI);
        Quaternion q2 = toQuat(Vector(1,0,0), 0.1*PI);
        test_operations(q1,q2);
    }
    {
        Quaternion q1 = toQuat(Vector(1,0,0), 0.99*PI);
        Quaternion q2 = toQuat(Vector(1,0,0), 1*PI);
        test_operations(q1,q2);
    }
    {
        test_diff(1E-5);
        test_diff(1E-1);
    }
    {
        test_exp_log();
    }
    {
        // exploring edge cases:
        Quaternion q1 = toQuat( Vector(0,0,1), 30.0/180.0*PI );
        Quaternion q2 = toQuat( Vector(0,0,1), 60.0/180.0*PI );
        Quaternion q3 = toQuat( Vector(0,0,1), 170.0/180.0*PI );
        test_diff(q1,q2,q3);
        q1 = toQuat( Vector(1,2,3), 30.0/180.0*PI );
        q2 = toQuat( Vector(1,2,3), 60.0/180.0*PI );
        q3 = toQuat( Vector(1,2,3), 210.0/180.0*PI );
        test_diff(q1,q2,q3);
        q1 = toQuat( Vector(1,2,3), 30.0/180.0*PI );
        q2 = toQuat( Vector(1,2,3), (210+0.001)/180.0*PI );
        q3 = toQuat( Vector(1,2,3), (210-0.001)/180.0*PI );
        test_diff(q1,q2,q3);
        q1 = toQuat( Vector(1,2,3), 30.0/180.0*PI );
        q2 = toQuat( Vector(1,2,3), 390.0/180.0*PI );
        q3 = toQuat( Vector(1,2,3), (390+0.001)/180.0*PI );
        test_diff(q1,q2,q3);
    }
    {
        Quaternion q1 = toQuat( Vector(1,2,3), 30.0/180.0*PI );
        Quaternion q2 = toQuat( Vector(1,2,3), (150.0)/180.0*PI );
        test_slerp(q1,q2);
        q1 = toQuat( Vector(1,2,3), 0.0/180.0*PI );
        q2 = toQuat( Vector(1,2,3), (150.0)/180.0*PI );
        test_slerp(q1,q2);
        q1 = toQuat( Vector(1,2,3), 0.0/180.0*PI );
        q2 = toQuat( Vector(1,2,3), (180.0)/180.0*PI );
        test_slerp(q1,q2);
        q1 = toQuat( Vector(1,2,3), 0.0/180.0*PI );
        q2 = toQuat( Vector(1,2,3), (181.0)/180.0*PI );
        test_slerp(q1,q2);
        q1 = toQuat( Vector(1,0,0), 10.0/180.0*PI );
        q2 = toQuat( Vector(0,1,0), (10.0)/180.0*PI );
        test_slerp(q1,q2);
    }
    {
        Rotation R = Rotation::EulerZYX(0.1,0.2,0.3);
        test_toQuat(R);
        Rotation R2 = Rotation::EulerZYX(30.0/180*PI,0,0);
        test_toQuat(R2);
    }
    /*
    cout << "===================================================" << endl;
    cout << "derivative of q " << endl;
    cout << "===================================================" << endl;
    {
        // d(q^-1) =  q^-1  dq q^-1
        double dt = 0.00001;
        Quaternion q(0.4,0.5,0.6,0.5);
        Quaternion dq(0.1,0.2,0.3,0.2);
        Quaternion q2 = q + dt*dq;
        Quaternion dqi_numeric = (inv(q2) - inv(q))/dt;
        Quaternion dqi = inv(q)*dq*inv(q);
        cout << "q " << q.w << " " << q.vec.transpose() << endl;
        cout << "dq " << dq.w << " " << dq.vec.transpose() << endl;
        cout << "numeric deriv q^{-1} " << dqi_numeric.w << " " << dqi_numeric.vec.transpose() << endl; 
        cout << "symbolic deriv q^{-1} " << dqi.w << " " << dqi.vec.transpose() << endl; 
    }
    **/
    test_exp_prop(); 
}
