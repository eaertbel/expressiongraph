#include <iostream>
#include <kdl/frames.hpp>
#include <kdl/expressiontree.hpp>
#include <kdl/conversions.hpp>
#include <Eigen/Cholesky>
using namespace KDL;

/**
// for now: only a version with scalar K
void fill(Expression<Frame>::Ptr arg, int row, const Frame& desired, const Twist& desired_dot, double K, Eigen::MatrixXd& J, Eigen::VectorXd& rhs)
{
    arg->value(); // obligatory
    rhs.block(row,1,6,1) = toEigen( K*diff(arg->value(),desired)  + desired_dot );
    int N = min(J.cols(),arg->number_of_derivatives());
    for (int i=0;i<N;++i) {
        J.block(row,i,6,1)= toEigen(arg->derivative(i));      
    }
}
*/

/**
 * @TODO: Bring EulerAngles to Expression Trees
 * @TODO: Inequality solver
 * @TODO: Multiple priority solver
 */
class SimpleSolver {
        Eigen::MatrixXd JTJ;
        Eigen::LLT<Eigen::MatrixXd> llt;  // for Cholesky decomposition
        Eigen::MatrixXd J;
        Eigen::VectorXd rhs;
        int row;
    public:
    SimpleSolver(int constraints, int nvar) : JTJ(nvar,nvar), llt(nvar), J(constraints,nvar), rhs(constraints) {
        row = 0;
    }
    
    void reset() {
        row = 0;
    }

    void add(Expression<double>::Ptr arg, double K,double desired=0.0, double desired_dot=0.0 ) {
        rhs(row) = K*(desired - arg->value()) + desired_dot;
        int N = min(J.cols(),arg->number_of_derivatives());
        for (int i=0;i<N;++i) {
            J(row,i) = arg->derivative(i);      
        }
        row+=1;
    }

    void solve(Eigen::VectorXd& qdot) {
        // a somewhat naieve DLS-solver:
        JTJ = J.transpose()*J;
        double lambda=0.001;
        for (int j=0;j<JTJ.rows();++j) {
            JTJ(j,j) += lambda*lambda;
        }
        llt.compute(JTJ);
        qdot = llt.solve(J.transpose()*rhs);
    }
};




int main(int argc,char* argv[]) {
    Eigen::MatrixXd Jacobian;
    using namespace KDL;
    using namespace std;
    using namespace Eigen;
    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;

    Chain chain;
    chain.addSegment(Segment("Segment 0", Joint("Joint 0", Joint::RotZ),Frame(Vector(0,0,L1))));
    chain.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotZ),Frame(Vector(0,0,L2))));
    chain.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotZ),Frame(Vector(0,0,L3))));
    chain.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotX),Frame(Vector(0,0,L4))));
    chain.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),Frame(Vector(0,0,0))));

    Expression<Frame>::Ptr kinchain = cached<Frame>(
         Constant(Frame(Vector(0,0,0.5))) * kinematic_chain( chain, 0 ) * Constant(Frame(Vector(0,0,0.3)))
    );
    Expression<Frame>::Ptr  line                  = Constant( Frame(Rotation::Identity(), Vector(0.2,0.2,0.0))); // z axis of frame corresponds to line.
    // constraint 1:
    Expression<Vector>::Ptr tmp                   = cached<Vector>( inv(line) * origin(kinchain) );
    Expression<double>::Ptr distance_to_line      = coord_x(tmp)*coord_x(tmp) + coord_y(tmp)*coord_y(tmp) - Constant(0.1*0.1);
    // constraint 2:
    Expression<double>::Ptr perpendicular_to_line = dot( unit_z(rotation(line)), unit_z(rotation(kinchain)));
    // constraint 3:
    Expression<double>::Ptr fixed_joint           = input(2) - Constant(0.5);

    // solving:
    int nvar = 7; 
    std::vector<double> joints(nvar);
    joints[0] = 0*M_PI*0.08;
    joints[1] = 1*M_PI*0.08;
    joints[2] = 0.5;
    joints[3] = 3*M_PI*0.08;
    joints[4] = 4*M_PI*0.08;
    joints[5] = 5*M_PI*0.08;
    joints[6] = 6*M_PI*0.08;   
    Eigen::VectorXd qdot(nvar);
    SimpleSolver solver(3,nvar);
    double dt=0.01; 
    double K = 4;
    for (double t=0;t<10;t+=dt) {
        solver.reset();
        distance_to_line->setInputValues(joints);
        solver.add(distance_to_line, K);
        perpendicular_to_line->setInputValues(joints);
        solver.add(perpendicular_to_line, K);
        fixed_joint->setInputValues(joints);
        solver.add(fixed_joint, K);
        
        solver.solve(qdot);

        // integration and print-out
        cout << distance_to_line->value() << "\t" 
             << perpendicular_to_line->value() << "\t"
             << fixed_joint->value() << "\t"; 
        for (int i=0;i<nvar;++i) {
            joints[i] +=  dt*qdot(i);
            cout << joints[i] << "\t";
        }
        cout << endl;
    }
}

