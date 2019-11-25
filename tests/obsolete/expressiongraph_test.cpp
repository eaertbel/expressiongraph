/*
 * expressiontree_test.cpp
 *
 *  Created on: Aug. 2012
 *  Author: Erwin Aertbelien 
 */


#include <kdl/expressiontree.hpp>
#include "expressiongraph_test.hpp"


using namespace KDL;


TEST(ExpressionTree, Scalars) {
        std::vector<int> ndx; 
        ndx.push_back(0);ndx.push_back(2);ndx.push_back(3);
        Expression<double>::Ptr a = random<double>(ndx);
        Expression<double>::Ptr b = random<double>(ndx);
        CHECK_ROT_WITH_NUM( -a );
        CHECK_ROT_WITH_NUM( sin(a) );
        CHECK_ROT_WITH_NUM( cos(a) );
        CHECK_ROT_WITH_NUM( tan(a) );
        CHECK_ROT_WITH_NUM( exp(a) );
        CHECK_ROT_WITH_NUM( asin(a) );
        CHECK_ROT_WITH_NUM( acos(a*a+Constant(0.01)) );
        CHECK_ROT_WITH_NUM( exp(a) );
        CHECK_ROT_WITH_NUM( log(a*a) );
        CHECK_ROT_WITH_NUM( sqr(a) );
        CHECK_ROT_WITH_NUM( sqrt(a*a+Constant(0.001)) );
        CHECK_ROT_WITH_NUM( abs(a) );
        CHECK_ROT_WITH_NUM( a + b );
        CHECK_ROT_WITH_NUM( a - b );
        CHECK_ROT_WITH_NUM( a * b );
        CHECK_ROT_WITH_NUM( a / b );
        CHECK_ROT_WITH_NUM( atan2(a,b)  );

        EXPECT_EQ_VALUES( a*b, b*a );
        EXPECT_EQ_VALUES( a+b, b+a );
        EXPECT_EQ_VALUES( a+Constant(0.0), a );
        EXPECT_EQ_VALUES( a*Constant(1.0), a );
        EXPECT_EQ_VALUES( a-b, -(b-a) );
        EXPECT_EQ_VALUES( a/b, a*(Constant(1.0)/b) );
        EXPECT_EQ_VALUES( a/b, (Constant(1.0)/b)*a );

        EXPECT_EQ_VALUES( a*a/Constant(10.0), acos(cos(a*a/Constant(10.0))) );
        EXPECT_EQ_VALUES( a/Constant(10.0), asin(sin(a/Constant(10.0))) );
        EXPECT_EQ_VALUES( a/Constant(10.0), atan(tan(a/Constant(10.0))) );
        EXPECT_EQ_VALUES( a*a, sqr(sqrt(a*a)) );
        EXPECT_EQ_VALUES( abs(a), sqrt(sqr(a)) );
        EXPECT_EQ_VALUES( atan2(sin(a),cos(a)) , a );
}

        
TEST(ExpressionTree, Vector) {
        std::vector<int> ndx; 
        ndx.push_back(1);ndx.push_back(2);ndx.push_back(3);
        Expression<Vector>::Ptr a = random<Vector>(ndx);
        Expression<Vector>::Ptr b = random<Vector>(ndx);
        Expression<double>::Ptr s = random<double>(ndx);
 
        CHECK_ROT_WITH_NUM( -a );
        CHECK_ROT_WITH_NUM( dot(a,b) );
        CHECK_ROT_WITH_NUM( a*b );
        CHECK_ROT_WITH_NUM( a-b );
        CHECK_ROT_WITH_NUM( a+b );
        CHECK_ROT_WITH_NUM( norm(a) );
        CHECK_ROT_WITH_NUM( squared_norm(a) );
        CHECK_ROT_WITH_NUM( a*s );
        CHECK_ROT_WITH_NUM( s*a );
        CHECK_ROT_WITH_NUM( coord_x(a) );
        CHECK_ROT_WITH_NUM( coord_y(a) );
        CHECK_ROT_WITH_NUM( coord_z(a) );
        CHECK_ROT_WITH_NUM( diff(a,b) );
 
        EXPECT_EQ_VALUES( a*b, -b*a );
        EXPECT_EQ_VALUES( dot(a,b), dot(b,a) );
        EXPECT_EQ_VALUES( a+b, b+a );
        EXPECT_EQ_VALUES( a-b, -(b-a) );
        EXPECT_EQ_VALUES( a*s, s*a );
        //EXPECT_EQ_VALUES( sqr(a), a*a );
        EXPECT_EQ_VALUES( squared_norm(a), dot(a,a) );
        EXPECT_EQ_VALUES( diff(a,b), -diff(b,a) );
        EXPECT_EQ_VALUES( a + Constant(Vector::Zero()), a );

}

TEST(ExpressionTree, Rotation) {
        // declare random variables to be used:
        std::vector<int> ndx; 
        ndx.push_back(1);
        ndx.push_back(2);ndx.push_back(3);
        Expression<Rotation>::Ptr a = random<Rotation>(ndx);
        Expression<Rotation>::Ptr b = random<Rotation>(ndx);
        Expression<Rotation>::Ptr I = Constant( Rotation::Identity() );
        Expression<Vector>::Ptr   v = random<Vector>(ndx);
        Expression<double>::Ptr   s = random<double>(ndx);
        Vector axis(1,2,3);
        random(axis);

        
        CHECK_ROT_WITH_NUM( rot(axis,s) );
        CHECK_ROT_WITH_NUM( rot_x(s) );
        CHECK_ROT_WITH_NUM( rot_y(s) );
        CHECK_ROT_WITH_NUM( rot_z(s) );
        CHECK_ROT_WITH_NUM( inv(a) );
        CHECK_ROT_WITH_NUM( a*b );
        CHECK_ROT_WITH_NUM( a*v );
        CHECK_ROT_WITH_NUM( unit_x(a) );
        CHECK_ROT_WITH_NUM( unit_y(a) );
        CHECK_ROT_WITH_NUM( unit_z(a) );
        
        EXPECT_EQ_VALUES( a*I, I*a );
        EXPECT_EQ_VALUES( a*inv(a), inv(a)*a );
        EXPECT_EQ_VALUES( a*inv(a), I );
        EXPECT_EQ_VALUES( (a*b)*v, a*(b*v) );
        EXPECT_EQ_VALUES( rot_x(s), rot(Vector(1,0,0),s) );
        EXPECT_EQ_VALUES( rot_y(s), rot(Vector(0,1,0),s) );
        EXPECT_EQ_VALUES( rot_z(s), rot(Vector(0,0,1),s) );
}

TEST(ExpressionTree, Frame) {
        // declare random variables to be used:
        std::vector<int> ndx; 
        ndx.push_back(0);
        ndx.push_back(1);ndx.push_back(2);
        Expression<Frame>::Ptr a = random<Frame>(ndx);
        Expression<Frame>::Ptr b = random<Frame>(ndx);
        Expression<Frame>::Ptr I = Constant( Frame::Identity() );
        Expression<Rotation>::Ptr R = random<Rotation>(ndx);
        Expression<Vector>::Ptr   v = random<Vector>(ndx);
        Expression<double>::Ptr   s = random<double>(ndx);

        
        CHECK_ROT_WITH_NUM( frame(R,v) );
        CHECK_ROT_WITH_NUM( frame(R) );
        CHECK_ROT_WITH_NUM( frame(v) );
        CHECK_ROT_WITH_NUM( inv(a) );
        CHECK_ROT_WITH_NUM( a*b );
        CHECK_ROT_WITH_NUM( a*v );
        CHECK_ROT_WITH_NUM( origin(a) );
        CHECK_ROT_WITH_NUM( rotation(a) );

        EXPECT_EQ_VALUES( a*I, I*a );
        EXPECT_EQ_VALUES( a*inv(a), inv(a)*a );
        EXPECT_EQ_VALUES( a*inv(a), I );
        EXPECT_EQ_VALUES( (a*b)*v, a*(b*v) );
}

TEST(ExpressionTree, Twist) {
        // declare random variables to be used:
        std::vector<int> ndx; 
        ndx.push_back(1);
        //ndx.push_back(2);ndx.push_back(3);
        Expression<Twist>::Ptr a = random<Twist>(ndx);
        Expression<Twist>::Ptr b = random<Twist>(ndx);
        Expression<Twist>::Ptr Z = Constant( Twist::Zero() );
        Expression<Frame>::Ptr F1 = random<Frame>(ndx);
        Expression<Frame>::Ptr F2 = random<Frame>(ndx);
        Expression<Rotation>::Ptr R1 = random<Rotation>(ndx);
        Expression<Rotation>::Ptr R2 = random<Rotation>(ndx);
        Expression<Vector>::Ptr   v1 = random<Vector>(ndx);
        Expression<Vector>::Ptr   v2 = random<Vector>(ndx);
        Expression<double>::Ptr   s = random<double>(ndx);

        
        CHECK_ROT_WITH_NUM( twist(v1,v2) );
        CHECK_ROT_WITH_NUM( transvel(a) );
        CHECK_ROT_WITH_NUM( rotvel(a) );
        CHECK_ROT_WITH_NUM( a+b );
        CHECK_ROT_WITH_NUM( a-b );
        CHECK_ROT_WITH_NUM( R1*a );
        CHECK_ROT_WITH_NUM( s*a );
        CHECK_ROT_WITH_NUM( a*s );
        CHECK_ROT_WITH_NUM( ref_point(a,v1) );

        EXPECT_EQ_VALUES( (R1*R2)*a, R1*(R2*a) );
        EXPECT_EQ_VALUES( a-b, -(b-a) );
        EXPECT_EQ_VALUES( a+a, Constant(2.0)*a );
        EXPECT_EQ_VALUES( a-a, Z );
        EXPECT_EQ_VALUES( a+Z, a );
        EXPECT_EQ_VALUES( Z+a, a );
        EXPECT_EQ_VALUES( transvel(twist(v1,v2)), v1 );
        EXPECT_EQ_VALUES( rotvel(twist(v1,v2)), v2 );
}

TEST(ExpressionTree, Wrench) {
        // declare random variables to be used:
        std::vector<int> ndx; 
        ndx.push_back(0);
        //ndx.push_back(2);ndx.push_back(3);
        Expression<Wrench>::Ptr a = random<Wrench>(ndx);
        Expression<Wrench>::Ptr b = random<Wrench>(ndx);
        Expression<Wrench>::Ptr Z = Constant( Wrench::Zero() );
        Expression<Frame>::Ptr F1 = random<Frame>(ndx);
        Expression<Frame>::Ptr F2 = random<Frame>(ndx);
        Expression<Rotation>::Ptr R1 = random<Rotation>(ndx);
        Expression<Rotation>::Ptr R2 = random<Rotation>(ndx);
        Expression<Vector>::Ptr   v1 = random<Vector>(ndx);
        Expression<Vector>::Ptr   v2 = random<Vector>(ndx);
        Expression<double>::Ptr   s = random<double>(ndx);

        
        CHECK_ROT_WITH_NUM( wrench(v1,v2) );
        CHECK_ROT_WITH_NUM( torque(a) );
        CHECK_ROT_WITH_NUM( force(a) );
        CHECK_ROT_WITH_NUM( a+b );
        CHECK_ROT_WITH_NUM( a-b );
        CHECK_ROT_WITH_NUM( R1*a );
        CHECK_ROT_WITH_NUM( s*a );
        CHECK_ROT_WITH_NUM( a*s );
        CHECK_ROT_WITH_NUM( ref_point(a,v1) );

        EXPECT_EQ_VALUES( (R1*R2)*a, R1*(R2*a) );
        EXPECT_EQ_VALUES( a-b, -(b-a) );
        EXPECT_EQ_VALUES( a+a, Constant(2.0)*a );
        EXPECT_EQ_VALUES( a-a, Z );
        EXPECT_EQ_VALUES( a+Z, a );
        EXPECT_EQ_VALUES( Z+a, a );
        EXPECT_EQ_VALUES( force(wrench(v1,v2)), v1 );
        EXPECT_EQ_VALUES( torque(wrench(v2,v2)), v2 );
}

class MonsterExpression : public ::testing::Test {
protected:
        std::vector<int> ndx; 
        Expression<Twist>::Ptr a;
        Expression<Wrench>::Ptr b;
        Expression<Wrench>::Ptr Z;
        Expression<Frame>::Ptr F;
        Expression<Rotation>::Ptr R1;
        Expression<Rotation>::Ptr R2;
        Expression<Vector>::Ptr   v1;
        Expression<Vector>::Ptr   v2;
        Expression<double>::Ptr   s;
        Expression<double>::Ptr   expr;



    virtual void SetUp() {
        ndx.push_back(0);
        //ndx.push_back(2);ndx.push_back(3);
        a  = random<Twist>(ndx);
        b  = random<Wrench>(ndx);
        Z  = Constant( Wrench::Zero() );
        s  = random<double>(ndx);
        Expression<double>::Ptr s1,s2,s3,s4;
        s1 = cached<double>( -sin(s)*cos(s)+tan(s)   );
        s2 = cached<double>( sin(s)+tan(s) )/s;
        s3 = s*cached<double>( s  );
        v1 = cached<Vector>(KDL::vector(s1,s2,s3));
        v2 = cached<Frame>(frame(rot_x(s3)*rot_z(s3)))*cached<Vector>(KDL::vector(s3,s2*s2,s1));
        expr = Constant(0.0001)*norm(v1*v2*dot(v1,v2));
    }
};

TEST_F(MonsterExpression, NumericalDerivative) {
        CHECK_ROT_WITH_NUM( expr );
}

TEST_F(MonsterExpression, DerivativeExpression) {
        CHECK_ROT_WITH_NUM( expr->derivativeExpression(0) );
        CHECK_ROT_WITH_NUM( expr->derivativeExpression(1) );
        CHECK_ROT_WITH_NUM( expr->derivativeExpression(2) );
        CHECK_ROT_WITH_NUM( expr->derivativeExpression(3) );
        CHECK_ROT_WITH_NUM( expr->derivativeExpression(4) );
}

TEST_F(MonsterExpression, ClonedExpression) {
        Expression<double>::Ptr expr2 = expr->clone();
        EXPECT_EQ_VALUES( expr, expr2 );
        // check whether they influence each other:
        double val1;
        setArbitraryInput<double>( expr );         
        val1=expr->value();
        setArbitraryInput<double>( expr2 );         
        EXPECT_NEAR( expr->value(), val1, 1E-8 );
}

TEST_F(MonsterExpression, DerivativeExpression2) {
        setArbitraryInput<double>( expr );         
        Expression<double>::Ptr a = expr->derivativeExpression(0);
        Expression<double>::Ptr b = expr->derivativeExpression(3);
        a->value();
        b->value();
        EXPECT_NEAR( a->derivative(3), b->derivative(0), 1E-8 );
}


TEST_F(MonsterExpression, ExpressionOptimizer) {
        Expression<double>::Ptr expr2 = expr->clone();
        std::vector<int> ndx;
        std::vector<double> values(3); 
        ndx.push_back(0);
        ndx.push_back(2);
        ndx.push_back(3);

        ExpressionOptimizer opt;
        opt.prepare(ndx);
        expr->addToOptimizer(opt);
        
        random(values[0]);
        random(values[1]);
        random(values[2]);
        opt.setInputValues(values);
        expr2->setInputValues(ndx,values);
        EXPECT_NEAR( expr->value(), expr2->value(), 1E-8 );
        EXPECT_NEAR( expr->derivative(0), expr2->derivative(0), 1E-8 );
        EXPECT_NEAR( expr->derivative(1), expr2->derivative(1), 1E-8 );
        EXPECT_NEAR( expr->derivative(2), expr2->derivative(2), 1E-8 );
        EXPECT_NEAR( expr->derivative(3), expr2->derivative(3), 1E-8 );

         
        random(values[0]);
        random(values[1]);
        random(values[2]);
        opt.setInputValues(values);
        expr2->setInputValues(ndx,values);
        EXPECT_NEAR( expr->value(), expr2->value(), 1E-8 );
        EXPECT_NEAR( expr->derivative(0), expr2->derivative(0), 1E-8 );
        EXPECT_NEAR( expr->derivative(1), expr2->derivative(1), 1E-8 );
        EXPECT_NEAR( expr->derivative(2), expr2->derivative(2), 1E-8 );
        EXPECT_NEAR( expr->derivative(3), expr2->derivative(3), 1E-8 );

  
        random(values[0]);
        random(values[1]);
        random(values[2]);
        opt.setInputValues(values);
        expr2->setInputValues(ndx,values);
        EXPECT_NEAR( expr->value(), expr2->value(), 1E-8 );
        EXPECT_NEAR( expr->derivative(0), expr2->derivative(0), 1E-8 );
        EXPECT_NEAR( expr->derivative(1), expr2->derivative(1), 1E-8 );
        EXPECT_NEAR( expr->derivative(2), expr2->derivative(2), 1E-8 );
        EXPECT_NEAR( expr->derivative(3), expr2->derivative(3), 1E-8 );


}

TEST(RotationalInputs, Simple) {
    Expression<Rotation>::Ptr e = inputRot(0);
    CHECK_ROT_WITH_NUM( e );
}

TEST(RotationalInputs, Simple2) {
    Expression<Rotation>::Ptr e = Constant(Rotation::EulerZYX(0.1,0.2,0.3))*inputRot(0);
    CHECK_ROT_WITH_NUM( e );
}

TEST(RotationalInputs, Simple3) {
    Expression<Rotation>::Ptr e = Constant(Rotation::EulerZYX(0.1,0.2,0.3))*inputRot(0)*Constant(Rotation::EulerZYX(0.1,0.2,0.3));
    CHECK_ROT_WITH_NUM( e );
}

TEST(RotationalInputs, Simple4) {
    Expression<Rotation>::Ptr e = inputRot(0)*inputRot(3);
    CHECK_ROT_WITH_NUM( e );
    Expression<Rotation>::Ptr e2 = inputRot(3)*inputRot(0)*inputRot(3);
    CHECK_ROT_WITH_NUM( e2 );
}

TEST(RotationalInputs, Simple5) {
    Expression<Rotation>::Ptr e = rot_x(input(0))*rot_y(input(1))*rot_z(input(2));
    CHECK_ROT_WITH_NUM( e );
}
TEST(RotationalInputs, NumericalDerivatives) {
    Expression<double>::Ptr e = dot( Constant(Vector(0,0,1)), inputRot(0)*KDL::vector(input(4),Constant(0.0),input(3))) ;
    CHECK_ROT_WITH_NUM( e );
}

TEST(VariableType, Optimizer) {
    // build up two equivalent expressions a and b
    // we will manually fill in the values for a such that it corresponds to b
    // and evaluate a and b in the same expression.
    std::vector<int> ndx;
    ndx.push_back(0); 
    ndx.push_back(2); 
    ndx.push_back(3); 
    VariableType<double>::Ptr a = Variable<double>(ndx); 

    Expression<double>::Ptr b = sin(input(0))*input(2) + input(3);

    // expression where b is used:
    Expression<double>::Ptr expr_b = sqr(b)+cached<double>(norm(KDL::vector(Constant(0.0),Constant(1.0),b)));
    // expression where a is used:
    Expression<double>::Ptr expr_a = sqr(a)+cached<double>(norm(KDL::vector(Constant(0.0),Constant(1.0),a)));

    ExpressionOptimizer opt;
    opt.prepare(ndx);
    expr_a->addToOptimizer(opt);

    std::vector<double> arg(3);


    random(arg[0]);
    random(arg[1]);
    random(arg[2]);
    expr_b->setInputValues(ndx,arg); // also the values on b are now set.

    a->setValue( b->value() ); 
    a->setJacobian(0,b->derivative(0));
    a->setJacobian(1,b->derivative(2));
    a->setJacobian(2,b->derivative(3));
    opt.setInputValues(arg);  // you should call these, to clear possible caches.

    EXPECT_NEAR( expr_a->value(), expr_b->value(), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(0), expr_b->derivative(0), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(1), expr_b->derivative(1), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(2), expr_b->derivative(2), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(3), expr_b->derivative(3), 1E-8 );


    random(arg[0]);
    random(arg[1]);
    random(arg[2]);
    expr_b->setInputValues(ndx,arg); // also the values on b are now set.

    a->setValue( b->value() ); 
    a->setJacobian(0,b->derivative(0));
    a->setJacobian(1,b->derivative(2));
    a->setJacobian(2,b->derivative(3));
    opt.setInputValues(arg);  // you should call these, to clear possible caches.

    EXPECT_NEAR( expr_a->value(), expr_b->value(), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(0), expr_b->derivative(0), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(1), expr_b->derivative(1), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(2), expr_b->derivative(2), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(3), expr_b->derivative(3), 1E-8 );


    random(arg[0]);
    random(arg[1]);
    random(arg[2]);
    expr_b->setInputValues(ndx,arg); // also the values on b are now set.

    a->setValue( b->value() ); 
    a->setJacobian(0,b->derivative(0));
    a->setJacobian(1,b->derivative(2));
    a->setJacobian(2,b->derivative(3));
    opt.setInputValues(arg);  // you should call these, to clear possible caches.

    EXPECT_NEAR( expr_a->value(), expr_b->value(), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(0), expr_b->derivative(0), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(1), expr_b->derivative(1), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(2), expr_b->derivative(2), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(3), expr_b->derivative(3), 1E-8 );



    }

TEST(VariableType, Consistency) {
    // build up two equivalent expressions a and b
    // we will manually fill in the values for a such that it corresponds to b
    // and evaluate a and b in the same expression.
    std::vector<int> ndx;
    ndx.push_back(0); 
    ndx.push_back(2); 
    ndx.push_back(3); 
    VariableType<double>::Ptr a = Variable<double>(ndx); 

    Expression<double>::Ptr b = sin(input(0))*input(2) + input(3);

    // expression where b is used:
    Expression<double>::Ptr expr_b = sqr(b)+cached<double>(norm(KDL::vector(Constant(0.0),Constant(1.0),b)));
    // expression where a is used:
    Expression<double>::Ptr expr_a = sqr(a)+cached<double>(norm(KDL::vector(Constant(0.0),Constant(1.0),a)));

    std::vector<double> arg(3);
    random(arg[0]);
    random(arg[1]);
    random(arg[2]);

    expr_b->setInputValues(ndx,arg); // also the values on b are now set.

    a->setValue( b->value() );
    a->setJacobian(0, b->derivative(0) );
    a->setJacobian(1, b->derivative(2) );
    a->setJacobian(2, b->derivative(3) );
    expr_a->setInputValues(ndx,arg); 
 
    EXPECT_NEAR( expr_a->value(), expr_b->value(), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(0), expr_b->derivative(0), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(1), expr_b->derivative(1), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(2), expr_b->derivative(2), 1E-8 );
    EXPECT_NEAR( expr_a->derivative(3), expr_b->derivative(3), 1E-8 );
}


// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
    epsilon=1E-12;
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}


