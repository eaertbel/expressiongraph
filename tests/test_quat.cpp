#include <kdl/expressiontree_quaternion.hpp>
#include <kdl/expressiontree.hpp>
#include <gtest/gtest.h>
#include <string.h>



namespace KDL {

int debug_output = 0;    // output additional debug values
int simple_values = 0;   // choose simple values, i.e. var. i has value i, otherwise random values are used.

/**
 * generate random entities for expression graphs of a given type.
 * \param [in] ndx, the random entities will depend on the variable numbers given by ndx.
 */
template <class T>
typename Expression<T>::Ptr random(std::vector<int>& ndx);


template<> Expression<double>::Ptr random<double>(std::vector<int>& ndx) {
    Expression<double>::Ptr result;
    if (ndx.size() > 0) {
        double a,b;
        random(a);
        random(b);
        result = Constant(a)*input(ndx[0])+Constant(b);
        for (size_t i=1;i<ndx.size();++i) {
            random(a);
            result = result + Constant(a)*input(ndx[i]);
        }
    }
    return result; 
}

template<> Expression<Vector>::Ptr random<Vector>(std::vector<int>& ndx) {
    return KDL::vector( random<double>(ndx), random<double>(ndx), random<double>(ndx) );
}

template<> Expression<Rotation>::Ptr random<Rotation>(std::vector<int>& ndx) {
    return rot_x( random<double>(ndx)) * rot_y(random<double>(ndx)) * rot_z(random<double>(ndx));
}

template<> Expression<Frame>::Ptr random<Frame>(std::vector<int>& ndx) {
    return frame( random<Rotation>(ndx), random<Vector>(ndx) );
}

template<> Expression<Twist>::Ptr random<Twist>(std::vector<int>& ndx) {
    return twist( random<Vector>(ndx), random<Vector>(ndx) );
}

template<> Expression<Wrench>::Ptr random<Wrench>(std::vector<int>& ndx) {
    return wrench( random<Vector>(ndx), random<Vector>(ndx) );
}

template<> Expression<Quaternion>::Ptr random<Quaternion>(std::vector<int>& ndx) {
    return quaternion( random<double>(ndx), random<Vector>(ndx) );
}

template <class T>
typename Expression<T>::Ptr testvar(int& ndx);


template<> Expression<double>::Ptr testvar<double>(int& ndx) {
    return input(ndx++);
}

template<> Expression<Vector>::Ptr testvar<Vector>(int& ndx) {
    return KDL::vector( testvar<double>(ndx), testvar<double>(ndx), testvar<double>(ndx) );
}

template<> Expression<Rotation>::Ptr testvar<Rotation>(int& ndx) {
    return rot_x( testvar<double>(ndx)) * rot_y(testvar<double>(ndx)) * rot_z(testvar<double>(ndx));
}

template<> Expression<Frame>::Ptr testvar<Frame>(int& ndx) {
    return frame( testvar<Rotation>(ndx), testvar<Vector>(ndx) );
}

template<> Expression<Twist>::Ptr testvar<Twist>(int& ndx) {
    return twist( testvar<Vector>(ndx), testvar<Vector>(ndx) );
}

template<> Expression<Wrench>::Ptr testvar<Wrench>(int& ndx) {
    return wrench( testvar<Vector>(ndx), testvar<Vector>(ndx) );
}

template<> Expression<Quaternion>::Ptr testvar<Quaternion>(int& ndx) {
    return quaternion( testvar<double>(ndx), testvar<Vector>(ndx) );
}


/******************************************************************************************************
 * CHECKING VALUE-TYPES
 ****************************************************************************************************/

bool Equal(Quaternion q1, const Quaternion& q2, double eps) {
    double c = dot(q1,q2);
    if (c<0) {
        q1=-q1;
    }
    if (fabs(q1.w-q2.w) > eps) return false;
    if (fabs(q1.vec[0]-q2.vec[0]) > eps) return false;
    if (fabs(q1.vec[1]-q2.vec[1]) > eps) return false;
    if (fabs(q1.vec[2]-q2.vec[2]) > eps) return false;
    return true;
}

::testing::AssertionResult checkVectorValues(const char* astr, const char* bstr, Vector q1, Vector q2) {
    double eps = 1E-8;
    if (debug_output) {
        std::cout << "Comparing " << q1 << std::endl;
        std::cout << "      and " << q2 << std::endl;   
    }
    if (!Equal(q1,q2,eps)) {
        return ::testing::AssertionFailure() 
           << "quaternions are not equal : \n"
           << astr << " : " << q1 << "\n"
           << bstr << " : " << q2 << "\n"
           << std::endl;
    }
    return ::testing::AssertionSuccess();
}

#define EXPECT_EQ_VECTOR( a, b ) \
    EXPECT_PRED_FORMAT2(checkVectorValues,a, b );



::testing::AssertionResult checkQuatValues(const char* astr, const char* bstr, Quaternion q1, Quaternion q2) {
    double eps = 1E-8;
    if (debug_output) {
        std::cout << "Comparing " << q1 << std::endl;
        std::cout << "      and " << q2 << std::endl;
    } 
    if (!Equal(q1,q2,eps)) {
        return ::testing::AssertionFailure() 
           << "quaternions are not equal : \n"
           << astr << " : " << q1 << "\n"
           << bstr << " : " << q2 << "\n"
           << std::endl;
    }
    return ::testing::AssertionSuccess();
}

#define EXPECT_EQ_QUAT( a, b ) \
    EXPECT_PRED_FORMAT2(checkQuatValues,a, b );

/******************************************************************************************************
 * CHECKING EXPRESSION TYPES 
 ****************************************************************************************************/


/**
 * compares the <b>value</b> of two expression graphs.  It fills in random values with setInputValue(..).
 * For each of the expressions the same value is filled in.  After that, the expressions are compared.
 */
template<class T>
::testing::AssertionResult         
EqualValues(const char* astr, const char* bstr, boost::shared_ptr< Expression<T> > a, boost::shared_ptr< Expression<T> > b) {
    double eps = 1E-4;
    int n1 = a->number_of_derivatives();
    int n2 = b->number_of_derivatives();
    // fill in arbitrary input values:
    for (int i=0;i<max(n1,n2);++i) {
        double arg;
        random(arg);
        a->setInputValue(i,arg);
        b->setInputValue(i,arg);
    }
    // compare:
    T va = a->value();
    T vb = b->value();
    if (debug_output) { 
        std::cout << "Comparing " << va << std::endl;
        std::cout << "      and " << vb << std::endl;
    }
    if (!Equal(va,vb,eps)) {
          std::stringstream os;
          os << "false because the value of the expressions is not equal  : \n"; 
          os << astr << " : \n";
          os << va << "\n";
          os << bstr << " : \n";
          os << vb;
          os << std::endl;
          return ::testing::AssertionFailure() << os.str();
    }
    for (int i=0;i<max(n1,n2);++i) {
        typedef typename AutoDiffTrait<T>::DerivType Td;
        Td ad = a->derivative(i);
        Td bd = b->derivative(i);
        if (debug_output) {
            std::cout << "      Comparing derivative : " << ad << std::endl;
            std::cout << "                       and : " << bd << std::endl;
        }
        if (!Equal(ad, bd,eps ) ) {
              std::stringstream os;
              os << "false because derivative " << i << " of the expressions is not equal  : \n"; 
              os << astr << " : \n";
              os << "value : "  << va << "\n";
              os << bstr << " : \n";
              os << "value : "  << vb << "\n";
              os << astr << " : \n";
              os << "derivative " << i << " : " << ad << "\n";
              os << bstr << " : \n";
              os << "derivative " << i << " : " << bd;
              os << std::endl;
              return ::testing::AssertionFailure() << os.str();
        }
    }
    return ::testing::AssertionSuccess();
}

/**
 * Specialization for quaternions.  To deal with double covering and derivatives:
 * To correctly interprete a comparision of derivatives, you also need to look at the
 * position-level variables.
 */
template<>
::testing::AssertionResult         
EqualValues(const char* astr, const char* bstr, boost::shared_ptr< Expression<Quaternion> > a, boost::shared_ptr< Expression<Quaternion> > b) {
    double eps = 1E-4;
    int n1 = a->number_of_derivatives();
    int n2 = b->number_of_derivatives();
    // fill in arbitrary input values:
    for (int i=0;i<max(n1,n2);++i) {
        double arg;
        random(arg);
        a->setInputValue(i,arg);
        b->setInputValue(i,arg);
    }
    // compare:
    Quaternion va = a->value();
    Quaternion vb = b->value();
    double c = dot(va,vb);
    if (c<0) {
        va = -va;        
    }
    if (debug_output) {
        std::cout << "Quaternion" << std::endl;
        std::cout << "Comparing " << va;
        if (c<0) {
            std::cout << "  (switched sign to take into account double covering of quaterions) ";
        }
        std::cout << std::endl;
        std::cout << "      and " << vb << std::endl;
    }
    if (!Equal(va,vb,eps)) {
          std::stringstream os;
          os << "false because the value of the expressions is not equal  : \n"; 
          os << astr << " : \n";
          os << va << "\n";
          os << bstr << " : \n";
          os << vb;
          os << std::endl;
          return ::testing::AssertionFailure() << os.str();
    }
    for (int i=0;i<max(n1,n2);++i) {
        Quaternion ad = a->derivative(i);
        Quaternion bd = b->derivative(i);
        if (c<0) {
            ad = -ad;
        }
        if (debug_output) {
            std::cout << "      Derivative of Quaternion" << std::endl;
            std::cout << "      Comparing " << ad;
            if (c<0) {
                std::cout << "  (switched sign to take into account double covering of quaterions) ";
            }
            std::cout << std::endl;
            std::cout << "            and " << bd << std::endl;
        }
        if (!Equal(ad, bd,eps ) ) {
              std::stringstream os;
              os << "false because derivative " << i << " of the expressions is not equal  : \n"; 
              os << astr << " : \n";
              os << ad << "\n";
              os << bstr << " : \n";
              os << bd;
              os << std::endl;
              return ::testing::AssertionFailure() << os.str();
        }
    }
    return ::testing::AssertionSuccess();
}


/**
 * An easy to use test to test the equality of the value of two expression graphs
 */
#define EXPECT_EQ_EXPR( a, b ) \
    EXPECT_PRED_FORMAT2(EqualValues, a, b );


double NumDiff(double v1, double v2, double dt) {
    return (v2-v1)/dt; 
} 

Vector NumDiff(const Vector& v1, const Vector& v2, double dt) {
    return (v2-v1)/dt; 
} 

Twist NumDiff(const Twist& v1, const Twist& v2, double dt) {
    return (v2-v1)/dt; 
} 

Wrench NumDiff(const Wrench& v1, const Wrench& v2, double dt) {
    return (v2-v1)/dt; 
} 

Quaternion NumDiff(const Quaternion& v1, const Quaternion& v2, double dt) {
    return (v2-v1)/dt; 
} 


// The KDL implementation of diff is not so good any more...
// my own based upon the quaternion routines.

Vector NumDiff(const Rotation& R1, const Rotation& R2, double dt) {
    Vector r = 2.0*logUnit( toQuat( R2*R1.Inverse() ) )/dt;
    return r;
} 

/** 
 * A predicate format to check the derivative using numerical differentiation:
 * Evaluation is for arbitrary input values, towards all relevant variables.
 *
 */
template <class T>
::testing::AssertionResult CheckNumerical(        
                                               const char* mstr,
                                               boost::shared_ptr< Expression<T> > m
                                               ) {
  double h   = 1E-9;
  double tol = 1E-4;
  typedef typename AutoDiffTrait<T>::DerivType Td;
  // generate arbitrary input for the expression:
  int n1 = m->number_of_derivatives();
  std::vector<double> arg(n1);
  for (int i=0;i<n1;++i) {
        if (simple_values) {
            arg[i]=i;
        } else {
            random(arg[i]);
        }
        m->setInputValue(i,arg[i]);
  }
  if (debug_output) {
         T value = m->value();
        std::cout << "Checking automatic differentiation and numerical differentiation" << std::endl;
        std::cout << "expression \n"; m->print(std::cout);
        std::cout << std::endl;
        std::cout << "value : " << value << std::endl;
  }
  for (int i=0;i<n1;++i) {
        //  compute numerical derivative towards variable i:
        m->setInputValue(i, arg[i]-h );
        T v1 = m->value();
        m->setInputValue(i, arg[i]+h );
        T v2 = m->value();
        Td numdiff = NumDiff(v1,v2,2*h);
        // automatic differentation:
        m->setInputValue(i,arg[i]);
        T value = m->value(); 
        Td autodiff = m->derivative(i);
        if (debug_output) {
            std::cout << "automatic derivarive towards " << i << " : " << numdiff << std::endl;
            std::cout << "numeric derivarive towards   " << i << " : " << autodiff << std::endl;
        }
        if (!Equal(numdiff,autodiff,tol)) {
            std::stringstream os;
            os << "check using numerical derivative failed for " << mstr << " : \n";  
            os << "expression \n"; m->print(os); os << "\n\n";
            os << "variables : ";
            for (int j=0;j<n1;++j) {
                os << j <<":"<<arg[j] <<" ";
            }
            os << "\n";
            os << "value                       : " << value << "\n";
            os << "derivative towards variable : " << i << "\n";
            os << "numerical differentiation   : " << numdiff << "\n"; 
            os << "automatic differentiation   : " << autodiff << "\n"; 
            return ::testing::AssertionFailure() << os.str();
        }
  }
  if (debug_output) {
        std::cout << std::endl;
  }
  return ::testing::AssertionSuccess();
}

#define CHECK_WITH_NUM( a ) \
    EXPECT_PRED_FORMAT1(CheckNumerical, a );

} // namespace KDL.

using namespace KDL;

Quaternion randomQuaternion() {
    double w,x,y,z;
    random(w);
    random(x);
    random(y);
    random(z);
    return Quaternion(w,x,y,z);
}

Vector randomVector() {
    double x,y,z;
    random(x);
    random(y);
    random(z);
    return Vector(x,y,z);
}

TEST(QuaternionValues, toRot) {
    Quaternion q1 = normalized(randomQuaternion()); 
    EXPECT_EQ_QUAT( q1, toQuat(toRot(q1)) );
    EXPECT_EQ_QUAT( -q1, toQuat(toRot(-q1)) );
}

TEST(QuaternionValues, comparisonWithMatrices) {
    Quaternion q1 = normalized(randomQuaternion()); 
    Quaternion q2 = normalized(randomQuaternion()); 
    EXPECT_EQ_QUAT( q1*q2, toQuat(toRot(q1)*toRot(q2)) );
    EXPECT_EQ_QUAT( Quaternion::Identity(), inv(q1)*q1 );
    Vector v = randomVector();
    EXPECT_EQ_VECTOR( apply(q1,v),  toRot(q1)*v );
}

TEST(QuaternionValues, normOps) {
    Quaternion q1 = randomQuaternion(); 
    EXPECT_EQ_QUAT( q1/norm(q1), normalized(q1) );
    EXPECT_EQ_QUAT( squarednorm(q1), norm(q1)*norm(q1) );
    EXPECT_EQ_QUAT( squarednorm(q1), dot(q1,q1) );
}


TEST(QuaternionValues, rotatingPoints) {
    Quaternion q1 = normalized(randomQuaternion()); 
    Vector v = randomVector();
    EXPECT_EQ_VECTOR( apply(q1,v),  (q1*(v*inv(q1))).vec );
    EXPECT_EQ_VECTOR( apply(q1,v),  (q1*(Quaternion(v)*inv(q1))).vec );
    EXPECT_NEAR( 0.0,  (q1*(Quaternion(v)*inv(q1))).w, QUAT_EPS );
    EXPECT_EQ_VECTOR( apply(q1,v),  toRot(q1)*v );
    EXPECT_EQ_VECTOR( apply(inv(q1),v),  (toRot(q1)).Inverse()*v );
    EXPECT_EQ_VECTOR( v,  apply(inv(q1),apply(q1,v)) );
}


TEST(QuaternionValues, testInverse) {
    Quaternion q1 = normalized(randomQuaternion()); 
    EXPECT_EQ_QUAT( Quaternion::Identity(), inv(q1)*q1 );
    EXPECT_EQ_QUAT( Quaternion::Identity(), q1*inv(q1) );
    EXPECT_EQ_QUAT( inv(q1), conj(q1)/dot(q1,q1) );
}

TEST(QuaternionValues, testExpLog) {
    Quaternion q1 = normalized(randomQuaternion()); 
    EXPECT_EQ_QUAT( q1, exp(log(q1)) );
    EXPECT_EQ_QUAT( q1, exp(logUnit(q1)) );
    Quaternion q2 = randomQuaternion(); 
    q2.w = 0; // a pure quaternion
    EXPECT_EQ_VECTOR( q2.vec, logUnit(exp(q2)) );
    EXPECT_EQ_QUAT( q2, log(exp(q2)) );
    EXPECT_EQ_QUAT( q2*exp(q2), exp(q2)*q2 );
}

TEST(QuaternionValues, operatorProperties) {
    Quaternion q1 = normalized(randomQuaternion()); 
    Quaternion q2 = normalized(randomQuaternion()); 
    Quaternion q3 = normalized(randomQuaternion()); 
    EXPECT_EQ_QUAT( q1*q2 + q1*q3  ,  q1*(q2+q3) );
    EXPECT_NEAR( dot(q1,q2) + dot(q1,q3)  ,  dot(q1,q2+q3), QUAT_EPS );
    EXPECT_EQ_QUAT( Quaternion::Identity(), inv(q1)*q1 );
}

TEST(QuaternionValues, axisAngle) {
    Quaternion q1 = normalized(randomQuaternion()); 
    Vector axis;
    double angle;
    axis_angle(q1,axis,angle);
    EXPECT_EQ_QUAT( q1  ,  toQuat(axis,angle) );
    EXPECT_EQ_QUAT( q1  ,  toQuat(axis*angle) );
    EXPECT_EQ_QUAT( q1  ,  exp(axis*angle/2) );
    EXPECT_EQ_VECTOR( axis*angle/2, logUnit(q1) );
    EXPECT_EQ_VECTOR( axis/norm(axis), KDL::axis(q1)/norm(KDL::axis(q1)) );
}

TEST(QuaternionValues, norm) {
    Quaternion q1 = randomQuaternion(); 
    Quaternion q2 = randomQuaternion(); 
    q2.w=0;
    EXPECT_EQ_QUAT( (q2*q2), squarednorm(q2)*Quaternion(1,0,0,0)) ;
    EXPECT_EQ_QUAT( squarednorm(q1), norm(q1)*norm(q1) );
    EXPECT_EQ_QUAT( squarednorm(q1), dot(q1,q1) );
}

TEST(QuaternionValues, operatorRelations) {
    Quaternion q1 = randomQuaternion(); 
    Quaternion q2 = randomQuaternion(); 
    q2.w = 0;
    EXPECT_EQ_QUAT( Quaternion::Zero(), q1-q1 );
    EXPECT_EQ_QUAT( 2*q1, q1+q1 );
    EXPECT_EQ_QUAT( q1, q1+q1+(-q1) );
    EXPECT_EQ_QUAT( -q1, (-1)*q1 );
    EXPECT_EQ_QUAT( -q1, q1*(-1) );
    EXPECT_EQ_QUAT( q1*2, q1+q1 );
    EXPECT_EQ_QUAT( q1*q2.vec, q1*Quaternion(q2.vec) );
    EXPECT_EQ_QUAT( q2.vec*q1, Quaternion(q2.vec)*q1 );
}

TEST(QuaternionValues, pow) {
    Quaternion q1 = normalized(randomQuaternion()); 
    Quaternion q2 = randomQuaternion(); 
    EXPECT_EQ_QUAT( pow(q2,2), q2*q2 );
    EXPECT_EQ_QUAT( powUnit(q1,2), q1*q1 );
    EXPECT_EQ_QUAT( pow(q2,3), q2*q2*q2 );
    EXPECT_EQ_QUAT( powUnit(q1,3), q1*q1*q1 );
    EXPECT_EQ_QUAT( pow(q2,-2), inv(q2*q2) );
    EXPECT_EQ_QUAT( powUnit(q1,3), q1*q1*q1 );
}

TEST(QuaternionValues, diff) {
    Vector a(1,2,3);
    a.Normalize();
    Quaternion q1 = toQuat( a, 20.0/180.0*PI); 
    Quaternion q2 = toQuat( a, 190.0/180.0*PI); 
    Vector d = diff(q1,q2);
    EXPECT_EQ_VECTOR( d, a*170.0/180.0*PI);
}



TEST(QuaternionCompareValueExpr,exp) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    EXPECT_EQ_QUAT( exp(q1)->value(), exp(q1->value()) );
} 

TEST(QuaternionCompareValueExpr,log) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    EXPECT_EQ_QUAT( log(q1)->value(), log(q1->value()) );
} 

TEST(QuaternionCompareValueExpr,logUnit) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized(random<Quaternion>(ndx));
    EXPECT_EQ_VECTOR( logUnit(q1)->value(), logUnit(q1->value()) );
} 

TEST(QuaternionCompareValueExpr,toRot) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized(random<Quaternion>(ndx));
    EXPECT_EQ( toRot(q1)->value(), toRot(q1->value()) );
} 

TEST(QuaternionCompareValueExpr,toQuat) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized(random<Quaternion>(ndx));
    EXPECT_EQ_QUAT( toQuat(toRot(q1))->value(), toQuat(toRot(q1->value())) );
} 

TEST(QuaternionCompareValueExpr,apply) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized(random<Quaternion>(ndx));
    std::vector<int> ndx2; 
    ndx2.push_back(5);
    ndx2.push_back(6);
    ndx2.push_back(3);
    Expression<Vector>::Ptr v = normalized(random<Vector>(ndx2));
    Vector v1 = apply(q1,v)->value();
    Vector v2 = apply(q1->value(), v->value() );
    EXPECT_EQ_VECTOR( v1, v2);
    EXPECT_EQ_VECTOR( apply(inv(q1),v)->value(), apply(inv(q1->value()), v->value() ));
} 

TEST(QuaternionCompareValueExpr,normOps) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    EXPECT_EQ_QUAT( (q1/norm(q1))->value(), q1->value() / norm( q1->value() ) );
    EXPECT_EQ_QUAT( squarednorm(q1)->value(), squarednorm(q1->value()) );
    EXPECT_EQ_QUAT( norm(q1)->value(), norm(q1->value()) );
    EXPECT_EQ_QUAT( dot(q1,q1)->value(), dot(q1->value(),q1->value()) );
} 
TEST(QuaternionCompareValueExpr,inverse) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    EXPECT_EQ_QUAT( inv(q1)->value(), inv(q1->value()) );
} 
TEST(QuaternionCompareValueExpr,pow) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    EXPECT_EQ_QUAT( pow(q1, Constant(2.0))->value(), q1->value()*q1->value() );
} 
TEST(QuaternionCompareValueExpr,operators) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    std::vector<int> ndx2; 
    ndx2.push_back(5);
    ndx2.push_back(6);
    ndx2.push_back(7);
    ndx2.push_back(8);
    Expression<Quaternion>::Ptr q2 = random<Quaternion>(ndx2);
 
    EXPECT_EQ_QUAT( (q1*q2)->value(), q1->value() * q2->value() );
    EXPECT_EQ_QUAT( (q1+q2)->value(), q1->value() + q2->value() );
    EXPECT_EQ_QUAT( (q1+Constant(-1.0)*q2)->value(), q1->value() + (-1)*q2->value() );
    EXPECT_EQ_QUAT( (q1-q2)->value(), q1->value()-q2->value() );
    EXPECT_EQ_QUAT( (-q2)->value(), -q2->value() );
} 
/*
TEST(QuaternionCompareValueExpr,axis) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    Vector axis;
    double angle;
    axis_angle(q1->value(),axis,angle);
    Expression<Vector>::Ptr a = logUnit(q1)*2; 
    EXPECT_EQ( a->value(), axis*angle  );
} 
*/
TEST(QuaternionCompareValueExpr,diff) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    std::vector<int> ndx2; 
    ndx2.push_back(5);
    ndx2.push_back(6);
    ndx2.push_back(7);
    ndx2.push_back(8);
    Expression<Quaternion>::Ptr q2 = random<Quaternion>(ndx2);
 
    EXPECT_EQ( (diff(q1,q2))->value(), diff(q1->value() , q2->value()) );
} 

TEST(DoubleNumDiff, Scalars) {
        int ndx=1;
        Expression<double>::Ptr a =  testvar<double>(ndx);
        Expression<double>::Ptr b = testvar<double>(ndx);
        CHECK_WITH_NUM( -a );
        CHECK_WITH_NUM( sin(a) );
        CHECK_WITH_NUM( cos(a) );
        CHECK_WITH_NUM( tan(a) );
        CHECK_WITH_NUM( exp(a) );
        CHECK_WITH_NUM( asin(a) );
        CHECK_WITH_NUM( acos(a*a+Constant(0.01)) );
        CHECK_WITH_NUM( exp(a) );
        CHECK_WITH_NUM( log(a*a) );
        CHECK_WITH_NUM( sqr(a) );
        CHECK_WITH_NUM( sqrt(a*a+Constant(0.001)) );
        CHECK_WITH_NUM( abs(a) );
        CHECK_WITH_NUM( a + b );
        CHECK_WITH_NUM( a - b );
        CHECK_WITH_NUM( a * b );
        CHECK_WITH_NUM( a / b );
        CHECK_WITH_NUM( atan2(a,b)  );
}



TEST(QuaternionNumDiff,operators) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = testvar<Quaternion>(ndx);
/*    Quaternion q1v(3,2,1,0);
    Quaternion q2v(7,6,5,4);
    Quaternion dq1v(0,0,0,0);
    Quaternion dq2v(0,0,0,1);
    std::cout << "q1 : " << q1v << std::endl; 
    std::cout << "q2 : " << q2v << std::endl; 
    std::cout << "q1*q2 : " << q1v * q2v << std::endl;
    std::cout << "dq1*q2 + q1*dq2 : " << dq1v*q2v + q1v*dq2v << std::endl;

    Expression<Quaternion>::Ptr qr = q1*q2;
    qr->setInputValue(0,0.0);
    qr->setInputValue(1,1.0);
    qr->setInputValue(2,2.0);
    qr->setInputValue(3,3.0);
    qr->setInputValue(4,4.0);
    qr->setInputValue(5,5.0);
    qr->setInputValue(6,6.0);
    qr->setInputValue(7,7.0);
    std::cout << "q1*q2 expr value: " << qr->value() << std::endl;
    std::cout << "q1*q2 expr derivative 0: " << qr->derivative(0) << std::endl;
    std::cout << "q1*q2 expr derivative 1: " << qr->derivative(1) << std::endl;
    std::cout << "q1*q2 expr derivative 2: " << qr->derivative(2) << std::endl;
    std::cout << "q1*q2 expr derivative 3: " << qr->derivative(3) << std::endl;
    std::cout << "q1*q2 expr derivative 4: " << qr->derivative(4) << std::endl;
    std::cout << "q1*q2 expr derivative 5: " << qr->derivative(5) << std::endl;
    std::cout << "q1*q2 expr derivative 6: " << qr->derivative(6) << std::endl;
    std::cout << "q1*q2 expr derivative 7: " << qr->derivative(7) << std::endl;
    double h=1E-9;
    double v=4.0;
    qr->setInputValue(4,v+h);
    Quaternion v1 = qr->value();
    std::cout << "q1*q2 numdiff expr v+h " << v1 << std::endl;
    qr->setInputValue(4,v-h);
    Quaternion v2 = qr->value();
    std::cout << "q1*q2 numdiff expr v-h " << v2 << std::endl;
    std::cout << "q1*q2 numdiff expr " << (v2-v1)/2.0/h << std::endl;
    */
    CHECK_WITH_NUM( q1 );
    CHECK_WITH_NUM( q2 );
    CHECK_WITH_NUM( -q1 );
    CHECK_WITH_NUM( q1+q2 );
    CHECK_WITH_NUM( q1+q2 );
    CHECK_WITH_NUM( q1-q2 );
    CHECK_WITH_NUM( q1*q2 );
    CHECK_WITH_NUM( -q1 );
}
TEST(QuaternionNumDiff,dot_conj) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = testvar<Quaternion>(ndx);
    CHECK_WITH_NUM( conj(q1) );
    CHECK_WITH_NUM( dot(q1,q2) );
    CHECK_WITH_NUM( dot(q1,conj(q1))-norm(q1)*norm(q1) );
}

TEST(QuaternionNumDiff, toQuat) {
    int ndx=0;
    Expression<double>::Ptr r = testvar<double>(ndx);
    Expression<double>::Ptr p = testvar<double>(ndx);
    Expression<double>::Ptr y = testvar<double>(ndx);
    Expression<Rotation>::Ptr R = rot_z(y)*rot_y(p)*rot_x(r);
    CHECK_WITH_NUM( toQuat(R) );
}

TEST(QuaternionNumDiff, toRot) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = normalized(testvar<Quaternion>(ndx));
    Expression<double>::Ptr r = testvar<double>(ndx);
    Expression<double>::Ptr p = testvar<double>(ndx);
    Expression<double>::Ptr y = testvar<double>(ndx);
    Expression<Rotation>::Ptr R = rot_z(y)*rot_y(p)*rot_x(r);
    CHECK_WITH_NUM( toRot(q1) );
    CHECK_WITH_NUM( toQuat(toRot(q1)) );
}


TEST(QuaternionNumDiff, exp) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<Vector>::Ptr v = testvar<Vector>(ndx);
    CHECK_WITH_NUM(exp(q1));
    CHECK_WITH_NUM(exp(v));
}

TEST(QuaternionNumDiff, normalize) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    CHECK_WITH_NUM(normalized(q1));
    CHECK_WITH_NUM(norm(q1));
    CHECK_WITH_NUM(squarednorm(q1));
}

TEST(QuaternionNumDiff, inv) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    CHECK_WITH_NUM(inv(q1));
}

TEST(QuaternionNumDiff, log) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    CHECK_WITH_NUM(logUnit(normalized(q1)));
    CHECK_WITH_NUM(log(q1));
}

TEST(QuaternionNumDiff, pow) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<double>::Ptr s = testvar<double>(ndx);
    CHECK_WITH_NUM(pow(q1, s));
    CHECK_WITH_NUM(powUnit(normalized(q1), s));
}


TEST(QuaternionNumDiff, apply) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<Vector>::Ptr v = testvar<Vector>(ndx);
    CHECK_WITH_NUM(apply(q1, v));
    CHECK_WITH_NUM(q1*v*conj(q1));
}

TEST(QuaternionNumDiff, diff) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = normalized(testvar<Quaternion>(ndx));
    Expression<Quaternion>::Ptr q2 = normalized(testvar<Quaternion>(ndx));
    CHECK_WITH_NUM(diff(q1, q2));
}
TEST(QuaternionNumDiff, toQuat2) {
    int ndx=0;
    Expression<Vector>::Ptr omega = testvar<Vector>(ndx);
    CHECK_WITH_NUM(toQuat(omega));
}
TEST(QuaternionNumDiff, toQuat3) {
    int ndx=0;
    Expression<Vector>::Ptr axis = testvar<Vector>(ndx);
    Expression<double>::Ptr angle = testvar<double>(ndx);
    CHECK_WITH_NUM(toQuat(axis,angle));
}
TEST(QuaternionNumDiff, axisAngle) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = normalized(testvar<Quaternion>(ndx));
    CHECK_WITH_NUM(axisAngle(q1));
}
TEST(QuaternionNumDiff, slerp) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = testvar<Quaternion>(ndx);
    Expression<double>::Ptr s = testvar<double>(ndx);
    CHECK_WITH_NUM(slerp(q1,q2,s));
}
TEST(QuaternionNumDiff, slerpUnit) {
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = normalized(testvar<Quaternion>(ndx));
    Expression<Quaternion>::Ptr q2 = normalized(testvar<Quaternion>(ndx));
    Expression<double>::Ptr s = testvar<double>(ndx);
    CHECK_WITH_NUM(slerpUnit(q1,q2,s));
}


TEST(QuaternionNumDiff, quaternion) {
    int ndx=0;
    Expression<Vector>::Ptr v = testvar<Vector>(ndx);
    Expression<double>::Ptr s = testvar<double>(ndx);
    CHECK_WITH_NUM(quaternion(v,s));
    CHECK_WITH_NUM(quaternion(s,v));
}
/*
TEST(QuaternionNumDiff,special_values) {
    int ndx=1;
    Expression<Quaternion>::Ptr q = testvar<Quaternion>(ndx);
    double h=1E-8;
    q->setInputVar(1,0.0); q->setInputVar(2,0.0); q->setInputVar(3,0.0); q->setInputVar(4,1.0);
    
}
*/


TEST(DoubleExpr, operators) {
        int ndx=0;
        Expression<double>::Ptr a =  testvar<double>(ndx);
        Expression<double>::Ptr b = testvar<double>(ndx);
        EXPECT_EQ_EXPR( a*b, b*a );
        EXPECT_EQ_EXPR( a+b, b+a );
        EXPECT_EQ_EXPR( a+Constant(0.0), a );
        EXPECT_EQ_EXPR( a*Constant(1.0), a );
        EXPECT_EQ_EXPR( a-b, -(b-a) );
        EXPECT_EQ_EXPR( a/b, a*(Constant(1.0)/b) );
        EXPECT_EQ_EXPR( a/b, (Constant(1.0)/b)*a );
        EXPECT_EQ_EXPR( a*a/Constant(10.0), acos(cos(a*a/Constant(10.0))) );
        EXPECT_EQ_EXPR( a/Constant(10.0), asin(sin(a/Constant(10.0))) );
        EXPECT_EQ_EXPR( a/Constant(10.0), atan(tan(a/Constant(10.0))) );
}

TEST(DoubleExpr, transcendental_functions) {
        int ndx=0;
        Expression<double>::Ptr a =  testvar<double>(ndx);
        Expression<double>::Ptr b = testvar<double>(ndx);
 
        EXPECT_EQ_EXPR( a*a, sqr(sqrt(a*a)) );
        EXPECT_EQ_EXPR( abs(a), sqrt(sqr(a)) );
        EXPECT_EQ_EXPR( a, log(exp(a)) );
        EXPECT_EQ_EXPR( a, log(exp(a)) );
        EXPECT_EQ_EXPR( atan2(sin(a),cos(a)) , a );
}

TEST(QuaternionExpr, mult_sum) {
    /*
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = random<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q3 = random<Quaternion>(ndx);
    */
    int ndx=0;
    Expression<Quaternion>::Ptr q1 = testvar<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = testvar<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q3 = testvar<Quaternion>(ndx);
 
    //EXPECT_EQ_EXPR( q1, toQuat(toRot(q1)) );
    EXPECT_EQ_EXPR( Constant(2.0)*q1, q1+q1 );
    EXPECT_EQ_EXPR( Constant(1.0)*q1, q1+q1-q1 );
    EXPECT_EQ_EXPR( Constant(1.0)*q1, q1+q1+Constant(-1.0)*q1 );
    EXPECT_EQ_EXPR( q1*(q2*q3), (q1*q2)*q3 );
}


TEST(QuaternionExpr, normOps) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = random<Quaternion>(ndx);
    EXPECT_EQ_EXPR( normalized(q1) , q1 / norm(q1) );
    EXPECT_EQ_EXPR( squarednorm(q1), dot(q1,q1 ));
    EXPECT_EQ_EXPR( squarednorm(q1), norm(q1)*norm(q1) );
    EXPECT_EQ_EXPR( Constant(1.0), norm( normalized(q1) * normalized(q1) ));
    EXPECT_EQ_EXPR( Constant(1.0), norm( normalized(q1)  ));
}

TEST(QuaternionExpr, toRot) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized( random<Quaternion>(ndx) );
    EXPECT_EQ_EXPR( q1, toQuat(toRot(q1)) );
    EXPECT_EQ_EXPR( -q1, -toQuat(toRot(q1)) );
}

TEST(QuaternionExpr, comparisonWithMatrices) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized( random<Quaternion>(ndx) );
    Expression<Quaternion>::Ptr q2 = normalized( random<Quaternion>(ndx) );
    std::vector<int> ndx2; 
    ndx2.push_back(5);
    ndx2.push_back(6);
    ndx2.push_back(2);
    Expression<Vector>::Ptr      p = random<Vector>(ndx2);
    EXPECT_EQ_EXPR( apply(q1,p), toRot(q1)*p );
    EXPECT_EQ_EXPR( q1*q2, toQuat( toRot(q1) * toRot(q2) ) );
    EXPECT_EQ_EXPR( Constant<Quaternion>(Quaternion::Identity()), inv(q1)*q1 );
    EXPECT_EQ_EXPR( Constant<Rotation>(Rotation::Identity()), toRot(inv(q1))*toRot(q1) );
}

TEST(QuaternionExpr, rotatingPoints) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized( random<Quaternion>(ndx) );
    std::vector<int> ndx2; 
    ndx2.push_back(5);
    ndx2.push_back(6);
    ndx2.push_back(2);
    Expression<Vector>::Ptr      v = random<Vector>(ndx2);
 
    EXPECT_EQ_EXPR( apply(q1,v),  vec(q1*(v*inv(q1))) );
    EXPECT_EQ_EXPR( apply(q1,v),  vec(q1*(v*inv(q1))) );
    EXPECT_EQ_EXPR( Constant(0.0),  w(q1*(v*inv(q1))) );
    EXPECT_EQ_EXPR( apply(q1,v),  toRot(q1)*v );
    EXPECT_EQ_EXPR( apply(inv(q1),v),  inv(toRot(q1))*v );
    EXPECT_EQ_EXPR( v,  apply(inv(q1),apply(q1,v)) );
}


TEST(QuaternionExpr, testInverse) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized( random<Quaternion>(ndx) );
    EXPECT_EQ_EXPR( Constant(Quaternion::Identity()), inv(q1)*q1 );
    EXPECT_EQ_EXPR( Constant(Quaternion::Identity()), q1*inv(q1) );
    EXPECT_EQ_EXPR( inv(q1), conj(q1)/dot(q1,q1) );
}

TEST(QuaternionExpr, testExpLog) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized( random<Quaternion>(ndx) );
    Expression<Quaternion>::Ptr q2 = quaternion(Constant(0.0),normalized( random<Vector>(ndx) ));
    EXPECT_EQ_EXPR( q1, exp(log(q1)) );
    EXPECT_EQ_EXPR( q1, exp(logUnit(q1)) );
    EXPECT_EQ_EXPR( vec(q2), logUnit(exp(q2)) );
    EXPECT_EQ_EXPR( q2, log(exp(q2)) );
    EXPECT_EQ_EXPR( q2*exp(q2), exp(q2)*q2 );
}

TEST(QuaternionExpr, operatorProperties) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized( random<Quaternion>(ndx) );
    Expression<Quaternion>::Ptr q2 = normalized( random<Quaternion>(ndx) );
    Expression<Quaternion>::Ptr q3 = normalized( random<Quaternion>(ndx) );
 
    EXPECT_EQ_EXPR( q1*q2 + q1*q3  ,  q1*(q2+q3) );
    EXPECT_EQ_EXPR( dot(q1,q2) + dot(q1,q3)  ,  dot(q1,q2+q3) );
    EXPECT_EQ_EXPR( Constant(Quaternion::Identity()), inv(q1)*q1 );
    EXPECT_EQ_EXPR( Constant(Quaternion::Identity()), inv(q1)*q1 );
}


TEST(QuaternionExpr, mult_keeps_norm) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    Expression<Quaternion>::Ptr q2 = random<Quaternion>(ndx);
 
    EXPECT_EQ_EXPR( norm(q1*q2) ,  norm(q1)*norm(q2) );
}


TEST(QuaternionExpr, pow) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    ndx.push_back(4);
    Expression<Quaternion>::Ptr q1 = normalized(random<Quaternion>(ndx));
    Expression<Quaternion>::Ptr q2 = random<Quaternion>(ndx);
 
    EXPECT_EQ_EXPR( exp(log(q2)*Constant(2.0)), q2*q2 );
    EXPECT_EQ_EXPR( pow(q2,Constant(2.0)), q2*q2 );
    EXPECT_EQ_EXPR( powUnit(q1,Constant(2.0)), q1*q1 );
    EXPECT_EQ_EXPR( pow(q2,Constant(3.0)), q2*q2*q2 );
    EXPECT_EQ_EXPR( powUnit(q1,Constant(3.0)), q1*q1*q1 );
    EXPECT_EQ_EXPR( pow(q2,Constant(-2.0)), inv(q2*q2) );
    EXPECT_EQ_EXPR( powUnit(q1,Constant(3.0)), q1*q1*q1 );
}

TEST(QuaternionExpr, diff) {
    std::vector<int> ndx; 
    ndx.push_back(1);
    ndx.push_back(2);
    ndx.push_back(3);
    Expression<Vector>::Ptr a = normalized(random<Vector>(ndx));
    Expression<Quaternion>::Ptr q1 = toQuat( a, Constant(20.0/180.0*PI)); 
    Expression<Quaternion>::Ptr q2 = toQuat( a, Constant(190.0/180.0*PI)); 
    Expression<Vector>::Ptr d = diff(q1,q2);
    EXPECT_EQ_EXPR( d, a*Constant(170.0/180.0*PI));
}




/*
TEST(QuaternionDeriv, conj_inv) {
        std::vector<int> ndx; 
        ndx.push_back(1);
        ndx.push_back(2);
        ndx.push_back(3);
        ndx.push_back(4);
        Expression<Quaternion>::Ptr q1 = random<Quaternion>(ndx);
    
        CHECK_WITH_NUM( conj(q1) );
}*/


/**        
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
 
        EXPECT_EQ_EXPR( a*b, -b*a );
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
**/

// Run all the tests that were declared with TEST()
int main(int argc, char **argv){
    debug_output=0;
    simple_values=0;
    testing::InitGoogleTest(&argc, argv);
    std::cout << "usage:"<< std::endl;
    std::cout << "   output    : print additional debugging info\n"
              << "   simpleval : choose simple test values, otherwise random\n"<<std::endl;;
    for (int i=1;i<argc;++i) {
        if (strcmp(argv[i],"output")==0) {
            std::cout << "additional debugging information will be printed" << std::endl;
            debug_output=1;
        }
        if (strcmp(argv[i],"simpleval")==0) {
            std::cout << "additional debugging information will be printed" << std::endl;
            simple_values=1;
        }
        
    }
    return RUN_ALL_TESTS();
}


