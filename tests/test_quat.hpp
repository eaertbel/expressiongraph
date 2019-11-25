#ifndef EXPRESSIONGRAPH_TEST_QUAT_HPP
#define EXPRESSIONGRAPH_TEST_QUAT_HPP

#include <kdl/expressiontree_quaternion.hpp>
#include <kdl/expressiontree.hpp>
#include <gtest/gtest.h>
#include <string.h>


namespace KDL {

extern bool debug_output;    // output additional debug values
extern bool simple_values;   // choose simple values, i.e. var. i has value i, otherwise random values are used.
extern bool check_deriv_expr;   // if 1, there will be testing of derivative expressions.

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
    if (fabs(q1.w-q2.w) > eps) return false;
    if (fabs(q1.vec[0]-q2.vec[0]) > eps) return false;
    if (fabs(q1.vec[1]-q2.vec[1]) > eps) return false;
    if (fabs(q1.vec[2]-q2.vec[2]) > eps) return false;
    return true;
}

bool EqualUnit(Quaternion q1, const Quaternion& q2, double eps) {
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

::testing::AssertionResult checkUnitQuatValues(const char* astr, const char* bstr, Quaternion q1, Quaternion q2) {
    double eps = 1E-8;
    if (debug_output) {
        std::cout << "Comparing " << q1 << std::endl;
        std::cout << "      and " << q2 << std::endl;
    } 
    if (!EqualUnit(q1,q2,eps)) {
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

#define EXPECT_EQ_UNITQUAT( a, b ) \
    EXPECT_PRED_FORMAT2(checkUnitQuatValues,a, b );


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
 * An easy to use test to test the equality of the value of two expression graphs
 */
#define EXPECT_EQ_EXPR( a, b ) \
    EXPECT_PRED_FORMAT2(EqualValues, a, b );


/**
 * Specialization for quaternions.  To deal with double covering and derivatives:
 * To correctly interprete a comparision of derivatives, you also need to look at the
 * position-level variables.
 */
::testing::AssertionResult         
EqualUnitValues(const char* astr, const char* bstr, boost::shared_ptr< Expression<Quaternion> > a, boost::shared_ptr< Expression<Quaternion> > b) {
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

#define EXPECT_EQ_UNITQUAT_EXPR( a, b ) \
    EXPECT_PRED_FORMAT2(EqualUnitValues, a, b );



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
        // value of the derivative expression
        if (check_deriv_expr) {
            typename Expression<Td>::Ptr deriv_expr = m->derivativeExpression(i);
            if (!deriv_expr) {
                std::stringstream os;
                os << "derivative expression not implemented for " << mstr << " : \n";  
                os << "expression \n"; m->print(os); os << "\n\n";
                return ::testing::AssertionFailure() << os.str();
            } else {
               // still retains the input values from the originating expression
                Td deriv_expr_value = deriv_expr->value();
                if (debug_output) {
                    std::cout << "value of derivative expr.    " << i << " : " << deriv_expr_value << std::endl;
                }
                if (!Equal(autodiff, deriv_expr_value, tol) ) {
                    std::stringstream os;
                    os << "check using the dervative value and the value of the derivative expression failed for " 
                       << mstr << " : \n";  
                    os << "expression \n"; m->print(os); os << "\n\n";
                    os << "variables : ";
                    for (int j=0;j<n1;++j) {
                        os << j <<":"<<arg[j] <<" ";
                    }
                    os << "\n";
                    os << "derivative towards variable : " << i << "\n";
                    os << "value of the derivative expr : " << deriv_expr_value << "\n"; 
                    os << "automatic differentiation   : " << autodiff << "\n"; 
                    return ::testing::AssertionFailure() << os.str();
                }
            }
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

#endif
