#ifndef EXPRESSIONGRAPH_TEST_HPP
#define EXPRESSIONGRAPH_TEST_HPP

#include <kdl/expressiontree.hpp>
#include <gtest/gtest.h>

namespace KDL {


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
            random(b);
            result = result + Constant(a)*input(ndx[i])+Constant(b);
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
        Td ad;
        Td bd;
        if (i<n1) {
            ad = a->derivative(i);
        } else {
            ad = AutoDiffTrait<T>::zeroDerivative();
        }
        if (i<n2) {
            bd = b->derivative(i);
        } else {
            bd = AutoDiffTrait<T>::zeroDerivative();
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
 * A predicate format to give back useful information in case an equal value test fails:
 *
template <class T>
::testing::AssertionResult EqualValues(        
                                               const char* mstr,
                                               const char* nstr,
                                               boost::shared_ptr< Expression<T> > m,
                                               boost::shared_ptr< Expression<T> > n) {
  if (EqualValue<T>(m,n,1E-4))
    return ::testing::AssertionSuccess();

  std::stringstream os;
  os << "false because the following are not equal : \n"; 
  os << mstr << " : \n";
  display<T>(os,m);
  os << nstr << " : \n";
  display<T>(os,n);
  os << std::endl;
  return ::testing::AssertionFailure() << os.str();
}
*/


/**
 * An easy to use test to test the equality of the value of two expression graphs
 */
#define EXPECT_EQ_VALUES( a, b ) \
    EXPECT_PRED_FORMAT2(EqualValues, a, b );





// the rotational velocity is expressed in the absolute frame, i.e.  \f$ \dot{R} = \omega \prod R \f$
template<typename ResultType>
typename AutoDiffTrait<ResultType>::DerivType
inline numerical_derivative_rot( typename Expression<ResultType>::Ptr expr, int towards_var, int component,const Rotation& value, double h=1E-7) {
    ResultType a,b;
    Rotation val;
    if (component==0) {
        val = Rotation::RotX(-h);
    } else if (component==1) {  
        val = Rotation::RotY(-h);
    } else {
        val = Rotation::RotZ(-h);
    }
    val = val*value;
    expr->setInputValue(towards_var,val);
    a = expr->value();

    if (component==0) {
        val = Rotation::RotX(h);
    } else if (component==1) {  
        val = Rotation::RotY(h);
    } else {
        val = Rotation::RotZ(h);
    }
    val = val*value;
 
    expr->setInputValue(towards_var,val);
    b = expr->value();

    return KDL::diff(a,b,1.0)/(2.0*h);
}


template <class T>
void setArbitraryInput(typename Expression<T>::Ptr e) {
    typedef std::set<int> Set;
    Set scalardep;
    Set rotdep;
    e->getScalarDependencies(scalardep);
    for (Set::iterator it=scalardep.begin();it!=scalardep.end();++it) {
        double arg;
        random(arg);
        e->setInputValue(*it, arg);
    }
    e->getRotDependencies(rotdep);
    for (Set::iterator it=rotdep.begin();it!=rotdep.end();++it) {
        Rotation arg;
        random(arg);
        e->setInputValue(*it, arg);
    }
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
        random(arg[i]);
        m->setInputValue(i,arg[i]);
  }
  m->value(); // needs to be called before using derivative(...)
  for (int i=0;i<n1;++i) {
        // check derivative towards variable i:
        Td d1 = numerical_derivative<T>(m,i,arg[i], h );
        std::cout << d1 << std::endl;
        Td d2 = m->derivative(i);
        std::cout << d2 << std::endl;
        if (!Equal(d1,d2,tol)) {
            std::stringstream os;
            os << "check using numerical derivative failed for " << mstr << " : \n";  
            os << "derivative towards variable " << i << "\n";
            os << "numerical differentiation:\n" << d1 << "\n"; 
            os << "automatic differentiation:\n" << d2 << "\n"; 
            return ::testing::AssertionFailure() << os.str();
        }
  }
  return ::testing::AssertionSuccess();
}

#define CHECK_WITH_NUM( a ) \
    EXPECT_PRED_FORMAT1(CheckNumerical, a );

/**
 * computing the numerical derivative for an expression tree
 * ( not using the derivative() function )
 * \param expr [in] expression to compute the numerical derivative.
 * \param towards_var [in]  variable number of the variable towards the derivative will be taken.
 * \param value [in]  value for the towards_var-th variable. 
 * \param h [in] interval over which to take the numerical derivative. 
 */
template<typename ResultType>
typename AutoDiffTrait<ResultType>::DerivType
inline numerical_derivative_print( typename Expression<ResultType>::Ptr expr, int towards_var, double value, double h=1E-7) {
    ResultType a,b;
    double val;
    val = value - h;
    expr->setInputValue(towards_var,val);
    a = expr->value();

    val = value + h;
    expr->setInputValue(towards_var,val);
    b = expr->value();
    std::cout << "f(x-h) " << a  << std::endl;
    std::cout << "f(x+h) " << b  << std::endl;
    std::cout << "diff " << KDL::diff(a,b,1.0)/2/h << std::endl;
    return KDL::diff(a,b,1.0)/(2.0*h);
}


/** 
 * A predicate format to check the derivative using numerical differentiation:
 * Evaluation is for arbitrary input values, towards all relevant variables.
 */
template <class T>
::testing::AssertionResult CheckNumericalRot(        
                                               const char* mstr,
                                               boost::shared_ptr< Expression<T> > e
                                               ) {
    typedef typename AutoDiffTrait<T>::DerivType Td;
    double h = 1E-8;
    double tol = 1E-4;
 
    typedef std::set<int> Set;
    Set scalardep;
    Set rotdep;
    e->getScalarDependencies(scalardep);
    std::vector<int>      scalar_ndx;
    std::vector<double>   scalar_value;
    std::vector<int>      rot_ndx;
    std::vector<Rotation> rot_value;
    for (Set::iterator it=scalardep.begin();it!=scalardep.end();++it) {
        double arg;
        random(arg);
        scalar_ndx.push_back(*it);
        scalar_value.push_back(arg);
        e->setInputValue(*it, arg);
    }
    e->getRotDependencies(rotdep);
    for (Set::iterator it=rotdep.begin();it!=rotdep.end();++it) {
        Rotation arg;
        random(arg);
        rot_ndx.push_back(*it);
        rot_value.push_back(arg);
        e->setInputValue(*it, arg);
    }
    e->value(); // needs to be called before using derivative(...)
    for (size_t i=0;i<scalar_ndx.size();++i) {
        // check derivative towards variable i:
        Td d1 = numerical_derivative<T>(e,scalar_ndx[i],scalar_value[i], h );
        Td d2 = e->derivative(scalar_ndx[i]);
        //std::cout << d1 << " compared to " << d2 << std::endl;
        if (!Equal(d1,d2,tol)) {
            std::stringstream os;
            os << "check using numerical derivative failed for " << mstr << " : \n";  
            os << "expression "; e->print(os); os << "\n";
            for (size_t j=0;j<scalar_ndx.size();++j) {
                os << "   var " << scalar_ndx[j] << " value=" << scalar_value[j] << "\n";
            }
            os << "value " << e->value() << "\n";
            Td d3 = numerical_derivative_print<T>(e,scalar_ndx[i],scalar_value[i], h );
            os << "derivative towards variable " << scalar_ndx[i] << "\n";
            os << "numerical differentiation:\n" << d1 << "\n"; 
            os << "automatic differentiation:\n" << d2 << "\n"; 
            return ::testing::AssertionFailure() << os.str();
        }
    }
    for (size_t i=0;i<rot_ndx.size();++i) {
        for (int c=0;c<3;++c) {
        // check derivative towards variable i:
            Td d1 = numerical_derivative_rot<T>(e,rot_ndx[i],c, rot_value[i],  h );
            Td d2 = e->derivative(rot_ndx[i]+c);
            if (!Equal(d1,d2,tol)) {
                std::stringstream os;
                os << "check using numerical derivative failed for " << mstr << " : \n";  
                os << "derivative towards rot variable " << scalar_ndx[i] << " and component " << c << "\n";
                os << "numerical differentiation:\n" << d1 << "\n"; 
                os << "automatic differentiation:\n" << d2 << "\n"; 
                return ::testing::AssertionFailure() << os.str();
            }
        }
    }
    return ::testing::AssertionSuccess();
}

#define CHECK_ROT_WITH_NUM( a ) \
    EXPECT_PRED_FORMAT1(CheckNumericalRot, a );




} // namespace KDL 
#endif
