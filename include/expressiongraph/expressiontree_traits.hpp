#ifndef KDL_EXPRESSIONTREE_TRAITS_HPP
#define KDL_EXPRESSIONTREE_TRAITS_HPP

/*
* expressiongraph library
* 
* Copyright 2014 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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
#ifdef __GNUC__
    #include <cxxabi.h>
    #include <malloc.h>
#endif
#include <Eigen/Core>
#include <expressiongraph/expressiontree_exceptions.hpp>



namespace KDL {

enum class ExpressionType {
    expression_double    = 0,
    expression_vector    = 1,
    expression_rotation  = 2,
    expression_frame     = 3,
    expression_twist     = 4,
    expression_wrench    = 5,
    expression_quaternion= 6,
    expression_matrix= 7,
    unknown   = 999 
};

// tollerance to consider a value equal.
// (e.g. during simplifications, equal to one)
extern double EG_EPS_EQUAL;


/**
 * Trait classes
 * Defines relationships between types and some utility methods
 * that can be used to write code about types in a more general way.
 */
template<typename T>
struct AutoDiffTrait {
    typedef T ValueType;
    typedef T DerivType;
    static ValueType zeroValue() {
        EG_ASSERT_MSG(0,"should never be called, only specialized versions should be called."); 
        return 0;
    }
    static DerivType zeroDerivative() {
        EG_ASSERT_MSG(0,"should never be called, only specialized versions should be called."); 
        return 0;
    }
    const static int size=0;
    const static ExpressionType expr_type = ExpressionType::unknown;
};

template<>
struct AutoDiffTrait<double> {
    typedef double ValueType;
    typedef double DerivType;
    static ValueType zeroValue() {
        return 0;
    }
    static DerivType zeroDerivative() {
        return 0;
    }
    const static int size=1;
    const static ExpressionType expr_type = ExpressionType::expression_double;
};

template<>
struct AutoDiffTrait<KDL::Vector>
{
    typedef KDL::Vector ValueType;
    typedef KDL::Vector DerivType;
    static ValueType zeroValue() {
        return KDL::Vector::Zero();
    }
    static DerivType zeroDerivative() {
        return KDL::Vector::Zero();
    }
    const static int size=3;
    const static ExpressionType expr_type = ExpressionType::expression_vector;
};

template<>
struct AutoDiffTrait <KDL::Rotation>
{
    typedef KDL::Rotation ValueType;
    typedef KDL::Vector DerivType;
    static ValueType zeroValue() {
        return KDL::Rotation::Identity();
    }
    static DerivType zeroDerivative() {
        return KDL::Vector::Zero();
    }
    const static int size=3;
    const static ExpressionType expr_type = ExpressionType::expression_rotation;
};

//Frame
template<>
struct AutoDiffTrait <KDL::Frame>
{
    typedef KDL::Frame ValueType;
    typedef KDL::Twist DerivType;
    static ValueType zeroValue() {
        return KDL::Frame::Identity();
    }
    static DerivType zeroDerivative() {
        return KDL::Twist::Zero();
    }
    const static int size=6;
    const static ExpressionType expr_type = ExpressionType::expression_frame;
};

//Wrench
template<>
struct AutoDiffTrait <KDL::Wrench>
{
    typedef KDL::Wrench ValueType;
    typedef KDL::Wrench DerivType;
    static ValueType zeroValue() {
        return KDL::Wrench::Zero();
    }
    static DerivType zeroDerivative() {
        return KDL::Wrench::Zero();
    }
    const static int size=6;
    const static ExpressionType expr_type = ExpressionType::expression_wrench;

};

//Twist
template<>
struct AutoDiffTrait <KDL::Twist>
{
    typedef KDL::Twist ValueType;
    typedef KDL::Twist DerivType;
    static ValueType zeroValue() {
        return KDL::Twist::Zero();
    }
    static DerivType zeroDerivative() {
        return KDL::Twist::Zero();
    }
    const static int size=6;
    const static ExpressionType expr_type = ExpressionType::expression_twist;
};

// Quaternion
template<>
    struct AutoDiffTrait<Quaternion> {
        typedef Quaternion ValueType;
        typedef Quaternion DerivType;
        static ValueType zeroValue() {
            return Quaternion( 0.0, Vector::Zero());
        }
        static DerivType zeroDerivative() {
            return Quaternion( 0.0, Vector::Zero());
        }
        const static int size=4;
        const static ExpressionType expr_type = ExpressionType::expression_quaternion;
    };


//Matrix
template<int n, int m>
struct AutoDiffTrait<Eigen::Matrix<double, n, m> >
{
    typedef Eigen::Matrix<double,n,m> ValueType;
    typedef Eigen::Matrix<double,n,m> DerivType;
    static ValueType zeroValue() {
        return Eigen::Matrix<double,n,m>::Zero();
    }
    static DerivType zeroDerivative() {
        return Eigen::Matrix<double,n,m>::Zero();
    }
    const static int size=n*m;
    const static ExpressionType expr_type = ExpressionType::expression_matrix;

};



template<typename T1, typename T2>
struct type_comparison {
    static T2 return_if_equal(const T1& arg) {
        return T2();
    }
};

template <typename T>
struct type_comparison<T,T> {
    static T return_if_equal(const T& arg) {
        return arg;
    }
};

#ifdef __GNUC__
    inline std::string demangle( const std::string& name) {
        int status;
        char* buf=abi::__cxa_demangle(name.c_str(), 0,0, &status);
        std::string fullname; 
        switch(status) {
            case 0:
                fullname=buf;
                break;
            case -1:
                fullname="MEMORY ALLOC FAILED";
                break;
            case -2:
                fullname="INVALID NAME";
                break;
            case -3:
                fullname="INVALID NAME";
                break;
        } 
        free(buf);
        return fullname;
    }
#else
    inline std::string demangle( const std::string& name) {
            return name;
    }
#endif


} // namespace KDL;
#endif
