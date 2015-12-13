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

namespace KDL {

/**
 * Trait classes
 * Defines relationships between types and some utility methods
 * that can be used to write code about types in a more general way.
 */
template<typename T>
struct AutoDiffTrait {
    typedef T ValueType;
    typedef T DerivType;
    static DerivType zeroDerivative() {
        assert(0 /* should never be called, only specialized versions should be called. */); 
        return 0;
    }
    const static int size=0;
};

template<>
struct AutoDiffTrait<double> {
    typedef double ValueType;
    typedef double DerivType;
    static DerivType zeroDerivative() {
        return 0;
    }
    const static int size=1;
};

template<>
struct AutoDiffTrait<KDL::Vector>
{
    typedef KDL::Vector ValueType;
    typedef KDL::Vector DerivType;
    static DerivType zeroDerivative() {
        return KDL::Vector::Zero();
    }
    const static int size=3;

};

template<>
struct AutoDiffTrait <KDL::Rotation>
{
    typedef KDL::Rotation ValueType;
    typedef KDL::Vector DerivType;
    static DerivType zeroDerivative() {
        return KDL::Vector::Zero();
    }
    const static int size=3;
};

//Frame
template<>
struct AutoDiffTrait <KDL::Frame>
{
    typedef KDL::Frame ValueType;
    typedef KDL::Twist DerivType;
    static DerivType zeroDerivative() {
        return KDL::Twist::Zero();
    }
    const static int size=6;

};

//Wrench
template<>
struct AutoDiffTrait <KDL::Wrench>
{
    typedef KDL::Wrench ValueType;
    typedef KDL::Wrench DerivType;
    static DerivType zeroDerivative() {
        return KDL::Wrench::Zero();
    }
    const static int size=6;

};

//Twist
template<>
struct AutoDiffTrait <KDL::Twist>
{
    typedef KDL::Twist ValueType;
    typedef KDL::Twist DerivType;
    static DerivType zeroDerivative() {
        return KDL::Twist::Zero();
    }
    const static int size=6;

};

//Matrix
template<int n, int m>
struct AutoDiffTrait<Eigen::Matrix<double, n, m> >
{
    typedef Eigen::Matrix<double,n,m> ValueType;
    typedef Eigen::Matrix<double,n,m> DerivType;
    static DerivType zeroDerivative() {
        return Eigen::Matrix<double,n,m>::Zero();
    }
    const static int size=n*m;

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
