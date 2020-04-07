/*
 * expressiontree_matrix.cpp
 *
 *
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

#ifndef KDL_EXPRESSIONTREE_MATRIX_HPP
#define KDL_EXPRESSIONTREE_MATRIX_HPP

#include <expressiongraph/expressiontree_expressions.hpp>
#include <Eigen/Core>


namespace KDL {

template< int n, int m>
typename Expression<Eigen::Matrix<double, n,m> >::Ptr 
addition( typename Expression< Eigen::Matrix<double,n,m> >::Ptr a1, 
            typename Expression< Eigen::Matrix<double,n,m> >::Ptr a2 );

template< int n, int m, int k>
typename Expression<Eigen::Matrix<double, n,k> >::Ptr 
multiply ( typename Expression< Eigen::Matrix<double,n,m> >::Ptr a1, 
            typename Expression< Eigen::Matrix<double,m,k> >::Ptr a2 );
	

template<int n, int m>
Expression<double>::Ptr get_element( int i,int j, Expression<KDL::Vector>::Ptr a);




	
template< int n, int m, int k>
class Matrix_Multiplication:
    public BinaryExpression<typename Eigen::Matrix<double,n,k>,  typename Eigen::Matrix<double,n,m>, typename Eigen::Matrix<double,m,k> >
{
    typedef BinaryExpression<typename Eigen::Matrix<double,n,k>, typename  Eigen::Matrix<double,n,m>,typename  Eigen::Matrix<double,m,k> > BinExpr;
    typename Eigen::Matrix<double,n,m> arg1value;
    typename Eigen::Matrix<double,m,k> arg2value;
    typename Eigen::Matrix<double,n,k> result;
public:
    Matrix_Multiplication(
			const   typename BinExpr::Argument1Expr::Ptr& arg1,
			const   typename BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("multiplication",arg1,arg2)
				{}
	virtual typename Eigen::Matrix<double,n,k> value() {
		arg1value = this->argument1->value();
		arg2value = this->argument2->value();
		result    = arg1value * arg2value;
        return result;
	}

	virtual typename Eigen::Matrix<double,n,k> derivative(int i){
        result = this->argument1->derivative(i) * arg2value;
        result += arg1value*this->argument2->derivative(i);
        return result;
	}

    virtual typename Expression<Eigen::Matrix<double,n,k> >::Ptr derivativeExpression(int i) {
        return addition<n,k>( multiply<n,m,k>( this->argument1, this->argument2->derivativeExpression(i) ),
                       multiply<n,m,k>( this->argument1->derivativeExpression(i), this->argument2));
    }

    virtual typename BinExpr::Ptr clone() {
            return boost::make_shared< Matrix_Multiplication<n,m,k> >( 
                this->argument1->clone(), 
                this->argument2->clone() 
            );
    }
};


/**
 * multiply two matrices, you will need to call this with 
 * explicit template parameters.
 */
template< int n, int m, int k>
inline typename Expression<Eigen::Matrix<double, n,k> >::Ptr 
multiply ( typename Expression< Eigen::Matrix<double,n,m> >::Ptr a1, 
            typename Expression< Eigen::Matrix<double,m,k> >::Ptr a2 ) {
	return boost::make_shared< Matrix_Multiplication<n,m,k>  >(a1,a2);
}

template< int n, int m>
class Matrix_Addition:
    public BinaryExpression<typename Eigen::Matrix<double,n,m>,  typename Eigen::Matrix<double,n,m>, typename Eigen::Matrix<double,n,m> >
{
    typedef BinaryExpression<typename Eigen::Matrix<double,n,m>, typename  Eigen::Matrix<double,n,m>,typename  Eigen::Matrix<double,n,m> > BinExpr;
    typename Eigen::Matrix<double,n,m> arg1value;
    typename Eigen::Matrix<double,n,m> arg2value;
    typename Eigen::Matrix<double,n,m> result;
public:
    Matrix_Addition(
			const   typename BinExpr::Argument1Expr::Ptr& arg1,
			const   typename BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("multiplication",arg1,arg2)
				{}
	virtual typename Eigen::Matrix<double,n,m> value() {
		arg1value = this->argument1->value();
		arg2value = this->argument2->value();
		result    = arg1value + arg2value;
        return result;
	}

	virtual typename Eigen::Matrix<double,n,m> derivative(int i){
        result = this->argument1->derivative(i) + this->argument2->derivative(i);
        return result;
	}

    virtual typename Expression<Eigen::Matrix<double,n,m> >::Ptr derivativeExpression(int i) {
        return addition<n,m>( this->argument1->derivativeExpression(i), this->argument2->derivativeExpression(i) );
    }

    virtual typename BinExpr::Ptr clone() {
            return boost::make_shared< Matrix_Addition<n,m> >( 
                this->argument1->clone(), 
                this->argument2->clone() 
            );
    }
};


/**
 * multiply two matrices, you will need to call this with 
 * explicit template parameters.
 */
template< int n, int m>
inline typename Expression<Eigen::Matrix<double, n,m> >::Ptr 
addition( typename Expression< Eigen::Matrix<double,n,m> >::Ptr a1, 
            typename Expression< Eigen::Matrix<double,n,m> >::Ptr a2 ) {
	return boost::make_shared< Matrix_Addition<n,m>  >(a1,a2);
}



//CoordX Vector
template<int n, int m>
class MatrixElement:
    public UnaryExpression<double, Eigen::Matrix<double,n,m> >
{
public:
    typedef UnaryExpression<double, Eigen::Matrix<double,n,m> > UnExpr;
    int i;
    int j;
public:
    MatrixElement(){}
    MatrixElement(const  typename UnExpr::ArgumentExpr::Ptr& arg, int _i, int _j):
                UnExpr("get_element",arg),
                i(_i),
                j(_j)
                {}

    virtual double value() {
    	return this->argument->value()(i,j);
    }

    virtual double derivative(int i) {
    	return this->argument->derivative(i)(i,j);
    }

    virtual typename Expression<double>::Ptr derivativeExpression(int c) {
        return get_element<n,m>(i,j,this->argument->derivativeExpression(c));
    }

    virtual  typename UnExpr::Ptr clone() {
        return boost::make_shared< MatrixElement >(this->argument->clone(), i, j);
    }
};

template<int n, int m>
inline Expression<double>::Ptr get_element( int i,int j, typename Expression<Eigen::Matrix<double,n,m> >::Ptr a) {
    if ((0<=i)&&(i<n)&&(0<=j)&&(j<m)) {
        return boost::make_shared<MatrixElement<n,m> >(a,i,j);
    } else {
        return Expression<double>::Ptr();
    }
}

} // end of namespace KDL
#endif
