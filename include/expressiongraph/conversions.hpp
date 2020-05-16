#ifndef KDL_CONVERSIONS_HPP
#define KDL_CONVERSIONS_HPP
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



#include <Eigen/Dense>
#include <kdl/frames.hpp>
#include <expressiongraph/expressiontree_exceptions.hpp>
/** @file conversions.hpp 
 *  @brief Conversions to/from KDL types
 */ 


namespace KDL {
    /** \addtogroup conversions Conversion to/from KDL Types. 
     *  @{
     */

    /** 
     * converts an Eigen vector to KDL::Vector
     */
    inline Vector toKDLVector( const Eigen::Vector3d& arg) {
        Vector tmp( arg(0), arg(1), arg(2) );
        return tmp;
    }


    /** 
     * converts an Eigen vector to KDL::Vector
     * The Eigen vector should have a length of 3.
     */
   inline  Vector toKDLVector( const Eigen::VectorXd& arg) {
        EG_ASSERT_MSG( arg.rows()==3,"Eigen::VectorXD should have 3 elements for conversion to KDL::Vector" );
        Vector tmp( arg(0), arg(1), arg(2) );
        return tmp;
    }

    /**
     * converts a KDL::Vector to Eigen vector (fixed size, allocated on stack, real-time)
     */
    inline Eigen::Vector3d toEigen( const KDL::Vector& arg ) {
        Eigen::Vector3d tmp;
        tmp(0)=arg.x();
        tmp(1)=arg.y();
        tmp(2)=arg.z();
        return tmp;
    }


    /** 
     * converts an Eigen vector to KDL::Twist 
     */
    inline KDL::Twist toKDLTwist( const Eigen::Matrix<double,6,1>& arg) {
        Twist tmp( Vector(arg(0), arg(1), arg(2)), Vector( arg(3),arg(4),arg(5)) );
        return tmp;
    }


    /** 
     * converts an Eigen vector to KDL::Twist 
     * The Eigen vector should have a length of 6.
     */
    inline KDL::Twist toKDLTwist( const Eigen::VectorXd& arg) {
        EG_ASSERT_MSG( arg.rows()==6,"Eigen::VectorXD should have 6 elements for conversion to KDL::Twist" );
        Twist tmp( Vector(arg(0), arg(1), arg(2)), Vector( arg(3),arg(4),arg(5)) );
        return tmp;
    }

    /**
     * converts a KDL::Twist to Eigen vector (fixed size, allocated on stack, real-time)
     */
    inline Eigen::Matrix<double,6,1> toEigen( const KDL::Twist& arg ) {
        Eigen::Matrix<double,6,1> tmp;
        tmp(0)=arg.vel.x();
        tmp(1)=arg.vel.y();
        tmp(2)=arg.vel.z();
        tmp(3)=arg.rot.x();
        tmp(4)=arg.rot.y();
        tmp(5)=arg.rot.z();
        return tmp;
    }
    
    /** 
     * converts an Eigen vector to KDL::Wrench 
     */
    inline KDL::Wrench toKDLWrench( const Eigen::Matrix<double,6,1>& arg) {
        Wrench tmp( Vector(arg(0), arg(1), arg(2)), Vector( arg(3),arg(4),arg(5)) );
        return tmp;
    }

    /** 
     * converts an Eigen vector to KDL::Wrench 
     * The Eigen vector should have a length of 6.
     */
    inline KDL::Wrench toKDLWrench( const Eigen::VectorXd& arg) {
        EG_ASSERT_MSG( arg.rows()==6,"Eigen::VectorXD should have 6 elements for conversion to KDL::Wrench" );
        Wrench tmp( Vector(arg(0), arg(1), arg(2)), Vector( arg(3),arg(4),arg(5)) );
        return tmp;
    }

    /**
     * converts a KDL::Wrench to Eigen vector (fixed size, allocated on stack, real-time)
     */
    inline Eigen::Matrix<double,6,1> toEigen( const KDL::Wrench& arg ) {
        Eigen::Matrix<double,6,1> tmp;
        tmp(0)=arg.force.x();
        tmp(1)=arg.force.y();
        tmp(2)=arg.force.z();
        tmp(3)=arg.torque.x();
        tmp(4)=arg.torque.y();
        tmp(5)=arg.torque.z();
        return tmp;
    }

    /**
     * converts Eigen matrix 3x3 to a KDL::Rotation
     */
    inline KDL::Rotation toKDLRotation(const Eigen::Matrix<double,3,3>& arg) {
        KDL::Rotation R( 
                        arg(0,0), arg(0,1), arg(0,2),
                        arg(1,0), arg(1,1), arg(1,2),
                        arg(2,0), arg(2,1), arg(2,2) 
        );
        return R;                  
    }
    
    /**
     * converts Eigen matrix to a KDL::Rotation
     * size should be 3 x 3
     */
    inline KDL::Rotation toKDLRotation(const Eigen::MatrixXd& arg) {
        EG_ASSERT_MSG( (arg.rows()==3) && (arg.cols()==3) , "Eigen MatrixXd should be 3x3 for conversion to KDL::Rotation"); 
        KDL::Rotation R( 
                        arg(0,0), arg(0,1), arg(0,2),
                        arg(1,0), arg(1,1), arg(1,2),
                        arg(2,0), arg(2,1), arg(2,2) 
        );
        return R; 
    }
   
    /**
     * converts a KDL::Rotation to an Eigen matrix
     */ 
    inline Eigen::Matrix<double,3,3> toEigen(const KDL::Rotation& arg) {
        Eigen::Matrix<double,3,3> m;
        m.col(0) = toEigen( arg.UnitX() );
        m.col(1) = toEigen( arg.UnitY() );
        m.col(2) = toEigen( arg.UnitZ() );
        return m;
    }

    /**
     * converts an Eigen Quaternion to a KDL::Rotation
     */
    inline KDL::Rotation toKDLRotation(const Eigen::Quaternion<double>& arg) {
        return toKDLRotation( arg.toRotationMatrix() );
    }
   
    /**
     * converts a KDL::Rotation to an Eigen Quaternion
     */ 
    inline Eigen::Quaternion<double> toEigenQuaternion(const KDL::Rotation& arg) {
        return Eigen::Quaternion<double>( toEigen(arg) );
    }

    /**
     * @}
     * end of group conversions
     */
};
#endif
