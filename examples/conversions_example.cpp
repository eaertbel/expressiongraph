/*
 * expressiontree_example.cpp
 *
 *  Created on: Sept., 2012
 *      Author: Erwin Aertbelien 
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

#include <kdl/frames.hpp>
#include <kdl/frames_io.hpp>
#include <expressiongraph/conversions.hpp>

/*
 * This small example show the conversion between Eigen types and KDL types
 *
 * The toKDLxxx() function specify the KDL type in the name because it is impossible
 * to determine the KDL type at compile time when using a dynamic eigen matrix.
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;
    cout << "Demonstration of simple conversion routines between KDL and Eigen types \n";
    Vector a(1.0,2.0,3.0);
    cout << "KDL Vector " << a << "\n";
    Eigen::Vector3d a_eigen;
    a_eigen = toEigen(a);
    cout << "Eigen representation " << a_eigen.transpose() << "\n";
    cout << "converted back to KDL " << toKDLVector(a_eigen) << "\n";
    cout << endl;

    Twist t(Vector(1.0,2.0,3.0),Vector(4.0,5.0,6.0));
    cout << "KDL Vector " << t << "\n";
    Eigen::Matrix<double,6,1> t_eigen;
    t_eigen = toEigen(t);
    cout << "Eigen representation " << t_eigen.transpose() << "\n";
    cout << "converted back to KDL " << toKDLTwist(t_eigen) << "\n";
    cout << endl;


    Wrench w(Vector(1.0,2.0,3.0),Vector(4.0,5.0,6.0));
    cout << "KDL Wrench " << t << "\n";
    Eigen::Matrix<double,6,1> w_eigen;
    w_eigen = toEigen(w);
    cout << "Eigen representation " << w_eigen.transpose() << "\n";
    cout << "converted back to KDL " << toKDLTwist(w_eigen) << "\n";

    Rotation R( Rotation::EulerZYX(30*deg2rad,45*deg2rad,0*deg2rad) );
    cout << "KDL Rotation " << R << "\n";
    Eigen::Matrix<double,3,3> R_eigen;
    R_eigen = toEigen( R );
    cout << "Eigen 3x3 matrix " << R_eigen << "\n";
    cout << "converted back to KDL " << toKDLRotation( R_eigen ) << "\n";
    Eigen::Quaternion<double> q;
    q = toEigenQuaternion( R );
    cout << "Eigen Quaternion " << q.coeffs() << "\n";
    cout << "converted back to KDL " << toKDLRotation( q ) << "\n"; 
    cout << endl;
	return 0;
}
