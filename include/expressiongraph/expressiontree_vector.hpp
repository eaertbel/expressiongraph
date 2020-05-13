/*
 * expressiontree_vector.cpp
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

#ifndef KDL_EXPRESSIONTREE_VECTOR_HPP
#define KDL_EXPRESSIONTREE_VECTOR_HPP

#include <expressiongraph/expressiontree_expressions.hpp>
#include <kdl/frames.hpp>

/*
 * Vector operations:
 *      - dot(Jacobian<Vector>, Jacobian<Vector>) returns Jacobian<Vector>(dot-product)
 *      - Jacobian<Vector> * Jacobian<Vector>   returns Jacobian<Vector>  (cross-product)
 *      - Jacobian<Vector> + Jacobian<Vector>   returns Jacobian<Vector>
 *      - Jacobian<Vector> - Jacobian<Vector>   returns Jacobian<Vector>
 *      -    - Jacobian<Vector>           returns Jacobian<Vector>
 *      - norm( Jacobian<Vector> )        returns Jacobian<double> representing the norm of the vector
 *      - Jacobian<Vector> * Jacobian<double>   returns Jacobian<Vector>
 *      - Jacobian<double> * Jacobian<Vector>   returns Jacobian<Vector>
 *      - CoordX( Jacobian<Vector> )      returns Jacobian<double>
 *      - CoordY( Jacobian<Vector> )      returns Jacobian<double>
 *      - CoordZ( Jacobian<Vector> )      returns Jacobian<double>
 *      - Diff( Jacobian<Vector>, Jacobian<Vector>) returns Jacobian<Vector>
 *
 * Constructor for Vector
 */

namespace KDL {
Expression<Vector>::Ptr vector( Expression<double>::Ptr a1, 
                                Expression<double>::Ptr a2,
                                Expression<double>::Ptr a3);




Expression<double>::Ptr dot( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2);

Expression<KDL::Vector>::Ptr operator*( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2);


Expression<KDL::Vector>::Ptr cross( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2);



Expression<KDL::Vector>::Ptr operator+( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2);


Expression<KDL::Vector>::Ptr operator-( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2);

Expression<KDL::Vector>::Ptr operator-( Expression<KDL::Vector>::Ptr a);

Expression<double>::Ptr squared_norm ( Expression<KDL::Vector>::Ptr a);

Expression<double>::Ptr norm ( Expression<KDL::Vector>::Ptr a);

Expression<KDL::Vector>::Ptr operator*( Expression<KDL::Vector>::Ptr a1, Expression<double>::Ptr a2);

Expression<KDL::Vector>::Ptr operator*( Expression<double>::Ptr a1, Expression<KDL::Vector>::Ptr a2);

Expression<double>::Ptr coord_x ( Expression<KDL::Vector>::Ptr a);

Expression<double>::Ptr coord_y ( Expression<KDL::Vector>::Ptr a);

Expression<double>::Ptr coord_z ( Expression<KDL::Vector>::Ptr a);

Expression<KDL::Vector>::Ptr diff ( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2 );

}; //namespace KDL
#endif
