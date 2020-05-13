/*
 * expressiontree_rotation.cpp
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

#ifndef KDL_EXPRESSIONTREE_ROTATION_HPP
#define KDL_EXPRESSIONTREE_ROTATION_HPP

#include <expressiongraph/expressiontree_expressions.hpp>
#include <expressiongraph/expressiontree_vector.hpp>
#include <expressiongraph/expressiontree_double.hpp>
#include <kdl/frames.hpp>


namespace KDL {
/*
 * Rotation
 *      - RotX(Expression<double>::Ptr)          					returns Expression<Rotation>::Ptr
 *      - RotY(Expression<double>::Ptr)          					returns Expression<Rotation>::Ptr
 *      - RotZ(Expression<double>::Ptr)          					returns Expression<Rotation>::Ptr
 *      - Inverse(Expression<Rotation>::Ptr)     					returns Expression<Rotation>::Ptr
 *      - Expression<Rotation>::Ptr*Expression<Rotation::Ptr> 		returns Expression<Rotation>::Ptr
 *      - Expression<Rotation>::Ptr*Expression<Vector::Ptr>			returns Expression<Vector>::Ptr
 *      - UnitX( Expression<Rotation>::Ptr)      					returns Expression<Vector>::Ptr
 *      - UnitY( Expression<Rotation>::Ptr)      					returns Expression<Vector>::Ptr
 *      - UnitZ( Expression<Rotation>::Ptr)      					returns Expression<Vector>::Ptr
 *      - construct_rotation_from_vectors                           returns Expression<Rotation>::Ptr
 */
// Rot Double

Expression<KDL::Rotation>::Ptr rot(const KDL::Vector& axis, Expression<double>::Ptr a);

Expression<KDL::Rotation>::Ptr rotVec(Expression<Vector>::Ptr a, Expression<double>::Ptr b);

Expression<KDL::Rotation>::Ptr rot_x( Expression<double>::Ptr a);

Expression<KDL::Rotation>::Ptr rot_y( Expression<double>::Ptr a);

Expression<KDL::Rotation>::Ptr rot_z( Expression<double>::Ptr a);

Expression<KDL::Rotation>::Ptr inv (Expression<KDL::Rotation>::Ptr a);

Expression<KDL::Rotation>::Ptr operator* ( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Rotation>::Ptr a2 );

Expression<KDL::Vector>::Ptr operator* ( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Vector>::Ptr a2 );

Expression<KDL::Vector>::Ptr unit_x ( Expression<KDL::Rotation>::Ptr a);

Expression<KDL::Vector>::Ptr unit_y ( Expression<KDL::Rotation>::Ptr a);

Expression<KDL::Vector>::Ptr unit_z ( Expression<KDL::Rotation>::Ptr a);

Expression<KDL::Rotation>::Ptr construct_rotation_from_vectors( Expression<KDL::Vector>::Ptr a, Expression<KDL::Vector>::Ptr b, Expression<KDL::Vector>::Ptr c);

Expression<KDL::Vector>::Ptr getRotVec( Expression<KDL::Rotation>::Ptr a);

Expression<KDL::Vector>::Ptr getRPY( Expression<KDL::Rotation>::Ptr a);


// just because the analog exists for rotation matrices, we define it
// also for rotation matrices:
Expression<Vector>::Ptr apply ( 
        Expression<Rotation>::Ptr F, Expression<Vector>::Ptr v);

} // end of namespace KDL
#endif
