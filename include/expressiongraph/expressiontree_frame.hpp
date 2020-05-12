/*
 * expressiontree_frame.hpp
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

#ifndef KDL_EXPRESSIONTREE_FRAME_HPP
#define KDL_EXPRESSIONTREE_FRAME_HPP

#include <expressiongraph/expressiontree_expressions.hpp>
#include <kdl/frames.hpp>
/*
 * \todo : complete documentation
 * Frame related expression tree classes.
 * Operations have the same semantics as the operations in frames.hpp,
 * (naming is similar, except for member functions that are replaced by functions) 
 *
 * Frame unary operations:
 *     - inv
 * Frame binary operations:
 *     - "*"
 *     -
 * Corresponding Frame unary Expression classes:
 *     -  Inverse_Frame
 * Corresponding Frame binary Expression classes:
 *     - Composition_FrameFrame
 *     - Composition_FrameVector
 *
 * 
 *      - Inverse(Jacobian<Frame>)        			returns Jacobian<Frame>
 *      - Jacobian<Frame> * Jacobian<Frame>     	returns Jacobian<Frame>
 *      - Jacobian<Frame> * Jacobian<Vector>    	returns Jacobian<Vector>
 *      - Origin( Jacobian<Frame> )       			returns Jacobian<Vector>
 *      - RotMat( Jacobian<Frame> )      			returns Jacobian<Rotation>
 *      change_coordinate_frame ...
 * Constructor:
 *    Frame_RotationVector -> frame(Expression<Rotation>::Ptr, Expression<Vector>::Ptr )
 *                         -> frame( Expression<Vector>::Ptr )
 *                         -> frame( Expression<Rotation>::Ptr )
 *                     
 */

namespace KDL {

Expression<Frame>::Ptr frame( Expression<Rotation>::Ptr a1, 
                                Expression<Vector>::Ptr a2);

Expression<Frame>::Ptr frame( Expression<Rotation>::Ptr a1 );

Expression<Frame>::Ptr frame(  Expression<Vector>::Ptr a2);

Expression<KDL::Frame>::Ptr inv ( Expression<KDL::Frame>::Ptr a);

Expression<KDL::Frame>::Ptr operator* ( Expression<KDL::Frame>::Ptr a1, Expression<KDL::Frame>::Ptr a2 );

Expression<KDL::Vector>::Ptr operator* ( Expression<KDL::Frame>::Ptr a1, Expression<KDL::Vector>::Ptr a2 );

Expression<KDL::Vector>::Ptr origin ( Expression<KDL::Frame>::Ptr a);

Expression<KDL::Rotation>::Ptr rotation( Expression<KDL::Frame>::Ptr a);

// just because the analog exists for rotation matrices, we define it
// also for rotation matrices:
Expression<Vector>::Ptr apply ( 
        Expression<Frame>::Ptr F, Expression<Vector>::Ptr v);

}; // namespace KDL
#endif
