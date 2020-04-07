/**
 * quat_io.hpp
 * expressiongraph library
 *
 * I/O for quaternions
 *
 * Copyright 2019 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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

#ifndef KDL_QUAT_IO_HPP
#define KDL_QUAT_IO_HPP

#include "quat.hpp"
#include <iostream>
namespace KDL {
inline std::ostream& operator << (std::ostream& os,const KDL::Quaternion& q) {
    os << "[ "<< q.w << " +  " << q.vec[0] << " i +  " << q.vec[1] << " j + " << q.vec[2] << " k  ]";
    return os;
}
} // namespace KDL
#endif 
