#ifndef EXPRESSIONGRAPH_INTEGRAL_51566a22_0d16_4ad0_b673_fa9a6d3cecec
#define EXPRESSIONGRAPH_INTEGRAL_51566a22_0d16_4ad0_b673_fa9a6d3cecec

/*
* expressiongraph library
* 
* Copyright 2020 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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

#include "expressiontree_expressions.hpp"
#include <string> 
#include <map>
#include <vector>

namespace KDL { 

    /**
     * constructs an expressiongraph node that integrates.
     * _integrand is a function definition that has one function parameter
     * _lower is an expression for the lower limit of the integral
     * _upper is an expression for the upper limit of the integral
     * _epsilon is the requested numerical accuracy
     * _maxRecDepth is the maximal recursion depth of the routine, 
     *
     */

    Expression<double>::Ptr make_Integral( 
                FunctionDefinition::Ptr _integrand, 
                Expression<double>::Ptr _lower, 
                Expression<double>::Ptr _upper, 
                double _epsilon, 
                double _maxstepsize, 
                int _maxRecDepth );

} // namespace KDL
#endif 
