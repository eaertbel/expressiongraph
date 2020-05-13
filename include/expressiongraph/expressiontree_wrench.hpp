/*
 * expressiontree_wrench.hpp
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

#ifndef KDL_EXPRESSIONTREE_WRENCH_HPP
#define KDL_EXPRESSIONTREE_WRENCH_HPP

#include "expressiontree_expressions.hpp"
#include <kdl/frames.hpp>

namespace KDL {

/*
 * Wrench operations:
 *     -Wrench                		returns Wrench
 *     Wrench+Wrench           		returns Wrench
 *     Wrench-Wrench           		returns Wrench
 *     Rotation*Wrench        		returns Wrench
 *     Wrench*double          		returns Wrench
 *     double*Wrench          		returns Wrench
 *     RefPoint(Wrench,Vector) 		returns Wrench
 *     Inverse(Stiffness, Wrench)	returns Twist
 */
Expression<Wrench>::Ptr wrench( Expression<KDL::Vector>::Ptr a, Expression<KDL::Vector>::Ptr b);

Expression<Vector>::Ptr force( Expression<Wrench>::Ptr a);

Expression<KDL::Vector>::Ptr torque( Expression<KDL::Wrench>::Ptr a);

Expression<KDL::Wrench>::Ptr operator-( Expression<KDL::Wrench>::Ptr a);

Expression<KDL::Wrench>::Ptr operator+( Expression<KDL::Wrench>::Ptr a1, Expression<KDL::Wrench>::Ptr a2);

Expression<KDL::Wrench>::Ptr operator-( Expression<KDL::Wrench>::Ptr a1, Expression<KDL::Wrench>::Ptr a2);

Expression<KDL::Wrench>::Ptr operator* ( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Wrench>::Ptr a2 );

Expression<KDL::Wrench>::Ptr operator* ( Expression<KDL::Wrench>::Ptr a1, Expression<double>::Ptr a2 );

Expression<KDL::Wrench>::Ptr operator* ( Expression<double>::Ptr a1, Expression<KDL::Wrench>::Ptr a2 );

Expression<KDL::Wrench>::Ptr ref_point ( Expression<KDL::Wrench>::Ptr a1, Expression<KDL::Vector>::Ptr a2 );

/*  **********************************************
//Inverse Stiffness Wrench
class Inverse_StiffnessWrench:
	public BinaryExpression<KDL::Twist, KDL::Stiffness, KDL::Wrench>
{
public:
	typedef BinaryExpression<KDL::Twist,KDL::Stiffness,KDL::Wrench> BinExpr;
private:
	KDL::Stiffness arg1value;
	KDL::Wrench arg2value;
public:
	Inverse_StiffnessWrench(
			const typename BinExpr::Argument1Expr::Ptr& arg1,
			const typename BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("inverse",arg1,arg2)
				{}

	virtual KDL::Twist value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value.Inverse(arg2value);
	}

	virtual KDL::Twist derivative(int i){ //TODO: Check implementation! The old implementation asserts false in the VV case.
		return arg1value.Inverse(argument2->derivative(i));
	}
    virtual Expression<Twist>::Ptr derivativeExpression(int i) {
        assert( 0 && "derivativeExpression for InverseStiffnessWrench NOT IMPLEMENTED (YET)");
    }


    virtual typename BinExpr::Ptr clone() {
        Expression<KDL::Twist>::Ptr expr(
            new Inverse_StiffnessWrench( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Twist>::Ptr inv ( Expression<KDL::Stiffness>::Ptr a1, Expression<KDL::Wrench>::Ptr a2 ) {
	Expression<KDL::Twist>::Ptr expr(new Inverse_StiffnessWrench( a1, a2 ));
	return expr;
}
**************************************************** */

} // end of namespace KDL

#endif
