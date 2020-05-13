/*
 * expressiontree_twist.hpp
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

#ifndef KDL_EXPRESSIONTREE_TWIST_HPP
#define KDL_EXPRESSIONTREE_TWIST_HPP

#include "expressiontree_expressions.hpp"
#include <kdl/frames.hpp>

namespace KDL {
/*
 * Twist operations:
 *     -Twist                	returns Twist
 *     Twist+Twist           	returns Twist
 *     Twist-Twist           	returns Twist
 *     Rotation*Twist        	returns Twist
 *     Twist*double          	returns Twist
 *     double*Twist          	returns Twist
 *     RefPoint(Twist,Vector) 	returns Twist
 *     Stiffness*Twist			returns Wrench
 *     vel() returns a Vector corresponding to translational velocity
 *     rot() returns a Vector corresponding to rotational velocity
 */


Expression<KDL::Twist>::Ptr twist( Expression<KDL::Vector>::Ptr a, Expression<KDL::Vector>::Ptr b);

Expression<KDL::Twist>::Ptr operator-( Expression<KDL::Twist>::Ptr a);

Expression<Vector>::Ptr transvel( Expression<Twist>::Ptr a);

Expression<KDL::Vector>::Ptr rotvel( Expression<KDL::Twist>::Ptr a);

Expression<KDL::Twist>::Ptr operator+( Expression<KDL::Twist>::Ptr a1, Expression<KDL::Twist>::Ptr a2);

Expression<KDL::Twist>::Ptr operator-( Expression<KDL::Twist>::Ptr a1, Expression<KDL::Twist>::Ptr a2);

Expression<KDL::Twist>::Ptr operator*( Expression<KDL::Rotation>::Ptr a1, Expression<KDL::Twist>::Ptr a2);

Expression<KDL::Twist>::Ptr operator*( Expression<KDL::Twist>::Ptr a1, Expression<double>::Ptr a2);

Expression<KDL::Twist>::Ptr operator*( Expression<double>::Ptr a1, Expression<KDL::Twist>::Ptr a2);

Expression<KDL::Twist>::Ptr ref_point ( Expression<KDL::Twist>::Ptr a1, Expression<KDL::Vector>::Ptr a2);


/**
 * /todo Finish composition with stiffness
 *
//Composition Stiffness, Twist
class Composition_StiffnessTwist:
	public BinaryExpression<KDL::Wrench, KDL::Stiffness, KDL::Twist>
{
public:
	typedef BinaryExpression<KDL::Wrench,KDL::Stiffness,KDL::Twist> BinExpr;
private:
	KDL::Stiffness arg1value;
	KDL::Twist arg2value;
public:
	Composition_StiffnessTwist(
			const typename BinExpr::Argument1Expr::Ptr& arg1,
			const typename BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("composition",arg1,arg2)
				{}

	virtual KDL::Wrench value() {
		arg1value = argument1->value();
		arg2value = argument2->value();
		return arg1value * arg2value;
	}

	virtual KDL::Wrench derivative(int i){
		return arg1value * argument2->derivative(i) + argument1->derivative(i)*arg2value;
	}

    virtual typename BinExpr::Ptr clone() {
        Expression<KDL::Wrench>::Ptr expr(
            new Composition_StiffnessTwist( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Wrench>::Ptr operator* ( Expression<KDL::Stiffness>::Ptr a1, Expression<KDL::Twist>::Ptr a2 ) {
	Expression<KDL::Wrench>::Ptr expr(new Composition_StiffnessTwist( a1, a2 ));
	return expr;
}
*/

} // end of namespace KDL

#endif
