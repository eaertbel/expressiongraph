/*
 * expressiontree_vector.cpp
 *
 *  Created on: Aug 7, 2012
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

#include <expressiongraph/expressiontree_double.hpp>
#include <expressiongraph/expressiontree_vector.hpp>
#include "util.hpp"
 
namespace KDL {

class Vector_DoubleDoubleDouble:
    public TernaryExpression<Vector,double,double,double>
{
public:

    Vector_DoubleDoubleDouble(){}
    Vector_DoubleDoubleDouble(  Expression<double>::Ptr a1, 
                                Expression<double>::Ptr a2,
                                Expression<double>::Ptr a3):
        TernaryExpression<Vector,double,double,double>("vector",a1,a2,a3) {}

    virtual Vector value() {
        return Vector(argument1->value(),argument2->value(),argument3->value());
    } 
    virtual Vector derivative(int i) {
        return Vector(
            argument1->derivative(i),
            argument2->derivative(i),
            argument3->derivative(i)
        );
    } 

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  Expression<Vector>::Ptr clone() {
         Expression<Vector>::Ptr expr(
            new Vector_DoubleDoubleDouble(argument1->clone(), argument2->clone(), argument3->clone())
        );
        return expr;
    } 
};
Expression<Vector>::Ptr vector( Expression<double>::Ptr a1, 
                                Expression<double>::Ptr a2,
                                Expression<double>::Ptr a3) {
    Expression<Vector>::Ptr expr(
        new Vector_DoubleDoubleDouble( a1, a2, a3 )
    );
    return expr;
}


//Dot vector vector
class Dot_VectorVector:
    public BinaryExpression<double, KDL::Vector, KDL::Vector>
{
public:
    typedef BinaryExpression<double, KDL::Vector, KDL::Vector> BinExpr;
    KDL::Vector arg1value;
    KDL::Vector arg2value;

    Dot_VectorVector(){}
    Dot_VectorVector(
                    const  BinExpr::Argument1Expr::Ptr& arg1,
                    const  BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("dot",arg1,arg2)
                    {}
    virtual double value() {
        arg1value = argument1->value();
        arg2value = argument2->value();
        return dot(arg1value,arg2value);
    }

    virtual double derivative(int i) {
        return dot(arg1value,argument2->derivative(i)) + dot(argument1->derivative(i),arg2value);
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Dot_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<double>::Ptr dot( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<double>::Ptr expr(
        new Dot_VectorVector( a1, a2 )
    );
    return expr;
}

//Composition Vector Vector
class CrossProduct_VectorVector:
    public BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector>
{
public:
    typedef BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector> BinExpr;
    KDL::Vector arg1value;
    KDL::Vector arg2value;
public:
    CrossProduct_VectorVector(){}
    CrossProduct_VectorVector(
                    const  BinExpr::Argument1Expr::Ptr& arg1,
                    const  BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("cross",arg1,arg2)
                    {}
    virtual KDL::Vector value() {
        arg1value = argument1->value();
        arg2value = argument2->value();
        return arg1value * arg2value;
    }

    virtual KDL::Vector derivative(int i) {
        return arg1value * argument2->derivative(i) + argument1->derivative(i) * arg2value;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new CrossProduct_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};
Expression<KDL::Vector>::Ptr operator*( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<KDL::Vector>::Ptr expr(
        new CrossProduct_VectorVector( a1, a2 )
    );
    return expr;
}

Expression<KDL::Vector>::Ptr cross( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<KDL::Vector>::Ptr expr(
        new CrossProduct_VectorVector( a1, a2 )
    );
    return expr;
}


//Addition Vector Vector
class Addition_VectorVector:
    public BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector>
{
public:
    typedef BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector> BinExpr;
public:
    Addition_VectorVector(){}
    Addition_VectorVector(
                    const  BinExpr::Argument1Expr::Ptr& arg1,
                    const  BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("add",arg1,arg2)
                    {}
    virtual KDL::Vector value() {
        return argument1->value() + argument2->value();
    }

    virtual KDL::Vector derivative(int i) {
        return argument1->derivative(i) + argument2->derivative(i);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Addition_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Vector>::Ptr operator+( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<KDL::Vector>::Ptr expr(
        new Addition_VectorVector( a1, a2 )
    );
    return expr;
}


//Subtraction Vector Vector
class Subtraction_VectorVector:
    public BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector>
{
public:
    typedef BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector> BinExpr;
public:
    Subtraction_VectorVector(){}
    Subtraction_VectorVector(
                    const  BinExpr::Argument1Expr::Ptr& arg1,
                    const  BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("sub",arg1,arg2)
                    {}

    virtual KDL::Vector value() {
        return argument1->value() - argument2->value();
    }

    virtual KDL::Vector derivative(int i) {
        return argument1->derivative(i) - argument2->derivative(i);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Subtraction_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Vector>::Ptr operator-( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<KDL::Vector>::Ptr expr(
        new Subtraction_VectorVector( a1, a2 )
    );
    return expr;
}

//Negate Vector
class Negate_Vector:
    public UnaryExpression<KDL::Vector, KDL::Vector>
{
public:
    typedef UnaryExpression<KDL::Vector, KDL::Vector> UnExpr;
public:
    Negate_Vector(){}
    Negate_Vector(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("negate",arg)
                {}

    virtual KDL::Vector value() {
    	return -argument->value();
    }

    virtual KDL::Vector derivative(int i) {
    	return -argument->derivative(i);
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Negate_Vector( argument->clone())
        );
        return expr;
    }
};
Expression<KDL::Vector>::Ptr operator-( Expression<KDL::Vector>::Ptr a) {
    Expression<KDL::Vector>::Ptr expr(
        new Negate_Vector( a )
    );
    return expr;
}

//Norm Vector
class SquaredNorm_Vector:
    public UnaryExpression<double, KDL::Vector>
{
public:
    typedef UnaryExpression<double, KDL::Vector> UnExpr;
    KDL::Vector val;

public:
    SquaredNorm_Vector(){}
    SquaredNorm_Vector(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("squared_norm",arg)
                {}

    virtual double value() {
    	val = argument->value();
    	return val.x()*val.x() + val.y()*val.y() + val.z()*val.z();
    }

    virtual double derivative(int i) {
        Vector vald = argument->derivative(i);
        return 2.0*val.x()*vald.x() + 2.0*val.y()*vald.y() + 2.0*val.z()*vald.z();
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new SquaredNorm_Vector( argument->clone())
        );
        return expr;
    }
};

Expression<double>::Ptr squared_norm ( Expression<KDL::Vector>::Ptr a) {
    Expression<double>::Ptr expr(
        new SquaredNorm_Vector( a )
    );
    return expr;
}

//Norm Vector
class Norm_Vector:
    public UnaryExpression<double, KDL::Vector>
{
public:
    typedef UnaryExpression<double, KDL::Vector> UnExpr;
    KDL::Vector val;
    double nval;
public:
    Norm_Vector(){}
    Norm_Vector(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("norm",arg)
                {}

    virtual double value() {
    	val = argument->value();
        nval = sqrt( val.x()*val.x() + val.y()*val.y() + val.z()*val.z()+1E-12);
    	return nval;
    }

    virtual double derivative(int i) {
    	// of course, when norm==0, there exists no derivative,...
        // derivative= (x*xdot+y*ydot+z*zdot)/n
        // where n=sqrt(x*x+y*y+z*z)
        // limit if norm->0: 
        //      x*xdot/n + y*ydot/n + z*zdot/n
        //      x*xdot/n with x->0 (because norm->0)  
        //      becomes sign(x)*xdot
        //   so in the limit:
        //      sign(x)*xdot + sign(y)*ydot + sign(z)*zdot 
        return dot(val, argument->derivative(i))/nval;
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Norm_Vector( argument->clone())
        );
        return expr;
    }
};
Expression<double>::Ptr norm ( Expression<KDL::Vector>::Ptr a) {
    Expression<double>::Ptr expr(
        new Norm_Vector( a )
    );
    return expr;
}

//Multiplication Vector Double
class Multiplication_VectorDouble:
    public BinaryExpression<KDL::Vector, KDL::Vector, double>
{
public:
    typedef BinaryExpression<KDL::Vector, KDL::Vector, double> BinExpr;
    KDL::Vector arg1value;
    double arg2value;

public:
    Multiplication_VectorDouble(){}
    Multiplication_VectorDouble(
                    const  BinExpr::Argument1Expr::Ptr& arg1,
                    const  BinExpr::Argument2Expr::Ptr& arg2):
                        BinExpr("scalarprod",arg1,arg2)
                    {}

    virtual KDL::Vector value() {
        arg1value = argument1->value();
        arg2value = argument2->value();
        return arg1value * arg2value;
    }

    virtual KDL::Vector derivative(int i) {
        return arg1value * argument2->derivative(i) + argument1->derivative(i) * arg2value;
    }

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Multiplication_VectorDouble( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Vector>::Ptr operator*( Expression<KDL::Vector>::Ptr a1, Expression<double>::Ptr a2) {
    Expression<KDL::Vector>::Ptr expr(
        new Multiplication_VectorDouble( a1, a2 )
    );
    return expr;
}

Expression<KDL::Vector>::Ptr operator*( Expression<double>::Ptr a1, Expression<KDL::Vector>::Ptr a2) {
    Expression<KDL::Vector>::Ptr expr(
        new Multiplication_VectorDouble( a2, a1 )
    );
    return expr;
}

//CoordX Vector
class CoordX_Vector:
    public UnaryExpression<double, KDL::Vector>
{
public:
    typedef UnaryExpression<double, KDL::Vector> UnExpr;
public:
    CoordX_Vector(){}
    CoordX_Vector(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("coord_x",arg)
                {}

    virtual double value() {
    	return argument->value()[0];
    }

    virtual double derivative(int i) {
    	return argument->derivative(i)[0];
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new CoordX_Vector( argument->clone())
        );
        return expr;
    }
};

Expression<double>::Ptr coord_x ( Expression<KDL::Vector>::Ptr a) {
    Expression<double>::Ptr expr(
        new CoordX_Vector(a)
    );
    return expr;
}



//CoordY Vector
class CoordY_Vector:
    public UnaryExpression<double, KDL::Vector>
{
public:
    typedef UnaryExpression<double, KDL::Vector> UnExpr;
public:
    CoordY_Vector(){}
    CoordY_Vector(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("coord_y",arg)
                {}

    virtual double value() {
    	return argument->value()[1];
    }

    virtual double derivative(int i) {
    	return argument->derivative(i)[1];
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new CoordY_Vector(argument->clone())
        );
        return expr;
    }
};
Expression<double>::Ptr coord_y ( Expression<KDL::Vector>::Ptr a) {
    Expression<double>::Ptr expr(
        new CoordY_Vector(a)
    );
    return expr;
}

//CoordZ Vector
class CoordZ_Vector:
    public UnaryExpression<double, KDL::Vector>
{
public:
    typedef UnaryExpression<double, KDL::Vector> UnExpr;
public:
    CoordZ_Vector(){}
    CoordZ_Vector(
                const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("coord_z",arg)
                {}

    virtual double value() {
    	return argument->value()[2];
    }

    virtual double derivative(int i) {
    	return argument->derivative(i)[2];
    }

    virtual Expression<double>::Ptr derivativeExpression(int i);

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new CoordZ_Vector( argument->clone())
        );
        return expr;
    }
};
Expression<double>::Ptr coord_z ( Expression<KDL::Vector>::Ptr a) {
    Expression<double>::Ptr expr(
        new CoordZ_Vector(a)
    );
    return expr;
}


//Diff Vector Vector
class Diff_VectorVector:
	public BinaryExpression<KDL::Vector, KDL::Vector, KDL::Vector>
{
public:
	typedef BinaryExpression<KDL::Vector,KDL::Vector,KDL::Vector> BinExpr;
public:
	Diff_VectorVector(){}
	Diff_VectorVector(
			const  BinExpr::Argument1Expr::Ptr& arg1,
			const  BinExpr::Argument2Expr::Ptr& arg2):
				BinExpr("diff",arg1,arg2)
				{}

	virtual KDL::Vector value() {
		return diff(argument1->value(),argument2->value());
	}

	virtual KDL::Vector derivative(int i){
		return argument2->derivative(i) - argument1->derivative(i);
	}

    virtual Expression<Vector>::Ptr derivativeExpression(int i);

    virtual  BinExpr::Ptr clone() {
        Expression<KDL::Vector>::Ptr expr(
            new Diff_VectorVector( argument1->clone(), argument2->clone())
        );
        return expr;
    }
};

Expression<KDL::Vector>::Ptr diff ( Expression<KDL::Vector>::Ptr a1, Expression<KDL::Vector>::Ptr a2 ) {
	Expression<KDL::Vector>::Ptr expr(new Diff_VectorVector( a1, a2 ));
	return expr;
}

//=====================================================================================
Expression<Vector>::Ptr Vector_DoubleDoubleDouble::derivativeExpression(int i) {
    bool d1 = isDepOn<double>(i,argument1);
    bool d2 = isDepOn<double>(i,argument2);
    bool d3 = isDepOn<double>(i,argument3);
    if (!d1 && !d2 && !d3) {
        return Constant(Vector::Zero());
    } else if (!d1 && !d2) {
        return vector(Constant(0.0), Constant(0.0), argument3->derivativeExpression(i) );
    } else if (!d1 && !d3) {
        return vector(Constant(0.0), argument2->derivativeExpression(i), Constant(0.0) );
    } else if (!d2 && !d3) {
        return vector(argument1->derivativeExpression(i), Constant(0.0), Constant(0.0) );
    } else if (!d1) {
        return vector(Constant(0.0), argument2->derivativeExpression(i), argument3->derivativeExpression(i) );
    } else if (!d2) {
        return vector(argument1->derivativeExpression(i), Constant(0.0), argument3->derivativeExpression(i) );
    } else if (!d3) {
        return vector(argument1->derivativeExpression(i), argument2->derivativeExpression(i), Constant(0.0) );
    } else  {
        return vector(argument1->derivativeExpression(i), argument2->derivativeExpression(i), argument3->derivativeExpression(i) );
    }
}

Expression<double>::Ptr Dot_VectorVector::derivativeExpression(int i) {
       int nr = getDep<Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant<double>(0.0);
        } if (nr==2) {
            return dot(argument1,argument2->derivativeExpression(i));
        } if (nr==3) {
            return dot(argument1->derivativeExpression(i),argument2);
        } else {
            return dot(argument1,argument2->derivativeExpression(i)) + 
               dot(argument1->derivativeExpression(i),argument2);
        }
}

Expression<Vector>::Ptr CrossProduct_VectorVector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } if (nr==2) {
            return argument1 * argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i) * argument2;
        } else {
            return argument1 * argument2->derivativeExpression(i) + 
                argument1->derivativeExpression(i) * argument2;
        }
}

Expression<Vector>::Ptr Addition_VectorVector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } if (nr==2) {
            return argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i);
        } else {
            return argument1->derivativeExpression(i) + argument2->derivativeExpression(i);
        }
}

Expression<Vector>::Ptr Subtraction_VectorVector::derivativeExpression(int i) {
    return argument1->derivativeExpression(i) - argument2->derivativeExpression(i);
}

Expression<Vector>::Ptr Negate_Vector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } else {
            return -argument->derivativeExpression(i);
        }
}

Expression<double>::Ptr Norm_Vector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<Vector>::Ptr val = cached<Vector>( argument );
            return dot( val, argument->derivativeExpression(i)  ) / norm(val) ;
        }
}

Expression<double>::Ptr SquaredNorm_Vector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
            Expression<Vector>::Ptr val = cached<Vector>( argument );
            return dot( val, argument->derivativeExpression(i)  )*Constant<double>(2.0);
        }
}

Expression<Vector>::Ptr Multiplication_VectorDouble::derivativeExpression(int i) {
        int nr = getDep2<Vector,double>(i,argument1,argument2);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } if (nr==2) {
            return argument1 * argument2->derivativeExpression(i);
        } if (nr==3) {
            return argument1->derivativeExpression(i) * argument2;
        } else {
            return argument1 * argument2->derivativeExpression(i) + argument1->derivativeExpression(i) * argument2;
        }
}

Expression<double>::Ptr CoordX_Vector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
    	    return coord_x(argument->derivativeExpression(i));
        }
}

Expression<double>::Ptr CoordY_Vector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
    	    return coord_y(argument->derivativeExpression(i));
        }
}

Expression<double>::Ptr CoordZ_Vector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
    	    return coord_z(argument->derivativeExpression(i));
        }
}

Expression<Vector>::Ptr Diff_VectorVector::derivativeExpression(int i) {
        int nr = getDep<Vector>(i,argument1,argument2);
        if (nr==1) {
            return Constant<Vector>(Vector::Zero());
        } if (nr==2) {
            return argument2->derivativeExpression(i);
        } if (nr==3) {
            return - argument1->derivativeExpression(i);
        } else {
            return argument2->derivativeExpression(i) - argument1->derivativeExpression(i);
        }
}


} // end of namespace KDL
