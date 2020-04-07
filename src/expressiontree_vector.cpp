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
