 /**
 * \file cloning.cpp
 * \brief example of cloning and Variable expressions as parameters
 *
 * \Author: Sept. 2016, Erwin Aertbelien 
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



/**
 *
 * This file demonstrates a technique that combines the use of Variable expressions and
 * cloning. The Variable expressions are cloned and are de facto constants after the clone
 * (because there is no access any more the change their value).  This behaviour can be exploited
 * for e.g. in fitting applications. 
 *
**/

#include <expressiongraph/expressiontree.hpp>
#include <expressiongraph/expressiontree_var.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <Eigen/SVD>
 
using namespace KDL;

void solve( std::vector<Expression<double>::Ptr> list,  int time_ndx, const std::vector<int>& ndx, Eigen::VectorXd& solution) {
    using namespace std;
    int m = list.size();
    int n = ndx.size();
    Eigen::VectorXd err(m);
    Eigen::MatrixXd Jac(m, n );
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, n );
    Eigen::VectorXd  d_sol(n);

    for (int it=0;it<10;++it) {
        for (int i=0;i<m;++i) {
            list[i]->setInputValues(ndx,solution);
            list[i]->setInputValue(time_ndx,0); // to force the variable expr to be read (that depend on time_ndx).
            err(i) = list[i]->value();
            for (int j=0;j<n;++j) {
                Jac(i,j) = list[i]->derivative(ndx[j]);
            }
        }
        //cout << "error " << err.transpose() << endl;
        //cout << "Jacobian " << endl << Jac << endl;
        svd.compute(Jac,Eigen::ComputeThinU| Eigen::ComputeThinV );
        d_sol = svd.solve(err);
        solution = solution - d_sol;
        cout << "solution at it="<<it<<"\t\t"<<solution.transpose() << endl;
    }
}


Expression<double>::Ptr  addMeasurement( 
            Expression<double>::Ptr model, VariableType<double>::Ptr& x, VariableType<double>::Ptr& y, 
            double x_val, double y_val ) {
    x->val = x_val;
    y->val = y_val;
    return model->clone();
}


int main() {
    using namespace std;
    // Building up the measurement model in function of a set of unknowns and some parameters x and y:
    VariableType<double>::Ptr x = Variable<double>(1,1);
    VariableType<double>::Ptr y = Variable<double>(1,1);
    Expression<double>::Ptr u1 = input(2);
    Expression<double>::Ptr u2 = input(3);  
    Expression<double>::Ptr model =  u1*x+u2 - y;

    // Enter measurements:
    std::vector<Expression<double>::Ptr> list;
    list.push_back( addMeasurement(model,x,y, 1, 2));
    list.push_back( addMeasurement(model,x,y, 1.5, 2.2));
    list.push_back( addMeasurement(model,x,y, 1.8, 2.3));
    list.push_back( addMeasurement(model,x,y, 2, 2.5));


    // Solve:
    std::vector<int> ndx(2); 
    ndx[0]=2;ndx[1]=3;
    Eigen::VectorXd solution(2);
    solution(0) = 1.4;
    solution(1) = 0;
    solve( list, 1, ndx, solution);

    return 0;
}

