#ifndef INTEGRATOR_HPP_51566a22_0d16_4ad0_b673_fa9a6d3cecec
#define INTEGRATOR_HPP_51566a22_0d16_4ad0_b673_fa9a6d3cecec

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
#include "expressiontree_function.hpp"
#include "expressiontree_var.hpp"
#include <string> 
#include <map>
#include <vector>
#include <functional>
#include <exception>


namespace KDL { 
        const int err_EDOM=1;
        const int err_ERANGE=2;

    using IntegrandType=std::function< void(double,double*) >;

    class IntegratorAdaptiveSimpson {
            double epsilon;
            int        maxRecDepth;
            int        diffRecDepth;
            int        err;

            int        vecsize;  // size of the vectors
            /** storage of temporaries **/
            double*    storage;
            int        currentp;       
            double     fa_value;
            double     fb_value;

            // allocation statistics:
            int maxalloc;  
            int alloclimit;

            double* alloc() {
                //std::cout << "alloc " << currentp << " limit " << alloclimit << std::endl;
                double* p = &storage[currentp];
                currentp+=vecsize;
                maxalloc=std::max( maxalloc, currentp);
                assert(currentp<alloclimit);
                return p;
            }
            void pop(int n=1) {
                currentp-=vecsize*n;
                assert(currentp>=0);
            }


            double* adaptiveSimpsonsAux(IntegrandType f, double a, double b, double eps,
                                      double* whole, double* fa, double* fb, double* fm, int rec);
        public:
            IntegratorAdaptiveSimpson( double _epsilon, int _minRecDepth, int _maxRecDepth, int _vecsize);
            int get_err();
            int get_maxalloc() { return maxalloc;}
            double get_fa_value() { return fa_value;}
            double get_fb_value() { return fb_value;}
            int integrate( IntegrandType f, double a,double b, double* result);
            ~IntegratorAdaptiveSimpson();
    };


    // we need FunctionEvaluation exposed for setInputValues
    //IntegrandType make_function( FunctionDefinition::Ptr func, const std::vector<int>& ndx ) {
    IntegrandType make_function( FunctionDefinition::Ptr func );

    // we need FunctionEvaluation exposed for setInputValues
    //IntegrandType make_function( FunctionDefinition::Ptr func, const std::vector<int>& ndx ) {
    IntegrandType make_function( FunctionDefinition::Ptr func, int& counter );

}// namespace KDL

#endif
 
