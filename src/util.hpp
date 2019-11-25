/*
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

#ifndef KDL_EXPRESSIONGRAPH_UTIL_HPP
#define KDL_EXPRESSIONGRAPH_UTIL_HPP

// private utilities for constant collapsing:

/*
        int nr = getDep<double>(i,argument1,argument2);
        if (nr==1) {
            return Constant<double>(0.0);
        } if (nr==2) {
        } if (nr==3) {
        } else {
        }

        int nr = getDep<double>(i,argument);
        if (nr==1) {
            return Constant<double>(0.0);
        } else {
        }


*/

namespace KDL {

template <typename T>
inline int getDep( int i, typename boost::shared_ptr<Expression<T> > argument1, typename boost::shared_ptr<Expression<T> > argument2) {
        //std::cout << "multiplication between " << std::endl;
        //argument1->print(std::cout);std::cout << std::endl;
        //argument2->print(std::cout);std::cout << std::endl;
        std::set<int> vset1;
        argument1->getDependencies(vset1);
        bool depend1 = vset1.count(i)>0;
        //std::cout << depend1 << std::endl;
        std::set<int> vset2;
        argument2->getDependencies(vset2);
        bool depend2 = vset2.count(i)>0;
        //std::cout << depend2 << std::endl;
        if (!depend1 && !depend2) {
            return 1;
        } else if (!depend1) {
            return 2;
        } else if (!depend2) {
            return 3;
        } else {
            return 4;
        }
}
template <typename T1,typename T2>
inline int getDep2( int i, 
              boost::shared_ptr<Expression<T1> > argument1, 
              boost::shared_ptr<Expression<T2> > argument2) {
        //std::cout << "multiplication between " << std::endl;
        //argument1->print(std::cout);std::cout << std::endl;
        //argument2->print(std::cout);std::cout << std::endl;
        std::set<int> vset1;
        argument1->getDependencies(vset1);
        bool depend1 = vset1.count(i)>0;
        //std::cout << depend1 << std::endl;
        std::set<int> vset2;
        argument2->getDependencies(vset2);
        bool depend2 = vset2.count(i)>0;
        //std::cout << depend2 << std::endl;
        if (!depend1 && !depend2) {
            return 1;
        } else if (!depend1) {
            return 2;
        } else if (!depend2) {
            return 3;
        } else {
            return 4;
        }
}
template<typename T>
bool isDepOn(int i, boost::shared_ptr<Expression<T> > a) {
        std::set<int> vset;
        a->getDependencies(vset);
        return vset.count(i)>0;
}

template <typename T>
inline int getDep( int i, boost::shared_ptr<Expression<T> > argument1) {
        std::set<int> vset1;
        argument1->getDependencies(vset1);
        bool depend1 = vset1.count(i)>0;
        if (!depend1) {
            return 1;
        } else {
            return 2;
        }
}


} // namespace KDL

#endif
