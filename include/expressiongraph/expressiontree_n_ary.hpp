/*
 * expressiontree_motprof_trap.hpp
 *
 *  Created on: Oct 15, 2013
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

#ifndef KDL_EXPRESSIONTREE_MOTPROF_TRAP_HPP
#define KDL_EXPRESSIONTREE_MOTPROF_TRAP_HPP

#include <expressiongraph/expressiontree_expressions.hpp>
#include <algorithm>

namespace KDL {

template< typename ResultType, typename T>
class N_aryExpression: public Expression<ResultType> {
public:
    typedef Expression<T>                       ArgumentExpr;
    typename std::vector<typename ArgumentExpr::Ptr> arguments;



    N_aryExpression() {}

    N_aryExpression( const std::string& name,
                      const std::vector<typename ArgumentExpr::Ptr>& _arguments
                      ): Expression<ResultType>(name) 
    {
        arguments.resize(_arguments.size);
        for (unsigned int i=0;i<_arguments.size();++i) {
            arguments[i] = checkConstant<T>(_arguments[i]); 
        }
    }

    virtual void setInputValues(const std::vector<double>& values) {
        for (unsigned int i=0;i<_arguments.size();++i) {
            arguments[i]->setInputValues(values); 
        }
    }

    virtual void setInputValue(int variable_number, double val) {
        for (unsigned int i=0;i<_arguments.size();++i) {
            arguments[i]->setInputValue(variable_number,val); 
        }
    }

    virtual void setInputValue(int variable_number, const Rotation& val) {
        for (unsigned int i=0;i<_arguments.size();++i) {
            arguments[i]->setInputValue(variable_number,val); 
        }
    }

    virtual int number_of_derivatives() {
        int n=0;
        for (unsigned int i=0;i<_arguments.size();++i) {
            n = std::max( n, arguments[i]->number_of_derivatives() ); 
        }
        return n;
 
    } 

    virtual typename Expression<Frame>::Ptr subExpression_Frame(const std::string& name) {
        typename Expression<Frame>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Frame(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    virtual typename Expression<Rotation>::Ptr subExpression_Rotation(const std::string& name) {
        typename Expression<Rotation>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Rotation(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    virtual typename Expression<Vector>::Ptr subExpression_Vector(const std::string& name) {
        typename Expression<Vector>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Vector(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    virtual typename Expression<Wrench>::Ptr subExpression_Wrench(const std::string& name) {
        typename Expression<Wrench>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Wrench(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    virtual typename Expression<Twist>::Ptr subExpression_Twist(const std::string& name) {
        typename Expression<Twist>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Twist(name);
            if (a) {
                return a;
            }
        }
        return a;

    }

    virtual typename Expression<double>::Ptr subExpression_Double(const std::string& name) {
        typename Expression<double>::Ptr a;
        for (unsigned int i=0;i<arguments.size();++i) {
            a = arguments[i]->subExpression_Double(name);
            if (a) {
                return a;
            }
        }
        return a;
    }

    virtual void addToOptimizer(ExpressionOptimizer& opt) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->addToOptimizer(opt);
        }
    }


    virtual void getDependencies(std::set<int>& varset) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->getDependencies(varset);
        }
    }

    virtual void getScalarDependencies(std::set<int>& varset) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->getScalarDependencies(varset);
        }
    }

    virtual void getRotDependencies(std::set<int>& varset) {
        for (unsigned int i=0;i<arguments.size();++i) {
            arguments[i]->getRotDependencies(varset);
        }
    }

    virtual void debug_printtree() {
        std::cout << Expression<ResultType>::name << "(";
        for (unsigned int i=0;i<arguments.size();++i) {
            if (i!=0) std::cout << ",";
            arguments[i]->debug_printtree();
        }
        std::cout << ")";
    }
    virtual void print(std::ostream& os) const {
        os << "" << Expression<ResultType>::name << "(";
        for (unsigned int i=0;i<arguments.size();++i) {
            if (i!=0) os << ",";
            arguments[i]->print(os);
        }
        os << ")";
    }

    virtual void write_dotfile_update(std::ostream& of, pnumber& thisnode) {
        thisnode=(size_t)this;
        of << "S"<<thisnode<<"[label=\"" << Expression<ResultType>::name << "\",shape=box,style=filled,fillcolor="
           << COLOR_OPERATION << ",color=black]\n";
        std::vector<pnumber> argnode( arguments.size());
        for (unsigned int=0;i<arguments.size();++i) {
            arguments[i]->write_dotfile_update(of,argnode[i]);
        }
        for (unsigned int=0;i<arguments.size();++i) {
            of << "S"<<thisnode<<" -> " << "S"<<argnode[i] << "\n"; //rev
        }
    }
    virtual void write_dotfile_init() {
        for (unsigned int=0;i<arguments.size();++i) {
            arguments[i]->write_dotfile_init();
        }
    } 
};

} // end of namespace KDL
#endif
