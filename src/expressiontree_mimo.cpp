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

#include <kdl/expressiontree_mimo.hpp>
namespace KDL {

MIMO::MIMO() {}

MIMO::MIMO( const std::string& _name):
            name(_name),
            cached(false)
{}


void MIMO::setInputValues(const std::vector<double>& values) {
        for (size_t i=0;i<inputDouble.size();++i) {
            inputDouble[i]->setInputValues(values);
        }
        for (size_t i=0;i<inputFrame.size();++i) {
            inputFrame[i]->setInputValues(values);
        }
        for (size_t i=0;i<inputTwist.size();++i) {
            inputTwist[i]->setInputValues(values);
        }
        cached = false;
    }

void MIMO::setInputValue(int variable_number, double val) {
        for (size_t i=0;i<inputDouble.size();++i) {
            inputDouble[i]->setInputValue(variable_number,val);
        }
        for (size_t i=0;i<inputFrame.size();++i) {
            inputFrame[i]->setInputValue(variable_number,val);
        }
        for (size_t i=0;i<inputTwist.size();++i) {
            inputTwist[i]->setInputValue(variable_number,val);
        }
        cached = false;
}

void MIMO::setInputValue(int variable_number, const Rotation& val) {
        for (size_t i=0;i<inputDouble.size();++i) {
            inputDouble[i]->setInputValue(variable_number,val);
        }
        for (size_t i=0;i<inputFrame.size();++i) {
            inputFrame[i]->setInputValue(variable_number,val);
        }
        for (size_t i=0;i<inputTwist.size();++i) {
            inputTwist[i]->setInputValue(variable_number,val);
        }
        cached = false;
}

int MIMO::number_of_derivatives() {
        int n=0;
        for (size_t i=0;i<inputDouble.size();++i) {
            int ni = inputDouble[i]->number_of_derivatives();
            n = n > ni ? n : ni;
        }
        for (size_t i=0;i<inputFrame.size();++i) {
            int ni = inputFrame[i]->number_of_derivatives();
            n = n > ni ? n : ni;
        }
        for (size_t i=0;i<inputTwist.size();++i) {
            int ni = inputTwist[i]->number_of_derivatives();
            n = n > ni ? n : ni;
        }
        return n;
} 

Expression<Frame>::Ptr MIMO::subExpression_Frame(const std::string& name) {
    for (size_t i=0;i<inputDouble.size();++i) {
        Expression<Frame>::Ptr tmp;
        tmp = inputDouble[i]->subExpression_Frame(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        Expression<Frame>::Ptr tmp;
        tmp = inputFrame[i]->subExpression_Frame(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        Expression<Frame>::Ptr tmp;
        tmp = inputTwist[i]->subExpression_Frame(name);
        if (tmp) {
            return tmp;
        }
    }
    return Expression<Frame>::Ptr();
}
Expression<Rotation>::Ptr MIMO::subExpression_Rotation(const std::string& name) {
    for (size_t i=0;i<inputDouble.size();++i) {
        Expression<Rotation>::Ptr tmp;
        tmp = inputDouble[i]->subExpression_Rotation(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        Expression<Rotation>::Ptr tmp;
        tmp = inputFrame[i]->subExpression_Rotation(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        Expression<Rotation>::Ptr tmp;
        tmp = inputTwist[i]->subExpression_Rotation(name);
        if (tmp) {
            return tmp;
        }
    }
    return Expression<Rotation>::Ptr();
}
Expression<Vector>::Ptr MIMO::subExpression_Vector(const std::string& name) {
    for (size_t i=0;i<inputDouble.size();++i) {
        Expression<Vector>::Ptr tmp;
        tmp = inputDouble[i]->subExpression_Vector(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        Expression<Vector>::Ptr tmp;
        tmp = inputFrame[i]->subExpression_Vector(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        Expression<Vector>::Ptr tmp;
        tmp = inputTwist[i]->subExpression_Vector(name);
        if (tmp) {
            return tmp;
        }
    }
 
    return Expression<Vector>::Ptr();
}
Expression<Wrench>::Ptr MIMO::subExpression_Wrench(const std::string& name) {
    for (size_t i=0;i<inputDouble.size();++i) {
        Expression<Wrench>::Ptr tmp;
        tmp = inputDouble[i]->subExpression_Wrench(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        Expression<Wrench>::Ptr tmp;
        tmp = inputFrame[i]->subExpression_Wrench(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        Expression<Wrench>::Ptr tmp;
        tmp = inputTwist[i]->subExpression_Wrench(name);
        if (tmp) {
            return tmp;
        }
    }
    return Expression<Wrench>::Ptr();
}
Expression<double>::Ptr MIMO::subExpression_Double(const std::string& name) {
    for (size_t i=0;i<inputDouble.size();++i) {
        Expression<double>::Ptr tmp;
        tmp = inputDouble[i]->subExpression_Double(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        Expression<double>::Ptr tmp;
        tmp = inputFrame[i]->subExpression_Double(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        Expression<double>::Ptr tmp;
        tmp = inputTwist[i]->subExpression_Double(name);
        if (tmp) {
            return tmp;
        }
    }
    return Expression<double>::Ptr();
}


Expression<Twist>::Ptr MIMO::subExpression_Twist(const std::string& name) {
    for (size_t i=0;i<inputDouble.size();++i) {
        Expression<Twist>::Ptr tmp;
        tmp = inputDouble[i]->subExpression_Twist(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        Expression<Twist>::Ptr tmp;
        tmp = inputFrame[i]->subExpression_Twist(name);
        if (tmp) {
            return tmp;
        }
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        Expression<Twist>::Ptr tmp;
        tmp = inputTwist[i]->subExpression_Twist(name);
        if (tmp) {
            return tmp;
        }
    }
    return Expression<Twist>::Ptr();
}

void MIMO::addToOptimizer(ExpressionOptimizer& opt) {
    CachedExpression::addToOptimizer(opt);
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->addToOptimizer(opt);
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->addToOptimizer(opt);
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->addToOptimizer(opt);
    }
}

void MIMO::getDependencies(std::set<int>& varset) {
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->getDependencies(varset);
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->getDependencies(varset);
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->getDependencies(varset);
    }
}

void MIMO::getScalarDependencies(std::set<int>& varset) {
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->getScalarDependencies(varset);
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->getScalarDependencies(varset);
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->getScalarDependencies(varset);
    }
}

void MIMO::getRotDependencies(std::set<int>& varset) {
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->getRotDependencies(varset);
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->getRotDependencies(varset);
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->getRotDependencies(varset);
    }
}

void MIMO::invalidate_cache() {
    cached=false;
}

void MIMO::debug_printtree() {
    std::cout << name << "(";
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->debug_printtree();
        std::cout << ",";
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->debug_printtree();
        std::cout << ",";
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->debug_printtree();
        std::cout << ",";
    }
    std::cout << ")";
}

void MIMO::print(std::ostream& os) const {
    os << name << "(";
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->debug_printtree();
        os << ",";
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->debug_printtree();
        os << ",";
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->debug_printtree();
        os << ",";
    }
    os << ")";
}

void MIMO::write_dotfile_update(std::ostream& of, pnumber& thisnode) {
    if (dot_already_written) {
        thisnode=(size_t)this;
        return;
    }
    dot_already_written=true;
    thisnode=(size_t)this;
    of << "S"<<thisnode<<"[label=\"" << name << "\",shape=box,style=filled,fillcolor="
       << COLOR_OPERATION << ",color=black]\n";
    pnumber argnode;
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->write_dotfile_update(of,argnode);
        of << "S"<<thisnode<<" -> " << "S"<<argnode<< "\n"; 
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->write_dotfile_update(of,argnode);
        of << "S"<<thisnode<<" -> " << "S"<<argnode<< "\n"; 
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->write_dotfile_update(of,argnode);
        of << "S"<<thisnode<<" -> " << "S"<<argnode<< "\n"; 
    }
}

void MIMO::write_dotfile_init() {
    dot_already_written=false;
    for (size_t i=0;i<inputDouble.size();++i) {
        inputDouble[i]->write_dotfile_init();
    }
    for (size_t i=0;i<inputFrame.size();++i) {
        inputFrame[i]->write_dotfile_init();
    }
    for (size_t i=0;i<inputTwist.size();++i) {
        inputTwist[i]->write_dotfile_init();
    }
} 

MIMO::Ptr MIMO::getClone(int count) {
    if (count >= (int)queue_of_clones.size()) {
        queue_of_clones.resize(count+1);
    }
    MIMO::Ptr tmp = queue_of_clones[count].lock();
    if (!tmp) {
        tmp =  this->clone();
        queue_of_clones[count] = tmp;
    }
    return tmp; 
}  

MIMO::~MIMO() {
}


} // namespace KDL

