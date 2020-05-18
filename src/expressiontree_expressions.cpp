/**
 * \file expressiontree_expressions.cpp
 * \brief example of using ExpressionTree in combination with KDL::Chain.
 *
 * \Author: Jan. 2013, Erwin Aertbelien 
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




#include <expressiongraph/expressiontree_expressions.hpp>
#include <algorithm>
#include <iterator>

using namespace std;

namespace KDL {

double EG_EPS_EQUAL = 1E-16;

void ExpressionOptimizer::prepare(const std::vector<int>& inputvarnr, const std::vector<int>& rotinputvarnr) {
    inputset.clear();   

    inputs.resize(inputvarnr.size());
    this->inputvarnr.resize(inputvarnr.size());
    for (size_t i=0;i<inputvarnr.size();++i) {
        this->inputvarnr[i] = inputvarnr[i];
        inputset.insert(inputvarnr[i]);
    }

    rotinputs.resize(rotinputvarnr.size());
    this->rotinputvarnr.resize(rotinputvarnr.size());
    for (size_t i=0;i<rotinputvarnr.size();++i) {
        this->rotinputvarnr[i] = rotinputvarnr[i];
        inputset.insert(rotinputvarnr[i]);
        inputset.insert(rotinputvarnr[i]+1);
        inputset.insert(rotinputvarnr[i]+2);
    }

    cached.clear();
    v_cached.clear();
}


void ExpressionOptimizer::prepare(const std::vector<int>& inputvarnr) {
    std::vector<int> rotinputvarnr(0);
    prepare(inputvarnr,rotinputvarnr);
}


void ExpressionOptimizer::addInput(InputType* obj) {
    for (size_t i=0;i<inputvarnr.size();++i) {
        if (inputvarnr[i]==obj->variable_number) {
            inputs[i].push_front(obj);
            //cout << "input " << obj->variable_number << " added \n";
            break;
        }
    }
}

void ExpressionOptimizer::addInput(InputRotationType* obj) {
    for (size_t i=0;i<rotinputvarnr.size();++i) {
        if (rotinputvarnr[i]==obj->variable_number) {
            rotinputs[i].push_front(obj);
            //cout << "input " << obj->variable_number << " added \n";
            break;
        }
    }
}

void ExpressionOptimizer::addCached(CachedExpression* obj) {
    InputSet dependency;
    obj->getDependencies( dependency );
    //cout << "cached dependency: ";
    //copy(dependency.begin(), dependency.end(), ostream_iterator<int>(cout, " "));
    //cout << "\n";
    for (InputSet::iterator it=dependency.begin();it!=dependency.end();++it) {
        if (inputset.find( *it )!= inputset.end() ) {
            if (cached.find(obj)==cached.end()) {
                cached.insert( obj );  
                v_cached.push_back( obj );
            }
            //cout << "cached " << obj->cached_name << " added \n";
            break;
        }
    }
}

void ExpressionOptimizer::setInputValues(const std::vector<double>& values) {
    assert( values.size() == inputvarnr.size() );
    for (size_t i=0; i< inputvarnr.size(); ++i) {
        //cout << "nr of input objects for "<< i << " variable in the list " << inputs[i].size() << endl;
        double value = values[i];
        for (ListInput::iterator it=inputs[i].begin(); it != inputs[i].end(); ++it ) {
            //cout << "set input value \n"<< endl;
            (*it)->val = value;
        }
    }

    // invalidate the appropriate caches:
    for (std::vector<CachedExpression*>::iterator   it=v_cached.begin();it!=v_cached.end();++it)  
        (*it)->invalidate_cache();
}

void ExpressionOptimizer::setInputValues(const std::vector<double>& values, const std::vector<Rotation>& rotvalues) {
    assert( values.size() == inputvarnr.size() );
    assert( rotvalues.size() == rotinputvarnr.size() );

    for (size_t i=0; i< inputvarnr.size(); ++i) {
        double value = values[i];
        //cout << "nr of input objects for "<< i << " variable in the list " << inputs[i].size() << endl;
        for (ListInput::iterator it=inputs[i].begin(); it != inputs[i].end(); ++it ) {
            //cout << "set input value \n"<< endl;
            (*it)->val = value;
        }
    }

    for (size_t i=0; i< rotinputvarnr.size(); ++i) {
        //cout << "nr of input objects for "<< i << " variable in the list " << inputs[i].size() << endl;
        const Rotation& value = rotvalues[i];
        for (ListRotInput::iterator it=rotinputs[i].begin(); it != rotinputs[i].end(); ++it ) {
            //cout << "set input value \n"<< endl;
            (*it)->val = value;
        }
    }

    // invalidate the appropriate caches:
    for (std::vector<CachedExpression*>::iterator   it=v_cached.begin();it!=v_cached.end();++it)  
        (*it)->invalidate_cache();
}


void ExpressionOptimizer::setInputValues(const Eigen::VectorXd& values) {
    assert( values.rows() == (int)inputvarnr.size() );
    for (size_t i=0; i< inputvarnr.size(); ++i) {
        double value = values[i];
        //cout << "nr of input objects for "<< i << " variable in the list " << inputs[i].size() << endl;
        for (ListInput::iterator it=inputs[i].begin(); it != inputs[i].end(); ++it ) {
            //cout << "set input value \n"<< endl;
            (*it)->val = value;
        }
    }

    // invalidate the appropriate caches:
    for (std::vector<CachedExpression*>::iterator   it=v_cached.begin();it!=v_cached.end();++it)  
        (*it)->invalidate_cache();
}


void ExpressionOptimizer::setInputValues(const Eigen::VectorXd& values, const std::vector<Rotation>& rotvalues) {
    assert( (size_t) values.size() == (size_t) inputvarnr.size() );
    assert( (size_t) rotvalues.size() == (size_t) rotinputvarnr.size() );

    for (size_t i=0; i< inputvarnr.size(); ++i) {
        double value = values[i];
        //cout << "nr of input objects for "<< i << " variable in the list " << inputs[i].size() << endl;
        for (ListInput::iterator it=inputs[i].begin(); it != inputs[i].end(); ++it ) {
            //cout << "set input value \n"<< endl;
            (*it)->val = value;
        }
    }

    for (size_t i=0; i< rotinputvarnr.size(); ++i) {
        const Rotation& value = rotvalues[i];
        //cout << "nr of input objects for "<< i << " variable in the list " << inputs[i].size() << endl;
        for (ListRotInput::iterator it=rotinputs[i].begin(); it != rotinputs[i].end(); ++it ) {
            //cout << "set input value \n"<< endl;
            (*it)->val = value; 
        }
    }


    // invalidate the appropriate caches:
    for (std::vector<CachedExpression*>::iterator   it=v_cached.begin();it!=v_cached.end();++it)  
        (*it)->invalidate_cache();
}





template class Expression<double>;
template class Expression<Vector>;
template class Expression<Frame>;
template class Expression<Rotation>;
template class Expression<Twist>;
template class Expression<Wrench>;
template class Expression<Quaternion>;


template class UnaryExpression<double,double>;
template class UnaryExpression<Vector,Vector>;
template class UnaryExpression<Frame,Frame>;
template class UnaryExpression<Rotation,Rotation>;
template class UnaryExpression<Twist,Twist>;
template class UnaryExpression<Wrench,Wrench>;
template class UnaryExpression<Quaternion,Quaternion>;
template class UnaryExpression<double,Quaternion>;
template class UnaryExpression<Vector,Quaternion>;
template class UnaryExpression<Quaternion,Vector>;
template class UnaryExpression<Quaternion,Rotation>;
template class UnaryExpression<Rotation,Quaternion>;

template class BinaryExpression<double,double,double>;
template class BinaryExpression<Vector,Vector,Vector>;
template class BinaryExpression<Rotation,Rotation,Rotation>;
template class BinaryExpression<Frame,Frame,Frame>;
template class BinaryExpression<Twist,Twist,Twist>;
template class BinaryExpression<Wrench,Wrench,Wrench>;
template class BinaryExpression<Quaternion,Quaternion,Quaternion>;

template class BinaryExpression<Quaternion,Quaternion,Vector>;
template class BinaryExpression<Quaternion,Vector,Quaternion>;
template class BinaryExpression<double,Quaternion,Quaternion>;
template class BinaryExpression<Vector,Quaternion,Vector>;
template class BinaryExpression<Quaternion,double,Vector>;



template class TernaryExpression<double,double,double,double>;
template class TernaryExpression<Vector,double,double,double>;
template class TernaryExpression<Rotation,double,double,double>;

template class FunctionType<double>;
template class FunctionType<Vector>;
template class FunctionType<Rotation>;
template class FunctionType<Frame>;
template class FunctionType<Twist>;
template class FunctionType<Wrench>;
template class FunctionType<Quaternion>;

template class ConstantType<double>;
template class ConstantType<Vector>;
template class ConstantType<Rotation>;
template class ConstantType<Frame>;
template class ConstantType<Twist>;
template class ConstantType<Wrench>;
template class ConstantType<Quaternion>;


template class CachedType<double>;
template class CachedType<Vector>;
template class CachedType<Rotation>;
template class CachedType<Frame>;
template class CachedType<Twist>;
template class CachedType<Wrench>;
template class CachedType<Quaternion>;


template class MakeConstantType<double>;
template class MakeConstantType<Vector>;
template class MakeConstantType<Rotation>;
template class MakeConstantType<Frame>;
template class MakeConstantType<Twist>;
template class MakeConstantType<Wrench>;
template class MakeConstantType<Quaternion>;

};// namespace KDL
