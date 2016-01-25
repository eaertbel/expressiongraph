/*
 * expressiontree_motprofiles.hpp
 *
 *  Created on: Jan, 2016
 *      Author: Erwin Aertbelien
 *
* expressiongraph library
* 
* Copyright 2016 Erwin Aertbelien - KU Leuven - Dep. of Mechanical Engineering
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

#ifndef KDL_EXPRESSIONTREE_MOTPROFILES_HPP
#define KDL_EXPRESSIONTREE_MOTPROFILES_HPP

#include <kdl/expressiontree_expressions.hpp>
#include <kdl/expressiontree_mimo.hpp>
#include <kdl/mptrap.hpp>
#include <stdexcept>

namespace KDL {

class MotionProfileTrapezoidal : 
    public MIMO {
public:
    int critical_output;  // index of the output that determines the duration
    std::vector< MPTrap > mp;
    double progrvar_value;
    typedef boost::shared_ptr<MotionProfileTrapezoidal> Ptr;
    MotionProfileTrapezoidal();

    /**
     * \brief declares an expression for the input variable for this motionprofile
     *
     * i.e. the progression variable, typically time or a time
     * related parameter
     */
    void setProgressExpression(const Expression<double>::Ptr& s); 

    /**
     * \brief gets an expression for the progression variable
     */
    Expression<double>::Ptr getProgressExpression();

    /**
     * \brief declares a given output motionprofile
     * All output motionprofiles will be synchronized, i.e. will start and end
     * at the same time.
     * \param startv : starting value expression
     * \param endv:    ending value expression
     * \param maxvel:  maximum velocity for this output
     * \param maxacc:  maximum acceleration for this output.
     */
    void addOutput( const Expression<double>::Ptr& startv, 
                    const Expression<double>::Ptr& endv, 
                    const Expression<double>::Ptr& maxvel, 
                    const Expression<double>::Ptr& maxacc);
    /**
     * \brief returns the number of declared outputs
     */
    int nrOfOutputs();
    
    /** 
     * \brief gets the starting value expression for output idx
     * \param idx output (from 0 to nrOfOutputs()-1 )
     */
    Expression<double>::Ptr getStartValue(int idx);

    /** 
     * \brief gets the ending value expression for output idx
     * \param idx output (from 0 to nrOfOutputs()-1 )
     */
    Expression<double>::Ptr getEndValue(int idx);

    /** 
     * \brief gets the maximum velocity  expression for output idx
     * \param idx output (from 0 to nrOfOutputs()-1 )
     */
    Expression<double>::Ptr getMaxVelocity(int idx);

    /** 
     * \brief gets the maximum acceleration expression for output idx
     * \param idx output (from 0 to nrOfOutputs()-1 )
     */
    Expression<double>::Ptr getMaxAcceleration(int idx);

    void compute();

    virtual MIMO::Ptr clone();
};


class MotionProfileTrapezoidalOutput : public MIMO_Output<double> {
    public:
        typedef boost::shared_ptr<MotionProfileTrapezoidalOutput> Ptr;
        int outputnr;
        int idx_base;

        MotionProfileTrapezoidalOutput(MIMO::Ptr m, int _outputnr);
        double value();
        double derivative(int i);
        MIMO_Output<double>::Ptr clone(); 
};

/**
 * \brief creates a MotionProfileTrapezoidal object
 * 
 * You will need to call a set of operations on this object:
 *  - setProgressExpression(...) to set the "time"-like variable in the motionprofile
 *  - addOutput( startval, endval, maxvel, maxacc) : all inputs are expressions.  The influence of changing maxvel
 *      and maxacc will not be propagated into the computations of the derivative.
 * 
 * After this, you can call get_output_profile(...) to get an expression for the different outputs you had
 * defined.      
 */
MotionProfileTrapezoidal::Ptr create_motionprofile_trapezoidal() {
    return boost::make_shared< MotionProfileTrapezoidal>();
}

/**
 * \brief gets an expression representing the motion profile for a given output
 * \param idx index of the output for which the expression is returned.
 */
Expression<double>::Ptr get_output_profile(MotionProfileTrapezoidal::Ptr& m,int output) {
    return boost::make_shared<MotionProfileTrapezoidalOutput>( m,output);
}

/**
 * \brief gets an expression representing the duration 
 */
Expression<double>::Ptr get_duration(MotionProfileTrapezoidal::Ptr& m) {
    return boost::make_shared<MotionProfileTrapezoidalOutput>( m,-1);
}





} // end of namespace KDL
#endif
