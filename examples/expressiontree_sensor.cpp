/*
 *  \file expressiontree_sensor.cpp
 *
 *  Created on: Dec., 2012
 *      Author: Erwin Aertbelien 
***************************************************************************/

#include <kdl/expressiontree.hpp>
#include <fstream>

#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace KDL;




class NoiseGen {
    boost::mt19937 rng;
    boost::normal_distribution<> nd;
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor; 
public:
    NoiseGen():nd(0,0.03), var_nor(rng,nd) {}
    double operator()() {
        return var_nor();
    }
};


/**
 * You can construct a expression graph node for a sensor
 * by definining a unary eg node and caching it and making
 * it depended on time.  Then, each time the time variable changes
 * (and only then) the expression is read out again...
 */
class Sensor:
    public UnaryExpression<double, double>
{
public:
    typedef UnaryExpression<double, double> UnExpr;
    double val;
    NoiseGen noise;
public:
    Sensor() {}

    Sensor(const  UnExpr::ArgumentExpr::Ptr& arg):
                UnExpr("sensor",arg)
        {}

    virtual double value() {
        val = 12.3+noise(); // read out sensor 
        return val;
    }

    virtual double derivative(int i) {
        return 0; 
    }

    virtual Expression<double>::Ptr derivativeExpression(int i) {
        return Constant(0.0);
    }

    virtual  UnExpr::Ptr clone() {
        Expression<double>::Ptr expr(
            new Sensor( argument->clone())
        );
        return expr;
    }
};

inline Expression<double>::Ptr getSensor() {
    // this assumes that the time variable has variable index 0.
    Expression<double>::Ptr expr(
        new Sensor( input(0) )
    );
    return cached<double>(expr);
}


/*
 * Example that demonstrates how to encapsulate a sensor.
 */
int main(int argc, char* argv[]) {
	using namespace KDL;
	using namespace std;

    Expression<double>::Ptr e = getSensor();
    std::cout << e->value() << std::endl;
    std::cout << e->value() << std::endl;
    e->setInputValue(0,0.1);
    std::cout << e->value() << std::endl;
    std::cout << e->value() << std::endl;
    e->setInputValue(0,0.2);
    std::cout << e->value() << std::endl;
}

