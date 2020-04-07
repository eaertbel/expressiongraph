/**
 * \file expressiontree_example6.cpp
 * \brief example of using ExpressionTree in combination with KDL::Chain.
 *
 * \Author: Sept. 2012, Erwin Aertbelien 
 *
***************************************************************************/


#include <expressiongraph/expressiontree.hpp>
#include <vector>

int main(int argc, char* argv[]) {
 	using namespace KDL;
	using namespace std;
    double L1=0.310;
    double L2=0.400;
    double L3=0.390;
    double L4=0.078;

    Chain chain;
    chain.addSegment(Segment("Segment 0", Joint("Joint 0", Joint::RotZ),Frame(Vector(0,0,L1))));
    chain.addSegment(Segment("Segment 1", Joint("Joint 1", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 2", Joint("Joint 2", Joint::RotZ),Frame(Vector(0,0,L2))));
    chain.addSegment(Segment("Segment 3", Joint("Joint 3", Joint::RotX),Frame(Vector(0,0,0))));
    chain.addSegment(Segment("Segment 4", Joint("Joint 4", Joint::RotZ),Frame(Vector(0,0,L3))));
    chain.addSegment(Segment("Segment 5", Joint("Joint 5", Joint::RotX),Frame(Vector(0,0,L4))));
    chain.addSegment(Segment("Segment 6", Joint("Joint 6", Joint::RotZ),Frame(Vector(0,0,0))));

    int robotvar;
    robotvar=0;
    std::string base_name="robot1";
    Expression<Frame>::Ptr robot1 = cached<Frame>(
        base_name +":ee",
        cached<Frame>(base_name+":wrist",
            cached<Frame>(base_name+":elbow",
                cached<Frame>(base_name+":shoulder",
                    frame( rot_z(input(robotvar)),  Constant(Vector(0,0,L1))) *
                    frame( rot_y(-input(robotvar+1)) ) *
                    frame( rot_z(input(robotvar+2) ) )
                ) *
                frame( rot_y(input(robotvar+3) ), Constant(Vector(0,0,L2)))
            ) *
            frame( rot_z(input(robotvar+4) ), Constant(Vector(0,0,L3)))
        ) *
        frame( rot_y(-input(robotvar+5) ) ) *
        frame( rot_z(input(robotvar+6) ), Constant(Vector(0,0,L4)) )
    );
    robotvar=7;
    base_name="robot2";
    Expression<Frame>::Ptr robot2 = cached<Frame>(
        base_name +":ee",
        cached<Frame>(base_name+":wrist",
            cached<Frame>(base_name+":elbow",
                cached<Frame>(base_name+":shoulder",
                    frame( rot_z(input(robotvar)),  Constant(Vector(0,0,L1))) *
                    frame( rot_y(-input(robotvar+1)) ) *
                    frame( rot_z(input(robotvar+2) ) )
                ) *
                frame( rot_y(input(robotvar+3) ), Constant(Vector(0,0,L2)))
            ) *
            frame( rot_z(input(robotvar+4) ), Constant(Vector(0,0,L3)))
        ) *
        frame( rot_y(-input(robotvar+5) ) ) *
        frame( rot_z(input(robotvar+6) ), Constant(Vector(0,0,L4)) )
    );


    // base at [0 0 1], tool at [0 0 0.3] wrt mounting plate: 
    // chain is copied when using kinematic_chain(...)
    Expression<Frame>::Ptr kinchain1 = cached<Frame>(
         Constant(Frame(Vector(0,0,1))) * robot1 * Constant(Frame(Vector(0,0,0.3)))
    );
    // base at [0 0 0], tool at [0 0 0.3] wrt mounting plate: 
    Expression<Frame>::Ptr kinchain2 = cached<Frame>(
         Constant(Frame(Vector(0,0,0))) * robot2 * Constant(Frame(Vector(0,0,0.3)))
    );
 
   
    // z-axes of the end effector of the 2 robots are perpendicular to each other: 
    Expression<double>::Ptr perpendicular = cached<double>(dot(unit_z( rotation( kinchain1 )), unit_z(rotation( kinchain2 ) )));

    std::vector<double> jval(14);
    for (int i=0;i<jval.size();++i) jval[i] = 0.05;
    cout << "numerical " << numerical_derivative<double>(perpendicular, 3,jval[3]) << "\n";


    for (int i=0;i<jval.size();++i) jval[i] = 0.05;
    perpendicular->setInputValues(jval);
    perpendicular->derivative(3);

    perpendicular->setInputValue(3,0.05);    	
    cout << perpendicular->value() << "\t";
    cout << perpendicular->derivative(3) << "\n";
    perpendicular->setInputValue(3,1.0);    	
    cout << perpendicular->value() << "\t";
    cout << perpendicular->derivative(3) << "\n";
   
    ExpressionOptimizer optimizer;
    std::vector<int>    variablelist;
    variablelist.push_back(3);
    optimizer.prepare( variablelist ); 
    perpendicular->addToOptimizer(optimizer);
    std::vector<double> values;
    values.resize(1);
    values[0] = 0.05;
    optimizer.setInputValues( values );
    cout << perpendicular->value() << "\t";
    cout  << perpendicular->derivative(3) << "\n";
    values[0] = 1.0;
    optimizer.setInputValues( values );
    cout << perpendicular->value() << "\t";
    cout  << perpendicular->derivative(3) << "\n";

    std::ofstream of("expressiontree_optimizer.dot");
    perpendicular->write_dotfile(of);

    return 0;
}
