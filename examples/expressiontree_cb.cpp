// example of a MIMO expressiongraph node:
// a multiplication by a 2x2 matrix.
#include <kdl/expressiontree.hpp>
#include <fstream>

using namespace KDL;


class MMult : public MIMO {
    public:
        double a,b,c,d;
        double y[2];
        double dydx[2][2];

        typedef boost::shared_ptr<MMult> Ptr;

        MMult(double _a, double _b, double _c, double _d, Expression<double>::Ptr x1, Expression<double>::Ptr x2):
            MIMO("MMult"),
            a(_a),b(_b),c(_c),d(_d) {
                inputDouble.push_back(x1);
                inputDouble.push_back(x2);
        }

        void compute() {
            if (cached) return;  
            y[0] = a*inputDouble[0]->value() + b*inputDouble[1]->value();
            y[1] = c*inputDouble[0]->value() + d*inputDouble[1]->value();
            dydx[0][0] = a;
            dydx[0][1] = b;
            dydx[1][0] = c;
            dydx[1][1] = d;
            cached=true;
        }

        virtual MIMO::Ptr clone() {
            MIMO::Ptr tmp( new MMult(a,b,c,d,inputDouble[0]->clone(),inputDouble[1]->clone()));
            return tmp;
        } 
};

MMult::Ptr create_MMult(double a, double b, double c, double d, Expression<double>::Ptr x1, Expression<double>::Ptr x2) {
    MMult::Ptr tmp( new MMult(a,b,c,d,x1,x2));
    return tmp;
}


class MMultOut : public MIMO_Output<double> {
    public:
        typedef boost::shared_ptr<MMultOut> Ptr;
        int outputnr;

        MMultOut(MIMO::Ptr m, int _outputnr):MIMO_Output<double>("MMult_Out",m), outputnr(_outputnr) {}
        
        double value() {
            MMult::Ptr p = boost::static_pointer_cast<MMult>(mimo);
            p->compute();
            return p->y[outputnr]; 
        }

        double derivative(int i){
            MMult::Ptr p = boost::static_pointer_cast<MMult>(mimo);
            p->compute();
            if ((0<=i)&&(i<=1)) {
                return p->dydx[outputnr][i];
            } else {
                return 0.0;
            }
        }
        MIMO_Output::Ptr clone() {
            MMultOut::Ptr tmp( new MMultOut(getMIMOClone(),outputnr));
            return tmp;
        }
};

MMultOut::Ptr create_MMult_Output( MMult::Ptr m, int i) {
    MMultOut::Ptr tmp( new MMultOut(m, i));
    return tmp;
}



int main(int argc,char* argv[]) {
    using namespace std;
    Expression<double>::Ptr x1 = input(0);
    Expression<double>::Ptr x2 = input(1);
    MMult::Ptr               m = create_MMult(1,2,3,4,x1,x2); 
    Expression<double>::Ptr y1 = create_MMult_Output(m,0);
    Expression<double>::Ptr y2 = create_MMult_Output(m,1);

    // display the graph:
    ofstream of("mimo.dot");
    write_dotfile_start(of);
    y1->write_dotfile_init();
    y2->write_dotfile_init();
    y1->write_dotfile_update(of);
    y2->write_dotfile_update(of);
    write_dotfile_end(of);
    of.close();

    // evaluate the graph, method 1 for two sets of values:
    y1->setInputValue(0,1.0);
    y2->setInputValue(0,1.0);
    y1->setInputValue(1,2.0);
    y2->setInputValue(1,2.0);
    cout << "method 1 with (x1=1,x2=2)" << endl;
    cout << "y1="<< y1->value() << endl;
    cout << "y2="<< y2->value() << endl;
    cout << "y1 derivatives = " << y1->derivative(0) << "     " << y1->derivative(1) << endl;
    cout << "y2 derivatives = " << y2->derivative(0) << "     " << y2->derivative(1) << endl;
    y1->setInputValue(0,0.0);
    y2->setInputValue(0,0.0);
    y1->setInputValue(1,1.0);
    y2->setInputValue(1,1.0);
    cout << "method 1 with (x1=0,x2=1)" << endl;
    cout << "y1="<< y1->value() << endl;
    cout << "y2="<< y2->value() << endl;
    cout << "y1 derivatives = " << y1->derivative(0) << "     " << y1->derivative(1) << endl;
    cout << "y2 derivatives = " << y2->derivative(0) << "     " << y2->derivative(1) << endl;

    // evaluate the graph, method 2, using expression optimizer:
    ExpressionOptimizer optimizer;
    std::vector<int>    variablelist;
    variablelist.push_back(0);
    variablelist.push_back(1);
    optimizer.prepare( variablelist ); 
    y1->addToOptimizer(optimizer);
    y2->addToOptimizer(optimizer);
    std::vector<double> values(2);
    values[0] = 1.0;
    values[1] = 2.0;
    optimizer.setInputValues( values );
    cout << "method 2 with (x1=1,x2=2)" << endl;
    cout << "y1="<< y1->value() << endl;
    cout << "y2="<< y2->value() << endl;
    cout << "y1 derivatives = " << y1->derivative(0) << "     " << y1->derivative(1) << endl;
    cout << "y2 derivatives = " << y2->derivative(0) << "     " << y2->derivative(1) << endl;
    values[0] = 0.0;
    values[1] = 1.0;
    optimizer.setInputValues( values );
    cout << "method 2 with (x1=0,x2=1)" << endl;
    cout << "y1="<< y1->value() << endl;
    cout << "y2="<< y2->value() << endl;
    cout << "y1 derivatives = " << y1->derivative(0) << "     " << y1->derivative(1) << endl;
    cout << "y2 derivatives = " << y2->derivative(0) << "     " << y2->derivative(1) << endl;

    // cloning: 
    Expression<double>::Ptr y1b = y1->clone();
    Expression<double>::Ptr y2b = y2->clone();

    // display the graph after cloning:
    ofstream of2("mimo_cloning.dot");
    write_dotfile_start(of2);
    y1->write_dotfile_init();
    y2->write_dotfile_init();
    y1b->write_dotfile_init();
    y2b->write_dotfile_init();
    y1->write_dotfile_update(of2);
    y2->write_dotfile_update(of2);
    y1b->write_dotfile_update(of2);
    y2b->write_dotfile_update(of2);
    write_dotfile_end(of2);
    of2.close();


    // evaluation after cloning
    ExpressionOptimizer optimizer2;
    optimizer2.prepare( variablelist ); 
    y1b->addToOptimizer(optimizer2);
    y2b->addToOptimizer(optimizer2);
    values[0] = 1.0;
    values[1] = 2.0;
    optimizer2.setInputValues( values );
    cout << "method 2(cloned) with (x1=1,x2=2)" << endl;
    cout << "y1b="<< y1b->value() << endl;
    cout << "y2b="<< y2b->value() << endl;
    cout << "y1 derivatives = " << y1b->derivative(0) << "     " << y1b->derivative(1) << endl;
    cout << "y2 derivatives = " << y2b->derivative(0) << "     " << y2b->derivative(1) << endl;
    values[0] = 0.0;
    values[1] = 1.0;
    optimizer2.setInputValues( values );
    cout << "method 2(cloned) with (x1=0,x2=1)" << endl;
    cout << "y1b="<< y1b->value() << endl;
    cout << "y2b="<< y2b->value() << endl;
    cout << "y1 derivatives = " << y1b->derivative(0) << "     " << y1b->derivative(1) << endl;
    cout << "y2 derivatives = " << y2b->derivative(0) << "     " << y2b->derivative(1) << endl;
 
    return 0;
}
