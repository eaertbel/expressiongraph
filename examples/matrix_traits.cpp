#include <kdl/expressiontree.hpp>
#include <kdl/expressiontree_matrix.hpp>
#include <Eigen/Core>
#include <iostream>
int main() {
    
   using namespace KDL;
   using namespace std;
   using namespace Eigen;


  typedef Matrix<double,4,4> Mat;
  typedef Matrix<double,3,3> MatA;
  typedef Matrix<double,3,2> MatB;

  cout << demangle(typeid( AutoDiffTrait<Mat>::ValueType ).name()) << endl;
  cout << demangle(typeid( AutoDiffTrait<Mat>::DerivType ).name()) << endl;
  cout << AutoDiffTrait<Mat>::zeroDerivative() << endl;
  cout << AutoDiffTrait<Mat>::size << endl;
  cout << "A="<< endl;
  MatA A;
  A << 1,2,3,
       4,5,2,
       -1,3,2;
  cout << A << endl;
  MatB B;
  cout << "B="<< endl;
  B << 3,4,
       2,1,
       5,2;
  cout << B << endl;
  Expression< MatA>::Ptr a   = Constant< MatA >( A );
  Expression< MatB>::Ptr b   = Constant< MatB >( B );
  cout << "A*B+B expression " << endl;
  Expression< MatB>::Ptr res = cached< MatB>( addition<3,2>( multiply<3,3,2>(a,b), b) );
  res->print(cout);
  cout << endl;
  cout << "value = " << endl;
  cout << res->value();
  cout << endl;
  cout << "element out of result " << endl;
  Expression<double>::Ptr e = get_element<3,2>(1,0, res);
  e->print(cout);
  cout << "value: " <<  endl;
  cout << e->value() << endl;
};
