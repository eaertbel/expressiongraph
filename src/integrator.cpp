#include <expressiongraph/expressiontree_expressions.hpp>
#include <expressiongraph/integrator.hpp>

using namespace KDL;

IntegratorAdaptiveSimpson::IntegratorAdaptiveSimpson( 
    double _epsilon, 
    int _minRecDepth, 
    int _maxRecDepth, 
    int _vecsize):
                epsilon(_epsilon), 
                diffRecDepth(_maxRecDepth-_minRecDepth), 
                maxRecDepth(_maxRecDepth), 
                storage( new double[_maxRecDepth*7*_vecsize]),
                currentp(0),
                vecsize(_vecsize),
                maxalloc(0), alloclimit(_maxRecDepth*7*_vecsize+1) 
    {
    }


IntegratorAdaptiveSimpson::~IntegratorAdaptiveSimpson() {
    delete[] storage;
}

int IntegratorAdaptiveSimpson::get_err(){
    return err;
}

//
// Vectorizing:
//     std::function< vec(double)> f
//     Integral variable type (scalar) : a,b,eps    + m,h,lm,rm
//     Integrand/Integral types (vector) : whole, fa,fb,fm  + flm,frm, left, right, delta
//     fabs(vec) ==> maximum of fabs(elements)
//Adaptive Simpson's Rule, Recursive Core 
//
// minRecDepth 
// maxRecDepth 
//
// depth = maxRecDepth-rec >= minRecDepth
// rec <= maxRecDepth-minRecDepth
//
// maxalloc = 7*depth*n
double* 
IntegratorAdaptiveSimpson::adaptiveSimpsonsAux(IntegrandType f, double a, double b, double eps,
                                  double* whole, double* fa, double* fb, double* fm, int rec) 
{
    double* result = alloc();
    double m   = (a + b)/2;
    double h   = (b - a)/2;
    double lm  = (a + m)/2;  
    double rm  = (m + b)/2;
    if ((eps/2 == eps) || (a == lm)) { err = err_EDOM; return whole; }
    double* flm=alloc(); f(lm,flm);
    double* frm=alloc(); f(rm,frm);
    double* left=alloc(); 
    double* right=alloc(); 
    double* delta=alloc();
    double h6 = h/6;
    double max_fabs=0;
    for (int i=0;i<vecsize;++i) {
        left[i]  = h6*(fa[i]+4*flm[i]+fm[i]);
        right[i] = h6*(fm[i]+4*frm[i]+fb[i]);
        delta[i] = left[i] + right[i] - whole[i];
        double fbs=fabs(delta[i]);
        if (fbs > max_fabs) max_fabs=fbs;
    } 
    if (rec<=diffRecDepth) {
        if ((rec <= 0 || max_fabs <= 15*eps)) {
            for (int i=0;i<vecsize;++i) {
                result[i] = left[i] + right[i] + (delta[i])/15; 
            }
            pop(5); // flm,frm,left,right,delta
            return result;
         }
         if (rec <= 0 && err != err_EDOM) err = err_ERANGE;  // depth limit too shallow
    }
    pop(1); // delta
    double* result1=adaptiveSimpsonsAux(f, a, m, eps/2, left,  fa, fm, flm, rec-1);
    double* result2=adaptiveSimpsonsAux(f, m, b, eps/2, right, fm, fb, frm, rec-1);
    for (int i=0;i<vecsize;++i) {
        result[i]=result1[i]+result2[i];
    }
    pop(6);  // result1, result2, flm, frm, left, right
    return result;
}

int IntegratorAdaptiveSimpson::integrate( IntegrandType f, double a,double b, double* result) {
    currentp = 0;
    err = 0;
    double h  = b - a;
    if (fabs(h)<= 1E-15 ) return 0;
    double* fa=alloc();   
    double* fb=alloc();   
    double* fm=alloc();
    double* S =alloc();
    f(a, fa ); fa_value = fa[0];
    f(b, fb ); fb_value = fb[0];
    f((a + b)/2, fm );
    double h6=h/6;
    for (int i=0;i<vecsize;++i) {
        S[i]  = h6*(fa[i] + 4*fm[i] + fb[i]);
    }
    double* retval = adaptiveSimpsonsAux(f, a, b, epsilon, S, fa, fb, fm, maxRecDepth);
    for (int i=0;i<vecsize;++i) {
        result[i] = retval[i];
    }
    pop(5);
    return err;
}

// we need FunctionEvaluation exposed for setInputValues
//IntegrandType make_function( FunctionDefinition::Ptr func, const std::vector<int>& ndx ) {
IntegrandType make_function( FunctionDefinition::Ptr func ) {
    auto var = Variable<double>({});
    auto funceval = new FunctionEvaluation<double>(func,{var});
    return [funceval,var](double arg,double* result) {
        var->setValue(arg); 
        result[0] = funceval->value();
        result[1] = funceval->value();
        //for (int i=0;i<ndx.size();++i) {
        //    result[i+1]=funceval->derivative(i);
        //}
        return;
    };
}

// we need FunctionEvaluation exposed for setInputValues
//IntegrandType make_function( FunctionDefinition::Ptr func, const std::vector<int>& ndx ) {
IntegrandType make_function( FunctionDefinition::Ptr func, int& counter ) {
    auto var = Variable<double>({});
    auto funceval = new FunctionEvaluation<double>(func,{var});
    return [funceval,var,&counter](double arg,double* result) {
        counter++;
        var->setValue(arg); 
        result[0] = funceval->value();
        result[1] = funceval->value();
        //for (int i=0;i<ndx.size();++i) {
        //    result[i+1]=funceval->derivative(i);
        //}
        return;
    };
}



