// Independent small-system oracle: general matrix exponential, not the GPU's
// eigendecomposition. Covers noncommuting pulses, mixed states and partial traces.
#include "calculation/gpu/gpu_backend.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <algorithm>
#include <complex>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using Z=std::complex<double>;
using M=Eigen::MatrixXcd;
using V=Eigen::VectorXcd;

static void requireClose(const M& actual, const M& expected, const char* what){
    double error=(actual-expected).cwiseAbs().maxCoeff();
    if (!actual.allFinite() || error>3e-11)
        throw std::runtime_error(std::string(what)+" error="+std::to_string(error));
}

static M traceBath(const M& full, int q, int b){
    M out=M::Zero(q,q);
    for (int i=0; i<q; ++i) for (int j=0; j<q; ++j)
        for (int k=0; k<b; ++k) out(i,j)+=full(i*b+k,j*b+k);
    return out;
}

static void runCase(int q, int b, bool vectorMode, int npulse, bool degenerate){
    int d=q*b;
    M a=M::Random(d,d);
    M h=(a+a.adjoint()).eval()/double(d);
    if (degenerate){
        // Complex, rotated eigenspaces with exactly repeated eigenvalues.
        M rotation=(Z(0,1)*h).exp();
        V diag(d);
        for (int k=0;k<d;++k) diag[k]=(k%2 ? 0.75 : -0.5);
        h=rotation*diag.asDiagonal()*rotation.adjoint();
    }
    V psi=V::Random(d); psi.normalize();
    M rho;
    if (vectorMode) rho=psi*psi.adjoint();
    else {
        M mixture=M::Random(d,d);
        rho=mixture*mixture.adjoint(); rho/=rho.trace().real();
    }
    std::vector<M> pulses, expanded;
    std::vector<Z> packed((size_t)npulse*q*q);
    for (int p=0;p<npulse;++p){
        M pa=M::Random(q,q);
        M pulse=(Z(0,-0.7)*(pa+pa.adjoint())).exp();
        pulses.push_back(pulse);
        std::copy(pulse.data(),pulse.data()+q*q,packed.data()+(size_t)p*q*q);
        M full=M::Zero(d,d);
        for (int i=0;i<q;++i) for (int j=0;j<q;++j)
            for (int k=0;k<b;++k) full(i*b+k,j*b+k)=pulse(i,j);
        expanded.push_back(full);
    }
    std::vector<double> fractions(npulse+1);
    double sum=0;
    for (int p=0;p<=npulse;++p) sum+=(fractions[p]=p+1.0);
    for (double& v: fractions) v/=sum;
    std::unique_ptr<GpuGcce,void(*)(GpuGcce*)> work(
        gpuGcceCreate(d,q,vectorMode,3,h.data(),vectorMode ? psi.data() : rho.data(),npulse,packed.data()),
        gpuGcceDestroy);
    // Reuse buffers across a full batch, a tail and another full batch. Request
    // full output only on alternating calls, exercising both output paths.
    for (int iteration=0;iteration<3;++iteration){
        int count=iteration==1 ? 1 : 3;
        double times[]={0.0,0.37,1.31};
        if (iteration==1) times[0]=0.91;
        std::vector<Z> reduced((size_t)count*q*q), full((size_t)count*d*d);
        gpuGcceEvolve(work.get(),count,times,fractions.data(),reduced.data(),iteration==2 ? nullptr : full.data());
        for (int t=0;t<count;++t){
            M u=M::Identity(d,d);
            for (int p=0;p<=npulse;++p){
                M free=(Z(0,-times[t]*fractions[p])*h).exp();
                u=(free*u).eval();
                if (p<npulse) u=(expanded[p]*u).eval();
            }
            M expected=u*rho*u.adjoint();
            M actual=Eigen::Map<const M>(reduced.data()+(size_t)t*q*q,q,q);
            requireClose(actual,traceBath(expected,q,b),"partial trace");
            if (std::abs(actual.trace()-Z(1,0))>3e-11) throw std::runtime_error("trace conservation");
            if (iteration!=2)
                requireClose(Eigen::Map<const M>(full.data()+(size_t)t*d*d,d,d),expected,"full density matrix");
        }
    }
    std::cout << "PASS q="<<q<<" bath="<<b<<" vector="<<vectorMode
              <<" pulses="<<npulse<<" degenerate="<<degenerate<<'\n';
}

int main(){
    char description[512]={};
    if (!gpuInitialize(0,1,description,sizeof(description))){
        std::cerr<<description<<'\n'; return 77;
    }
    std::cout<<description<<'\n';
    try {
        for (bool vectorMode : {true,false}){
            runCase(3,1,vectorMode,0,false);
            runCase(2,3,vectorMode,1,false);
            runCase(3,4,vectorMode,3,true);
            runCase(1,1,vectorMode,1,true);
        }
    }catch (const std::exception& e){
        std::cerr<<"FAIL: "<<e.what()<<'\n'; gpuFinalize(); return 1;
    }
    gpuFinalize();
}
