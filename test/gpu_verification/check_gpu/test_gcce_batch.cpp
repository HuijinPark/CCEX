// Independent general-exponential checks for heterogeneous Hamiltonians in a
// homogeneous-dimension GPU cluster batch. Full output is checked through the
// sequential API between batches to exercise buffer swaps and reuse as well.
#include "calculation/gpu/gpu_backend.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
using Z=std::complex<double>;using M=Eigen::MatrixXcd;using V=Eigen::VectorXcd;
static void check(int q,int b,bool vectorMode,int pulses,bool refine=true){
    const int d=q*b,n=5,batch=3;
    std::vector<M> h,rho;
    std::vector<std::vector<M>> expanded;
    std::vector<std::unique_ptr<GpuGcce,void(*)(GpuGcce*)>> owned;
    std::vector<GpuGcce*> works;
    for(int c=0;c<n;++c){
        M a=M::Random(d,d);h.push_back((a+a.adjoint()).eval()/double(d));
        V psi=V::Random(d);psi.normalize();
        if(vectorMode)rho.push_back(psi*psi.adjoint());
        else {M mix=M::Random(d,d);M r=mix*mix.adjoint();r/=r.trace().real();rho.push_back(r);}
        std::vector<Z> packed((size_t)pulses*q*q);expanded.emplace_back();
        for(int p=0;p<pulses;++p){
            M a=M::Random(q,q),pulse=(Z(0,-.4)*(a+a.adjoint())).exp();
            std::copy(pulse.data(),pulse.data()+q*q,packed.data()+(size_t)p*q*q);
            M pe=M::Zero(d,d);
            for(int i=0;i<q;++i)for(int j=0;j<q;++j)for(int k=0;k<b;++k)pe(i*b+k,j*b+k)=pulse(i,j);
            expanded.back().push_back(pe);
        }
        owned.emplace_back(gpuGcceCreate(d,q,vectorMode,batch,h.back().data(),
            vectorMode ? psi.data() : rho.back().data(),pulses,packed.data(),refine),gpuGcceDestroy);
        works.push_back(owned.back().get());
    }
    std::vector<double> delays(pulses+1);double sum=0;
    for(int p=0;p<=pulses;++p)sum+=(delays[p]=p+1.0);
    for(double& x:delays)x/=sum;
    double maximum=0;
    for(int count:{3,1,2,3}){
        double times[]={.07,.619,1.327};std::vector<Z> reduced((size_t)n*count*q*q);
        gpuGcceEvolveBatch(works.data(),n,count,times,delays.data(),reduced.data());
        const auto before=gpuMemoryStats();
        std::vector<Z> repeated(reduced.size());
        gpuGcceEvolveBatch(works.data(),n,count,times,delays.data(),repeated.data());
        const auto after=gpuMemoryStats();
        if(before.pointerUploads!=after.pointerUploads || before.allocations!=after.allocations || reduced!=repeated)
            throw std::runtime_error("gCCE unchanged batch did not reuse tables/output exactly");
        for(int c=0;c<n;++c){
            std::vector<Z> single((size_t)count*q*q),full((size_t)count*d*d);
            gpuGcceEvolve(works[c],count,times,delays.data(),single.data(),full.data());
            for(int t=0;t<count;++t){
                M u=M::Identity(d,d);
                for(int p=0;p<=pulses;++p){u=(((Z(0,-times[t]*delays[p])*h[c]).exp())*u).eval();if(p<pulses)u=(expanded[c][p]*u).eval();}
                M expected=u*rho[c]*u.adjoint();M partial=M::Zero(q,q);
                for(int i=0;i<q;++i)for(int j=0;j<q;++j)for(int k=0;k<b;++k)partial(i,j)+=expected(i*b+k,j*b+k);
                M actual=Eigen::Map<const M>(reduced.data()+((size_t)c*count+t)*q*q,q,q);
                M sequential=Eigen::Map<const M>(single.data()+(size_t)t*q*q,q,q);
                M fullActual=Eigen::Map<const M>(full.data()+(size_t)t*d*d,d,d);
                double error=std::max({(actual-partial).cwiseAbs().maxCoeff(),
                    (sequential-partial).cwiseAbs().maxCoeff(),(fullActual-expected).cwiseAbs().maxCoeff()});
                if(!actual.allFinite() || !std::isfinite(error) || error>3e-11)throw std::runtime_error("gCCE cluster batch mismatch");
                maximum=std::max(maximum,error);
            }
        }
        std::swap(works[1],works[2]);std::swap(h[1],h[2]);
        std::swap(rho[1],rho[2]);std::swap(expanded[1],expanded[2]);
    }
    std::cout<<"PASS gCCE batch q="<<q<<" bath="<<b<<" vector="<<vectorMode<<" pulses="<<pulses<<" refine="<<refine<<" max_error="<<maximum<<'\n';
}
int main(){
    char description[512];if(!gpuInitialize(0,1,description,sizeof(description))){std::cerr<<description<<'\n';return 77;}
    try {for(bool v:{false,true}){for(int p:{0,1,3,4})check(3,2,v,p);check(1,1,v,2);check(2,3,v,1);check(3,2,v,3,false);}}
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';gpuFinalize();return 1;}
    gpuFinalize();
}
