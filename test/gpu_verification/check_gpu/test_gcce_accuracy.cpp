// A six-dimensional NV/13C Hamiltonian with a large qubit splitting and small
// bath energies. Compare both refined double paths with a long-double spectral
// reference, including the near-zero coherences amplified by cluster powers.
#include "calculation/gpu/gpu_backend.h"
#include "calculation/gpu/eigen_refinement.h"
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
using Z=std::complex<double>;
using LZ=std::complex<long double>;
using M=Eigen::MatrixXcd;
using LM=Eigen::Matrix<LZ,Eigen::Dynamic,Eigen::Dynamic>;
int main(int argc,char**argv){
    if(std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits){
        std::cerr<<"SKIP: reference needs extended long double precision\n";return 77;
    }
    std::ifstream input(argc>1 ? argv[1] : "test/gpu_verification/check_gpu/data/gcce_scale_separation.txt");
    int d=0,q=0,np=0;input>>d>>q>>np;
    if(d<1 || q<1 || d%q || np!=1)throw std::runtime_error("invalid reference fixture");
    M h(d,d),rho(d,d),pulse(q,q);
    for(M* m:{&h,&rho,&pulse})for(int k=0;k<m->size();++k){
        double re,im;input>>re>>im;m->data()[k]=Z(re,im);
    }
    if(!input)throw std::runtime_error("incomplete reference fixture");
    int b=d/q;M expanded=M::Zero(d,d);
    for(int i=0;i<q;++i)for(int j=0;j<q;++j)for(int k=0;k<b;++k)expanded(i*b+k,j*b+k)=pulse(i,j);
    Eigen::SelfAdjointEigenSolver<M> cpu(h);
    M vectors=cpu.eigenvectors();Eigen::VectorXd values=cpu.eigenvalues(),correction;
    refineEigensystem(h,vectors,values,correction);
    const LM lh=h.cast<LZ>(),lrho=rho.cast<LZ>(),lpulse=expanded.cast<LZ>();
    Eigen::SelfAdjointEigenSolver<LM> reference(lh);
    char description[512];if(!gpuInitialize(0,1,description,sizeof(description))){std::cerr<<description<<'\n';return 77;}
    double largestCpu=0,largestGpu=0,largestLegacy=0;
    try {
        std::unique_ptr<GpuGcce,void(*)(GpuGcce*)> work(gpuGcceCreate(d,q,false,1,h.data(),rho.data(),np,pulse.data()),gpuGcceDestroy);
        double time=0,dt=(float).08,delays[]={.5,.5};
        for(int step=0;step<30;++step,time+=dt){
            auto phases=reference.eigenvalues().cast<LZ>().eval();
            for(int k=0;k<d;++k)phases[k]=std::exp(LZ(0,-reference.eigenvalues()[k]*(static_cast<long double>(time)*.5L)));
            LM free=reference.eigenvectors()*phases.asDiagonal()*reference.eigenvectors().adjoint();
            LM u=free*lpulse*free,expected=u*lrho*u.adjoint();
            Eigen::VectorXcd phase;refinedPhases(values,correction,time*.5,phase);
            M cf=vectors*phase.asDiagonal()*vectors.adjoint();
            M cu=cf*expanded*cf,actual=cu*rho*cu.adjoint();
            for(int k=0;k<d;++k)phase[k]=std::exp(Z(0,-cpu.eigenvalues()[k]*(time*.5)));
            M oldFree=cpu.eigenvectors()*phase.asDiagonal()*cpu.eigenvectors().adjoint();
            M oldU=oldFree*expanded*oldFree,old=oldU*rho*oldU.adjoint();
            M gpu(q,q);gpuGcceEvolve(work.get(),1,&time,delays,gpu.data(),nullptr);
            for(int i=0;i<q;++i)for(int j=0;j<q;++j){
                LZ want=0,got=0,legacy=0;
                for(int k=0;k<b;++k){want+=expected(i*b+k,j*b+k);got+=LZ(actual(i*b+k,j*b+k));legacy+=LZ(old(i*b+k,j*b+k));}
                largestCpu=std::max(largestCpu,double(std::abs(got-want)));
                largestGpu=std::max(largestGpu,double(std::abs(LZ(gpu(i,j))-want)));
                largestLegacy=std::max(largestLegacy,double(std::abs(legacy-want)));
            }
        }
        std::cout<<std::setprecision(12)<<"max_abs_error cpu="<<largestCpu<<" cuda="<<largestGpu<<" legacy_cpu="<<largestLegacy<<'\n';
        if(!std::isfinite(largestCpu) || !std::isfinite(largestGpu) || largestCpu>2e-11 || largestGpu>2e-11)
            throw std::runtime_error("scale-separated Hamiltonian accuracy regression");
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';gpuFinalize();return 1;}
    gpuFinalize();std::cout<<"PASS scale-separated gCCE accuracy\n";
}
