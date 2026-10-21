// The parallel GPU refinement (A3) on clusters larger than the accuracy fixtures: an NV spin-1 with a
// large splitting and six spin-1/2 nuclei in two groups of equal Zeeman energy (degenerate bath levels,
// dim 192), an odd dimension (15, random bath with a degenerate pair) and the smallest ones (1, 2).
// The GPU Hahn echo is compared with a long-double spectral reference, next to the CPU refined path.
#include "calculation/gpu/gpu_backend.h"
#include "calculation/gpu/eigen_refinement.h"
#include <Eigen/Dense>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
using Z=std::complex<double>;
using LZ=std::complex<long double>;
using M=Eigen::MatrixXcd;
using LM=Eigen::Matrix<LZ,Eigen::Dynamic,Eigen::Dynamic>;

static M kron(const M& a,const M& b){
    M out(a.rows()*b.rows(),a.cols()*b.cols());
    for(int i=0;i<a.rows();++i)for(int j=0;j<a.cols();++j)out.block(i*b.rows(),j*b.cols(),b.rows(),b.cols())=a(i,j)*b;
    return out;
}
static M site(const M& op,int k,int n){
    M out=M::Identity(1,1);
    for(int i=0;i<n;++i)out=kron(out,i==k ? op : M::Identity(2,2));
    return out;
}

// Maximum over 30 Hahn times of |reduced density error| for the GPU and for the CPU refined path.
static void check(const std::string& name,const M& h,int q,const M& rho,const M& pulse){
    const int d=h.rows(),b=d/q;
    M expanded=M::Zero(d,d);
    for(int i=0;i<q;++i)for(int j=0;j<q;++j)for(int k=0;k<b;++k)expanded(i*b+k,j*b+k)=pulse(i,j);
    Eigen::SelfAdjointEigenSolver<M> cpu(h);
    M vectors=cpu.eigenvectors();Eigen::VectorXd values=cpu.eigenvalues(),correction;
    refineEigensystem(h,vectors,values,correction);
    Eigen::SelfAdjointEigenSolver<LM> reference(h.cast<LZ>());
    if(cpu.info()!=Eigen::Success || reference.info()!=Eigen::Success)throw std::runtime_error(name+": eigensolver failed");
    const LM lrho=rho.cast<LZ>(),lpulse=expanded.cast<LZ>();
    std::unique_ptr<GpuGcce,void(*)(GpuGcce*)> work(gpuGcceCreate(d,q,false,1,h.data(),rho.data(),1,pulse.data()),gpuGcceDestroy);
    double largestCpu=0,largestGpu=0,time=0;const double dt=(float).08,delays[]={.5,.5};
    for(int step=0;step<30;++step,time+=dt){
        Eigen::Matrix<LZ,Eigen::Dynamic,1> phases(d);
        for(int k=0;k<d;++k)phases[k]=std::exp(LZ(0,-reference.eigenvalues()[k]*(static_cast<long double>(time)*.5L)));
        LM free=reference.eigenvectors()*phases.asDiagonal()*reference.eigenvectors().adjoint();
        LM u=free*lpulse*free,expected=u*lrho*u.adjoint();
        Eigen::VectorXcd phase;refinedPhases(values,correction,time*.5,phase);
        M cf=vectors*phase.asDiagonal()*vectors.adjoint();
        M cu=cf*expanded*cf,actual=cu*rho*cu.adjoint();
        M gpu(q,q);gpuGcceEvolve(work.get(),1,&time,delays,gpu.data(),nullptr);
        for(int i=0;i<q;++i)for(int j=0;j<q;++j){
            LZ want=0,got=0;
            for(int k=0;k<b;++k){want+=expected(i*b+k,j*b+k);got+=LZ(actual(i*b+k,j*b+k));}
            const double cpuError=double(std::abs(got-want)),gpuError=double(std::abs(LZ(gpu(i,j))-want));
            // std::max keeps the old maximum over a NaN, so a non-finite error fails here (review R1)
            if(!std::isfinite(cpuError) || !std::isfinite(gpuError))throw std::runtime_error(name+": non-finite result");
            largestCpu=std::max(largestCpu,cpuError);largestGpu=std::max(largestGpu,gpuError);
        }
    }
    std::cout<<std::setprecision(6)<<name<<" dim="<<d<<" max_abs_error cpu="<<largestCpu<<" cuda="<<largestGpu<<'\n';
    // fixed before the first run: the GPU within 2e-11 or within twice the CPU refined path
    if(largestGpu>std::max(2e-11,2*largestCpu))
        throw std::runtime_error(name+": parallel GPU refinement accuracy regression");
}

int main(){
    if(std::numeric_limits<long double>::digits<=std::numeric_limits<double>::digits){
        std::cerr<<"SKIP: reference needs extended long double precision\n";return 77;
    }
    char description[512];if(!gpuInitialize(0,1,description,sizeof(description))){std::cerr<<description<<'\n';return 77;}
    try {
        const Z i1(0,1);
        M sx(2,2),sy(2,2),sz(2,2);sx<<0,.5,.5,0;sy<<0,-.5*i1,.5*i1,0;sz<<.5,0,0,-.5;
        {   // NV (ms = -1, 0, +1) and six nuclei: dim 3 x 64
            const int n=6;const int b=1<<n;
            M bath=M::Zero(b,b),hyperfine=M::Zero(b,b);
            for(int k=0;k<n;++k){
                bath+=(k<3 ? 6.7e3 : 1.1e4)*site(sz,k,n);
                hyperfine+=(40.0+35*k)*site(sz,k,n)+(9.0-2*k)*site(sx,k,n);
                for(int l=0;l<k;++l)bath+=(0.8+0.3*k-0.2*l)*(site(sx,k,n)*site(sx,l,n)+site(sy,k,n)*site(sy,l,n)-2.*site(sz,k,n)*site(sz,l,n));
            }
            M qubit=M::Zero(3,3),qz=M::Zero(3,3);qubit(0,0)=1.79e7;qubit(2,2)=2.31e7;qz(0,0)=-1;qz(2,2)=1;
            M h=kron(qubit,M::Identity(b,b))+kron(M::Identity(3,3),bath)+kron(qz,hyperfine);
            M psi=M::Zero(3,1);psi(0,0)=psi(1,0)=std::sqrt(.5);
            M rho=kron(psi*psi.adjoint(),M::Identity(b,b)/double(b));
            M pulse=M::Zero(3,3);pulse(0,1)=pulse(1,0)=-i1;pulse(2,2)=1;
            check("nv_six_nuclei",h,3,rho,pulse);
        }
        {   // odd dimension: qubit 3 x random bath 5 with a degenerate pair
            std::srand(7);M a=M::Random(5,5);M bath=(a+a.adjoint())*20.;
            Eigen::SelfAdjointEigenSolver<M> s(bath);Eigen::VectorXd e=s.eigenvalues();e[1]=e[0];
            bath=s.eigenvectors()*e.cast<Z>().asDiagonal()*s.eigenvectors().adjoint();
            M qubit=M::Zero(3,3);qubit(0,0)=1.79e7;qubit(2,2)=2.31e7;
            M c=M::Random(5,5);M coupling=(c+c.adjoint())*3.;M qz=M::Zero(3,3);qz(0,0)=-1;qz(2,2)=1;
            M h=kron(qubit,M::Identity(5,5))+kron(M::Identity(3,3),bath)+kron(qz,coupling);
            M psi=M::Zero(3,1);psi(0,0)=psi(1,0)=std::sqrt(.5);
            M rho=kron(psi*psi.adjoint(),M::Identity(5,5)/5.);
            M pulse=M::Zero(3,3);pulse(0,1)=pulse(1,0)=-i1;pulse(2,2)=1;
            check("odd_dimension",h,3,rho,pulse);
        }
        {   // the smallest dimensions
            M h1(1,1);h1<<2.5e3;M one=M::Identity(1,1);check("dim1",h1,1,one,one);
            M h2(2,2);h2<<1.2e7,3.,3.,-1.2e7;M rho2(2,2);rho2<<.5,.5,.5,.5;M x(2,2);x<<0,1,1,0;check("dim2",h2,2,rho2,x);
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';gpuFinalize();return 1;}
    gpuFinalize();std::cout<<"PASS parallel GPU refinement on larger, odd and smallest dimensions\n";
}
