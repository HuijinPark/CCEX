#include "calculation/hamiltonian_check.h"
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
using Eigen::Index;
using Eigen::MatrixXcd;
typedef std::complex<double> cd;

// The expression finiteAndHermitian replaced, kept here as the reference verdict.
static bool reference(const MatrixXcd& h){
    return !(!h.allFinite() || (h.adjoint()-h).cwiseAbs().maxCoeff()>1e-12*h.cwiseAbs().maxCoeff());
}
static int cases=0;
static void same(const std::string& name,const MatrixXcd& h,int expected=-1){
    const bool old=reference(h),now=finiteAndHermitian(h);
    ++cases;
    if(old!=now || (expected>=0 && now!=bool(expected)))
        throw std::runtime_error(name+": reference "+std::to_string(old)+", one pass "+std::to_string(now)+
                                 (expected>=0 ? ", expected "+std::to_string(expected) : ""));
}
static MatrixXcd hermitian(Index n,bool sparse){
    MatrixXcd h=MatrixXcd::Random(n,n);
    h=(h+h.adjoint()).eval();
    // what a spin Hamiltonian looks like: the diagonal and the single-flip and flip-flop elements
    if(sparse)for(Index j=0;j<n;++j)for(Index i=0;i<n;++i){const Index x=i^j;if(x && (x&(x-1)) && ((x&(x-1))&((x&(x-1))-1)))h(i,j)=0;}
    return h;
}
template<class F> static double seconds(F f,const MatrixXcd& h,bool& verdict){
    const auto start=std::chrono::steady_clock::now();
    verdict=f(h);
    return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
int main(int argc,char** argv){
    const double inf=std::numeric_limits<double>::infinity(),nan=std::numeric_limits<double>::quiet_NaN();
    std::srand(20260930);
    same("1x1 real",MatrixXcd::Constant(1,1,cd(3,0)),1);
    same("1x1 zero",MatrixXcd::Zero(1,1),1);
    same("1x1 imaginary",MatrixXcd::Constant(1,1,cd(3,1)),0);
    same("1x1 nan",MatrixXcd::Constant(1,1,cd(nan,0)),0);
    // the tile is 64: sizes below, at, across and well past it, and not a multiple of it
    for(Index n:{2,3,63,64,65,127,128,129,200}) for(int sparse=0;sparse<2;++sparse){
        const std::string tag=" n="+std::to_string(n)+(sparse ? " sparse" : " dense");
        const MatrixXcd h=hermitian(n,sparse);
        same("hermitian"+tag,h,1);
        same("all zero"+tag,MatrixXcd::Zero(n,n),1);
        const double top=h.cwiseAbs().maxCoeff();
        // where the perturbation sits: corners, the tile seams, the diagonal, both triangles
        const Index spots[][2]={{0,0},{n-1,n-1},{n/2,n/2},{0,n-1},{n-1,0},{0,1},{1,0},{n/2,n/3},{n/3,n/2},
                                {std::min<Index>(n-1,63),std::min<Index>(n-1,64)},{std::min<Index>(n-1,64),std::min<Index>(n-1,63)}};
        for(const auto& s:spots){
            const std::string at=tag+" ("+std::to_string(s[0])+","+std::to_string(s[1])+")";
            // a bump of size e at one entry moves max|h-h^dagger| to e (an imaginary one on the diagonal to 2e):
            // well below, just below, just above and well above 1e-12*max|h|
            for(double factor:{0.01,0.4,0.49,0.51,0.9,0.999,1.001,1.1,2.0,100.0}){
                MatrixXcd g=h;g(s[0],s[1])+=cd(0,factor*1e-12*top);
                same("imaginary bump "+std::to_string(factor)+at,g);
                g=h;g(s[0],s[1])+=cd(factor*1e-12*top,0);
                same("real bump "+std::to_string(factor)+at,g);
            }
            for(const cd bad:{cd(nan,0),cd(0,nan),cd(inf,0),cd(0,-inf),cd(nan,inf)}){
                MatrixXcd g=h;g(s[0],s[1])=bad;
                same("non-finite"+at,g,0);
                g=MatrixXcd::Zero(n,n);g(s[0],s[1])=bad;
                same("non-finite in zeros"+at,g,0);
            }
            MatrixXcd g=MatrixXcd::Zero(n,n);g(s[0],s[1])=cd(1e-300,0);
            same("one entry in zeros"+at,g,s[0]==s[1]);
        }
        // the threshold itself: max|h|=1 and one entry below the diagonal that is the difference, stored
        // as it is (added to 0.5 the three would round to one number): one step under 1e-12, 1e-12, one over
        if(n>=2){
            const double skews[]={std::nextafter(1e-12,0.0),1e-12,std::nextafter(1e-12,1.0)};
            for(int k=0;k<3;++k){
                MatrixXcd g=MatrixXcd::Zero(n,n);g(0,0)=1;g(1,0)=cd(skews[k],0);
                same("threshold"+tag+" step "+std::to_string(k-1),g,k<2);
                g(1,0)=0;g(0,1)=cd(0,skews[k]);
                same("threshold above the diagonal"+tag+" step "+std::to_string(k-1),g,k<2);
            }
        }
    }
    // max|h| sits below the diagonal only (its partner is two steps smaller), and the difference of
    // another pair walks across 1e-12*max|h| one step at a time: the scale has to come from both triangles
    { double skew=1e-12;
      for(int step=0;step<5;++step,skew=std::nextafter(skew,1.0)){
          MatrixXcd g=MatrixXcd::Zero(3,3);
          g(0,1)=1;g(1,0)=std::nextafter(std::nextafter(1.0,2.0),2.0);
          g(0,2)=cd(0.1,0);g(2,0)=cd(0.1,skew);
          same("scale below the diagonal, step "+std::to_string(step),g);
      } }
    // entries so large that |h| overflows while every entry is finite
    { MatrixXcd g=MatrixXcd::Zero(3,3);g(0,1)=cd(1.5e308,1.5e308);g(1,0)=cd(1.5e308,-1.5e308);same("overflowing modulus",g);
      g(1,0)=cd(1.5e308,1.5e308);same("overflowing modulus, not hermitian",g); }
    std::cout<<"PASS one-pass Hamiltonian check gives the reference verdict in "<<cases<<" cases\n";
    // timing, on request: test_hamiltonian_check.out DIM [DIM...]
    for(int a=1;a<argc;++a) for(int sparse=0;sparse<2;++sparse){
        const Index n=std::atoi(argv[a]);
        const MatrixXcd h=hermitian(n,sparse);
        bool old,now;
        const double told=seconds(reference,h,old),tnow=seconds(finiteAndHermitian,h,now);
        if(old!=now)throw std::runtime_error("timing matrix: verdicts differ");
        std::cout<<"dim "<<n<<(sparse ? " sparse" : " dense ")<<"  reference "<<told<<" s  one pass "<<tnow<<" s\n";
    }
}
