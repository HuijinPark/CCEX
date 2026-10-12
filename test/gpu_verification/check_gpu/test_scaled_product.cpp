#include "calculation/scaled_product.h"
#include <cstring>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <stdexcept>
void close(std::complex<double> actual,std::complex<double> expected){
    if(!std::isfinite(std::abs(actual)) || std::abs(actual-expected)>1e-12)
        throw std::runtime_error("scaled product mismatch");
}
int main(){
    ScaledComplex product;
    const ScaledComplex small(std::complex<double>(1e-250,2e-250));
    for(int i=0;i<40;++i)product=product*small;
    for(int i=0;i<40;++i)product=product*small.inverse();
    close(product.value(),{1,0});
    const ScaledComplex complex(std::complex<double>(.37,-.89));
    close((complex.power(2000)*complex.power(-2000)).value(),{1,0});
    close(ScaledComplex(std::complex<double>(0,0)).power(0).value(),{1,0});
    close((complex*ScaledComplex(std::complex<double>(0,0))).value(),{0,0});
    ScaledComplex unequal;unequal.real=.5L;unequal.imag=std::ldexp(.5L,-1200);unequal.exponent=1201;
    if(!std::isinf(unequal.value().real()) || unequal.value().imag()!=1)throw std::runtime_error("component exponents not preserved");
    std::cout<<"PASS scaled products retain finite results across underflow and overflow intermediates\n";

    // improvements.md A9 (after review 20261007g): assembleFactor, the step calculate() uses, against the old
    // always-normalized step, compared after EVERY step on normalized copies (components and exponent, signed zeros and
    // NaN classes included), over short deterministic sequences of edge values and a long random one.
    auto sameLD=[](long double a,long double b){
        if(std::isnan(a)||std::isnan(b))return std::isnan(a)&&std::isnan(b)&&std::signbit(a)==std::signbit(b);
        return a==b && std::signbit(a)==std::signbit(b);
    };
    auto same=[&](ScaledComplex a,ScaledComplex b){
        a.normalize();b.normalize();
        const bool finite=std::isfinite(a.real)&&std::isfinite(a.imag)&&a.real!=0;
        return sameLD(a.real,b.real)&&sameLD(a.imag,b.imag)&&(!finite || a.exponent==b.exponent);
    };
    int steps=0,shortcut=0,gcceShortcut=0;
    auto run=[&](const char* name,const std::vector<std::pair<std::complex<double>,int>>& seq,std::complex<double> zero){
        const ScaledComplex inverseZero=ScaledComplex(zero).inverse();
        const bool unitZero=zero.real()==1.0 && zero.imag()==0.0 && !std::signbit(zero.imag());
        ScaledComplex oldWith,oldWithout,newWith,newWithout;
        for(size_t k=0;k<seq.size();++k){
            const auto [z,p]=seq[k];
            const auto factor=(ScaledComplex(z)*ScaledComplex(zero).inverse()).power(p);   // the old step
            oldWith=oldWith*factor; if(p>0)oldWithout=oldWithout*factor;
            assembleFactor(newWith,newWithout,z,inverseZero,p,unitZero,nullptr);
            if(p==1&&unitZero&&z.real()!=0&&z.imag()!=0&&std::isfinite(z.real())&&std::isfinite(z.imag()))++shortcut;
            if(p==1&&!unitZero&&z.real()!=0&&z.imag()!=0&&std::isfinite(z.real())&&std::isfinite(z.imag())){   // A24's path
                const long double r=z.real()*inverseZero.real-z.imag()*inverseZero.imag,i=z.real()*inverseZero.imag+z.imag()*inverseZero.real;
                if(r!=0&&i!=0&&std::isfinite(r)&&std::isfinite(i))++gcceShortcut;
            }
            ++steps;
            // the final reconstruction: rank 0 multiplies by L0, the output is value()
            const auto ow=(oldWith*ScaledComplex(zero)).value(), nw=([&]{auto x=newWith;x.normalize();return x*ScaledComplex(zero);}()).value();
            if(!same(oldWith,newWith)||!same(oldWithout,newWithout)||std::memcmp(&ow,&nw,sizeof ow)!=0){
                std::cerr<<name<<" step "<<k<<": raw "<<z<<" power "<<p<<" old "<<ow<<" new "<<nw<<"\n";
                throw std::runtime_error("assembly step differs");
            }
        }
    };
    const double dmin=std::numeric_limits<double>::denorm_min(), big=1e308;
    const std::complex<double> one(1,0);
    // signed zeros, negative reals with zero imaginary parts, zeros with nonzero partners, every power kind
    for(std::complex<double> z: {std::complex<double>(0.0,0.0),{-0.0,0.0},{0.0,-0.0},{-0.0,-0.0},{-1,-0.0},{-1,0.0},{1,-0.0},
                                 {-0.0,1},{0.0,-1},{-0.4152891829,-0.0}})
        for(int p: {1,2,-1,-3,0}){
            run("signed zero",{{std::polar(0.9,0.3),1},{z,p},{std::polar(0.8,-0.2),1},{z,1}},one);
            run("signed zero first",{{z,p},{-1.0,1},{z,1}},one);
        }
    // subnormal and huge doubles, then compensating factors (both guard boundaries crossed and recovered)
    run("denorm",{{{dmin,dmin},1},{{dmin,0.5},1},{{1/dmin*0.5,0.25},1},{{0.5,dmin},-2},{{1e-300,1e-300},1}},one);
    run("huge",{{{big,big},1},{{big,-big},1},{{1e-308,1e-308},1},{{big,big},-1},{{2e-308,1e-308},1}},one);
    {std::vector<std::pair<std::complex<double>,int>> seq;
     for(int i=0;i<1002;++i)seq.push_back({{0.5,0.25},1});
     for(int i=0;i<1002;++i)seq.push_back({{2,1},1});
     for(int i=0;i<3;++i)seq.push_back({{0.5,0.25},-400});
     run("boundaries",seq,one);}
    run("non-finite",{{{INFINITY,1},1},{{0.5,0.5},1}},one);
    run("nan",{{{NAN,1},1},{{0.5,0.5},1},{{1,NAN},2}},one);
    // a zeroth cluster that is not 1 (gCCE): the shortcut is off, the cached inverse is used
    run("nonunit L0",{{std::polar(0.9,0.1),1},{std::polar(0.7,-0.3),-2},{std::polar(0.99,0.0),1}},std::complex<double>(0.95,0.05));
    run("L0 1-0i",{{std::polar(0.9,0.1),1},{{-1,-0.0},1}},std::complex<double>(1,-0.0));
    // improvements.md A24: L0 != 1 (gCCE). Every sequence above again under zeroth clusters that are not 1: complex, real
    // of either sign with a signed-zero imaginary part, purely imaginary, tiny and huge (the cached inverse then carries a
    // large exponent), subnormal, and invalid ones (0, inf, nan) where the old expression must be kept.
    {std::vector<std::complex<double>> zeros={{0.95,0.05},{0.3,-0.9},{0.8,0.0},{0.8,-0.0},{-0.7,0.0},{-0.7,-0.0},{0.0,0.6},{-0.0,-0.6},
                                              {1e-300,2e-300},{3e300,-1e300},{dmin,0.5},{0.5,dmin},{dmin,dmin},{1e308,1e308},
                                              {0.0,0.0},{-0.0,0.0},{INFINITY,0.5},{0.5,NAN}};
     for(auto zero: zeros){
         for(std::complex<double> z: {std::complex<double>(0.0,0.0),{-0.0,0.0},{0.0,-0.0},{-0.0,-0.0},{-1,-0.0},{-1,0.0},{1,-0.0},
                                      {-0.0,1},{0.0,-1},{-0.4152891829,-0.0}})
             for(int p: {1,2,-1,-3,0}){
                 run("L0 != 1, signed zero",{{std::polar(0.9,0.3),1},{z,p},{std::polar(0.8,-0.2),1},{z,1}},zero);
                 run("L0 != 1, signed zero first",{{z,p},{-1.0,1},{z,1}},zero);
             }
         // factors whose product with inverse(L0) cancels to zero in one component or is purely real / imaginary
         run("L0 != 1, cancelling",{{zero,1},{-zero,1},{zero*std::complex<double>(0,1),1},{std::conj(zero),1},{{0.5,0.5},1}},zero);
         run("L0 != 1, denorm",{{{dmin,dmin},1},{{dmin,0.5},1},{{1/dmin*0.5,0.25},1},{{0.5,dmin},-2},{{1e-300,1e-300},1}},zero);
         run("L0 != 1, huge",{{{big,big},1},{{big,-big},1},{{1e-308,1e-308},1},{{big,big},-1},{{2e-308,1e-308},1}},zero);
         run("L0 != 1, non-finite",{{{INFINITY,1},1},{{0.5,0.5},1},{{NAN,1},1},{{0.5,0.5},1},{{1,NAN},2}},zero);
     }
     // L0 = (0.5, 0.5): inverse(L0) is exact, so these raw values cancel to an exact zero in one component of the factor
     for(auto zero: {std::complex<double>(0.5,0.5),{0.5,-0.5},{-0.5,0.5}})
         run("L0 != 1, exact cancellation",{{{1,-1},1},{{1,1},1},{{-1,1},1},{{-1,-1},1},{{0.25,0.25},1},{{1,-1},-1}},zero);
     std::vector<std::pair<std::complex<double>,int>> seq;   // both guard boundaries crossed and recovered, L0 != 1
     for(int i=0;i<1002;++i)seq.push_back({{0.5,0.25},1});
     for(int i=0;i<1002;++i)seq.push_back({{2,1},1});
     for(int i=0;i<3;++i)seq.push_back({{0.5,0.25},-400});
     for(auto zero: {std::complex<double>(0.95,0.05),{1e-300,2e-300},{3e300,-1e300}})run("L0 != 1, boundaries",seq,zero);}
    // a long random sequence near 1 with occasional negative powers
    std::mt19937 g(20261007);
    std::uniform_real_distribution<double> mag(0.97,1.0),phase(-0.05,0.05),pick(0,1);
    {std::vector<std::pair<std::complex<double>,int>> seq;
     for(int c=0;c<20000;++c)seq.push_back({std::polar(mag(g),phase(g)),pick(g)<0.01 ? -int(pick(g)*3000) : (pick(g)<0.02 ? 2 : 1)});
     run("random",seq,one);
     // the same under L0 != 1 (gCCE: L0 a decaying, rotating qubit-only coherence), raw values near L0
     for(auto zero: {std::polar(0.93,0.4),std::polar(1e-200,-2.0),std::polar(1.0,3.0)}){
         std::vector<std::pair<std::complex<double>,int>> g2;
         for(auto [z,p]: seq)g2.push_back({z*zero,p});
         run("random, L0 != 1",g2,zero);
     }}
    std::cout<<"PASS the assembly step equals the old normalized step after every step ("<<steps<<" steps, "<<shortcut
             <<" through the cCCE shortcut, "<<gcceShortcut<<" through the gCCE one; signed zeros, subnormal and huge doubles, guard crossings, non-finite values, L0 != 1)\n";
    if(ScaledComplex::lazySafe && gcceShortcut<60000)throw std::runtime_error("the gCCE shortcut was hardly exercised");
    if(!ScaledComplex::lazySafe)std::cout<<"note: long double exponent range too narrow here, deferred normalization is off\n";
}
