#ifndef CCEX_SCALED_PRODUCT_H
#define CCEX_SCALED_PRODUCT_H
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <limits>

// Keep a product's binary exponent separate until the final observable is
// formed. A representable final result can otherwise be lost when a local MPI
// product underflows before later inverse-cluster factors cancel its scale.
struct ScaledComplex {
    long double real=1,imag=0;
    long long exponent=0;
    ScaledComplex()=default;
    explicit ScaledComplex(std::complex<double> z):real(z.real()),imag(z.imag()){normalize();}
    void normalize(){
        const long double scale=std::max(std::abs(real),std::abs(imag));
        if(!std::isfinite(scale) || scale==0){exponent=0;return;}
        int shift=0;std::frexp(scale,&shift);
        real=std::ldexp(real,-shift);imag=std::ldexp(imag,-shift);exponent+=shift;
    }
    ScaledComplex inverse() const {
        ScaledComplex out;
        const long double norm=real*real+imag*imag;
        out.real=real/norm;out.imag=-imag/norm;out.exponent=-exponent;out.normalize();return out;
    }
    ScaledComplex operator*(const ScaledComplex& other) const {
        ScaledComplex out;
        out.real=real*other.real-imag*other.imag;
        out.imag=real*other.imag+imag*other.real;
        out.exponent=exponent+other.exponent;out.normalize();return out;
    }
    // Deferred normalization (improvements.md A9) needs a long double exponent range much wider than double's: a raw
    // double (down to its subnormals, up to ~2^1024) times an accumulator kept in [2^-1000, 2^1000] must stay normal.
    // That holds for the x87 80-bit and IEEE quad formats; elsewhere the always-normalized path is used.
    static constexpr bool lazySafe=std::numeric_limits<long double>::min_exponent<=-16000 &&
                                   std::numeric_limits<long double>::max_exponent>=16000;
    // A value taken as it is, not normalized (the exponent stays 0) - only where lazySafe.
    static ScaledComplex raw(std::complex<double> z){
        if(!lazySafe)return ScaledComplex(z);
        ScaledComplex x;x.real=z.real();x.imag=z.imag();x.exponent=0;return x;
    }
    // In-place product. Where lazySafe it normalizes only when the larger component leaves [2^-1000, 2^1000] or is zero
    // or not finite (then, as normalize does, the exponent is reset). Rescaling by a power of two commutes with rounding
    // while the products stay normal, which this range keeps for the larger components; a minor component more than
    // ~2^-15000 below the larger one can round differently (far below what a double output can show).
    void multiplyLazy(const ScaledComplex& other){
        if(!lazySafe){*this=*this*other;return;}
        const long double r=real*other.real-imag*other.imag, i=real*other.imag+imag*other.real;
        real=r;imag=i;exponent+=other.exponent;
        const long double scale=std::max(std::abs(r),std::abs(i));
        if(!(scale>=0x1p-1000L && scale<=0x1p1000L))normalize();
    }
    ScaledComplex power(int power) const {
        long long n=power;ScaledComplex base=*this,out;
        if(n<0){base=base.inverse();n=-n;}
        while(n){if(n&1)out=out*base;n>>=1;if(n)base=base*base;}
        return out;
    }
    std::complex<double> value() const {
        auto component=[&](long double x){
            if(!std::isfinite(x) || x==0)return double(x);
            int shift=0;const long double fraction=std::frexp(x,&shift);
            const long long total=exponent+shift;
            if(total>std::numeric_limits<double>::max_exponent+2)
                return std::copysign(std::numeric_limits<double>::infinity(),double(x));
            if(total<std::numeric_limits<double>::min_exponent-std::numeric_limits<double>::digits-2)
                return std::copysign(0.0,double(x));
            return double(std::ldexp(fraction,int(total)));
        };
        return {component(real),component(imag)};
    }
};

// One step of the CCE product assembly in calculate() (improvements.md A9), shared with the unit test so that the test
// exercises this code: multiply the factor (raw / L0)^power into the product with and, for power > 0, without the
// negative exponents. unitZero says that L0 is exactly (1, +0) at every time point of the state. The shortcut takes the
// cluster value as the factor only for power 1 and finite components that are both nonzero: there the old
// (raw * inverse(L0)).power(1) is exact and equal to raw, signed zeros included; zeros, non-finite values and every
// other case keep the old expression. Returns the factor's value when saved is requested (savemode cluster_coherence).
// L0 != 1 (gCCE, improvements.md A24), power 1, where lazySafe: the factor is raw * inverse(L0) formed from the raw value
// as it is, not normalized. The old expression normalizes raw first (an exact power-of-two scaling), multiplies, normalizes
// and multiplies by 1 in power(1). With finite raw and inverse(L0) every product here is a normal long double or zero
// (raw is a double, the cached inverse is normalized and its minor part is not below ~2^-2200 of its major part), so
// the unscaled product is the old one times a power of two, component by component, the signs of zeros included; when
// both of its components are finite and nonzero the old power(1) leaves them as they are (1*x - 0*y = x), so the factor
// has the old value and the deferred normalization of multiplyLazy gives the old product. A non-finite inverse(L0) gives
// a non-finite factor and so the old expression too. (An exact zero here is +0, as the old one would be; the test keeps
// the argument short.)
inline void assembleFactor(ScaledComplex& with,ScaledComplex& without,std::complex<double> raw,
                           const ScaledComplex& inverseZero,int power,bool unitZero,std::complex<double>* saved){
    if(power==1 && !saved && raw.real()!=0 && raw.imag()!=0 &&
       std::isfinite(raw.real()) && std::isfinite(raw.imag())){
        if(unitZero){
            const ScaledComplex factor=ScaledComplex::raw(raw);
            with.multiplyLazy(factor);without.multiplyLazy(factor);
            return;
        }
        if(ScaledComplex::lazySafe){
            ScaledComplex factor;
            factor.real=raw.real()*inverseZero.real-raw.imag()*inverseZero.imag;
            factor.imag=raw.real()*inverseZero.imag+raw.imag()*inverseZero.real;
            factor.exponent=inverseZero.exponent;
            if(factor.real!=0 && factor.imag!=0 && std::isfinite(factor.real) && std::isfinite(factor.imag)){
                with.multiplyLazy(factor);without.multiplyLazy(factor);
                return;
            }
        }
    }
    const auto factor=(ScaledComplex(raw)*inverseZero).power(power);
    with.multiplyLazy(factor);
    if(power>0)without.multiplyLazy(factor);
    if(saved)*saved=factor.value();
}

// improvements.md A34: count power-1 cluster values raw[0], raw[stride], ... of one time point multiplied, in that order,
// into one product - the same operations as count calls of assembleFactor(product, other, raw[k*stride], inverseZero, 1,
// unitZero, nullptr) for this product, so the same bits. calculate() gathers a block of clusters first so that the run is
// long, and the product is held in local variables over the run instead of being stored and reloaded at every factor
// (x87 long double stores are slow). In the rare cases (a normalization, a value off assembleFactor's shortcuts) the
// product goes back to `product`, is handled there out of line, and is taken up again; the locals' addresses never
// escape, so the compiler can keep them in registers. The cCCE and gCCE runs are separate functions because x87 register
// allocation is fragile: in one function, a change to one loop slowed the other twofold (icpx).
//
// Guard-free runs: multiplyLazy normalizes when the larger component of the product leaves [2^-1000, 2^1000]. Before a
// run, a scan of the doubles bounds every factor's size by powers of two; when those bounds keep the product (and every
// partial product of the run) at least 2^10 inside that window, the guard cannot fire, and the loop leaves it out - the
// same operations as with it. Sizes: for a component maximum a = f 2^e (frexp, f in [1/2, 1)) a complex number lies in
// [2^(e-1), 2^(e+1/2)); a product's larger component M and modulus |p| satisfy M <= |p| < 2 M; each multiplication
// changes |p| by |factor| times (1 +- 2^-60) at most (x87 rounding, 64-bit significands, all values normal).
namespace scaled_run {
__attribute__((noinline,cold)) inline void normalize(ScaledComplex& product){product.normalize();}
__attribute__((noinline,cold)) inline void other(ScaledComplex& product,std::complex<double> raw,
                                                  const ScaledComplex& inverseZero,bool unitZero){
    ScaledComplex unused;
    assembleFactor(product,unused,raw,inverseZero,1,unitZero,nullptr);
}
// both components finite and nonzero: assembleFactor's condition for its shortcuts
inline bool regular(std::complex<double> z){
    return z.real()!=0 && z.imag()!=0 && std::isfinite(z.real()) && std::isfinite(z.imag());
}
// e with |x| = f 2^e, f in [1/2, 1) (as frexp) for a finite nonzero double, read from its bits; a subnormal (and zero, not
// used) gives -1100, below every subnormal's e, which only makes the bounds below wider
inline int sizeExponent(double x){
    unsigned long long bits;std::memcpy(&bits,&x,sizeof bits);
    const int biased=int((bits>>52)&0x7ff);
    return biased==0 ? -1100 : biased-1022;
}
// the binary exponent e of max(|re|, |im|) = f 2^e, f in [1/2, 1); false if that is zero or not finite
inline bool scaleExponent(long double re,long double im,int& e){
    const long double m=std::max(std::abs(re),std::abs(im));
    if(!(m>0 && m<=std::numeric_limits<long double>::max()))return false;
    std::frexp(m,&e);return true;
}
// the guard cannot fire when the product starts at M0 = f 2^m and the factors' size exponents add up, over any prefix of
// the run, to no less than low and no more than high: then 2^(m-2+low) <= M <= 2^(m+2+high), with 2^10 to spare
inline bool guardFree(int m,long long low,long long high){return m-2+low>=-990 && m+2+high<=990;}
// multiplyLazy by (fr, fi, fe) of the product held in re, im, exponent
#define CCEX_SCALED_RUN_LAZY(FR,FI,FE) do{ \
        const long double r_=re*(FR)-im*(FI), i_=re*(FI)+im*(FR); \
        re=r_;im=i_;exponent+=(FE); \
        const long double scale_=std::max(std::abs(r_),std::abs(i_)); \
        if(__builtin_expect(!(scale_>=0x1p-1000L && scale_<=0x1p1000L),0)){ \
            product.real=re;product.imag=im;product.exponent=exponent;normalize(product); \
            re=product.real;im=product.imag;exponent=product.exponent; \
        } \
    }while(0)
#define CCEX_SCALED_RUN_OTHER(Z,UNIT) do{ \
        product.real=re;product.imag=im;product.exponent=exponent;other(product,Z,inverseZero,UNIT); \
        re=product.real;im=product.imag;exponent=product.exponent; \
    }while(0)
// The scan of a run, once for the two products of a time point: whether every value is regular (the usual case) and
// bounds on the factors' sizes, as sums of size exponents that hold over any prefix of the run.
// L0 = 1 (cCCE): the factor is ScaledComplex::raw(z), exponent 0; |z| in [2^(e-1), 2^(e+1/2)), so the size exponents add up
// to at least the sum of min(e-1, 0) and at most the sum of max(e+1, 0).
// L0 != 1 (gCCE), see runGeneral: |factor| in [2^(e-2), 2^(e+1)), the sums of min(e-2, 0) and max(e+2, 0).
struct Scan { bool allRegular=true; long long low=0,high=0; };
inline Scan scan(const std::complex<double>* raw,size_t stride,int count,bool unitZero){
    Scan out;
    const int margin=unitZero ? 1 : 2;
    for(int k=0;k<count;++k){
        const std::complex<double> z=raw[size_t(k)*stride];
        out.allRegular&=regular(z);
        const int e=std::max(sizeExponent(z.real()),sizeExponent(z.imag()));
        out.low+=std::min(e-margin,0);out.high+=std::max(e+margin,0);
    }
    return out;
}
// L0 = 1 (cCCE)
__attribute__((noinline)) inline void runUnit(ScaledComplex& product,const std::complex<double>* raw,size_t stride,int count,
                                              const ScaledComplex& inverseZero,const Scan& sc){
    const bool allRegular=sc.allRegular;
    const long long low=sc.low,high=sc.high;
    long double re=product.real,im=product.imag;
    long long exponent=product.exponent;
    int m=0;
    if(allRegular && count<=4096 && scaleExponent(re,im,m) && guardFree(m,low,high)){
        for(int k=0;k<count;++k){
            const std::complex<double> z=raw[size_t(k)*stride];
            const long double fr=z.real(),fi=z.imag();
            const long double r=re*fr-im*fi, i=re*fi+im*fr;
            re=r;im=i;
        }
    }else if(allRegular){
        for(int k=0;k<count;++k){
            const std::complex<double> z=raw[size_t(k)*stride];
            const long double fr=z.real(),fi=z.imag();
            CCEX_SCALED_RUN_LAZY(fr,fi,0);
        }
    }else for(int k=0;k<count;++k){
        const std::complex<double> z=raw[size_t(k)*stride];
        if(regular(z)){const long double fr=z.real(),fi=z.imag();CCEX_SCALED_RUN_LAZY(fr,fi,0);}
        else CCEX_SCALED_RUN_OTHER(z,true);
    }
    product.real=re;product.imag=im;product.exponent=exponent;
}
// L0 != 1 (gCCE): the factor is z * inverse(L0), formed without normalization (A24). inverse(L0) is normalized: where its
// components are at most 1 in size, z * inverse(L0) of a finite double z is at most 2 * DBL_MAX and so finite, and
// assembleFactor's isfinite test of the factor can be left out (an x87 isfinite goes through memory). Where in addition
// its larger component is at least 1/2, |inverse(L0)| is in [1/2, 2^(1/2)) and |factor| in [2^(e-2), 2^(e+1)) (the scan's
// bounds). A factor with a zero component still goes through assembleFactor, and the run then continues with the guard.
__attribute__((noinline)) inline void runGeneral(ScaledComplex& product,const std::complex<double>* raw,size_t stride,int count,
                                                 const ScaledComplex& inverseZero,const Scan& sc){
    const long double zr=inverseZero.real,zi=inverseZero.imag;
    const long long ze=inverseZero.exponent;
    long double re=product.real,im=product.imag;
    long long exponent=product.exponent;
    int k=0;
    if(std::abs(zr)<=1 && std::abs(zi)<=1){
        const bool allRegular=sc.allRegular && std::max(std::abs(zr),std::abs(zi))>=0.5L;
        const long long low=sc.low,high=sc.high;
        int m=0;
        if(allRegular && count<=4096 && scaleExponent(re,im,m) && guardFree(m,low,high)){
            for(;k<count;++k){
                const std::complex<double> z=raw[size_t(k)*stride];
                const long double fr=z.real()*zr-z.imag()*zi, fi=z.real()*zi+z.imag()*zr;
                if(__builtin_expect(!(fr!=0 && fi!=0),0))break;
                const long double r=re*fr-im*fi, i=re*fi+im*fr;
                re=r;im=i;exponent+=ze;
            }
        }
        for(;k<count;++k){
            const std::complex<double> z=raw[size_t(k)*stride];
            if(__builtin_expect(regular(z),1)){
                const long double fr=z.real()*zr-z.imag()*zi, fi=z.real()*zi+z.imag()*zr;
                if(__builtin_expect(fr!=0 && fi!=0,1)){CCEX_SCALED_RUN_LAZY(fr,fi,ze);continue;}
            }
            CCEX_SCALED_RUN_OTHER(z,false);
        }
    }else for(;k<count;++k){
        const std::complex<double> z=raw[size_t(k)*stride];
        if(regular(z)){
            const long double fr=z.real()*zr-z.imag()*zi, fi=z.real()*zi+z.imag()*zr;
            if(fr!=0 && fi!=0 && std::isfinite(fr) && std::isfinite(fi)){CCEX_SCALED_RUN_LAZY(fr,fi,ze);continue;}
        }
        CCEX_SCALED_RUN_OTHER(z,false);
    }
    product.real=re;product.imag=im;product.exponent=exponent;
}
#undef CCEX_SCALED_RUN_LAZY
#undef CCEX_SCALED_RUN_OTHER
}
// the run's scan, computed once for both products of a time point
inline scaled_run::Scan assembleScan(const std::complex<double>* raw,size_t stride,int count,bool unitZero){
    return scaled_run::scan(raw,stride,count,unitZero);
}
inline void assembleRun(ScaledComplex& product,const std::complex<double>* raw,size_t stride,int count,
                        const ScaledComplex& inverseZero,bool unitZero,const scaled_run::Scan& scan){
    if(!ScaledComplex::lazySafe){
        for(int k=0;k<count;++k)scaled_run::other(product,raw[size_t(k)*stride],inverseZero,unitZero);
    }else if(unitZero)scaled_run::runUnit(product,raw,stride,count,inverseZero,scan);
    else scaled_run::runGeneral(product,raw,stride,count,inverseZero,scan);
}
inline void assembleRun(ScaledComplex& product,const std::complex<double>* raw,size_t stride,int count,
                        const ScaledComplex& inverseZero,bool unitZero){
    assembleRun(product,raw,stride,count,inverseZero,unitZero,assembleScan(raw,stride,count,unitZero));
}
#endif
