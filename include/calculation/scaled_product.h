#ifndef CCEX_SCALED_PRODUCT_H
#define CCEX_SCALED_PRODUCT_H
#include <algorithm>
#include <cmath>
#include <complex>
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
// other case keep the old expression. Returns the factor's value when saved is requested (savemode all).
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
#endif
