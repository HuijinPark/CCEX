#ifndef CCEX_COMPENSATED_H
#define CCEX_COMPENSATED_H
#include <cmath>
#ifdef __CUDACC__
#define CCEX_HD __host__ __device__
#else
#define CCEX_HD
#endif

// Two-double arithmetic for the short, cancellation-sensitive Rayleigh sums.
// Hamiltonians, eigenvectors and evolution buffers remain complex double.
struct DoublePair { double high,low; };
CCEX_HD inline DoublePair pairAdd(DoublePair a,DoublePair b){
    const double sum=a.high+b.high,back=sum-a.high;
    const double error=((a.high-(sum-back))+(b.high-back))+a.low+b.low;
    const double high=sum+error;
    return {high,error-(high-sum)};
}
CCEX_HD inline DoublePair pairProduct(double a,double b){
    const double high=a*b;return {high,::fma(a,b,-high)};
}
CCEX_HD inline DoublePair pairScale(DoublePair a,double b){
    return pairAdd(pairProduct(a.high,b),{a.low*b,0});
}
CCEX_HD inline DoublePair rayleighPair(double ar,double ai,double hr,double hi,double br,double bi){
    const auto real=pairAdd(pairProduct(ar,br),pairProduct(ai,bi));
    const auto imag=pairAdd(pairProduct(ar,bi),pairProduct(-ai,br));
    return pairAdd(pairScale(real,hr),pairScale(imag,-hi));
}
CCEX_HD inline DoublePair pairDifference(DoublePair a,DoublePair b){
    return pairAdd(a,{-b.high,-b.low});
}
#undef CCEX_HD
#endif
