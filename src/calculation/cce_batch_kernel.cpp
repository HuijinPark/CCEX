/**
 * @file cce_batch_kernel.cpp
 * @brief The batched time evolution of small cCCE clusters on the CPU (compiled without contraction).
 */
// A30: the batched time evolution of small cCCE clusters (cce_batch.h). Compiled with -ffp-contract=off (both builds):
// every fused multiply-add here is written out, nothing else is fused. No Eigen in this file: an Eigen function
// instantiated here with these flags could be the copy the linker keeps for the scalar path too.
//
// Per cluster the operations of the scalar path in its order (simulator_cce.cpp evolveSingleState and the ensemble Hahn
// branch, simulator_cce_recursive.cpp), as the build computes them:
//   - Eigen's complex packet product pmul(a, b): re = ar br - ai bi, im = ar bi + ai br; unfused for Packet2cd (AVX,
//     the size-2 objects everywhere and the size-4 ones with AVX2), fused (_mm512_fmaddsub_pd: re = fma(ar,br,-(ai bi)),
//     im = fma(ar,bi,ai br)) for Packet4cd, which AVX-512 uses for the size-4 objects. A product with a conjugate is
//     pmul(a, pconj(b)).
//   - colTimes (exact): Eigen's column-major gemv: per row, from +0, (A(r,j) s_j) + c in order of j, then (c alpha) + 0
//     with alpha = (1,0) as a packet. recursive: the fixed-size lazy product, A(r,0) s_0, then (A(r,k) s_k) + c.
//   - rowTimes: Eigen's row-major gemv: AVX2 sums the even and the odd j in two lanes from +0 and adds them; one Packet4cd
//     (AVX-512, D 4) adds each product to +0 and sums (0+2)+(1+3); AVX-512 with D 2 has no full packet and sums the
//     scalar products onto +0 in order. Then res = 0 + alpha*cc as std::complex.
//   - the final row * col: (p0+p2)+(p1+p3) for D 4 (two Packet2cd or one Packet4cd), p0+p1 for D 2 (exact with AVX-512
//     as a scalar sum of scalar products, the same bits).
//   - sin and cos as scalar libm calls (a vector libm would change bits).
//   - the scalar std::complex products of the scalar path (the recursive phase step, alpha*cc, the AVX-512 D 2 row gemv,
//     the ensemble's M Dd and the cast gemv): the compiler fuses them or not site by site (table below).
//   - the ensemble Hahn echo (both the exact branch and A18's recursive one): X = (M Dd) M^dagger coefficient by
//     coefficient in order, |X|^2 = fma(re,re,im im), t_c = sum_k |X(k,c)|^2 conj(Uu_k) in order then 0 + alpha t,
//     (t Uu) in pairs, and / (D,0) as libgcc's __divdc3 does it for finite values (others through that division itself);
//     with AVX-512 the products of X, of (t Uu) and of the recursive phase step are Packet4cd's (fused) for D 4 (found on
//     idm by trying the combinations against the scalar path, A30 step 2).
// The step-1 study and the tables: (internal notes)
#include "../../include/calculation/cce_batch.h"
#include <cfloat>
#include <cmath>

// The builds whose arithmetic is reproduced. Anything else - another compiler or version, fast math, an instruction set
// not listed - has no kernel and every cluster goes through the scalar path.
#if !defined(__FAST_MATH__) && defined(__AVX2__) && defined(__FMA__)
#  if defined(__INTEL_LLVM_COMPILER)
#    if __INTEL_LLVM_COMPILER == 20230200 && defined(__AVX512F__)
#      define CCEX_BATCH_VARIANT 3          // icpx 2023.2, AVX-512 (idm4)
#    elif __INTEL_LLVM_COMPILER == 20260101 && !defined(__AVX512F__)
#      define CCEX_BATCH_VARIANT 2          // icpx 2026.1, AVX2 (wsl-icpx)
#    endif
#  elif defined(__GNUC__) && !defined(__clang__) && !defined(__INTEL_COMPILER) && __GNUC__ == 13 && !defined(__AVX512F__)
#    define CCEX_BATCH_VARIANT 1            // g++ 13, AVX2 (wsl-gnu)
#  endif
#endif
#ifndef CCEX_BATCH_VARIANT
#  define CCEX_BATCH_VARIANT 0
#endif

#if CCEX_BATCH_VARIANT == 0

namespace ccebatch {
int lanes(){ return 0; }
const char* variant(){ return "none"; }
bool covers(const Plan&, int){ return false; }
void evolve(const Plan&, int, const Pack&, std::complex<double>* const*){}
}

#else

#include <immintrin.h>
// No contraction here even when the build forgets -ffp-contract=off (a Makefile generated before A30 lacks the rule: g++
// then fused the kernel's own products and 96 % of the outputs moved, found 2026-10-13 in check_rdip_ladder's pCCE cases)
#if defined(__clang__)
#  pragma clang fp contract(off)
#else
#  pragma GCC optimize("fp-contract=off")
#endif

namespace {
#if CCEX_BATCH_VARIANT == 3
typedef __m512d V; constexpr int B = 8;
inline V vadd(V a, V b){ return _mm512_add_pd(a,b); }
inline V vsub(V a, V b){ return _mm512_sub_pd(a,b); }
inline V vmul(V a, V b){ return _mm512_mul_pd(a,b); }
inline V vdiv(V a, V b){ return _mm512_div_pd(a,b); }
inline V vfmadd(V a, V b, V c){ return _mm512_fmadd_pd(a,b,c); }    // a b + c, one rounding
inline V vfmsub(V a, V b, V c){ return _mm512_fmsub_pd(a,b,c); }    // a b - c, one rounding
inline V vset1(double x){ return _mm512_set1_pd(x); }
inline V vneg(V a){ return _mm512_castsi512_pd(_mm512_xor_si512(_mm512_castpd_si512(a), _mm512_set1_epi64(0x8000000000000000LL))); }
inline V vload(const double* p){ return _mm512_load_pd(p); }
inline void vstore(double* p, V a){ _mm512_store_pd(p,a); }
constexpr bool fusedPacket4 = true;     // Packet4cd
#else
typedef __m256d V; constexpr int B = 4;
inline V vadd(V a, V b){ return _mm256_add_pd(a,b); }
inline V vsub(V a, V b){ return _mm256_sub_pd(a,b); }
inline V vmul(V a, V b){ return _mm256_mul_pd(a,b); }
inline V vdiv(V a, V b){ return _mm256_div_pd(a,b); }
inline V vfmadd(V a, V b, V c){ return _mm256_fmadd_pd(a,b,c); }
inline V vfmsub(V a, V b, V c){ return _mm256_fmsub_pd(a,b,c); }
inline V vset1(double x){ return _mm256_set1_pd(x); }
inline V vneg(V a){ return _mm256_xor_pd(a, _mm256_set1_pd(-0.0)); }
inline V vload(const double* p){ return _mm256_load_pd(p); }
inline void vstore(double* p, V a){ _mm256_store_pd(p,a); }
constexpr bool fusedPacket4 = false;
#endif
constexpr int phaseAnchor = 16;   // simulator_cce.cpp, simulator_cce_recursive.cpp

struct CV { V re, im; };
inline CV cadd(CV a, CV b){ return { vadd(a.re,b.re), vadd(a.im,b.im) }; }
inline CV cconj(CV b){ return { b.re, vneg(b.im) }; }
inline CV czero(){ return { vset1(0.0), vset1(0.0) }; }
// Eigen's packet product: unfused (Packet2cd's addsub) or fused (Packet4cd's fmaddsub)
template<bool F> inline CV pmul(CV a, CV b){
    if (F) return { vfmsub(a.re,b.re,vmul(a.im,b.im)), vfmadd(a.re,b.im,vmul(a.im,b.re)) };
    return { vsub(vmul(a.re,b.re), vmul(a.im,b.im)), vadd(vmul(a.re,b.im), vmul(a.im,b.re)) };
}
// A scalar std::complex product a*b of the scalar path as the compiler made it: unfused (0), or re = fma(ar,br,-(ai bi))
// and im = fma(ar,bi,ai br) (1) or fma(ai,br,ar bi) (3).
inline CV smul(int v, CV a, CV b){
    if (v == 0) return pmul<false>(a,b);
    const V re = vfmsub(a.re,b.re,vmul(a.im,b.im));
    const V im = v == 1 ? vfmadd(a.re,b.im,vmul(a.im,b.re)) : vfmadd(a.im,b.re,vmul(a.re,b.im));
    return { re, im };
}
// Which of these each site of the scalar path is (the compilers decide it site by site):
//   siteUp/siteDown: evolveSingleState's up[k] *= stepUp[k], down[k] *= stepDown[k] (fixed size, unrolled);
//   siteAlpha: the gemv's res[i] += alpha*cc0; siteRow: the row-major gemv's cc0 += cj.pmul(lhs(i,j), b0) (AVX-512, D 2);
//   siteEns: the ensemble's M(r,k) Dd(k) and the cast gemv's |X(k,c)|^2 conj(Uu_k).
// g++ 13: read from its object code (A30 step 1); icpx 2026.1 (AVX2): none fused (step 1); icpx 2023.2 (AVX-512): the
// row-major gemv's product is fused (variant 3), the phase step is not (A30 step 2, idm).
#if CCEX_BATCH_VARIANT == 1
template<int D> constexpr int siteUp(int k){ return D == 2 ? 1 : (k % 2 == 0 ? 3 : 1); }
template<int D> constexpr int siteDown(int k){ return k == D-1 ? 3 : 1; }
constexpr int siteAlpha = 1, siteRow = 0, siteEns = 1;
#elif CCEX_BATCH_VARIANT == 2
template<int D> constexpr int siteUp(int){ return 0; }
template<int D> constexpr int siteDown(int){ return 0; }
constexpr int siteAlpha = 0, siteRow = 0, siteEns = 0;
#else
template<int D> constexpr int siteUp(int){ return 0; }
template<int D> constexpr int siteDown(int){ return 0; }
constexpr int siteAlpha = 0, siteRow = 3, siteEns = 0;
#endif
// res = 0 + alpha*cc, alpha = (1,0) as std::complex (the gemv's res[i] += alpha*cc0)
inline CV alphaInto(CV cc){
    const CV alpha{ vset1(1.0), vset1(0.0) };
    return cadd(czero(), smul(siteAlpha,alpha,cc));
}

inline CV load(const double (*re)[ccebatch::maxLanes], const double (*im)[ccebatch::maxLanes], int e){
    return { vload(re[e]), vload(im[e]) };
}
// exp(-i E tau) per lane: scalar cos and sin, as the scalar path forms doublec(cos(a), -sin(a))
template<int D> inline void phases(const double (*E)[ccebatch::maxLanes], double tau, CV (&u)[D]){
    alignas(64) double re[D][B], im[D][B];
    for (int k=0; k<D; ++k)
        for (int l=0; l<B; ++l){
            const double a = E[k][l]*tau;
            re[k][l] = std::cos(a); im[k][l] = -std::sin(a);
        }
    for (int k=0; k<D; ++k) u[k] = { vload(re[k]), vload(im[k]) };
}
inline void storeOut(CV v, std::complex<double>* const* out, int i){
    alignas(64) double re[B], im[B];
    vstore(re,v.re); vstore(im,v.im);
    for (int l=0; l<B; ++l) if (out[l]) out[l][i] = std::complex<double>(re[l],im[l]);
}

// col = A s: exact, Eigen's column-major gemv (colTimes); recursive, the fixed-size lazy product
template<int D, bool F> inline void colGemv(const CV (&A)[D*D], const CV (&s)[D], CV (&col)[D]){
    const CV alpha{ vset1(1.0), vset1(0.0) };
    for (int r=0; r<D; ++r){
        CV c = czero();
        for (int j=0; j<D; ++j) c = cadd(pmul<F>(A[r*D+j],s[j]), c);
        col[r] = cadd(pmul<F>(c,alpha), czero());
    }
}
template<int D, bool F> inline void colLazy(const CV (&A)[D*D], const CV (&s)[D], CV (&col)[D]){
    for (int r=0; r<D; ++r){
        CV c = pmul<F>(A[r*D],s[0]);
        for (int k=1; k<D; ++k) c = cadd(pmul<F>(A[r*D+k],s[k]), c);
        col[r] = c;
    }
}
// out = in A (rowTimes, Eigen's row-major gemv over A's storage: out_i = sum_j A(j,i) in_j)
template<int D> inline void rowGemv(const CV (&in)[D], const CV (&A)[D*D], CV (&out)[D]){
    for (int i=0; i<D; ++i){
        CV cc;
        if (B == 8 && D == 4){              // one Packet4cd: each product + 0, then (0+2)+(1+3)
            CV q[4];
            for (int j=0; j<4; ++j) q[j] = cadd(pmul<true>(A[j*D+i],in[j]), czero());
            cc = cadd(cadd(q[0],q[2]), cadd(q[1],q[3]));
        } else if (B == 8){                  // D 2: no full packet; the scalar products onto +0 in order
            cc = czero();
            for (int j=0; j<D; ++j) cc = cadd(cc, smul(siteRow,A[j*D+i],in[j]));
        } else {                             // Packet2cd: the even and the odd j from +0, then added
            CV l0 = czero(), l1 = czero();
            for (int j=0; j<D; j+=2){ l0 = cadd(pmul<false>(A[j*D+i],in[j]),l0); l1 = cadd(pmul<false>(A[(j+1)*D+i],in[j+1]),l1); }
            cc = cadd(l0,l1);
        }
        out[i] = alphaInto(cc);
    }
}
// row * col: (p0+p2)+(p1+p3), or p0+p1
template<int D, bool F> inline CV dot(const CV (&r)[D], const CV (&c)[D]){
    if (D == 4) return cadd(cadd(pmul<F>(r[0],c[0]),pmul<F>(r[2],c[2])), cadd(pmul<F>(r[1],c[1]),pmul<F>(r[3],c[3])));
    return cadd(pmul<F>(r[0],c[0]),pmul<F>(r[1],c[1]));
}

// One state (evolveSingleState)
template<int D, bool REC>
void stateKernel(const ccebatch::Plan& pl, const ccebatch::Pack& P, std::complex<double>* const* out){
    constexpr bool F = D == 4 && fusedPacket4;
    CV M[D*D], Madj[D*D], r0[D], c0[D];
    for (int e=0; e<D*D; ++e){ M[e] = load(P.Mre,P.Mim,e); Madj[e] = load(P.Are,P.Aim,e); }
    for (int k=0; k<D; ++k){ r0[k] = load(P.r0re,P.r0im,k); c0[k] = load(P.c0re,P.c0im,k); }
    CV up[ccebatch::maxIntervals][D], down[ccebatch::maxIntervals][D];
    CV stepUp[REC ? ccebatch::maxIntervals : 1][D], stepDown[REC ? ccebatch::maxIntervals : 1][D];
    CV row[D], rowNext[D], col[D], sc[D];
    double tfree = 0.0;
    for (int i=0; i<pl.nstep; ++i){
        for (int p=0; p<pl.intervals; ++p){
            const int same = pl.same[p];
            CV (&Uu)[D] = up[same];
            CV (&Ud)[D] = down[same];
            if (same == p){
                if (!REC || i % phaseAnchor == 0){
                    const double tau = tfree * pl.fraction[p];
                    phases<D>(P.Eu,tau,Uu); phases<D>(P.Ed,tau,Ud);
                    if (REC && i == 0){
                        const double dtau = pl.deltat * pl.fraction[p];
                        phases<D>(P.Eu,dtau,stepUp[same]); phases<D>(P.Ed,dtau,stepDown[same]);
                    }
                } else if (REC){
                    for (int k=0; k<D; ++k){
                        Uu[k] = smul(siteUp<D>(k),Uu[k],stepUp[same][k]);
                        Ud[k] = smul(siteDown<D>(k),Ud[k],stepDown[same][k]);
                    }
                }
            }
            if (p == 0){
                for (int k=0; k<D; ++k){ row[k] = pmul<F>(r0[k],cconj(Ud[k])); sc[k] = pmul<F>(Uu[k],c0[k]); }
                if (REC) colLazy<D,F>(Madj,sc,col); else colGemv<D,F>(Madj,sc,col);
            } else if (p % 2 == 1){
                rowGemv<D>(row,Madj,rowNext);
                for (int k=0; k<D; ++k){ row[k] = pmul<F>(rowNext[k],cconj(Uu[k])); sc[k] = pmul<F>(Ud[k],col[k]); }
                if (REC) colLazy<D,F>(M,sc,col); else colGemv<D,F>(M,sc,col);
            } else {
                rowGemv<D>(row,M,rowNext);
                for (int k=0; k<D; ++k){ row[k] = pmul<F>(rowNext[k],cconj(Ud[k])); sc[k] = pmul<F>(Uu[k],col[k]); }
                if (REC) colLazy<D,F>(Madj,sc,col); else colGemv<D,F>(Madj,sc,col);
            }
        }
        storeOut(dot<D,F>(row,col),out,i);
        tfree += pl.deltat;
    }
}

// The ensemble Hahn echo: Tr[X Du X^dagger Du^*]/D with X = M Dd M^dagger (the exact branch of calCoherenceCceInto and
// A18's recursive one, tau = t/2)
template<int D, bool REC>
void ensembleKernel(const ccebatch::Plan& pl, const ccebatch::Pack& P, std::complex<double>* const* out){
    constexpr bool F = D == 4 && fusedPacket4;      // the dynamic-size products of 4 (Packet4cd with AVX-512)
    CV M[D*D], Uu[D], Ud[D], sUu[D], sUd[D], A[D*D], t[D];
    V X2[D*D];
    for (int e=0; e<D*D; ++e) M[e] = load(P.Mre,P.Mim,e);
    const V zero = vset1(0.0), dim = vset1(double(D));
    double tfree = 0.0;
    for (int i=0; i<pl.nstep; ++i){
        if (!REC || i % phaseAnchor == 0){
            phases<D>(P.Eu,tfree/2.0,Uu); phases<D>(P.Ed,tfree/2.0,Ud);
            if (REC && i == 0){ phases<D>(P.Eu,pl.deltat/2.0,sUu); phases<D>(P.Ed,pl.deltat/2.0,sUd); }
        } else {
            for (int k=0; k<D; ++k){ Uu[k] = pmul<F>(Uu[k],sUu[k]); Ud[k] = pmul<F>(Ud[k],sUd[k]); }
        }
        for (int r=0; r<D; ++r) for (int k=0; k<D; ++k) A[r*D+k] = smul(siteEns,M[r*D+k],Ud[k]);
        for (int r=0; r<D; ++r) for (int c=0; c<D; ++c){
            CV x = pmul<F>(A[r*D],cconj(M[c*D]));
            for (int k=1; k<D; ++k) x = cadd(pmul<F>(A[r*D+k],cconj(M[c*D+k])),x);
            X2[r*D+c] = vfmadd(x.re,x.re,vmul(x.im,x.im));
        }
        for (int c=0; c<D; ++c){
            CV s = smul(siteEns,CV{X2[c],zero},cconj(Uu[0]));
            for (int k=1; k<D; ++k) s = cadd(s,smul(siteEns,CV{X2[k*D+c],zero},cconj(Uu[k])));
            t[c] = alphaInto(s);
        }
        const CV tu = dot<D,F>(t,Uu);
        // / (D,0) as __divdc3 for finite values: x = (a + 0*b)/D, y = (b - 0*a)/D
        const CV o{ vdiv(vadd(tu.re,vmul(zero,tu.im)),dim), vdiv(vsub(tu.im,vmul(zero,tu.re)),dim) };
        storeOut(o,out,i);
        alignas(64) double ar[B], br[B];
        vstore(ar,tu.re); vstore(br,tu.im);
        for (int l=0; l<B; ++l)
            if (out[l] && (!std::isfinite(ar[l]) || !std::isfinite(br[l]))){
                volatile double d = D, z = 0.0;      // the library division, as the scalar path calls it
                out[l][i] = std::complex<double>(ar[l],br[l]) / std::complex<double>(d,z);
            }
        tfree += pl.deltat;
    }
}

template<int D>
void evolveD(const ccebatch::Plan& pl, const ccebatch::Pack& P, std::complex<double>* const* out){
    if (pl.ensemble){ if (pl.recursive) ensembleKernel<D,true>(pl,P,out); else ensembleKernel<D,false>(pl,P,out); }
    else            { if (pl.recursive) stateKernel<D,true>(pl,P,out);    else stateKernel<D,false>(pl,P,out); }
}
}

namespace ccebatch {
int lanes(){ return B; }
const char* variant(){
    return CCEX_BATCH_VARIANT == 1 ? "g++ 13 AVX2" : CCEX_BATCH_VARIANT == 2 ? "icpx 2026.1 AVX2" : "icpx 2023.2 AVX-512";
}
bool covers(const Plan& plan, int D){
    if (D != 2 && D != 4) return false;
    if (plan.ensemble) return true;               // the caller gives the ensemble only for the Hahn echo
    if (plan.intervals < 1 || plan.intervals > maxIntervals) return false;
    for (int p=0; p<plan.intervals; ++p) if (plan.same[p] < 0 || plan.same[p] > p) return false;
    return true;
}
void evolve(const Plan& plan, int D, const Pack& pack, std::complex<double>* const* out){
    if (D == 4) evolveD<4>(plan,pack,out); else evolveD<2>(plan,pack,out);
}
}
#endif
