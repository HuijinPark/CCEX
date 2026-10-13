#ifndef CCEX_SMALL_CLUSTER_MATH_H
#define CCEX_SMALL_CLUSTER_MATH_H
// A26 step 1: the cCCE coherence of one small cluster (bath dimension D = 2 or 4), written once for the batched GPU kernel
// (one thread per cluster, small_cluster_gpu.cu) and for the host (the CPU check of the same arithmetic, test_small_cluster).
// No Eigen, no cuComplex: plain structs, so that nvcc and the host compiler see the same code.
//
// What it computes is calCoherenceCceInto's (simulator_cce.cpp) for one cluster:
//   eigensystems  Halpha = Vu diag(Eu) Vu^dagger, Hbeta = Vd diag(Ed) Vd^dagger (cyclic complex Jacobi; the lower triangle
//                 is read, as Eigen's SelfAdjointEigenSolver does), M = Vu^dagger Vd;
//   a single state psi: L(t) = row . col with row = psi^dagger Vd D_d0^* [ M^dagger D_u1^* | M D_d2^* ] ...,
//                 col = ... [ M D_d1 | M^dagger D_u2 ] M^dagger D_u0 Vu^dagger psi  (evolveSingleState);
//   the ensemble: the same chain averaged over the basis states e_k, i.e. Tr[...]/D - which is what both CPU ensemble
//                 branches give (the CPMG/Hahn/Ramsey products after their even-pulse conjugation, and the Manual one).
// Interval p lasts t * fraction[p]; intervals with the same fraction share their phases (distinct slot same[p]).
// Phases exp(-i E tau): exact = sin and cos at every step; recursive = so every 16 steps and in between the previous
// step's times exp(-i E fraction deltat), as evolution_update = recursive does on the CPU (phaseAnchor 16).
//
// Every loop over the matrix indices has a compile-time trip count and is unrolled, so the matrices stay in registers.
// A26 step 0 found the device eigenvalues wrong for D = 4 when the Jacobi loops were not unrolled (the matrix then lives
// in local memory with run-time indices): A26 step 1 traced it (job 680) to nvcc 12.8's NVVM optimizer (cicc -O3) - after
// the first rotation the device reads the diagonal from before it; the same source is right with -G, with -Xcicc -O0,
// and on the host. The standing GPU check (test_small_cluster_cuda) compares this code on the card with CPU-dumped clusters.
#include <cmath>

#if defined(__CUDACC__)
#define SC_HD __host__ __device__ __forceinline__
#define SC_UNROLL _Pragma("unroll")
#else
#define SC_HD inline
#define SC_UNROLL
#endif

namespace smallcluster {

constexpr int maxIntervals = 65;   // npulse + 1
constexpr int phaseAnchor = 16;    // evolution_update = recursive: sin and cos every 16 steps (simulator_cce.cpp)
constexpr int jacobiSweeps = 20;
// status of a cluster: 0 computed; 1 the Jacobi iteration did not bring the off-diagonal down, or an eigenvalue or
// eigenvector is not finite; 2 a coherence is not finite (the host then computes the cluster on the CPU in both cases)
enum { statusOk = 0, statusNotConverged = 1, statusNotFinite = 2 };

struct C { double re, im; };
SC_HD C cmake(double re, double im){ C z; z.re=re; z.im=im; return z; }
SC_HD C cadd(C a, C b){ return cmake(a.re+b.re, a.im+b.im); }
SC_HD C csub(C a, C b){ return cmake(a.re-b.re, a.im-b.im); }
SC_HD C cmul(C a, C b){ return cmake(a.re*b.re-a.im*b.im, a.re*b.im+a.im*b.re); }
SC_HD C cconj(C a){ return cmake(a.re, -a.im); }
SC_HD C cscale(C a, double s){ return cmake(a.re*s, a.im*s); }
SC_HD double cabs2(C a){ return a.re*a.re+a.im*a.im; }

// The time plan of one calculation: interval fractions, which distinct slot each interval reads, the time grid.
struct Plan {
    int intervals;                    // npulse + 1
    int distinct;                     // distinct fractions (1, 2 or 4 slots are compiled)
    int nstep;
    double deltat;
    double fraction[maxIntervals];    // of the distinct slot: fraction[same[p]] is interval p's
    int same[maxIntervals];
};

// The phase slots compiled for a plan with this many distinct fractions (1, 2 or 4); 0: none fits, the CPU computes.
inline int slotsFor(int distinct){ return distinct<=1 ? 1 : distinct<=2 ? 2 : distinct<=4 ? 4 : 0; }

// exp(-i E tau) = (cos, -sin)
SC_HD C phase(double angle){
#if defined(__CUDA_ARCH__)
    double s, c; sincos(angle, &s, &c); return cmake(c, -s);
#else
    return cmake(std::cos(angle), -std::sin(angle));
#endif
}

// Eigen-decomposition of the Hermitian D x D matrix h (column-major, the lower triangle read): h = V diag(E) V^dagger.
// The rotations run on h scaled by a power of two (exact) so that its largest component is in [0.5, 1): a finite h of any
// size - 1e200 or 1e-200 - neither overflows nor underflows in them (review 20261012-233058-452278 P1-01/P3-02; the
// magnitudes squared of the unscaled h did both and returned wrong eigenvalues as converged); magnitudes are hypot.
// Returns false when the off-diagonal did not come down to 1e-12 of the largest eigenvalue, or when h, an eigenvalue
// (after scaling back) or an eigenvector component is not finite.
template<int D> SC_HD bool jacobi(const C* h, C (&V)[D][D], double (&E)[D]){
    C A[D][D];
    double largest = 0;
    bool finite = true;
    SC_UNROLL for(int j=0;j<D;++j) SC_UNROLL for(int i=0;i<D;++i){
        A[i][j] = i>=j ? h[i+j*D] : cconj(h[j+i*D]);
        finite = finite && std::isfinite(A[i][j].re) && std::isfinite(A[i][j].im);
        largest = std::fmax(largest, std::fmax(std::fabs(A[i][j].re), std::fabs(A[i][j].im)));
    }
    SC_UNROLL for(int i=0;i<D;++i) A[i][i].im = 0;
    SC_UNROLL for(int i=0;i<D;++i) SC_UNROLL for(int j=0;j<D;++j) V[i][j] = cmake(i==j ? 1.0 : 0.0, 0.0);
    if(!finite){
        SC_UNROLL for(int i=0;i<D;++i) E[i] = 0;
        return false;
    }
    int exponent = 0;
    if(largest > 0){
        std::frexp(largest, &exponent);               // largest = f 2^exponent, f in [0.5, 1)
        SC_UNROLL for(int i=0;i<D;++i) SC_UNROLL for(int j=0;j<D;++j)
            A[i][j] = cmake(std::ldexp(A[i][j].re, -exponent), std::ldexp(A[i][j].im, -exponent));
    }
    for(int sweep=0; sweep<jacobiSweeps; ++sweep){
        bool rotated = false;
        SC_UNROLL for(int p=0;p<D-1;++p) SC_UNROLL for(int q=p+1;q<D;++q){
            const C off = A[p][q];
            const double m = std::hypot(off.re, off.im);
            const double ap = A[p][p].re, aq = A[q][q].re;
            if(m > 1e-17*(std::fabs(ap)+std::fabs(aq)) && m > 0){
                rotated = true;
                // the rotation that zeroes A[p][q] (that of the step-0 prototype, checked against Eigen to 1e-13)
                const double g = aq-ap;
                const double t = g==0 ? 1.0 : std::copysign(2*m/(std::fabs(g)+std::hypot(g,2*m)), g);
                const double c = 1/std::sqrt(1+t*t);
                const C s = cscale(off, t*c/m), sc = cconj(s);
                SC_UNROLL for(int k=0;k<D;++k){     // columns: A J, V J
                    const C a=A[k][p], b=A[k][q];
                    A[k][p] = csub(cscale(a,c), cmul(sc,b)); A[k][q] = cadd(cmul(s,a), cscale(b,c));
                    const C va=V[k][p], vb=V[k][q];
                    V[k][p] = csub(cscale(va,c), cmul(sc,vb)); V[k][q] = cadd(cmul(s,va), cscale(vb,c));
                }
                SC_UNROLL for(int k=0;k<D;++k){     // rows: J^dagger A
                    const C a=A[p][k], b=A[q][k];
                    A[p][k] = csub(cscale(a,c), cmul(s,b)); A[q][k] = cadd(cmul(sc,a), cscale(b,c));
                }
                A[p][q] = cmake(0,0); A[q][p] = cmake(0,0);
            }
        }
        if(!rotated) break;
    }
    // the test on the scaled matrix; "!(x <= y)" and the explicit checks so that a NaN fails it (fmax drops a NaN)
    double off = 0, diag = 0;
    SC_UNROLL for(int i=0;i<D;++i){
        finite = finite && std::isfinite(A[i][i].re);
        diag = std::fmax(diag, std::fabs(A[i][i].re));
        SC_UNROLL for(int j=0;j<D;++j){
            finite = finite && std::isfinite(V[i][j].re) && std::isfinite(V[i][j].im);
            if(i!=j){
                const double a = std::hypot(A[i][j].re, A[i][j].im);
                finite = finite && std::isfinite(a);
                off = std::fmax(off, a);
            }
        }
    }
    SC_UNROLL for(int i=0;i<D;++i){
        E[i] = std::ldexp(A[i][i].re, exponent);
        finite = finite && std::isfinite(E[i]);
    }
    return finite && off <= 1e-12*diag;
}

// Phases of the distinct slots at step i (tfree = the step's time, as the CPU sums it)
template<int D,int NF,bool REC> SC_HD void phases(const Plan& plan, int i, double tfree, const double (&Eu)[D], const double (&Ed)[D],
                                                C (&up)[NF][D], C (&down)[NF][D], C (&stepUp)[NF][D], C (&stepDown)[NF][D]){
    SC_UNROLL for(int f=0;f<NF;++f){
        if(f >= plan.distinct) break;
        const double frac = plan.fraction[f];
        if(!REC || i % phaseAnchor == 0){
            const double tau = tfree*frac;
            SC_UNROLL for(int k=0;k<D;++k){ up[f][k] = phase(Eu[k]*tau); down[f][k] = phase(Ed[k]*tau); }
            if(REC && i == 0){
                const double dtau = plan.deltat*frac;
                SC_UNROLL for(int k=0;k<D;++k){ stepUp[f][k] = phase(Eu[k]*dtau); stepDown[f][k] = phase(Ed[k]*dtau); }
            }
        } else {
            SC_UNROLL for(int k=0;k<D;++k){ up[f][k] = cmul(up[f][k], stepUp[f][k]); down[f][k] = cmul(down[f][k], stepDown[f][k]); }
        }
    }
}

// The interval chain for one start (r0 = psi^dagger Vd, c0 = Vu^dagger psi): row . col
template<int D,int NF> SC_HD C chain(const Plan& plan, const C (&M)[D][D], const C (&r0)[D], const C (&c0)[D],
                                    const C (&up)[NF][D], const C (&down)[NF][D]){
    C row[D], col[D], tmp[D], u[D], d[D];
    for(int p=0; p<plan.intervals; ++p){
        const int f = plan.same[p];
        SC_UNROLL for(int ff=0; ff<NF; ++ff) if(ff==f){ SC_UNROLL for(int k=0;k<D;++k){ u[k]=up[ff][k]; d[k]=down[ff][k]; } }
        if(p==0){
            SC_UNROLL for(int k=0;k<D;++k){ row[k] = cmul(r0[k], cconj(d[k])); tmp[k] = cmul(u[k], c0[k]); }
            SC_UNROLL for(int i=0;i<D;++i){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(cconj(M[k][i]), tmp[k])); col[i]=s; }
        } else if(p%2==1){
            // row <- (row M^dagger) .* conj(Du); col <- M (Dd .* col)
            SC_UNROLL for(int j=0;j<D;++j){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(row[k], cconj(M[j][k]))); tmp[j]=s; }
            SC_UNROLL for(int j=0;j<D;++j) row[j] = cmul(tmp[j], cconj(u[j]));
            SC_UNROLL for(int k=0;k<D;++k) tmp[k] = cmul(d[k], col[k]);
            SC_UNROLL for(int i=0;i<D;++i){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(M[i][k], tmp[k])); col[i]=s; }
        } else {
            // row <- (row M) .* conj(Dd); col <- M^dagger (Du .* col)
            SC_UNROLL for(int j=0;j<D;++j){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(row[k], M[k][j])); tmp[j]=s; }
            SC_UNROLL for(int j=0;j<D;++j) row[j] = cmul(tmp[j], cconj(d[j]));
            SC_UNROLL for(int k=0;k<D;++k) tmp[k] = cmul(u[k], col[k]);
            SC_UNROLL for(int i=0;i<D;++i){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(cconj(M[k][i]), tmp[k])); col[i]=s; }
        }
    }
    C L = cmake(0,0);
    SC_UNROLL for(int k=0;k<D;++k) L = cadd(L, cmul(row[k], col[k]));
    return L;
}

// One cluster. in: Halpha (D*D), Hbeta (D*D), then psi (D) for a single state; out: plan.nstep coherences.
// Returns the status.
template<int D,bool ENS,int NF,bool REC> SC_HD int cluster(const Plan& plan, const C* in, C* out){
    C Vu[D][D], Vd[D][D], M[D][D];
    double Eu[D], Ed[D];
    const bool okUp = jacobi<D>(in, Vu, Eu);
    const bool okDown = jacobi<D>(in+D*D, Vd, Ed);
    if(!okUp || !okDown){
        for(int t=0;t<plan.nstep;++t) out[t] = cmake(0,0);
        return statusNotConverged;
    }
    SC_UNROLL for(int i=0;i<D;++i) SC_UNROLL for(int j=0;j<D;++j){
        C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(cconj(Vu[k][i]), Vd[k][j])); M[i][j]=s;
    }
    C r0[D], c0[D];                         // psi^dagger Vd and Vu^dagger psi (a single state)
    if(!ENS){
        const C* psi = in+2*D*D;
        SC_UNROLL for(int j=0;j<D;++j){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(cconj(psi[k]), Vd[k][j])); r0[j]=s; }
        SC_UNROLL for(int i=0;i<D;++i){ C s=cmake(0,0); SC_UNROLL for(int k=0;k<D;++k) s=cadd(s, cmul(cconj(Vu[k][i]), psi[k])); c0[i]=s; }
    }
    C up[NF][D], down[NF][D], stepUp[NF][D], stepDown[NF][D];
    double tfree = 0.0;
    for(int i=0;i<plan.nstep;++i){
        phases<D,NF,REC>(plan, i, tfree, Eu, Ed, up, down, stepUp, stepDown);
        C L;
        if(!ENS){
            L = chain<D,NF>(plan, M, r0, c0, up, down);
        } else if(plan.intervals == 2 && plan.same[0] == plan.same[1]){
            // Hahn echo: Tr[X^dagger Du^* X Du]/D with X = M Dd M^dagger, i.e. sum_ab |X_ab|^2 conj(Du_a) Du_b / D -
            // the identity of the CPU's ensemble Hahn branch (D^3 per step instead of D chains)
            const C (&u)[D] = up[0];        // slots are numbered by first appearance: interval 0 has slot 0
            const C (&d)[D] = down[0];
            C sum = cmake(0,0);
            SC_UNROLL for(int a=0;a<D;++a){
                C rowSum = cmake(0,0);
                SC_UNROLL for(int b=0;b<D;++b){
                    C x = cmake(0,0);
                    SC_UNROLL for(int k=0;k<D;++k) x = cadd(x, cmul(cmul(M[a][k], d[k]), cconj(M[b][k])));
                    rowSum = cadd(rowSum, cscale(u[b], cabs2(x)));
                }
                sum = cadd(sum, cmul(cconj(u[a]), rowSum));
            }
            L = cscale(sum, 1.0/D);
        } else {
            L = cmake(0,0);
            SC_UNROLL for(int k=0;k<D;++k){
                C rk[D], ck[D];   // e_k^dagger Vd and Vu^dagger e_k
                SC_UNROLL for(int j=0;j<D;++j){ rk[j] = Vd[k][j]; ck[j] = cconj(Vu[k][j]); }
                L = cadd(L, chain<D,NF>(plan, M, rk, ck, up, down));
            }
            L = cscale(L, 1.0/D);
        }
        if(!std::isfinite(L.re) || !std::isfinite(L.im)){
            for(int t=0;t<plan.nstep;++t) out[t] = cmake(0,0);
            return statusNotFinite;
        }
        out[i] = L;
        tfree += plan.deltat;
    }
    return statusOk;
}

} // namespace smallcluster
#endif
