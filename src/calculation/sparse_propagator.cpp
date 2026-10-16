/**
 * @file sparse_propagator.cpp
 * @brief propagator = sparse: the sparse cluster Hamiltonian and the Krylov expansion of exp(-i H tau).
 */
#include "../../include/calculation/sparse_propagator.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/calculation/small_hamiltonian.h"
#include <algorithm>
#include <cmath>
#include <utility>

// A2 (improvements.md), propagator = sparse. See sparse_propagator.h.

bool sparsePropagatorFor(const Config* cnf, long dim){
    return !strcasecmp(cnf->propagator,"sparse") && dim >= cnf->sparse_min_dim;
}

////////////////////////////////////////////////////////////////////////
// The Hamiltonian in compressed rows
////////////////////////////////////////////////////////////////////////

// y = (H - shift) x; the diagonal entry minus the shift is formed before it multiplies x (exact when the two are close),
// so the shift costs no precision where the state is
void CsrMatrix::apply(const doublec* x, doublec* y, double shift) const {
    for (long r=0; r<n; r++){
        doublec sum(0.0,0.0);
        const long d = diag[r];
        for (long k=rowPtr[r]; k<d; k++) sum += val[k] * x[col[k]];
        sum += (val[d] - shift) * x[r];
        for (long k=d+1; k<rowPtr[r+1]; k++) sum += val[k] * x[col[k]];
        y[r] = sum;
    }
}

namespace {
void finishRow(CsrMatrix& out, long r, std::vector<std::pair<long,doublec>>& row, double& normInf){
    // the terms of one entry summed in the order they were added (stable sort), exact zeros dropped except on the
    // diagonal (kept for the shift of apply)
    row.push_back({r, doublec(0.0,0.0)});
    std::stable_sort(row.begin(), row.end(), [](const std::pair<long,doublec>& a, const std::pair<long,doublec>& b){
        return a.first < b.first; });
    double rowSum = 0.0;
    for (size_t k=0; k<row.size();){
        const long c = row[k].first;
        doublec v(0.0,0.0);
        for (; k<row.size() && row[k].first == c; k++) v += row[k].second;
        if (v.real() == 0.0 && v.imag() == 0.0 && c != r) continue;
        if (c == r) out.diag.push_back(long(out.col.size()));
        out.col.push_back(int(c)); out.val.push_back(v);
        rowSum += std::abs(v);
    }
    out.rowPtr.push_back(long(out.col.size()));
    normInf = std::max(normInf, rowSum);
    row.clear();
}
}

void csrFromLocal(const LocalHamiltonian& H, CsrMatrix& out){
    const int nsite = int(H.dims.size());
    std::vector<long> stride(nsite, 1);
    long n = 1;
    for (int s=nsite-1; s>=0; s--){ stride[s] = n; n *= H.dims[s]; }
    // the nonzero entries of each local row of each term
    struct Entry { int localCol; doublec v; };
    std::vector<std::vector<std::vector<Entry>>> rows(H.terms.size());
    for (size_t t=0; t<H.terms.size(); t++){
        const MatrixXcd& op = H.terms[t].op;
        rows[t].resize(op.rows());
        for (int a=0; a<op.rows(); a++)
            for (int b=0; b<op.cols(); b++)
                if (op(a,b).real() != 0.0 || op(a,b).imag() != 0.0) rows[t][a].push_back({b, op(a,b)});
    }
    out.n = n; out.rowPtr.assign(1, 0); out.diag.clear(); out.col.clear(); out.val.clear(); out.normInf = 0.0;
    std::vector<std::pair<long,doublec>> row;
    for (long r=0; r<n; r++){
        for (size_t t=0; t<H.terms.size(); t++){
            const LocalTerm& term = H.terms[t];
            const int s0 = term.site[0];
            const int d0 = H.dims[s0];
            const int x0 = int((r/stride[s0]) % d0);
            if (term.nsite == 1){
                const long base = r - x0*stride[s0];
                for (const Entry& e: rows[t][x0]) row.push_back({base + e.localCol*stride[s0], e.v});
            } else {
                const int s1 = term.site[1], d1 = H.dims[s1];
                const int x1 = int((r/stride[s1]) % d1);
                const long base = r - x0*stride[s0] - x1*stride[s1];
                for (const Entry& e: rows[t][x0*d1 + x1])
                    row.push_back({base + (e.localCol/d1)*stride[s0] + (e.localCol%d1)*stride[s1], e.v});
            }
        }
        finishRow(out, r, row, out.normInf);
    }
}

void csrFromDense(const MatrixXcd& H, CsrMatrix& out){
    out.n = H.rows(); out.rowPtr.assign(1, 0); out.diag.clear(); out.col.clear(); out.val.clear(); out.normInf = 0.0;
    std::vector<std::pair<long,doublec>> row;
    for (long r=0; r<out.n; r++){
        for (long c=0; c<H.cols(); c++) row.push_back({c, H(r,c)});
        finishRow(out, r, row, out.normInf);
    }
}

double csrExpectation(const CsrMatrix& H, const Eigen::VectorXcd& v){
    Eigen::VectorXcd w(H.n);
    H.apply(v.data(), w.data());
    return v.dot(w).real() / v.squaredNorm();
}

bool csrFiniteHermitian(const CsrMatrix& H, double& asym, double& scale){
    asym = 0.0; scale = 0.0;
    bool finite = true;
    for (long r=0; r<H.n; r++)
        for (long k=H.rowPtr[r]; k<H.rowPtr[r+1]; k++){
            const doublec v = H.val[k];
            if (!std::isfinite(v.real()) || !std::isfinite(v.imag())) finite = false;
            scale = std::max(scale, std::abs(v));
            // the entry (c, r) of row c, 0 when absent
            const long c = H.col[k];
            const int* first = H.col.data() + H.rowPtr[c];
            const int* last  = H.col.data() + H.rowPtr[c+1];
            const int* hit = std::lower_bound(first, last, int(r));
            const doublec w = (hit != last && *hit == r) ? H.val[hit - H.col.data()] : doublec(0.0,0.0);
            asym = std::max(asym, std::abs(v - std::conj(w)));
        }
    return finite && !(scale > 0.0 && asym > 1e-12*scale);
}

namespace {
void checkHamiltonian(const CsrMatrix& H, const char* what){
    double asym, scale;
    if (!csrFiniteHermitian(H, asym, scale)){
        fprintf(stderr,"Error : execution.propagator=sparse : the %s Hamiltonian (dim %ld) is not finite and Hermitian "
                       "(max|H-H^dag| = %.3e, max|H| = %.3e)\n", what, H.n, asym, scale);
        ccexFatal();
    }
}
}

////////////////////////////////////////////////////////////////////////
// exp(-i H t) v by Lanczos
////////////////////////////////////////////////////////////////////////
//
// One substep of length tau from v (|v| = beta): the Lanczos basis V_k of the Krylov space of H and v (fully
// reorthogonalized), T_k = V_k^dag H V_k tridiagonal, and
//     exp(-i H tau) v  ~  beta V_k exp(-i T_k tau) e_1,
// exact when the space is invariant (breakdown). The error is estimated by the leading term of the remainder
// (Saad 1992; Hochbruck, Lubich and Selhofer 1998):
//     err ~ beta b_k tau |e_k^T phi_1(-i T_k tau) e_1|,   phi_1(z) = (e^z - 1)/z,
// with b_k the norm of the next Lanczos vector. The basis does not depend on tau, so a substep that fails its share of
// the tolerance (tol |v| tau / t) is retried shorter on the same basis. The interval is covered by as many substeps as
// that takes; the step after an accepted one grows with the margin left.
void krylovExpv(const CsrMatrix& H, double t, Eigen::VectorXcd& v, double tol, KrylovStats& stats, double shift){
    if (t == 0.0) return;
    const double beta0 = v.norm();
    if (beta0 == 0.0) return;
    const double sign = t > 0.0 ? 1.0 : -1.0;
    const double total = std::abs(t);
    const long n = H.n;
    // Krylov dimension 30; a space of at most 64 is taken whole (the expansion is then exact in one step)
    const int mmax = int(n <= 64 ? n : 30);
    static MatrixXcd V;
    static Eigen::VectorXcd w;
    if (V.rows() != n || V.cols() < mmax+1) V.resize(n, mmax+1);
    w.resize(n);
    std::vector<double> alpha(mmax), beta(mmax);
    double done = 0.0, next = total;
    while (done < total){
        const double remaining = total - done;
        const double b0 = v.norm();
        V.col(0) = v / b0;
        int k = 0;
        bool invariant = false;
        for (int j=0; j<mmax; j++){
            H.apply(V.col(j).data(), w.data(), shift);
            stats.matvecs++;
            const double hv = w.norm();
            alpha[j] = V.col(j).dot(w).real();
            w -= alpha[j] * V.col(j);
            if (j > 0) w -= beta[j-1] * V.col(j-1);
            // Gram-Schmidt against the whole basis, repeated while it removes most of what is left (Daniel, Gragg, Kaufman
            // and Stewart 1976): near an invariant subspace the new vector is mostly rounding, and one pass would leave it
            // far from orthogonal to the basis
            double before = w.norm();
            for (int pass=0; pass<4; pass++){
                const Eigen::VectorXcd c = V.leftCols(j+1).adjoint() * w;
                w.noalias() -= V.leftCols(j+1) * c;
                const double after = w.norm();
                if (after > 0.7 * before) break;
                before = after;
            }
            beta[j] = w.norm();
            k = j+1;
            // Invariant: the whole space, or a remainder at the level of the rounding of H V_j (it is then rounding
            // alone, and normalizing it would put a vector that is not orthogonal to the basis into it: measured on an
            // h-BN gCCE cluster, a remainder of 1e-72 at k = 31 of 63 gave a 5e-6 error)
            if (k == n || beta[j] <= 1e-14 * hv){ invariant = true; break; }
            V.col(j+1) = w / beta[j];
        }
        // The small problem in long double: an eigenvalue of T near 0 (the state's own energy after the shift) next to
        // others of the size of the Zeeman energies would otherwise carry an absolute error of eps |T|, which the phase
        // multiplies by tau.
        typedef Eigen::Matrix<long double,Eigen::Dynamic,Eigen::Dynamic> MatrixL;
        typedef std::complex<long double> complexL;
        MatrixL T = MatrixL::Zero(k,k);
        for (int j=0; j<k; j++){
            T(j,j) = alpha[j];
            if (j+1 < k){ T(j,j+1) = beta[j]; T(j+1,j) = beta[j]; }
        }
        Eigen::SelfAdjointEigenSolver<MatrixL> es(T);
        const MatrixL& Q = es.eigenvectors();
        const Eigen::Matrix<long double,Eigen::Dynamic,1>& lam = es.eigenvalues();
        const double bk = invariant ? 0.0 : beta[k-1];
        double tau = std::min(next, remaining);
        double err = 0.0, allowed = 0.0;
        for (;;){
            allowed = tol * beta0 * tau / total;
            err = 0.0;
            if (bk > 0.0){
                complexL s(0.0L,0.0L);
                for (int l=0; l<k; l++){
                    const long double theta = sign * lam(l) * (long double)tau;
                    // phi_1(-i theta) = (exp(-i theta) - 1) / (-i theta)
                    complexL phi;
                    if (std::abs(theta) < 1e-5L) phi = complexL(1.0L - theta*theta/6.0L, -theta/2.0L);
                    else phi = (complexL(std::cos(theta), -std::sin(theta)) - 1.0L) / complexL(0.0L, -theta);
                    s += Q(k-1,l) * Q(0,l) * phi;
                }
                err = double(b0 * bk * tau * std::abs(s));
            }
            if (err <= allowed) break;
            stats.rejected++;
            tau *= std::min(0.9, std::max(0.1, 0.9 * std::pow(allowed/err, 1.0/k)));
            if (tau < total * 1e-13){
                fprintf(stderr,"Error : execution.propagator=sparse : the Krylov expansion does not converge (dim %ld, t = %.6e, "
                               "tol = %.3e, ||H||_inf = %.6e)\n", n, t, tol, H.normInf);
                ccexFatal();
            }
        }
        // v = b0 V_k Q exp(-i lam tau) Q^T e_1
        std::vector<complexL> phase(k);
        for (int l=0; l<k; l++){
            const long double theta = sign * lam(l) * (long double)tau;
            phase[l] = Q(0,l) * complexL(std::cos(theta), -std::sin(theta));
        }
        Eigen::VectorXcd y(k);
        for (int j=0; j<k; j++){
            complexL sum(0.0L,0.0L);
            for (int l=0; l<k; l++) sum += Q(j,l) * phase[l];
            y(j) = doublec(double(sum.real()), double(sum.imag()));
        }
        v.noalias() = b0 * (V.leftCols(k) * y);
        stats.substeps++;
        done = (tau == remaining) ? total : done + tau;
        next = err > 0.0 ? tau * std::min(5.0, std::max(1.0, 0.9 * std::pow(allowed/err, 1.0/k))) : remaining;
    }
}

////////////////////////////////////////////////////////////////////////
// cCCE, one state
////////////////////////////////////////////////////////////////////////
//
// calCoherenceCceInto's single-state branch computes psi^dag (UdownDagger Uup) psi, with interval k of the sequence
// evolved by H_alpha in the forward chain and H_beta in the backward one for even k, the other way round for odd k:
//     L(t) = <phi_b | phi_a>,  phi_a = ... e^{-i H_b tau_1} e^{-i H_a tau_0} psi,  phi_b = ... e^{-i H_a tau_1} e^{-i H_b tau_0} psi,
// tau_k = t * sequence[k][2]. The first interval's states are carried from one time point to the next (tau_0 grows by
// sequence[0][2] deltat per point); the others start from them at every point.
void calCoherenceCceSparse(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, MatrixXcd* result, double start){
    if (strcasecmp(Config_getQuantity(cnf),"coherence")!=0){
        fprintf(stderr,"Error : calCoherenceCce : quantity is not coherence\n");
        exit(1);
    }
    LocalHamiltonian La, Lb;
    cceClusterTerms(qa, ba, cnf, La, Lb);
    CsrMatrix Ha, Hb;
    csrFromLocal(La, Ha);
    csrFromLocal(Lb, Hb);
    La.terms.clear(); Lb.terms.clear();
    checkHamiltonian(Ha, "cCCE alpha");
    checkHamiltonian(Hb, "cCCE beta");
    const double buildEnd = cnf->profile ? MPI_Wtime() : 0.0;

    const Eigen::VectorXcd psi0 = BathArray_Psi0(ba);
    // H - c for both Hamiltonians, c the state's energy: the Lanczos coefficients stay small numbers instead of
    // differences of Zeeman energies, and exp(-i c t) is the same on both sides of <phi_b|phi_a> and cancels exactly
    const double shift = csrExpectation(Ha, psi0);
    const int npulse = Pulse_getNpulse(pulse);
    double** sequence = Pulse_getSequence(pulse);
    const int nstep = Config_getNstep(cnf);
    const double deltat = (double)Config_getDeltat(cnf);
    const double tol = cnf->sparse_tol;
    KrylovStats stats;
    Eigen::VectorXcd firstA = psi0, firstB = psi0, phiA, phiB;
    double tau0Done = 0.0, tfree = 0.0;
    const double setupEnd = cnf->profile ? MPI_Wtime() : 0.0;
    for (int i=0; i<nstep; i++){
        const double tau0 = tfree * sequence[0][2];
        if (tau0 != tau0Done){
            krylovExpv(Ha, tau0 - tau0Done, firstA, tol, stats, shift);
            krylovExpv(Hb, tau0 - tau0Done, firstB, tol, stats, shift);
            tau0Done = tau0;
        }
        phiA = firstA; phiB = firstB;
        for (int k=1; k<=npulse; k++){
            const double tau = tfree * sequence[k][2];
            krylovExpv(k%2==0 ? Ha : Hb, tau, phiA, tol, stats, shift);
            krylovExpv(k%2==0 ? Hb : Ha, tau, phiB, tol, stats, shift);
        }
        result[i].resize(1,1);
        result[i](0,0) = phiB.dot(phiA);
        tfree += deltat;
    }
    if (cnf->profile){
        computeProfile(cnf,"cce","sparse",int(Ha.n),buildEnd-start,setupEnd-buildEnd,MPI_Wtime()-setupEnd);
        printf("[sparse] rank=%d method=cce dim=%ld nnz=%ld matvecs=%lld substeps=%lld rejected=%lld\n",
               rank, Ha.n, long(Ha.val.size()), stats.matvecs, stats.substeps, stats.rejected);
    }
}

////////////////////////////////////////////////////////////////////////
// gCCE, evolution = vector
////////////////////////////////////////////////////////////////////////

namespace {
// out = (U (x) 1_bdim) in
void applyPulseBlocks(const MatrixXcd& U, int bdim, const Eigen::VectorXcd& in, Eigen::VectorXcd& out){
    const int qdim = int(U.rows());
    out.setZero(in.size());
    for (int a=0; a<qdim; a++)
        for (int b=0; b<qdim; b++)
            if (U(a,b).real() != 0.0 || U(a,b).imag() != 0.0)
                out.segment(long(a)*bdim, bdim) += U(a,b) * in.segment(long(b)*bdim, bdim);
}
}

// exp(-i H tau) v, v split along the eigenvectors e_a of the qubit Hamiltonian: v = sum_a e_a (x) <e_a|v>. Each part
// lies near one qubit level (the levels are far apart next to the bath's energies: the zero-field splitting, the electron
// Zeeman energy) and is expanded around its own energy c_a, exp(-i H tau) v_a = exp(-i c_a tau) exp(-i (H - c_a) tau) v_a,
// the phase in long double. One expansion of the whole v carried the level splitting in every Lanczos coefficient, with
// an absolute rounding of eps times it (3e-8 against eigen on the NV-13C gCCE2 check case).
void GcceSparseWork::evolve(double tau, Eigen::VectorXcd& v, double tol){
    if (tau == 0.0) return;
    sum.setZero(v.size());
    for (int a=0; a<qdim; a++){
        part.setZero(bdim);
        for (int q=0; q<qdim; q++) part += std::conj(levels(q,a)) * v.segment(long(q)*bdim, bdim);
        if (part.squaredNorm() == 0.0) continue;
        piece.resize(v.size());
        for (int q=0; q<qdim; q++) piece.segment(long(q)*bdim, bdim) = levels(q,a) * part;
        const double c = csrExpectation(H, piece);
        krylovExpv(H, tau, piece, tol, stats, c);
        const long double theta = (long double)c * (long double)tau;
        sum += doublec(double(std::cos(theta)), -double(std::sin(theta))) * piece;
    }
    v.swap(sum);
}

// Utot psi0 = Ufree_N P_{N-1} ... P_0 Ufree_0 psi0 (propagateStateGcce); Ufree_0 psi0 carried from the previous time point
const Eigen::VectorXcd& propagateStateGcceSparse(GcceSparseWork& w, const Eigen::VectorXcd& psi0, Pulse* pulse,
                                                 const MatrixXcd* Upulses, double tfree, double tol){
    double** sequence = Pulse_getSequence(pulse);
    const int npulse = Pulse_getNpulse(pulse);
    if (!w.started){ w.first = psi0; w.firstTau = 0.0; w.started = true; }
    const double tau0 = tfree * sequence[0][2];
    if (tau0 != w.firstTau){
        w.evolve(tau0 - w.firstTau, w.first, tol);
        w.firstTau = tau0;
    }
    w.psi = w.first;
    for (int ipulse=1; ipulse<npulse+1; ipulse++){
        applyPulseBlocks(Upulses[ipulse-1], w.bdim, w.psi, w.scratch);
        w.psi.swap(w.scratch);
        w.evolve(tfree * sequence[ipulse][2], w.psi, tol);
    }
    return w.psi;
}
