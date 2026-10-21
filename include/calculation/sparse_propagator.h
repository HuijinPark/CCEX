/**
 * @file sparse_propagator.h
 * @brief propagator = sparse: exp(-i H tau) on a state vector by a Krylov expansion, without diagonalizing H.
 */
#ifndef __CCEX_SPARSE_PROPAGATOR_H_
#define __CCEX_SPARSE_PROPAGATOR_H_
#include "calculation/simulator.h"

/**
 * propagator = sparse: exp(-i H tau) applied to a state vector without diagonalizing H. The cluster
 * Hamiltonian is assembled as a sparse (CSR) matrix from its one- and two-spin terms - never as a dense matrix - and the
 * action of the exponential is a Lanczos (Krylov) expansion with an a-posteriori error estimate per substep
 * (sparse_propagator.cpp). Paths: the single-state cCCE (nstate > 0) and gCCE with evolution = vector, for clusters of
 * dimension >= sparse_min_dim; smaller clusters keep the eigen path. main.cpp refuses every other combination
 * (Config_validatePropagator).
 */

/** One term of a cluster Hamiltonian acting on one site or two (op on the space of those sites, kron order). */
struct LocalTerm { int nsite; int site[2]; MatrixXcd op; };
/** A Hamiltonian as a sum of local terms on sites of the given dimensions (site 0 is the most significant digit). */
struct LocalHamiltonian { std::vector<int> dims; std::vector<LocalTerm> terms; };

/** The two conditional Hamiltonians of the cCCE (bath + ms_a,b * secular hyperfine) as local terms (small_hamiltonian.cpp). */
void cceClusterTerms(QubitArray* qa, BathArray* ba, Config* cnf, LocalHamiltonian& Halpha, LocalHamiltonian& Hbeta);
/** The gCCE Htot (one qubit, nspin > 0) as local terms, site 0 the qubit; Hq as in gcceClusterHamiltonian. */
void gcceClusterTerms(QubitArray* qa, BathArray* ba, Config* cnf, LocalHamiltonian& Htot, MatrixXcd& Hq);

/** A Hermitian matrix in compressed rows. */
struct CsrMatrix {
    long n = 0;
    std::vector<long> rowPtr, diag;   // diag: the position of each row's diagonal entry (always stored)
    std::vector<int> col;
    std::vector<doublec> val;
    double normInf = 0.0;          // largest row sum of |H_ij|
    void apply(const doublec* x, doublec* y, double shift = 0.0) const;   // y = (H - shift) x
};
void csrFromLocal(const LocalHamiltonian& H, CsrMatrix& out);
void csrFromDense(const MatrixXcd& H, CsrMatrix& out);
/** <v|H|v> / <v|v>. */
double csrExpectation(const CsrMatrix& H, const Eigen::VectorXcd& v);
/** Finite and Hermitian to 1e-12 of the largest entry (the test of validateGcceHamiltonian); else asym and scale set. */
bool csrFiniteHermitian(const CsrMatrix& H, double& asym, double& scale);

/** Counters of the Krylov expansions of one cluster (matrix-vector products, accepted substeps, rejected trials). */
struct KrylovStats { long long matvecs = 0, substeps = 0, rejected = 0; };
/** v <- exp(-i (H - shift) t) v; the estimated error of the whole interval at most tol * |v|. */
void krylovExpv(const CsrMatrix& H, double t, Eigen::VectorXcd& v, double tol, KrylovStats& stats, double shift = 0.0);

/** Whether a cluster of this Hilbert dimension takes the sparse path. */
bool sparsePropagatorFor(const Config* cnf, long dim);
/** calCoherenceCceInto for a cluster on the sparse path (one state, coherence). */
void calCoherenceCceSparse(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, MatrixXcd* result, double start);

/** gCCE evolution = vector on the sparse path: the Hamiltonian of one cluster and the state exp(-i H tau_0) psi0 of the
 *  first interval, carried from one time point to the next. */
struct GcceSparseWork {
    CsrMatrix H;
    int qdim = 0, bdim = 0;
    MatrixXcd levels;                       // eigenvectors of the qubit Hamiltonian (columns)
    Eigen::VectorXcd first, psi, scratch, sum, part, piece;
    double firstTau = 0.0;
    bool started = false;
    KrylovStats stats;
    /** v <- exp(-i H tau) v, one expansion per qubit level (sparse_propagator.cpp). */
    void evolve(double tau, Eigen::VectorXcd& v, double tol);
};
/** |psi(tfree)> of the pulse sequence (propagateStateGcce on the sparse path); calls come with increasing tfree. The
 *  pulses Upulses (qdim x qdim) are applied block by block over the bath. */
const Eigen::VectorXcd& propagateStateGcceSparse(GcceSparseWork& w, const Eigen::VectorXcd& psi0, Pulse* pulse,
                                                 const MatrixXcd* Upulses, double tfree, double tol);

#endif
