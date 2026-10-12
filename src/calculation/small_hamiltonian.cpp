#include "../../include/calculation/small_hamiltonian.h"
#include "../../include/calculation/hamiltonian.h"
#include <cstring>
#include <memory>

// A25 (improvements.md): the cluster Hamiltonians without the allocations and repeated work of the general functions.
// Every value is formed by the expression simulator_hamiltonian.cpp / bath.cpp / hamiltonian.cpp / qubit.cpp use for it,
// on the same Eigen types, in the same order. What changes:
//  - storage: buffers kept from one cluster to the next instead of new matrices;
//  - reuse of pure functions: the spin operators per dimension, the tensor of a pair (read up to three times per pair),
//    ms of the qubit states, the qubit Hamiltonian and the pulse propagators while their inputs are the same bits;
//  - terms that are all zeros are not added. A sum that starts at +0 is never -0 (x + y is -0 only when both are), so
//    adding +0 or -0 to it leaves every entry as it is. A term is all zeros when its coefficient is +0 or -0 and the
//    operator it multiplies is finite.
// The products the general code rounds in a function of another file before adding them (mean-field and Overhauser
// sums) are rounded the same way here: product() is not inlined, so the compiler cannot fuse them into the addition.
namespace {

// the spin operators (I, X, Y, Z) of one dimension: getPauliOperators, formed once
std::vector<std::unique_ptr<MatrixXcd[]>> pauliByDim;
MatrixXcd* pauli(int d){
    if ((int)pauliByDim.size() <= d) pauliByDim.resize(d+1);
    if (!pauliByDim[d]) pauliByDim[d].reset(getPauliOperators(d));
    return pauliByDim[d].get();
}
std::vector<MatrixXcd> identityByDim;
const MatrixXcd& identity(int d){
    if ((int)identityByDim.size() <= d) identityByDim.resize(d+1);
    if (identityByDim[d].rows() != d) identityByDim[d] = MatrixXcd::Identity(d,d);
    return identityByDim[d];
}

// buffers
std::vector<MatrixXcd*> bathSigmas, allSigmas, qubitSigmas;
MatrixXcd zeeman, detuning, quad, disorder, single, pair, hqb, small, expanded, hqExpand, hbExpand, bath;
std::vector<MatrixXcd> tensors;      // BathArray_int_i_j(i,j), i < j, at i*nspin + j

__attribute__((noinline)) double product(double a, float b){ return a * b; }

bool sameBits(const MatrixXcd& a, const MatrixXcd& b){
    return a.rows()==b.rows() && a.cols()==b.cols() &&
           std::memcmp(a.data(), b.data(), sizeof(doublec)*size_t(a.size())) == 0;
}
template<class T> bool sameBits(const T& a, const T& b){ return std::memcmp(&a, &b, sizeof(T)) == 0; }
inline bool zero(const doublec& c){ return c.real() == 0.0 && c.imag() == 0.0; }

// calHamiltonianSingleInt(Vector, P): H = 0; H += Vector(0,i) * P[i+1]
void singleInt(MatrixXcd& H, const doublec v[3], const MatrixXcd* P){
    H = MatrixXcd::Zero(P[0].rows(), P[0].cols());
    for (int i=0; i<3; i++) if (!zero(v[i])) H += v[i] * P[i+1];
}
// calHamiltonianSelfInt(P, T): H = 0; H += T(i,j) * P[i+1] * P[j+1]
void selfInt(MatrixXcd& H, const MatrixXcd* P, const MatrixXcd& T){
    H = MatrixXcd::Zero(P[0].rows(), P[0].cols());
    for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
            if (!zero(T(i,j))) H += T(i,j) * P[i+1] * P[j+1];
}

// addHamiltonianHeteroInt(H, P, T, nspin, i, j). Entry (xu, yv) of its two-spin operator is
//   sum_ab T_ab (S_i^a(x,y) S_j^b(u,v))   in the order ab = 00, 01, ..., 22,
// started from the first term. Here the products S^a S^b are formed once per pair of dimensions, the zero ones dropped,
// and the sum started from +0: the same nonzero terms in the same order, and a zero entry is added to a sum that is
// never -0. Each S^a S^b of spin operators is real or imaginary (X, Z real, Y imaginary), so in T_ab times it one of the
// two cross products is a zero and each part is a single rounded product, fused into the sum or not - the same value.
// Not used (the general function is) when a product has two nonzero parts or T is not finite (0 * inf).
struct PairTerm { int ab; doublec c; };
struct PairTable { int di = 0, dj = 0; bool usable = true; std::vector<std::vector<PairTerm>> entries; };
std::vector<PairTable> pairTables;
const PairTable& pairTable(const MatrixXcd* Si, const MatrixXcd* Sj){
    const int di = Si[0].rows(), dj = Sj[0].rows();
    for (const PairTable& t: pairTables) if (t.di == di && t.dj == dj) return t;
    PairTable t; t.di = di; t.dj = dj; t.entries.resize(size_t(di*dj)*size_t(di*dj));
    for (int x=0; x<di; x++) for (int y=0; y<di; y++) for (int u=0; u<dj; u++) for (int v=0; v<dj; v++)
        for (int a=0; a<3; a++) for (int b=0; b<3; b++){
            const doublec c = Si[a+1](x,y) * Sj[b+1](u,v);
            if (c.real() != 0.0 && c.imag() != 0.0) t.usable = false;
            if (!zero(c)) t.entries[size_t(x*dj+u)*(di*dj) + (y*dj+v)].push_back({a*3+b, c});
        }
    pairTables.push_back(std::move(t));
    return pairTables.back();
}
void addPair(MatrixXcd& H, MatrixXcd** P, const MatrixXcd& T, int nSpin, int iSpin, int jSpin){
    const PairTable& table = pairTable(P[iSpin], P[jSpin]);
    if (!table.usable || !T.allFinite()){ addHamiltonianHeteroInt(H, P, T, nSpin, iSpin, jSpin); return; }
    const int di = table.di, dj = table.dj, n = di*dj;
    small.resize(n, n);
    for (int row=0; row<n; row++)
        for (int col=0; col<n; col++){
            doublec value(0.0, 0.0);
            for (const PairTerm& term: table.entries[size_t(row)*n + col]) value += T(term.ab/3, term.ab%3) * term.c;
            small(row, col) = value;
        }
    long long dim = 1, si = 1, sj = 1;
    for (int k=0; k<nSpin; k++) dim *= P[k][0].rows();
    for (int k=iSpin+1; k<nSpin; k++) si *= P[k][0].rows();
    for (int k=jSpin+1; k<nSpin; k++) sj *= P[k][0].rows();
    for (long long r=0; r<dim; r++){
        const int x = int((r/si)%di), u = int((r/sj)%dj);
        const long long base = r - x*si - u*sj;
        for (int y=0; y<di; y++) for (int v=0; v<dj; v++)
            H(r, base + y*si + v*sj) += small(x*dj+u, y*dj+v);
    }
}

// expandHamiltonian(Pmats, Hi, nSpin, iSpin) into H
void expandInto(MatrixXcd& H, MatrixXcd** Pmats, const MatrixXcd& Hi, int nSpin, int iSpin){
    long long dim = 1, stride = 1;
    for (int k=0; k<nSpin; k++) dim *= Pmats[k][0].rows();
    for (int k=iSpin+1; k<nSpin; k++) stride *= Pmats[k][0].rows();
    const int d = Hi.rows();
    H = MatrixXcd::Zero(dim, dim);
    for (long long r=0; r<dim; r++){
        const int a = int((r/stride)%d);
        const long long base = r - a*stride;
        for (int b=0; b<d; b++) H(r, base + b*stride) = Hi(a,b);
    }
}
// kron(A, B) (utilities.cpp) into C
void kronInto(MatrixXcd& C, const MatrixXcd& A, const MatrixXcd& B){
    C.resize(A.rows()*B.rows(), A.cols()*B.cols());
    for (int i=0; i<A.rows(); i++)
        for (int j=0; j<A.cols(); j++)
            C.block(i*B.rows(),j*B.cols(),B.rows(),B.cols()) = A(i,j)*B;
}

MatrixXcd** spinOperators(BathArray* ba){
    const int nspin = BathArray_getNspin(ba);
    bathSigmas.resize(nspin);
    for (int ib=0; ib<nspin; ib++) bathSigmas[ib] = pauli(BathArray_dimBath_i(ba,ib));
    return bathSigmas.data();
}

// the pair tensors of the cluster: BathArray_int_i_j, which the bath Hamiltonian reads up to three times per pair
void pairTensors(BathArray* ba){
    const int nspin = BathArray_getNspin(ba);
    if ((int)tensors.size() < nspin*nspin) tensors.resize(nspin*nspin);
    for (int i=0; i<nspin; i++)
        for (int j=i+1; j<nspin; j++)
            tensors[i*nspin+j] = BathArray_int_i_j(ba,i,j);
}
const MatrixXcd& tensorOf(int nspin, int i, int j){ return tensors[i*nspin+j]; }

// BathArray_DisorderHamil(ba, sigmas, ib, rm_overlap = true) into H
void disorderInto(MatrixXcd& H, BathArray* ba, MatrixXcd** sigmas, int ib){
    const int nspin = BathArray_getNspin(ba);
    double value = BathArray_getBath_i_disorder(ba,ib);
    double overlap = 0.0;
    double xy[2] = {ba->bath[ib]->disorder_xy[0], ba->bath[ib]->disorder_xy[1]};
    const bool vector = ba->bath[ib]->disorder_vector;
    for (int jb=0; jb<nspin; jb++){
        // BathArray_getBath_i_disorder_j: J_zz m_j (0 for jb == ib)
        const float mj = BathArray_getBath_i_state(ba,jb);
        if (ib < jb)       overlap += product(tensorOf(nspin,ib,jb)(2,2).real(), mj);
        else if (ib == jb) overlap += 0.0;
        else               overlap += product(tensorOf(nspin,jb,ib)(2,2).real(), mj);
        if (vector && ba->bath[jb]->disorder_vector){
            // BathArray_getBath_i_disorderXY_j: J(0,2), J(1,2) times m_j, J of (ib, jb) (transposed when ib > jb); 0 for jb == ib
            double f[2] = {0.0, 0.0};
            if (ib < jb){
                const MatrixXcd& J = tensorOf(nspin,ib,jb);
                f[0] = product(J(0,2).real(), mj); f[1] = product(J(1,2).real(), mj);
            } else if (ib > jb){
                const MatrixXcd& J = tensorOf(nspin,jb,ib);
                f[0] = product(J(2,0).real(), mj); f[1] = product(J(2,1).real(), mj);
            }
            xy[0] -= f[0]; xy[1] -= f[1];
        }
    }
    value = value - overlap;
    const doublec v[3] = {doublec(xy[0]), doublec(xy[1]), doublec(value)};
    singleInt(H, v, sigmas[ib]);
}

// HamilBath(ba, sigmas, cnf) into Hb
void bathInto(MatrixXcd& Hb, BathArray* ba, MatrixXcd** sigmas, Config* cnf){
    const int nspin = BathArray_getNspin(ba);
    const int bdim = BathArray_dim(ba);
    float* bfield = Config_getBfield(cnf);
    const int nstate = Config_getNstate(cnf);
    Hb = MatrixXcd::Zero(bdim,bdim);
    pair = MatrixXcd::Zero(bdim,bdim);
    for (int ib=0; ib<nspin; ib++){
        const int d = BathArray_dimBath_i(ba,ib);
        const double gyro = BathArray_getBath_i_gyro(ba,ib);
        const doublec vz[3] = {doublec(-1.0 * gyro * bfield[0]), doublec(-1.0 * gyro * bfield[1]), doublec(-1.0 * gyro * bfield[2])};
        singleInt(zeeman, vz, sigmas[ib]);
        const doublec vd[3] = {doublec(0.0), doublec(0.0), doublec(BathArray_getBath_i_detuning(ba,ib))};
        singleInt(detuning, vd, sigmas[ib]);
        selfInt(quad, sigmas[ib], ba->bath[ib]->quad);
        if (nstate > 0) disorderInto(disorder, ba, sigmas, ib);
        else            disorder = MatrixXcd::Zero(d,d);
        single = zeeman + detuning + quad + disorder;
        addExpandedHamiltonian(Hb, sigmas, single, nspin, ib);
    }
    for (int i=0; i<nspin; i++)
        for (int j=i+1; j<nspin; j++){
            const MatrixXcd& tensor = tensorOf(nspin,i,j);
            if (!tensor.isZero(FLT_EPSILON)) addPair(pair, sigmas, tensor, nspin, i, j);
        }
    Hb += pair;
}

// ms of psia, psib (findZbasisSubLevel), kept while the two states are the same bits
MatrixXcd psiaSeen, psibSeen;
float msaSeen = 0.0f, msbSeen = 0.0f;
void qubitLevels(QubitArray* qa, float& msa, float& msb){
    if (!sameBits(qa->psia, psiaSeen) || !sameBits(qa->psib, psibSeen)){
        msaSeen = findZbasisSubLevel(qa->psia);
        msbSeen = findZbasisSubLevel(qa->psib);
        psiaSeen = qa->psia; psibSeen = qa->psib;
    }
    msa = msaSeen; msb = msbSeen;
}

// QubitArray_TotalHamil for one qubit, kept while what it reads (gyro, detuning, Overhauser field, zfs tensor, field,
// dimension) is the same bits; the Overhauser field of the qubit changes with the bath state
struct QubitKey { double gyro, detuning, overhaus; float field[3]; int dim; };
bool totalKept = false;
QubitKey totalKey;
MatrixXcd totalZfs, totalHq;
const MatrixXcd& qubitTotal(QubitArray* qa, MatrixXcd** qsigmas, float* bfield){
    QubitKey key;
    std::memset(&key, 0, sizeof(key));
    key.gyro = QubitArray_getQubit_i_gyro(qa,0);
    key.detuning = QubitArray_getQubit_i_detuning(qa,0);
    key.overhaus = QubitArray_getQubit_i_overhaus(qa,0);
    for (int k=0; k<3; k++) key.field[k] = bfield[k];
    key.dim = QubitArray_dimQubit_i(qa,0);
    if (!totalKept || !sameBits(key, totalKey) || !sameBits(qa->intmap[0][0], totalZfs)){
        totalHq = QubitArray_TotalHamil(qa, qsigmas, bfield);
        totalKey = key; totalZfs = qa->intmap[0][0]; totalKept = true;
    }
    return totalHq;
}
}

void cceClusterHamiltonians(QubitArray* qa, BathArray* ba, Config* cnf, MatrixXcd& Halpha, MatrixXcd& Hbeta){
    MatrixXcd** sigmas = spinOperators(ba);
    pairTensors(ba);
    bathInto(bath, ba, sigmas, cnf);
    if (Config_getKnight(cnf) || Config_getHfmedi(cnf)){
        // the general function for these two terms (not in the production inputs)
        std::unique_ptr<MatrixXcd[]> interaction(HamilQubitBathSecularApp(qa,ba,sigmas,cnf));
        Halpha = bath + interaction[0];
        Hbeta  = bath + interaction[1];
        return;
    }
    // HamilQubitBathSecularApp without knight and hfmedi: ms * (sum of Azx X + Azy Y + Azz Z on each spin) + 0
    if (QubitArray_getNqubit(qa) > 1){
        fprintf(stderr,"Error : nqubit > 1, current nqubit = %d\n",QubitArray_getNqubit(qa));
        ccexFatal();
    }
    const int nspin = BathArray_getNspin(ba);
    const int bdim = BathArray_dim(ba);
    hqb = MatrixXcd::Zero(bdim,bdim);
    for (int ib=0; ib<nspin; ib++){
        const MatrixXcd& tensor = ba->bath[ib]->hypf[0];
        const doublec v[3] = {doublec(tensor(2,0).real()), doublec(tensor(2,1).real()), doublec(tensor(2,2).real())};
        singleInt(small, v, sigmas[ib]);
        addExpandedHamiltonian(hqb, sigmas, small, nspin, ib);
    }
    float msa, msb;
    qubitLevels(qa, msa, msb);
    const double ms[2] = {msa, msb};
    // cceHamiltonians: interaction[k] = ((ms[k] * Hqb).array() + 0).matrix(), then bath + interaction[k] in place
    Halpha = ((ms[0] * hqb).array() + doublec(0.0,0.0)).matrix();
    Hbeta  = ((ms[1] * hqb).array() + doublec(0.0,0.0)).matrix();
    Halpha = bath + Halpha;
    Hbeta  = bath + Hbeta;
}

void gcceClusterHamiltonian(QubitArray* qa, BathArray* ba, Config* cnf, MatrixXcd& Htot, MatrixXcd& Hq){
    const int nqubit = QubitArray_getNqubit(qa);
    const int nspin = BathArray_getNspin(ba);
    const int qdim = QubitArray_dim(qa);
    const int bdim = BathArray_dim(ba);
    const int totdim = qdim*bdim;
    qubitSigmas.resize(nqubit);
    for (int iq=0; iq<nqubit; iq++) qubitSigmas[iq] = pauli(QubitArray_dimQubit_i(qa,iq));
    MatrixXcd** qsigmas = qubitSigmas.data();

    // HamilQubit: the qubit Hamiltonian less the Overhauser fields the cluster's spins put in it
    if (nqubit == 1) Hq = qubitTotal(qa, qsigmas, Config_getBfield(cnf));
    else             Hq = QubitArray_TotalHamil(qa, qsigmas, Config_getBfield(cnf));
    for (int iq=0; iq<nqubit; iq++){
        double overhaus = 0.0;
        for (int ib=0; ib<nspin; ib++)   // BathArray_getBath_i_overhaus_j: A_zz m
            overhaus += product(ba->bath[ib]->hypf[iq](2,2).real(), BathArray_getBath_i_state(ba,ib));
        const doublec v[3] = {doublec(0.0), doublec(0.0), doublec(overhaus)};
        singleInt(small, v, qsigmas[iq]);
        expandInto(expanded, qsigmas, small, nqubit, iq);
        Hq -= expanded;
    }
    kronInto(hqExpand, Hq, identity(bdim));

    // HamilBath
    MatrixXcd** bsigmas = spinOperators(ba);
    pairTensors(ba);
    bathInto(bath, ba, bsigmas, cnf);
    kronInto(hbExpand, identity(qdim), bath);

    // HamilQubitBath: each qubit-spin term added where it has entries (Hqb starts at +0 and is never -0, so the +0 the
    // full matrix of calHamiltonianHeteroInt added elsewhere changed nothing)
    allSigmas.resize(nqubit+nspin);
    for (int iq=0; iq<nqubit; iq++) allSigmas[iq] = qsigmas[iq];
    for (int ib=0; ib<nspin; ib++) allSigmas[nqubit+ib] = bsigmas[ib];
    hqb = MatrixXcd::Zero(totdim,totdim);
    for (int iq=0; iq<nqubit; iq++)
        for (int ib=0; ib<nspin; ib++)
            addPair(hqb, allSigmas.data(), ba->bath[ib]->hypf[iq], nqubit+nspin, iq, ib+nqubit);

    Htot = hqExpand + hbExpand + hqb;
}

// calUpulses for pulseiter false: exp(-i sigma angle / 2) per pulse, sigma from psia and psib, kept while psia, psib, the
// angles and the axes are the same bits
namespace {
MatrixXcd pulsePsia, pulsePsib;
std::vector<double> pulseAngles;
std::vector<char> pulseAxes;
std::vector<MatrixXcd> pulseKept;
bool pulsesKept = false;
}
void cachedUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse){
    const int npulse = Pulse_getNpulse(pulse);
    const double* angles = Pulse_getPulseAngles(pulse);
    const char* axes = Pulse_getPulseAxes(pulse);
    const bool same = pulsesKept && (int)pulseKept.size() == npulse && sameBits(qa->psia, pulsePsia) &&
                      sameBits(qa->psib, pulsePsib) &&
                      (npulse == 0 || (std::memcmp(angles, pulseAngles.data(), sizeof(double)*npulse) == 0 &&
                                       std::memcmp(axes, pulseAxes.data(), size_t(npulse)) == 0));
    if (!same){
        pulseKept.assign(npulse, MatrixXcd());
        calUpulses(pulseKept.data(), qa, pulse);
        pulsePsia = qa->psia; pulsePsib = qa->psib;
        pulseAngles.assign(angles, angles+npulse);
        pulseAxes.assign(axes, axes+npulse);
        pulsesKept = true;
    }
    for (int i=0; i<npulse; i++) Upulses[i] = pulseKept[i];
}

bool gcceHamiltonianFinite(const MatrixXcd& H){ return H.allFinite(); }

// validateGcceHamiltonian's asym = max|H - H^dag|, scale = max|H| decide; the squared moduli first: when they put asym
// far below 1e-12 scale the answer is known without the square roots, otherwise (and at the ends of the double range)
// the two values are computed as before.
bool gcceHamiltonianHermitian(const MatrixXcd& H, double& asym, double& scale){
    const Eigen::Index n = H.rows();
    double asym2 = 0.0, scale2 = 0.0;
    for (Eigen::Index j=0; j<n; j++)
        for (Eigen::Index i=0; i<n; i++){
            scale2 = std::max(scale2, std::norm(H(i,j)));
            asym2 = std::max(asym2, std::norm(H(i,j) - std::conj(H(j,i))));
        }
    if (scale2 > 1e-200 && scale2 < 1e200 && asym2 <= 0.25e-24 * scale2) return true;
    asym  = (H - H.adjoint()).cwiseAbs().maxCoeff();
    scale = H.cwiseAbs().maxCoeff();
    return !(scale > 0.0 && asym > 1e-12*scale);
}
