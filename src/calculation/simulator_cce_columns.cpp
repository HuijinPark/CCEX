/**
 * @file simulator_cce_columns.cpp
 * @brief The one-state cCCE time evolution with a block of time points as matrix columns (column_evolution).
 */
// A62 (2026-10-15; the idea of origin/cce-fast-mkl_hyeonsu cf05bfa3d): the one-state cCCE time evolution of cceEvolve with
// a block of time points as matrix columns. evolveSingleState (simulator_cce.cpp) carries a row vector and a column vector
// through the intervals of the sequence for one time point at a time - two matrix-vector products per interval, each reading
// M or M^dagger (16 D^2 bytes) again; here `columns` time points go together, so each interval is one matrix-matrix product
// (columns x D times D x D on the row side, D x D times D x columns on the column side) and M is read once per block.
// The same factors and the same number of multiplications; the summation order differs, so the last bits do too.
// A file of its own, so that simulator_cce.cpp compiles as before.
#include "../../include/calculation/simulator.h"
#include <algorithm>
#include <cmath>

namespace {
constexpr int phaseAnchor = 16;   // as simulator_cce.cpp (A16 b): recursive computes sin and cos at steps 0, 16, 32, ...

// C = A B: the library's zgemm in a build with MKL (A46), Eigen's product otherwise
inline void product(const MatrixXcd& A, const MatrixXcd& B, MatrixXcd& C){
#ifdef CCEX_BLAS_DIRECT
    ccexBlas::product(A,B,C);
#else
    C.noalias() = A * B;
#endif
}

template<bool Recursive>
void evolveColumns(const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown, double** sequence,
                   const int* sequenceIndices, int npulse, int nstep, double deltat, int columns,
                   MatrixXcd* result, CceWorkspace& ws){
    const int bdim = int(ws.M.rows());
    const MatrixXcd& M = ws.M;
    const MatrixXcd& Madj = ws.Madj;
    // the phases of the current step per distinct interval, laid out as evolveSingleState lays them out (up, down and with
    // Recursive their one-step factors), and computed by the same expressions in the same order of steps
    constexpr int slots = Recursive ? 4 : 2;
    std::vector<doublec>& buf = ws.phases;
    if (buf.size() < size_t(slots*bdim*(npulse+1))) buf.resize(size_t(slots*bdim*(npulse+1)));
    ws.blockUp.resize(size_t(npulse+1)); ws.blockDown.resize(size_t(npulse+1));
    double tfree = 0.0;
    for (int first = 0; first < nstep; first += columns){
        const int count = std::min(columns, nstep-first);
        // this block's phases: column j is step first + j
        for (int j = 0; j < count; j++){
            const int i = first + j;
            for (int ipulse = 0; ipulse < npulse+1; ipulse++){
                if (sequenceIndices[ipulse] != ipulse) continue;
                doublec* up = buf.data()+size_t(slots*bdim*ipulse);
                doublec* down = up+bdim;
                doublec* stepUp = Recursive ? up+2*bdim : nullptr;
                doublec* stepDown = Recursive ? up+3*bdim : nullptr;
                if (!Recursive || i % phaseAnchor == 0){
                    const double tau = tfree * sequence[ipulse][2];
                    for (int k=0; k<bdim; k++){
                        const double a_up = valUp(k).real()*tau, a_down = valDown(k).real()*tau;
                        up[k]   = doublec(std::cos(a_up),   -std::sin(a_up));
                        down[k] = doublec(std::cos(a_down), -std::sin(a_down));
                    }
                    if (Recursive && i == 0){
                        const double dtau = deltat * sequence[ipulse][2];
                        for (int k=0; k<bdim; k++){
                            const double a_up = valUp(k).real()*dtau, a_down = valDown(k).real()*dtau;
                            stepUp[k]   = doublec(std::cos(a_up),   -std::sin(a_up));
                            stepDown[k] = doublec(std::cos(a_down), -std::sin(a_down));
                        }
                    }
                } else {
                    for (int k=0; k<bdim; k++){ up[k] *= stepUp[k]; down[k] *= stepDown[k]; }
                }
                MatrixXcd& Pu = ws.blockUp[ipulse];
                MatrixXcd& Pd = ws.blockDown[ipulse];
                if (j == 0){ Pu.resize(bdim,count); Pd.resize(bdim,count); }
                Pu.col(j) = Eigen::Map<const Eigen::VectorXcd>(up,bdim);
                Pd.col(j) = Eigen::Map<const Eigen::VectorXcd>(down,bdim);
            }
            tfree += deltat;
        }
        // the chain of evolveSingleState, one row (of R, count x D) and one column (of C, D x count) per time point
        MatrixXcd& R = ws.blockRow;  MatrixXcd& Rnext = ws.blockRowNext;
        MatrixXcd& C = ws.blockCol;  MatrixXcd& Cnext = ws.blockColNext;
        for (int ipulse = 0; ipulse < npulse+1; ipulse++){
            const int same = sequenceIndices[ipulse];
            const MatrixXcd& Pu = ws.blockUp[same];
            const MatrixXcd& Pd = ws.blockDown[same];
            if (ipulse == 0){
                R = (Pd.adjoint().array().rowwise() * ws.rowStart.array()).matrix();   // psi^dagger M_down D_d^*
                Cnext = (Pu.array().colwise() * ws.colStart.array()).matrix();
                product(Madj,Cnext,C);                                                // M^dagger D_u M_up^dagger psi
            } else if (ipulse%2 == 1){
                product(R,Madj,Rnext); R = (Rnext.array() * Pu.adjoint().array()).matrix();   // * M^dagger D_u^*
                Cnext = (Pd.array() * C.array()).matrix();
                product(M,Cnext,C);                                                           // M D_d *
            } else {
                product(R,M,Rnext);    R = (Rnext.array() * Pd.adjoint().array()).matrix();   // * M D_d^*
                Cnext = (Pu.array() * C.array()).matrix();
                product(Madj,Cnext,C);                                                        // M^dagger D_u *
            }
        }
        for (int j = 0; j < count; j++){
            result[first+j].resize(1,1);
            result[first+j](0,0) = (R.row(j) * C.col(j))(0,0);
        }
    }
}
}

void cceEvolveColumns(const Eigen::VectorXcd& valUp, const Eigen::VectorXcd& valDown, double** sequence,
                      const int* sequenceIndices, int npulse, int nstep, double deltat, bool recursive, int columns,
                      MatrixXcd* result, CceWorkspace& ws){
    if (recursive) evolveColumns<true>(valUp,valDown,sequence,sequenceIndices,npulse,nstep,deltat,columns,result,ws);
    else           evolveColumns<false>(valUp,valDown,sequence,sequenceIndices,npulse,nstep,deltat,columns,result,ws);
}
