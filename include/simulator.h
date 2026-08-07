#ifndef __CCEX_SIMULATOR_H_
#define __CCEX_SIMULATOR_H_

#include <unsupported/Eigen/MatrixFunctions>
#include "utilities.h"
#include "qubit.h"
#include "bath.h"
#include "defect.h"
#include "output.h"
#include "pulse.h"
#include "cluster.h"
#include "general.h"

//MatrixXcd cal_Upulse_withDetuning(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma, MatrixXcd sigma_z);
void calUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse);
void calTDUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse, MatrixXcd Hq);
void calculate(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cc, Pulse* pulse, Cluster* cluster, Output* output, int*** localClusters);

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, DefectArray* dfa, int nqubit);

MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Output* op);
MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls);

/**
 * Scratch for the gCCE propagator, owned by the caller and reused across the nstep
 * loop. Htot does not depend on the time step -- calCoherenceGcce builds it once and
 * the loop only varies tau -- so one eigendecomposition covers every tau:
 *     exp(-i*Htot*tau) = M * diag(exp(-i*lambda_k*tau)) * M^dag
 * That replaces a matrix exponential per step with a handful of scalar exponentials.
 * This is what simulator_cce.cpp has always done; gCCE was calling Eigen's general
 * .exp(), which runs a complex Schur decomposition and does not know Htot is Hermitian.
 */
struct GcceWork {
    bool useEigen;                       // false => legacy per-step matrix exponential
    MatrixXcd M, Madj;                   // eigenvectors of Htot, and their adjoint
    Eigen::VectorXd  evals;              // eigenvalues (real -- Htot is Hermitian)
    Eigen::VectorXd  evals_tau;          // scratch
    Eigen::VectorXcd phase;              // scratch
    MatrixXcd scratch;                   // scratch
    std::vector<MatrixXcd> Ufrees;       // per-pulse free-evolution propagators

    // evolution=vector only (see Config_setEvolution): |psi(t)>, the eigenbasis
    // intermediate of one Ufree application, the phase factor of each distinct tau, and
    // the pulse propagators already expanded over the bath.
    Eigen::VectorXcd psi, psiscratch;
    std::vector<Eigen::VectorXcd> phases;
    std::vector<MatrixXcd> UpulsesExpanded;
};

/** Builds everything in GcceWork that is fixed for the cluster. Call once per cluster. */
void prepareGcceWork(GcceWork* w, const MatrixXcd& Htot, bool useEigen);
MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pls, double tfree, MatrixXcd* Upulses, GcceWork* w);

/** evolution=vector: applies the same pulse sequence as calPropagatorGcce, but to the
 *  state vector psi0 instead of assembling the propagator. Result in w->psi.
 *  prepareGcceVectorWork must have been called for this cluster first. */
void prepareGcceVectorWork(GcceWork* w, MatrixXcd* Upulses, int npulse, int bdim);
const Eigen::VectorXcd& propagateStateGcce(const Eigen::VectorXcd& psi0, Pulse* pls, double tfree, GcceWork* w);

MatrixXcd HamilQubit(QubitArray* qa, BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilBath(BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilQubitBath(QubitArray* qa, BathArray* ba, MatrixXcd** qsigmas, MatrixXcd** bsigmas, Config* cnf);
MatrixXcd* HamilQubitBathSecularApp(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas, Config* cnf);
double calQubit_Energy_difference(QubitArray* qa, MatrixXcd Hq);
MatrixXcd cal_TDUpulse(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, MatrixXcd ox, MatrixXcd oy, char axis, double pulse_angle, double Ediff, double detuning_factor, double pulse_time, int Nstep);


MatrixXcd** gatherSigmas(MatrixXcd** qsigmas, MatrixXcd** bsigmas, int nqubit, int nspin);
MatrixXcd** gatherSigmas_singlequbit(MatrixXcd* qsigma, MatrixXcd** bsigmas, int nqubit, int nspin);
int close_state_index(MatrixXcd state, MatrixXcd eigenVectors);

#endif // __CCEX_SIMULATOR_H_
