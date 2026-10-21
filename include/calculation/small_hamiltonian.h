/**
 * @file small_hamiltonian.h
 * @brief The cluster Hamiltonians of the CPU cCCE and gCCE paths without per-call allocations.
 */
#ifndef __CCEX_SMALL_HAMILTONIAN_H_
#define __CCEX_SMALL_HAMILTONIAN_H_
#include "calculation/simulator.h"

/**
 * The cluster Hamiltonians of the CPU cCCE and gCCE calculations, built with the arithmetic of HamilBath,
 * HamilQubitBathSecularApp, HamilQubit and HamilQubitBath (simulator_hamiltonian.cpp) - the same Eigen expressions on the
 * same types in the same order - without their per-call allocations: the spin operators are kept per dimension, the
 * pieces go into buffers kept from one cluster to the next, the dipolar tensor of a pair is formed once instead of
 * three times. The results are the same bits (small_hamiltonian.cpp says why at each step).
 */
/** Halpha, Hbeta of calCoherenceCceInto (bath + ms_a,b * secular hyperfine); written over. */
void cceClusterHamiltonians(QubitArray* qa, BathArray* ba, Config* cnf, MatrixXcd& Halpha, MatrixXcd& Hbeta);
/** Htot (qubit (x) 1 + 1 (x) bath + qubit-bath) and Hq of calCoherenceGcceInto, nspin > 0; written over. */
void gcceClusterHamiltonian(QubitArray* qa, BathArray* ba, Config* cnf, MatrixXcd& Htot, MatrixXcd& Hq);
/** calUpulses for pulseiter false, the matrices formed again only when psia, psib, the angles or the axes change. */
void cachedUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse);
/** The check of validateGcceHamiltonian (finite, Hermitian to 1e-12 of the largest entry) without its temporaries. */
bool gcceHamiltonianFinite(const MatrixXcd& H);
bool gcceHamiltonianHermitian(const MatrixXcd& H, double& asym, double& scale);

#endif
