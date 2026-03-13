#include "../include/hamiltonian.h"
#include "../include/simulator.h"
#include "../include/superoperator.h"
#include <iostream>

/**
 * @brief Qubit incoherent superoperator for ME-CCE (secular approximation)
 * @details For qubit jump operator C_q, the contribution to d(rho_tilde)/dt
 *          where rho_tilde = <alpha|rho|beta> is:
 *          [ <a|C_q|a><b|C_q|b>* - 1/2<a|C†C|a> - 1/2<b|C†C|b> ] * rho_tilde
 *          which is a scalar × I in superoperator space.
 *          Use bath_name = "center" in jump_operators to target the qubit.
 */
static MatrixXcd qubit_incoherent_superoperator_mecce(QubitArray* qa, JumpOperatorArray* joa, int N2){
    MatrixXcd L_qubit = MatrixXcd::Zero(N2, N2);
    if (joa == NULL || joa->njump == 0) return L_qubit;

    int qdim = QubitArray_dim(qa);
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    for (int ij = 0; ij < joa->njump; ij++){
        bool isQubit = (strcasecmp(joa->ops[ij]->bath_name, "qubit") == 0 ||
                        QubitArray_getQubitIdx_fromName(qa, joa->ops[ij]->bath_name) >= 0);
        if (!isQubit) continue;

        const char* op_type = joa->ops[ij]->op_type;
        double rate = joa->ops[ij]->rate;

        MatrixXcd C_q   = get_jump_operator_matrix(op_type, qdim, rate);
        MatrixXcd CdagC = C_q.adjoint() * C_q;

        doublec caa      = (psia.adjoint() * C_q   * psia)(0,0);
        doublec cbb_conj = std::conj((psib.adjoint() * C_q   * psib)(0,0));
        doublec cdcaa    = (psia.adjoint() * CdagC * psia)(0,0);
        doublec cdcbb    = (psib.adjoint() * CdagC * psib)(0,0);

        doublec coeff = caa * cbb_conj - 0.5 * cdcaa - 0.5 * cdcbb;
        L_qubit += coeff * MatrixXcd::Identity(N2, N2);
    }
    return L_qubit;
}

MatrixXcd *calCoherenceMecce(QubitArray *qa, BathArray *ba, Config *cnf,
                             Pulse *pls, JumpOperatorArray *joa) {

  int nspin = BathArray_getNspin(ba);
  int bdim = BathArray_dim(ba);

  int nqubit = QubitArray_getNqubit(qa);
  if (nqubit > 1) {
    fprintf(stderr,
            "Error : calCoherenceMecce : nqubit > 1, current nqubit = %d\n",
            nqubit);
    exit(1);
  }

  ////////////////////////////////
  // Hamiltonian
  ////////////////////////////////

  // Pauli matrices
  MatrixXcd **bsigmas = BathArray_PauliOperators(ba);

  // Bath Hamiltonian
  MatrixXcd Hb = HamilBath(ba, bsigmas, cnf);

  // Qubit - Bath Interaction Hamiltonian (projected)
  MatrixXcd *Hqb = HamilQubitBathSecularApp(qa, ba, bsigmas, cnf);

  // Total projected Hamiltonians
  MatrixXcd Halpha = Hb + Hqb[0];
  MatrixXcd Hbeta = Hb + Hqb[1];

  ////////////////////////////////
  // Build Lindbladian superoperator
  ////////////////////////////////

  // Coherent part: projected superoperator
  MatrixXcd L_coh = projected_coherent_superoperator(Halpha, Hbeta);

  // Incoherent part: bath spin dissipators
  MatrixXcd L_incoh = incoherent_superoperator(ba, bsigmas, joa, bdim);

  // Incoherent part: qubit (center) dissipators — scalar × I contribution
  MatrixXcd L_qubit = qubit_incoherent_superoperator_mecce(qa, joa, bdim * bdim);

  // Total Lindbladian
  MatrixXcd L = L_coh + L_incoh + L_qubit;

  // Flipped Lindbladian (swap alpha/beta) for CPMG
  // L_qubit is symmetric under alpha<->beta swap, so it stays the same
  MatrixXcd L_flipped =
      projected_coherent_superoperator(Hbeta, Halpha) + L_incoh + L_qubit;

  ////////////////////////////////
  // Initial state
  ////////////////////////////////
  int nstate = Config_getNstate(cnf);
  bool isEnsemble = (nstate == 0);

  MatrixXcd rho0;
  if (isEnsemble) {
    rho0 = MatrixXcd::Identity(bdim, bdim) / doublec(bdim, 0.0);
  } else {
    MatrixXcd psi0 = BathArray_Psi0(ba);
    rho0 = psi0 * psi0.adjoint();
  }
  VectorXcd rho0_vec = mat_to_vec(rho0);

  ////////////////////////////////
  // Propagation
  ////////////////////////////////

  int nstep = Config_getNstep(cnf);
  double deltat = (double)Config_getDeltat(cnf);
  MatrixXcd *result = new MatrixXcd[nstep];

  int npulse = Pulse_getNpulse(pls);
  double tfree = 0.0;

  for (int i = 0; i < nstep; i++) {

    MatrixXcd super_prop;

    if (npulse == 0) {
      // Ramsey: simple propagation
      super_prop = simple_incoherent_propagator(tfree, L);

    } else {
      // CPMG: use propagate_superpropagators
      double tau = tfree / (2.0 * npulse);

      MatrixXcd U_before = simple_incoherent_propagator(tau, L);
      MatrixXcd U_after = simple_incoherent_propagator(tau, L_flipped);

      super_prop = propagate_superpropagators(U_before, U_after, npulse);
    }

    // Evolve density matrix
    VectorXcd rho_t_vec = super_prop * rho0_vec;
    MatrixXcd rho_t = vec_to_mat(rho_t_vec);

    // Extract coherence: Tr(rho_t)
    // For ensemble: rho_0 = I/bdim already contains the 1/bdim factor
    // so Tr(rho_t) directly gives the coherence
    result[i] = MatrixXcd::Zero(1, 1);
    result[i](0, 0) = rho_t.trace();

    tfree += deltat;
  }

  // Free
  for (int i = 0; i < nspin; i++) {
    delete[] bsigmas[i];
  }
  delete[] bsigmas;
  delete[] Hqb;

  return result;
}
