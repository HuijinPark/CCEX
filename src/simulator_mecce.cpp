#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include "../include/superoperator.h"
#include <iostream>

MatrixXcd* calCoherenceMecce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, JumpOperatorArray* joa){

    int nspin = BathArray_getNspin(ba);
    int bdim = BathArray_dim(ba);

    int nqubit = QubitArray_getNqubit(qa);
    if (nqubit > 1){
        fprintf(stderr, "Error : calCoherenceMecce : nqubit > 1, current nqubit = %d\n", nqubit);
        exit(1);
    }

    ////////////////////////////////
    // Hamiltonian
    ////////////////////////////////

    // Pauli matrices
    MatrixXcd** bsigmas = BathArray_PauliOperators(ba);

    // Bath Hamiltonian
    MatrixXcd Hb = HamilBath(ba, bsigmas, cnf);

    // Qubit - Bath Interaction Hamiltonian (projected)
    MatrixXcd* Hqb = HamilQubitBathSecularApp(qa, ba, bsigmas, cnf);

    // Total projected Hamiltonians
    MatrixXcd Halpha = Hb + Hqb[0];
    MatrixXcd Hbeta  = Hb + Hqb[1];

    ////////////////////////////////
    // Build Lindbladian superoperator
    ////////////////////////////////

    // Coherent part: projected superoperator
    MatrixXcd L_coh = projected_coherent_superoperator(Halpha, Hbeta);

    // Incoherent part: dissipators from jump operators
    MatrixXcd L_incoh = incoherent_superoperator(ba, bsigmas, joa, bdim);

    // Total Lindbladian
    MatrixXcd L = L_coh + L_incoh;

    // Flipped Lindbladian (swap alpha/beta) for CPMG
    MatrixXcd L_flipped = projected_coherent_superoperator(Hbeta, Halpha) + L_incoh;

    ////////////////////////////////
    // Initial state
    ////////////////////////////////
    int nstate = Config_getNstate(cnf);
    bool isEnsemble = (nstate == 0);

    MatrixXcd rho0;
    if (isEnsemble){
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
    MatrixXcd* result = new MatrixXcd[nstep];

    int npulse = Pulse_getNpulse(pls);
    double tfree = 0.0;

    for (int i = 0; i < nstep; i++){

        MatrixXcd super_prop;

        if (npulse == 0){
            // Ramsey: simple propagation
            super_prop = simple_incoherent_propagator(tfree, L);

        } else {
            // CPMG: use propagate_superpropagators
            double tau = tfree / (2.0 * npulse);

            MatrixXcd U_before = simple_incoherent_propagator(tau, L);
            MatrixXcd U_after  = simple_incoherent_propagator(tau, L_flipped);

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
    for (int i = 0; i < nspin; i++){
        delete[] bsigmas[i];
    }
    delete[] bsigmas;
    delete[] Hqb;

    return result;
}
