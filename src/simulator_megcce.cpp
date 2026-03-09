#include <iomanip>
#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include "../include/superoperator.h"
#include "../include/memory.h"

/* ================================================================
 * ME-gCCE: Master Equation generalized Cluster Correlation Expansion
 *
 * Combines gCCE (full Hamiltonian, no secular approximation) with
 * Lindblad master equation for dissipative dynamics.
 *
 * Key differences from gCCE:
 *   - Uses superoperator formalism instead of unitary evolution
 *   - Handles dissipation via Lindblad jump operators
 *   - Propagator: exp(L*t) instead of exp(-iHt)
 *
 * Key differences from ME-CCE:
 *   - Uses full Hamiltonian (not projected/secular)
 *   - Works in qubit x bath product space
 *   - Supports multi-qubit, arbitrary pulses, density matrix output
 * ================================================================ */

/**
 * @brief Build incoherent superoperator for bath jump operators in qubit x bath space
 * @details Jump operators act on bath spins: C_full = I_qubit x C_bath
 */
static MatrixXcd incoherent_superoperator_megcce(BathArray* ba, MatrixXcd** bsigmas,
                                                  JumpOperatorArray* joa, int bdim, int qdim){
    int totdim = qdim * bdim;
    int N2 = totdim * totdim;
    MatrixXcd L_incoh = MatrixXcd::Zero(N2, N2);

    if (joa == NULL || joa->njump == 0){
        return L_incoh;
    }

    int nspin = BathArray_getNspin(ba);
    MatrixXcd Iq = MatrixXcd::Identity(qdim, qdim);

    for (int ij = 0; ij < joa->njump; ij++){
        const char* target_name = joa->ops[ij]->bath_name;
        const char* op_type = joa->ops[ij]->op_type;
        double rate = joa->ops[ij]->rate;

        for (int ib = 0; ib < nspin; ib++){
            char* spin_name = BathArray_getBath_i_name(ba, ib);

            if (strcasecmp(spin_name, target_name) == 0){
                int dim_i = BathArray_dimBath_i(ba, ib);

                // Build single-spin jump operator matrix
                MatrixXcd C_single = get_jump_operator_matrix(op_type, dim_i, rate);

                // Expand to full bath Hilbert space
                MatrixXcd C_bath = expand_operator(C_single, bsigmas, nspin, ib);

                // Expand to full qubit x bath space: I_q x C_bath
                MatrixXcd C_full = kron(Iq, C_bath);

                // Add the collapse superoperator
                L_incoh += collapse_superoperator_single(C_full);
            }
        }
    }

    return L_incoh;
}

/**
 * @brief Compute ME-gCCE super-propagator for a given free evolution time
 * @details Builds: S_total = S_free(tau_{N+1}) * S_pulse_N * ... * S_pulse_1 * S_free(tau_1)
 *          where S_free = exp(L * tau), S_pulse = op_to_supop(U, U†)
 */
static MatrixXcd calSuperPropagatorMegcce(QubitArray* qa, BathArray* ba, const MatrixXcd& L,
                                           Pulse* pulse, double tfree,
                                           MatrixXcd* Upulses, MatrixXcd* bUpulses){

    int nspin = BathArray_getNspin(ba);
    int qdim = QubitArray_dim(qa);
    int total_npulse = Pulse_getTotNpulse(pulse);
    double** total_sequence = Pulse_getTotSequence(pulse);
    int* total_sequence_indices = Pulse_getTotSequenceIndices(pulse);

    int N2 = L.rows(); // superoperator dimension (totdim^2)

    MatrixXcd Stotal;
    MatrixXcd* Sfrees = new MatrixXcd[total_npulse + 1];

    for (int ipulse = 0; ipulse < total_npulse + 1; ipulse++){

        double tau = tfree * total_sequence[ipulse][2];
        int sameTauIndex = total_sequence_indices[ipulse];

        // Build pulse superoperator
        MatrixXcd Spulse;
        if (ipulse < total_npulse){
            // Combine qubit and bath pulse operators
            MatrixXcd Upulse = Upulses[ipulse];
            MatrixXcd bUpulse;
            if (nspin > 0){
                bUpulse = bUpulses[ipulse];
            } else {
                bUpulse = MatrixXcd::Identity(1, 1);
            }
            MatrixXcd Upulse_full = kron(Upulse, bUpulse);

            // Convert to superoperator: S[rho] = U * rho * U†
            Spulse = op_to_supop(Upulse_full, Upulse_full.adjoint());
        } else {
            // Identity superoperator for last segment
            Spulse = MatrixXcd::Identity(N2, N2);
        }

        // Free evolution superoperator
        MatrixXcd Sfree;
        if (sameTauIndex == ipulse){
            // Compute: exp(L * tau)
            Sfree = (L * tau).exp();
        } else {
            Sfree = Sfrees[sameTauIndex];
        }

        // Compose total super-propagator
        if (ipulse == 0 && ipulse != total_npulse){
            Stotal = Spulse * Sfree;
        } else if (ipulse == 0 && ipulse == total_npulse){
            Stotal = Sfree;
        } else if (ipulse == total_npulse){
            Stotal = Sfree * Stotal;
        } else {
            Stotal = Spulse * Sfree * Stotal;
        }

        Sfrees[ipulse] = Sfree;
    }

    delete[] Sfrees;
    return Stotal;
}

/**
 * @brief ME-gCCE coherence/dm calculation
 */
MatrixXcd* calCoherenceMegcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, Output* op, JumpOperatorArray* joa){

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int ntotspin = nqubit + nspin;

    int qdim = QubitArray_dim(qa);
    int bdim = BathArray_dim(ba);
    int totdim = qdim * bdim;

    ////////////////////////////////
    // Hamiltonian (same as gCCE)
    ////////////////////////////////
    MatrixXcd Htot = MatrixXcd::Zero(totdim, totdim);

    // Pauli matrices of Qubit
    MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);

    // Qubit Hamiltonian
    MatrixXcd Hq = HamilQubit(qa, ba, qsigmas, cnf);
    MatrixXcd Hq_expand = kron(Hq, MatrixXcd::Identity(bdim, bdim));

    // Qubit's pulse operators (same as gCCE)
    MatrixXcd* Upulses = new MatrixXcd[pulse->total_npulse];
    if (pulse->detuning_factor == 1.0){
        calUpulses(Upulses, qa, pulse);
    } else {
        calTDUpulses(Upulses, qa, pulse, Hq);
    }

    MatrixXcd* bUpulses = new MatrixXcd[pulse->total_npulse];
    MatrixXcd** bsigmas = NULL;

    if (nspin > 0){
        // Bath's pulse operators
        calbUpulses(bUpulses, ba, pulse);

        // Bath Hamiltonian
        bsigmas = BathArray_PauliOperators(ba);
        MatrixXcd Hb = HamilBath(ba, bsigmas, cnf);
        MatrixXcd Hb_expand = kron(MatrixXcd::Identity(qdim, qdim), Hb);

        // Qubit-Bath Interaction Hamiltonian (full, not secular)
        MatrixXcd Hqb = HamilQubitBath(qa, ba, qsigmas, bsigmas, cnf);

        // Total Hamiltonian
        Htot = Hq_expand + Hb_expand + Hqb;
    } else {
        Htot = Hq;
    }

    ////////////////////////////////
    // Build Lindbladian superoperator
    ////////////////////////////////

    // Coherent part: L_coh = -i(H x I - I x H^T)
    MatrixXcd L_coh = coherent_superoperator(Htot);

    // Incoherent part: dissipators from jump operators on bath spins
    MatrixXcd L_incoh;
    if (nspin > 0 && joa != NULL && joa->njump > 0){
        L_incoh = incoherent_superoperator_megcce(ba, bsigmas, joa, bdim, qdim);
    } else {
        int N2 = totdim * totdim;
        L_incoh = MatrixXcd::Zero(N2, N2);
    }

    // Total Lindbladian
    MatrixXcd L = L_coh + L_incoh;

    ////////////////////////////////
    // Density matrix
    ////////////////////////////////
    int nstate = Config_getNstate(cnf);
    bool isEnsemble = (nstate == 0);

    MatrixXcd qrho0 = QubitArray_Rho0(qa);
    MatrixXcd brho0 = BathArray_Rho0(ba, isEnsemble);
    MatrixXcd rho0 = kron(qrho0, brho0);
    VectorXcd rho0_vec = mat_to_vec(rho0);

    ////////////////////////////////
    // Propagation
    ////////////////////////////////
    char* quantity = Config_getQuantity(cnf);
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    int nstep = Config_getNstep(cnf);
    double deltat = (double)Config_getDeltat(cnf);
    double tfree = 0.0;
    MatrixXcd* result = new MatrixXcd[nstep];

    MatrixXcd* result_tot = nullptr;
    if (strcasecmp(op->savemode, "allfull") == 0){
        result_tot = new MatrixXcd[nstep];
    }

    for (int i = 0; i < nstep; i++){

        // Compute super-propagator for this time step
        MatrixXcd Sprop = calSuperPropagatorMegcce(qa, ba, L, pulse, tfree, Upulses, bUpulses);

        // Evolve density matrix
        VectorXcd rhot_vec = Sprop * rho0_vec;
        MatrixXcd rhot = vec_to_mat(rhot_vec);

        // Save total density matrix if requested
        if (strcasecmp(op->savemode, "allfull") == 0){
            result_tot[i] = rhot;
        }

        // Partial trace over bath
        MatrixXcd reducedRhot = rhot;
        for (int ib = nspin - 1; ib >= 0; ib--){
            int bdim_i = BathArray_dimBath_i(ba, ib);
            reducedRhot = partialtrace(reducedRhot, bdim_i, bdim_i);
        }

        // Extract result
        if (strcasecmp(quantity, "coherence") == 0){
            result[i] = (psia.adjoint() * reducedRhot * psib);
        } else if (strcasecmp(quantity, "dm") == 0){
            result[i] = reducedRhot;
        } else {
            fprintf(stderr, "Error : calCoherenceMegcce : Quantity is neither coherence or dm\n");
            exit(1);
        }

        tfree += deltat;
    }

    // Save total density matrix
    if (strcasecmp(op->savemode, "allfull") == 0){
        int* cluster = allocInt1d(nspin);
        for (int i = 0; i < nspin; i++){
            cluster[i] = i;
        }
        Output_save_all(op, result_tot, nstep, deltat, cluster, nspin + nqubit, 0);
        freeInt1d(&cluster);
    }

    // Free
    if (bsigmas != NULL){
        for (int i = 0; i < nspin; i++){
            delete[] bsigmas[i];
        }
        delete[] bsigmas;
    }
    for (int i = 0; i < nqubit; i++){
        delete[] qsigmas[i];
    }
    delete[] qsigmas;
    delete[] Upulses;
    delete[] bUpulses;
    if (result_tot != nullptr){
        delete[] result_tot;
    }

    return result;
}
