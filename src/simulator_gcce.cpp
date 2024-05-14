#include "../include/simulator.h"
#include "../include/hamiltonian.h"


MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls){

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int ntotspin = nqubit + nspin;

    int qdim = QubitArray_dim(qa);
    int bdim = BathArray_dim(ba);
    int totdim = qdim*bdim;

    ////////////////////////////////
    // Hamiltonian
    ////////////////////////////////

    // Declare the total Hamiltonian
    MatrixXcd Htot = MatrixXcd::Zero(totdim,totdim);

    // Pauli matrices of Qubit
    MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);

    // Qubit Hamiltonian
    MatrixXcd Hq = HamilQubit(qa,ba,qsigmas,cnf);
    MatrixXcd Hq_expand = kron(Hq,MatrixXcd::Identity(bdim,bdim));

    // Bath Hamiltonian
    if (nspin > 0){
        MatrixXcd** bsigmas = BathArray_PauliOperators(ba);
        MatrixXcd Hb = HamilBath(ba,bsigmas,cnf);
        MatrixXcd Hb_expand = kron(MatrixXcd::Identity(qdim,qdim),Hb);

        // Qubit - Bath Interaction Hamiltonian
        MatrixXcd Hqb = HamilQubitBath(qa,ba,qsigmas,bsigmas,cnf);

        // Total Hamiltonian
        Htot = Hq_expand + Hb_expand + Hqb;


        for (int i=0; i<nspin; i++){
            delete[] bsigmas[i];
        }
        delete[] bsigmas;

    }else{
        Htot = Hq;
    }

    ////////////////////////////////
    // Density matrix
    ////////////////////////////////

    // Check if this is the ensemble calculation
    int nstate = Config_getNstate(cnf);
    bool isEnsemble = false;
    if (nstate == 0){isEnsemble = true;}

    MatrixXcd qrho0 = QubitArray_Rho0(qa);
    MatrixXcd brho0 = BathArray_Rho0(ba, isEnsemble);
    MatrixXcd rho0 = kron(qrho0,brho0);

    ////////////////////////////////
    // Propagation
    ////////////////////////////////

    // Quantity
    char* quantity = Config_getQuantity(cnf);

    // alpha, beta
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    // Alloc the result variable
    int nstep = Config_getNstep(cnf);
    double deltat = (double)Config_getDeltat(cnf);
    MatrixXcd* result = new MatrixXcd[nstep];

    double tfree = 0.0;

    for (int i=0; i<nstep; i++){

        MatrixXcd Utot = calPropagatorGcce(qa, Htot, pls, tfree);

        // Density matrix for time
        MatrixXcd rhot = Utot * rho0 * Utot.adjoint();

        // Trace for tha bath state
        MatrixXcd reducedRhot = rhot;
        for (int ib=nspin-1; ib>=0; ib--){
            int bdim_i = BathArray_dimBath_i(ba, ib);
            reducedRhot = partialtrace(reducedRhot, bdim_i, bdim_i);
        }

        // Get the phase for qubit' two states
        if (strcasecmp(quantity,"coherence")==0){ // Coherence
            result[i] = (psia.adjoint() * reducedRhot * psib);

        }else if (strcasecmp(quantity,"dm")==0){ // Density matrix
            result[i] = reducedRhot;

        }else{
            fprintf(stderr,"Error : Quantity is neither coherence or dm\n");
            exit(1);
        }

        // Update tFree
        tfree += deltat;
    }

    // free
    for (int i=0; i<nqubit; i++){
        delete[] qsigmas[i];
    }
    delete[] qsigmas;

    return result;
}



MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pls, double tfree){

    int npulse = Pulse_getNpulse(pls);
    double** sequence = Pulse_getSequence(pls);
    bool pulseiter = Pulse_getPulseiter(pls);
    int qdim = QubitArray_dim(qa);
    
    /////////////////////////////////////////////////////////////
    // Propagator for pulse
    MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);

    if (!pulseiter){
        MatrixXcd* sigmas_general = QubitArray_PauliOperator_fromPsiaPsib(qa);
        // Note here, i'm only considering the following pulse situation:
        //  Upulse = exp(-i*sigma_x*pi/2)
        MatrixXcd UpulseExponent = ((-1.0) * doublec(0.0,1.0) * sigmas_general[1] * M_PI/2.0);
        Upulse = UpulseExponent.exp();
    }else{
        int nqubit = QubitArray_getNqubit(qa);
        for (int iq=0; iq<nqubit; iq++){

            MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa,iq);
            MatrixXcd beta = QubitArray_getQubit_i_beta(qa,iq);

            if (alpha.rows() == 0 && alpha.cols() == 0){
                fprintf(stderr,"Error : calPropagatorGcce : alpha or bata is not set\n");
                fprintf(stderr,"If pulseiter turned on, you have to give alpha,beta\n");
                exit(1);
            }

            /////////////////////////////////////////////////////////////
            // Generate the general Pauli operators for i-th qubit
            MatrixXcd sigmaExpanded_x;
            MatrixXcd* sigmas = getGeneralPauliOperators(alpha,beta);
            /////////////////////////////////////////////////////////////
            MatrixXcd sigma;
            for (int jq=0; jq<iq; jq++){
                if (iq==jq){
                    sigma = sigmas[1];
                }else{
                    int jqdim = QubitArray_dimQubit_i(qa,jq);
                    sigma = MatrixXcd::Identity(jqdim,jqdim);
                }

                if (jq==0){
                    sigmaExpanded_x = sigma;
                }else{
                    sigmaExpanded_x = kron(sigmaExpanded_x,sigma);
                }
            }
            /////////////////////////////////////////////////////////////
            // -i*sigma_x_i
            MatrixXcd Upulse_iq = ((-1.0) * doublec(0.0,1.0) * sigmaExpanded_x);
            Upulse = Upulse * Upulse_iq;
        }
    }
    /////////////////////////////////////////////////////////////

    // Propagator for pulse : Expansion for bath spins : Upulse x II
    int bdim = Htot.rows() / Upulse.rows();
    Upulse = kron(Upulse,MatrixXcd::Identity(bdim,bdim));

    // Propagator for total
    MatrixXcd Utotal;
    MatrixXcd* Ufrees = new MatrixXcd[npulse+1];

    // Propagator(total) =  U(tauN+1) (U_pulseN) ... (U_pulse1) U(tau1)
    // Pulse Index  0  1  2  3  4  5  6  7  8 
    // Pulse (8#)    __|__|__|__|__|__|__|__|__
    // Pulse delay  0  tau                    tfree
    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        
        // Propagator for free evolution
        MatrixXcd Ufree;

        double tau = tfree * sequence[ipulse][2];
        int sameTauIndex = (int)sequence[ipulse][3];

        // Calculate operators (Ufree)
        if (sameTauIndex == ipulse){ // if Ufree haven't been calculated
            // U = exp(-iHtau)
            Ufree = ((-1.0) * doublec(0.0,1.0) * Htot * tau).exp();
        }
        else{ // if Ufree have been calculated
            // get previously calculated Ufree for the same tau
            Ufree = Ufrees[sameTauIndex];
        }

        // Calculate operators (Upulse)
        if (ipulse==0 && ipulse != npulse){
            Utotal = Upulse * Ufree;
        }
        else if (ipulse==0 && ipulse == npulse){
            Utotal = Ufree;
        }
        else if (ipulse==npulse){
            Utotal = Ufree * Utotal;
        }
        else {
            Utotal = Upulse * Ufree * Utotal;
        }
        Ufrees[ipulse] = Ufree;
    }
    
    // free Ufrees
    delete[] Ufrees;

    return Utotal;
}