#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include "../include/memory.h"

MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pulse, Output* op){

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
    int nstep     = Config_getNstep(cnf);
    double deltat = (double)Config_getDeltat(cnf);
    double tfree  = 0.0;
    MatrixXcd* result = new MatrixXcd[nstep];

    // Obtain full density matrix of a cluster for a certain purpose
    MatrixXcd* result_tot = nullptr; // Density matrix before trace out (total density matrix)
    
    if (strcasecmp(op->savemode,"allfull")==0){
        result_tot = new MatrixXcd[nstep]; //!
    }

    for (int i=0; i<nstep; i++){
        
        MatrixXcd Utot = calPropagatorGcce(qa, Htot, pulse, tfree);
        // Density matrix for time
        MatrixXcd rhot = Utot * rho0 * Utot.adjoint();

        // Save the density matrix
        if (strcasecmp(op->savemode,"allfull")==0){
            result_tot[i] = rhot;
        }
    
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

    // Save the density matrix
    if (strcasecmp(op->savemode,"allfull")==0){
        int* cluster = allocInt1d(nspin);
        for (int i = 0; i<nspin; i++){
            cluster[i] = i;
        }
        Output_save_all(op,result_tot,nstep,deltat,cluster,nspin+nqubit,0);
        freeInt1d(&cluster);
    }

    // free
    for (int i=0; i<nqubit; i++){
        delete[] qsigmas[i];
    }
    delete[] qsigmas;

    if (result_tot != nullptr){
        delete[] result_tot;
    }

    return result;
}

MatrixXcd getPauliOperator(MatrixXcd* sigmas, char axis) {
    switch (axis) {
        case 'I': return sigmas[0];
        case 'X': return sigmas[1];
        case 'Y': return sigmas[2];
        case 'Z': return sigmas[3];
        default:
            throw std::invalid_argument("Invalid axis character");
    }
}
double getAngle(double degree) {
    return degree * (M_PI / 180.0);  // Radian
}

MatrixXcd cal_Upulse_when_pulseiter_isFalse(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma){
    MatrixXcd exponent = ((-1.0) * doublec(0.0,1.0) * sigma * angle / 2.0); 
    MatrixXcd Upulse = exponent.exp();
    return Upulse;
}

MatrixXcd cal_Upulse_when_pulseiter_isTrue(QubitArray* qa, Pulse* pulse, double angle, int ipulse) {
    int nqubit = QubitArray_getNqubit(qa);
    MatrixXcd U_total = MatrixXcd::Identity(1, 1); // 초기 유니타리 (1x1 단위행렬)

    for (int iq = 0; iq < nqubit; iq++) {
        MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa, iq);
        MatrixXcd beta  = QubitArray_getQubit_i_beta(qa, iq);

        if (alpha.rows() == 0 || alpha.cols() == 0 || beta.rows() == 0 || beta.cols() == 0) {
            fprintf(stderr, "[Error] Function - cal_Upulse_when_pulseiter_isTrue: alpha or beta not set for qubit %d\n", iq);
            exit(EXIT_FAILURE);
        }

        MatrixXcd* sigmas = getGeneralPauliOperators(alpha, beta); // sigma[0]=I, [1]=X, [2]=Y, [3]=Z
        MatrixXcd  sigma = getPauliOperator(sigmas, (pulse->pulse_axes[ipulse]));

        // Expand (tensor product)
        MatrixXcd sigma_expanded = MatrixXcd::Identity(1, 1);
        for (int jq = 0; jq < nqubit; jq++) {
            int dim = QubitArray_dimQubit_i(qa, jq);
            MatrixXcd part;
            if (jq == iq) {
                part = sigma;
            } else {
                part = MatrixXcd::Identity(dim, dim);
            }
            sigma_expanded = kron(sigma_expanded, part);
        }

        MatrixXcd exponent = ((-1.0) * doublec(0.0,1.0) * angle * sigma_expanded / 2.0);
        MatrixXcd U_iq = exponent.exp();

        U_total = U_iq * U_total;

        delete[] sigmas;
    }

    return U_total;
}

MatrixXcd cal_Upulse_when_pulseiter_isTrue_whenOnlyXPulse(QubitArray* qa, Pulse* pulse, double angle){

    int nqubit = QubitArray_getNqubit(qa);
    int  qdim  = QubitArray_dim(qa);
    MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);

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
    return Upulse;
}

MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pulse, double tfree){

    double** sequence = Pulse_getSequence(pulse);
    int     npulse           = Pulse_getNpulse(pulse);
    double* pulse_angles     = Pulse_getPulseAngles(pulse);
    char*   pulse_axes       = Pulse_getPulseAxes(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);
    bool pulseiter           = Pulse_getPulseiter(pulse);
    int  qdim                = QubitArray_dim(qa);
    
    /////////////////////////////////////////////////////////////
    // Propagator for pulse
    //  Upulse = exp(-i*sigma_x*pi/2) -> Pi-Pulse
    MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);

    /////////////////////////////////////////////////////////////

    // Propagator for pulse : Expansion for bath spins : Upulse x II
    int bdim = Htot.rows() / Upulse.rows();
    //Upulse = kron(Upulse,MatrixXcd::Identity(bdim,bdim));

    // Propagator for total
    MatrixXcd Utotal;
    MatrixXcd* Ufrees = new MatrixXcd[npulse+1];

    // Propagator(total) =  U(tauN+1) (U_pulseN) ... (U_pulse1) U(tau1)
    // Pulse Index  0    1    2    3    4    5    6    7    8    |
    // Pulse (8#)   |____|____|____|____|____|____|____|____|____|
    //              |tau1|tau2|           .......           |tau9|
    //              |seq1|seq2|           .......           |seq9|
    // Pulse delay  0                                          tfree
    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        
        // Propagator for free evolution
        double       tau = tfree * sequence[ipulse][2];
        double     angle = getAngle(pulse_angles[ipulse]);  // Radian
        int sameTauIndex = sequence_indices[ipulse];

        // Calculate Pulse Operator
        if (!pulseiter){
            MatrixXcd* sigmas = QubitArray_PauliOperator_fromPsiaPsib(qa);
            MatrixXcd  sigma = getPauliOperator(sigmas, (pulse_axes[ipulse]));
            Upulse = cal_Upulse_when_pulseiter_isFalse(qa, pulse, angle, sigma);
            delete[] sigmas;
        } else{
            if ((strcasecmp(pulse->pulsename, "CPMG") == 0) || (strcasecmp(pulse->pulsename, "HahnEcho") == 0) || (strcasecmp(pulse->pulsename, "Ramsey") == 0)){
                Upulse = cal_Upulse_when_pulseiter_isTrue_whenOnlyXPulse(qa, pulse, angle);
            }else{
                Upulse = cal_Upulse_when_pulseiter_isTrue(qa, pulse, angle, ipulse);
            }
        }

        Upulse = kron(Upulse,MatrixXcd::Identity(bdim,bdim));
        // Pulse Operator Up = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t

        // Propagator for free evolution
        MatrixXcd Ufree;
        // Calculate free evolution operators (Ufree)
        if (sameTauIndex == ipulse){ // if Ufree haven't been calculated
            // U = exp(-iHtau)
            Ufree = ((-1.0) * doublec(0.0,1.0) * Htot * tau).exp();
        }
        else{ // if Ufree have been calculated
            // get previously calculated Ufree for the same tau
            Ufree = Ufrees[sameTauIndex];
        }

        // =================================== //
        //     Calculate U_total operator      //
        // =================================== //
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

//////////////////////////////// 
// Huijin's Code !! //
//////////////////////////////// 
//MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pulse, double tfree){
//
//    int npulse = Pulse_getNpulse(pulse);
//    double** sequence = Pulse_getSequence(pulse);
//    bool pulseiter = Pulse_getPulseiter(pulse);
//    int qdim = QubitArray_dim(qa);
//    
//    /////////////////////////////////////////////////////////////
//    // Propagator for pulse
//    MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);
//
//    if (!pulseiter){
//        MatrixXcd* sigmas_general = QubitArray_PauliOperator_fromPsiaPsib(qa);
//        // Note here, i'm only considering the following pulse situation:
//        //  Upulse = exp(-i*sigma_x*pi/2)
//        MatrixXcd UpulseExponent = ((-1.0) * doublec(0.0,1.0) * sigmas_general[1] * M_PI/2.0);
//        Upulse = UpulseExponent.exp();
//    }else{
//        int nqubit = QubitArray_getNqubit(qa);
//        for (int iq=0; iq<nqubit; iq++){
//
//            MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa,iq);
//            MatrixXcd beta = QubitArray_getQubit_i_beta(qa,iq);
//
//            if (alpha.rows() == 0 && alpha.cols() == 0){
//                fprintf(stderr,"Error : calPropagatorGcce : alpha or bata is not set\n");
//                fprintf(stderr,"If pulseiter turned on, you have to give alpha,beta\n");
//                exit(1);
//            }
//
//            /////////////////////////////////////////////////////////////
//            // Generate the general Pauli operators for i-th qubit
//            MatrixXcd sigmaExpanded_x;
//            MatrixXcd* sigmas = getGeneralPauliOperators(alpha,beta);
//            /////////////////////////////////////////////////////////////
//            MatrixXcd sigma;
//            for (int jq=0; jq<iq; jq++){
//                if (iq==jq){
//                    sigma = sigmas[1];
//                }else{
//                    int jqdim = QubitArray_dimQubit_i(qa,jq);
//                    sigma = MatrixXcd::Identity(jqdim,jqdim);
//                }
//
//                if (jq==0){
//                    sigmaExpanded_x = sigma;
//                }else{
//                    sigmaExpanded_x = kron(sigmaExpanded_x,sigma);
//                }
//            }
//            /////////////////////////////////////////////////////////////
//            // -i*sigma_x_i
//            MatrixXcd Upulse_iq = ((-1.0) * doublec(0.0,1.0) * sigmaExpanded_x);
//            Upulse = Upulse * Upulse_iq;
//        }
//    }
//    /////////////////////////////////////////////////////////////
//
//    // Propagator for pulse : Expansion for bath spins : Upulse x II
//    int bdim = Htot.rows() / Upulse.rows();
//    Upulse = kron(Upulse,MatrixXcd::Identity(bdim,bdim));
//
//    // Propagator for total
//    MatrixXcd Utotal;
//    MatrixXcd* Ufrees = new MatrixXcd[npulse+1];
//
//    // Propagator(total) =  U(tauN+1) (U_pulseN) ... (U_pulse1) U(tau1)
//    // Pulse Index  0  1  2  3  4  5  6  7  8 
//    // Pulse (8#)    __|__|__|__|__|__|__|__|__
//    // Pulse delay  0  tau                    tfree
//    for (int ipulse=0; ipulse<npulse+1; ipulse++){
//        
//        // Propagator for free evolution
//        MatrixXcd Ufree;
//
//        double tau = tfree * sequence[ipulse][2];
//        int sameTauIndex = (int)sequence[ipulse][3];
//
//        // Calculate operators (Ufree)
//        if (sameTauIndex == ipulse){ // if Ufree haven't been calculated
//            // U = exp(-iHtau)
//            Ufree = ((-1.0) * doublec(0.0,1.0) * Htot * tau).exp();
//        }
//        else{ // if Ufree have been calculated
//            // get previously calculated Ufree for the same tau
//            Ufree = Ufrees[sameTauIndex];
//        }
//
//        // Calculate operators (Upulse)
//        if (ipulse==0 && ipulse != npulse){
//            Utotal = Upulse * Ufree;
//        }
//        else if (ipulse==0 && ipulse == npulse){
//            Utotal = Ufree;
//        }
//        else if (ipulse==npulse){
//            Utotal = Ufree * Utotal;
//        }
//        else {
//            Utotal = Upulse * Ufree * Utotal;
//        }
//        Ufrees[ipulse] = Ufree;
//    }
//    
//    // free Ufrees
//    delete[] Ufrees;
//
//    return Utotal;
//}
