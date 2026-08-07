#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include "../include/memory.h"

#define sPI 3.14159265358979323846

// exp(-i*x) for a real vector x. simulator_cce.cpp carries the same helper; both are
// kept file-local rather than shared through a header so that adding this one cannot
// perturb the CCE path's codegen, which is verified byte-for-byte against the reference.
static inline void phaseVectorGcce(const Eigen::VectorXd& x, Eigen::VectorXcd& out){
    const int n = (int)x.size();
    out.resize(n);
    for (int k=0; k<n; ++k){
        out(k) = doublec(std::cos(x(k)), -std::sin(x(k)));
    }
}

// exp(-i*Htot*tau) applied to w->psi, in place, without ever forming the propagator:
//     exp(-i*H*tau)|v> = M ( exp(-i*lambda*tau) .* (M^dag |v>) )
// Two matrix-vector products and one elementwise scaling, reading the decomposition
// prepareGcceWork already computed. `phase` is exp(-i*lambda*tau) for this tau.
// evolution=vector only.
static inline void applyUfreeGcce(GcceWork* w, const Eigen::VectorXcd& phase){
    w->psiscratch.noalias() = w->Madj * w->psi;
    w->psiscratch.array()  *= phase.array();
    w->psi.noalias()        = w->M * w->psiscratch;
}


void prepareGcceWork(GcceWork* w, const MatrixXcd& Htot, bool useEigen){

    w->useEigen = useEigen;

    // propagator=expm keeps the original per-step matrix exponential, so there is no
    // decomposition to prepare and no Hermiticity requirement to enforce.
    if (!useEigen) return;

    // The eigendecomposition below stands in for exp(-i*Htot*tau) only if Htot really
    // is Hermitian. It is a Hamiltonian, so it must be -- check rather than assume,
    // because SelfAdjointEigenSolver reads one triangle and would silently accept a
    // non-Hermitian matrix by ignoring the other half.
    double asym  = (Htot - Htot.adjoint()).cwiseAbs().maxCoeff();
    double scale = Htot.cwiseAbs().maxCoeff();
    if (scale > 0.0 && asym > 1e-12*scale){
        fprintf(stderr,"Error : prepareGcceWork : Htot is not Hermitian "
                       "(max|H-H^dag| = %.3e, max|H| = %.3e)\n", asym, scale);
        exit(1);
    }

    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(Htot);
    if (es.info() != Eigen::Success){
        fprintf(stderr,"Error : prepareGcceWork : eigendecomposition of Htot failed\n");
        exit(1);
    }
    w->evals = es.eigenvalues();
    w->M     = es.eigenvectors();
    w->Madj  = w->M.adjoint();
}

Matrix3cd Sx_spin1() {
    Matrix3cd Sx = Matrix3cd::Zero();
    double factor = 1.0 / sqrt(2.0);
    //double factor = 1.0;

    Sx(0,1) = factor;
    Sx(1,0) = factor;
    Sx(1,2) = factor;
    Sx(2,1) = factor;

    return Sx;
}

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
    //std::cout << "Hq = \n" << Hq << std::endl;
    //exit(1);
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

    // evolution=vector propagates |psi> instead of rho, which needs rho0 to BE a pure
    // state. QubitArray_Rho0 always is (psi0*psi0^dag) and BathArray_Rho0 is too except
    // for the ensemble average, so the requirement is nstate > 0. Refuse rather than
    // silently falling back -- a run that quietly used a different algorithm than the
    // one it was asked for is worse than a run that stops.
    const bool useVector = (strcasecmp(Config_getEvolution(cnf),"vector")==0);
    if (useVector){
        if (isEnsemble){
            fprintf(stderr,"Error : evolution=vector needs a pure rho0, but nstate=0 "
                           "(ensemble) makes the bath density matrix mixed. "
                           "Use nstate > 0, or evolution=matrix.\n");
            exit(1);
        }
        if (strcasecmp(Config_getPropagator(cnf),"eigen")!=0){
            fprintf(stderr,"Error : evolution=vector requires propagator=eigen "
                           "(with expm the per-step matrix exponential dominates and "
                           "there is nothing to gain).\n");
            exit(1);
        }
    }

    MatrixXcd rho0;
    Eigen::VectorXcd psi0;
    if (useVector){
        // BathArray_Psi0 leaves its result empty when nspin==0; the bath factor is the
        // 1x1 identity there, so the qubit state is already the whole state.
        psi0 = (nspin > 0) ? Eigen::VectorXcd(kron(QubitArray_getPsi0(qa),BathArray_Psi0(ba)))
                           : Eigen::VectorXcd(QubitArray_getPsi0(qa));
    }else{
        MatrixXcd qrho0 = QubitArray_Rho0(qa);
        MatrixXcd brho0 = BathArray_Rho0(ba, isEnsemble);
        rho0 = kron(qrho0,brho0);
    }

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

    MatrixXcd* result_tot = nullptr; // Density matrix before trace out (total density matrix)
    if (strcasecmp(op->savemode,"allfull")==0){
        result_tot = new MatrixXcd[nstep]; //!
    }

    MatrixXcd* Upulses = new MatrixXcd[pulse->npulse];
    if (pulse->detuning_factor == 1.0){
        calUpulses(Upulses, qa, pulse);
        //calTDUpulses(Upulses, qa, pulse, Hq);
    } else{
        calTDUpulses(Upulses, qa, pulse, Hq);
    }

    // Everything fixed for this cluster: the eigendecomposition that replaces a
    // matrix exponential on every step.
    GcceWork gw;
    prepareGcceWork(&gw, Htot, strcasecmp(Config_getPropagator(cnf),"eigen")==0);
    if (useVector){ prepareGcceVectorWork(&gw, Upulses, pulse->npulse, bdim); }

    MatrixXcd reducedRhot;
    for (int i=0; i<nstep; i++){

        if (useVector){
            // |psi(t)> = Utot |psi0>, without ever forming Utot.
            const Eigen::VectorXcd& psit = propagateStateGcce(psi0, pulse, tfree, &gw);

            // Tracing the bath out of |psi><psi| is a reshape. kron() lays the state out
            // as psi(iq*bdim + ib), so reading it as a qdim x bdim matrix Psi gives
            //   (Tr_bath |psi><psi|)_{ij} = sum_b psi(i,b) conj(psi(j,b)) = (Psi Psi^dag)_{ij}
            // which is exactly what the partialtrace chain below computes term by term.
            Eigen::Map<const Eigen::Matrix<doublec,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >
                Psi(psit.data(), qdim, bdim);
            reducedRhot.noalias() = Psi * Psi.adjoint();

            // allfull wants the full qubit+bath density matrix; from a pure state that is
            // an outer product, O(n^2), cheaper than the O(n^3) the matrix path pays.
            if (strcasecmp(op->savemode,"allfull")==0){
                result_tot[i] = psit * psit.adjoint();
            }
        }else{
            MatrixXcd Utot = calPropagatorGcce(qa, Htot, pulse, tfree, Upulses, &gw);
            // Density matrix for time
            MatrixXcd rhot = Utot * rho0 * Utot.adjoint();

            // Save the density matrix
            if (strcasecmp(op->savemode,"allfull")==0){
                result_tot[i] = rhot;
            }

            // Trace for tha bath state
            reducedRhot = rhot;
            for (int ib=nspin-1; ib>=0; ib--){
                int bdim_i = BathArray_dimBath_i(ba, ib);
                reducedRhot = partialtrace(reducedRhot, bdim_i, bdim_i);
            }
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
    delete[] Upulses;
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

double getPhase(char axis) {
    double PI = 3.14159265358979323846;
    switch (axis) {
        case 'X': return 0.0;
        case 'Y': return PI/2;
        default:
            throw std::invalid_argument("Invalid axis character (For Time-dependent pulse operator)");
    }
}

double getAngle(double degree) {
    return degree * (M_PI / 180.0);  // Radian
} 
/**
 * evolution=vector. calPropagatorGcce expands kron(Upulses[i], I_bdim) again on every
 * step; nothing in it depends on the step, so the vector path expands once per cluster.
 */
void prepareGcceVectorWork(GcceWork* w, MatrixXcd* Upulses, int npulse, int bdim){

    w->UpulsesExpanded.resize(npulse);
    for (int ipulse=0; ipulse<npulse; ipulse++){
        w->UpulsesExpanded[ipulse] = kron(Upulses[ipulse],MatrixXcd::Identity(bdim,bdim));
    }
}


/**
 * evolution=vector. The same operator as calPropagatorGcce, applied to a state instead
 * of assembled. That function builds, for npulse pulses,
 *     Utotal = Ufree_N (Upulse_{N-1} Ufree_{N-1}) ... (Upulse_0 Ufree_0)
 * so applying it to a state runs right to left: Ufree_0 first, then alternating the
 * pulse and the next Ufree. Every application is a matrix-VECTOR product, O(n^2), where
 * assembling Utotal costs a chain of O(n^3) matrix-matrix products.
 *
 * Ufree(tau) is never formed either -- see applyUfreeGcce. Requires propagator=eigen;
 * Config_setEvolution documents why, and calCoherenceGcce enforces it.
 */
const Eigen::VectorXcd& propagateStateGcce(const Eigen::VectorXcd& psi0, Pulse* pulse, double tfree, GcceWork* w){

    double** sequence        = Pulse_getSequence(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);
    int     npulse           = Pulse_getNpulse(pulse);

    // A pulse sequence repeats its delays -- CPMG runs npulse-1 identical inner taus --
    // and the matrix path has always exploited that through sequence_indices. Do the
    // same here: exp(-i*lambda*tau) costs one sincos per eigenvalue, so recomputing it
    // for a tau already seen is the most expensive thing to get wrong. An identical tau
    // gives an identical phase vector, so reusing it is bit-for-bit the same as
    // recomputing it.
    if ((int)w->phases.size() < npulse+1){ w->phases.resize(npulse+1); }
    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        if (sequence_indices[ipulse] == ipulse){
            w->evals_tau = w->evals * (tfree * sequence[ipulse][2]);
            phaseVectorGcce(w->evals_tau, w->phases[ipulse]);
        }
    }

    w->psi = psi0;
    applyUfreeGcce(w, w->phases[sequence_indices[0]]);

    for (int ipulse=1; ipulse<npulse+1; ipulse++){
        w->psiscratch.noalias() = w->UpulsesExpanded[ipulse-1] * w->psi;
        w->psi.swap(w->psiscratch);
        applyUfreeGcce(w, w->phases[sequence_indices[ipulse]]);
    }

    return w->psi;
}


MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pulse, double tfree, MatrixXcd* Upulses, GcceWork* w){

    double** sequence        = Pulse_getSequence(pulse);
    int*    sequence_indices = Pulse_getSequenceIndices(pulse);
    int     qdim             = QubitArray_dim(qa);
    int     npulse           = Pulse_getNpulse(pulse);

    /////////////////////////////////////////////////////////////
    // Propagator for total
    MatrixXcd Utotal;
    // Reused across the nstep loop instead of new/delete on every call.
    if ((int)w->Ufrees.size() < npulse+1){ w->Ufrees.resize(npulse+1); }

    // Propagator(total) =  U(tauN+1) (U_pulseN) ... (U_pulse1) U(tau1)
    // Pulse Index  0    1    2    3    4    5    6    7    8    |
    // Pulse (8#)   |____|____|____|____|____|____|____|____|____|
    //              |tau1|tau2|           .......           |tau9|
    //              |seq1|seq2|           .......           |seq9|
    // Pulse delay  0                                          tfree
    for (int ipulse=0; ipulse<npulse+1; ipulse++){
        
        // Propagator for free evolution
        double       tau = tfree * sequence[ipulse][2];
        int sameTauIndex = sequence_indices[ipulse];

        MatrixXcd Upulse;
        if (ipulse < npulse){
            Upulse = Upulses[ipulse];
        } else {
            Upulse = MatrixXcd::Identity(qdim, qdim);
        }
        int bdim = Htot.rows() / Upulse.rows();
        Upulse   = kron(Upulse, MatrixXcd::Identity(bdim, bdim));
        
        // Calculate free evolution operators (Ufree)
        if (sameTauIndex == ipulse){ // if Ufree haven't been calculated
            if (w->useEigen){
                // U = exp(-iHtau), from the decomposition prepareGcceWork already did:
                //   exp(-i*H*tau) = M * diag(exp(-i*lambda_k*tau)) * M^dag
                // Htot does not change with the step, so one decomposition serves every
                // tau instead of a matrix exponential per step.
                w->evals_tau = w->evals * tau;
                phaseVectorGcce(w->evals_tau, w->phase);
                w->scratch.noalias()        = w->M * w->phase.asDiagonal();
                w->Ufrees[ipulse].noalias() = w->scratch * w->Madj;
            }else{
                // propagator=expm : Eigen's general matrix exponential, as CCEX has
                // always done it. Kept for reproducing numbers generated before the
                // eigendecomposition path existed.
                w->Ufrees[ipulse] = ((-1.0) * doublec(0.0,1.0) * Htot * tau).exp();
            }
        }
        else{ // if Ufree have been calculated
            // get previously calculated Ufree for the same tau
            w->Ufrees[ipulse] = w->Ufrees[sameTauIndex];
        }
        const MatrixXcd& Ufree = w->Ufrees[ipulse];

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
    }
    

    return Utotal;
}

MatrixXcd cal_Upulse_when_pulseiter_isFalse(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma){
    MatrixXcd exponent = ((-1.0) * doublec(0.0,1.0) * sigma * angle / 2.0); 
    MatrixXcd Upulse = exponent.exp();
    return Upulse;
}

MatrixXcd cal_TDUpulse_when_pulseiter_isTrue(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, double angle, int ipulse) {

    int nqubit             = QubitArray_getNqubit(qa);
    double detuning_factor = Pulse_getPulseDetuningFactor(pulse);
    char axis              = pulse->pulse_axes[ipulse];
    double Ediff           = calQubit_Energy_difference(qa, Hq);
    std::cout << "Ediff : " << Ediff << std::endl;
    std::cout << "Qubit Hamiltonian : \n" << Hq  << std::endl;

    MatrixXcd U_total = MatrixXcd::Identity(1, 1); // 초기 유니타리

    // 전체 유니타리 구성
    for (int iq = 0; iq < nqubit; iq++) {
        MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa, iq);
        MatrixXcd beta  = QubitArray_getQubit_i_beta(qa, iq);
        std::cout << "alpha: \n" << alpha << std::endl;
        std::cout << "beta : \n" << beta << std::endl;
        exit(1);
        if (alpha.size() == 0 || beta.size() == 0) {
            fprintf(stderr, "[Error] alpha or beta not set for qubit %d\n", iq);
            exit(EXIT_FAILURE);
        }

        MatrixXcd* sigmas = getGeneralPauliOperators(alpha, beta);
        MatrixXcd ox = getPauliOperator(sigmas, 'X');
        MatrixXcd oy = getPauliOperator(sigmas, 'Y');
        MatrixXcd tempTDUpulse = cal_TDUpulse(qa, Hq, pulse, ox, oy, axis, angle, Ediff, detuning_factor, pulse->pulse_time, 10000);
        U_total = kron(U_total, tempTDUpulse);

        delete[] sigmas;
    }

    return U_total;
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
        // MatrixXcd Upulse_iq = ((-1.0) * doublec(0.0,1.0) * sigmaExpanded_x); // Before 25.08.08
        MatrixXcd exponent_iq = ((-1.0) * doublec(0.0,1.0) * sigmaExpanded_x * M_PI / 2.0); // Changed 25.08.08
        MatrixXcd Upulse_iq = exponent_iq.exp(); // Changed 25.08.08
        Upulse = Upulse * Upulse_iq;
    }
    return Upulse;
}

void calUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse){

    // Pulse Operator Up  = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t
    int     npulse        = Pulse_getNpulse(pulse);
    double* pulse_angles  = Pulse_getPulseAngles(pulse);
    char*   pulse_axes    = Pulse_getPulseAxes(pulse);
    bool    pulseiter     = Pulse_getPulseiter(pulse);
    int     qdim          = QubitArray_dim(qa);

    /////////////////////////////////////////////////////////////
    // Propagator for pulse
    //  Upulse = exp(-i*sigma_x*pi/2) -> Pi-Pulse
    for (int ipulse=0; ipulse<npulse; ipulse++){
        MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);
        double angle = getAngle(pulse_angles[ipulse]);  // Radian
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
        Upulses[ipulse] = Upulse;
    }
    //return Upulses;
}

//MatrixXcd cal_Upulse_withDetuning(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma, MatrixXcd sigma_z){
//
//    printf("Function: cal_Upulse_withDetuning \n");
//    // Pulse Operator Up  = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t
//    double detuning_factor = Pulse_getPulseDetuningFactor(pulse);
//    double ptime = (pulse->pulse_time) / 1e6;
//    //std::cout << "Detun/Omega: " << detuning / angle * ptime << std::endl;
//
//    MatrixXcd factor = (detuning * ptime * sigma_z) + (angle * sigma);
//    MatrixXcd U_detuning = ((-1.0) * doublec(0.0,1.0) * factor / 2.0).exp(); 
//    return U_detuning;
//}

double calQubit_Energy_difference(QubitArray* qa, MatrixXcd Hq){
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    double Ediff = 0.0;
    int alphaIdx;
    int betaIdx;
    //std::cout << "psia: " << psia << std::endl;
    //std::cout << "psib: " << psib << std::endl;
    for (int i = 0; i < psia.rows(); ++i) {
        if (psia(i, 0) != std::complex<double>(0.0, 0.0)) {
            alphaIdx = i;
            //std::cout << "non-zero element idx for pisa: " << i << std::endl;
        }
        if (psib(i, 0) != std::complex<double>(0.0, 0.0)) {
            betaIdx = i;
            //std::cout << "non-zero element idx for psib: " << i <<std::endl;
        }
    }
    Ediff = std::abs(Hq(alphaIdx, alphaIdx).real() - Hq(betaIdx, betaIdx).real());
    //std::cout << "Ediff: " << Ediff << std::endl;
    return Ediff;
}

void calTDUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse, MatrixXcd Hq){

    // Pulse Operator Up  = exp [ -(i/2) * sigma * theta ] where theta = (Rabi_Freq) * t
    int     npulse        = Pulse_getNpulse(pulse);
    double* pulse_angles  = Pulse_getPulseAngles(pulse);
    char*   pulse_axes    = Pulse_getPulseAxes(pulse);
    bool    pulseiter     = Pulse_getPulseiter(pulse);
    int     qdim          = QubitArray_dim(qa);
    double detuning_factor = Pulse_getPulseDetuningFactor(pulse);

    /////////////////////////////////////////////////////////////
    // Calculate time-dependent pulse operator. 
    //printf("Calculating pulse energy level difference! \n");
    double Ediff = calQubit_Energy_difference(qa, Hq);

    for (int ipulse=0; ipulse<npulse; ipulse++){
        MatrixXcd Upulse = MatrixXcd::Identity(qdim,qdim);
        double angle = getAngle(pulse_angles[ipulse]);  // Radian
        if (!pulseiter){
            MatrixXcd* sigmas = QubitArray_PauliOperator_fromPsiaPsib(qa);
            MatrixXcd  ox = getPauliOperator(sigmas, 'X');
            MatrixXcd  oy = getPauliOperator(sigmas, 'Y');
            char axis     = pulse_axes[ipulse];
            // printf(" || Calculating Time-dependent pulse operator || \n");
            Upulse = cal_TDUpulse(qa, Hq, pulse, ox, oy, axis, angle, Ediff, detuning_factor, pulse->pulse_time, 10000);
            delete[] sigmas;

        } else{
            printf("There's code developed for pulse iter time-dependent pulse operator!\n");
            Upulse = cal_Upulse_when_pulseiter_isTrue(qa, pulse, angle, ipulse);
            exit(1);
        }
        Upulses[ipulse] = Upulse;
    }
}

MatrixXcd cal_TDUpulse(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, MatrixXcd ox, MatrixXcd oy, char axis, double pulse_angle, double Ediff, double detuning_factor, double pulse_time, int Nstep) {

    //printf("Function: cal_TDUpulse\n");

    // Omega(Rabi_frequency) = Angle * pulse_duration_time //
    //double Omega = pulse_angle / pulse->pulse_time * 1e6;
    
    double ptime = pulse_time / 1e6;
    double Omega_Rabi = pulse_angle / ptime;

    // detuning = Ediff - omega // 
    double omega    = Ediff * detuning_factor;
    double detuning = Ediff - omega; // unit : rad*kHz
    double phase    = getPhase(axis);
    double Omega    = std::sqrt(Omega_Rabi*Omega_Rabi + detuning*detuning);

    const int qdim = ox.rows();
    double dt = ptime / static_cast<double>(Nstep);

     printf("Qubit Hamiltonian : \n");
     std::cout << Hq << std::endl;
     printf("\n");
     std::cout << "Pulse angle: " << pulse_angle << "  rad"  << std::endl;
     std::cout << "Pulse Time : " << pulse_time << "  ns" << std::endl;
     //std::cout << "dt         : " << dt*1e6 << "  ns" << std::endl;
     std::cout << "w0         : " << Ediff*0.001  << "  rad*MHz" << std::endl;
     std::cout << "Rabi       : " << Omega_Rabi*0.001  << "  rad*MHz" << std::endl;
     std::cout << "w          : " << omega     *0.001  << "  rad*MHz" << std::endl;
     std::cout << "Omega      : " << Omega     *0.001  << "  rad*MHz" << std::endl;
     std::cout << "detuning   : " << detuning  *0.001  << "  rad*MHz" << std::endl;
     //std::cout << "phase      : " << phase                      << std::endl;

    MatrixXcd TDUpulse = MatrixXcd::Identity(qdim, qdim);
    //MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);
    for (int n = 0; n < Nstep; ++n) {
        double t = n * dt;
        //MatrixXcd H = Hq + Omega * cos(omega * t) * ox / 2.0;
        MatrixXcd H = Hq + Omega * ( (cos(omega*t + phase)*ox) + (sin(omega*t + phase)*oy) ) / 2.0;
        MatrixXcd U = ((-1.0) * doublec(0.0, 1.0) * H * dt).exp();
        TDUpulse = U * TDUpulse;  
    }
    std::cout << "TDUpuslse : \n" << TDUpulse << std::endl;
    //exit(1);

    return TDUpulse;
}
