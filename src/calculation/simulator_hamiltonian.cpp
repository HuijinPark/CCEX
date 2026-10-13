#include "../../include/calculation/simulator.h"
#include "../../include/calculation/hamiltonian.h"

MatrixXcd HamilQubit(QubitArray* qa, BathArray* ba, MatrixXcd** sigmas, Config* cnf){

    int nspin = BathArray_getNspin(ba);
    int nqubit = QubitArray_getNqubit(qa);
    int qdim = QubitArray_dim(qa);
    float* bfield = Config_getBfield(cnf);

    // Qubit Hamiltonian
    MatrixXcd Hq = MatrixXcd::Zero(qdim,qdim);
    
    Hq = QubitArray_TotalHamil(qa, sigmas, bfield);
    // Subtract the double counted Overhauser fields
    MatrixXcd Hq_i_overlapped = MatrixXcd::Zero(qdim,qdim);
    for (int iq = 0; iq < nqubit; iq++){
        double overhaus = 0.0;
        for (int ib = 0; ib < nspin; ib++){
            overhaus += BathArray_getBath_i_overhaus_j(ba,ib,iq); // if ensemble calculation, overhaus = 0.0 (ms = 0.0)
        }
        MatrixXcd vecOverhaus = calOverhauserVector(overhaus);
        MatrixXcd Hq_i_overlapped_tmp = calHamiltonianSingleInt(vecOverhaus,sigmas[iq]);
        MatrixXcd Hq_i_overlapped = expandHamiltonian(sigmas, Hq_i_overlapped_tmp, nqubit, iq);
        Hq -= Hq_i_overlapped;
    }

    return Hq;
}

MatrixXcd HamilBath(BathArray* ba, MatrixXcd** sigmas, Config* cnf){

    int nspin = BathArray_getNspin(ba);
    int bdim = BathArray_dim(ba);
    float* bfield = Config_getBfield(cnf);
    int nstate = Config_getNstate(cnf);

    // Bath Hamiltonian. Each term is added into its sum at the entries it has
    // (addExpandedHamiltonian, BathArray_addInteractionHamil) instead of being built
    // as a bdim x bdim matrix and added whole: the sums are the same bit for bit, and
    // for 14 spin-1/2 the build went from 161 s to the time of a few passes (2026-09-30).
    MatrixXcd Hb_single  = MatrixXcd::Zero(bdim,bdim); // bath spin's total single spin Hamiltonian
    MatrixXcd Hb_pair    = MatrixXcd::Zero(bdim,bdim); // bath spin's total paired spin Hamiltonian

    ///////////////////////////////////////////////////////////
    // Single correlation 

    for (int ib=0; ib<nspin; ib++){
        ///////////////////////////////////////////////////////////
        // Zeeman Hamiltonian     : -1.0 * gamma * vec{B0} * vec{Spin Operator}
        // Detuning Hamiltonian   : [0,0,detuning] * vec{Spin Operator}
        // Quad / ZFS Hamiltonian : vec{Spin Operator} 
        //                            * tensor{quad/zfs} * vec{Spin Operator}
        // Disorder Hamiltonian   : [0,0,disorder] * vec{Spin Operator} (Only single app.)
        ///////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////
        // Declare i-th bath Hamiltonian for each interaction
        int bdim_i = BathArray_dimBath_i(ba,ib); // i-th spin dimension
        MatrixXcd zeros_bdim_i  = MatrixXcd::Zero(bdim_i,bdim_i);
        MatrixXcd Hbi           = zeros_bdim_i; // single   Hamiltonian of i-th bath spin 
        MatrixXcd Hbi_Zeeman    = zeros_bdim_i; // Zeeman   Hamiltonian
        MatrixXcd Hbi_Detuning  = zeros_bdim_i; // Detuning Hamiltonian
        MatrixXcd Hbi_Quad      = zeros_bdim_i; // Quad/ZFS Hamiltonian
        MatrixXcd Hbi_Disorder  = zeros_bdim_i; // Disorder Hamiltonian (Only single app)
        bool      Hbi_Disorder_rm_overlap = true; // Remove overlapped disorder
                                                  // Subtract the double counted disorders

        // Calculate each Hamiltonian
        Hbi_Zeeman   = BathArray_ZeemanHamil(ba,sigmas,ib,bfield);
        Hbi_Detuning = BathArray_DetuningHamil(ba,sigmas,ib);
        Hbi_Quad     = BathArray_QuadHamil(ba,sigmas,ib);

        if (nstate > 0){    // Single-sample app.
            Hbi_Disorder = BathArray_DisorderHamil(ba,sigmas,ib, \
                                                    Hbi_Disorder_rm_overlap);
        }else{ // Ensemble app.
            Hbi_Disorder = zeros_bdim_i; 
        }
        ///////////////////////////////////////////////////////////

        ///////////////////////////////////////////////////////////
        // Sum all single interaction of i-th spin
        Hbi = Hbi_Zeeman + Hbi_Detuning + Hbi_Quad + Hbi_Disorder;
        ///////////////////////////////////////////////////////////
        char mesg[500];
        //sprintf(mesg,"Spin[%d] : %s",ib,"ZE");
        //printInlineMatrixXcd(mesg,Hbi_Zeeman);
        //sprintf(mesg,"Spin[%d] : %s",ib,"Detun");
        //printInlineMatrixXcd(mesg,Hbi_Detuning);
        //sprintf(mesg,"Spin[%d] : %s",ib,"Quad");
        //printInlineMatrixXcd(mesg,Hbi_Quad);
        //sprintf(mesg,"Spin[%d] : %s",ib,"Disd");
        //printInlineMatrixXcd(mesg,Hbi_Disorder);

        ///////////////////////////////////////////////////////////
        // Expand dimension
        addExpandedHamiltonian(Hb_single, sigmas, Hbi, nspin, ib);
        ///////////////////////////////////////////////////////////
    }

    ///////////////////////////////////////////////////////////
    // Pair correlation
    for (int i=0; i<nspin; i++){
        for (int j=i+1; j<nspin; j++){
            BathArray_addInteractionHamil(Hb_pair,ba,sigmas,i,j);

            //char mesg[500];
            //sprintf(mesg,"Spin[%d] - Spin[%d] (%s) ",i,j,"Pair");
            //printInlineMatrixXcd(mesg,Hbij);
        }
    }
    ///////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////
    // Bath Hamiltonian : Hb_single + Hb_pair (total bath spin Hamiltonian)
    Hb_single += Hb_pair;
    ///////////////////////////////////////////////////////////

    return Hb_single;
}

MatrixXcd HamilQubitBath(QubitArray* qa, BathArray* ba, MatrixXcd** qsigmas, MatrixXcd** bsigmas, Config* cnf){

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int ntotspin = nqubit + nspin;
    int qdim = QubitArray_dim(qa);
    int bdim = BathArray_dim(ba);
    int totdim = qdim*bdim;

    // Pauli matrices
    MatrixXcd** sigmas = gatherSigmas(qsigmas,bsigmas,nqubit,nspin);

    // Qubit - Bath Interaction Hamiltonian
    MatrixXcd Hqb = MatrixXcd::Zero(totdim,totdim);
    
    for (int iq = 0; iq < nqubit; iq++){
        for (int ib = 0; ib < nspin; ib++){
            MatrixXcd tensor = BathArray_getBath_i_hypf_j(ba,ib,iq);
            MatrixXcd Hqbij = calHamiltonianHeteroInt(sigmas, tensor, ntotspin, iq, ib+nqubit);
            Hqb += Hqbij;
        }
    }

    // free
    for (int i=0; i<ntotspin; i++){
        delete[] sigmas[i];
    }
    delete[] sigmas;

    return Hqb;
}

// Hyperfine-mediated term of one branch (complete second order in the hyperfine coupling):
//   H_med^(alpha) = sum_{beta != alpha} E_ab E_ab^dagger / (w_a - w_b),  E_ab = sum_i (S_ab . A_i) . I_i,
// S_ab = <alpha|S|beta> with the adjoint on the bra; alpha is the eigenstate of Hq closest to psi and beta runs over
// every other eigenstate. E_ba = E_ab^dagger (A real, I Hermitian), so each term is Hermitian by construction, and
// expanding E E^dagger gives every ordered pair (i, j) and the i = j terms.
static MatrixXcd hfmediBranch(const MatrixXcd& psi, const MatrixXcd& eigenVectors, const VectorXcd& eigenValues,
                              MatrixXcd** qsigmas, BathArray* ba, MatrixXcd** bsigmas, int nspin, int bdim, int iq){
    const int ai = close_state_index(psi, eigenVectors);
    const VectorXcd ea = eigenVectors.col(ai);
    const double w_a = eigenValues[ai].real();
    // degenerate: within 1e-12 of the spread of the levels, which an energy shift (c I added to the ZFS) leaves as it is
    // (P1-03 of review 20261012-233058-452278; it was 1e-12 of the largest |level|)
    const double degenerate = 1e-12 * (eigenValues.real().maxCoeff() - eigenValues.real().minCoeff());
    // E = sum_i (s . A_i) . I_i
    auto coupling = [&](const MatrixXcd& s){
        MatrixXcd E = MatrixXcd::Zero(bdim,bdim);
        for (int ib = 0; ib < nspin; ib++){
            MatrixXcd A = BathArray_getBath_i_hypf_j(ba,ib,iq);    // 3x3, row = qubit component
            MatrixXcd Hi = calHamiltonianSingleInt(s * A, bsigmas[ib]);
            addExpandedHamiltonian(E, bsigmas, Hi, nspin, ib);
        }
        return E;
    };
    MatrixXcd H = MatrixXcd::Zero(bdim,bdim);
    for (int k = 0; k < eigenVectors.cols(); k++){
        if (k == ai) continue;                                    // beta != alpha
        const VectorXcd eb = eigenVectors.col(k);
        const double w_b = eigenValues[k].real();
        MatrixXcd s_ab(1,3);                                      // <alpha|S|beta>
        for (int c = 0; c < 3; c++) s_ab(0,c) = (ea.adjoint() * qsigmas[iq][c+1] * eb)(0,0);
        const MatrixXcd E_ab = coupling(s_ab);
        // beta degenerate with alpha (review 20261012-194049-af8e42 R1): decided on the coupling itself (P1-02 of review
        // 20261012-233058-452278) - with E_ab = 0 the term is 0 (S_ab = 0, e.g. +1 and -1 of a spin-1 at zero field, or a
        // hyperfine that vanishes), anything else has no second-order correction. Zero: below 1e-12 of the coupling a
        // unit S would have along each axis
        if (std::abs(w_a - w_b) <= degenerate){
            double scale = 0;
            for (int c = 0; c < 3; c++){
                MatrixXcd unit = MatrixXcd::Zero(1,3);
                unit(0,c) = 1.0;
                scale += coupling(unit).norm();
            }
            if (E_ab.norm() <= 1e-12 * scale) continue;
            fprintf(stderr,"Error : hfmedi: qubit level %d (energy %.10g) is degenerate with level %d (energy %.10g), "
                    "which it couples to through the hyperfine interaction; the second-order hyperfine-mediated term does "
                    "not exist there (a nonzero field or zero-field splitting lifts the degeneracy)\n", ai, w_a, k, w_b);
            ccexFatal();
        }
        H += (E_ab * E_ab.adjoint()) / (w_a - w_b);
    }
    return H;
}

MatrixXcd* HamilQubitBathSecularApp(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas, Config* cnf){

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int bdim = BathArray_dim(ba);

    // Qubit - Bath Interaction Hamiltonian
    MatrixXcd Hqb = MatrixXcd::Zero(bdim,bdim);

    // nqubit > 1, err
    if (nqubit > 1){
        fprintf(stderr,"Error : nqubit > 1, current nqubit = %d\n",nqubit);
        ccexFatal();   // reached per cluster by one rank: every rank stops (review 20261008-225534-82c666 R5)
    }

    // SzAzzIz + SzAzxIx + SzAzyIy
    // Here, Sz is defined for alpha,beta
    int iq = 0;

    MatrixXcd Hqbij;
    for (int ib = 0; ib < nspin; ib++){
        MatrixXcd tensor = BathArray_getBath_i_hypf_j(ba,ib,iq);
        MatrixXcd vector = MatrixXcd::Zero(3,1);
        vector(0,0) = tensor(2,0).real(); // Azx
        vector(1,0) = tensor(2,1).real(); // Azy
        vector(2,0) = tensor(2,2).real(); // Azz

        MatrixXcd Hqbij_tmp = calHamiltonianSingleInt(vector,bsigmas[ib]);
        char mesg[500];
        //sprintf(mesg,"Qubit - Spin[%d]", ib);
        //printInlineMatrixXcd(mesg,Hqbij_tmp);
 
        addExpandedHamiltonian(Hqb, bsigmas, Hqbij_tmp, nspin, ib);
    }

     //////////////////////////////////////////////////////////////////////////////////////
    // Knight shift for transverse magnetic field
    // (bdim x bdim only when the term is on; off, it is the +0 it contributes below)
    MatrixXcd Hqb_knight_a;
    MatrixXcd Hqb_knight_b;

    // (3|alpha_eig|-2)/D * (gamma_e/ gamma_n) * mat(Axx Axy Axz Ayx Ayy Ayz 0 0 0)
    if (Config_getKnight(cnf)){
        Hqb_knight_a = MatrixXcd::Zero(bdim,bdim);
        Hqb_knight_b = MatrixXcd::Zero(bdim,bdim);

        // Delta (D)
        MatrixXcd zfs_tensor = QubitArray_getIntmap_i_j(qa,iq,iq);
        double Delta = zfs_tensor(2,2).real() * 3.0 / 2.0;

        ////////////////////////////////////////////
        // Get alpha/beta_abs
        //
        // Qubit effective two states
        MatrixXcd psia = QubitArray_getPsia(qa);
        MatrixXcd psib = QubitArray_getPsib(qa);

        // Pauli matrices of Qubit
        MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);

        // Eigenstates close to alpha and beta (alpha_eig, beta_eig)
        int qdim = QubitArray_dim(qa);
        MatrixXcd Hq = QubitArray_TotalHamil(qa,qsigmas,Config_getBfield(cnf));

        Eigen::SelfAdjointEigenSolver<MatrixXcd> es(qdim);
        es.compute(Hq);

        MatrixXcd eigenVectors = es.eigenvectors();
        VectorXcd eigenValues = es.eigenvalues();

        // The closest eigen state from the alpha state
        int ai = close_state_index(psia, eigenVectors); // Find <ei|alpha> > condidence (0.95)
        int bi = close_state_index(psib, eigenVectors); // Find <ei|alpha> > condidence (0.95)

        VectorXcd alpha_eig = eigenVectors.col(ai);
        VectorXcd beta_eig  = eigenVectors.col(bi);

        // |m_s| of each branch as <S_z^2> (the norm of the eigenvector, used before, is always 1)
        MatrixXcd Sz2 = qsigmas[iq][3] * qsigmas[iq][3];
        double alpha_abs = (alpha_eig.adjoint() * Sz2 * alpha_eig)(0,0).real();
        double beta_abs  = ( beta_eig.adjoint() * Sz2 *  beta_eig)(0,0).real();
        ////////////////////////////////////////////

        // get gyros
        double qgyro = QubitArray_getQubit_i_gyro(qa,iq);

        // Bfield
        float* bfield = Config_getBfield(cnf);

        // Main part of simulation 
        for (int ib = 0; ib < nspin; ib++){
            double bgyro = BathArray_getBath_i_gyro(ba,ib);
            MatrixXcd tensor_iq_ib = BathArray_getBath_i_hypf_j(ba,ib,iq);
            MatrixXcd g_tensor_tmp = tensor_iq_ib;
            g_tensor_tmp(2,0) = doublec(0.0,0.0);
            g_tensor_tmp(2,1) = doublec(0.0,0.0);
            g_tensor_tmp(2,2) = doublec(0.0,0.0);
            MatrixXcd g_tensor_alpha = ((3.0*alpha_abs - 2.0)/Delta) * (qgyro / bgyro) * g_tensor_tmp;
            MatrixXcd g_tensor_beta  = ((3.0*beta_abs - 2.0)/Delta) * (qgyro / bgyro) * g_tensor_tmp;

            // Vector K = -gamma_n B * g_tensor (1x3): the nuclear Zeeman vector, not the bare field
            MatrixXcd zeeman_vec = calZeemanVector(bgyro, bfield);
            MatrixXcd Bvec_gten_alpha = zeeman_vec * g_tensor_alpha;
            MatrixXcd Bvec_gten_beta = zeeman_vec * g_tensor_beta;

            MatrixXcd Hqb_knight_a_tmp = calHamiltonianSingleInt(Bvec_gten_alpha, bsigmas[ib]); 
            MatrixXcd Hqb_knight_b_tmp = calHamiltonianSingleInt(Bvec_gten_beta , bsigmas[ib]); 

            addExpandedHamiltonian(Hqb_knight_a, bsigmas, Hqb_knight_a_tmp, nspin, ib);
            addExpandedHamiltonian(Hqb_knight_b, bsigmas, Hqb_knight_b_tmp, nspin, ib);

        } 

        // free
        for (int i=0; i<nqubit; i++){
            delete[] qsigmas[i]; 
        }
        delete[] qsigmas;
    }



    //////////////////////////////////////////////////////////////////////////////////////
    // HF mediated term (second order in the hyperfine coupling)
    MatrixXcd Hqb_medi_a;
    MatrixXcd Hqb_medi_b;

    if (Config_getHfmedi(cnf)){
        // Qubit effective two states
        MatrixXcd psia = QubitArray_getPsia(qa);
        MatrixXcd psib = QubitArray_getPsib(qa);

        // Pauli matrices of Qubit
        MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);

        // Eigenstates of the qubit Hamiltonian: every one but the branch's own is an intermediate state
        int qdim = QubitArray_dim(qa);
        MatrixXcd Hq = QubitArray_TotalHamil(qa,qsigmas,Config_getBfield(cnf));

        Eigen::SelfAdjointEigenSolver<MatrixXcd> es(qdim);
        es.compute(Hq);

        MatrixXcd eigenVectors = es.eigenvectors();
        VectorXcd eigenValues = es.eigenvalues();

        // Second-order (virtual flip of the qubit) correction, each branch on its own.
        // For the equation, see N. Zhao et al., Nat. Nanotechnol. 6, 242 (2011), SI Eq. 9,
        // summed over every ordered pair and i = j (complete second order, as PyCCE bath_mediated).
        Hqb_medi_a = hfmediBranch(psia, eigenVectors, eigenValues, qsigmas, ba, bsigmas, nspin, bdim, iq);
        Hqb_medi_b = hfmediBranch(psib, eigenVectors, eigenValues, qsigmas, ba, bsigmas, nspin, bdim, iq);

        // free
        for (int i=0; i<nqubit; i++){
            delete[] qsigmas[i]; 
        }
        delete[] qsigmas;
    }
    //////////////////////////////////////////////////////////////////////////////////////

    // Get qubit sublevel
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);
    double ms_a = findZbasisSubLevel(psia);
    double ms_b = findZbasisSubLevel(psib);

    // Calculate all Hqb : ms * Hqb + Hqb_medi + Hqb_knight, in that order.
    // A term that is off used to be a zero matrix. Adding +0 only turns a -0 of
    // ms * Hqb into +0, so the first one is kept as a scalar +0 and the second,
    // which then changes nothing, is left out - the same bits without the matrices.
    MatrixXcd* Hqb_PsiaPsib = new MatrixXcd[2];
    const double ms[2] = {ms_a, ms_b};
    const MatrixXcd* medi[2] = {&Hqb_medi_a, &Hqb_medi_b};
    const MatrixXcd* knight[2] = {&Hqb_knight_a, &Hqb_knight_b};
    for (int k=0; k<2; k++){
        if (Config_getHfmedi(cnf)){
            Hqb_PsiaPsib[k] = ms[k] * Hqb + *medi[k];
        }else{
            Hqb_PsiaPsib[k] = ((ms[k] * Hqb).array() + doublec(0.0,0.0)).matrix();
        }
        if (Config_getKnight(cnf)){
            Hqb_PsiaPsib[k] += *knight[k];
        }
    }

    return Hqb_PsiaPsib;
}

MatrixXcd** gatherSigmas(MatrixXcd** qsigmas, MatrixXcd** bsigmas, int nqubit, int nspin){

    int ntotspin = nqubit+nspin;

    MatrixXcd** sigmas = new MatrixXcd*[nqubit+nspin];

    for (int i=0; i<ntotspin; i++){
        sigmas[i] = new MatrixXcd[4];

        if (i<nqubit){
            sigmas[i][0] = qsigmas[i][0];
            sigmas[i][1] = qsigmas[i][1];
            sigmas[i][2] = qsigmas[i][2];
            sigmas[i][3] = qsigmas[i][3];
        }else{
            sigmas[i][0] = bsigmas[i-nqubit][0];
            sigmas[i][1] = bsigmas[i-nqubit][1];
            sigmas[i][2] = bsigmas[i-nqubit][2];
            sigmas[i][3] = bsigmas[i-nqubit][3];
        }
    }

    return sigmas;
}

MatrixXcd** gatherSigmas_singlequbit(MatrixXcd* qsigma, MatrixXcd** bsigmas, int nqubit, int nspin){

    if (nqubit > 1){
        fprintf(stderr,"Error gatherSigmas_singlequbit : nqubit > 1, current nqubit = %d\n",nqubit);
        ccexFatal();
    }

    int ntotspin = nqubit+nspin;

    MatrixXcd** sigmas = new MatrixXcd*[nqubit+nspin];

    for (int i=0; i<ntotspin; i++){
        sigmas[i] = new MatrixXcd[4];

        if (i<nqubit){
            sigmas[i][0] = qsigma[0];
            sigmas[i][1] = qsigma[1];
            sigmas[i][2] = qsigma[2];
            sigmas[i][3] = qsigma[3];
        }else{
            sigmas[i][0] = bsigmas[i-nqubit][0];
            sigmas[i][1] = bsigmas[i-nqubit][1];
            sigmas[i][2] = bsigmas[i-nqubit][2];
            sigmas[i][3] = bsigmas[i-nqubit][3];
        }
    }

    return sigmas;
} 

int close_state_index(MatrixXcd state, MatrixXcd eigenVectors){

    //double level_confidence = 0.95; //Original from pycce
    double level_confidence = 0.51;

    VectorXcd ev_state = eigenVectors.adjoint() * state;

    for (int i = 0; i < ev_state.size(); ++i) {

        double fidelity = std::norm(ev_state[i]);  // |⟨eigen_i | state⟩|^2

        if (fidelity > level_confidence){
            return i;
        }
    }

    fprintf(stderr,"Error close_state_index (HF medi): no eigenstate holds more than %.2f of the qubit state\n", level_confidence);
    ccexFatal();   // only the rank with this cluster gets here: the others would wait in MPI_Reduce (review 82c666 R5)
}
