#include "../include/simulator.h"
#include "../include/hamiltonian.h"

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

    // Bath Hamiltonian
    MatrixXcd zeros_bdim = MatrixXcd::Zero(bdim,bdim);
    MatrixXcd Hb         = zeros_bdim; // total bath spin Hamiltonian
    MatrixXcd Hb_single  = zeros_bdim; // bath spin's total single spin Hamiltonian
    MatrixXcd Hb_pair    = zeros_bdim; // bath spin's total paired spin Hamiltonian

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
        Hb_single += expandHamiltonian(sigmas, Hbi, nspin, ib);
        ///////////////////////////////////////////////////////////
    }

    ///////////////////////////////////////////////////////////
    // Pair correlation
    for (int i=0; i<nspin; i++){
        for (int j=i+1; j<nspin; j++){
            MatrixXcd Hbij = BathArray_InteractionHamil(ba,sigmas,i,j);
            Hb_pair += Hbij;

            //char mesg[500];
            //sprintf(mesg,"Spin[%d] - Spin[%d] (%s) ",i,j,"Pair");
            //printInlineMatrixXcd(mesg,Hbij);
        }
    }
    ///////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////
    // Bath Hamiltonian
    Hb = Hb_single + Hb_pair;
    ///////////////////////////////////////////////////////////

    return Hb;
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

MatrixXcd* HamilQubitBathSecularApp(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas, Config* cnf){

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int bdim = BathArray_dim(ba);

    // Qubit - Bath Interaction Hamiltonian
    MatrixXcd Hqb = MatrixXcd::Zero(bdim,bdim);

    // nqubit > 1, err
    if (nqubit > 1){
        fprintf(stderr,"Error : nqubit > 1, current nqubit = %d\n",nqubit);
        exit(1);
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
        MatrixXcd Hqbij = expandHamiltonian(bsigmas, Hqbij_tmp, nspin, ib);
        char mesg[500];
        //sprintf(mesg,"Qubit - Spin[%d]", ib);
        //printInlineMatrixXcd(mesg,Hqbij);
 
        Hqb += Hqbij;
    }

     //////////////////////////////////////////////////////////////////////////////////////
    // Knight shift for transverse magnetic field
    MatrixXcd Hqb_knight_a = MatrixXcd::Zero(bdim,bdim);
    MatrixXcd Hqb_knight_b = MatrixXcd::Zero(bdim,bdim);

    // (3|alpha_eig|-2)/D * (gamma_e/ gamma_n) * mat(Axx Axy Axz Ayx Ayy Ayz 0 0 0)
    if (Config_getKnight(cnf)){

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

        double alpha_abs = (alpha_eig.adjoint() * alpha_eig)(0,0).real();
        double beta_abs  = ( beta_eig.adjoint() *  beta_eig)(0,0).real();
        ////////////////////////////////////////////

        // get gyros
        double qgyro = QubitArray_getQubit_i_gyro(qa,iq);

        // Bfield
        float* bfield = Config_getBfield(cnf);
        MatrixXcd bfield_vector = MatrixXcd::Zero(1,3);
        bfield_vector(0,0) = (double)bfield[0];
        bfield_vector(0,1) = (double)bfield[1];
        bfield_vector(0,2) = (double)bfield[2];

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

            // Vector K = B * g_tensor (1x3)
            MatrixXcd Bvec_gten_alpha = bfield_vector * g_tensor_alpha;
            MatrixXcd Bvec_gten_beta = bfield_vector * g_tensor_beta;

            MatrixXcd Hqb_knight_a_tmp = calHamiltonianSingleInt(Bvec_gten_alpha, bsigmas[ib]); 
            MatrixXcd Hqb_knight_b_tmp = calHamiltonianSingleInt(Bvec_gten_beta , bsigmas[ib]); 

            Hqb_knight_a += expandHamiltonian(bsigmas, Hqb_knight_a_tmp, nspin, ib);
            Hqb_knight_b += expandHamiltonian(bsigmas, Hqb_knight_b_tmp, nspin, ib);

        } 

        // free
        for (int i=0; i<nqubit; i++){
            delete[] qsigmas[i]; 
        }
        delete[] qsigmas;
    }



    //////////////////////////////////////////////////////////////////////////////////////
    // HF mediated term (Pair correlation)
    MatrixXcd Hqb_medi_a = MatrixXcd::Zero(bdim,bdim);
    MatrixXcd Hqb_medi_b = MatrixXcd::Zero(bdim,bdim);

    if (Config_getHfmedi(cnf)){

        // Qubit effective two states
        MatrixXcd psia = QubitArray_getPsia(qa);
        MatrixXcd psib = QubitArray_getPsib(qa);

        // Pauli matrices of Qubit
        MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);

        // Projected qubit operator
        // sigma_ab = (sab_x, sab_y, sab_z)
        // sigma_ba = (sba_x, sba_y, sba_z) //sab is not an operator
        MatrixXcd sigma_ab = projOperators(qsigmas[iq],psia,psib); // <psia|S|psib>
        MatrixXcd sigma_ba = projOperators(qsigmas[iq],psib,psia); // <psia|S|psib>

        // Eigenvalues for the eigenstates close to alpha and beta
        int qdim = QubitArray_dim(qa);
        MatrixXcd Hq = QubitArray_TotalHamil(qa,qsigmas,Config_getBfield(cnf));

        Eigen::SelfAdjointEigenSolver<MatrixXcd> es(qdim);
        es.compute(Hq);

        MatrixXcd eigenVectors = es.eigenvectors();
        VectorXcd eigenValues = es.eigenvalues();

        // The closest eigen state from the alpha state
        int ai = close_state_index(psia, eigenVectors); // Find <ei|alpha> > condidence (0.95)
        int bi = close_state_index(psib, eigenVectors); // Find <ei|alpha> > condidence (0.95)

        double omega_alpha = eigenValues[ai].real();
        double omega_beta  = eigenValues[bi].real();

        // Main part of HF medi (Virtual flip of the qubit) simulation.
        // For the equation, see N. Zhao, Nature Nanotechnology, 6, 242-246 (2011)
        for (int ib=0; ib<nspin; ib++){
            for (int jb=ib+1; jb<nspin; jb++){
                MatrixXcd tensor_iq_ib = BathArray_getBath_i_hypf_j(ba,ib,iq);
                MatrixXcd tensor_iq_jb = BathArray_getBath_i_hypf_j(ba,jb,iq);

                // Main result
                MatrixXcd* Iib_AiAj_Ijb_ab = calHamiltonianHeteroInt_o2(sigma_ab, sigma_ba, bsigmas, tensor_iq_ib, tensor_iq_jb, nspin, ib, jb);

                Hqb_medi_a += Iib_AiAj_Ijb_ab[0]/(omega_alpha - omega_beta );
                Hqb_medi_b += Iib_AiAj_Ijb_ab[1]/(omega_beta  - omega_alpha);

                delete[] Iib_AiAj_Ijb_ab;
            }
        }

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

    // Calculate all Hqb
    MatrixXcd* Hqb_PsiaPsib = new MatrixXcd[2];
    Hqb_PsiaPsib[0] = ms_a * Hqb + Hqb_medi_a + Hqb_knight_a;
    Hqb_PsiaPsib[1] = ms_b * Hqb + Hqb_medi_b + Hqb_knight_b;

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
        exit(1);
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

    fprintf(stderr,"Error close_state_index (HF medi) \n");
    exit(1);
}
