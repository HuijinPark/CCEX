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

    //! hyperfine mediated term
    //! if (hfmedi){

    //! }

    // get qubit sublevel
    // MatrixXcd* qsigmas_general = QubitArray_PauliOperator_fromPsiaPsib(qa);
    // MatrixXcd qsz = qsigmas_general[3];
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    double ms_a = findZbasisSubLevel(psia);
    double ms_b = findZbasisSubLevel(psib);

    // Hqb
    MatrixXcd* Hqb_PsiaPsib = new MatrixXcd[2];
    Hqb_PsiaPsib[0] = ms_a * Hqb;
    Hqb_PsiaPsib[1] = ms_b * Hqb;

    return Hqb_PsiaPsib;
}
