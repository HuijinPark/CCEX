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

    // Bath Hamiltonian
    MatrixXcd Hb = MatrixXcd::Zero(bdim,bdim);
    Hb = BathArray_TotalHamil(ba, sigmas, bfield);

    // Subtract the double counted disorders
    if (Config_getNstate(cnf)>0){
        MatrixXcd Hb_i_overlapped = MatrixXcd::Zero(bdim,bdim);
        for (int ib = 0; ib < nspin; ib++){
            double disorder = 0.0;
            for (int jb = 0; jb < nspin; jb++){
                disorder += BathArray_getBath_i_disorder_j(ba,ib,jb); // if ensemble calculation, disorder = 0.0
            }
            MatrixXcd vecDisorder = calDetuningVector(disorder);
            MatrixXcd Hb_i_overlapped_tmp = calHamiltonianSingleInt(vecDisorder,sigmas[ib]);
            MatrixXcd Hb_i_overlapped = expandHamiltonian(sigmas, Hb_i_overlapped_tmp, nspin, ib);
            Hb -= Hb_i_overlapped;
        }
    }

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
        // std::cout << "H_S_I [" << ib << "] = " <<  Hqbij << std::endl;
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
