#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include <iostream>
MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls){
    
    double stime = MPI_Wtime();
    double etime = 0.0;

    int nspin = BathArray_getNspin(ba);
    int bdim = BathArray_dim(ba);

    int nqubit = QubitArray_getNqubit(qa);
    if (nqubit > 1){
        fprintf(stderr,"Error : calCoherenceCce : nqubit > 1, current nqubit = %d\n",nqubit);
        exit(1);
    }

    ////////////////////////////////
    // Hamiltonian
    ////////////////////////////////

    // Declare the total Hamiltonian
    MatrixXcd Halpha = MatrixXcd::Zero(bdim,bdim);
    MatrixXcd Hbeta = MatrixXcd::Zero(bdim,bdim);

    // Pauli matrices
    MatrixXcd** bsigmas = BathArray_PauliOperators(ba);

    // Bath Hamiltonian
    MatrixXcd Hb = HamilBath(ba,bsigmas,cnf);
    //std::cout << "Hb" << Hb << std::endl;

    // Qubit - Bath Interaction Hamiltonian Hqb[0] = alpha * Hqb, Hqb[1] = beta * Hqb
    // Hqb[0] = alpha * (Sum_i Hqb[i]), Hqb[1] = beta * (Sum_i Hqb[i])
    MatrixXcd* Hqb = HamilQubitBathSecularApp(qa,ba,bsigmas,cnf); 
    //std::cout << "Hqb[0] (alpha)" << Hqb[0] << std::endl;
    //std::cout << "Hqb[1] (beta)" << Hqb[1] << std::endl;

    // Total Hamiltonian
    Halpha = Hb + Hqb[0];
    Hbeta = Hb + Hqb[1];

    //printInlineMatrixXcd("Halpha", Halpha);
    //printInlineMatrixXcd("Hbeta", Hbeta  );

    //etime = MPI_Wtime();
    //printf("          Wall time step1 = %.5f s\n", stime - etime);

    ////////////////////////////////
    // Diagonalize the Hamiltonian
    ////////////////////////////////
    MatrixXcd M_up (bdim,bdim);
    MatrixXcd M_down(bdim,bdim);
    MatrixXcd M(bdim,bdim);
    Eigen::SelfAdjointEigenSolver<MatrixXcd> es_up(bdim), es_down(bdim);
    es_up.compute(Halpha);
    es_down.compute(Hbeta);
    M_up   = es_up.eigenvectors();
    M_down = es_down.eigenvectors();
    M      = M_up.inverse()*M_down;
    Eigen::VectorXcd es_up_Val=es_up.eigenvalues();
    Eigen::VectorXcd es_down_Val=es_down.eigenvalues();

    //std::cout << "es_up" << es_up_Val  << std::endl;
    //std::cout << "\n\nes_down" << es_down_Val << std::endl;
    //etime = MPI_Wtime();
    //printf("          Wall time step2 = %.5f s\n", stime - etime);


    ////////////////////////////////
    // Initial state of bath
    ////////////////////////////////

    // Check if this is the ensemble calculation
    int nstate = Config_getNstate(cnf);
    bool isEnsemble = false;
    if (nstate == 0){isEnsemble = true;}

    MatrixXcd psi0;
    if (!isEnsemble){
        psi0 = BathArray_Psi0(ba);
    }

    ////////////////////////////////
    // Propagation
    ////////////////////////////////

    // Quantity
    char* quantity = Config_getQuantity(cnf);
    if (strcasecmp(quantity,"coherence")!=0){
        fprintf(stderr,"Error : calCoherenceCce : quantity is not coherence\n");
        exit(1);
    }

    // Npulse
    int npulse = Pulse_getNpulse(pls);
    double** sequence = Pulse_getSequence(pls);

    // Alloc the result variable
    int nstep = Config_getNstep(cnf);
    double deltat = (double)Config_getDeltat(cnf);
    MatrixXcd* result = new MatrixXcd[nstep];

    //etime = MPI_Wtime();
    //printf("          Wall time step3 = %.5f s\n", stime - etime);


    double tfree = 0.0;

    if (isEnsemble){

        MatrixXcd U_up; MatrixXcd U_down;
        MatrixXcd U_up2; MatrixXcd U_down2;

        MatrixXcd L_mat; 
        MatrixXcd L_a_even; MatrixXcd L_b_even;
        MatrixXcd L_a_odd; MatrixXcd L_b_odd;
        MatrixXcd L_a2_even; MatrixXcd L_b2_even;
        MatrixXcd L_a2_odd; MatrixXcd L_b2_odd;
        MatrixXcd sum_L_a; MatrixXcd sum_L_b;

        for (int i = 0; i < nstep; i++) {
            //U_up and U down: Propagator for the spin-up Hamiltonian
            //Note that the free evolution time between the pulses 
            // is t_free/2.0 for the Hahn-echo sequence.
            U_up   = MatrixXcd::Zero(bdim,bdim);
            U_down = MatrixXcd::Zero(bdim,bdim);
            U_up2   = MatrixXcd::Zero(bdim,bdim);
            U_down2 = MatrixXcd::Zero(bdim,bdim);

            //T2_star (Ramsey)
            if(npulse==0){

                Eigen::VectorXcd es_up_Val_tau      =es_up_Val*tfree;
                Eigen::VectorXcd es_down_Val_tau    =es_down_Val*tfree;
                Eigen::VectorXcd U_up               =es_up_Val_tau.array().cos() + doublec(0,-1)*es_up_Val_tau.array().sin();
                Eigen::VectorXcd U_down             =es_down_Val_tau.array().cos() + doublec(0,-1)*es_down_Val_tau.array().sin();

                // Coherence function
                //L_mat     =  M * U_down * M.adjoint()* U_up.adjoint();
                L_mat     =  (M.array().rowwise() * U_down.transpose().eval().array() ).matrix() 
                            * (M.adjoint().array().rowwise() * U_up.adjoint().array()).matrix();



            }else{
            //T2 (CPMG sequence)
            //if you want to use the Hahn-echo, set the N=1
                float tau = tfree/(2*npulse);

                //define diagonal matrix to vector 
                Eigen::VectorXcd es_up_Val_tau      = es_up_Val*tau;
                Eigen::VectorXcd es_down_Val_tau    = es_down_Val*tau;
                Eigen::VectorXcd U_up               = es_up_Val_tau.array().cos() + doublec(0,-1)*es_up_Val_tau.array().sin();
                Eigen::VectorXcd U_down             = es_down_Val_tau.array().cos() + doublec(0,-1)*es_down_Val_tau.array().sin();

                Eigen::VectorXcd es_up_Val_tau2     = es_up_Val*2*tau;
                Eigen::VectorXcd es_down_Val_tau2   = es_down_Val*2*tau;
                Eigen::VectorXcd U_up2              = es_up_Val_tau2.array().cos() + doublec(0,-1)*es_up_Val_tau2.array().sin();
                Eigen::VectorXcd U_down2            = es_down_Val_tau2.array().cos() + doublec(0,-1)*es_down_Val_tau2.array().sin();
                

                //pulse multiple
                L_a_even = M.array().rowwise()           * U_down.transpose().eval().array();   //M * U_down
                L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array();       //M.adjoint * U_up.adjoint

                if(npulse % 2 != 0){
                    L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().eval().array();  //M.adjoint() * U_up
                    L_b_odd = M.array().rowwise()           * U_down.adjoint().array();  //M * U_down.adjoint
                }

                switch (npulse){
                default :
                case 3:
                    L_a2_even = M.array().rowwise()           * U_down2.transpose().eval().array();  //M * U_down2
                    L_b2_even = M.adjoint().array().rowwise() * U_up2.adjoint().array();      //M.adjoint * U_up2.adjoint
                case 2:
                    L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().eval().array();     //M.adjoint * U_up2
                    L_b2_odd = M.array().rowwise()           * U_down2.adjoint().array();     //M * U_down2.adjoint
                    break;
                case 1:
                    break;
                }

                sum_L_a = L_a_even; // M * U_down
                sum_L_b = L_b_even; // M.adjoint * U_up.adjoint


                int loopN=1;
                for (int nn = 1; nn <= npulse; nn ++){
                    if ( nn == npulse){
                        if ( nn % 2==0){  //even number pulse
                            //L_a = L_a * M * U_down;
                            //L_b = M.adjoint() * U_up.adjoint() * L_b;

                            //L_a_even = M.array().rowwise()           * U_down.transpose().eval().array(); //M * U_down
                            //L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array(); //M.adjoint * U_up.adjoint

                            sum_L_a = sum_L_a * L_a_even;
                            sum_L_b = L_b_even * sum_L_b;

                        }else{ //odd number pulse
                            //L_a = L_a * M.adjoint() * U_up;
                            //L_b = M * U_down.adjoint() * L_b;

                            //L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().eval().array(); //M.adjoint() * U_up
                            //L_b_odd = M.array().rowwise()           * U_down.adjoint().array(); //M * U_down.adjoint
                            
                            sum_L_a = sum_L_a * L_a_odd;
                            sum_L_b = L_b_odd * sum_L_b;

                        }
                    }else{
                        if (loopN % 2 != 0){ //odd number term
                            //L_a = L_a * M.adjoint()* U_up2;
                            //L_b = M * U_down2.adjoint() * L_b;

                            //L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().eval().array(); //M.adjoint * U_up2
                            //L_b2_odd = M.array().rowwise()           * U_down2.adjoint().array(); //M * U_down2.adjoint

                            sum_L_a = sum_L_a * L_a2_odd;
                            sum_L_b = L_b2_odd * sum_L_b;

                        }else{ //even number term 
                            //L_a = L_a * M * U_down2;
                            //L_b = M.adjoint() * U_up2.adjoint() * L_b;

                            //L_a2_even = M.array().rowwise()           * U_down2.transpose().eval().array(); //M * U_down2
                            //L_b2_even = M.adjoint().array().rowwise() * U_up2.adjoint().array();     //M.adjoint * U_up2.adjoint

                            sum_L_a = sum_L_a * L_a2_even;
                            sum_L_b = L_b2_even * sum_L_b;

                        }
                        loopN++;
                    }
                }
                L_mat = sum_L_a * sum_L_b;
            }

            //////////////////////////////////////////////////////////////
            result[i] = MatrixXcd::Zero(1,1);
            result[i](0,0) = L_mat.trace() / doublec(bdim*1.0,0.0);
            tfree += deltat;
        }
    }else{

         //#pragma omp parallel for
         for (int i = 0; i < nstep; i++) {
            //double stime_a = MPI_Wtime();

            // Total U operator for each ms state
            // L = < J | (Udown)^dagger (Udup) | J >
            MatrixXcd UdownDagger(bdim,bdim);
            MatrixXcd Uup(bdim,bdim);

            // Temporary storage of time evolution operators
            //  for each interval of pulses
            MatrixXcd* Mdagger_UupTauDagger = new MatrixXcd[npulse+1]; 
            MatrixXcd* M_UdownTauDagger = new MatrixXcd[npulse+1]; 
            MatrixXcd* Mdagger_UupTau = new MatrixXcd[npulse+1]; 
            MatrixXcd* M_UdownTau = new MatrixXcd[npulse+1]; 

            // U operator related variables
            Eigen::VectorXcd es_up_tau;    Eigen::VectorXcd UupTau;
            Eigen::VectorXcd es_down_tau;  Eigen::VectorXcd UdownTau;

            // Assign the U operators for each sequence
            // and calculate total U operators
            for (int ipulse=0; ipulse<npulse+1; ipulse++){

                double tau = tfree * sequence[ipulse][2];
                int SameTauIndex = (int)sequence[ipulse][3];

                //printf("[%d] tau = %.10lf\n",ipulse,tau);

                // Calculate operators 
                if (SameTauIndex == ipulse){
                    // U = exp(-iHt)
                    // U = M(Ud)Mdagger
                    // Ud = exp(-i *Hd * t)
                    // ev_up_tau = Hd_up * tau
                    // ev_down_tau = Hd_down * tau
                    es_up_tau   =es_up_Val*tau;
                    es_down_tau =es_down_Val*tau;

                    // UupTau = Ud_up
                    // UdownTau = Ud_down
                    UupTau = es_up_tau.array().cos() 
                        - doublec(0.0,1.0)*es_up_tau.array().sin();
                    UdownTau = es_down_tau.array().cos() 
                            - doublec(0.0,1.0)*es_down_tau.array().sin();
    
                    //std::cout << "Eigenvalues\n" << std::endl;
                    //std::cout << es_up_tau << std::endl;
                    //std::cout << es_down_tau << std::endl;
                    //std::cout << "Uoperators\n" << std::endl;
                    //std::cout << UupTau << std::endl;
                    //std::cout << UdownTau << std::endl;

                    //Mdagger_UupTauDagger = M_dagger * UupTauDagger
                    //M_UdownTauDagger = M * UdownTauipulse,Dagger
                    //Mdagger_UupTau = M_dagger * UupTau
                    //M_UdownTau = M * UdownTau
                    Mdagger_UupTauDagger[ipulse] = (M.adjoint().array().rowwise() * UupTau.adjoint().array()).matrix();
                    M_UdownTauDagger[ipulse]     = (M.array().rowwise()           * UdownTau.adjoint().array()).matrix();
                    Mdagger_UupTau[ipulse]       = (M.adjoint().array().rowwise() * UupTau.transpose().eval().array()).matrix();
                    M_UdownTau[ipulse]           = (M.array().rowwise()           * UdownTau.transpose().eval().array()).matrix();
                }
                else{
                    Mdagger_UupTauDagger[ipulse] = Mdagger_UupTauDagger[SameTauIndex];
                    M_UdownTauDagger[ipulse]     = M_UdownTauDagger[SameTauIndex]    ;
                    Mdagger_UupTau[ipulse]       = Mdagger_UupTau[SameTauIndex]      ;
                    M_UdownTau[ipulse]           = M_UdownTau[SameTauIndex]          ;
                }

                if (ipulse==0){
                    // M_down * U_down_dagger 
                    UdownDagger = (M_down.array().rowwise() 
                                * UdownTau.adjoint().array()).matrix();
                    // M_down_dagger * M_up * U_up * M_up_dagger 
                    Uup         = (M.adjoint().array().rowwise() * UupTau.transpose().eval().array()).matrix() 
                                * (M_up.adjoint());

                    //std::cout << "M_down * U_down_dagger" << std::endl;
                    //std::cout << UdownDagger << std::endl;
                    //std::cout << "M * U_up * M_up_dagger" << std::endl;
                    //std::cout << Uup << std::endl;
                }
                else if (ipulse!=0 && ipulse%2==1){
                    // odd
                    UdownDagger = UdownDagger            * Mdagger_UupTauDagger[ipulse];
                    Uup         = M_UdownTau[ipulse]     * Uup;
                    //std::cout << "M_dagger * U_up_dagger" << std::endl;
                    //std::cout << Mdagger_UupTauDagger[ipulse] << std::endl;
                    //std::cout << "M * U_down" << std::endl;
                    //std::cout << M_UdownTau[ipulse] << std::endl;
                }
                else if (ipulse!=0 && ipulse%2==0){
                    //even
                    UdownDagger = UdownDagger            * M_UdownTauDagger[ipulse];
                    Uup         = Mdagger_UupTau[ipulse] * Uup;
                }
                else{printf("Error calsingle.cpp\n");exit(1);}
            }

            //std::cout << "UdownDagger" << UdownDagger << std::endl;
            //std::cout << "Uup" << Uup << std::endl;
            //std::cout << psi0 << std::endl;
            result[i] = (psi0.adjoint() * UdownDagger * Uup * psi0);

            //printf("%.10lf  %.10lf %+.10lf j \n",tfree,result[i](0,0).real(),result[i](0,0).imag());
            tfree += deltat;                      

            delete[] Mdagger_UupTauDagger;
            delete[] M_UdownTauDagger;
            delete[] Mdagger_UupTau;
            delete[] M_UdownTau;

            //etime = MPI_Wtime();
            //printf("          Wall time step4-%d# = %.5f s\n",i, stime_a - etime);
        }
    }

    //etime = MPI_Wtime();
    //printf("          Wall time step done = %.5f s\n",stime - etime);

    // free
    for (int i=0; i<nspin; i++){
        delete[] bsigmas[i];
    }
    delete[] bsigmas;

    delete[] Hqb;

    return result;
}

