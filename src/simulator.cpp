#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include "../include/memory.h"

void calculate(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Cluster* cls, Output* op){


    ////////////////////////////////
    // MPI : distribute the clusters
    ////////////////////////////////
    int order = Cluster_getOrder(cls);
    int*** clusters = Cluster_getClusinfo(cls);

    // int*** localclusters = Cluster_getClusinfo_MPI(order,clusters); // MPI
    int*** localclusters = clusters; // No MPI

    ////////////////////////////////
    // Check is GCCE or CCE
    ////////////////////////////////
    int iter_0th = Cluster_getClusinfo_iter(cls,0,0);
    char* method = Config_getMethod(cnf);
    bool isGCCE = false;
    if (iter_0th != 0){isGCCE = true;}

    if (rank==0){
        printLineSection();
        if (isGCCE){
            if (strcasecmp(method,"gcce")!=0){
                fprintf(stderr,"Error : Method is not GCCE\n");
                exit(1);
            }
            printTitle("gCCE calculation");
        }else{  
            if (strcasecmp(method,"cce")!=0){
                fprintf(stderr,"Error : Method is not CCE\n");
                exit(1);
            }
            printTitle("CCE calculation");
        }
    }

    ////////////////////////////////
    // Save mode
    ////////////////////////////////
    char* savemode = Output_getSavemode(op);

    ////////////////////////////////
    // Main calculation loop
    ////////////////////////////////

    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    int nstate = Config_getNstate(cnf);
    int nstep = Config_getNstep(cnf);
    float deltat = Config_getDeltat(cnf);

    for (int istate=0; istate<=nstate; istate++){

        ////////////////////////////////
        // Update QubitArray, BathArray 
        ////////////////////////////////

        // Ensemble app. doesn't generate random states

        // Single app. generates random states
        if (istate>0){
    
            // set Bath random states
            BathArray_setBathStatesRandom(ba);
            //! Defect_setDefectSpinStateRandom(defect); // if defect is used

            // Set Bath disorders
            BathArray_setBathDisorders(ba);
            //! Disorder due to defect spins

            // Set OverHauser fields
            for (int iq=0; iq<nqubit; iq++){
                double overhaus = BathArray_getOverhaus(ba,iq);
                //!  BathArray_getOverhaus_sub(ba,iq); if apprx. is used 
                //! overhaus due to defect spins
                QubitArray_setQubit_i_overhaus(qa,overhaus,iq);
            }


        }

        // Set the initial state and psia, psib
        if (istate==0 && istate < nstate){
            ; // Single, no cal. for istate=0
        }else{

            // Calculate the qubit Hamiltonian
            float* bfield = Config_getBfield(cnf);
            int* alphaidx = QubitArray_get_alphaidx(qa);
            if (alphaidx!=NULL){
                QubitArray_setPsiaPsib_fromIdx(qa,bfield);
            }

            // Qubit initial state
            MatrixXcd psi0 = QubitArray_getPsi0(qa);
            if (psi0.rows()==0 && psi0.cols()==0){
                QubitArray_setPsi0_fromPsiaPsib(qa);
            }

            ////////////////////////////////
            //  Calculate the coherence
            ////////////////////////////////

            ////////////////////////////////
            // Initialize the results
            MatrixXcd* result_wD = new MatrixXcd[nstep];
            MatrixXcd* result_nD = new MatrixXcd[nstep];

            MatrixXcd* result_0th = new MatrixXcd[nstep];
            MatrixXcd* result_nth = new MatrixXcd[nstep];

            MatrixXcd* result_0th_inv = new MatrixXcd[nstep];

            int dim = 0;
            if (isGCCE){
                dim = QubitArray_dim(qa);
            }else{
                dim = 1;
            }          

            for (int istep=0; istep<nstep; istep++){
                result_wD[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_nD[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_0th[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_0th_inv[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_nth[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
            }

            ////////////////////////////////
            // Zero-th order 
            double steptime_sta = MPI_Wtime();
            double steptime_end = 0.0;
            char message[500];

            if (isGCCE){ // gcce
                // Zero-th order
                if (rank==0){
                    printLine();
                    printSubTitle("Calculate 0-th cluster...\n");
                }
                result_0th = calCoherenceGcce(qa,NULL,cnf,pls);
                for (int istep=0; istep<nstep; istep++){
                    result_0th_inv[istep] = powMatrixXcdElementWise(result_0th[istep],-1);
                }
                if (rank==0){
                    steptime_end = MPI_Wtime();
                    // make message
                    sprintf(message,"Wall time : %.2f[s]\n",steptime_end-steptime_sta);
                    printSubTitle(message);
                    printLine(); printf("\n");
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);

            ////////////////////////////////
            // Higher orders
            for (int n=1; n<=order; n++){

                // Message
                if (rank==0){
                    sprintf(message,"Calculate %d-th cluster...\n",n);
                    printLine();
                    printSubTitle(message);
                }

                // Number of clusters
                int ncluster = localclusters[n][0][0]-1;

                // Create BathArray for the cluster (defect spins can be added)
                for (int ic=1; ic<ncluster+1; ic++){

                    // Get only spin index part from cluster
                    int cluster[n];
                    copyInt1dPart(cluster,localclusters[n][ic],1,n);

                    // Get the power
                    int power = localclusters[n][ic][0];

                    // Create the bath array for the cluster
                    BathArray* ba_cluster = createBathArray(cluster,nspin,ba, nqubit);

                    // Calculate the coherence
                    if (isGCCE){ // gcce
                        result_nth = calCoherenceGcce(qa, ba_cluster, cnf, pls);
                    }else{ // cce
                        result_nth = calCoherenceCce(qa, ba_cluster, cnf, pls);
                    }

                    // Update the total result
                    // zeroth cluster : result_nth = L_C_(k=n) * ( 1 / L_C_(k=0) )
                    // power : result_nth = ( L_C_(k=n) / L_C_(k=0) ) ^ power
                    // total result : result = result * ( L_C_(k=n) / L_C_(k=0) ) ^ power
                    for (int istep=0; istep<nstep; istep++){
                        result_nth[istep] = mulMatrixXcdElementWise(result_nth[istep], result_0th_inv[istep]);
                        result_nth[istep] = powMatrixXcdElementWise(result_nth[istep], power);

                        if (power>0){
                            result_nD[istep] = mulMatrixXcdElementWise(result_nD[istep], result_nth[istep]);
                        }

                        result_wD[istep] = mulMatrixXcdElementWise(result_wD[istep], result_nth[istep]);
                    }

                    if (strcasecmp(savemode,"all")==0){
                        Output_setOutfile_all(op,cluster,n,istate);
                        Output_save_all(op,result_nth,nstep,deltat);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);

                    if (rank==0){
                        steptime_end = MPI_Wtime();
                        if (ic % 142){
                            sprintf(message," At rank 0, %d cluster (%d%%) is computed at %.2f s\n",ic, int(ic*100/ncluster), steptime_end-steptime_sta);
                            printSubTitle(message);
                        }else if (ic == ncluster){
                            sprintf(message," At rank 0, %d cluster (%d%%) is computed at %.2f s\n",ic, int(ic*100/ncluster), steptime_end-steptime_sta);
                            printSubTitle(message);
                        }
                    }
                }

                if (rank==0){
                    steptime_end = MPI_Wtime();
                    // make message
                    sprintf(message,"Wall time : %.2f[s]\n",steptime_end-steptime_sta);
                    printSubTitle(message);
                    printLine(); printf("\n");
                }
            }

            if (rank==0){
                // zeroth cluster  L_C_(k=0) * PI_C_(k=1){ L_C_(k=1) } * PI_C_(k=2){ L_C_(k=2) }  .. 
                for (int istep=0; istep<nstep; istep++){
                    result_nD[istep] = mulMatrixXcdElementWise(result_nD[istep], result_0th[istep]);
                    result_wD[istep] = mulMatrixXcdElementWise(result_wD[istep], result_0th[istep]);
                }           
            }
    
            delete[] result_0th;
            delete[] result_0th_inv;
            delete[] result_nth;

            MPI_Barrier(MPI_COMM_WORLD);

            ////////////////////////////////
            // Gather the results (MPI)
            ////////////////////////////////

            // Initialize the result_all
            MatrixXcd* result_wD_all = new MatrixXcd[nstep];
            MatrixXcd* result_nD_all = new MatrixXcd[nstep];
            for (int istep=0; istep<nstep; istep++){
                result_wD_all[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_nD_all[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // Reduce the results
            int err;
            for (int istep=0; istep<nstep; istep++){
                err = MPI_Reduce(result_wD[istep].data(),result_wD_all[istep].data(),dim*dim,MPI_DOUBLE_COMPLEX,MPI_PROD,0,MPI_COMM_WORLD);
                err = MPI_Reduce(result_nD[istep].data(),result_nD_all[istep].data(),dim*dim,MPI_DOUBLE_COMPLEX,MPI_PROD,0,MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
            }

            if (rank !=0){
                delete[] result_wD_all;
                delete[] result_nD_all;
            }

            delete[] result_wD;
            delete[] result_nD;

            MPI_Barrier(MPI_COMM_WORLD);

            ////////////////////////////////
            // Normalize the results
            ////////////////////////////////
            if (rank == 0){

                // quantity
                char* quantity = Config_getQuantity(cnf);

                // get normalization factor
                MatrixXcd normalization;
                MatrixXcd normalization_inv;
                MatrixXcd rho0 = QubitArray_Rho0(qa);

                if ( strcasecmp(quantity,"coherence")==0 && isGCCE){
                    MatrixXcd psia = QubitArray_getPsia(qa);
                    MatrixXcd psib = QubitArray_getPsib(qa);
                    normalization = psia.adjoint() * rho0 * psib;
                    normalization_inv = powMatrixXcdElementWise(normalization,-1); 
                }
                else if ( strcasecmp(quantity,"dm")==0 && isGCCE ){
                    normalization = rho0;
                    normalization_inv = powMatrixXcdElementWise(normalization,-1); 
                }else{
                    normalization = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                    normalization_inv = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                }
                
                // do normalization : Ltot / normalization
                for (int i=0; i<nstep; i++){
                    result_wD_all[i] = mulMatrixXcdElementWise(result_wD_all[i], normalization_inv);
                    result_nD_all[i] = mulMatrixXcdElementWise(result_nD_all[i], normalization_inv);
                }
            }

            ////////////////////////////////
            // Save the results
            ////////////////////////////////
            if (rank==0){
                if (strcasecmp(savemode,"normal")==0){
                    Output_setOutfile_normal(op,istate);
                    Output_save(op,result_wD_all,result_nD_all,nstep,deltat);
                }else if (strcasecmp(savemode,"all")==0){
                    Output_setOutfile_normal(op,istate);
                    Output_save(op,result_wD_all,result_nD_all,nstep,deltat);
                }else if (strcasecmp(savemode,"info")==0){
                    char* bathfile = Config_getBathfiles_i(cnf,0);
                    MatrixXcd Atensor = BathArray_getBath_i_hypf_j(ba,0,0);
                    double Azx = Atensor(2,0).real();
                    double Azz = Atensor(2,2).real();
                    Output_setOutfile_info(op,istate);
                    Output_save_info(op,result_wD_all,result_nD_all,nstep,deltat, Azx, Azz, bathfile);
                }
            }

            ////////////////////////////////
            // Free the results
            ////////////////////////////////
            if (rank==0){
                delete[] result_wD_all;
                delete[] result_nD_all;
            }
        }
    }

}

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, int nqubit){
    
    BathArray* ba_cluster = BathArray_init();

    // Allocate the bath spins 
    BathArray_setNspin(ba_cluster,nspin);
    BathArray_allocBath(ba_cluster,nqubit);
    
    int bdim = BathArray_dim(ba);

    // Create the bath array for the cluster
    int ispin;
    BathSpin* bs_i;
    for (int ic=0; ic<nspin; ic++){
        ispin = cluster[ic];
        bs_i  = BathArray_getBath_i(ba,ispin);

        BathArray_setBath_i(ba_cluster, bs_i, ic, nqubit);
    }

    // Set the properties
    int nspecies = BathArray_getProp_nspecies(ba);
    BathArray_setProp_nspecies(ba_cluster,nspecies);

    BathArray_allocProp(ba_cluster);
    char** names = BathArray_getProp_names(ba);
    double* gyros = BathArray_getProp_gyros(ba);
    float* spins = BathArray_getProp_spins(ba);
    for (int i=0; i<nspecies; i++){
        BathArray_setProp_names_i(ba_cluster,names[i],i);
        BathArray_setProp_gyros_i(ba_cluster,gyros[i],i);
        BathArray_setProp_spins_i(ba_cluster,spins[i],i);       
    }
    return ba_cluster;
}

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

    // Pauli matrices
    MatrixXcd** qsigmas = QubitArray_PauliOperators(qa);
    MatrixXcd** bsigmas = BathArray_PauliOperators(ba);

    // Qubit Hamiltonian
    MatrixXcd Hq = HamilQubit(qa,ba,qsigmas,cnf);
    MatrixXcd Hq_expand = kron(Hq,MatrixXcd::Identity(bdim,bdim));

    // Bath Hamiltonian
    if (nspin > 0){
        MatrixXcd Hb = HamilBath(ba,bsigmas,cnf);
        MatrixXcd Hb_expand = kron(MatrixXcd::Identity(qdim,qdim),Hb);

        // Qubit - Bath Interaction Hamiltonian
        MatrixXcd Hqb = HamilQubitBath(qa,ba,qsigmas,bsigmas,cnf);

        // Total Hamiltonian
        Htot = Hq_expand + Hb_expand + Hqb;

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

    for (int i=0; i<nspin; i++){
        delete[] bsigmas[i];
    }
    delete[] bsigmas;

    return result;
}

MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pls, double tfree){

    int npulse = Pulse_getNpulse(pls);
    double** sequence = Pulse_getSequence(pls);

    // Propagator for pulse
    MatrixXcd* sigmas_general = QubitArray_PauliOperator_fromPsiaPsib(qa);
    MatrixXcd UpulseExponent = ((-1.0) * doublec(0.0,1.0) * sigmas_general[1] * M_PI/2.0);
    MatrixXcd Upulse = UpulseExponent.exp();

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

MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls){

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

    // Qubit - Bath Interaction Hamiltonian Hqb[0] = alpha * Hqb, Hqb[1] = beta * Hqb
    MatrixXcd* Hqb = HamilQubitBathSecularApp(qa,ba,bsigmas,cnf);

    // Total Hamiltonian
    Halpha = Hb + Hqb[0];
    Hbeta = Hb + Hqb[1];

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

    ////////////////////////////////
    // Initial state of bath
    ////////////////////////////////

    // Check if this is the ensemble calculation
    int nstate = Config_getNstate(cnf);
    bool isEnsemble = false;
    if (nstate == 0){isEnsemble = true;}

    MatrixXcd psi0 = BathArray_Psi0(ba);

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
                L_mat     =  (M.array().rowwise() * U_down.transpose().array() ).matrix() 
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
                L_a_even = M.array().rowwise()           * U_down.transpose().array();   //M * U_down
                L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array();       //M.adjoint * U_up.adjoint

                if(npulse % 2 != 0){
                    L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().array();  //M.adjoint() * U_up
                    L_b_odd = M.array().rowwise()           * U_down.adjoint().array();  //M * U_down.adjoint
                }

                switch (npulse){
                default :
                case 3:
                    L_a2_even = M.array().rowwise()           * U_down2.transpose().array();  //M * U_down2
                    L_b2_even = M.adjoint().array().rowwise() * U_up2.adjoint().array();      //M.adjoint * U_up2.adjoint
                case 2:
                    L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().array();     //M.adjoint * U_up2
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

                            //L_a_even = M.array().rowwise()           * U_down.transpose().array(); //M * U_down
                            //L_b_even = M.adjoint().array().rowwise() * U_up.adjoint().array(); //M.adjoint * U_up.adjoint

                            sum_L_a = sum_L_a * L_a_even;
                            sum_L_b = L_b_even * sum_L_b;

                        }else{ //odd number pulse
                            //L_a = L_a * M.adjoint() * U_up;
                            //L_b = M * U_down.adjoint() * L_b;

                            //L_a_odd = M.adjoint().array().rowwise() * U_up.transpose().array(); //M.adjoint() * U_up
                            //L_b_odd = M.array().rowwise()           * U_down.adjoint().array(); //M * U_down.adjoint
                            
                            sum_L_a = sum_L_a * L_a_odd;
                            sum_L_b = L_b_odd * sum_L_b;

                        }
                    }else{
                        if (loopN % 2 != 0){ //odd number term
                            //L_a = L_a * M.adjoint()* U_up2;
                            //L_b = M * U_down2.adjoint() * L_b;

                            //L_a2_odd = M.adjoint().array().rowwise() * U_up2.transpose().array(); //M.adjoint * U_up2
                            //L_b2_odd = M.array().rowwise()           * U_down2.adjoint().array(); //M * U_down2.adjoint

                            sum_L_a = sum_L_a * L_a2_odd;
                            sum_L_b = L_b2_odd * sum_L_b;

                        }else{ //even number term 
                            //L_a = L_a * M * U_down2;
                            //L_b = M.adjoint() * U_up2.adjoint() * L_b;

                            //L_a2_even = M.array().rowwise()           * U_down2.transpose().array(); //M * U_down2
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

         for (int i = 0; i < nstep; i++) {

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
                    //M_UdownTauDagger = M * UdownTauDagger
                    //Mdagger_UupTau = M_dagger * UupTau
                    //M_UdownTau = M * UdownTau
                    Mdagger_UupTauDagger[ipulse] = (M.adjoint().array().rowwise() * UupTau.adjoint().array()).matrix();
                    M_UdownTauDagger[ipulse]     = (M.array().rowwise()           * UdownTau.adjoint().array()).matrix();
                    Mdagger_UupTau[ipulse]       = (M.adjoint().array().rowwise() * UupTau.transpose().array()).matrix();
                    M_UdownTau[ipulse]           = (M.array().rowwise()           * UdownTau.transpose().array()).matrix();
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
                    Uup         = (M.adjoint().array().rowwise() * UupTau.transpose().array()).matrix() 
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

            //std::cout << UdownDagger * Uup << std::endl;

            result[i] = (psi0.adjoint() * UdownDagger * Uup * psi0);
        
            tfree += deltat;                      

            delete[] Mdagger_UupTauDagger;
            delete[] M_UdownTauDagger;
            delete[] Mdagger_UupTau;
            delete[] M_UdownTau;
        }
    }

    // free
    for (int i=0; i<nspin; i++){
        delete[] bsigmas[i];
    }
    delete[] bsigmas;

    delete[] Hqb;

    return result;
}


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

        Hqb += Hqbij;
    }

    //! hyperfine mediated term
    //! if (hfmedi){

    //! }

    // get qubit sublevel
    MatrixXcd* qsigmas_general = QubitArray_PauliOperator_fromPsiaPsib(qa);
    MatrixXcd qsz = qsigmas_general[3];
    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);
    double ms_a = (psia.adjoint() * qsz * psia)(0,0).real();
    double ms_b = (psib.adjoint() * qsz * psib)(0,0).real();

    // Hqb
    MatrixXcd* Hqb_PsiaPsib = new MatrixXcd[2];
    Hqb_PsiaPsib[0] = ms_a * Hqb;
    Hqb_PsiaPsib[1] = ms_b * Hqb;

    return Hqb_PsiaPsib;
}