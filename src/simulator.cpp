#include "../include/simulator.h"
#include "../include/hamiltonian.h"
#include "../include/memory.h"
#include "../include/reader.h"

void calculate(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, Cluster* cls, Output* op){

    char message[500];

    if (rank==0){
        printLineSection();
        printTitle("MPI distribution for the clusters");
    }
    ////////////////////////////////
    // MPI : distribute the clusters
    ////////////////////////////////
    int order = Cluster_getOrder(cls);
    int*** clusters = Cluster_getClusinfo(cls);
    int*** localclusters = MPI_getLocalClusters(order,clusters);
    // int*** localclusters = clusters; // No MPI
    // reportClusinfo(localclusters,order);
    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////
    // Check is GCCE or CCE
    ////////////////////////////////
    int iter_0th = Cluster_getClusinfo_iter(cls,0,0);
    char* method = Config_getMethod(cnf);
    char* quantity = Config_getQuantity(cnf);

    bool isGCCE = false;
    if (iter_0th != 0){isGCCE = true;}

    if (rank==0){ 
        printLineSection();
        if (isGCCE){
            if (strcasecmp(method,"gcce")!=0){
                fprintf(stderr,"Error : Method is not gCCE\n");
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
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0){
        printLineSection();
        printTitle("Main calculation loop");

        // Print the information
        if (nstate == 0){
            printLine();
            printMessage("Ensemble calculation");
            printLine();
            printf("\n");
        }else{
            printLine();
            printMessage("Single calculation");
            printLine();
            printf("\n");
        }

    }

    for (int istate=0; istate<=nstate; istate++){

        ////////////////////////////////
        // Update QubitArray, BathArray 
        ////////////////////////////////
        if (rank==0 && nstate!=0 && istate!=0){
            sprintf(message,"Iteration %d : Update QubitArray, BathArray",istate);
            printMessage(message);
            printLine();
        }
        
        // Ensemble app. doesn't generate random states

        // Single app. generates random states
        if (istate>0){

            // Randomize bath states
            // BathArray_setBathStatesRandom(ba);
            setBathStates(ba,cnf,istate); // Read or random set
            if (rank==0){
                BathArray_reportBath_states(ba);
            }

            // Set Subbath states
            if (DefectArray_getNdefect(dfa)>0){
                setSubbathStates(dfa,ba,cnf,istate); // Read or random set
                if (rank==0){
                    DefectArray_reportSubbath_states(dfa);
                }
            }
            

            // Calculate bath disorders
            BathArray_setBathDisorders(ba);
            if (rank==0){
                BathArray_reportBath_disorders(ba);
            }
            
            // Calculate overhauser fields of qubits
            bool isOvh = QubitArray_getOverhaus(qa);
            if (isOvh){
                for (int iq=0; iq<nqubit; iq++){
                    double overhaus = BathArray_getOverhaus(ba,iq);
                    QubitArray_setQubit_i_overhaus(qa,overhaus,iq);
                }
            }
            if (rank==0){
                QubitArray_reportQubit_overhaus(qa);
            }

            if (DefectArray_getNdefect(dfa)>0){
                // Update BathArray from DefectArray
                updateMainSpins_fromDefectArray(dfa,ba);
                updateDisorder_main_sub(dfa,ba);
                updateDisorder_sub_sub(dfa);
                updateOverhaus_qubit_sub(dfa,qa);
                if (rank==0){
                    BathArray_reportBath_disorders(ba);
                    DefectArray_reportSubbath_disorders(dfa);
                    QubitArray_reportQubit_overhaus(qa);
                }
            }
            

        }

        // Set the initial state and psia, psib
        if (istate==0 && istate < nstate){
            ; // Single, no cal. for istate=0
        }else{
            
            // Calculate the qubit Hamiltonian
            float* bfield = Config_getBfield(cnf);
            int* alphaidx = QubitArray_get_alphaidx(qa);
            int* betaidx = QubitArray_get_betaidx(qa);
            if (alphaidx!=NULL && betaidx!=NULL){
                QubitArray_setPsiaPsib_fromIdx(qa,bfield);
            }else if (alphaidx==NULL && betaidx==NULL){
                QubitArray_setPsiaPsib_fromQubit(qa);
            }else{
                fprintf(stderr,"Error : calculate,, Invalid qubit index (cannot set the psia,psib)\n");
                exit(1);
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

            MatrixXcd* result_0th; //= new MatrixXcd[nstep];
            MatrixXcd* result_nth; //= new MatrixXcd[nstep];

            MatrixXcd* result_0th_inv = new MatrixXcd[nstep];

            int dim = 0;
            if (isGCCE && strcasecmp(quantity,"dm")==0){
                dim = QubitArray_dim(qa);
            }else{
                dim = 1;
            }          

            for (int istep=0; istep<nstep; istep++){
                result_wD[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_nD[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                // result_0th[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_0th_inv[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                // result_nth[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
            }

            ////////////////////////////////
            // Zero-th order 
            double steptime_sta = MPI_Wtime();
            double steptime_end = 0.0;

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
            }else{
                result_0th = new MatrixXcd[nstep];
                for (int istep=0; istep<nstep; istep++){
                    result_0th[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
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
                for (int ic=1; ic<=ncluster; ic++){

                    // Get only spin index part from cluster
                    int cluster[n];
                    copyInt1dPart(cluster,localclusters[n][ic],1,n);

                    // Get the power
                    int power = localclusters[n][ic][0];

                    // Create the bath array for the cluster
                    BathArray* ba_cluster = createBathArray(cluster,n,ba,dfa,nqubit);
                    // BathArray_report(ba_cluster);

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
                        Output_save_all(op,result_nth,nstep,deltat,cluster,n,istate);
                    }

                    if (rank==0){
                        steptime_end = MPI_Wtime();
                        if (ic % 142){
                            sprintf(message,"At rank 0, %d cluster (%d%%) is computed at %.2f s",ic, int(ic*100/ncluster), steptime_end-steptime_sta);
                            printMessage(message);

                        }else if (ic == ncluster){
                            sprintf(message,"At rank 0, %d cluster (%d%%) is computed at %.2f s",ic, int(ic*100/ncluster), steptime_end-steptime_sta);
                            printMessage(message);
                        }
                    }
                    delete[] result_nth;
                }

                if (rank==0){
                    steptime_end = MPI_Wtime();
                    // make message
                    printf("\n");
                    sprintf(message,"# Wall time : %.2f[s]\n",steptime_end-steptime_sta);
                    printMessage(message);
                    printLine(); printf("\n");
                }
            }


            MPI_Barrier(MPI_COMM_WORLD);
            if (rank==0){
                printSubTitle("Product local results...");
                // zeroth cluster  L_C_(k=0) * PI_C_(k=1){ L_C_(k=1) } * PI_C_(k=2){ L_C_(k=2) }  .. 
                for (int istep=0; istep<nstep; istep++){
                    result_nD[istep] = mulMatrixXcdElementWise(result_nD[istep], result_0th[istep]);
                    result_wD[istep] = mulMatrixXcdElementWise(result_wD[istep], result_0th[istep]);
                }           
            }
    
            delete[] result_0th;
            delete[] result_0th_inv;
            

            MPI_Barrier(MPI_COMM_WORLD);

            ////////////////////////////////
            // Gather the results (MPI)
            ////////////////////////////////
            if (rank==0){
                printSubTitle("Gather the results (MPI)...");
            }
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
            }

            MPI_Barrier(MPI_COMM_WORLD);

            if (rank !=0){
                delete[] result_wD_all;
                delete[] result_nD_all;
            }

            delete[] result_wD;
            delete[] result_nD;

            ////////////////////////////////
            // Normalize the results
            ////////////////////////////////
            if (rank == 0){

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
            
            MPI_Barrier(MPI_COMM_WORLD);
            ////////////////////////////////
            // Save the results
            ////////////////////////////////
            if (rank==0){

                printSubTitle("Save the results...");
                
                if (strcasecmp(savemode,"normal")==0){
                    Output_save(op,result_wD_all,result_nD_all,nstep,deltat,istate);
                }else if (strcasecmp(savemode,"all")==0){
                    Output_save(op,result_wD_all,result_nD_all,nstep,deltat,istate);
                }else if (strcasecmp(savemode,"info")==0){
                    char* bathfile = Config_getBathfiles_i(cnf,0);
                    MatrixXcd Atensor = BathArray_getBath_i_hypf_j(ba,0,0);
                    double Azx = Atensor(2,0).real();
                    double Azz = Atensor(2,2).real();
                    Output_save_info(op,result_wD_all,result_nD_all,nstep,deltat,istate, Azx, Azz, bathfile);
                }
            }

            ////////////////////////////////
            // Free the results
            ////////////////////////////////
            if (rank==0){
                delete[] result_wD_all;
                delete[] result_nD_all;
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

}

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, DefectArray* dfa, int nqubit){

    BathArray* ba_cluster = BathArray_init();

    // Allocate the bath spins 
    BathArray_setNspin(ba_cluster,nspin);
    BathArray_allocBath(ba_cluster,nspin);

    int bdim = BathArray_dim(ba);

    // Create the bath array for the cluster
    for (int ic=0; ic<nspin; ic++){
        int ibs = cluster[ic];
        BathSpin* bs  = BathArray_getBath_i(ba,ibs);
        BathArray_setBath_i(ba_cluster, bs, ic, nqubit);
    }

    // // Set the properties
    // int nspecies = BathArray_getProp_nspecies(ba);
    // BathArray_setProp_nspecies(ba_cluster,nspecies);
    // BathArray_allocProp(ba_cluster);
    // char** names = BathArray_getProp_names(ba);
    // double* gyros = BathArray_getProp_gyros(ba);
    // float* spins = BathArray_getProp_spins(ba);
    // for (int i=0; i<nspecies; i++){
    //     BathArray_setProp_names_i(ba_cluster,names[i],i);
    //     BathArray_setProp_gyros_i(ba_cluster,gyros[i],i);
    //     BathArray_setProp_spins_i(ba_cluster,spins[i],i);       
    // }
    
        
    // Check if spins in the cluster is electronic defect
    int addedspin = 0;
    for (int i=0; i<nspin; i++){
        int ibs = cluster[i];
        char* name = BathArray_getBath_i_name(ba,ibs);
        int idf = DefectArray_findDefectIndex(dfa,name);

        if (idf<0){
            ; // doesn't have additional spin information
        }else{
            // has additional spin information
            int naddspin = DefectArray_getDefect_idf_naddspin(dfa,idf);
            for (int isp=0; isp<naddspin; isp++){
                
                // Reallocate the bath array
                int nspin_old = BathArray_getNspin(ba_cluster);
                int nspin_new = nspin_old + 1;
                BathArray_setNspin(ba_cluster,nspin_new);
                BathArray_reallocBath(ba_cluster,nspin_old,nspin_new,nqubit);

                // Add the additional spin to BathArray_for_cluster (ba_cluster)
                BathSpin* bs_sub = DefectArray_getSubbath_i_isp(dfa,ibs,isp);
                BathArray_setBath_i(ba_cluster,bs_sub,nspin_new-1,nqubit);

                // Add the additional spin to the cluster
                addedspin++;
            }
        }
    }

    if (rank==0 && addedspin>0 && verbosity){
        printStructElementInt("Added spin",addedspin);
    }

    return ba_cluster;
}
