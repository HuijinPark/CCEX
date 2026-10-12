#include "../../include/calculation/gpu/gpu_batch_memory.h"
#include "../../include/calculation/simulator.h"
#include "../../include/calculation/hamiltonian.h"
#include "../../include/common/memory.h"
#include "../../include/input/reader.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/calculation/scaled_product.h"
#include "../../include/auto_convergence/rdip_ladder.h"
#include "../../include/clustering/bundle/bundle_work.h"
#include "../../include/clustering/bundle/bundle_connection_check.h"
#include "../../include/clustering/bundle/bundle_contributions.h"
#include <cstddef>

static void reduceScaledProduct(void* input,void* output,int* count,MPI_Datatype*){
    auto* a=static_cast<ScaledComplex*>(input);
    auto* b=static_cast<ScaledComplex*>(output);
    for(int i=0;i<*count;++i)b[i]=a[i]*b[i];
}

void calculate(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, Cluster* cls, Output* op, int*** localclusters){

    char message[500];

    ////////////////////////////////
    // Get parameters
    ////////////////////////////////

    // Method-related parameters
    char* method = Config_getMethod(cnf);
    char* quantity = Config_getQuantity(cnf);
    int nstate = Config_getNstate(cnf);

    // Calculation condition parameters
    int nstep = Config_getNstep(cnf);
    double deltat = Config_getDeltat(cnf);

    // The number of Bath and Qubit spins
    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    
    // Output parameters
    char* savemode = Output_getSavemode(op);

    // Power of the 0th cluster
    int order = Cluster_getOrder(cls);
    int iter_0th = Cluster_getClusinfo_iter(cls,0,0);
    bool isGCCE = false;
    
    ////////////////////////////////
    // Print the calculation method
    ////////////////////////////////

    // gCCE or conventional CCE
    if (strcasecmp(method,"gcce")==0 && iter_0th==0){
        // gCCE method requires the 0th cluster
        fprintf(stderr,"Error : calculate,, gCCE method requires the 0th cluster\n");
        exit(1);
    }else if (strcasecmp(method,"gcce")==0 && iter_0th>0){
        isGCCE = true;
    }

    // The Overhauser term reaches the Hamiltonian only through QubitArray_SingleHamil
    // -> HamilQubit (simulator_hamiltonian.cpp:4), and calCoherenceGcce is its only
    // caller. Under CCE the qubit is frozen into the classical labels alphams/betams, so
    // its own Hamiltonian never enters -- the field is still parsed, stored and printed,
    // which reads exactly like it took effect.
    if (!isGCCE && QubitArray_getOverhaus(qa) && rank==0){
        fprintf(stderr,"Warning: \"overhaus\" has no effect with method=%s.\n"
                       "         CCE projects the qubit out, so the qubit Hamiltonian --\n"
                       "         where the Overhauser field lives -- is never built.\n"
                       "         Use method=gCCE for mean-field.\n", method);
    }


    MPI_Barrier(MPI_COMM_WORLD);

    // profile: per-rank wall time of the parts of a state outside the cluster calls (GPU plan A8). Only measured with
    // "profile": true; printed once by rank 0, rank by rank, after the last state.
    enum {PH_PREP, PH_CREATE, PH_CLUSTER, PH_ASSEMBLE, PH_FREE, PH_WAIT, PH_REDUCE, PH_SAVE, PH_N};
    const bool phases=cnf->profile;
    double phase[PH_N]={0};
    int phaseStates=0;
    double mark=0;
    auto lap=[&](int p){ if(phases){ const double t=MPI_Wtime(); phase[p]+=t-mark; mark=t; } };

    rdipladder::Ladder ladder(cnf,cls,ba,nstep);     // off unless rdip_ladder is given

    ////////////////////////////////
    // Main calculation loop
    ////////////////////////////////
    DisorderCache disorderCache;                          // A13: this rank's couplings, kept over the states
    CceWorkspace cceWork;                                 // A14: the cCCE storage reused from cluster to cluster
    std::unique_ptr<MatrixXcd[]> cceResult(new MatrixXcd[Config_getNstep(cnf)]);
    Eigen::RowVectorXcd gcceAlphaRho;                      // A14: psia^dagger rho of the gCCE coherence, reused
    const bool checkDisorder = getenv("CCEX_CHECK_DISORDER") && strcmp(getenv("CCEX_CHECK_DISORDER"),"0");
    for (int istate=0; istate<=nstate; istate++){

        // wall time
        double time_start_step = MPI_Wtime();
        mark=time_start_step;

        if (rank==0 && nstate!=0 && istate!=0){
            sprintf(message,"%s calculation (%s) : Iteration # %d",method,"Single-sample approach",istate); 
            printTitle(message);
        }else if (rank==0 && nstate==0 && istate==0){
            sprintf(message,"%s calculation (%s)",method,"Ensemble approach"); 
            printTitle(message);
        }

        /////////////////////////////////////////////////////////
        // Update QubitArray, BathArray, DefectArray 
        /////////////////////////////////////////////////////////
        // Generate random states in single-sample approach
        if (istate>0){

            ////////////////////////////////
            // Randomize bath states
            ////////////////////////////////
            // Set Bath states
            if (rank==0){
                printSubTitle("Randomize main-bath states...");
            }
            const bool restoredCheckState=restoreAdaptiveCheckState(ba,dfa,istate);
            if(!restoredCheckState)setBathStates(ba,cnf,istate); // Read or random set
            if (rank==0){
                BathArray_reportBath_states(ba);
            }

            // Set Subbath states
            if (DefectArray_getNdefect(dfa)>0){
                if (rank==0){
                    printSubTitle("Randomize sub-bath states...");
                }
                if(!restoredCheckState)setSubbathStates(dfa,ba,cnf,istate); // Read or random set
                if (rank==0){
                    DefectArray_reportSubbath_states(dfa);
                }
            }
            
            ////////////////////////////////
            // Calculate bath disorders
            ////////////////////////////////

            // Set the bath disorders: the spins split over the ranks, the couplings kept per rank (improvements.md A13).
            // CCEX_CHECK_DISORDER=1 also runs the one-rank-does-all function and stops on any bit that differs.
            BathArray_setBathDisordersRanks(ba, cnf->vector_mean_field, disorderCache);
            if (checkDisorder){
                const int ns = BathArray_getNspin(ba);
                std::vector<double> got(3*size_t(ns));
                for (int i=0; i<ns; i++){
                    got[3*i] = BathArray_getBath_i_disorder(ba,i);
                    got[3*i+1] = ba->bath[i]->disorder_xy[0]; got[3*i+2] = ba->bath[i]->disorder_xy[1];
                }
                BathArray_setBathDisorders(ba, cnf->vector_mean_field);
                for (int i=0; i<ns; i++){
                    const double want[3] = {BathArray_getBath_i_disorder(ba,i), ba->bath[i]->disorder_xy[0], ba->bath[i]->disorder_xy[1]};
                    if (memcmp(want, &got[3*i], sizeof want) != 0){
                        fprintf(stderr,"Error: CCEX_CHECK_DISORDER: rank %d, state %d, spin %d: %.17g %.17g %.17g (split) vs %.17g %.17g %.17g\n",
                                rank,istate,i,got[3*i],got[3*i+1],got[3*i+2],want[0],want[1],want[2]);
                        ccexFatal();
                    }
                }
                if (rank==0) printf("      CCEX_CHECK_DISORDER: state %d, %d spins, the same bits\n",istate,ns);
            }
            
            // Calculate overhauser fields of qubits
            bool isOvh = QubitArray_getOverhaus(qa);
            if (isOvh){
                for (int iq=0; iq<nqubit; iq++){
                    double overhaus = BathArray_getOverhaus(ba,iq);
                    QubitArray_setQubit_i_overhaus(qa,overhaus,iq);
                }
            }

            if (DefectArray_getNdefect(dfa)>0){
                // Update BathArray from DefectArray
                updateDisorder_sub_sub(dfa);
                updateDisorder_main_sub(dfa,ba);
                if (isOvh){
                    updateOverhaus_qubit_sub(dfa,qa);
                }
            }

            if (rank==0){
                QubitArray_reportQubit_overhaus(qa);
                BathArray_reportBath_disorders(ba);
                if (DefectArray_getNdefect(dfa)>0){
                    DefectArray_reportSubbath_disorders(dfa);
                }
            }
        }


        if (!(nstate !=0 && istate==0)){
            // Calculating condition : 
            //  nstate !=0 and istate!=0 : single sample approach
            //  nstate !=0 and istate==0 : X
            //  nstate ==0 and istate==0 : ensemble approach
            //  nstate ==0 and istate!=0 : X



            /////////////////////////////////////////////////////////
            // Set initial state and alpha, beta of qubits
            /////////////////////////////////////////////////////////

            // Get psia, psib as adiabatic states
            MatrixXcd psia = QubitArray_getPsia(qa);
            MatrixXcd psib = QubitArray_getPsib(qa);

            float* bfield = Config_getBfield(cnf);
            int* alphaidx = QubitArray_get_alphaidx(qa);
            int* betaidx = QubitArray_get_betaidx(qa);

            if (psia.rows()==0 && psia.cols()==0 && psib.rows()==0 && psib.cols()==0){
                if (alphaidx!=NULL && betaidx!=NULL){
                    if (rank==0){
                        printSubTitle("Set psia, psib from qubit index...");
                    }
                    QubitArray_setPsiaPsib_fromIdx(qa,bfield);
                    QubitArray_checkIdxAgainstMs(qa);
                }else if (alphaidx==NULL && betaidx==NULL){
                    if (rank==0){
                        printSubTitle("Set psia, psib from qubit alpha, beta...");
                    }
                    QubitArray_setPsiaPsib_fromQubit(qa);
                }else{
                    fprintf(stderr,"Error : calculate,, Invalid qubit index (cannot set the psia,psib)\n");
                    exit(1);
                }
            }else{
                if (rank==0){
                    printSubTitle("Set psia, psib from input...\n");
                }
            }

            if (rank==0){
                QubitArray_reportPsiaPsib(qa);
            }
            
            // Qubit initial state
            MatrixXcd psi0 = QubitArray_getPsi0(qa);
            if (psi0.rows()==0 && psi0.cols()==0){
                if (rank==0){
                    printSubTitle("Set psi0 from psia, psib...");
                }
                QubitArray_setPsi0_fromPsiaPsib(qa);
            }else{
                if (rank==0){
                    printSubTitle("Set psi0 from input...");
                }
            }

            // Report
            if (rank==0){
                QubitArray_reportPsi0(qa);
            }

            /////////////////////////////////////////////////////////
            //  Calculate the coherence
            /////////////////////////////////////////////////////////

            if(adaptiveCheckState(qa,ba,dfa,cnf,pls,op,istate)) {
                MPI_Barrier(MPI_COMM_WORLD);
                continue;
            }
            if(adaptiveWork(cnf).mode=="probe") {
                // no barrier between probe states (A1): a rank goes on to its probes of the next state, the node
                // memory budget of adaptiveProbe keeps the large ones from all running at once
                adaptiveProbe(qa,ba,dfa,cnf,pls,op,istate);
                continue;
            }
            if(adaptiveWork(cnf).mode=="factors") {
                // One state's culprit factors (adaptive-regroup-stage3.md 1절). The states before it were prepared
                // all the same, so its bath state and mean fields are the run's; after it nothing is wanted.
                if(istate!=adaptiveWork(cnf).state)continue;
                adaptiveFactors(qa,ba,dfa,cnf,pls,op,istate);
                break;
            }

            ////////////////////////////////
            // Initialize the results
            ClusterContributions contributions(cnf,cls,ba);
            MatrixXcd* result_wD = new MatrixXcd[nstep];
            MatrixXcd* result_nD = new MatrixXcd[nstep];

            MatrixXcd* result_0th; //= new MatrixXcd[nstep];
            MatrixXcd* result_nth; //= new MatrixXcd[nstep];


            int dim = 0;
            // Preserve a separate exponent through powers, local products and
            // MPI. DSJ/ITB variants retain their existing execution path.
            const bool stableProduct=!strcasecmp(quantity,"coherence") &&
                (isGCCE || !strcasecmp(method,"cce") || !strcasecmp(method,"pcce"));
            std::vector<ScaledComplex> stableWith(stableProduct ? nstep : 0);
            std::vector<ScaledComplex> stableWithout(stableProduct ? nstep : 0);
            // mode "except" (adaptive-regroup-stage3b.md 1절): the run's own product is untouched; beside it P_K, the
            // product of every term not in the except list, and the raw curves of the excepted terms as probes.
            const bool exceptMode=adaptiveWork(cnf).mode=="except";
            std::vector<ScaledComplex> keptWith(exceptMode ? nstep : 0),keptWithout(exceptMode ? nstep : 0);
            std::unique_ptr<std::ofstream> exceptCurves;
            if(exceptMode)exceptCurves=adaptiveExceptCurves(cnf,istate);
            long long exceptMatched=0;
            bool reportedInvalidProduct=false;
            double dmDropped=0.0;   // the largest cluster value in an element the free evolution keeps at zero (dm only)
            if (isGCCE && strcasecmp(quantity,"dm")==0){
                dim = QubitArray_dim(qa);
            }else{
                dim = 1;
            }          

            for (int istep=0; istep<nstep; istep++){
                result_wD[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                result_nD[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
            }
            ////////////////////////////////

            ////////////////////////////////
            // Zero-th order 
            double steptime_sta = MPI_Wtime();
            double steptime_end = 0.0;

            if (isGCCE){ // gcce
                // Zero-th order
                if (rank==0){
                    printf("\n");
                    printSubTitle("Calculate 0-th cluster...");
                }
                result_0th = calCoherenceGcce(qa,NULL,cnf,pls,op);
                if (rank==0){
                    steptime_end = MPI_Wtime();
                    // make message
                    sprintf(message,"# Wall time : %.2f[s]",steptime_end-steptime_sta);
                    printMessage(message);
                    printLine(); printf("\n");
                }
            }else{
                result_0th = new MatrixXcd[nstep];
                for (int istep=0; istep<nstep; istep++){
                    result_0th[istep] = MatrixXcd::Constant(dim, dim, doublec(1.0,0.0));
                }
            }

            // Product assembly (improvements.md A9): the savemode test and the inverse of the zeroth cluster are taken
            // once per state, not per cluster and time point.
            const bool saveAll=!strcasecmp(savemode,"all");
            std::vector<ScaledComplex> inverseZero;
            bool unitZero=false;
            if(stableProduct){
                inverseZero.resize(nstep);unitZero=true;
                for(int istep=0;istep<nstep;istep++){
                    const auto zero=result_0th[istep](0,0);
                    inverseZero[istep]=ScaledComplex(zero).inverse();
                    unitZero=unitZero && zero.real()==1.0 && zero.imag()==0.0 && !std::signbit(zero.imag());
                }
            }

            ladder.reset();
            lap(PH_PREP);
            MPI_Barrier(MPI_COMM_WORLD);
            lap(PH_WAIT);
            if(phases)++phaseStates;
            ////////////////////////////////

            ////////////////////////////////
            // Higher orders
            for (int n=1; n<=order; n++){

                // Message
                if (rank==0){
                    printf("\n");
                    sprintf(message,"Calculate %d-th cluster...",n);
                    printSubTitle(message);
                }

                // Number of clusters
                int ncluster = localclusters[n][0][0]-1;
                std::vector<BathArray*> batchBaths;
                std::vector<MatrixXcd*> batchResults;
                int batchFirst=0;
                // Create BathArray for the cluster (defect spins can be added)
                for (int ic=1; ic<=ncluster; ic++){

                    // Get only spin index part from cluster
                    int cluster[n];
                    copyInt1dPart(cluster,localclusters[n][ic],1,n);

                    // Get the power
                    int power = localclusters[n][ic][0];

                    // Create the bath array for the cluster
                    BathArray* ba_cluster = batchResults.empty() ?
                        createBathArray(cluster,n,ba,dfa,nqubit) : batchBaths[ic-batchFirst];
                    lap(PH_CREATE);
                    // BathArray_report(ba_cluster);

#ifdef CCEX_USE_CUDA
                    const int qubitFactor=isGCCE ? QubitArray_dim(qa) : 1;
                    if(batchResults.empty() && (!isGCCE || strcasecmp(savemode,"allfull")) &&
                       cnf->gpu_batch_clusters>1 && computeUseCuda(cnf,qubitFactor*BathArray_dim(ba_cluster))){
                        const int d=qubitFactor*BathArray_dim(ba_cluster);
                        const int timeBatch=std::min(cnf->gpu_batch_steps,std::max(1,nstep));
                        // A conservative buffer estimate limits bundleing. A single
                        // cluster still follows the existing allocation/error path.
                        const bool vectorState=isGCCE ? !strcasecmp(cnf->evolution,"vector") : cnf->nstate>0;
                        const long double bytes=gpuClusterMemoryEstimate(d,qubitFactor,timeBatch,pls->npulse,isGCCE,vectorState);
                        const long double budget=1048576.0L*cnf->gpu_cluster_batch_mb;
                        batchBaths.push_back(ba_cluster);
                        while(batchBaths.size()<(size_t)cnf->gpu_batch_clusters &&
                              (size_t)ic+batchBaths.size()<=(size_t)ncluster &&
                              bytes*(batchBaths.size()+1)<=budget){
                            const int next=ic+int(batchBaths.size());
                            BathArray* candidate=createBathArray(localclusters[n][next]+1,n,ba,dfa,nqubit);
                            if(qubitFactor*BathArray_dim(candidate)!=d){BathArray_freeAll(candidate);break;}
                            batchBaths.push_back(candidate);
                        }
                        if(batchBaths.size()>1){
                            batchFirst=ic;
                            batchResults=isGCCE ? calCoherenceGcceCudaBatch(qa,batchBaths,cnf,pls,op) :
                                                  calCoherenceCceCudaBatch(qa,batchBaths,cnf,pls);
                        }else batchBaths.clear();
                    }
#endif

                    // Calculate the coherence
                    bool ownedResult = true;                  // false: the reused cCCE buffer, not deleted here
                    if(!batchResults.empty()){
                        result_nth=batchResults[ic-batchFirst];
                    }else if (isGCCE){ // gcce: into the reused buffer (A14)
                        calCoherenceGcceInto(qa, ba_cluster, cnf, pls, op, cceResult.get(), &gcceAlphaRho);
                        result_nth = cceResult.get();
                        ownedResult = false;
                    }else{ // cce: into this calculation's buffer and workspace, reused from cluster to cluster (A14)
                        calCoherenceCceInto(qa, ba_cluster, cnf, pls, cceResult.get(), cceWork);
                        result_nth = cceResult.get();
                        ownedResult = false;
                    }
                    lap(PH_CLUSTER);

                    int exceptIndex=-1;
                    if(exceptMode && (exceptIndex=adaptiveExceptIndex(cluster,n))>=0) {
                        adaptiveExceptValues(*exceptCurves,exceptIndex,result_nth,nstep);
                        ++exceptMatched;
                    }

                    if(cnf->cluster_contributions)
                        contributions.capture(cluster,n,(isGCCE ? QubitArray_dim(qa) : 1)*BathArray_dim(ba_cluster),result_nth);

                    // Update the total result
                    // zeroth cluster : result_nth = L_C_(k=n) * ( 1 / L_C_(k=0) )
                    // power : result_nth = ( L_C_(k=n) / L_C_(k=0) ) ^ power
                    // total result : result = result * ( L_C_(k=n) / L_C_(k=0) ) ^ power
                    ladder.add(cluster,n,power,result_nth,inverseZero,unitZero);
                    for (int istep=0; istep<nstep; istep++){
                        if((cnf->cluster_contributions || ladder.on) && power==0)continue;
                        if(stableProduct){
                            const auto raw=result_nth[istep](0,0),zero=result_0th[istep](0,0);
                            const bool invalid=!std::isfinite(raw.real()) || !std::isfinite(raw.imag()) ||
                                !std::isfinite(zero.real()) || !std::isfinite(zero.imag()) || zero==doublec(0,0) ||
                                (power<0 && raw==doublec(0,0));
                            if(invalid && !reportedInvalidProduct) {
                                fprintf(stderr,"[numerical-invalid] rank=%d state=%d step=%d spin_order=%d local_cluster=%d power=%d raw=(%.17g,%.17g) L0=(%.17g,%.17g)\n",
                                    rank,istate,istep,n,ic,power,raw.real(),raw.imag(),zero.real(),zero.imag());
                                reportedInvalidProduct=true;
                            }
                            std::complex<double> saved;
                            assembleFactor(stableWith[istep],stableWithout[istep],raw,inverseZero[istep],power,unitZero,
                                           saveAll ? &saved : nullptr);
                            if(exceptMode && exceptIndex<0)
                                assembleFactor(keptWith[istep],keptWithout[istep],raw,inverseZero[istep],power,unitZero,nullptr);
                            if(saveAll)result_nth[istep](0,0)=saved;
                            continue;
                        }
                        // rho_C / rho_0, element by element; an element rho_0 keeps at zero gets 1 (divMatrixXcdElementWise)
                        result_nth[istep] = divMatrixXcdElementWise(result_nth[istep], result_0th[istep], &dmDropped);
                        result_nth[istep] = powMatrixXcdElementWise(result_nth[istep], power);

                        if (power>0){
                            result_nD[istep] = mulMatrixXcdElementWise(result_nD[istep], result_nth[istep]);
                        }

                        result_wD[istep] = mulMatrixXcdElementWise(result_wD[istep], result_nth[istep]);
                    }

                    lap(PH_ASSEMBLE);
                    if (strcasecmp(savemode,"all")==0){
                        Output_save_all(op,result_nth,nstep,deltat,cluster,n,istate);
                    }

                    if (rank==0){
                        steptime_end = MPI_Wtime();
                        if (ic % std::max(1, 5000/nprocess) == 0){
                            sprintf(message,"At rank 0, %d-th cluster (%d%%) is computed at %.2f s",ic, int(ic*100/ncluster), steptime_end-steptime_sta);
                            printMessage(message);

                        }else if (ic == ncluster){
                            sprintf(message,"At rank 0, %d-th cluster (%d%%) is computed at %.2f s",ic, int(ic*100/ncluster), steptime_end-steptime_sta);
                            printMessage(message);
                        }
                    }
                    if (ownedResult) delete[] result_nth;

                    // createBathArray builds ba_cluster fresh for every cluster and it was
                    // never released -- one BathArray, its BathSpin objects, their hypf
                    // arrays and every MatrixXcd inside, leaked per (state, config,
                    // cluster). Safe to free here: BathArray_setBath_i copies name,
                    // scalars, quad, hypf_sub and each hypf[j] BY VALUE (Eigen's
                    // operator= deep-copies), so ba_cluster aliases nothing in the parent
                    // bath, and result_nth is an independent new[] that has already been
                    // consumed above.
                    BathArray_freeAll(ba_cluster);
                    if(!batchResults.empty() && ic-batchFirst+1==(int)batchResults.size()){
                        batchResults.clear();batchBaths.clear();
                    }
                    lap(PH_FREE);
                }

                if (rank==0){
                    steptime_end = MPI_Wtime();
                    // make message
                    printf("\n");
                    sprintf(message,"# Wall time : %.2f[s]",steptime_end-steptime_sta);
                    printMessage(message);
                    printLine(); printf("\n");
                }
            }


            for(auto& x:stableWith)x.normalize();       // the lazy products are normalized before they are combined (A9)
            for(auto& x:stableWithout)x.normalize();
            for(auto& x:keptWith)x.normalize();
            if(exceptMode) {
                exceptCurves->flush();
                if(!*exceptCurves){fprintf(stderr,"Error: adaptive work (rank %d): writing except curves failed\n",rank);ccexFatal();}
                exceptCurves.reset();
                // P_K as defined needs L0 = 1 (cCCE); the except list is then the whole difference (stage3b 1절)
                if(!unitZero){fprintf(stderr,"Error: adaptive work (rank %d): except mode needs L0 = 1\n",rank);ccexFatal();}
                adaptiveExceptCheck(exceptMatched);
            }
            ladder.normalize();
            lap(PH_ASSEMBLE);
            MPI_Barrier(MPI_COMM_WORLD);
            lap(PH_WAIT);
            if (rank==0){
                printSubTitle("Product local results...");
                // zeroth cluster  L_C_(k=0) * PI_C_(k=1){ L_C_(k=1) } * PI_C_(k=2){ L_C_(k=2) }  .. 
                for (int istep=0; istep<nstep; istep++){
                    if(stableProduct){
                        const ScaledComplex zero(result_0th[istep](0,0));
                        stableWith[istep]=stableWith[istep]*zero;
                        stableWithout[istep]=stableWithout[istep]*zero;
                        if(exceptMode)keptWith[istep]=keptWith[istep]*zero;
                        for(int k=0;k<ladder.K && ladder.on;++k){   // the rungs take L0 at the same place as the run
                            ladder.with[k][istep]=ladder.with[k][istep]*zero;
                            ladder.without[k][istep]=ladder.without[k][istep]*zero;
                        }
                        continue;
                    }
                    result_nD[istep] = mulMatrixXcdElementWise(result_nD[istep], result_0th[istep]);
                    result_wD[istep] = mulMatrixXcdElementWise(result_wD[istep], result_0th[istep]);
                }           
            }

            if (dim > 1){
                double dmDroppedAll=0.0;
                MPI_Reduce(&dmDropped,&dmDroppedAll,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
                if (rank==0 && dmDroppedAll > 0.0){
                    printf("[dm] elements the free evolution keeps at zero are left at zero; clusters put values up to %.3e into them, which the elementwise product cannot carry (state %d)\n",dmDroppedAll,istate);
                }
            }
            contributions.finish(result_0th,istate,op);
            delete[] result_0th;
            

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
            if(stableProduct){
                int lengths[]={2,1};
                MPI_Aint offsets[]={offsetof(ScaledComplex,real),offsetof(ScaledComplex,exponent)};
                MPI_Datatype fields[]={MPI_LONG_DOUBLE,MPI_LONG_LONG},packed,type;
                MPI_Type_create_struct(2,lengths,offsets,fields,&packed);
                MPI_Type_create_resized(packed,0,sizeof(ScaledComplex),&type);
                MPI_Type_commit(&type);MPI_Type_free(&packed);
                MPI_Op multiply;MPI_Op_create(reduceScaledProduct,1,&multiply);
                std::vector<ScaledComplex> with(nstep),without(nstep);
                MPI_Reduce(stableWith.data(),with.data(),nstep,type,multiply,0,MPI_COMM_WORLD);
                MPI_Reduce(stableWithout.data(),without.data(),nstep,type,multiply,0,MPI_COMM_WORLD);
                if(exceptMode) {
                    std::vector<ScaledComplex> kept(nstep);
                    MPI_Reduce(keptWith.data(),kept.data(),nstep,type,multiply,0,MPI_COMM_WORLD);
                    adaptiveWriteExcept(cnf,kept,istate);
                }
                for(int k=0;k<ladder.K && ladder.on;++k){   // rdip ladder: the same reduction and normalization per rung
                    std::vector<ScaledComplex> lw(nstep),ln(nstep);
                    MPI_Reduce(ladder.with[k].data(),lw.data(),nstep,type,multiply,0,MPI_COMM_WORLD);
                    MPI_Reduce(ladder.without[k].data(),ln.data(),nstep,type,multiply,0,MPI_COMM_WORLD);
                    if(rank==0){
                        MatrixXcd* wD=new MatrixXcd[nstep];MatrixXcd* nD=new MatrixXcd[nstep];
                        // as the run's own product below: gCCE divides by its L(t=0)
                        const auto normalization=isGCCE ? lw[0].inverse() : ScaledComplex();
                        for(int i=0;i<nstep;++i){
                            ScaledComplex w=lw[i],n=ln[i];
                            if(isGCCE){w=w*normalization;n=n*normalization;}
                            wD[i]=MatrixXcd::Constant(1,1,w.value());nD[i]=MatrixXcd::Constant(1,1,n.value());
                        }
                        rdipladder::save(op,ladder.r[k],wD,nD,nstep,deltat,istate);
                        delete[] wD;delete[] nD;
                    }
                }
                MPI_Op_free(&multiply);MPI_Type_free(&type);
                if(rank==0){
                    // Only gCCE has the historical L(t=0) normalization.
                    const auto normalization=isGCCE ? with[0].inverse() : ScaledComplex();
                    for(int i=0;i<nstep;++i){
                        result_wD_all[i](0,0)=(with[i]*normalization).value();
                        result_nD_all[i](0,0)=(without[i]*normalization).value();
                    }
                }
            }else for (int istep=0; istep<nstep; istep++){
                err = MPI_Reduce(result_wD[istep].data(),result_wD_all[istep].data(),dim*dim,MPI_DOUBLE_COMPLEX,MPI_PROD,0,MPI_COMM_WORLD);
                err = MPI_Reduce(result_nD[istep].data(),result_nD_all[istep].data(),dim*dim,MPI_DOUBLE_COMPLEX,MPI_PROD,0,MPI_COMM_WORLD);
            }

            MPI_Barrier(MPI_COMM_WORLD);
            lap(PH_REDUCE);

            if (rank !=0){
                delete[] result_wD_all;
                delete[] result_nD_all;
            }

            delete[] result_wD;
            delete[] result_nD;

            ////////////////////////////////
            // Normalize the results
            ////////////////////////////////
            if (rank == 0 && !stableProduct){
                std::complex<double> normalization;
                if ( strcasecmp(quantity,"coherence")==0 && isGCCE){
                    normalization = result_wD_all[0].coeff(0);
                }else{
                    normalization = std::complex<double>(1, 0);
                }
                // do normalization : Ltot / normalization
                for (int i=0; i<nstep; i++){
                    result_wD_all[i] /= normalization;
                    result_nD_all[i] /= normalization;
                }
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            ////////////////////////////////
            // Save the results
            ////////////////////////////////
            if (rank==0){

                printSubTitle("Save the results...");
                adaptiveWriteCurve(cnf,result_wD_all,istate);
                
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
                // print wall clock time
                steptime_end = MPI_Wtime();
                printLine();
                printf("          Wall time = %.5f s\n", steptime_end - time_start_step);
                printLine();
                printf("\n\n");
            }

            ////////////////////////////////
            // Free the results
            ////////////////////////////////
            if (rank==0){
                delete[] result_wD_all;
                delete[] result_nD_all;
            }
            MPI_Barrier(MPI_COMM_WORLD);
            lap(PH_SAVE);
        }
    }
    // every rank, after its last probe state: the node budget is freed collectively (A1)
    if(adaptiveWork(cnf).mode=="probe" || adaptiveWork(cnf).mode=="factors") adaptiveProbeFinish();
    if(phases && phaseStates>0){
        // Gathered and printed by rank 0 in rank order, like computeProfileFlush, so records of ranks do not interleave.
        std::vector<double> all(rank==0 ? PH_N*nprocess : 0);
        MPI_Gather(phase,PH_N,MPI_DOUBLE,all.data(),PH_N,MPI_DOUBLE,0,MPI_COMM_WORLD);
        if(rank==0)for(int r=0;r<nprocess;++r){
            const double* p=all.data()+PH_N*r;
            printf("[phase] rank=%d states=%d prep_s=%.6f create_s=%.6f cluster_s=%.6f assemble_s=%.6f free_s=%.6f wait_s=%.6f "
                   "reduce_s=%.6f save_s=%.6f\n",r,phaseStates,p[PH_PREP],p[PH_CREATE],p[PH_CLUSTER],p[PH_ASSEMBLE],p[PH_FREE],
                   p[PH_WAIT],p[PH_REDUCE],p[PH_SAVE]);
        }
        fflush(stdout);
    }
}

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, DefectArray* dfa, int nqubit){

    // Initialize  ba_cluster
    BathArray* ba_cluster = BathArray_init();
    int nspin_cluster = 0;

    // Create the bath array for the cluster
    int addedspin = 0;
    for (int ic=0; ic<nspin; ic++){

        // Get BathSpin within the clusters
        int ibs = cluster[ic];
        BathSpin* bs  = BathArray_getBath_i(ba,ibs);

        // Allocate the bathArray
        nspin_cluster+=1;
        BathArray_setNspin(ba_cluster, nspin_cluster);

        if (nspin_cluster==1){ 
            BathArray_allocBath(ba_cluster,nqubit); 
        }
        else{                  
            BathArray_reallocBath(ba_cluster,nspin_cluster-1,nspin_cluster,nqubit); 
        }
        
        // Set BathSpin in cluster batharray
        BathArray_setBath_i(ba_cluster, bs, nspin_cluster-1, nqubit);

        //if (rank==0 && verbosity){
        //    BathArray_reportBath(ba_cluster);
        //}

        // Add sub-spins in cluster
        char* name = BathArray_getBath_i_name(ba,ibs);
        int idf = DefectArray_findDefectIndex(dfa,name);
        // If the bath is paramagnetic defect bath,
        //  additional spin might need to be added
        if (idf>=0){
            // has additional spin information
            int naddspin = DefectArray_getDefect_idf_naddspin(dfa,idf);
            for (int isp=0; isp<naddspin; isp++){
                
                // Reallocate the bath array
                nspin_cluster +=1;
                BathArray_setNspin(ba_cluster,nspin_cluster);
                BathArray_reallocBath(ba_cluster,nspin_cluster-1,nspin_cluster,nqubit);

                // Add the additional spin to BathArray_for_cluster (ba_cluster)
                BathSpin* bs_sub = DefectArray_getSubbath_i_isp(dfa,ibs,isp);

                // Set Addtional bath spin in cluster batharray
                BathArray_setBath_i(ba_cluster,bs_sub,nspin_cluster-1,nqubit);

                // Add the additional spin to the cluster
                addedspin++;
            }
        }
    }

    if (rank==0 && addedspin>0 && verbosity){
        printLineSection();
        BathArray_reportBath(ba_cluster);
        printStructElementInt(" * Added spin",addedspin);
        printf("\n");

        BathArray_reportBath_states(ba_cluster);
        BathArray_reportBath_detunings(ba_cluster);
        BathArray_reportBath_disorders(ba_cluster);
        BathArray_reportBath_hypf(ba_cluster,nqubit);
        BathArray_reportBath_quad(ba_cluster);
        BathArray_reportBath_hypf_sub(ba_cluster);
    }

    return ba_cluster;
}
