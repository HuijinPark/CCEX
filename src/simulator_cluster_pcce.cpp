#include "../include/cluster.h"
#include "../include/bath.h"
#include "../include/cluster_pcce.h"
#include "../include/memory.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "mpi.h"

void simulator_cluster_partition(BathSpin** bath        , Partition_info* pinfo             , Cluster* cls,
                                 Point*     best_centers, int*            best_assigned_idx               ){
    
    clock_t start     = clock();
    char* inputfile   = NULL;
    //////////////////////////////////////////////////////////
    ////            !!  Setting MPI process  !!             //
    //////////////////////////////////////////////////////////
    int max_trial       = cls      -> max_trial;
    int trials_per_proc = max_trial / nprocess;
    int start_trial     = rank      * trials_per_proc;
    int end_trial       = (rank == nprocess - 1) ? max_trial : start_trial + trials_per_proc;

    //for (int i=0; i<nprocess; i++){
    //    if (rank == i){
    //        printf("Rank: %d\n", rank);
    //        printf("trials_per_proc: %d\n", trials_per_proc);
    //        printf("start_trial: %d\n", start_trial);
    //        printf("end_trial: %d\n", end_trial);
    //        print_BD("=", 55);printf("\n");
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    //}
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    ////        !!  Setting the information array !!        //
    //////////////////////////////////////////////////////////
    int pcce_nspin = pinfo -> pcce_nspin;
    int ncenter    = pinfo -> ncenter;

    Point* local_best_spins      = (Point*) malloc(sizeof(Point) * pcce_nspin);
    Point* local_best_centers    = (Point*) malloc(sizeof(Point) * ncenter  );
    int* local_best_assigned_idx = (int*)   calloc(pcce_nspin, sizeof(int));
    int* local_best_shuffle_idx  = (int*)   calloc(pcce_nspin, sizeof(int));

    double **local_distances = allocDouble2d(pcce_nspin, ncenter);
    
    //////////////////////////////////////////////////////////
    ////                     !!  Run !!                     //
    //////////////////////////////////////////////////////////
    double local_best_inertia =  INFINITY;
    double local_best_sil     = -INFINITY;
    int    local_best_trial   =  0;
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    for (int trial = start_trial; trial < end_trial; trial++){

        Point* local_spins          = (Point*) malloc(sizeof(Point) * pcce_nspin);
        Point* local_centers        = (Point*) malloc(sizeof(Point) * ncenter  );
        int* local_assigned_idx     = (int*)   calloc(pcce_nspin, sizeof(int));
        int* local_shuffle_idx      = (int*)   calloc(pcce_nspin, sizeof(int));
        copy_bath_to_point_type(bath, local_spins, pcce_nspin);

        for (int i=0; i<pcce_nspin; i++){ 
            local_shuffle_idx[i] = i;
        }

        // Shuffle "spins" & initialize the "local_centers". // 
        shuffle_spins(local_spins, local_shuffle_idx, pcce_nspin);
        choose_initial_method(cls->kmeans_pp, local_spins, local_centers, pcce_nspin, ncenter);

        // Iteration info. //
        int iteration       = 0;
        bool converged      = false;
        clock_t trial_start = clock();

        // Start the k-means optimization in local process. //
        while (!converged && iteration < cls->max_iter) {
        
            calc_dist_to_centers(local_spins, local_centers, local_distances, pcce_nspin, ncenter);
            assignCentroids(local_assigned_idx, local_distances, pcce_nspin, ncenter);
            converged = updateCentroids(local_spins, local_assigned_idx, local_centers, pcce_nspin, ncenter);
            iteration++;

            if (cls->iter_detail == true){
                print_kmeans_converge_info(trial, iteration, ncenter, &local_centers);
            }
        }

        // Find the local optimization centers configuration. // 
        double local_inertia = calculateInertia(local_spins, local_centers, local_assigned_idx, pcce_nspin, ncenter);
        if (local_inertia < local_best_inertia){

            double local_sil = silhouetteCoefficient(local_spins, local_assigned_idx, local_centers, pcce_nspin, ncenter);
            //printf("\n\n\nLOCAL & BEST: %.2lf & %.2lf\n\n\n", local_sil, local_best_sil);
            local_best_trial   = trial;
            local_best_inertia = local_inertia;
            local_best_sil     = local_sil;

            copy_best_centers(local_centers, local_best_centers, ncenter);
            copy_best_assigned(local_assigned_idx, local_best_assigned_idx, pcce_nspin);
            copy_best_idxList(local_shuffle_idx, local_best_shuffle_idx, pcce_nspin);
            copy_best_spins(local_spins, local_best_spins, pcce_nspin);
        }
        //else if (local_inertia < local_best_inertia) {
        else if (local_inertia == local_best_inertia) {
            double local_sil = silhouetteCoefficient(local_spins, local_assigned_idx, local_centers, pcce_nspin, ncenter);
            if (local_sil > local_best_sil){
                local_best_trial   = trial;
                local_best_inertia = local_inertia;
                local_best_sil     = local_sil;

                copy_best_centers(local_centers, local_best_centers, ncenter);
                copy_best_assigned(local_assigned_idx, local_best_assigned_idx, pcce_nspin);
                copy_best_spins(local_spins, local_best_spins, pcce_nspin);
                copy_best_idxList(local_shuffle_idx, local_best_shuffle_idx, pcce_nspin);
                //printf("\n\n\nLOCAL & BEST: %.2lf & %.2lf\n\n\n", local_sil, local_best_sil);
            }
        }

        ////////////////////////////////////////////////////////
        clock_t trial_end = clock();
        for (int i=0; i < nprocess; i++){
            if (rank==i){
                print_optimize_kmeans_info(trial, iteration, local_inertia, local_best_inertia, local_best_sil, rank);
                printf("\n||---------------- Optimization Time ----------------||\n\n");
                printf("                 "); print_time(trial_start, trial_end);
                print_BD("=", 55); printf("\n");
                fflush(stdout);
            }
        MPI_Barrier(MPI_COMM_WORLD);
        }
        ////////////////////////////////////////////////////////
        free(local_spins);
        free(local_centers);
        free(local_assigned_idx);
        free(local_shuffle_idx);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //for (int i=0; i<nprocess; i++){
    //    if (rank == i){
    //        printf("\n||-------------------  Rank : %d | Best Local -------------------||\n\n", rank);
    //        print_BD("=", 55);
    //        print_centers(ncenter, local_best_centers, rank);
    //        print_assign_idx(local_best_assigned_idx, pcce_nspin, rank);
    //        print_spins(local_best_spins, pcce_nspin, rank);
    //        print_BD("=", 55);
    //        fflush(stdout);
    //    }
    //    MPI_Barrier(MPI_COMM_WORLD);
    //}
    //free_2Darr_double(&local_distances, pcce_nspin);
    freeDouble2d(&local_distances, pcce_nspin);

    ////////////////////////////////////////////////////////
    //  !!  Setting the "All & Global best" varibles  !!  //
    ////////////////////////////////////////////////////////
    double*  all_best_inertia = (double*)calloc(nprocess, sizeof(double));
    double*  all_best_sil     = (double*)calloc(nprocess, sizeof(double));
    int*     all_best_trial   = (int*)   calloc(nprocess, sizeof(int));

    int    best_rank           =  0;
    int    global_best_trial   =  0;
    double global_best_inertia =  INFINITY;
    double global_best_sil     = -INFINITY;

    ////////////////////////////////////////////////////////
    //   Gathering local best parameters in each procs.   //
    ////////////////////////////////////////////////////////
    MPI_Gather(&local_best_inertia, 1, MPI_DOUBLE, all_best_inertia, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_best_sil    , 1, MPI_DOUBLE, all_best_sil    , 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_best_trial  , 1, MPI_INT   , all_best_trial  , 1, MPI_INT   , 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ////////////////////////////////////////////////////////
    //           !!  Find the best rank index !!          //
    ////////////////////////////////////////////////////////
    if (rank == 0) {
        for (int i = 0; i < nprocess; i++){
            if (all_best_inertia[i] < global_best_inertia || 
                all_best_inertia[i] == global_best_inertia && all_best_sil[i] < global_best_sil) {

                best_rank = i;
                update_global_best(i, all_best_inertia, all_best_sil, all_best_trial, &global_best_inertia, &global_best_sil, &global_best_trial);
            }
        
        }
        free(all_best_inertia);
        free(all_best_sil);
        free(all_best_trial);

        printf("\n"); print_BD("=", 55); printf("\n");
        printf("Best rank: %d, Iner.: %.2lf, Sil.: %.2lf, Trial: %d\n", best_rank, global_best_inertia, global_best_sil, global_best_trial);
        printf("\n"); print_BD("=", 55); printf("\n");
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Bcast(&best_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    ////////////////////////////////////////////////////////
    //     !!  Send & Receive the best rank resuts !!     //
    ////////////////////////////////////////////////////////
    // Chnage the indecies for "local_spins => ba->bath"
    Point* global_best_centers              = local_best_centers;

    Point* global_best_inverse_spins        = (Point*)malloc(pcce_nspin*sizeof(Point));
    int*   global_best_inverse_assigned_idx = (int*)calloc(pcce_nspin, sizeof(int));
    int*   global_best_inverse_shuffle_idx  = (int*)calloc(pcce_nspin, sizeof(int)); 
    changeIdx_inverse_shuffle_idx(local_best_shuffle_idx, global_best_inverse_shuffle_idx, pcce_nspin);
    
    if (rank == best_rank){
        for (int i=0; i<pcce_nspin; i++){
            int inverse_idx                               = local_best_shuffle_idx[i];
            global_best_inverse_assigned_idx[inverse_idx] = local_best_assigned_idx[i];
            global_best_inverse_spins[inverse_idx]        = local_best_spins[i];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(global_best_inverse_spins       , pcce_nspin* 3, MPI_DOUBLE, best_rank,  MPI_COMM_WORLD);
    MPI_Bcast(global_best_centers             , ncenter  * 3 , MPI_DOUBLE, best_rank,  MPI_COMM_WORLD);
    MPI_Bcast(global_best_inverse_assigned_idx, pcce_nspin   , MPI_INT   , best_rank,  MPI_COMM_WORLD);

    for (int i=0; i<ncenter; i++){
        best_centers[i].x = global_best_centers[i].x; 
        best_centers[i].y = global_best_centers[i].y; 
        best_centers[i].z = global_best_centers[i].z; 
    }
    for (int i=0; i<pcce_nspin; i++){
        best_assigned_idx[i] = global_best_inverse_assigned_idx[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0){
        printf("\n||-------------------------------------------------------||\n");
        printf("\n||----------------! FINISH OPTIMIZATION !----------------||\n");
        printf("\n||-------------------------------------------------------||\n\n");
        printf("\n||-----------------! Center Information !----------------||\n\n");
        print_centers(ncenter, best_centers, rank);
        printf("\n||--------------! Assigned Index for Center !------------||\n");
        print_assign_idx(best_assigned_idx, pcce_nspin, rank);
        print_BD("=", 55);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    ////      !!  Close the dynamical allocation !!         //
    //////////////////////////////////////////////////////////
    FILE *c_savefile; 
    FILE *s_savefile; 
    
    if (rank == 0){
        //write_centers(&c_savefile, &(ba->centers), ncenter);
        //write_spins(&s_savefile, &global_best_inverse_spins, &(ba->assigned_idx), pcce_nspin);
        //write_centers(&c_savefile, &best_centers, ncenter);
        //write_spins(&s_savefile, &global_best_inverse_spins, &best_assigned_idx, pcce_nspin);
        print_best_centers(&best_centers, ncenter);
        print_best_spins(&global_best_inverse_spins, &best_assigned_idx, pcce_nspin);
    }

    free(global_best_centers);
    free(global_best_inverse_spins);
    free(global_best_inverse_assigned_idx);
    free(global_best_inverse_shuffle_idx);

    free(local_best_spins);
    free(local_best_assigned_idx);

    ////////////////////////////////////////////////////////
    //                   !!  Finish  !!                   //
    ////////////////////////////////////////////////////////
    clock_t end = clock();

    if (rank==0) {
        printf("\n||--------------- Total Execution Time ---------------||\n\n");
        printf("                 "); print_time(start, end);
        print_BD("=", 55); printf("\n");
    }   
}
