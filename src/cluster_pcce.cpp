#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "../include/cluster_hash.h"
#include "../include/cluster_pcce.h"
#include "../include/memory.h"
#include "../include/qubit.h"
#include "../include/cluster.h"
#include "../include/bath.h"
#include "../include/hamiltonian.h"


void clusterizePcce(Cluster* cls, BathArray* ba, QubitArray* qa, Config* config){

    ////////////////////////////////////////////////////
    // About pCCE method :
    //  pCCE(N,K) indicates a clustering method of 
    //  the N-clusters having K number of bath spins
    ////////////////////////////////////////////////////
    //
    ////////////////////////////////////////////////////
    // pCCE related information : 
    //   sK, ncenter, pcce_nspin, rest_nspin 
    Partition_info pinfo;

    ////////////////////////////////////////////////////
    // Cut the # of spins (pcce grouping)
    // Total number of spins : ba.nspin -> NK
    ////////////////////////////////////////////////////
    set_spinfinite(ba, qa, cls->sK, &pinfo, rank);
    
    ////////////////////////////////////////////////////
    // Find the center position & assigned index 
    // (Find a representative position of a cluster 
    //  to be included k spins)
    //  Arg return : best_centers, best_assigned_idx
    ////////////////////////////////////////////////////
    int ncenter              = pinfo.ncenter; 
    int pcce_nspin           = pinfo.pcce_nspin;
    Point* best_centers      = (Point*) calloc(ncenter, sizeof(Point));
    int*   best_assigned_idx = (int*)   calloc(pcce_nspin, sizeof(int));
    simulator_cluster_partition(ba->bath, &pinfo, cls, best_centers, best_assigned_idx);
    
    ////////////////////////////////////////////////////
    // Clusterize for the center
    ////////////////////////////////////////////////////
    int order = Cluster_getOrder(cls);

    Cluster_setClusinfo_0th(cls);

    if (order == 1){ // For pCCE, + k=1 (order = n)
        
        if (rank==0){printf("\n\t Clustering 1st order ... \n");}
    
        int cluster[1] = {0};
        int iter = 1; // clusinfo[N][i][0]
        // center >> cls->info for o1
        for (int centeridx = 0; centeridx < ncenter; centeridx++){
            cluster[0] = centeridx; 
            Cluster_setClusinfo_addcluster(cls, order, iter, cluster);
        }
    }
    else if (order > 1){ 
    
        int** cmap    = NULL;  // connectivity map : cmap[nspin][nspin] = 1 if connected else 0
        float** stmap = NULL;  // strength map : stmap[nspin][nspin] = strength
        int** spmap   = NULL;  // sparse map : spmap[nspin][ncol] = n connected spin + 1d array of connected spins index
    
        float rdip    = Config_getRdip(config);
        float rdipcut = Config_getRdipcut(config);
    
        // center >> cls->info for on
        if (strcasecmp(cls->method, "pcce")==0){
            BathArray_connectivity_pcce(&cmap, &stmap, best_centers, rdip, rdipcut, ncenter);
            makeSparsemap(&spmap, cmap, ncenter);
            // main
            clusterizeHash(cls, ncenter, spmap, stmap); // clusterize + solve tilde
        }
        else{
            printf("Error: clusterize: method is not defined\n");
            exit(1);
        }
        freeInt2d(   &cmap, ncenter);
        freeInt2d(  &spmap, ncenter);
        freeFloat2d(&stmap, ncenter); 
    
    }else{
        fprintf(stderr,"Error: Cluster_clusterize: order(%d) is not defined\n",order);
        exit(1);
    }
    
    // ========================== //
    // Setting cls_pcce !! <= cls //
    // ========================== //
    int** cs2dArr = Cluster_setCenterIdx_spinIdx_2dArr(ncenter, cls->sK, pcce_nspin, best_assigned_idx);
    if (rank == 0){
        for (int tempi=0; tempi<ncenter; tempi++){
            printf("cs2dArr[%d] ", tempi);
            for (int tempj=0; tempj<cls->sK; tempj++){
                printf("%d ", cs2dArr[tempi][tempj]);
            }
            printf("\n");
        }
    }
    printf("\n");
    
    int*** pcceClusInfo = convert_centerIdx_to_spinIdx(cls, cls->order, cls->sK, cs2dArr);
    freeInt2d(&cs2dArr, ncenter);
    
    int spinOrder = (cls->order * cls->sK);
    Cluster_freeClusinfo(cls);
    cls->clusinfo = pcceClusInfo;
    Cluster_setOrder(cls, spinOrder);
    //printf("! ================= !\n");
    //printf("! -- cls print --!!\n\n");
    //reportClusinfo(cls->clusinfo, 4);
    //printf("\n! ================= !\n");
}

// Find the Center position & assigned index 
// arg return : center position(best_centers) and related info (pinfo, best_assigned_idx).
void simulator_cluster_partition(BathSpin** bath        , Partition_info* pinfo             , Cluster* cls,
                                 Point*     best_centers, int*            best_assigned_idx               ){
    
    clock_t start     = clock();
    char* inputfile   = NULL;
    //////////////////////////////////////////////////////////
    ////            !!  Setting MPI process  !!             //
    //////////////////////////////////////////////////////////
    int max_trial       = cls      -> max_trial;
    int trials_per_proc = max_trial / nprocess;
    if (trials_per_proc == 0){trials_per_proc = 1;}

    //int start_trial     = rank      * trials_per_proc;
    //int end_trial       = (rank == nprocess - 1) ? max_trial : start_trial + trials_per_proc;

    for (int i=0; i<nprocess; i++){
        if (rank == i){
            printf(" || Rank: %d\n", rank);
            printf(" || trials_per_proc: %d\n", trials_per_proc);
            //printf("start_trial: %d\n", start_trial);
            //printf("end_trial: %d\n", end_trial);
            print_BD("=", 55);printf("\n");
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //////////////////////////////////////////////////////////
    ////        !!  Setting the information array !!        //
    //////////////////////////////////////////////////////////
    int pcce_nspin = pinfo -> pcce_nspin;
    int ncenter    = pinfo -> ncenter;

    Point* local_best_spins      = (Point*) calloc(pcce_nspin, sizeof(Point));
    Point* local_best_centers    = (Point*) calloc(ncenter,    sizeof(Point));
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
    for (int trial = 0; trial < trials_per_proc; trial++){

        Point* local_spins          = (Point*) calloc(pcce_nspin, sizeof(Point));
        Point* local_centers        = (Point*) calloc(ncenter,    sizeof(Point));
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
        //else if (local_inertia < local_best_inertia) 
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

        //printf("\n"); print_BD("=", 55); printf("\n");
        //printf("Best rank: %d, Iner.: %.2lf, Sil.: %.2lf, Trial: %d\n", best_rank, global_best_inertia, global_best_sil, global_best_trial);
        //printf("\n"); print_BD("=", 55); printf("\n");
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

    //if (rank == 0){
    //    printf("\n||-------------------------------------------------------||\n");
    //    printf("\n||----------------! FINISH OPTIMIZATION !----------------||\n");
    //    printf("\n||-------------------------------------------------------||\n\n");
    //    printf("\n||-----------------! Center Information !----------------||\n\n");
    //    print_centers(ncenter, best_centers, rank);
    //    printf("\n||--------------! Assigned Index for Center !------------||\n");
    //    print_assign_idx(best_assigned_idx, pcce_nspin, rank);
    //    print_BD("=", 55);
    //}
    //MPI_Barrier(MPI_COMM_WORLD);

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

///////////////////////////////////////////////////
void update_global_best(int index, double *inertia, double *sil, 
                        int *trial, double *global_best_inertia, double *global_best_sil, int *global_best_trial) {

    *global_best_inertia = inertia[index];
    *global_best_sil     = sil[index];
    *global_best_trial   = trial[index];
}

///////////////////////////////////////////////
//       Util. & Setting the spin number     // 
///////////////////////////////////////////////
double calc_dist(double* arr1, double* arr2){

    double dx = arr1[0] - arr2[0];
    double dy = arr1[1] - arr2[1];
    double dz = arr1[2] - arr2[2];

    double dist = sqrt(dx*dx + dy*dy + dz*dz);
    return dist;
}

int compare_dist(const void *a, const void *b){
    BathSpin* s1 = *(BathSpin**)a;
    BathSpin* s2 = *(BathSpin**)b;
    if (s1->min_dist < s2->min_dist) return -1;
    if (s1->min_dist > s2->min_dist) return  1;
    return 0;
}

void print_pcce_info(int bath_nspin, int nqubit, int ncenter, int pcce_nspin, int rest_nspin){
    print_BD("-", 30);
    printf("Total Bath Spin: %d\n", bath_nspin);
    printf("Qubit #: %d\n", nqubit);
    printf("ncenter: %d\n", ncenter);
    printf("pcce_nspin: %d\n", pcce_nspin);
    printf("rest_nspin: %d\n", rest_nspin);
    print_BD("-", 30);
}

void shrink_restspins(BathArray* ba, int rest_nspin){

    if (rest_nspin > ba->nspin){
        printf("Rest_nspin exceeds the size of the BathArray ba.\n");
    }
    
    ba->nspin -= rest_nspin;
    BathSpin** temp = (BathSpin**)realloc(ba->bath, ba->nspin * sizeof(BathSpin*));
    if (ba->bath == NULL && ba->nspin>0){
        perror("Failed to realloc ba->bath\n");
        exit(EXIT_FAILURE);
    }
    ba->bath = temp;
}

void set_spinfinite(BathArray* ba, QubitArray* qa, int sK, Partition_info* pinfo, int rank){
    
    int ncenter    = (ba->nspin / sK);
    int pcce_nspin = (  ncenter * sK);
    int rest_nspin = (ba->nspin % sK);

    pinfo->sK         = sK;
    pinfo->ncenter    = ncenter;
    pinfo->pcce_nspin = pcce_nspin;
    pinfo->rest_nspin = rest_nspin;
    
    if (rank == 0){
        print_pcce_info(ba->nspin, qa->nqubit, ncenter, pcce_nspin, rest_nspin);
    }

    for (int spinIdx=0; spinIdx<(ba->nspin); spinIdx++){

        double min_dist = -1;
        for (int qubitIdx=0; qubitIdx<(qa->nqubit); qubitIdx++){

            //double dist = dist((qa->qubit[qubitIdx]->xyz), ba->bath[spinIdx]->xyz);
            double dist = calc_dist((qa->qubit[qubitIdx]->xyz), ba->bath[spinIdx]->xyz);
            if (min_dist == -1 || dist < min_dist){
                min_dist = dist;
            }
        }
        ba->bath[spinIdx]->min_dist = min_dist;
    }

    qsort(ba->bath, ba->nspin, sizeof(BathSpin*), compare_dist);
    shrink_restspins(ba, rest_nspin);
}

///////////////////////////////////////////////
//       Constrained k-means function.       // 
///////////////////////////////////////////////
int if_contain(int* arr, int arr_size, int value) {
    for (int k=0; k < arr_size;  k++){
        if (arr[k] == value) {
            return 1;
        }
    }
    return 0;           
}

// Inititialize the centers
void initializeCentroids(Point* spins, Point* centers, int ncenter, int pcce_nspin) {
    
    // Make and initialize the Random number index list // 
    int count = 0;
    int* RandomList = (int*) calloc(ncenter, sizeof(int));
    for (int i=0; i<ncenter; i++){
        RandomList[i] = (pcce_nspin+100);
    }
    
    // Random vales => srand //
    // srand(time(NULL)+i*i);
    for (int count=0; count <ncenter;) {
        int random_index = rand() % pcce_nspin;
        //printf("%d", random_index);
    
        if (!if_contain(RandomList, count, random_index)) {
            RandomList[count] = random_index;
            centers[count].x = spins[random_index].x;
            centers[count].y = spins[random_index].y;
            centers[count].z = spins[random_index].z;
            count++;
        }
    }
    //printf("\n\n Finish making the random number index list!!\n\n");
    free(RandomList);
}

void choose_initial_method(bool kmeans_pp, Point* spins, Point* centers, int pcce_nspin, int ncenter){
    //printf("||------------------- Start Iter. -------------------||\n\n");
    if (kmeans_pp == true){
        //printf(" - Centroids Initialization: K-Means++ \n\n");
        //printf(" || Centroids Initialization: K-Means++ ||\n\n");
        //print_BD("=", 55);
        kMeansPlusPlus(spins, centers, pcce_nspin, ncenter);
    }
    else {
        //printf(" - Centroids Initialization: Random \n");
        //printf(" || Centroids Initialization: Random ||\n\n");
        //print_BD("=", 55);
        initializeCentroids(spins, centers, ncenter, pcce_nspin);
    }
    printf("\n");
}

void kMeansPlusPlus(Point *spins, Point *centers, int pcce_nspin, int ncenter) {

    int index = rand() % pcce_nspin;  // 무작위 인덱스
    
    centers[0].x = spins[index].x;
    centers[0].y = spins[index].y;
    centers[0].z = spins[index].z;
    //centers[0] = spins[index];  // 첫 번째 중심점 선택

    double *dists = (double*) malloc(pcce_nspin * sizeof(double));

    for (int i = 1; i < ncenter; i++) {
        double sum = 0.0;
        // 각 포인트에 대한 최소 거리를 제곱하여 계산
        for (int j = 0; j < pcce_nspin; j++) {
            double d = distance(spins[j], centers[0]);  // 첫 번째 중심점과의 거리
            //double d = distance_SC(spins[j], centers[0]);  // 첫 번째 중심점과의 거리
            for (int c = 1; c < i; c++) {
                double tmp = distance(spins[j], centers[c]);
                //double tmp = distance_SC(spins[j], centers[c]);
                if (tmp < d) d = tmp;
            }
            dists[j] = d * d;
            sum += dists[j];
        }

        // 확률적으로 다음 중심점 선택
        double r = ((double)rand() / (RAND_MAX)) * sum;
        sum = 0;
        for (int j = 0; j < pcce_nspin; j++) {
            sum += dists[j];
            if (sum >= r) {
                //centers[i] = spins[j];
                centers[i].x = spins[j].x;
                centers[i].y = spins[j].y;
                centers[i].z = spins[j].z;
                break;
            }
        }
    }
    free(dists);
}
void shuffle_spins(Point *array, int* shuffle_idxs, int pcce_nspin) {
    // 배열의 마지막 요소부터 시작하여 첫 번째 요소 이전까지 반복
    for (int i = pcce_nspin - 1; i > 0; i--) {
        // 0부터 i까지의 무작위 인덱스 선택
        int j = rand() % (i + 1);

        Point temp = array[i];
        array[i]   = array[j];
        array[j]   = temp;

        int temp_idx    = shuffle_idxs[i];
        shuffle_idxs[i] = shuffle_idxs[j];
        shuffle_idxs[j] = temp_idx;

    }
}

void changeIdx_inverse_shuffle_idx(int* shuffle_idx, int* inverse_shuffle_idx, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        int j = shuffle_idx[i];
        inverse_shuffle_idx[j] = i;
    }
}

void assignCentroids(int* assigned_cluster_idx, double** distances, int pcce_nspin, int ncenter) {
    int spinsPerCenter = pcce_nspin / ncenter;
    int *centroidCounts = (int *)malloc(ncenter * sizeof(int));
    for (int i = 0; i < ncenter; i++) {
        centroidCounts[i] = 0;
    }

    for (int i = 0; i < pcce_nspin; i++) {
        int minIndex = -1;
        double minDistance = HUGE_VAL;
        for (int j = 0; j < ncenter; j++) {
            if (distances[i][j] < minDistance && centroidCounts[j] < spinsPerCenter) {
                minDistance = distances[i][j];
                minIndex = j;
            }
        }
        if (minIndex == -1) {
            // It may happen that no centroid can take more points, distribute remaining points
            for (int j = 0; j < ncenter; j++) {
                if (centroidCounts[j] < spinsPerCenter) {
                    minIndex = j;
                    break;
                }
            }
        }
        assigned_cluster_idx[i] = minIndex;
        centroidCounts[minIndex]++;
    }
    free(centroidCounts);
}

void calc_dist_to_centers(Point* spins, Point* centers, double** distances, int pcce_nspin, int ncenter) {
    for (int i = 0; i < pcce_nspin; i++) {
        for (int j = 0; j < ncenter; j++) {
            double dx = spins[i].x - centers[j].x;
            double dy = spins[i].y - centers[j].y;
            double dz = spins[i].z - centers[j].z;
            //dists[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            distances[i][j] = (dx * dx + dy * dy + dz * dz);
            //printf("|| DISTANCE[%d][%d] =%lf || \n", i, j, (dx * dx + dy * dy + dz * dz));
        }
    }
}

bool updateCentroids(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter) {

    bool converged = true;
    int *counts = (int *)malloc(ncenter * sizeof(int));
    Point *sums = (Point *)malloc(ncenter * sizeof(Point));
    
    // Initialize //
    for (int i = 0; i < ncenter; i++) {
        sums[i].x = sums[i].y = sums[i].z = 0;
        counts[i] = 0;
    }

    for (int i = 0; i < pcce_nspin; i++) {
        int cid = assignments[i];
        sums[cid].x += spins[i].x;
        sums[cid].y += spins[i].y;
        sums[cid].z += spins[i].z;
        counts[cid]++;
    }

    double TOLERANCE = 0.0001;
    for (int i = 0; i < ncenter; i++) {
        if (counts[i] > 0) {
            double newX = sums[i].x / counts[i];
            double newY = sums[i].y / counts[i];
            double newZ = sums[i].z / counts[i];
            // Check if centroid moved significantly
            if (fabs(newX - centers[i].x) > TOLERANCE || fabs(newY - centers[i].y) > TOLERANCE || fabs(newZ - centers[i].z) > TOLERANCE){
                converged = false;
            }
            centers[i].x = newX;
            centers[i].y = newY;
            centers[i].z = newZ;
        }
    }
    free(counts);
    free(sums);
    return converged;
}

double calculateInertia(Point *spins, Point *centers, int *assignments, int pcce_nspin, int ncenter) {
    double inertia = 0.0;
    for (int i = 0; i < pcce_nspin; i++) {
        int cid = assignments[i];
        double dx = spins[i].x - centers[cid].x;
        double dy = spins[i].y - centers[cid].y;
        double dz = spins[i].z - centers[cid].z;
        inertia += dx * dx + dy * dy + dz * dz;
    }
    return inertia;
}

double silhouetteCoefficient(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter) {

    double *a = (double*)malloc(pcce_nspin * sizeof(double)); // 클러스터 내 거리
    double *b = (double*)malloc(pcce_nspin * sizeof(double)); // 클러스터 간 거리
    double totalSilhouette = 0.0; // 전체 실루엣 계수 합

    // 모든 포인트에 대한 클러스터 내 거리 a(i) 초기화
    for (int i = 0; i < pcce_nspin; i++) {
        double sumDist = 0.0;
        int clusterSize = 0;
        for (int j = 0; j < pcce_nspin; j++) {
            if (assignments[i] == assignments[j]) {
                sumDist += distance(spins[i], spins[j]);
                clusterSize++;
            }
        }
        a[i] = sumDist / (clusterSize - 1); // 자기 자신을 제외
    }

    // 모든 포인트에 대한 클러스터 간 최소 거리 b(i) 계산
    for (int i = 0; i < pcce_nspin; i++) {
        double minDist = FLT_MAX;
        for (int c = 0; c < ncenter; c++) {
            if (c != assignments[i]) { // 다른 클러스터
                double clusterDist = 0.0;
                int clusterSize = 0;
                for (int j = 0; j < pcce_nspin; j++) {
                    if (assignments[j] == c) {
                        clusterDist += distance(spins[i], spins[j]);
                        clusterSize++;
                    }
                }
                double avgDist = clusterDist / clusterSize;
                if (avgDist < minDist) {
                minDist = avgDist;
                }
            }
        }
        b[i] = minDist;
    }
    
    // 각 포인트에 대한 실루엣 계수 계산
    for (int i = 0; i < pcce_nspin; i++) {
        double si = (b[i] - a[i]) / (a[i] > b[i] ? a[i] : b[i]);
        totalSilhouette += si;
    }
    double averageSilhouette = totalSilhouette / pcce_nspin;
    free(a);
    free(b);

    return averageSilhouette;
}

double distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

/////////////////////////////////////////////////////////////////////////////
// Function: Copy now array --> best array //
void copy_best_centers(Point* centers, Point* best_centers, int ncenter){
    for (int i=0; i<ncenter; i++){
        best_centers[i] = centers[i];
    }
}

void copy_best_assigned(int* assigned, int* best_assigned, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        best_assigned[i] = assigned[i];
    }
}

void copy_best_idxList(int* idxList, int* best_idxList, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        best_idxList[i] = idxList[i];
    }
}

void copy_best_spins(Point* spins, Point* best_spins, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        best_spins[i] = spins[i];
    }
}

void copy_bath_to_point_type(BathSpin** bath, Point* spins, int pcce_nspin){
    for (int i=0; i<pcce_nspin; i++){
        spins[i].x = bath[i]->xyz[0];
        spins[i].y = bath[i]->xyz[1];
        spins[i].z = bath[i]->xyz[2];
    }
}


///////////////////////////////////////////////
//            Print & Write Function         // 
///////////////////////////////////////////////
void print_BD(char* str, int num){
    printf("\n");
    for (int b=0; b<num; b++){
        printf(str);
    }
    printf("\n");
}

///////////////////////////////////
void print_centers(int ncenter, Point* centers, int rank){
    print_BD("-", 55);
    printf("  - Rank: %d\n\n", rank);
    for (int i=0; i < ncenter; i++){
        printf(" Point %5d: %10.3lf %10.3lf %10.3lf\n", i, centers[i].x, centers[i].y, centers[i].z);
    }
}

///////////////////////////////////
void print_kmeans_converge_info(int trial, int iteration, int ncenter, Point** centers){
    print_BD("-", 55);
    printf(" || Iter. %d - (%d) ||\n", trial, iteration);
    for (int i=0; i < ncenter; i++){
        printf(" Centroid %d: %10.5lf %10.5lf %10.5lf\n", i, (*centers)[i].x, (*centers)[i].y, (*centers)[i].z);
    }
}

///////////////////////////////////
void print_optimize_kmeans_info(int trial, int iteration, double inertia, double best_inertia, double best_sil, int rank){

    //print_BD("=", 55); 
    printf("\n||------------- Optimization Infomation -------------||\n\n");
    //printf(" ||  Iteration  infomation  ||\n\n");
    printf("    - Rank: %d \n\n", rank);
    printf("    - %d trials & %d iterations.\n\n", trial, iteration);
    printf("    - Iner. & Best Iner. : %.2lf & %.2lf\n\n", inertia, best_inertia);
    printf("    - Best Sil.  : %.2lf\n\n", best_sil);
    //print_BD("=", 55); 
}

///////////////////////////////////
void print_assign_idx(int *assign_idx, int pcce_nspin, int rank) {
    printf("\n || Point assignments to centers ||\n\n");
    printf(" - Rank: %d\n\n", rank);
    for (int i = 0; i < pcce_nspin ; i++) {
        printf(" Point %5d -> Center %4d\n", i, assign_idx[i]);
    }
    printf("\n\n");
}

void print_spins(Point *spins, int pcce_nspin, int rank){
    //printf("shuffle!!\n");
    printf(" - Rank: %d\n\n", rank);
    for (int i=0; i<pcce_nspin; i++){
        printf("Point %5d - %10.5lf %10.5lf %10.5lf\n", i, spins[i].x, spins[i].y, spins[i].z);
    }
}

void print_time(clock_t start, clock_t end){
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    int hours      = (int) cpu_time_used / 3600;
    int minutes    = ((int) cpu_time_used % 3600) / 60;
    double seconds = cpu_time_used - (hours * 3600) - (minutes*60);
    printf("    %d : %d : %.3f \n\n", hours, minutes, seconds);
}

///////////////////////////////////
void write_centers(FILE **c_savefile, Point **centers, int ncenter){

    *c_savefile = fopen("./centers.txt", "w");
    //*c_savefile = fopen("./out/centers.txt", "w");
    if (*c_savefile == NULL){
        printf("File could not be opend\n");
    }
    fprintf(*c_savefile, " x   y   z \n");
    for (int p=0; p<ncenter; p++){
        fprintf(*c_savefile, " %10.5lf %10.5lf %10.5lf\n", (*centers)[p].x, (*centers)[p].y, (*centers)[p].z);
    }
    fclose(*c_savefile);
}

///////////////////////////////////
void write_spins(FILE **s_savefile, Point **spins, int** assigned_idx, int pcce_nspin){

    *s_savefile = fopen("./spins.txt", "w");
    if (*s_savefile == NULL){
        printf("File could not be opend\n");
    }
    fprintf(*s_savefile, " x   y   z   center Idx\n");
    for (int p=0; p<pcce_nspin; p++){
        fprintf(*s_savefile, " %10.5lf %10.5lf %10.5lf %10d\n", (*spins)[p].x, (*spins)[p].y, (*spins)[p].z, (*assigned_idx)[p]);
    }
    fclose(*s_savefile);
}

///////////////////////////////////
void print_best_spins(Point **spins, int** assigned_idx, int pcce_nspin){

    printf(" ===============\n");
    printf(" ==== Spins ====\n");
    printf(" ===============\n");
    printf(" x   y   z   center Idx\n");
    for (int p=0; p<pcce_nspin; p++){
        printf(" %10.5lf %10.5lf %10.5lf %10d\n", (*spins)[p].x, (*spins)[p].y, (*spins)[p].z, (*assigned_idx)[p]);
    }
    printf(" ===============\n");
} 


void print_best_centers(Point **centers, int ncenter){

    printf(" ===========\n");
    printf(" = Centers =\n");
    printf(" ===========\n");
    printf(" x   y   z \n");
    for (int p=0; p<ncenter; p++){
        printf(" %10.5lf %10.5lf %10.5lf\n", (*centers)[p].x, (*centers)[p].y, (*centers)[p].z);
    }
    printf(" ===========\n");
}
///////////////////////////////////
void BathArray_connectivity_pcce(int*** cmap, float*** stmap, Point* centers, float rdip, float rdipcut, int ncenter){
     
    //Point* centers = ba->centers;
    // cmap  [ncenter][ncenter] = 1 if connected else 0
    // stmap [ncenter][ncenter] = strength of the connection

    // Connectivity Map and Strength Map
    *cmap  =   allocInt2d(ncenter, ncenter);
    *stmap = allocFloat2d(ncenter, ncenter);

    // Find connectivity Map and strength Map
    double r = 0.0;
    float strength = 0.0;
    MatrixXcd tensor;

    for (int sp1=0; sp1<ncenter; sp1++){
        for (int sp2=sp1+1; sp2<ncenter; sp2++){
            
            Point c1 = centers[sp1];  double cxyz1[3] = {c1.x, c1.y, c1.z};
            Point c2 = centers[sp2];  double cxyz2[3] = {c2.x, c2.y, c2.z};
        
            r        = dist(cxyz1, cxyz2);
            tensor   = calPointDipoleTensor(cxyz1, cxyz2, 10000, 10000);
            strength = tensor(2,2).real();

            //printf("r = %10.5lf\n", r);
            if ( r < rdip && r > rdipcut){
                (*cmap)[sp1][sp2] = 1;
                (*cmap)[sp2][sp1] = 1;
            }

            (*stmap)[sp1][sp2] = fabs(strength);
            (*stmap)[sp2][sp1] = fabs(strength);
        }
    }
}

///////////////////////////////////
int** Cluster_setCenterIdx_spinIdx_2dArr(int ncenter, int sK, int pcce_nspin, int* best_assigned_idx){

    int** cs2dArr = allocInt2d(ncenter, sK);
    int*  assignNumArr = allocInt1d(ncenter);

    for (int i=0; i<pcce_nspin; i++){
        int centerIdx = best_assigned_idx[i];
        int assignNum = assignNumArr[centerIdx];

        cs2dArr[centerIdx][assignNum] = i;
        assignNumArr[centerIdx] += 1;

        if (assignNumArr[centerIdx] > sK){
            perror("Too many spins are assigned in center!!\n");
            exit(EXIT_FAILURE);
        }
    } 
    //for (int i=0; i<ncenter; i++){
    //    for (int j=0; j<sK; j++){
    //        printf("cs2dArr[%d][%d] = %d\n", i, j, cs2dArr[i][j]);
    //    }
    //}
    //printf("! ======================= !\n\n");
    free(assignNumArr);
    return cs2dArr;
}

void print_clusinfo_sorder_i(int*** pcceClusInfo, int sorder, int n_cInfo){
    for (int n_cinfo=0; n_cinfo<n_cInfo; n_cinfo++){
        printf("pcceClusInfo[%d][%d] = %d\n", sorder, n_cinfo, pcceClusInfo[sorder][n_cinfo]);
    }
    printf("=======\n");
}

void print_clusinfo_sorder_i_j(int*** pcceClusInfo, int sorder, int n_cinfo, int cIdx){
    //printf("=======\n");
    //printf("sorder == %d\n", sorder);
    //printf("sorder / sK = 0\n");
    //printf("=======\n");
    for (int cidx=0; cidx<cIdx; cidx++){
        printf("pcceClusInfo[%d][%d][%d] = %d\n", sorder, n_cinfo, cidx, pcceClusInfo[sorder][n_cinfo][cidx]);
    }
}

int*** convert_centerIdx_to_spinIdx(Cluster* cls, int cOrder, int sK, int** cs2dArr){
    //int*** pcceClusInfo;

    int sOrder              = cOrder * sK; // (sOrder = spin Order)  &  (cOrder = center Order <== cce.input file) //
    int*** pcceClusInfo     = (int***)allocArray1d(sOrder+1, sizeof(int**)); 

    // When cOrder = 0 // 
    pcceClusInfo[0]       = (int**)allocArray1d(1, sizeof(int*));
    pcceClusInfo[0][0]    = (int* )allocArray1d(1, sizeof(int));
    pcceClusInfo[0][0][0] = 1;
    
    // When cOrder > 0 // 
    for (int sorder = 1; sorder <= sOrder; sorder++){

        if (sorder % sK==0){
            int corder                 = int(sorder/sK);
            int n_cInfo                = cls->clusinfo[corder][0][0] - 1;               // n_clnfo: The number of "center clusinfo"'s n-th order cluster //
            pcceClusInfo[sorder]       = (int**)allocArray1d(n_cInfo+1, sizeof(int*));
            print_clusinfo_sorder_i(pcceClusInfo, sorder, n_cInfo+1);

            if (n_cInfo != 0) {
                pcceClusInfo[sorder][0]    = (int* )allocArray1d(1, sizeof(int));
                pcceClusInfo[sorder][0][0] = (n_cInfo+1);
                //printf("pcceClusInfo[%d][0][0] = %d\n", sorder, pcceClusInfo[sorder][0][0]);

                for (int n_cinfo=1; n_cinfo<(n_cInfo+1); n_cinfo++){
                    //printf("!=========================!\n");
                    //printf("sorder & n_cinfo & n_cInfo: %d & %d & %d \n\n", sorder, n_cinfo, n_cInfo);
                    pcceClusInfo[sorder][n_cinfo] = (int* )allocArray1d((sOrder+1), sizeof(int));

                    int bd = 1;
                    pcceClusInfo[sorder][n_cinfo][0] = cls->clusinfo[corder][n_cinfo][0];
                    //printf("cls      -> clustinfo[%d][%d][0] = %d\n"  , corder, n_cinfo, cls      -> clusinfo[corder][n_cinfo][0]);
                    //printf("pcceClusInfo[%d][%d][0] = %d\n\n", sorder, n_cinfo, pcceClusInfo[sorder][n_cinfo][0]);
                    for (int cidx=1; cidx<corder+1; cidx++){
                        int centerIdx = cls->clusinfo[corder][n_cinfo][cidx];
                        //printf("corder, n_cinfo, cidx: %d, %d, %d\n", corder, n_cinfo, cidx);
                        for (int sk=0; sk<sK; sk++){
                            pcceClusInfo[sorder][n_cinfo][bd++] = cs2dArr[centerIdx][sk];
                        }
                    }
                    //print_clusinfo_sorder_i_j(pcceClusInfo, sorder, n_cinfo, sOrder+1);
                    //printf("\n!=========================!\n");
                    //printf("AFTER SETTING nsij!!\n");
                    //printf("---------------------\n");
                }
            }

            else {
                pcceClusInfo[sorder][0]      = (int* )allocArray1d(1, sizeof(int));
                pcceClusInfo[sorder][0][0]   = (n_cInfo+1);
            }
        }

        else {
            pcceClusInfo[sorder]       = (int**)allocArray1d(1,sizeof(int*));
            pcceClusInfo[sorder][0]    = (int*)allocArray1d(1,sizeof(int));
            pcceClusInfo[sorder][0][0] = 1;

            //printf("=======\n");
            //printf("sorder == %d\n", sorder);
            //printf("sorder !~ sK != 0\n");
            //printf("pcceClusInfo[sorder][0][0] = %d\n", pcceClusInfo[sorder][0][0]);
            //printf("=======\n");
        }
    }
    printf("!! FINSIH: Function \"Convert center indices => spin indices. \" !!\n");
    return pcceClusInfo;
}
