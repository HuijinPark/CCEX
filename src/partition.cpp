#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "../include/partition.h"
#include "../include/qubit.h"
#include "../include/cluster.h"
#include "../include/bath.h"

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

void print_partition_info(int bath_nspin, int nqubit, int npartition, int partition_nspin, int rest_nspin){
    print_BD("-", 30);
    printf("Total Bath spin: %d\n", bath_nspin);
    printf("Qubit #: %d\n", nqubit);
    printf("npartition: %d\n", npartition);
    printf("partition_nspin: %d\n", partition_nspin);
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

void set_spinfinite(BathArray* ba, QubitArray* qa, int sK, Partition_info* pinfo){
    
    int npartition      = (ba->nspin / sK);
    int partition_nspin = (npartition * sK);
    int rest_nspin      = (ba->nspin % sK);

    pinfo->sK = sK;
    pinfo->npartition = npartition;
    pinfo->partition_nspin = partition_nspin;
    pinfo->rest_nspin = rest_nspin;
    
    print_partition_info(ba->nspin, qa->nqubit, npartition, partition_nspin, rest_nspin);

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
void initializeCentroids(Point* spins, Point* centers, int ClusterNum, int TotalSpinNum) {
    
    // Make and initialize the Random number index list // 
    int count = 0;
    int* RandomList = (int*) calloc(ClusterNum, sizeof(int));
    for (int i=0; i<ClusterNum; i++){
        RandomList[i] = (TotalSpinNum+100);
    }
    
    // Random vales => srand //
    // srand(time(NULL)+i*i);
    for (int count=0; count <ClusterNum;) {
        int random_index = rand() % TotalSpinNum;
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

void choose_initial_method(bool kmeans_pp, Point* spins, Point* centers, int TotalSpinNum, int ClusterNum){
    //printf("||------------------- Start Iter. -------------------||\n\n");
    if (kmeans_pp == true){
        //printf(" - Centroids Initialization: K-Means++ \n\n");
        //printf(" || Centroids Initialization: K-Means++ ||\n\n");
        //print_BD("=", 55);
        kMeansPlusPlus(spins, centers, TotalSpinNum, ClusterNum);
    }
    else {
        //printf(" - Centroids Initialization: Random \n");
        //printf(" || Centroids Initialization: Random ||\n\n");
        //print_BD("=", 55);
        initializeCentroids(spins, centers, ClusterNum, TotalSpinNum);
    }
    printf("\n");
}

void kMeansPlusPlus(Point *spins, Point *centers, int TotalSpinNum, int ClusterNum) {

    int index = rand() % TotalSpinNum;  // 무작위 인덱스
    
    centers[0].x = spins[index].x;
    centers[0].y = spins[index].y;
    centers[0].z = spins[index].z;
    //centers[0] = spins[index];  // 첫 번째 중심점 선택

    double *dists = (double*) malloc(TotalSpinNum * sizeof(double));

    for (int i = 1; i < ClusterNum; i++) {
        double sum = 0.0;
        // 각 포인트에 대한 최소 거리를 제곱하여 계산
        for (int j = 0; j < TotalSpinNum; j++) {
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
        for (int j = 0; j < TotalSpinNum; j++) {
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
void shuffle_spins(Point *array, int* shuffle_idxs, int TotalSpinNum) {
    // 배열의 마지막 요소부터 시작하여 첫 번째 요소 이전까지 반복
    for (int i = TotalSpinNum - 1; i > 0; i--) {
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

void changeIdx_inverse_shuffle_idx(int* shuffle_idx, int* inverse_shuffle_idx, int TotalSpinNum){
    for (int i=0; i<TotalSpinNum; i++){
        int j = shuffle_idx[i];
        inverse_shuffle_idx[j] = i;
    }
}

void assignCentroids(int* assigned_cluster_idx, double** distances, int TotalSpinNum, int ClusterNum) {
    int spinsPerCenter = TotalSpinNum / ClusterNum;
    int *centroidCounts = (int *)malloc(ClusterNum * sizeof(int));
    for (int i = 0; i < ClusterNum; i++) {
        centroidCounts[i] = 0;
    }

    for (int i = 0; i < TotalSpinNum; i++) {
        int minIndex = -1;
        double minDistance = HUGE_VAL;
        for (int j = 0; j < ClusterNum; j++) {
            if (distances[i][j] < minDistance && centroidCounts[j] < spinsPerCenter) {
                minDistance = distances[i][j];
                minIndex = j;
            }
        }
        if (minIndex == -1) {
            // It may happen that no centroid can take more points, distribute remaining points
            for (int j = 0; j < ClusterNum; j++) {
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

void calc_dist_to_centers(Point* spins, Point* centers, double** distances, int TotalSpinNum, int ClusterNum) {
    for (int i = 0; i < TotalSpinNum; i++) {
        for (int j = 0; j < ClusterNum; j++) {
            double dx = spins[i].x - centers[j].x;
            double dy = spins[i].y - centers[j].y;
            double dz = spins[i].z - centers[j].z;
            //dists[i][j] = sqrt(dx * dx + dy * dy + dz * dz);
            distances[i][j] = (dx * dx + dy * dy + dz * dz);
            //printf("|| DISTANCE[%d][%d] =%lf || \n", i, j, (dx * dx + dy * dy + dz * dz));
        }
    }
}

bool updateCentroids(Point *spins, int *assignments, Point *centers, int TotalSpinNum, int ClusterNum) {

    bool converged = true;
    int *counts = (int *)malloc(ClusterNum * sizeof(int));
    Point *sums = (Point *)malloc(ClusterNum * sizeof(Point));
    
    // Initialize //
    for (int i = 0; i < ClusterNum; i++) {
        sums[i].x = sums[i].y = sums[i].z = 0;
        counts[i] = 0;
    }

    for (int i = 0; i < TotalSpinNum; i++) {
        int cid = assignments[i];
        sums[cid].x += spins[i].x;
        sums[cid].y += spins[i].y;
        sums[cid].z += spins[i].z;
        counts[cid]++;
    }

    double TOLERANCE = 0.0001;
    for (int i = 0; i < ClusterNum; i++) {
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

double calculateInertia(Point *spins, Point *centers, int *assignments, int TotalSpinNum, int ClusterNum) {
    double inertia = 0.0;
    for (int i = 0; i < TotalSpinNum; i++) {
        int cid = assignments[i];
        double dx = spins[i].x - centers[cid].x;
        double dy = spins[i].y - centers[cid].y;
        double dz = spins[i].z - centers[cid].z;
        inertia += dx * dx + dy * dy + dz * dz;
    }
    return inertia;
}

double silhouetteCoefficient(Point *spins, int *assignments, Point *centers, int TotalSpinNum, int ClusterNum) {

    double *a = (double*)malloc(TotalSpinNum * sizeof(double)); // 클러스터 내 거리
    double *b = (double*)malloc(TotalSpinNum * sizeof(double)); // 클러스터 간 거리
    double totalSilhouette = 0.0; // 전체 실루엣 계수 합

    // 모든 포인트에 대한 클러스터 내 거리 a(i) 초기화
    for (int i = 0; i < TotalSpinNum; i++) {
        double sumDist = 0.0;
        int clusterSize = 0;
        for (int j = 0; j < TotalSpinNum; j++) {
            if (assignments[i] == assignments[j]) {
                sumDist += distance(spins[i], spins[j]);
                clusterSize++;
            }
        }
        a[i] = sumDist / (clusterSize - 1); // 자기 자신을 제외
    }

    // 모든 포인트에 대한 클러스터 간 최소 거리 b(i) 계산
    for (int i = 0; i < TotalSpinNum; i++) {
        double minDist = FLT_MAX;
        for (int c = 0; c < ClusterNum; c++) {
            if (c != assignments[i]) { // 다른 클러스터
                double clusterDist = 0.0;
                int clusterSize = 0;
                for (int j = 0; j < TotalSpinNum; j++) {
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
    for (int i = 0; i < TotalSpinNum; i++) {
        double si = (b[i] - a[i]) / (a[i] > b[i] ? a[i] : b[i]);
        totalSilhouette += si;
    }
    double averageSilhouette = totalSilhouette / TotalSpinNum;
    free(a);
    free(b);

    return averageSilhouette;
}

double distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y) + (p1.z - p2.z) * (p1.z - p2.z));
}

/////////////////////////////////////////////////////////////////////////////
// Function: Copy now array --> best array //
void copy_best_centers(Point* centers, Point* best_centers, int ClusterNum){
    for (int i=0; i<ClusterNum; i++){
        best_centers[i] = centers[i];
    }
}

void copy_best_assigned(int* assigned, int* best_assigned, int TotalSpinNum){
    for (int i=0; i<TotalSpinNum; i++){
        best_assigned[i] = assigned[i];
    }
}

void copy_best_idxList(int* idxList, int* best_idxList, int TotalSpinNum){
    for (int i=0; i<TotalSpinNum; i++){
        best_idxList[i] = idxList[i];
    }
}

void copy_best_spins(Point* spins, Point* best_spins, int TotalSpinNum){
    for (int i=0; i<TotalSpinNum; i++){
        best_spins[i] = spins[i];
    }
}

void copy_best_config(Point*  centers     , Point* best_centers  , 
                      Point*  spins       , Point*   best_spins  , 
                      int*    assigned    , int*    best_assigned, 
                      int     TotalSpinNum, int     ClusterNum    ) {

    for (int i=0; i<ClusterNum; i++){
        best_centers[i] = centers[i];
    }
    for (int i=0; i<TotalSpinNum; i++){
        best_assigned[i] = assigned[i];
    }
    for (int i=0; i<TotalSpinNum; i++){
        best_spins[i] = spins[i];
    }

}

void copy_bath_to_point_type(BathSpin** bath, Point* spins, int TotalSpinNum){
    for (int i=0; i<TotalSpinNum; i++){
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
void print_calc_detail(InputData inputdata){

    print_BD("=", 55);
    printf("\n||---------------- Calculation Detail ----------------||\n\n");
    printf("kmeans_pp   = %s\n", inputdata.kmeans_pp ? "true" : "false");
    printf("iter_detail = %s\n", inputdata.iter_detail ? "true" : "false");
    printf("rBath       = %lf\n", inputdata.rBath);
    printf("sK          = %d\n",  inputdata.sK);
    printf("bathfile    = %s\n",  inputdata.bathfile);
    printf("outfile     = %s\n",  inputdata.outfile);
    printf("max_trial   = %d\n",  inputdata.max_trial);
    printf("max_iter.   = %d\n",  inputdata.max_iteration);
    print_BD("=", 55);

}

///////////////////////////////////
void print_Partition_info(int ClusterNum, int restNum, int TotalSpinNum, int sK, int spinNum){
    //printf("\n");
    print_BD("-", 55);
    printf(" || Cluster Infomation: || \n\n");
    printf("   - Point # (in rBath ) = %d\n", spinNum);      printf("\n");
    printf("   - Point # (calc. )    = %d\n", TotalSpinNum); printf("\n");
    printf("   - Total Point #       = %d\n", TotalSpinNum); printf("\n");
    printf("   - sK                 = %d\n", sK);           printf("\n");
    printf("   - Cluster    #       = %d\n", ClusterNum);   printf("\n");
    printf("   - restNum    #       = %d\n", restNum);      //printf("\n");
    print_BD("-", 55);
}

///////////////////////////////////
void print_centers(int ClusterNum, Point* centers, int rank){
    print_BD("-", 55);
    printf("  - Rank: %d\n\n", rank);
    for (int i=0; i < ClusterNum; i++){
        printf(" Point %5d: %10.3lf %10.3lf %10.3lf\n", i, centers[i].x, centers[i].y, centers[i].z);
    }
}

///////////////////////////////////
void print_kmeans_converge_info(int trial, int iteration, int ClusterNum, Point** centers){
    print_BD("-", 55);
    printf(" || Iter. %d - (%d) ||\n", trial, iteration);
    for (int i=0; i < ClusterNum; i++){
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
void print_assign_idx(int *assign_idx, int TotalSpinNum, int rank) {
    printf("\n || Point assignments to centers ||\n\n");
    printf(" - Rank: %d\n\n", rank);
    for (int i = 0; i < TotalSpinNum ; i++) {
        printf(" Point %5d -> Center %4d\n", i, assign_idx[i]);
    }
    printf("\n\n");
}

void print_spins(Point *spins, int TotalSpinNum, int rank){
    //printf("shuffle!!\n");
    printf(" - Rank: %d\n\n", rank);
    for (int i=0; i<TotalSpinNum; i++){
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
void write_centers(FILE **c_savefile, Point **centers, int ClusterNum){

    printf("WRITE CENTERS!!\n");
    *c_savefile = fopen("./out/centers.txt", "w");
    if (*c_savefile == NULL){
        printf("File could not be opend\n");
    }
    fprintf(*c_savefile, " x   y   z \n");
    for (int p=0; p<ClusterNum; p++){
        fprintf(*c_savefile, " %10.5lf %10.5lf %10.5lf\n", (*centers)[p].x, (*centers)[p].y, (*centers)[p].z);
    }
    fclose(*c_savefile);
}

///////////////////////////////////
void write_spins(FILE **s_savefile, Point **spins, int** assigned_idx, int TotalSpinNum){

    *s_savefile = fopen("./out/spins.txt", "w");
    if (*s_savefile == NULL){
        printf("File could not be opend\n");
    }
    fprintf(*s_savefile, " x   y   z   cidx\n");
    for (int p=0; p<TotalSpinNum; p++){
        fprintf(*s_savefile, " %10.5lf %10.5lf %10.5lf %10d\n", (*spins)[p].x, (*spins)[p].y, (*spins)[p].z, (*assigned_idx)[p]);
    }
    fclose(*s_savefile);
}

