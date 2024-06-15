#ifndef __CCEX_CLUSTER_PCCE_H_
#define __CCEX_CLUSTER_PCCE_H_

#include <stdio.h>
#include <time.h>

#include "bath.h"
#include "mpi.h"
#include "cluster.h"
#include "general.h"
#include "qubit.h"


//typedef struct{
//    double x, y, z;
//} Point;

typedef struct{
    bool kmeans_pp;
    bool iter_detail;
    double rBath;
    int sK;
    int max_trial;
    int max_iteration;
    char bathfile[100];
    char outfile[100];
} InputData;

typedef struct{
    int sK;
    int npartition;
    int partition_nspin;
    int rest_nspin;
} Partition_info;

void update_global_best(int index, double *inertia, double *sil, int *trial, double *global_best_inertia, double *global_best_sil, int *global_best_trial);
void print_partition_info(int bath_nspin, int nqubit, int npartition, int partition_nspin, int rest_nspin);
void shrink_restspins(BathArray* ba, int rest_nspin);
void set_spinfinite(BathArray* ba, QubitArray* qa, int sK, Partition_info* pinfo);
void initializeCentroids(Point* spins, Point* centers, int ClusterNum, int TotalSpinNum) ;
void choose_initial_method(bool kmeans_pp, Point* spins, Point* centers, int TotalSpinNum, int ClusterNum);
void kMeansPlusPlus(Point *spins, Point *centers, int TotalSpinNum, int ClusterNum) ;
void shuffle_spins(Point *array, int* shuffle_idxs, int TotalSpinNum) ;
void changeIdx_inverse_shuffle_idx(int* shuffle_idx, int* inverse_shuffle_idx, int TotalSpinNum);
void assignCentroids(int* assigned_cluster_idx, double** distances, int TotalSpinNum, int ClusterNum) ;
void calc_dist_to_centers(Point* spins, Point* centers, double** distances, int TotalSpinNum, int ClusterNum) ;
void copy_best_centers(Point* centers, Point* best_centers, int ClusterNum);
void copy_best_assigned(int* assigned, int* best_assigned, int TotalSpinNum);
void copy_best_idxList(int* idxList, int* best_idxList, int TotalSpinNum);
void copy_best_spins(Point* spins, Point* best_spins, int TotalSpinNum);

void copy_best_config(Point*  centers     , Point* best_centers  , 
                      Point*  spins       , Point*   best_spins  , 
                      int*    assigned    , int*    best_assigned, 
                      int     TotalSpinNum, int     ClusterNum    );

void copy_bath_to_point_type(BathSpin** bath, Point* spins, int TotalSpinNum);
void print_BD(char* str, int num);
void print_calc_detail(InputData inputdata);
void print_Partition_info(int ClusterNum, int restNum, int TotalSpinNum, int sK, int spinNum);
void print_centers(int ClusterNum, Point* centers, int rank);
void print_kmeans_converge_info(int trial, int iteration, int ClusterNum, Point** centers);
void print_optimize_kmeans_info(int trial, int iteration, double inertia, double best_inertia, double best_sil, int rank);
void print_assign_idx(int *assign_idx, int TotalSpinNum, int rank) ;
void print_spins(Point *spins, int TotalSpinNum, int rank);
void print_time(clock_t start, clock_t end);
void write_centers(FILE **c_savefile, Point **centers, int ClusterNum);
void write_spins(FILE **s_savefile, Point **spins, int** assigned_idx, int TotalSpinNum);
//void simulator_cluster_partition(BathArray* ba, Partition_info* pinfo, Cluster* cls);
void simulator_cluster_partition(BathSpin** bath, Partition_info* pinfo, Cluster* cls, Point* best_centers, int* best_assigned_idx);

double silhouetteCoefficient(Point *spins, int *assignments, Point *centers, int TotalSpinNum, int ClusterNum);
double calculateInertia(Point *spins, Point *centers, int *assignments, int TotalSpinNum, int ClusterNum);
double calc_dist(double* arr1, double* arr2);
double distance(Point p1, Point p2);
int compare_dist(const void *a, const void *b);
int if_contain(int* arr, int arr_size, int value) ;
bool updateCentroids(Point *spins, int *assignments, Point *centers, int TotalSpinNum, int ClusterNum) ;

#endif // __CCEX_CLUSTER_PCCE_H_
