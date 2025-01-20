#ifndef __CCEX_CLUSTER_PCCE_H_
#define __CCEX_CLUSTER_PCCE_H_

#include <stdio.h>
#include <time.h>

#include "bath.h"
#include "mpi.h"
#include "cluster.h"
#include "general.h"
#include "qubit.h"


typedef struct{
    double x, y, z;
} Point;


typedef struct{
    int sK;
    int ncenter;
    int pcce_nspin;
    int rest_nspin;
} Partition_info;

void clusterizePcce(Cluster* cls, BathArray* ba, QubitArray* qa, Config* config);
void simulator_cluster_partition(BathSpin** bath, Partition_info* pinfo, Cluster* cls, Point* best_centers, int* best_assigned_idx);
void update_global_best(int index, double *inertia, double *sil, int *trial, double *global_best_inertia, double *global_best_sil, int *global_best_trial);
void print_pcce_info(int bath_nspin, int nqubit, int ncenter, int pcce_nspin, int rest_nspin);
void shrink_restspins(BathArray* ba, int rest_nspin);
void set_spinfinite(BathArray* ba, QubitArray* qa, int sK, Partition_info* pinfo, int rank);
void initializeCentroids(Point* spins, Point* centers, int ncenter, int pcce_nspin) ;
void choose_initial_method(bool kmeans_pp, Point* spins, Point* centers, int pcce_nspin, int ncenter);
void kMeansPlusPlus(Point *spins, Point *centers, int pcce_nspin, int ncenter) ;
void shuffle_spins(Point *array, int* shuffle_idxs, int pcce_nspin) ;
void changeIdx_inverse_shuffle_idx(int* shuffle_idx, int* inverse_shuffle_idx, int pcce_nspin);
void assignCentroids(int* assigned_cluster_idx, double** distances, int pcce_nspin, int ncenter) ;
void calc_dist_to_centers(Point* spins, Point* centers, double** distances, int pcce_nspin, int ncenter) ;
void copy_best_centers(Point* centers, Point* best_centers, int ncenter);
void copy_best_assigned(int* assigned, int* best_assigned, int pcce_nspin);
void copy_best_idxList(int* idxList, int* best_idxList, int pcce_nspin);
void copy_best_spins(Point* spins, Point* best_spins, int pcce_nspin);
void copy_bath_to_point_type(BathSpin** bath, Point* spins, int pcce_nspin);
void print_BD(char* str, int num);
void print_centers(int ncenter, Point* centers, int rank);
void print_kmeans_converge_info(int trial, int iteration, int ncenter, Point** centers);
void print_optimize_kmeans_info(int trial, int iteration, double inertia, double best_inertia, double best_sil, int rank);
void print_assign_idx(int *assign_idx, int pcce_nspin, int rank) ;
void print_spins(Point *spins, int pcce_nspin, int rank);
void print_time(clock_t start, clock_t end);
void write_centers(FILE **c_savefile, Point **centers, int ncenter);
void write_spins(FILE **s_savefile, Point **spins, int** assigned_idx, int pcce_nspin);
void print_best_spins(Point **spins, int** assigned_idx, int pcce_nspin);
void print_best_centers(Point **centers, int ncenter);

double silhouetteCoefficient(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter);
double calculateInertia(Point *spins, Point *centers, int *assignments, int pcce_nspin, int ncenter);
double calc_dist(double* arr1, double* arr2);
double distance(Point p1, Point p2);
int compare_dist(const void *a, const void *b);
int if_contain(int* arr, int arr_size, int value) ;
bool updateCentroids(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter) ;


void BathArray_connectivity_pcce(int*** cmap, float*** stmap, Point* centers, float rdip, float rdipcut, int pcce_nspin);
int** Cluster_setCenterIdx_spinIdx_2dArr(int ncenter, int sK, int pcce_nspin, int* best_assigned_idx);
void print_clusinfo_sorder_i(int*** pcceClusInfo, int sorder, int n_cInfo);
void print_clusinfo_sorder_i_j(int*** pcceClusInfo, int sorder, int n_cinfo, int cIdx);
int*** convert_centerIdx_to_spinIdx(Cluster* cls, int cOrder, int sK, int** cs2dArr);
#endif // __CCEX_CLUSTER_PCCE_H_
