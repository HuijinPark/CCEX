#ifndef __CCEX_CLUSTER_PCCE_H_
#define __CCEX_CLUSTER_PCCE_H_

#include <stdio.h>
#include <time.h>

#include "physical/bath.h"
#include "mpi.h"
#include "clustering/cluster.h"
#include "input/general.h"
#include <cstdint>
#include <random>
#include "clustering/pcce_cycle.h"

// One independent, reproducible stream per global trial. It never consumes the
// C RNG used for bath states. Explicit sampling avoids library-specific uniform
// distribution implementations; the mt19937_64 engine itself is standardized.
struct PcceRandom {
    std::mt19937_64 engine;
    PcceRandom(int seed, int trial)
        : engine((uint64_t(uint32_t(seed))<<32) | uint32_t(trial)) {}
    int index(int n){
        uint64_t bound=uint64_t(n), threshold=-bound%bound, value;
        do { value=engine(); } while (value<threshold);
        return int(value%bound);
    }
    double unit(){ return (engine()>>11)*0x1.0p-53; }
};

typedef struct{
    double x, y, z;
} Point;

struct PcceTrialWorkspace {
    std::vector<int> counts, active;
    // Scratch of assignCentroidsGrid: centers bucketed by cell, the open centers.
    std::vector<int> cell_start, cell_open, cell_items, center_cell, center_slot, open, open_slot;
    PcceCycleTracker cycle;
    PcceTrialWorkspace(int nspin,int ncenter)
        : counts(ncenter), active(ncenter),
          cycle(size_t(ncenter)*sizeof(Point),size_t(nspin)*sizeof(int)) {}
};

struct PcceTrialResult {
    int iterations=0, executed=0, skipped=0, cycle_period=0;
    bool converged=false;
};

PcceTrialResult pcceOptimizeTrial(Point* spins,Point* centers,int* assigned,
                                  int nspin,int ncenter,int max_iter,PcceTrialWorkspace& workspace);


typedef struct{
    int sK;
    int ncenter;
    int pcce_nspin;
    int rest_nspin;
} Partition_info;

void clusterizePcce(Cluster* cls, BathArray* ba, Config* config);
void simulator_cluster_partition(BathSpin** bath, Partition_info* pinfo, Cluster* cls, Point* best_centers, int* best_assigned_idx, const Config* config);
void shrink_restspins(BathArray* ba, int rest_nspin);
void update_global_best(int index, double *inertia, double *sil, int *trial, double *global_best_inertia, double *global_best_sil, int *global_best_trial);
bool pcce_better_partition(double inertia, double silhouette, double best_inertia, double best_silhouette, bool deterministic=true);
void print_pcce_info(int bath_nspin, int nqubit, int ncenter, int pcce_nspin, int rest_nspin);
void initializeCentroids(Point* spins, Point* centers, int ncenter, int pcce_nspin, PcceRandom* rng=nullptr);
void choose_initial_method(bool kmeans_pp, Point* spins, Point* centers, int pcce_nspin, int ncenter, PcceRandom* rng=nullptr);
void kMeansPlusPlus(Point *spins, Point *centers, int pcce_nspin, int ncenter, PcceRandom* rng=nullptr);
void shuffle_spins(Point *array, int* shuffle_idxs, int pcce_nspin, PcceRandom* rng=nullptr);
void changeIdx_inverse_shuffle_idx(int* shuffle_idx, int* inverse_shuffle_idx, int pcce_nspin);
void assignCentroids(int* assigned_cluster_idx, double** distances, int pcce_nspin, int ncenter) ;
// Scratch arrays each contain ncenter integers and may be reused between calls.
void assignCentroidsFromPoints(int* assigned, const Point* spins, const Point* centers,
                               int pcce_nspin, int ncenter, int* counts, int* active);
// The same assignment as assignCentroidsFromPoints (A22), looking only at the centers
// of nearby grid cells; inputs it does not cover go to assignCentroidsFromPoints.
void assignCentroidsGrid(int* assigned, const Point* spins, const Point* centers,
                         int pcce_nspin, int ncenter, PcceTrialWorkspace& workspace);
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
//int compare_dist(const void *a, const void *b);
int if_contain(int* arr, int arr_size, int value) ;
bool updateCentroids(Point *spins, int *assignments, Point *centers, int pcce_nspin, int ncenter) ;


void BathArray_connectivity_pcce(int*** cmap, float*** stmap, Point* centers, float rdip, float rdipcut, int pcce_nspin);
void BathArray_connectivity_pcce_spins(int*** cmap, float*** stmap, BathArray* ba, int* assigned, int ncenter, float rdip, float rdipcut);
int** Cluster_setCenterIdx_spinIdx_2dArr(int ncenter, int sK, int pcce_nspin, int* best_assigned_idx);
void print_clusinfo_sorder_i(int*** pcceClusInfo, int sorder, int n_cInfo);
void print_clusinfo_sorder_i_j(int*** pcceClusInfo, int sorder, int n_cinfo, int cIdx);
int*** convert_centerIdx_to_spinIdx(Cluster* cls, int cOrder, int sK, int** cs2dArr);
#endif // __CCEX_CLUSTER_PCCE_H_
