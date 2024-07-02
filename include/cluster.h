#ifndef __CCEX_CLUSTER_H_
#define __CCEX_CLUSTER_H_

#include "general.h"
#include "bath.h"
#include "utilities.h"

/**
 * @struct ClusterConfig
 * @brief  This structure contains clustering algorithm-related parameters.
 * @details Each parameter is read from the input file or options.
 *          We will refine the parameters later to use in the simulation.
 *          Those parameters will be converted into "struct Cluster" 
 *          (See directory "clustering").
*/
typedef struct {

    /**
     * @brief The order of CCE (default : 0)
     * @details Maximum nuclear spins included in a cluster (order >= 0)
    */
    int order;

    /** pcce

    */
    int  sK;
    int  max_trial;
    int  max_iter;

    bool kmeans_pp;
    bool iter_detail;

    /**
     * @brief Clustering algorithm : cce (default) | gcce | dsj | itb | dsjitb | kmeans
     * @details The following is the description of each algorithm.
     * 
     *  - cce    : Conventional CCE method *using hash algorithm
     *  - gcce   : Generalized CCE method *using hash algorithm
     *  - dsj    : Disjoint clustering algorithm
     *             , meaning that it consider the disjointed cluster
     *  - itb    : Inter-bathcluster clustering algorithm
     *             , meaning that it consider the pair interaction between disjointed bath cluster
     *  - dsjitb : Disjoint + inter-bathcluster clustering algorithm
     *             , meaning that it consider the disjointed group 
     *             and the interaction between disjointed clusters.
     *             (dsjitb with order = 2 gives the same result as normal CCE2)
     * 
    */
    char method[MAX_CHARARRAY_LENGTH];

    /**
     * @brief The number of clusters for each order "k"
     * @details 
     * 
     *  nk is useful when you calculate CCE with order >= 3
     *  At order >=3, Creating clusters is based on rdip, a ton of cluster will be formed (at least 1e6).      
     *  This value would consider nk # clusters for each order after sorting this clusters by coupling-strength.
     *  ( NOTE : this parameter is only available at "clusalgo = hash" )
     * 
     * Input file format :
     *   nk = k:n (0<=k<=order, 0<=n<=Maximum number of clusters)
     *   e.g. nk = { 1:1000, 2:100, 3:10, 4:5 }
     *        nk = { 2:100000, 3:100000 }
     *        nk = { 4:0 } 
     *        (if n=0, it means we will include all clusters at order 4)
     * 
     * Array format: 
     *   nk[order + 1] = {0,} (default, including all clusters for each order)
     *   nk[i=0] = The number of clusters at order 0 (e.i. nk[0] = 1 always)
     *   nk[i>i] = The number of clusters at order i
    */
    int* nk; 

    /**
     * @brief Include all sub-clusters of the highest order : on | off
     * @details 
     * 
     *  When calcuting the coherence in electron spin bath, 
     *  CCE with order >= 3 gives lots of diverged point in coherence fucntion.
     *  If turn on this option, 
     *  it will include all sub-clusters of the highest order 
     *  and give the more robust result.
     *  (default : off)
    */
    bool addsubclus; 

    /**
     * @brief The clusters for each order "k"
     * @details 
     * 
     * clusters[i][j][k]
     * i : i-th order
     * j : (if j==0) The number of clusters at order i
     *     (if j!=0) j-th clusters at order i
     * k : (if k==0) The number that you will have to 
     *                multiply(+)/divide(-) in the coherence calcualtion
     *     (if k!=0) The index of the spin in the cluster
     * 
     * e.g.
     * clusters[i=0][j=0][k=0] :
     *      The number how many you should calculate in 0-th order
     * clusters[i=1][j=0][k=0] :
     *       The number of 1st order clusters
     * clusters[i=1][j!=0][k=0] :
     *       The number how many you should calculate 
     *       the coherence for 1st order and j-th cluster
     * clusters[i=1][j!=0][k!=0] :
     *       The spin index in j-th cluster at 1st order
     * 
    */
    int*** clusinfo; 

} Cluster;

/* High Level --------------------------------------------------------*/



// init
Cluster* Cluster_init();
void Cluster_clusterize(Cluster* cls, BathArray* ba, QubitArray* qa, Config* config);
void Cluster_clusterize_pcce(Cluster* cls, BathArray* ba, QubitArray* qa, Config* config, int rank);

// free
void Cluster_freeAll(Cluster* cls);

// MPI
// int*** Clsuter_getLocalClusters_MPI(Cluster* cls);

// find clusters
// void Cluster_find_clusinfo_hash(Cluster* cls, BathArray* batharr, Config* config); // normal (hash) algorithm
// void Cluster_find_clusinfo_dsj(Cluster* cls, BathArray* batharr, Config* config); // disjoint algorithm
// void Cluster_find_clusinfo_itb(Cluster* cls, BathArray* batharr, Config* config); // inter-bathcluster algorithm
// void Cluster_find_clusinfo_dsjitb(Cluster* cls, BathArray* batharr, Config* config); // disjoint + inter-bathcluster algorithm

// report
void Cluster_report(Cluster* cls);
void Cluster_reportNk(Cluster* cls);
void Cluster_reportClusinfo(Cluster* cls);
void reportClusinfo(int*** clusinfo, int order);

/* Low Level --------------------------------------------------------*/

// clusinfo controller

void Cluster_allocClusinfo(Cluster* cls, int order); // alloc clusinfo[order+1] : 0 ~ order
void Cluster_freeClusinfo(Cluster* cls);

int Cluster_setClusinfo_addcluster(Cluster* cls, int order, int iter, int* cluster); // add new cluster, return the index of the cluster
void Cluster_setClusinfo_chgcluster(Cluster* cls, int order, int* cluster, int ic); // change cluster at ic-th cluster
void Cluster_setClusinfo_chgiter(Cluster* cls, int order, int iter, int ic); // change iter at ic-th cluster 

//void Cluster_setSk(Cluster* cls, int sk); //
//void Cluster_getSk(Cluster* cls); // output



// what's iter :     
// The number how many you will multiply/divide 
// in the coherence calculation for 0 th order

int Cluster_getClusinfo_ncluster(Cluster* cls, int order); // get the number of clusters at order
int* Cluster_getClusinfo_itercluster(Cluster* cls, int order, int ic); // get ic-th iter+cluster at order [0] = iter, [1~] = cluster
int Cluster_getClusinfo_iter(Cluster* cls, int order, int ic); // get ic-th iter+cluster at order [0] = iter, [1~] = cluster
int* Cluster_getClusinfo_cluster_copy(Cluster* cls, int order, int ic); // get ic-th cluster at order (copy)
int*** Cluster_getClusinfo(Cluster* cls);

// nk controller
void Cluster_allocNk(Cluster* cls);
void Cluster_freeNk(Cluster* qa);
int Cluster_getNk_order(Cluster* cls, int i); // get the number of clusters at ith-orderr at ith-order
int* Cluster_getNk(Cluster* cls);
int Cluster_getSk(Cluster* cls);
int Cluster_getMax_iter(Cluster* cls);
int Cluster_getMax_trial(Cluster* cls);
bool Cluster_getKmeans_pp(Cluster* cls);
bool Cluster_getIter_detail(Cluster* cls);

// others controller
char* Cluster_getMethod(Cluster* cls);
int Cluster_getOrder(Cluster* cls);
bool Cluster_getAddsubclus(Cluster* cls);

// set 
void Cluster_setSk(Cluster* cls,          int sK);
void Cluster_setMax_trial(Cluster* cls,   int max_trial);
void Cluster_setMax_iter(Cluster* cls,    int max_iter);
void Cluster_setKmeans_pp(Cluster* cls,   bool kmeans_pp);
void Cluster_setIter_detail(Cluster* cls, bool iter_detail);
void Cluster_setOrder(Cluster* cls, int order);
void Cluster_setMethod(Cluster* cls, char* method);
void Cluster_setAddsubclus(Cluster* cls, bool addsubclus);
void Cluster_setNk(Cluster* cls, int* nk); // set the number of clusters at ith-order


#endif // __CCEX_CLUSTER_H_
