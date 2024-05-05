#ifndef __CCEX_CLUSTER_HASH_H_
#define __CCEX_CLUSTER_HASH_H_

#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <time.h>
#include "../zlib/uthash/src/uthash.h"

#include "cluster.h"

typedef struct {
    // id : clustered spin indeices(char)
    // spins : clustered spin
    // strength : coupling strength
    const char *id;
    int *spins;
    float strength;
    // solve tilde 
    //  [0][0][0] : None
    //  [i][0][0] : The number of subclusters
    //  [i][j][0] : each subcluster j, length = order
    //  [i][j][k] : k-th spin in j-th subcluster 
    //int*** subclusters_connected;
    //int*** subclusters_notConnected;
    
    // count : The number that you have to multi. or div.
    int count;
    
    UT_hash_handle hh;   
} Property; 

typedef struct {
    //N : cluster order
    //prop : cluster properties
    const char* N;
    Property *prop; 
    UT_hash_handle hh;   
} HashCluster;

void clusterizeHash(Cluster* CCE, int nspin, int** spmap, float** stmap);
void convertClusinfoToHash(HashCluster** hashClusters, Cluster* CCE);
void solveTilde(HashCluster** hashcluster, Cluster* CCE, int nspin);

void makeHashClusterO1(HashCluster** hashclusters, int nspin);
void makeHashClusterO2(HashCluster** hashclusters, int nspin, int** spmap, float** stmap);
void makeHashClusterOn(HashCluster** hashclusters, int order, int nspin, int** spmap, float** stmap, int* nks); //nks = CCE->nk

void freeHashCluster(HashCluster** hashClusters, int order);

int addCluster(HashCluster** hashClusters, int order, const char* id,int* spins, float strength,int count);
int by_strength(const Property* a, const Property* b);
HashCluster* findCluster(HashCluster* hashClusters, int order);

int setMaxLengthStr(int nSpin, int order);

void addSpin(int** newcluster, int* oldcluster, int oldn, int spin);
void typeStr(char** destination, int* cluster, int order, int nSpin);

int binarySearch(int* arr, int low, int high, int target);

float minStrength(float currentStrength, float** strengthMap, int* oldcluster, int oldn, int newspin);
float addAllStrength(float currentStrength, float** strengthMap, int* oldcluster, int oldn, int newspin);

void printClusters(HashCluster* hashClusters);
void printProperties(Property* hashProperties);

int* parseClusterIdToIntArray(const char* id, int* count);
int countDigits(int number);
int search2dArr(int** arr, int low, int high, int left, int  right, int* target);

void updateNk(int** Nk, int order, HashCluster* hashClusters);
void addSubClusters(HashCluster** hashclusters, int nspin, float** stmap, int n);
void generateCombinations(int* arr, int*** data, int* tempCombination, int start, int end, int index, int r, int* nCombination);
float addAllStrengthForAllSpins(float** strengthMap, int* cluster, int order);

#endif // __CCEX_CLUSTER_HASH_H_
