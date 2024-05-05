#include "../include/cluster.h"
#include "../include/cluster_hash.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include <stdlib.h>

bool verbosity = true;
int rank = 0;

int main(){
    Cluster* cls = Cluster_init();

    int order = 5;
    char method[MAX_CHARARRAY_LENGTH] = "hash";
    int* nk = allocInt1d(order+1);
    for (int i=0; i<order+1; i++){
        // random
        nk[i] = rand()%100;
    }
    bool addsubclus = true;

    ////////////////////////////////////////////
    // Check set
    printf("     Set the Cluster structure\n");
    printf("\n    >> Cluster properties\n");
    Cluster_setOrder(cls,order);
    Cluster_setMethod(cls,method);
    Cluster_setAddsubclus(cls,addsubclus);
    Cluster_report(cls);

    Cluster_allocNk(cls);
    Cluster_setNk(cls,nk);
    Cluster_reportNk(cls);

    ////////////////////////////////////////////
    
    ////////////////////////////////////////////
    // Clusterize without BathArray
    ////////////////////////////////////////////

    // Random generator to make cmap, stmap
    int nspin = 10;
    int** cmap = allocInt2d(nspin,nspin);
    float** stmap = allocFloat2d(nspin,nspin);

    for (int i=0; i<nspin; i++){
        for (int j=0; j<nspin; j++){
            cmap[i][j] = rand()%2;
            if (cmap[i][j] == 1){
                stmap[i][j] = rand()%100;
            }
        }
    }
    
    // Make Sparsemap
    int** spmap = NULL;
    makeSparsemap(&spmap,cmap,nspin);

    // Hash
    printf("\n    >> Clusterizing with %s-algorithm\n\n", Cluster_getMethod(cls));
    
    order = Cluster_getOrder(cls);

    Cluster_allocClusinfo(cls,order);
    clusterizeHash(cls, nspin, spmap, stmap);

    Cluster_reportClusinfo(cls);
    Cluster_reportNk(cls);
    Cluster_freeClusinfo(cls);

    ////////////////////////////////////////////
    // Clusterize
    ////////////////////////////////////////////
    
    // printf("\n    >> Clusterizing algorithm\n\n");

    // Cluster_setMethod(cluster,"cce");
    // printf("        \"%s\" algorithm \n",Cluster_getMethod(cluster));
    // Cluster_clusterize_hash(cluster); //set 0th [0][0][0] = 0
    // Cluster_reportClusinfo(cluster);
    // Cluster_freeClusinfo(cluster);

    // Cluster_setMethod(cluster,"gcce");
    // printf("        \"%s\" algorithm \n",Cluster_getMethod(cluster));
    // Cluster_clusterize_hash(cluster); // set 0th [0][0][0] != 0
    // Cluster_reportClusinfo(cluster);
    // Cluster_freeClusinfo(cluster);

    // Cluster_setMethod(cluster,"dsj");
    // printf("        \"%s\" algorithm \n",Cluster_getMethod(cluster));
    // Cluster_clusterize_dsjitb(cluster); // ignore itb
    // Cluster_reportClusinfo(cluster);
    // Cluster_freeClusinfo(cluster);

    // Cluster_setMethod(cluster,"itb");
    // printf("        \"%s\" algorithm \n",Cluster_getMethod(cluster));
    // Cluster_clusterize_dsjitb(cluster); // ignore dsj
    // Cluster_reportClusinfo(cluster);
    // Cluster_freeClusinfo(cluster);

    // Cluster_setMethod(cluster,"dsjitb");
    // printf("        \"%s\" algorithm \n",Cluster_getMethod(cluster));
    // Cluster_clusterize_dsjitb(cluster);
    // Cluster_reportClusinfo(cluster);
    // Cluster_freeClusinfo(cluster);

    // Cluster_setMethod(cluster,"pcce");
    // printf("        \"%s\" algorithm \n",Cluster_getMethod(cluster));
    // Cluster_clusterize_pcce_kmeans(cluster);
    // Cluster_reportClusinfo(cluster);
    // Cluster_freeClusinfo(cluster);
    ////////////////////////////////////////////

    ////////////////////////////////////////////
    // Free
    ////////////////////////////////////////////
    Cluster_freeAll(cls);
    ////////////////////////////////////////////
    return 0;
}