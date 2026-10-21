// Partition invariants and numerical edge cases; no MPI runtime is needed.
#include "clustering/cluster_pcce.h"
#include "clustering/cluster_hash.h"
#include "common/memory.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

int rank=0, nprocess=1, seed=19;
bool verbosity=false;

int main(int argc, char** argv){
    if (argc>1 && !strcmp(argv[1],"snapshot")){
        // Seeds, centers, and the following RNG value detect changes to both
        // the partition initialization and subsequent sampled bath states.
        std::vector<Point> spins(96), centers(96);
        for (int i=0;i<96;++i) spins[i]={double((i*17)%97),double((i*29)%89),double((i*43)%83)};
        for (int count : {1,2,8,48,96}) for (int initialSeed : {1,19,123,12345}){
            srand(initialSeed);
            kMeansPlusPlus(spins.data(),centers.data(),96,count);
            printf("%d %d",count,initialSeed);
            for (int i=0;i<count;++i) printf(" %.17g %.17g %.17g",centers[i].x,centers[i].y,centers[i].z);
            printf(" %d\n",rand());
        }
        return 0;
    }
    Point spins[]={{0,0,0},{2,0,0},{10,0,0},{12,0,0}};
    Point centers[]={{1,0,0},{11,0,0}};
    int pairs[]={0,0,1,1}, singletons[]={0,1,2,3}, one[]={0,0,0,0};
    double normal=silhouetteCoefficient(spins,pairs,centers,4,2);
    double singleton=silhouetteCoefficient(spins,singletons,spins,4,4);
    double only=silhouetteCoefficient(spins,one,centers,4,1);
    Point coincident[4]={};
    double zero=silhouetteCoefficient(coincident,pairs,centers,4,2);
    printf("silhouette normal=%.17g singleton=%.17g one_group=%.17g coincident=%.17g\n",normal,singleton,only,zero);
    fflush(stdout);
    assert(fabs(normal-0.7979797979797979)<1e-14);
    assert(singleton==0 && only==0 && zero==0);
    assert(pcce_better_partition(1,0.1,2,0.9));
    assert(pcce_better_partition(1,0.9,1,0.1));
    assert(!pcce_better_partition(1,0.1,1,0.9));
    assert(!pcce_better_partition(1,0.9,1,0.9));

    assert(pcce_better_partition(1,0.1,1,0.9,false));
    assert(!pcce_better_partition(1,0.9,1,0.1,false));
    assert(pcce_better_partition(1,0.9,2,0.1,false));

    // An absent highest order must not prevent lower-order cluster products.
    for (int connected : {0,1}){
        Cluster* cls=Cluster_init();
        Cluster_setOrder(cls,4);
        Cluster_allocNk(cls);
        Cluster_allocClusinfo(cls,4);
        cls->addsubclus=true;
        int** sparse=allocInt2d(3,3);
        float** strengths=allocFloat2d(3,3);
        for (int i=0;i<3;++i){
            sparse[i][0]=1;
            for (int j=0;j<3;++j) if (connected && i!=j){
                sparse[i][sparse[i][0]++]=j;
                strengths[i][j]=1;
            }
        }
        clusterizeHash(cls,3,sparse,strengths);
        assert(Cluster_getClusinfo_ncluster(cls,4)==0);
        if (connected){
            assert(Cluster_getClusinfo_ncluster(cls,3)==1);
            assert(Cluster_getClusinfo_iter(cls,3,1)==1);
            assert(Cluster_getClusinfo_ncluster(cls,1)==0);
        }else{
            assert(Cluster_getClusinfo_ncluster(cls,1)==3);
            assert(Cluster_getClusinfo_ncluster(cls,2)==0);
        }
        freeInt2d(&sparse,3);
        freeFloat2d(&strengths,3);
        Cluster_freeAll(cls);
    }
    puts("PASS partition numeric edge cases and empty cluster orders");
}
