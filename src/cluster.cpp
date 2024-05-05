#include "../include/cluster.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include "../include/cluster_hash.h"
// #include "mpi.h"

// init
// clusterize
// free
// report

/* High Level --------------------------------------------------------*/

// init
Cluster* Cluster_init(){
    Cluster* cls = (Cluster*)allocArray1d(1,sizeof(Cluster));
    Cluster_setOrder(cls, 0);
    Cluster_setAddsubclus(cls, false);
    Cluster_setNk(cls, NULL);
    cls->method[0] = '\0';
    cls->clusinfo = NULL;
    
    return cls;
}

int*** Clsuter_getLocalClusters_MPI(Cluster* cls){
    
    int order = Cluster_getOrder(cls);
    int*** clusters = Cluster_getClusinfo(cls);

    MPI_Request req; 
    MPI_Status status;

    MPI_Barrier(MPI_COMM_WORLD);
    // get the number of clusters for each order
    int MPI_size[order+1];

    for (int n=0; n<order+1; n++){
        if (n==0){
            // The number of cluster = 1#
            MPI_size[n] = 1;
        }
        else{
            // The number of cluster
            MPI_size[n] = clusters[n][0][0];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    // get sendcount, ista, iend for each order and each rank
    int MPI_sendcount[order+1][nprocess];
    int MPI_ista[order+1][nprocess];
    int MPI_iend[order+1][nprocess];

    for (int irank=0; irank<nprocess; irank++){
        for (int n=0; n<order+1; n++){
            if (n==0){
                MPI_sendcount[n][irank] = 1;
                MPI_ista[n][irank] = 0;
                MPI_iend[n][irank] = 0;
            }
            else{
                int ista, iend;
                para_range(2, MPI_size[n], nprocess, irank, &(ista) ,&(iend));
                MPI_Barrier(MPI_COMM_WORLD);
                //if (rank==0){
                //    printf("rank%d, size(ncluster)=%d, ista=%d , iend=%d\n",irank,MPI_size[n],ista,iend);
                //    printf("rank%d, sendcount%d \n",irank,iend-ista+1);
                //}
                // cluster : ista - 1 <= i < iend
                // 9# , 0 , 1 ... , 8 , rank = 0 .. 7
                // rank==0~7 then, sendcount = 1 ista=2, iend=1
                // else, sendcount = 0
                MPI_sendcount[n][irank] = iend - ista + 1;
                MPI_ista[n][irank] = ista - 1;
                MPI_iend[n][irank] = iend;
                
                // case : ista - 1 = iend
                // case : ista - 1 > iend
                if (ista-1 >= iend){
                    MPI_iend[n][irank] = MPI_ista[n][irank];
                    MPI_sendcount[n][irank] = 0;
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
    }

    //if (rank==0){
    //    printf("nprocess : %d\n",nprocess);
    //    for (int n=0; n<order+1; n++){
    //        printf("size[%d] : %d\n",n,MPI_size[n]);
    //        for (int ir=0; ir<nprocess; ir++){
    //            printf("ista[%d][%d] : %d\n",n,ir,MPI_ista[n][ir]);
    //            printf("iend[%d][%d] : %d\n",n,ir,MPI_iend[n][ir]);
    //            printf("sendcount[%d][%d] : %d\n",n,ir,MPI_sendcount[n][ir]);
    //        }
    //    }
    //}

    MPI_Barrier(MPI_COMM_WORLD);


    // make local clusters for each rank
    int*** localClusters = (int***)calloc(order+1,sizeof(int**));

    // zeroth cluster 
    localClusters[0] = (int**)calloc(1,sizeof(int*));
    localClusters[0][0] = (int*)calloc(1,sizeof(int));
    localClusters[0][0][0] = clusters[0][0][0]; // = 1
    if (rank!=0){
        localClusters[0][0][0] = 0;
    } 
    MPI_Barrier(MPI_COMM_WORLD);

    // n > 0 clusters
    for (int n=1; n<order+1; n++){

        int size = MPI_sendcount[n][rank] + 1;
        localClusters[n] = (int**)allocArray2d(size,n+1,sizeof(int));
        localClusters[n][0][0] = size;

        int rootClusterista = MPI_ista[n][rank];
        int rootClusteriend = MPI_iend[n][rank];
        int iroot = rootClusterista - 1;
        for (int i=1; i<size; i++){
            iroot = rootClusterista + i - 1;
            //printf("rank[%d] : iroot = %d\n",rank,iroot);
            for (int j=0; j<n+1; j++){
                localClusters[n][i][j] = clusters[n][iroot][j];
            }
        }

        if (iroot != rootClusteriend-1 ){
            printf("rank[%d] : iroot = %d, rootClusteriend = %d\n",rank,iroot,rootClusteriend);
            perror("iroot != rootClusteriend");
            exit(1);
        }        
    }
    MPI_Barrier(MPI_COMM_WORLD);

    return localClusters;
}

void Cluster_clusterize(Cluster* cls, BathArray* ba, Config* config){

    int order = Cluster_getOrder(cls);
    int nspin = BathArray_getNspin(ba);

    Cluster_allocClusinfo(cls, order); // clusinfo[order+1] : 0 ~ order, clusinfo[i][0][0] = 1

    if (order==0){
        if (rank==0){printf("\n\t Clustering 0th order ... \n");}
    }
    else if (order == 1){
        
        if (rank==0){printf("\n\t Clustering 1st order ... \n");}

        int cluster[1] = {0};
        int iter = 1;
        for (int spinidx = 0; spinidx < nspin; spinidx++){
            cluster[0] = spinidx;
            Cluster_setClusinfo_addcluster(cls, order, iter, cluster);
        }

    }
    else if (order > 1){ 

        int** cmap = NULL; // connectivity map : cmap[nspin][nspin] = 1 if connected else 0
        float** stmap = NULL; // strength map : stmap[nspin][nspin] = strength
        int** spmap = NULL; // sparse map : spmap[nspin][ncol] = n connected spin + 1d array of connected spins index

        float rdip = Config_getRdip(config);
        float rdipcut = Config_getRdipcut(config);

        // Connectivity Map
        BathArray_connectivity(&cmap, &stmap, ba, rdip, rdipcut);
        makeSparsemap(&spmap, cmap, nspin);

        if (strcasecmp(cls->method, "cce")==0 || strcasecmp(cls->method, "gcce")==0){
            clusterizeHash(cls, nspin, spmap, stmap); // clusterize + solve tilde
        }
        // else if (strcmp(cls->method, "dsj")==0){
        //     clusterizeDsj(cls, cmap, stmap, spmap, ba, config);
        // }
        // else if (strcmp(cls->method, "itb")==0){
        //     // clusterizeItb(cls, cmap, stmap, spmap, ba, config);
        // }
        // else if (strcmp(cls->method, "dsjitb")==0){
        //     // clusterizeDsjitb(cls, cmap, stmap, spmap, ba, config);
        // }
        // else if (strcmp(cls->method, "pcce")==0){
        //     // clusterizePcce(cls, cmap, stmap, spmap, ba, config);
        else{
            printf("Error: clusterize: method is not defined\n");
            exit(1);
        }

        freeInt2d(cmap,nspin);
        freeInt2d(spmap,nspin);
        freeFloat2d(stmap,nspin);

    }else{
        fprintf(stderr,"Error: Cluster_clusterize: order(%d) is not defined\n",order);
        exit(1);
    }

    char* method = Cluster_getMethod(cls);
    if (strcasecmp(method, "gcce")!=0 and strcasecmp(method, "cce")!=0){
        // cls -> hash
        HashCluster* hashClusters = NULL;
        convertClusinfoToHash(&hashClusters, cls);
        Cluster_freeClusinfo(cls);

        // solve tilde, hash -> cls
        solveTilde(&hashClusters, cls, nspin);
        freeHashCluster(&hashClusters, order);
    }

    if (strcasecmp(method, "cce")){
        Cluster_setClusinfo_chgiter(cls, 0, 1, 0);
    }
}


// free
void Cluster_freeAll(Cluster* cls){
    freeArray1d(cls);
}

/* Low Level --------------------------------------------------------*/

// report
void Cluster_report(Cluster* cls){

    printLineSection();
    printTitle("Structure Cluster");

    printSubTitle("General properties");
    printStructElementInt("order",Cluster_getOrder(cls));
    printStructElementChar("method",cls->method);
    printStructElementBool("addsubclus",Cluster_getAddsubclus(cls));
    
    if (cls->nk!=NULL){
        printSubTitle("The number of clusters for each order");
        Cluster_reportNk(cls);
    }
    
    if (cls->clusinfo!=NULL){
        printSubTitle("Clustering information");
        Cluster_reportClusinfo(cls);
    }

    printf("\n");
    printLineSection();
}

void Cluster_reportNk(Cluster* cls){
    int order = Cluster_getOrder(cls);
    printStructElementInt1dIdx("nk",Cluster_getNk(cls),order+1);
}

void Cluster_reportClusinfo(Cluster* cls){

    int order = Cluster_getOrder(cls);
    printf("      %-18s:   \n\n","clusinfo");

    // 0-th order
    int iter_0th = Cluster_getClusinfo_itercluster(cls,0,0)[0];
    int ncluster_0th = Cluster_getClusinfo_ncluster(cls,0);
    printf("%10s 0-th order ... \n"," ");
    printf("%10s Cluster[0][%5d] : %3d (Iter)\n"," ",0,iter_0th);
    printf("\n%10s Total 0-Cluster # :  %d # \n\n\n"," ",ncluster_0th);

    for (int i = 1; i <= order; i++){
        printf("\n%10s %d-th order ... \n"," ",i);
        int ncluster = Cluster_getClusinfo_ncluster(cls, i);

        for (int j = 1; j<= ncluster; j++){

            if (verbosity || (j<5 || j>ncluster-4)){
                
                int* itercluster = Cluster_getClusinfo_itercluster(cls,i,j);
                int iter = itercluster[0];

                printf("%10s Cluster[%d][%5d] : %3d ["," ",i,j,iter);
                for (int k = 1; k <= i; k++){
                    printf("%-3d",itercluster[k]);
                    if (k != i){ printf(", ");} else { printf(" ]\n");}
                }
            }
            
            if (!verbosity && j==5){
                printf("%10s      : \n"," ");
            }
        }
        printf("\n%10s Total %d-Cluster # :  %d # \n\n\n"," ",i,ncluster);
    }
}

/* Low Level --------------------------------------------------------*/

// alloc
void Cluster_allocNk(Cluster* cls){ // 0:0-th ... k:k-th order 
    int length = cls->order + 1;
    cls->nk = allocInt1d(length);
}

void Cluster_allocClusinfo(Cluster* CCE, int order){
    CCE->clusinfo = (int***)allocArray1d(order+1,sizeof(int**));

    for (int i=0; i<=order; i++){
        CCE->clusinfo[i] = CCE->clusinfo[i] = (int**)allocArray1d(1,sizeof(int*));
        CCE->clusinfo[i][0] = (int*)allocArray1d(1,sizeof(int));
        CCE->clusinfo[i][0][0] = 1;
    }

}

int Cluster_setClusinfo_addcluster(Cluster* CCE, int order, int iter, int* cluster){

    if (order==0){
        printf("Error: Cluster_setClusinfo_addcluster: at 0-th order, you cannot add cluster\n");
        exit(1);
    }

    if (CCE->clusinfo==NULL){
        printf("Error: Cluster_setClusinfo_addcluster: clusinfo is not allocated\n");
        exit(1);
    }

    if (iter != 1 ){
        printf("Warning: Cluster_setClusinfo_addcluster: iter is not 1 when you add new cluster ... \n");
    }

    int ncluster_old = CCE->clusinfo[order][0][0];
    int ncluster_new = ncluster_old + 1;
    int ic = ncluster_new - 1; // new cluster index

    // set the number of cluster at order
    CCE->clusinfo[order][0][0] = ncluster_new;

    // realloc
    reallocInt2d(&(CCE->clusinfo[order]),ncluster_old,ncluster_new, order+1);

    // set new cluster
    copyInt1d(&(CCE->clusinfo[order][ic][1]),cluster,order+1);

    // set iter at 0-th element (It is the number how many it will calculated)
    CCE->clusinfo[order][ic][0] = iter;

    return ic;
}

void Cluster_setClusinfo_chgcluster(Cluster* CCE, int order, int* cluster, int ic){
    if (order==0){
        printf("Error : Cluster_setClusinfo_chgcluster : at 0-th order, you cannot change cluster information\n");
        printf("For the 0th order, you can only change iter\n");
        exit(1);
    }
    copyInt1d(&(CCE->clusinfo[order][ic][1]),cluster,order);
}

void Cluster_setClusinfo_chgiter(Cluster* CCE, int order, int iter, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_setClusinfo_chgiter : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    CCE->clusinfo[order][ic][0] = iter;
}


int Cluster_getClusinfo_ncluster(Cluster* CCE, int order){
    if (order==0){
        return 1;
    }
    return CCE->clusinfo[order][0][0]-1;
}

int* Cluster_getClusinfo_itercluster(Cluster* CCE, int order, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_getClusinfo_cluster : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    return CCE->clusinfo[order][ic];
}

int Cluster_getClusinfo_iter(Cluster* CCE, int order, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_getClusinfo_cluster : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    return CCE->clusinfo[order][ic][0];
}

int* Cluster_getClusinfo_cluster_copy(Cluster* CCE, int order, int ic){
    if (order==0){
        printf("Error : Cluster_getClusinfo_cluster_copy : at 0-th order, you cannot get cluster information\n");
        exit(1);
    }
    int* cluster = allocInt1d(order);
    copyInt1d(cluster,&(CCE->clusinfo[order][ic][1]),order);
    return cluster;
}

int*** Cluster_getClusinfo(Cluster* CCE){
    return CCE->clusinfo;
}

// get 

char* Cluster_getMethod(Cluster* CCE){
    return CCE->method;
}

int Cluster_getOrder(Cluster* cls){
    return cls->order;
}

bool Cluster_getAddsubclus(Cluster* cls){
    return cls->addsubclus;
}

int* Cluster_getNk(Cluster* cls){
    return cls->nk;
}

int Cluster_getNk_order(Cluster* cls, int i){
    if (i<cls->order && i>=0){
        return cls->nk[i];
    }else{
        printf("Error: Cluster_getNk_order: i is out of range\n");
        exit(1);
    }
}



// set 
void Cluster_setOrder(Cluster* cls, int order){
    cls->order = order;
}

void Cluster_setMethod(Cluster* cls, char* method){
    strcpy(cls->method, method);
}
void Cluster_setAddsubclus(Cluster* cls, bool addsubclus){
    cls->addsubclus = addsubclus;
}
void Cluster_setNk(Cluster* cls, int* nk){
    copyInt1d(cls->nk,nk,cls->order+1);
}

// free
void Cluster_freeNk(Cluster* cls){
    freeInt1d(cls->nk);
}

void Cluster_freeClusinfo(Cluster* cls){
    
    int order = cls->order;

    if (cls->clusinfo!=NULL){
        for (int i=0; i<=order; i++){

            int ncluster = Cluster_getClusinfo_ncluster(cls,i);

            if (i==0){
                freeArray2d((void**)cls->clusinfo[i],1);
            }
            else{
                freeArray2d((void**)cls->clusinfo[i],ncluster);
            }
        }
        free(cls->clusinfo);
        cls->clusinfo = NULL;
    }
}