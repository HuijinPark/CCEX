#include "../include/cluster.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include "../include/cluster_hash.h"
#include "../include/cluster_pcce.h"
#include "../include/bath.h"

/* High Level --------------------------------------------------------*/

// init
Cluster* Cluster_init(){
    Cluster* cls = (Cluster*)allocArray1d(1,sizeof(Cluster));
    Cluster_setOrder(cls, 0);
    Cluster_setAddsubclus(cls, false);
    cls->nk = NULL;
    cls->method[0] = '\0';
    cls->clusinfo = NULL;
    
    return cls;
}

void Cluster_clusterize_pcce(Cluster* cls, BathArray* ba, QubitArray* qa, Config* config, int rank){
//void Cluster_clusterize_pcce(Cluster* cls_pcce, Cluster* cls, BathArray* ba, QubitArray* qa, Config* config, int rank){

    if (strcasecmp(cls->method, "pcce")==0){
        Partition_info pinfo;
        // Cut the # of spins (pcce grouping) //
        set_spinfinite(ba, qa, cls->sK, &pinfo, rank);

        int ncenter              = pinfo.ncenter;
        int pcce_nspin           = pinfo.pcce_nspin;
        Point* best_centers      = (Point*) malloc(sizeof(Point) * ncenter);
        int*   best_assigned_idx = (int*)   malloc(sizeof(int)   * pcce_nspin);
        
        // Find the Center position & assigned index // 
        // bath >> center 
        simulator_cluster_partition(ba->bath, &pinfo, cls, best_centers, best_assigned_idx);

        int order = Cluster_getOrder(cls);
        //int nspin = BathArray_getNspin(ba);

        // cls->info alloc => cls_pcce 
        Cluster_allocClusinfo(cls, order); // clusinfo[order+1] : 0 ~ order, clusinfo[i][0][0] = 1

        if (order==0){
            if (rank==0){printf("\n\t Clustering 0th order ... \n");}
        }
        else if (order == 1){ // For pCCE, + k=1 (order = n)
            
            if (rank==0){printf("\n\t Clustering 1st order ... \n");}

            int cluster[1] = {0};
            int iter = 1;
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
        //if (rank == 0){
        //    for (int tempi=0; tempi<ncenter; tempi++){
        //        for (int tempj=0; tempj<cls->sK; tempj++){
        //            printf("cs2dArr[%d][%d]=%d\n", tempi, tempj, cs2dArr[tempi][tempj]);
        //        }
        //    }
        //}
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
}

void Cluster_clusterize(Cluster* cls, BathArray* ba, QubitArray* qa, Config* config){

    int order = Cluster_getOrder(cls);
    int nspin = BathArray_getNspin(ba);

    Cluster_allocClusinfo(cls, order); // clusinfo[order+1] : 0 ~ order, clusinfo[i][0][0] = 1

    if (order==0){
        if (rank==0){printf("\n\t Clustering 0th order ... \n");}
    }
    else if (order == 1){ // For pCCE, + k=1 (order = n)
        
        if (rank==0){printf("\n\t Clustering 1st order ... \n");}

        int cluster[1] = {0};
        int iter = 1;
        for (int spinidx = 0; spinidx < nspin; spinidx++){
            cluster[0] = spinidx;
            Cluster_setClusinfo_addcluster(cls, order, iter, cluster);
        }

    }
    else if (order > 1){ 

        int** cmap    = NULL; // connectivity map : cmap[nspin][nspin] = 1 if connected else 0
        float** stmap = NULL; // strength map : stmap[nspin][nspin] = strength
        int** spmap   = NULL; // sparse map : spmap[nspin][ncol] = n connected spin + 1d array of connected spins index

        //!!if (strcasecmp(cls->method, "pcce")==0){
        //!!    double** centerpositions = NULL; // start from 0
        //!!    int centerarraylength = 0;
        //!!    int* spinidx = NULL ; // start from 0
        //!!    getCenters(&centerpositions, &centerarraylength, &spinidx, ba);
        //!!    int BathArray_getNspin(ba)
        //!!}

        float rdip    = Config_getRdip(config);
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

        freeInt2d(&cmap,nspin);
        freeInt2d(&spmap,nspin);
        freeFloat2d(&stmap,nspin); 

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

    if (strcasecmp(method, "gcce")!=0){
        Cluster_setClusinfo_chgiter(cls, 0, 0, 0);
    }
}


void change_clusinfo_for_pcce(int*** clusinfo_pcce, int*** clusinfo_old){
    
    reportClusinfo(clusinfo_old, 2);
}

// free
void Cluster_freeAll(Cluster* cls){
    freeArray1d((void**)&cls);
}

/* Low Level --------------------------------------------------------*/

// report
void Cluster_report(Cluster* cls){

    printTitle("Structure Cluster");

    printSubTitle("General properties");
    printStructElementInt("order",Cluster_getOrder(cls));
    printStructElementChar("method",cls->method);
    printStructElementBool("addsubclus",Cluster_getAddsubclus(cls));
    
    if (strcasecmp(Cluster_getMethod(cls),"pcce") == 0){
        printStructElementInt("sK",Cluster_getSk(cls));
        printStructElementInt("max_iter",Cluster_getMax_iter(cls));
        printStructElementInt("max_trial",Cluster_getMax_trial(cls));
        printStructElementBool("kmeans_pp",Cluster_getKmeans_pp(cls));
        printStructElementBool("iter_detail",Cluster_getIter_detail(cls));
    }

    if (cls->nk!=NULL){
        printSubTitle("The number of clusters for each order");
        Cluster_reportNk(cls);
    }
    
    if (cls->clusinfo!=NULL){
        printSubTitle("Clustering information");
        Cluster_reportClusinfo(cls);
    }

    printf("\n");
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

void reportClusinfo(int*** clusinfo, int order){
 
    // 0-th order
    int iter_0th = clusinfo[0][0][0];
    int ncluster_0th = 1;
    printf("%10s 0-th order ... \n"," ");
    printf("%10s Cluster[0][%5d] : %3d (Iter)\n"," ",0,iter_0th);
    printf("\n%10s Total 0-Cluster # :  %d # \n\n\n"," ",ncluster_0th);

    for (int i = 1; i <= order; i++){
        printf("\n%10s %d-th order ... \n"," ",i);
        int ncluster = clusinfo[i][0][0]-1;

        for (int j = 1; j<= ncluster; j++){

            if (verbosity || (j<5 || j>ncluster-4)){
                
                int* itercluster = &(clusinfo[i][j][0]);
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

void Cluster_allocClusinfo(Cluster* cls, int order){
    cls->clusinfo = (int***)allocArray1d(order+1,sizeof(int**));

    for (int i=0; i<=order; i++){
        cls->clusinfo[i] = cls->clusinfo[i] = (int**)allocArray1d(1,sizeof(int*));
        cls->clusinfo[i][0] = (int*)allocArray1d(1,sizeof(int));
        cls->clusinfo[i][0][0] = 1;
    }

}

int Cluster_setClusinfo_addcluster(Cluster* cls, int order, int iter, int* cluster){

    if (order==0){
        printf("Error: Cluster_setClusinfo_addcluster: at 0-th order, you cannot add cluster\n");
        exit(1);
    }

    if (cls->clusinfo==NULL){
        printf("Error: Cluster_setClusinfo_addcluster: clusinfo is not allocated\n");
        exit(1);
    }

    if (iter != 1 ){
        printf("Warning: Cluster_setClusinfo_addcluster: iter is not 1 when you add new cluster ... \n");
    }

    int ncluster_old = cls->clusinfo[order][0][0];
    int ncluster_new = ncluster_old + 1;
    int ic = ncluster_new - 1; // new cluster index

    // set the number of cluster at order
    cls->clusinfo[order][0][0] = ncluster_new;

    // realloc
    reallocInt2d(&(cls->clusinfo[order]),ncluster_old,ncluster_new, order+1);

    // set new cluster
    copyInt1d(&(cls->clusinfo[order][ic][1]),cluster,order+1);

    // set iter at 0-th element (It is the number how many it will calculated)
    cls->clusinfo[order][ic][0] = iter;

    return ic;
}

void Cluster_setClusinfo_chgcluster(Cluster* cls, int order, int* cluster, int ic){
    if (order==0){
        printf("Error : Cluster_setClusinfo_chgcluster : at 0-th order, you cannot change cluster information\n");
        printf("For the 0th order, you can only change iter\n");
        exit(1);
    }
    copyInt1d(&(cls->clusinfo[order][ic][1]),cluster,order);
}

void Cluster_setClusinfo_chgiter(Cluster* cls, int order, int iter, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_setClusinfo_chgiter : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    cls->clusinfo[order][ic][0] = iter;
}


int Cluster_getClusinfo_ncluster(Cluster* cls, int order){
    if (order==0){
        return 1;
    }
    return cls->clusinfo[order][0][0]-1;
}

int* Cluster_getClusinfo_itercluster(Cluster* cls, int order, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_getClusinfo_cluster : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    return cls->clusinfo[order][ic];
}

int Cluster_getClusinfo_iter(Cluster* cls, int order, int ic){
    if (order!=0 && ic==0){
        printf("Error : Cluster_getClusinfo_cluster : order!=0, ic=0 is not the cluster index\n");
        exit(1);
    }
    return cls->clusinfo[order][ic][0];
}

int* Cluster_getClusinfo_cluster_copy(Cluster* cls, int order, int ic){
    if (order==0){
        printf("Error : Cluster_getClusinfo_cluster_copy : at 0-th order, you cannot get cluster information\n");
        exit(1);
    }
    int* cluster = allocInt1d(order);
    copyInt1d(cluster,&(cls->clusinfo[order][ic][1]),order);
    return cluster;
}

int*** Cluster_getClusinfo(Cluster* cls){
    return cls->clusinfo;
}

// get 

char* Cluster_getMethod(Cluster* cls){
    return cls->method;
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

int Cluster_getSk(Cluster* cls){
    return cls->sK;
}

int Cluster_getMax_iter(Cluster* cls){
    return cls->max_iter;
}

int Cluster_getMax_trial(Cluster* cls){
    return cls->max_trial;
}

bool Cluster_getKmeans_pp(Cluster* cls){
    return cls->kmeans_pp;
}

bool Cluster_getIter_detail(Cluster* cls){
    return cls->iter_detail;
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

void Cluster_setSk(Cluster* cls, int sK){
    cls->sK = sK;
}

void Cluster_setMax_trial(Cluster* cls, int max_trial){
    cls->max_trial = max_trial;
}

void Cluster_setMax_iter(Cluster* cls, int max_iter){
    cls->max_iter = max_iter;
}

void Cluster_setKmeans_pp(Cluster* cls, bool kmeans_pp){
    cls->kmeans_pp= kmeans_pp;
}

void Cluster_setIter_detail(Cluster* cls, bool iter_detail){
    cls->iter_detail= iter_detail;
}




// free
void Cluster_freeNk(Cluster* cls){
    freeInt1d(&(cls->nk));
}

void Cluster_freeClusinfo(Cluster* cls){
    
    int order = cls->order;

    if (cls->clusinfo!=NULL){
        for (int i=0; i<=order; i++){

            int ncluster = Cluster_getClusinfo_ncluster(cls,i);

            if (i==0){
                freeArray2d((void***)&(cls->clusinfo[i]),1);
            }
            else{
                freeArray2d((void***)&(cls->clusinfo[i]),ncluster);
            }
        }
        free(cls->clusinfo);
        cls->clusinfo = NULL;
    }
}


