#include <string.h>  /* strcpy */
#include <stdlib.h>  /* malloc */
#include <stdio.h>   /* printf */
#include <time.h>
#include <math.h>
#include <string.h>
#include "../zlib/uthash/src/uthash.h"
#include "../include/cluster_hash.h"

extern int rank;

void clusterizeHash(Cluster* CCE, int nspin, int** spmap, float** stmap){

    int order = Cluster_getOrder(CCE);
    int* nks = Cluster_getNk(CCE);
    bool addsubclus = Cluster_getAddsubclus(CCE);

    //////////////////////////////////////////////////
    // Clusterize
    //////////////////////////////////////////////////
    if (rank==0){printf("\n\t CCE1 ... \n");}
    HashCluster *hashClusters = NULL;
    makeHashClusterO1(&hashClusters,nspin);
    updateNk(&(nks),order,hashClusters);
    
    if (rank==0){printf("\n\t CCE2 ... \n");}
    makeHashClusterO2(&hashClusters,nspin,spmap,stmap);
    updateNk(&(nks),order,hashClusters);
    
    // Higher order clusters
    makeHashClusterOn(&hashClusters,order,nspin,spmap,stmap,nks);
    updateNk(&(nks),order,hashClusters);
    if (rank==0 && verbosity){
        printClusters(hashClusters);
    }
    //////////////////////////////////////////////////
    // Add all sub clusters
    //////////////////////////////////////////////////
    if (addsubclus) {
        for (int n=order; n>2; n--){
            if (rank==0){printf("\n\t Add subclusters ... \n");}
            addSubClusters(&hashClusters,nspin,stmap,n);
        }
        updateNk(&(nks),order,hashClusters);
    }
    //////////////////////////////////////////////////
    // Solve tilde
    // Here, we'll find subclusters, 
    // and the number how many you will do the calculation.
    // Additionally, we will move the hashcluster information
    // to the CCE->clusinfo
    solveTilde(&hashClusters,CCE,nspin);
    //////////////////////////////////////////////////
    freeHashCluster(&hashClusters,order);
}

void freeHashCluster(HashCluster** hashClusters, int order){

    for (int n=1; n<order+1; n++){
        HashCluster *existingClusters = findCluster(*hashClusters,n);
        Property *existingCluster, *tmp;
        if (existingClusters != NULL){
            HASH_ITER(hh,existingClusters->prop, existingCluster, tmp)
            {
                free(existingCluster->spins);
                free((char*)existingCluster->id);
                HASH_DEL(existingClusters->prop,existingCluster);
            }
        }
        free(existingClusters->prop);
        free((char*)existingClusters->N);
        HASH_DEL(*hashClusters,existingClusters);
    }
}

void convertClusinfoToHash(HashCluster** hashClusters, Cluster* CCE){

    int order = Cluster_getOrder(CCE);
    int nspin = Cluster_getClusinfo_ncluster(CCE,1);

    for (int n=1; n<order+1; n++){

        int ncluster = Cluster_getClusinfo_ncluster(CCE,n);
        int length = ncluster+1;

        for (int i=1; i<length; i++){
            
            int iter = Cluster_getClusinfo_iter(CCE,n,i);
            int* cluster = Cluster_getClusinfo_cluster_copy(CCE,n,i);
            float strength = 0.0; // Not used (because we don't do sort)
            
            if (iter != 1){
                printf("Error,convertClusinfoToHash, iter is not 1 (It already has been processed the solveTilde part)\n");
            }

            // make id
            int maxlenId = setMaxLengthStr(nspin,n);
            char* id = (char*)calloc(maxlenId,sizeof(char));
            typeStr(&id,cluster,n,nspin);

            // add cluster
            int res = addCluster(hashClusters,n,id,cluster,strength,iter);
            if (!res) {
                printf("error,convertClusinfoToHash, addcluster"); 
                exit(1);
            }
        }
    }
}


int addCluster(HashCluster** hashClusters, int order, const char* id,int* spins, float strength,int count)
{
    char* orderStr = (char*)calloc(countDigits(order)+1,sizeof(char));
    sprintf(orderStr,"%d",order);
    /* Check overlapping of key */
    HashCluster *c;
    HASH_FIND_STR(*hashClusters, orderStr, c);
    
    Property *p;
    if (c != NULL){
        HASH_FIND_STR(c->prop, id, p);
    }

    /* Cluster is not found*/
    if (c == NULL || p == NULL) {
        /* Build new Hashcluster */
        if (c==NULL){
            c = (HashCluster*)malloc(sizeof(HashCluster));
            if (c==NULL){
                exit(-1);
            }
            /* Value assignment 1st */
            c->N = orderStr;
            c->prop = NULL;
            HASH_ADD_KEYPTR(hh, *hashClusters, c->N, strlen(c->N), c);
        }
        p = (Property*)malloc(sizeof(Property));
        ///* Value assignment 2nd */
        p->id = id; 
        p->spins = spins;
        p->strength = strength;
        p->count = count;
        /* Add new Hashcluster*/
        HASH_ADD_KEYPTR(hh, c->prop, p->id, strlen(p->id), p);
    }
    else{ 
        free(orderStr);
        // The cluster already exist, so we cannot add
        return 0;
    }
    // The cluster is successfully added. 
    return 1;
}
int setMaxLengthStr(int nSpin, int order){
    int nSpinDigit = countDigits(nSpin);
    /*Additional length 1 for null*/
    return (nSpinDigit*order) + (order-1) + (1); 
}

void makeHashClusterO1(HashCluster** hashclusters, int nspin){

    int n = 1;
    int res = 0;
    for (int i=0; i<nspin; i++)
    {
        /////////////////////////////////////
        int* single = (int*)calloc(1,sizeof(int));
        single[0] = i;
        /////////////////////////////////////
        float strength = 0.0;
        /////////////////////////////////////
        int maxlenId = setMaxLengthStr(nspin,n);
        char* id = (char*)calloc(maxlenId,sizeof(char));
        typeStr(&id,single,n,nspin);
        /////////////////////////////////////
        int count = 1;
        /////////////////////////////////////
        res = addCluster(hashclusters,n,id,single,strength,count);
        if (!res) {
            printf("error,makeHashClusterO1, addcluster"); 
            exit(1);
        }
    }
}


void makeHashClusterO2(HashCluster** hashclusters, int nspin, int** spmap, float** stmap){

    int n = 2;
    int res=0;
    for (int row=0; row<nspin; row++)
    {
        int nConnectedSpins = spmap[row][0];
    
        for (int icol=1; icol<nConnectedSpins; icol++)
        {
            int col = spmap[row][icol];  
            int* pair = (int*)calloc(2,sizeof(int));
            pair[0] = row;
            pair[1] = col;
            if (row<col){
                QuickSort(&pair,0,1);
                float strength = stmap[row][col];
                /////////////////////////////////////
                int maxlenId = setMaxLengthStr(nspin,n);
                char* id = (char*)calloc(maxlenId,sizeof(char));
                typeStr(&id,pair,n,nspin);
                /////////////////////////////////////
                int count = 1;
                /////////////////////////////////////
                res = addCluster(hashclusters,n,id,pair,strength,count);
                if (!res) {
                    printf("error,makeHashClusterO2, addcluster"); 
                    exit(1);
                }
            }
        }
    }
    HASH_SORT(findCluster(*hashclusters,2)->prop,by_strength);
}


void makeHashClusterOn(HashCluster** hashclusters, int order, int nspin, int** spmap, float** stmap, int* nks){

    int res = 0;

    for (int n=3; n<order+1; n++){
        if (rank==0){printf("\n\t CCE%d ... \n",n);}
        int o = n-1;
        HashCluster *existingClusters = findCluster(*hashclusters,o);
//        printf("order : %s\n", existingClusters ? existingClusters->N : "unknown\n");
        Property *existingCluster, *tmp;
        if (existingClusters != NULL){
            HASH_ITER(hh,existingClusters->prop, existingCluster, tmp)
            {
               // take one cluster
//                printf("id              : %s\n",existingCluster->id);
//                printf("strength        : %f\n",existingCluster->strength);
//                printf("current cluster : %s\n",existingCluster->id);
                int* cluster = existingCluster->spins;
                float strength = existingCluster->strength;

                int nPossibleSpins = 1;
                int* possibleSpins = (int*)calloc(nPossibleSpins,sizeof(int));
                possibleSpins[0] = nPossibleSpins;

                int nNewClusters = 1;
                int** newClusters = (int**)calloc(nNewClusters,sizeof(int*));
                newClusters[0] = (int*)calloc(n,sizeof(int));
                newClusters[0][0] = nNewClusters;

                // take one spin in a cluster to find connected spins
                //clock_t startTimeFindClusters, endTimeFindClusters;
                //startTimeFindClusters = clock();

                for (int i=0; i<o; i++){
                    
                    int spin = cluster[i];
                    int nConnectedSpins = spmap[spin][0];

                    for (int p=1; p<nConnectedSpins; p++){
                        
                        int possibleSpin = spmap[spin][p];

                        // overlap check to existing cluster  
                        int IsExistInCurrentCluster = binarySearch(cluster,0,o-1,possibleSpin);
                        // overlap check to PossibleSpin list
                        int IsExistInPossibleSpins = binarySearch(possibleSpins,1,nPossibleSpins-1,possibleSpin);

//                        printf("possible spin : %2d , (judge : %2d,%2d => ",possibleSpin,IsExistInCurrentCluster,IsExistInPossibleSpins);
                        if (IsExistInCurrentCluster == -1 && IsExistInPossibleSpins == -1){
//                            printf(" possible )\n");
                            // Get possible spin list
                            possibleSpins = (int*)realloc(possibleSpins,sizeof(int)*(nPossibleSpins+1));
                            possibleSpins[nPossibleSpins] = possibleSpin;
                            // Get potential new cluster
                            int* newClusterTmp = (int*)calloc(n,sizeof(int));

                            ////clock_t startTimeTmp;
                            ////startTimeTmp = clock();

                            addSpin(&newClusterTmp, cluster, o, possibleSpin);
                            QuickSort(&(newClusterTmp),0,n-1);
                            // Check the tempoaray cluster exist in the created cluster list
                            int IsExistInNewClusters = search2dArr(newClusters,1,nNewClusters-1,0,n-1,newClusterTmp);
                            if (IsExistInCurrentCluster == -1 ){
                                // Set NewCluster
                                newClusters = (int**)realloc(newClusters,sizeof(int*)*(nNewClusters+1));
                                newClusters[nNewClusters] = newClusterTmp;
                                // strength
                                //float newStrength = minStrength(strength, stmap, cluster, o, possibleSpin);
                                float newStrength = addAllStrength(strength, stmap, cluster, o, possibleSpin);
                                int maxlenId = setMaxLengthStr(nspin,n);
                                char* id = (char*)calloc(maxlenId,sizeof(char));
                                typeStr(&id,newClusters[nNewClusters],n,nspin);
                                // count
                                int count = 1;
                                // add
                                res = addCluster(hashclusters,n,id,newClusters[nNewClusters],newStrength,count);
                                nNewClusters++;
                            }
                            else{ free(newClusterTmp); }
                            nPossibleSpins++;
                        }
                        else{
                        }
                    }
                }
                possibleSpins[0] = nPossibleSpins;
                newClusters[0][0] = nNewClusters;

                free(possibleSpins);
            }
            HASH_SORT(findCluster(*hashclusters,n)->prop,by_strength);
            //HASH_DELETE(hh2,potentialClusters,strongCluster);
            //printProperties(findCluster(*hashclusters,n)->prop); 
            /////////////////////////////////////////////////////////
            //
            //  Cut the number of clusters up to Nk
            //
            /////////////////////////////////////////////////////////
            Property *madeClusters, *madeCluster, *tmp;
            madeClusters = findCluster(*hashclusters,n)->prop;
            /////////////////////////////////////////////////////////
            int nk = nks[n];
            int nCluster = 0;
            nCluster = HASH_COUNT(madeClusters);
            if (nk == 0){
                nk = nCluster;
            }
            /////////////////////////////////////////////////////////
            if (nCluster > nk ){
                if (rank==0){printf("\t\t apply Nk=%d ... \n",nk);}
                nCluster = 0;
                HASH_ITER(hh, madeClusters, madeCluster, tmp)
                {
                    nCluster++;
                    if (nCluster > nk){
                        HASH_DEL(madeClusters,madeCluster);
                    }
                }
            }
            else if (nCluster == nk){;}
            else{
                if (n == order){
                    char message[500];
                    sprintf(message,"Warning, during making clusters of %d order, made clusters(%d) is smaller than given nCluster(%d)",n,nCluster,nk);
                    printMessage(message);
                }
                else{;}
            }
            /////////////////////////////////////////////////////////
            //printProperties(findCluster(*hashclusters,n)->prop); 
        }
    }
}

void typeStr(char** destination, int* cluster, int order, int nSpin){
    int nSpinDigit = countDigits(nSpin);
    for (int i=0; i<order; i++){
        char str[nSpinDigit];
        sprintf(str,"%d", cluster[i]);
        strcat(*destination,str);
        if (i!=order-1){
            strcat(*destination,",");
        }
    }
}

int by_strength(const Property* a, const Property* b)
{   
    // descending order for stregnth (ex. 100 - 50 - 0)
    if (a->strength == b->strength) return 0;
    return a->strength < b->strength ? 1 : -1;
}

HashCluster* findCluster(HashCluster* hashClusters, int order){

    char* orderStr = (char*)calloc(countDigits(order)+1,sizeof(char));
    sprintf(orderStr,"%d",order);
    /* Check overlapping of key */
    HashCluster *c;
    HASH_FIND_STR(hashClusters, orderStr, c);
    
    free(orderStr);
    return c;
}

Property* findProperty(HashCluster* hashClusters, int order, const char* id){

    HashCluster *c = findCluster(hashClusters,order);
    
    Property *p;
    if (c != NULL){
        HASH_FIND_STR(c->prop, id, p);
    }

    return p;
}

int binarySearch(int* arr, int low, int high, int target) {

    while (low <= high) {
        int mid = low + (high - low) / 2;

        if (arr[mid] == target)
            return mid;

        if (arr[mid] < target)
            low = mid + 1;
        else
            high = mid - 1;
    }

    return -1;  
} 

void addSpin(int** newcluster, int* oldcluster, int oldn, int spin){

    for (int i=0; i<oldn; i++){ 
        (*newcluster)[i] = oldcluster[i];
    }
    (*newcluster)[oldn] = spin;
}

float minStrength(float currentStrength, float** strengthMap, int* oldcluster, int oldn, int newspin){

    float newStrength = currentStrength;

    for (int i=0; i<oldn; i++){ 
        int spin = oldcluster[i];
        if (newStrength > strengthMap[spin][newspin]){
            newStrength = strengthMap[spin][newspin];
        }
    }

    return newStrength;
}

void printProperties(Property* hashProperties)
{
    Property *item2, *tmp2;
    int nCluster = 0;
    int count = 0;
    HASH_ITER(hh, hashProperties, item2, tmp2) {
        printf("\t\t$items[%15s] : ", item2->id);

        int* clusterFromId = parseClusterIdToIntArray(item2->id,&count);
        for (int i =0; i<count; i++){
            printf("%4d ",item2->spins[i]);
        }
        printf(" (strength : %12.5f)\n",item2->strength);
        nCluster++;
        free(clusterFromId);
    }
    printf("\n");
    printf("\t\tTotal %d-Cluster # :  %d\n",count,nCluster);
    printf("\t\t_________________________________________________________\n");
    printf("\n\n");
    
}

void printClusters(HashCluster* hashClusters)
{
    HashCluster *item1, *tmp1;
    Property *item2, *tmp2;
    HASH_ITER(hh, hashClusters, item1, tmp1) {
        int nCluster = 0;
        HASH_ITER(hh, item1->prop, item2, tmp2) {
            printf("$items[%3s][%10s] : ", item1->N, item2->id);

            int count = 0;
            int* clusterFromId = parseClusterIdToIntArray(item2->id,&count);
            for (int i =0; i<count; i++){
                printf("%3d ",item2->spins[i]);
            }
            printf(" (strength : %10.5f)\n",item2->strength);
            nCluster++;
            free(clusterFromId);
        }
        printf("\n");
        printf("Total %s-Cluster # :  %d\n",item1->N,nCluster);
        printf("_________________________________________________________\n");
        printf("\n\n");
    }
}

// Parse a string into an array of integers
int* parseClusterIdToIntArray(const char* id, int* count){

    // Destination
    int* intArray;

    // Copy the input string
    char* strCopy = strdup(id);
    
    // Convert each token to an integer and store it in the array
    char* token = strtok(strCopy, ",");
    int index = 0;
    while (token != NULL) {
        // Allocation
        if (index == 0){
            intArray = (int*)calloc(index+1,sizeof(int));    
        }
        else{
            intArray = (int*)realloc(intArray,sizeof(int)*(index+1));
        }
        // Give value at index "and then" add 1 to index
        intArray[index] = atoi(token);
        token = strtok(NULL, ",");
        index++;
    }

    // Give the length of int array
    *count = index;

    // Set the count of the resulting array and free the memory
    free(strCopy);
    return intArray;
}

 
int countDigits(int number) {
    int count = 0;

    // Handle the special case when the number is 0
    if (number == 0) {
        return 1;
    }

    // For negative numbers, convert to positive
    if (number < 0) {
        number = -number;
    }

    while (number != 0) {
        // Increment the count for each digit
        count++;

        // Remove the last digit by dividing by 10
        number = number / 10;
    }

    return count;
}

int search2dArr(int** arr, int low, int high, int left, int  right, int* target) {
   
    for(int i=low; i<=high; i++){
        int count = 0;
        for(int j=left; j<=right; j++){
            if (arr[i][j] == target[j]){
                count++;
            }
        }
        if (count == right-left+1){return i;} 
    } 
    return -1;  
}

void updateNk(int** Nk, int order, HashCluster* hashClusters){

    unsigned int nk;

    HashCluster *item, *tmp;
    HASH_ITER(hh, hashClusters, item, tmp){
        int k = atoi(item->N);
        nk = HASH_COUNT(item->prop);
        (*Nk)[k] = nk;
    }

    if (rank==0){
        printf("\t\t update \"Nk\"\n\t\t ");
        for (int i=0; i<=order; i++){
            printf("%d:%d  ",i,(*Nk)[i]);
        }
        printf("\n\n");
    }

}

float addAllStrength(float currentStrength, float** strengthMap, int* oldcluster, int oldn, int newspin){

    float addedStrength = currentStrength;
    for (int i=0; i<oldn; i++){ 
        int spin = oldcluster[i];
        addedStrength+=strengthMap[spin][newspin];
    }
    return addedStrength;
}

float addAllStrengthForAllSpins(float** strengthMap, int* cluster, int order){

    float addedStrength = 0.0; 
    for (int i=0; i<order; i++){ 
        for (int j=i+1; j<order; j++){ 
            addedStrength+=strengthMap[cluster[i]][cluster[j]];
        }
    }

    return addedStrength;
}

// Function to generate combinations and store them in data
void generateCombinations(int* arr, int*** data, int* tempCombination, int start, int end, int index, int r, int* nCombination) {
    // If combination size is r, copy it to data
    if (index == r) {
        //printf("nCombination = %d\n",*nCombination);
        for (int i = 0; i < r; i++) {
            (*data)[*nCombination][i] = tempCombination[i];
            //printf("%d ",(*data)[*nCombination][i]);
        }
        //printf("\n");
        *nCombination += 1; 
        return;
    }

    // If we reach the end of the array or the desired combination size, return
    if (start == end) {
        return;
    }

    // Include the current element in the combination
    tempCombination[index] = arr[start];
    generateCombinations(arr, data, tempCombination, start + 1, end, index + 1, r, nCombination);

    // Exclude the current element from the combination
    generateCombinations(arr, data, tempCombination, start + 1, end, index, r, nCombination);
}

void addSubClusters(HashCluster** hashclusters, int nspin, float** stmap, int n){

    int res = 0;

    HashCluster *existingClusters = findCluster(*hashclusters,n);
    Property *existingCluster, *tmp;
    if (existingClusters != NULL){
        HASH_ITER(hh,existingClusters->prop, existingCluster, tmp)
        {
            int* cluster = existingCluster->spins;
            //printf("cluster :  %s\n",existingCluster->id);
            for (int r = n-1; r >=2; r--){
                ///////////////////////////////////////////
                // Find subclusters
                ///////////////////////////////////////////
                // Calculate the number of combinations (n choose r)
                int numCombinations = 1;
                for (int i = 0; i < r; i++) {
                    numCombinations *= (n - i);
                    numCombinations /= (i + 1);
                }
                // Allocate memory for data
                int** subclusters = (int**)calloc(numCombinations,sizeof(int*));
                for (int i = 0; i < numCombinations; i++) {
                    subclusters[i] = (int*)calloc(r,sizeof(int));
                }
                int* tempCombination = (int*)calloc(r,sizeof(int)); // Temporary storage for combinations
                int nCombination = 0 ;
                generateCombinations(cluster, &subclusters, tempCombination, 0, n, 0, r, &nCombination);
                // Free memory for tempCombination
                free(tempCombination); 
                ///////////////////////////////////////////
                // Add subclusters
                ///////////////////////////////////////////
                for (nCombination=0; nCombination<numCombinations; nCombination++){
                    QuickSort(&(subclusters[nCombination]),0,r-1);
                    //printf("subcluster[%2d] : ",nCombination);
                    //for ( int k =0; k<r; k++){
                    //    printf("  %d",subclusters[nCombination][k]);
                    //}
                    //printf("\n");
                    // strength
                    float newStrength = addAllStrengthForAllSpins(stmap, subclusters[nCombination], r);
                    int maxlenId = setMaxLengthStr(nspin,r);
                    char* id = (char*)calloc(maxlenId,sizeof(char));
                    typeStr(&id,subclusters[nCombination],r,nspin);
                    // count
                    int count = 1;
                    // add
                    res = addCluster(hashclusters,r,id,subclusters[nCombination],newStrength,count);
                    if (!res){
                        free(id);
                        free(subclusters[nCombination]);
                    }
                }
            }
        }
    }
}

void solveTilde(HashCluster** hashclusters,Cluster* CCE, int nspin){
 
    for (int n=CCE->order; n>=1; n--){
        HashCluster *existingClusters = findCluster(*hashclusters,n);
        Property *existingCluster, *tmp;

        int nClusInfo = 1;

        HASH_ITER(hh,existingClusters->prop, existingCluster, tmp)
        {
            int* cluster = existingCluster->spins;
            //printf("MAIN CLUSTER : %s \n",existingCluster->id);
            int cluster_count = existingCluster->count;
            if (cluster_count != 0){
                //////////////////////////////////////////
                // The cluster that you will find its sub-clusters
                // allocate memory for clusinfo
                CCE->clusinfo[n] = (int**)realloc(CCE->clusinfo[n],(nClusInfo+1)*sizeof(int*));
                CCE->clusinfo[n][nClusInfo] = (int*)calloc(n+1,sizeof(int));
                // value assignment
                CCE->clusinfo[n][nClusInfo][0] = cluster_count;
                for (int i=1; i<n+1; i++){
                    // Here i'm transfering the hashcluster to clusinfo
                    // cluster[n] = {ispin, jspin, kspin, ..., nspin}
                    // For -1 of cluster[i-1] is the following:
                    //   the clusinfo index start from 1 but cluster index start from 0
                    // For -1 of cluster[i-1] "-1" is the following:
                    //   -1 for actual spin index(BathArray) is starting from 0
                    CCE->clusinfo[n][nClusInfo][i] = cluster[i-1]; 
                    
                }
                // The number of clusters
                CCE->clusinfo[n][0][0] = nClusInfo+1;
                //////////////////////////////////////////
                // The 0-th cluster (qubit)
                // L0 (L1' L2' L3' L4') (L12' L13' L24') (L123' L124') (' = tilde)
                // L123' = L123 / ( (L0) (L1' L2' L3') (L12' L13' L23')) -> L0 count = -1
                // L124' = L124 / ( (L0) (L1' L2' L4') (L12' L14' L24')) -> L0 count = -1
                // L0^? (L1'^-1 L2'^-1 ) (L12'^-1 L14'^-1 L23'^-1) L123 L124
                // L23'^-1 = (L23 / ( (L0) (L2' L3')))^-1 -> L0 count = 1
                // L14'^-1 = (L14 / ( (L0) (L1' L4')))^-1 -> L0 count = 1
                // L12'^-1 = (L12 / ( (L0) (L1' L2')))^-1 -> L0 count = 1
                // L0^? (L1'^1 L2'^1 L3'^1 L4'^1 ) L23 L14 L12 L123 L124
                // L0^? (L1 L2 L3 L4 ) L23 L14 L12 L123 L124 -> L0 count = -1, -1, -1, -1
                // Total L0 count = -2
                // CCE->clusinfo[0][0][0] += -1 * cluster_count;
                nClusInfo++;
            }
            //////////////////////////////////////////
            for (int r = n-1; r >=1; r--){
                ///////////////////////////////////////////
                // Find subclusters
                ///////////////////////////////////////////
                // Calculate the number of combinations (n choose r)
                int numCombinations = 1;
                for (int i = 0; i < r; i++) {
                    numCombinations *= (n - i);
                    numCombinations /= (i + 1);
                }
                // Allocate memory for data
                int** subclusters = (int**)malloc(numCombinations * sizeof(int*));
                for (int i = 0; i < numCombinations; i++) {
                    subclusters[i] = (int*)malloc(r * sizeof(int));
                }
                int* tempCombination = (int*)malloc(r * sizeof(int)); // Temporary storage for combinations
                int nCombination = 0 ;
                generateCombinations(cluster, &subclusters, tempCombination, 0, n, 0, r, &nCombination);
                // print subclusters
                // Free memory for tempCombination
                free(tempCombination); 
                ///////////////////////////////////////////
                // Add subclusters
                ///////////////////////////////////////////
                for (nCombination=0; nCombination<numCombinations; nCombination++){
                    // sort
                    QuickSort(&(subclusters[nCombination]),0,r-1);
                    // id 
                    int maxlenId = setMaxLengthStr(nspin,r);
                    char* id = (char*)calloc(maxlenId,sizeof(char));
                    typeStr(&id,subclusters[nCombination],r,nspin);
                    // count
                    Property* subcluster;
                    subcluster = findProperty(*hashclusters,r,id);
                    //printf("subcluster : %s \n",id);
                    if (subcluster!=NULL){
                        subcluster->count = subcluster->count - cluster_count;
                    }
                    else{
                        if (CCE->addsubclus){
                            printf("Error, subcluster = null\n");
                            exit(1);
                        }
                        else{;
                            int count = 1;
                            int res = 0;
                            float strength = 0.0; //addAllStrengthForAllSpins(stmap, subclusters[nCombination], r);
                            count = count - cluster_count;
                            res = addCluster(hashclusters,r,id,subclusters[nCombination],strength,count);
                        }
                    }
                    //free(id); 
                    //free(subclusters[nCombination]);
                }
                //free(subclusters);
            }           
        } 
    }
}