#include "../include/utilities.h"
#include "../include/memory.h"
#include <float.h>  // FLT_EPSILON
#include "mpi.h"

/* utils ---------------------------------------------------*/
MatrixXcd Double1dToMatrixXcd(double* val, int n){
    MatrixXcd mat = MatrixXcd::Zero(n,1);
    for (int i=0; i<n; i++){
        mat(i,0) = doublec(val[i],0.0);
    }
    return mat;
}

double* MatrixXcdToDouble1d(MatrixXcd mat){
    double* val = allocDouble1d(mat.rows());
    for (int i=0; i<mat.rows(); i++){
        val[i] = mat(i,0).real();
    }
    return val;
}

/* math functions -----------------------------------------------*/

/////Calculate the parameter///////////////////////////////////
double dist(double spin1[],double spin2[]){
	double dist = 0;
	dist = sqrt(pow(spin1[0] - spin2[0],2) + pow(spin1[1] - spin2[1],2) + pow(spin1[2] - spin2[2],2));
	return dist;
}

/////polar angle from the vector_m to vector_n/////////////////
double cosTheta(double spin1[], double spin2[],double dist){ 
 //Note that spin_m[0] = gyromagnetic ratio, spin_m[0 ~ 2] --> x, y, z coordinates
    return (spin2[2] - spin1[2]) / dist;
}

double sinTheta(double spin1[], double spin2[], double dist){
    double tmp_mn = sqrt(pow(spin2[0] - spin1[0], 2) + pow(spin2[1] - spin1[1], 2)); // x^2 + y^2
    return tmp_mn / dist;
}

/////azimuthal angle from the vector_m to vector_n/////////////
double cosPhi(double spin1[], double spin2[]){
    double tmp_mn = sqrt(pow(spin2[0] - spin1[0], 2) + pow(spin2[1] - spin1[1], 2)); // x^2 + y^2
    if (tmp_mn == 0) return 1.0;
    else return (spin2[0] - spin1[0]) / tmp_mn;
}
double sinPhi(double spin1[], double spin2[]){
    double tmp_mn = sqrt(pow(spin2[0] - spin1[0], 2) + pow(spin2[1] - spin1[1], 2)); // x^2 + y^2
    if (tmp_mn == 0) return 0.0;
    else return (spin2[1] - spin1[1]) / tmp_mn;
}

///////////////////////////////////////////////////////////////
// Physics functions
///////////////////////////////////////////////////////////////

// Obtain possible substates of the spin
float* substates(float S){
    float* ms=NULL;
    int n = (int)(2*S+1);
    ms = allocFloat1d(n);
    for (int i=0; i<n; i++){
        ms[i] = S-i;
    }
    return ms;
}

// Obtain the spinor array
MatrixXcd getSpinor(float S, float ms){

    if (ms > S || ms < -S){
        perror("Error(spinorArray): ms is not satisfying -S <= ms <= S");
        exit(EXIT_FAILURE);
    }

    int n = (int)(2*S+1);
    MatrixXcd spinor = MatrixXcd::Zero(n,1);
    //e.g. S=1
    //    i = 0, 1, 2 (ms = S-i)
    //    i = 0 -> ms = 1 -> arr = [ 1 0 0] (else = 0)
    //    i = 1 -> ms = 0 -> arr = [ 0 1 0] (else = 0)
    //    i = 2 -> ms = -1 -> arr = [ 0 0 1] (else = 0)
    bool chkms = false;
    for (int i=0; i<n; i++){
        float mstmp = S-i;
        if (mstmp == ms){
            spinor(i,0) = doublec(1.0,0.0);
            chkms = true;
        }
        else{
            spinor(i,0) = doublec(0.0,0.0);
        }
    }

    if (chkms){
        // std::cout << "spinor" << std::endl;
        // std::cout << spinor << std::endl;
        return spinor;
    }
    else{
        perror("Error(spinorArray): ms is not the sub level of S");
        exit(EXIT_FAILURE);
    }
}


MatrixXcd kron(MatrixXcd A, MatrixXcd B){
    MatrixXcd C(A.rows()*B.rows(),A.cols()*B.cols());
    for (int i=0; i<A.rows(); i++){
        for (int j=0; j<A.cols(); j++){
            C.block(i*B.rows(),j*B.cols(),B.rows(),B.cols()) = A(i,j)*B;
        }
    }
    return C;
}


/*! 
 * @brief Calculate the norm of the matrix
 * @details Spinor Properties : normalized set
 * @param[in] m gsl_matrix_complex
 * @return norm of the matrix
 */
double calNorm(MatrixXcd m){
    int row = m.rows();
    int col = m.cols();

    // Compare only vector
    if (col!=1){
        perror("Error(norm): the matrix is not a vector");
        exit(EXIT_FAILURE);
    }

    // Find norm^2
    doublec norm2 = (m.adjoint() * m)(0,0);
    
    // check if the norm is real number
    if (fabs(norm2.imag()) > FLT_EPSILON){
        perror("Error(norm): the norm is not a real number");
        exit(EXIT_FAILURE);
    }

    return sqrt(norm2.real());
}


int normalize(MatrixXcd* m){

    double norm = calNorm(*m);

    if (fabsf((float)norm -1.0f) <= FLT_EPSILON){
        return 0; // already normalized
    }
    else{
        *m = *m / norm;
        return 1; // give normalized matrix
    }

}

float findZbasisSubLevel(MatrixXcd spinor){
    int n = spinor.rows();
    float S = (float)(n-1)/2.0;
    float ms = 0.0;
    
    bool chkms = false;
    for (int i=0; i<n; i++){
        double spinoriReal = spinor(i,0).real();
        if (fabs(spinoriReal - 1.0) <= FLT_EPSILON){
            ms = S-i;
            if (!chkms){
                chkms = true;
            }
            else{
                perror("Error(findZbasisSubLevel): ms is not unique");
                exit(EXIT_FAILURE);
            }
        }
    }

    if (chkms){
        return ms;
    }
    else{
        perror("Error(findMs): ms is not the sub level of S");
        exit(EXIT_FAILURE);
    }
}

MatrixXcd partialtrace(MatrixXcd Mij, int dimrow, int dimcol){

    // Mij : Full Matrix
    // Bkl : Block Matrix
    // Rkl : Reduced Matrix

    //        | B00 B01 B02 | 
    // Mij =  | B10 B11 B12 | 
    //        | B20 B21 B22 | 
    MatrixXcd Bkl = MatrixXcd::Zero(dimrow,dimcol);

    //        | r00 r01 r02 |
    // Rkl =  | r10 r11 r12 | = PartialTrace[Mij]
    //        | r20 r21 r22 |
    MatrixXcd Rkl = MatrixXcd::Zero(Mij.rows()/dimrow,Mij.cols()/dimcol);
    
    // Get partial trace of M after finding block matrix
    int irow = 0;
    int icol = 0;
    for (int k=0; k<Mij.rows(); k+=dimrow){
        icol=0;
        for (int l=0; l<Mij.cols(); l+=dimcol){
            // Find block matrix Bkl
            // Bkl,pq = Mij 
            // * i = (a)*dimrow + p < (a+1) * dimrow
            // * j = (b)*dimcol + q < (b+1) * dimcol
            Bkl = Mij.block(k,l,dimrow,dimcol);
            // Trace of Bkl 
            Rkl(irow,icol) = Bkl.trace();
            icol++;
        }
        irow++;
    }
    return Rkl;
}

MatrixXcd powMatrixXcdElementWise(MatrixXcd a, int n){
    // NOTE!!! 
    // here the input power is only allowing "int"
    // Because i found float n gives wrong value in this case.
    // have to use int n
    // but n<0 is okay to use
    int nrow = a.rows();
    int ncol = a.cols();
    MatrixXcd res(nrow,ncol);
 
    for (int i=0; i<nrow; i++){
        for (int j=0; j<ncol; j++){
            res(i,j) = std::pow(a(i,j),n);
        }
    }
    return res;
}

MatrixXcd mulMatrixXcdElementWise(MatrixXcd a, MatrixXcd b){
    return a.array() * b.array();
}

/* Easy print ---------------------------------------------------*/
void printInlineMatrixXcd(char* key, MatrixXcd mat){
    printf("      %-18s:   [ ", key);

    if (mat.rows() * mat.cols() > 9){
        int count = 0;
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                std::complex<double> z = mat(i, j);
                printf("%3.2fj%+-3.2f", z.real(), z.imag());
                if (i != mat.rows() - 1 || j != mat.cols() - 1) {
                    printf(", ");
                }
                count++;
            }
            if (count%9 == 0 && count != (mat.rows() * mat.cols())){
                printf("\n%30s"," ");
            }          
        }
    }else{
        for (int i = 0; i < mat.rows(); ++i) {
            for (int j = 0; j < mat.cols(); ++j) {
                std::complex<double> z = mat(i, j);
                printf("%3.2fj%+-3.2f", z.real(), z.imag());
                if (i != mat.rows() - 1 || j != mat.cols() - 1) {
                    printf(", ");
                }
            }
        }
    }
    printf(" ]\n");
}

void printStructElementChar(char* key,char* val){
    printf("      %-18s:   %-21s\n", key, val);
}

void printStructElementChar2d(char* key, char** val, int n){
    printf("      %-18s:   [ ", key);
    for (int i = 0; i < n; i++){
        printf("%-10s", val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
}

void printStructElementInt(char* key, int val){
    printf("      %-18s:   %-21d\n", key, val);
}

void printStructElementInt1dIdx(char* key, int* val, int n){
    printf("      %-18s:   [ ", key);
    for (int i = 0; i < n; i++){
        printf("%3d : %-5d",i,val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
}

void printStructElementFloat(char* key, float val){
    printf("      %-18s:   %-21.6g\n", key, val);
}

void printStructElementFloat1d(char* key, float* val, int n){
    printf("      %-18s:   [ ", key);
    for (int i = 0; i < n; i++){
        printf("%-10.2f", val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
}

void printStructElementDouble(char* key, double val){
    printf("      %-18s:   %-21.3g\n", key, val);
}

void printStructElementDouble1d(char* key, double* val, int n){
    printf("      %-18s:   [ ", key);
    for (int i = 0; i < n; i++){
        printf("%-10.2f", val[i]);
        if (i != n - 1){
            printf(", ");
        } else {
            printf(" ]\n");
        }
    }
}

void printStructElementBool(char* key, bool val){
    printf("      %-18s:   %-21s\n", key, val ? "true" : "false");
}

void printLine(){
    printf("      -------------------------------------------------\n");
}

void printLineSection(){
    printf("\n    ===============================================================\n\n");
}

void printTitle(char* title){
    printf("    < %s > \n\n",title);
}

void printSubTitle(char* title){
    printf("\n    >> %s\n\n",title);
}

void printMessage(char* title){
    printf("        %s\n",title);
}

/* Find index ---------------------------------------------------*/
int findIndexInt(int* array, int ista, int iend, int val){
    for (int i = ista; i <= iend; i++){
        if (array[i] == val){
            return i;
        }
    }
    return -1;
}

int findIndexCharFix(char array[][MAX_CHARARRAY_LENGTH], int ista, int iend, char* val){
    for (int i = ista; i <= iend; i++){
        if (strcasecmp(array[i],val) == 0){
            return i;
        }
    }
    return -1;
}

int findIndexChar(char** array, int ista, int iend, char* val){
    for (int i = ista; i <= iend; i++){
        if (strcasecmp(array[i],val) == 0){
            return i;
        }
    }
    return -1;
}

/* Quick sort ---------------------------------------------------*/

void QuickSort(int** d_Array, int left, int right){
    if(left<=right){
        int pivot = Partition(d_Array, left, right);
        QuickSort(d_Array, left, pivot-1);
        QuickSort(d_Array, pivot+1, right);
    }
}

int Partition(int** d_Array, int left, int right){
    int pivot = (*d_Array)[left];
    int low = left +1;
    int high = right;
    while(low<=high){
        while(low<=right && pivot >= (*d_Array)[low]){
            low++;
        }
        while(high>=(left+1) && pivot <= (*d_Array)[high]){
            high--;
        }
        if (low <= high){
            Swap(d_Array, low, high);
        }
    }
    Swap(d_Array, left, high);
    return high;
}

void Swap(int** d_Array, int a, int b){
    int temp = (*d_Array)[a];
    (*d_Array)[a] = (*d_Array)[b];
    (*d_Array)[b] = temp;
}

/* Type checker ---------------------------------------------------*/

//Check the string is the double
int isStringDouble(char *s){
    //0--> no number, 1--> is number
    size_t size = strlen(s);
    if (size == 0){return 0;} //no number
    if (size == 1){if(!isdigit(s[0])){return 0;}}

    bool CheckConti=false;
    for (int i=0;i<(int) size; i++){
        if(s[i] == '.' || s[i] =='-' || s[i] == '+' || isspace(s[i])){
            if (CheckConti==true){
                return 0; //if special string is continue, the string is not number
            }else{
                CheckConti=true;continue;
            }
        } // ok 
        if(s[i] < '0' || s[i] > '9') {return 0;} //not number
        CheckConti=false; //for check that continue special string
    }

    return 1; //this string is number
}

/* MPI ---------------------------------------------------*/
void para_range(int n1,int n2, int nprocs, int myrank, int*ista, int *iend)
{
    int iwork1, iwork2;
    iwork1 = (n2-n1+1)/nprocs;
    iwork2 = (n2-n1+1) % nprocs;
    *ista = myrank*iwork1 + n1 + min(myrank, iwork2);
    *iend = *ista + iwork1 - 1;
    if(iwork2 > myrank) *iend = *iend + 1;
}

int min(int x, int y)
{
    int v;
    if (x>=y) v = y;
    else v = x;
    return v;
}


int*** MPI_getLocalClusters(int order, int*** clusters){

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

// MatrixXcd* MPI_reduceLocalResult(int nstep, MatrixXcd* local){

//     int dim = local[0].rows();

//     // Initialize the result variable   
//     MatrixXcd* result = new MatrixXcd[nstep];
//     for (int istep=0; istep<nstep; istep++){
//         result[istep] = MatrixXcd::Constant(dim,dim,doublec(1.0,0.0));
//     }
//     MPI_Barrier(MPI_COMM_WORLD);

//     // Local data reduce to root

//     if (rank==0){printf("Start MPI_Reduce ... \n");}
//     int err;
//     for (int istep=0; istep<nstep; istep++){
//         //if (rank==0){printf("Start rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         //if (rank==10){printf("Start rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         err = MPI_Reduce(local[istep].data(),result[istep].data(),dim*dim,MPI_DOUBLE_COMPLEX,MPI_PROD,0,MPI_COMM_WORLD);
//         //if (rank==0){printf("End rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         //if (rank==10){printf("End rank %d : MPI_Reduce istep = %d\n ... ",rank,istep);}
//         MPI_Barrier(MPI_COMM_WORLD);
//     }

//     if (rank==0 && err == MPI_SUCCESS){printf("Succeed MPI_Reduce\n ... ");}
//     MPI_Barrier(MPI_COMM_WORLD);

//     if (rank!=0){
//         delete[] result;
//         return NULL;
//     }else{
//         return result;
//     }

//     //return result;
// } 