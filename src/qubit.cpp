#include "../include/qubit.h"
#include "../include/utilities.h"
#include "../include/memory.h"
#include "../include/hamiltonian.h"
#include <errno.h>
#include <float.h>

/* High level --------------------------------------------------------*/

// init
QubitArray* QubitArray_init(){

    QubitArray* qa = (QubitArray*)allocArray1d(1,sizeof(QubitArray));

    QubitArray_set_alphaidx(qa,NULL);
    QubitArray_set_betaidx(qa,NULL); 

    QubitArray_setNqubit(qa,0);
    QubitArray_setOverhaus(qa,false);
    QubitArray_setQubit(qa,NULL);

    qa->intmap = NULL;

    // MatrixXcd psia; do not need to initialize
    // MatrixXcd psib; do not need to initialize
    // MatrixXcd psi0; do not need to initialize
    
    return qa;

}

// free
void QubitArray_freeAll(QubitArray* qa){
    QubitArray_freeQubit(qa);
    QubitArray_freeIntmap(qa);
    QubitArray_free_alphaidx_betaidx(qa);
    freeArray1d(qa);
}

// phyiscal properties
int QubitArray_dim(QubitArray *qa){
    int dim = 1;
    int nqubit = QubitArray_getNqubit(qa);
    for (int i=0; i<nqubit; i++){
        dim *= QubitArray_dimQubit_i(qa,i);
    }
    return dim;
}

int QubitArray_dimQubit_i(QubitArray *qa, int i){
    float S = QubitArray_getQubit_i_spin(qa,i);
    return (int)(2*S+1);
}

double QubitArray_mindist(double* xyz, QubitArray* qa){

    // mininum distance from the qubit array
    int nqubit = QubitArray_getNqubit(qa);
    double* q0_xyz = QubitArray_getQubit_i_xyz(qa,0);
    
    double r_min = dist(xyz, q0_xyz);

    for (int i=0; i<nqubit; i++){
        double* qi_xyz = QubitArray_getQubit_i_xyz(qa,i);
        double r = dist(xyz, qi_xyz);
        if (r < r_min){
            r_min = r;
        }
    }
    return r_min;
}

// From Qubit[i]->alpha, beta, set Psia, Psib by krnoecker product
void QubitArray_setPsiaPsib_fromQubit(QubitArray* qa){

    int nqubit = QubitArray_getNqubit(qa);
    MatrixXcd psia;
    MatrixXcd psib;

    for (int i=0; i<nqubit; i++){
        MatrixXcd alpha_i = QubitArray_getQubit_i_alpha(qa,i);
        MatrixXcd beta_i = QubitArray_getQubit_i_beta(qa,i);
        if (i == 0){
            psia = alpha_i;
            psib = beta_i;
        }else{
            psia = kron(psia,alpha_i);
            psib = kron(psib,beta_i);
        }
    }
    QubitArray_setPsia(qa,psia);
    QubitArray_setPsib(qa,psib);
}

void QubitArray_setPsiaPsib_fromIdx(QubitArray* qa, float* bfield){

    int* alphaidx = QubitArray_get_alphaidx(qa);
    int* betaidx = QubitArray_get_betaidx(qa);

    if (alphaidx == NULL && betaidx == NULL){
        //error message
        fprintf(stderr,"Error: QubitArray_setPsiaPsib_fromIdx: alphaidx and betaidx are not set\n");
        exit(EXIT_FAILURE);
    }

    // Hamiltonian
    MatrixXcd** sigmas = QubitArray_PauliOperators(qa);
    MatrixXcd Hq = QubitArray_TotalHamil(qa,sigmas,bfield);

    for (int i=0; i<QubitArray_getNqubit(qa); i++){
        delete[] sigmas[i];
    }
    delete[] sigmas;
    
    // Eigenvalues and Eigenvectors
    int dim = QubitArray_dim(qa);

    Eigen::SelfAdjointEigenSolver<MatrixXcd> es(dim);
    es.compute(Hq);

    MatrixXcd eigenVectors = es.eigenvectors();
    VectorXcd eigenValues = es.eigenvalues();

    // Sort in desending order
    int* idx = getIndexInOrder(eigenValues);
    VectorXcd sortedEigenValues = sortEigenValues(eigenValues,idx);
    MatrixXcd sortedEigenVectors = sortEigenVectors(eigenVectors, idx);

    // Get Psia, Psib
    MatrixXcd psia = sortedEigenVectors.col(*alphaidx);
    MatrixXcd psib = sortedEigenVectors.col(*betaidx);

    // Normalize
    int res;
    res = normalize(&psia);
    res = normalize(&psib);

    // Set Psia, Psib
    QubitArray_setPsia(qa,psia);
    QubitArray_setPsib(qa,psib);
}

// From Psia & Psib, set psi0 by adding them
void QubitArray_setPsi0_fromPsiaPsib(QubitArray* qa){
    
    MatrixXcd psi0 = QubitArray_getPsi0(qa);
    if (psi0.rows()!=0 || psi0.cols()!=0){
        // warning message
        fprintf(stderr,"Warning: QubitArray_setPsi0_fromPsiaPsib: psi0 is already set\n");
    }

    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    psi0 = psia + psib;
    int res = normalize(&psi0);

    QubitArray_setPsi0(qa,psi0);
}

// utils
int QubitArray_getQubitIdx_fromName(QubitArray* qa, const char* name){
    int nqubit = QubitArray_getNqubit(qa);
    for (int i=0; i<nqubit; i++){
        if (strcmp(QubitArray_getQubit_i_name(qa,i),name) == 0){
            return i;
        }
    }
    return -1;
}

MatrixXcd** QubitArray_PauliOperators(QubitArray* qa){

    int nqubit = QubitArray_getNqubit(qa);
    int dimQA = QubitArray_dim(qa);

    // Pauli operators
    MatrixXcd** sigmas = new MatrixXcd*[nqubit];
    for (int iq=0; iq<nqubit; iq++){
        int dimQ = QubitArray_dimQubit_i(qa,iq);
        sigmas[iq] = getPauliOperators(dimQ);
    }

    return sigmas;
}

MatrixXcd* QubitArray_PauliOperator_fromPsiaPsib(QubitArray* qa){

    MatrixXcd psia = QubitArray_getPsia(qa);
    MatrixXcd psib = QubitArray_getPsib(qa);

    // sigmas:
    // sigmas[0] = I, sigmas[1] = X, sigmas[2] = Y, sigmas[3] = Z
    MatrixXcd* sigmas = getGeneralPauliOperators(psia,psib);

    return sigmas;
}

MatrixXcd QubitArray_TotalHamil(QubitArray* qa, MatrixXcd** sigmas, float* bfield){

    int nqubit = QubitArray_getNqubit(qa);
    int dimQA = QubitArray_dim(qa);

    // Qubit Hamiltonian
    MatrixXcd HqTotal = MatrixXcd::Zero(dimQA,dimQA);

    // Single Hamiltonian
    for (int iq=0; iq<nqubit; iq++){
        MatrixXcd Hqi = QubitArray_SingleHamil(qa,sigmas,iq,bfield);
        MatrixXcd HqiExpand = expandHamiltonian(sigmas, Hqi, nqubit, iq);
        HqTotal += HqiExpand;
    }

    // Interaction Hamiltonian
    for (int i=0; i<nqubit; i++){
        MatrixXcd Hqij;
        for (int j=i+1; j<nqubit; j++){
            Hqij = QubitArray_InteractionHamil(qa,sigmas,i,j);
            HqTotal += Hqij;
        }
    }
    
    return HqTotal;
}

MatrixXcd QubitArray_SingleHamil(QubitArray* qa, MatrixXcd** sigmas, int iq, float* bfield){

    // Single Hamiltonian
    MatrixXcd Hqi_Zeeman = QubitArray_ZeemanHamil(qa,sigmas,iq,bfield);    
    MatrixXcd Hqi_Detuning = QubitArray_DetuningHamil(qa,sigmas,iq);
    MatrixXcd Hqi_Overhaus = QubitArray_OverhausHamil(qa,sigmas,iq);
    MatrixXcd Hqi_ZFS = QubitArray_ZFSHamil(qa,sigmas,iq);

    return Hqi_Zeeman + Hqi_Detuning + Hqi_Overhaus + Hqi_ZFS;
}

MatrixXcd QubitArray_ZeemanHamil(QubitArray* qa, MatrixXcd** sigmas, int iq, float* bfield){

    int nqubit = QubitArray_getNqubit(qa);

    if (iq >= nqubit){
        fprintf(stderr,"Error: QubitArray_ZeemanHamil: iq = %d is out of range\n",iq);
        exit(EXIT_FAILURE);
    }
    
    float qgyro = QubitArray_getQubit_i_gyro(qa,iq);
    MatrixXcd vecZeeman = calZeemanVector(qgyro,bfield);

    return calHamiltonianSingleInt(vecZeeman,sigmas[iq]);
}

MatrixXcd QubitArray_DetuningHamil(QubitArray* qa, MatrixXcd** sigmas, int iq){

    int nqubit = QubitArray_getNqubit(qa);

    if (iq >= nqubit){
        fprintf(stderr,"Error: QubitArray_DetuningHamil: iq = %d is out of range\n",iq);
        exit(EXIT_FAILURE);
    }

    double detun = QubitArray_getQubit_i_detuning(qa,iq);
    MatrixXcd vecDetun = calDetuningVector(detun);

    return calHamiltonianSingleInt(vecDetun,sigmas[iq]);
}

MatrixXcd QubitArray_OverhausHamil(QubitArray* qa, MatrixXcd** sigmas, int iq){

    int nqubit = QubitArray_getNqubit(qa);

    if (iq >= nqubit){
        fprintf(stderr,"Error: QubitArray_OverhausHamil: iq = %d is out of range\n",iq);
        exit(EXIT_FAILURE);
    }

    double overhaus = QubitArray_getQubit_i_overhaus(qa,iq);
    MatrixXcd vecOverhaus = calOverhauserVector(overhaus);

    return calHamiltonianSingleInt(vecOverhaus,sigmas[iq]);
}

MatrixXcd QubitArray_ZFSHamil(QubitArray* qa, MatrixXcd** sigmas, int iq){

    int nqubit = QubitArray_getNqubit(qa);

    if (iq >= nqubit){
        fprintf(stderr,"Error: QubitArray_ZFSHamil: iq = %d is out of range\n",iq);
        exit(EXIT_FAILURE);
    }

    MatrixXcd zfs = QubitArray_getIntmap_i_j(qa,iq,iq);
    MatrixXcd HqSelf = calHamiltonianSelfInt(sigmas[iq],zfs);

    return HqSelf;
}

MatrixXcd QubitArray_InteractionHamil(QubitArray* qa, MatrixXcd** sigmas, int iq, int jq){

    int nqubit = QubitArray_getNqubit(qa);
    int dimQA = QubitArray_dim(qa);

    // Qubit Hamiltonian
    MatrixXcd Hqij = MatrixXcd::Zero(dimQA,dimQA);

    if (iq<jq){
        MatrixXcd tensor = QubitArray_getIntmap_i_j(qa,iq,jq);
        bool isEmpty = tensor.isZero(FLT_EPSILON);

        if (!isEmpty){
            Hqij = calHamiltonianHeteroInt(sigmas, tensor, nqubit, iq, jq);
        }
    }else{
        fprintf(stderr,"Error: QubitArray_InteractionHamil: iq,jq = %d,%d is out of range or iq>jq\n",iq,jq);
        exit(EXIT_FAILURE);
    }

    return Hqij;
}

// density matrix
MatrixXcd   QubitArray_Rho0(QubitArray* qa){
    MatrixXcd psi0 = QubitArray_getPsi0(qa);
    return psi0 * psi0.adjoint();
}

/* -------------------------------------------------------------------*/
/*                                                                    */
/* Low level                                                          */
/*                                                                    */
/* -------------------------------------------------------------------*/

////////////////////////////////////////////////////////////////
// QubitArray properties
////////////////////////////////////////////////////////////////
// alloc

void QubitArray_allocQubit(QubitArray* qa){
    // alloc as much as nqubit
    int nqubit = QubitArray_getNqubit(qa);
    qa->qubit = (Qubit**)allocArray2d(nqubit,1,sizeof(Qubit));
}

void QubitArray_allocIntmap(QubitArray* qa){
    // alloc as much as nqubit*nqubit
    int nqubit = QubitArray_getNqubit(qa);
    qa->intmap = new MatrixXcd*[nqubit];
    for (int i=0; i<nqubit; i++){
        qa->intmap[i] = new MatrixXcd[nqubit];
        for (int j=0; j<nqubit; j++){
            qa->intmap[i][j] = MatrixXcd::Zero(3,3);
        }
    }
}

void QubitArray_alloc_alphaidx_betaidx(QubitArray* qa){
    qa->_alphaidx = (int*)allocArray1d(1,sizeof(int));
    qa->_betaidx = (int*)allocArray1d(1,sizeof(int));
}

// set
void QubitArray_setNqubit(QubitArray* qa, const int nqubit){
    qa->nqubit = nqubit;
}

void QubitArray_setOverhaus(QubitArray* qa, const bool overhaus){
    qa->overhaus = overhaus;
}

void QubitArray_set_alphaidx(QubitArray* qa, const int* alphaidx){
    if (alphaidx == NULL){
        qa->_alphaidx = NULL;
    }else{
        *(qa->_alphaidx) = *alphaidx;
    }
}

void QubitArray_setQubit(QubitArray* qa, Qubit** qubits){
    qa->qubit = qubits;
}

void QubitArray_set_betaidx(QubitArray* qa, const int* betaidx){
    if (betaidx == NULL){
        qa->_betaidx = NULL;
    }else{
        *(qa->_betaidx) = *betaidx;
    }
}

void QubitArray_setIntmap_i_j(QubitArray* qa, const MatrixXcd tensor, int i, int j){
    int nqubit = QubitArray_getNqubit(qa);
    if (i >= nqubit || j >= nqubit || i > j || i < 0 || j < 0){
        fprintf(stderr,"Error: QubitArray_setIntmap_i_j: i,j = %d,%d is out of range or i>j\n",i,j);
        exit(EXIT_FAILURE);
    }
    qa->intmap[i][j] = tensor;
}

void QubitArray_setPsia(QubitArray* qa, const MatrixXcd psia){ // auto normalized
    qa->psia = psia.normalized();
}

void QubitArray_setPsib(QubitArray* qa, const MatrixXcd psib){ // auto normalized
    qa->psib = psib.normalized();
}

void QubitArray_setPsi0(QubitArray* qa, const MatrixXcd psi0){ // auto normalized
    qa->psi0 = psi0.normalized();
}


// get
int QubitArray_getNqubit(const QubitArray* qa){
    return qa->nqubit;
}

bool QubitArray_getOverhaus(const QubitArray* qa){
    return qa->overhaus;
}

int* QubitArray_get_alphaidx(const QubitArray* qa){
    return qa->_alphaidx;
}

int* QubitArray_get_betaidx(const QubitArray* qa){
    return qa->_betaidx;
}

MatrixXcd** QubitArray_getIntmap(const QubitArray* qa){
    return qa->intmap;
}

MatrixXcd QubitArray_getIntmap_i_j(const QubitArray* qa, int i, int j){
    return qa->intmap[i][j];
}

MatrixXcd QubitArray_getPsia(const QubitArray* qa){
    return qa->psia;
}

MatrixXcd QubitArray_getPsib(const QubitArray* qa){
    return qa->psib;
}

MatrixXcd QubitArray_getPsi0(const QubitArray* qa){
    return qa->psi0;
}

// free

void QubitArray_freeQubit(QubitArray* qa){
    if (qa->qubit == NULL){
        return;
    }
    int nqubit = QubitArray_getNqubit(qa);
    freeArray2d((void**)qa->qubit,nqubit);
}

void QubitArray_freeIntmap(QubitArray* qa){
    if (qa->intmap == NULL){
        return;
    }
    int nqubit = QubitArray_getNqubit(qa);
    for (int i=0; i<nqubit; i++){
        delete[] qa->intmap[i];
    }
    delete[] qa->intmap;
}

void QubitArray_free_alphaidx_betaidx(QubitArray* qa){
    freeArray1d(qa->_alphaidx);
    freeArray1d(qa->_betaidx);
}

////////////////////////////////////////////////////////////////
// i-th qubit properties
////////////////////////////////////////////////////////////////

//set 
void QubitArray_setQubit_i_name(QubitArray* qa, const char* name, int i){
    strcpy(qa->qubit[i]->name,name);
}

void QubitArray_setQubit_i_spin(QubitArray* qa, const float spin, int i){
    qa->qubit[i]->spin = spin;
}

void QubitArray_setQubit_i_gyro(QubitArray* qa, const double gyro, int i){
    qa->qubit[i]->gyro = gyro;
}

void QubitArray_setQubit_i_xyz(QubitArray* qa, const double* xyz, int i){
    qa->qubit[i]->xyz[0] = xyz[0];
    qa->qubit[i]->xyz[1] = xyz[1];
    qa->qubit[i]->xyz[2] = xyz[2];
}

void QubitArray_setQubit_i_detuning(QubitArray* qa, const double detuning, int i){
    qa->qubit[i]->detuning = detuning;
}

void QubitArray_setQubit_i_overhaus(QubitArray* qa, const double overhaus, int i){
    qa->qubit[i]->overhaus = overhaus;
}

void QubitArray_setQubit_i_alpha(QubitArray* qa, const MatrixXcd alpha, int i){
    qa->qubit[i]->alpha = alpha;
}

void QubitArray_setQubit_i_beta(QubitArray* qa, const MatrixXcd beta, int i){
    qa->qubit[i]->beta = beta;
}

void QubitArray_setQubit_i_alpha_fromMs(QubitArray* qa, const float ms, int i){
    float S = QubitArray_getQubit_i_spin(qa,i);
    QubitArray_setQubit_i_alpha(qa,getSpinor(S,ms),i);
}

void QubitArray_setQubit_i_beta_fromMs(QubitArray* qa, const float ms, int i){
    float S = QubitArray_getQubit_i_spin(qa,i);
    QubitArray_setQubit_i_beta(qa,getSpinor(S,ms),i);
}

// get
Qubit* QubitArray_getQubit_i(const QubitArray* qa, int i){
    return qa->qubit[i];
}

char* QubitArray_getQubit_i_name(const QubitArray* qa, int i){
    return qa->qubit[i]->name;
}

float QubitArray_getQubit_i_spin(const QubitArray* qa, int i){
    return qa->qubit[i]->spin;
}

double QubitArray_getQubit_i_gyro(const QubitArray* qa, int i){
    return qa->qubit[i]->gyro;
}

double* QubitArray_getQubit_i_xyz(const QubitArray* qa, int i){
    return qa->qubit[i]->xyz;
}

MatrixXcd QubitArray_getQubit_i_alpha(const QubitArray* qa, int i){
    return qa->qubit[i]->alpha;
}

MatrixXcd QubitArray_getQubit_i_beta(const QubitArray* qa, int i){
    return qa->qubit[i]->beta;
}

double QubitArray_getQubit_i_detuning(const QubitArray* qa, int i){
    return qa->qubit[i]->detuning;
}

double QubitArray_getQubit_i_overhaus(const QubitArray* qa, int i){
    return qa->qubit[i]->overhaus;
}


/* Report ---------------------------------------------------------*/
void QubitArray_reportQubit_i(QubitArray* qa, int i){

    printStructElementChar("qubit name",QubitArray_getQubit_i_name(qa,i));
    printStructElementFloat("spin",QubitArray_getQubit_i_spin(qa,i));
    printStructElementDouble("gyro (radkHz/G)",QubitArray_getQubit_i_gyro(qa,i));
    printStructElementDouble("detuning (kHz)",QubitArray_getQubit_i_detuning(qa,i));
    printStructElementDouble("overhaus (kHz)",QubitArray_getQubit_i_overhaus(qa,i));
    printStructElementDouble1d("xyz (A)",QubitArray_getQubit_i_xyz(qa,i),3);

    MatrixXcd alpha = QubitArray_getQubit_i_alpha(qa,i);
    printInlineMatrixXcd("alpha",alpha);

    MatrixXcd beta = QubitArray_getQubit_i_beta(qa,i);
    printInlineMatrixXcd("beta",beta);
}

void QubitArray_reportIntmap(QubitArray* qa){
    if (qa->intmap == NULL){
        return;
    }
    int nqubit = QubitArray_getNqubit(qa);
    for (int i=0; i<nqubit; i++){
        for (int j=i; j<nqubit; j++){
            char key[20] = "\0";
            snprintf(key,20,"%s[%d][%d]","intmap",i,j);
            printInlineMatrixXcd(key,QubitArray_getIntmap_i_j(qa,i,j));
        }
    }
}

void QubitArray_reportPsiaPsib(QubitArray* qa){
    printInlineMatrixXcd("psia",QubitArray_getPsia(qa));
    printInlineMatrixXcd("psib",QubitArray_getPsib(qa));
}

void QubitArray_reportPsi0(QubitArray* qa){
    printInlineMatrixXcd("psi0",QubitArray_getPsi0(qa));
}


void QubitArray_report(QubitArray* qa){

    printLineSection();
    printTitle("Structure QubitArray");
    
    printSubTitle("General properties");
    printStructElementBool("overhaus",QubitArray_getOverhaus(qa));
    
    printSubTitle("Qubit properties");
    int nqubit = QubitArray_getNqubit(qa);
    printStructElementInt("nqubit (#)",QubitArray_getNqubit(qa));
    printLine();

    for (int i=0; i<nqubit; i++){
        printf("\n");
        QubitArray_reportQubit_i(qa,i);
    }
    printf("\n");
    printLine();
    
    printSubTitle("Interactions (kHz)");
    QubitArray_reportIntmap(qa);
    printf("\n");

    printSubTitle(" Qubit total alpha, beta");
    QubitArray_reportPsiaPsib(qa);
    printf("\n");

    printSubTitle("Qubit initial state");
    QubitArray_reportPsi0(qa);
    printf("\n");

    printSubTitle("Qubit alpha,beta index");
    if (QubitArray_get_alphaidx(qa) != NULL){
        printStructElementInt("alphaidx",*(QubitArray_get_alphaidx(qa)));
        printStructElementInt("betaidx",*(QubitArray_get_betaidx(qa)));
    }else{ 
        printStructElementChar("alphaidx","NULL");
        printStructElementChar("betaidx","NULL");
    }
    printf("\n");
    printLineSection();
}

