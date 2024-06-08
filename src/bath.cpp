#include "../include/bath.h"
#include "../include/memory.h"
#include "../include/hamiltonian.h"
#include <float.h> // FLT_EPSILON

/* High level functions --------------------------------------------*/

void BathArray_setBathHypfs(BathArray* ba, QubitArray* qa){
    int nqubit = QubitArray_getNqubit(qa);
    int nspin = BathArray_getNspin(ba);
    
    for (int iqubit=0; iqubit<nqubit; iqubit++){

        double qgyro = QubitArray_getQubit_i_gyro(qa,iqubit);
        double* qxyz = QubitArray_getQubit_i_xyz(qa,iqubit);

        for (int ispin=0; ispin<nspin; ispin++){

            double spgyro = BathArray_getBath_i_gyro(ba,ispin);
            double* spxyz = BathArray_getBath_i_xyz(ba,ispin);
            MatrixXcd A = BathArray_getBath_i_hypf_j(ba,ispin,iqubit);
            if (A.rows() == 0 && A.cols() == 0){
                MatrixXcd Adip = calPointDipoleTensor(qxyz, spxyz, qgyro, spgyro); // radkHz            
                BathArray_setBath_i_hypf_j(ba, Adip, ispin, iqubit);    
            }
        }
    }
}

void BathArray_setBathStatesRandom(BathArray* ba){
    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        float S = BathArray_getBath_i_spin(ba,i);
        int dim = (int)(2*S+1);
        int r = rand() % dim; 
        float ms = S - r;
        BathArray_setBath_i_state(ba,ms,i);
    }
}

double BathArray_getOverhaus(BathArray* ba, int iq){

    int nspin = BathArray_getNspin(ba);
    double overhaus = 0.0;

    for (int ib=0; ib<nspin; ib++){
        overhaus += BathArray_getBath_i_overhaus_j(ba,ib,iq);
    }

    return overhaus;
}

double BathArray_getBath_i_overhaus_j(BathArray* ba, int isp, int iq){
    // 𝐻_(𝑜ℎ,𝑖𝑎)^𝒥=⟨𝑚_𝑎^𝒥 | 𝑆_(𝑧,𝑖) 𝐴_(𝑧𝑧,𝑖𝑎) 𝐼_(𝑧,𝑎) |𝑚_𝑎^𝒥 ⟩=𝑚_𝑎^𝒥 𝐴_(𝑧𝑧,𝑖𝑎) 𝑆_(𝑧,𝑖)
    // Here we'll get overhaus = 𝑚_𝑎^𝒥 𝐴_(𝑧𝑧,𝑖𝑎)
    float state = BathArray_getBath_i_state(ba,isp);
    MatrixXcd hypf = BathArray_getBath_i_hypf_j(ba,isp,iq);
    return hypf(2,2).real() * state;
}

void BathArray_setBathDisorders(BathArray* ba){
    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        double disorder = 0.0;
        for (int j=0; j<nspin; j++){
            if (i!=j){
                disorder += BathArray_getBath_i_disorder_j(ba,i,j);
            }
        }
        BathArray_setBath_i_disorder(ba,disorder,i);
    }
}

double BathArray_getBath_i_disorder_j(BathArray* ba, int isp, int jsp){

    // disorder for ith spin : 
    // 𝐻_(𝑚𝑓,𝑘𝑎)^𝒥 = ⟨𝑚_𝑎^𝒥 | 𝐼_(𝑧,𝑘) 𝐽_(𝑧𝑧,𝑘𝑎) 𝐼_(𝑧,𝑎) |𝑚_𝑎^𝒥 ⟩ = 𝑚_𝑎^𝒥 𝐽_(𝑧𝑧,𝑘𝑎) 𝐼_(𝑧,𝑘)
    // disorder for spin k = 𝑚_𝑎^𝒥 𝐽_(𝑧𝑧,𝑘𝑎)

    // 𝑚_𝑎^𝒥
    float mj = BathArray_getBath_i_state(ba,jsp); 

    if (isp < jsp){
        // 𝐽_(𝑧𝑧,𝑘𝑎)
        double Jzz = BathArray_int_i_j(ba,isp,jsp)(2,2).real(); 
        // disorder for spin k = 𝑚_𝑎^𝒥 𝐽_(𝑧𝑧,𝑘𝑎)
        return Jzz * mj;

    }else if (isp == jsp){
        return 0.0;

    }else{ // (isp > jsp)
        double Jzz = BathArray_int_i_j(ba,jsp,isp)(2,2).real();
        return Jzz * mj;
    }
}

// spin pairs
void BathArray_connectivity(int*** cmap, float*** stmap, BathArray* ba, float rdip, float rdipcut){
     
    // cmap [nspin][nspin] = 1 if connected else 0
    // stmap [nspin][nspin] = strength of the connection

    int nspin = BathArray_getNspin(ba);
    
    // Connectivity Map and Strength Map
    *cmap = allocInt2d(nspin,nspin);
    *stmap = allocFloat2d(nspin,nspin);

    // Find connectivity Map and strength Map
    double r = 0.0;
    float strength = 0.0;
    MatrixXcd tensor;

    for (int sp1=0; sp1<nspin; sp1++){
        for (int sp2=sp1+1; sp2<nspin; sp2++){
            
            r      = BathArray_dist_i_j(ba,sp1,sp2);
            tensor = BathArray_int_i_j(ba,sp1,sp2);
            strength = tensor(2,2).real();

            if ( r < rdip && r > rdipcut){
                (*cmap)[sp1][sp2] = 1;
                (*cmap)[sp2][sp1] = 1;
            }

            (*stmap)[sp1][sp2] = fabs(strength);
            (*stmap)[sp2][sp1] = fabs(strength);
        }
    }
}

void makeSparsemap(int*** spmap, int** cmap, int nspin){

    // spmap[nspin][ncol]
    // ncol = connected spin number + 1
    //   [i][0] : connected spin #
    //   [i][j] : connected spin index
    *spmap = (int**)calloc(nspin,sizeof(int*));

    for (int row=0; row<nspin; row++){
        (*spmap)[row] = (int*)calloc(1,sizeof(int));
        (*spmap)[row][0] = 1;
    }

    for (int row=0; row<nspin; row++){
        int colLength = 1;
        for (int col=0; col<nspin; col++){
            reallocInt1d(&((*spmap)[row]),colLength+1);
            if (cmap[row][col]!=0){
                (*spmap)[row][colLength] = col;
                colLength++;
            }
        }
        (*spmap)[row][0] = colLength;
    }
}

int BathArray_dimBath_i(BathArray* ba, int i){

    if (ba==NULL){
        fprintf(stderr,"Error(BathArray_dimBath_i): ba is NULL\n");
        exit(EXIT_FAILURE);
    }

    float s = BathArray_getBath_i_spin(ba,i);
    return int((2.0*s)+1.0);
}

int  BathSpin_dim(BathSpin* bs){

    if (bs==NULL){
        return 0;
    }

    float s = bs->spin;
    return int((2.0*s)+1.0);
}

int BathArray_dim(BathArray* ba){

    if (ba==NULL){
        return 1;
    }

    int nspin = BathArray_getNspin(ba);
    int dim = 1;
    for (int i=0; i<nspin; i++){
        dim *= BathArray_dimBath_i(ba,i);
    }
    return dim;
}

// spin interaction tensor between two bath spins
MatrixXcd BathArray_int_i_j(BathArray* ba, int i, int j){

    if (i>=j){
        fprintf(stderr,"Error(BathArray_Int_i_j): i and j should be i<j\n");
        exit(EXIT_FAILURE);
    }

    ///////////////////////////////////////////////////////////
    // Check if i-th, j-th spin are sharing the same mainspin
    ///////////////////////////////////////////////////////////
    int ib_mainspidx = BathArray_getBath_i_mainspidx(ba,i);
    int jb_mainspidx = BathArray_getBath_i_mainspidx(ba,j);
    
    if (ib_mainspidx==jb_mainspidx && (ib_mainspidx != -1 && jb_mainspidx != -1)){
        // Check mainspin or subspin
        bool ib_is_subspin = false;
        bool jb_is_subspin = false;
        char* ib_name = BathArray_getBath_i_name(ba,i);
        char* jb_name = BathArray_getBath_i_name(ba,j);


        // if subspin, then the name would be "main_type"
        char ib_name_[MAX_CHARARRAY_LENGTH] = "";
        char jb_name_[MAX_CHARARRAY_LENGTH] = "";
        sprintf(ib_name_,"%s_",ib_name);
        sprintf(jb_name_,"%s_",jb_name);

        // Check if ib is subspin
        if (strstr(ib_name,jb_name_) != NULL){
            ib_is_subspin = true;
        }

        // Check if jb is subspin
        if (strstr(jb_name,ib_name_) != NULL){
            jb_is_subspin = true;
        }
        
        if (ib_is_subspin && jb_is_subspin){
            ; // sub - sub > point-dipole

        }else if (ib_is_subspin && !jb_is_subspin){

            // sub - main

            // Check summation of hypf_sub
            MatrixXcd ib_hypf_sub = BathArray_getBath_i_hypf_sub(ba,i);

            // if the summation of hypf_sub is not zero, then point-dipole
            if (!ib_hypf_sub.isZero(FLT_EPSILON)){
                return ib_hypf_sub;
            }
            
        }else if (!ib_is_subspin && jb_is_subspin){

            // main - sub

            // Check summation of hypf_sub
            MatrixXcd jb_hypf_sub = BathArray_getBath_i_hypf_sub(ba,j);

            // if the summation of hypf_sub is not zero, then point-dipole
            if (!jb_hypf_sub.isZero(FLT_EPSILON)){
                return jb_hypf_sub;
            }

        }else{
            // main - main
            fprintf(stderr,"Error(BathArray_int_i_j): i=%d and j=%d are both main spins\n",i,j);
            fprintf(stderr,"Error(BathArray_int_i_j): They should have different mainspidx\n");
            exit(EXIT_FAILURE);
        }
    }
    ///////////////////////////////////////////////////////////
    
    // Both are the main spins
    double* xyz1 = BathArray_getBath_i_xyz(ba,i);
    double* xyz2 = BathArray_getBath_i_xyz(ba,j);
    double  gyro1 = BathArray_getBath_i_gyro(ba,i);
    double  gyro2 = BathArray_getBath_i_gyro(ba,j);

    MatrixXcd pdtensor = calPointDipoleTensor(xyz1,xyz2,gyro1,gyro2);
    return pdtensor;    
}

double BathArray_dist_i_j(BathArray* ba, int i, int j){
    double* xyz1 = BathArray_getBath_i_xyz(ba,i);
    double* xyz2 = BathArray_getBath_i_xyz(ba,j);
    return dist(xyz1,xyz2);

}

// Hamiltonian
MatrixXcd BathArray_TotalHamil(BathArray* ba, MatrixXcd** sigmas, float* bfield){

    int nspin = BathArray_getNspin(ba);
    int dimBA = BathArray_dim(ba);

    // Bath Hamiltonian
    MatrixXcd HbTotal = MatrixXcd::Zero(dimBA,dimBA);

    // Single Hamiltonian
    for (int ib=0; ib<nspin; ib++){
        MatrixXcd Hbi = BathArray_SingleHamil(ba,sigmas,ib,bfield);
        MatrixXcd HbiExpand = expandHamiltonian(sigmas, Hbi, nspin, ib);
        HbTotal += HbiExpand;
    }

    // Interaction Hamiltonian
    for (int i=0; i<nspin; i++){
        MatrixXcd Hbij;
        for (int j=i+1; j<nspin; j++){
            Hbij = BathArray_InteractionHamil(ba,sigmas,i,j);
            HbTotal += Hbij;
        }
    }

    return HbTotal;
}

MatrixXcd BathArray_SingleHamil(BathArray* ba, MatrixXcd** sigmas, int ib, float* bfield){

    // Single Hamiltonian
    MatrixXcd Hbi_Zeeman = BathArray_ZeemanHamil(ba,sigmas,ib,bfield);
    MatrixXcd Hbi_Detuning = BathArray_DetuningHamil(ba,sigmas,ib);
    MatrixXcd Hbi_Disorder = BathArray_DisorderHamil(ba,sigmas,ib);
    MatrixXcd Hbi_Quad = BathArray_QuadHamil(ba,sigmas,ib);

    return Hbi_Zeeman + Hbi_Detuning + Hbi_Disorder + Hbi_Quad;
}

MatrixXcd BathArray_ZeemanHamil(BathArray* ba, MatrixXcd** sigmas, int ib, float* bfield){

    int nspin = BathArray_getNspin(ba);

    if (ib >= nspin){
        fprintf(stderr,"Error: BathArray_ZeemanHamil: ib = %d is out of range\n",ib);
        exit(EXIT_FAILURE);
    }

    double spgyro = BathArray_getBath_i_gyro(ba,ib);
    MatrixXcd vecZeeman = calZeemanVector(spgyro,bfield);

    return calHamiltonianSingleInt(vecZeeman,sigmas[ib]);
}

MatrixXcd BathArray_DetuningHamil(BathArray* ba, MatrixXcd** sigmas, int ib){

    int nspin = BathArray_getNspin(ba);

    if (ib >= nspin){
        fprintf(stderr,"Error: BathArray_DetuningHamil: ib = %d is out of range\n",ib);
        exit(EXIT_FAILURE);
    }

    double detuning = BathArray_getBath_i_detuning(ba,ib);
    MatrixXcd vecDetuning = calDetuningVector(detuning);
    return calHamiltonianSingleInt(vecDetuning,sigmas[ib]);
}

MatrixXcd BathArray_DisorderHamil(BathArray* ba, MatrixXcd** sigmas, int ib){

    int nspin = BathArray_getNspin(ba);

    if (ib >= nspin){
        fprintf(stderr,"Error: BathArray_DisorderHamil: ib = %d is out of range\n",ib);
        exit(EXIT_FAILURE);
    }

    double disorder = BathArray_getBath_i_disorder(ba,ib);
    MatrixXcd vecDisorder = calDetuningVector(disorder);

    return calHamiltonianSingleInt(vecDisorder,sigmas[ib]);

}

MatrixXcd BathArray_QuadHamil(BathArray* ba, MatrixXcd** sigmas, int ib){

    int nspin = BathArray_getNspin(ba);

    if (ib >= nspin){
        fprintf(stderr,"Error: BathArray_QuadHamil: ib = %d is out of range\n",ib);
        exit(EXIT_FAILURE);
    }

    MatrixXcd quad = BathArray_getBath_i_quad(ba,ib);

    return calHamiltonianSelfInt(sigmas[ib],quad);
}

MatrixXcd BathArray_InteractionHamil(BathArray* ba, MatrixXcd** sigmas, int ib, int jb){

    int nspin = BathArray_getNspin(ba);
    int dimBA = BathArray_dim(ba);

    // Bath Hamiltonian
    MatrixXcd Hbij = MatrixXcd::Zero(dimBA,dimBA);

    if (ib<jb){
        MatrixXcd tensor = BathArray_int_i_j(ba,ib,jb);
        bool isEmpty = tensor.isZero(FLT_EPSILON);

        if (!isEmpty){
            Hbij = calHamiltonianHeteroInt(sigmas, tensor, nspin, ib, jb);
        }
    }else{
        fprintf(stderr,"Error: BathArray_InteractionHamil: ib,jb = %d,%d is out of range or ib>jb\n",ib,jb);
        exit(EXIT_FAILURE);
    }

    return Hbij;

}

MatrixXcd** BathArray_PauliOperators(BathArray* ba){

    int nspin = BathArray_getNspin(ba);
    int dimBA = BathArray_dim(ba);

    // Pauli operators
    MatrixXcd** sigmas = new MatrixXcd*[nspin];

    for (int ib=0; ib<nspin; ib++){
        int dimB = BathArray_dimBath_i(ba,ib);
        sigmas[ib] = getPauliOperators(dimB);
    }

    return sigmas;
}

// density matrix
MatrixXcd BathArray_Rho0(BathArray* ba, bool isEnsemble){

    int dimBA = BathArray_dim(ba);

    if (isEnsemble || dimBA == 1){
        return MatrixXcd::Identity(dimBA,dimBA) * (1.0/(double)dimBA);
    }else{
        MatrixXcd psi0 = BathArray_Psi0(ba);
        return psi0 * psi0.adjoint();
    }

}

MatrixXcd BathArray_Psi0(BathArray* ba){

    int nspin = BathArray_getNspin(ba);
    int dimBA = BathArray_dim(ba);

    MatrixXcd psi0;
    for (int ib=0; ib<nspin; ib++){
        float S = BathArray_getBath_i_spin(ba,ib);
        float ms = BathArray_getBath_i_state(ba,ib);

        if (ib==0){
            psi0 = getSpinor(S,ms);
        }else{
            psi0 = kron(psi0,getSpinor(S,ms));
        }        
    }

    return psi0;
}


/* Low level functions --------------------------------------------*/

BathArray*  BathArray_init(){
    BathArray* ba = (BathArray*)allocArray1d(1,sizeof(BathArray));
    BathArray_setNspin(ba,0);
    BathArray_setProp_nspecies(ba,0);
    
    ba->bath = NULL;
    ba->prop_names = NULL;
    ba->prop_gyros = NULL;
    ba->prop_spins = NULL;
    return ba;
}

void BathArray_report(BathArray* ba){
    BathArray_reportSpinProperties(ba);
    BathArray_reportBath(ba);
}

void BathArray_reportSpinProperties(BathArray* ba){
    printStructElementInt("nspecies",BathArray_getProp_nspecies(ba));
    printStructElementChar2d("names",BathArray_getProp_names(ba),BathArray_getProp_nspecies(ba));
    printStructElementDouble1d("gyros",BathArray_getProp_gyros(ba),BathArray_getProp_nspecies(ba));
    printStructElementFloat1d("spins",BathArray_getProp_spins(ba),BathArray_getProp_nspecies(ba));
}

void BathArray_reportBath(BathArray* ba){

    int nspin = BathArray_getNspin(ba);
    printStructElementInt("nspin (#)",nspin);
    printLine();

    for (int i=0; i<nspin; i++){
        
        if (verbosity || (i<3 || i>nspin-3)){ 
            BathArray_reportBath_i_props(ba,i);
        }

        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
    printf("\n");
    printLine();
}


void BathArray_reportBath_i_props(BathArray* ba, int i){
    char* name = BathArray_getBath_i_name(ba,i);
    float spin = BathArray_getBath_i_spin(ba,i);
    double gyro = BathArray_getBath_i_gyro(ba,i);
    double* xyz = BathArray_getBath_i_xyz(ba,i);    
    printf("      [%3d] %5s %7.3lf %7.3lf %7.3lf   ( S = %2.1f, gyro = %7.3lf )\n",i,name,xyz[0],xyz[1],xyz[2],spin,gyro);   
}

void BathArray_reportBath_states(BathArray* ba){

    int nspin = BathArray_getNspin(ba);
    printStructElementInt("nbathspin (#)",nspin);
    printLine();

    for (int i=0; i<nspin; i++){
        if (verbosity || (i<3 || i>nspin-3)){ 
            float state = BathArray_getBath_i_state(ba,i);

            char message[100];
            sprintf(message,"bath[%d].state ",i);
            printStructElementFloat(message,state);
        }
        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
    printLine();
    printf("\n");
}

void BathArray_reportBath_detunings(BathArray* ba){

    printSubTitle("Bath Detunings");
    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        if (verbosity || (i<3 || i>nspin-3)){ 
            double detuning = BathArray_getBath_i_detuning(ba,i);
            printf("      [%3d] detuning = %7.3g\n",i,detuning);
        }
        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
    printf("\n");
}

void BathArray_reportBath_disorders(BathArray* ba){

    printSubTitle("Bath Disorders (radkHz)");
    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        if (verbosity || (i<3 || i>nspin-3)){ 
            char message [100];
            snprintf(message,100,"Bath[%d].disorder",i);
            printStructElementDouble(message,BathArray_getBath_i_disorder(ba,i));
        }
        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
    printf("\n");
}

void BathArray_reportBath_hypf(BathArray* ba, int nqubit){

    printSubTitle("Bath Hyperfine tensors (radkHz) ");

    printf("      hypf[ib][iq] (ib : bath spin index, iq : qubit index)\n");

    int nspin = BathArray_getNspin(ba);

    for (int iq=0; iq<nqubit; iq++){
        for (int i=0; i<nspin; i++){
            if (verbosity || (i<3 || i>nspin-3)){              
                char key[100];
                sprintf(key,"hypf[%d][%d]",i,iq);
                MatrixXcd hypf = BathArray_getBath_i_hypf_j(ba,i,iq);
                printInlineMatrixXcd(key,hypf);
            }
            if (!verbosity && i==3){
                printf("         :\n\n");
            }
        }
        printf("\n");
    }

    printLineSection();
}

void BathArray_reportBath_quad(BathArray* ba){

    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        if (verbosity || (i<3 || i>nspin-3)){ 
            char key[100];
            snprintf(key,100,"quad[%d]",i);
            MatrixXcd quad = BathArray_getBath_i_quad(ba,i);
            printInlineMatrixXcd(key,quad);
        }
        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
}

void BathArray_reportBath_hypf_sub(BathArray* ba){

    int nspin = BathArray_getNspin(ba);
    for (int i=0; i<nspin; i++){
        if (verbosity || (i<3 || i>nspin-3)){ 
            MatrixXcd hypf_sub = BathArray_getBath_i_hypf_sub(ba,i);
            int mainspidx = BathArray_getBath_i_mainspidx(ba,i);

            char key[100];
            snprintf(key,100,"hypf_sub[%d] (mainsp = %d)",i, mainspidx);
            printInlineMatrixXcd(key,hypf_sub);
        }
        if (!verbosity && i==3){
            printf("         :\n");
        }
    }
}

//alloc

void BathArray_allocBath(BathArray* ba, int nqubit){
    int nspin = BathArray_getNspin(ba);
    ba->bath = (BathSpin**)allocArray2d(nspin,1,sizeof(BathSpin));
    for (int i=0; i<nspin; i++){
        BathArray_allocBath_i_hypf(ba,i,nqubit);
        BathArray_setBath_i_mainspidx(ba,-1,i);

        BathArray_setBath_i_name(ba,"",i);
        BathArray_setBath_i_spin(ba,0.0,i);
        BathArray_setBath_i_gyro(ba,0.0,i);
        BathArray_setBath_i_xyz(ba,(double[]){0.0,0.0,0.0},i);
        BathArray_setBath_i_state(ba,0.0,i);
        BathArray_setBath_i_detuning(ba,0.0,i);
        BathArray_setBath_i_disorder(ba,0.0,i);
        ba->bath[i]->quad = MatrixXcd::Zero(3,3);
        //BathArray_setBath_i_quad(ba,MatrixXcd::Zero(3,3),i);
    }
}

void BathArray_reallocBath(BathArray* ba, int old_length, int new_length, int nqubit){
    if (new_length <= old_length){
        fprintf(stderr,"Error(BathArray_reallocBath): it cannot be new_length <= old_length\n");
        exit(EXIT_FAILURE);
    }
    
    reallocArray2d((void***)&(ba->bath),old_length,new_length,1,sizeof(BathSpin));

    for (int i=old_length; i<new_length; i++){
        BathArray_allocBath_i_hypf(ba,i,nqubit);
        BathArray_setBath_i_mainspidx(ba,-1,i);

        BathArray_setBath_i_name(ba,"",i);
        BathArray_setBath_i_spin(ba,0.0,i);
        BathArray_setBath_i_gyro(ba,0.0,i);
        BathArray_setBath_i_xyz(ba,(double[]){0.0,0.0,0.0},i);
        BathArray_setBath_i_state(ba,0.0,i);
        BathArray_setBath_i_detuning(ba,0.0,i);
        BathArray_setBath_i_disorder(ba,0.0,i);
        ba->bath[i]->quad = MatrixXcd::Zero(3,3);
        //BathArray_setBath_i_quad(ba,MatrixXcd::Zero(3,3),i);
    }
}

void BathArray_allocProp(BathArray* ba){
    int nspecies = BathArray_getProp_nspecies(ba);
    ba->prop_names = allocChar2d(nspecies,MAX_CHARARRAY_LENGTH);
    ba->prop_gyros = allocDouble1d(nspecies);
    ba->prop_spins = allocFloat1d(nspecies);
}

void BathArray_reallocProp(BathArray* ba, int nspin_old, int nspin_new){
    if (nspin_new <= nspin_old){
        fprintf(stderr,"Error(BathArray_realloc_Prop): it cannot be nspin_new <= nspin_old\n");
        exit(EXIT_FAILURE);
    }
    reallocChar2d(&(ba->prop_names),nspin_old,nspin_new,MAX_CHARARRAY_LENGTH);
    reallocDouble1d(&(ba->prop_gyros),nspin_new);
    reallocFloat1d(&(ba->prop_spins),nspin_new);
}

void BathArray_allocBath_i_hypf(BathArray* ba, int i, int nqubit){
    ba->bath[i]->hypf = new MatrixXcd[nqubit];
}

//sets
void BathArray_setNspin(BathArray* ba, const int nspin){
    ba->nspin = nspin;
}

void BathArray_setBath_i(BathArray* ba, const BathSpin* bath, int i, int nqubit){

    int nspin = BathArray_getNspin(ba);
    if (bath==NULL || i >= nspin){
        fprintf(stderr,"Error(BathArray_setBath_i): bath couldn't be NULL\n");
        exit(EXIT_FAILURE);   
    }
    BathArray_setBath_i_name(ba,bath->name,i);
    BathArray_setBath_i_spin(ba,bath->spin,i);
    BathArray_setBath_i_gyro(ba,bath->gyro,i);
    BathArray_setBath_i_xyz(ba,bath->xyz,i);
    BathArray_setBath_i_state(ba,bath->state,i);
    BathArray_setBath_i_detuning(ba,bath->detuning,i);
    BathArray_setBath_i_disorder(ba,bath->disorder,i);

    BathArray_setBath_i_hypf_sub(ba,bath->hypf_sub,i);
    BathArray_setBath_i_mainspidx(ba,bath->mainspidx,i);
    ba->bath[i]->quad = bath->quad;
    //BathArray_setBath_i_quad(ba,bath->quad,i);
    for (int j=0; j<nqubit; j++){
        BathArray_setBath_i_hypf_j(ba,bath->hypf[j],i,j);
    }
}

void BathArray_setBath_i_name(BathArray* ba, const char* name, int i){
    strcpy(ba->bath[i]->name,name);
}

void BathArray_setBath_i_spin(BathArray* ba, const float spin, int i){
    ba->bath[i]->spin = spin;
}

void BathArray_setBath_i_gyro(BathArray* ba, const double gyro, int i){
    ba->bath[i]->gyro = gyro;
}

void BathArray_setBath_i_xyz(BathArray* ba, const double* xyz, int i){
    ba->bath[i]->xyz[0] = xyz[0];
    ba->bath[i]->xyz[1] = xyz[1];
    ba->bath[i]->xyz[2] = xyz[2];
}

void BathArray_setBath_i_state(BathArray* ba, const float state, int i){

    // float S = BathArray_getBath_i_spin(ba,i);

    // if (!isSubLevel(S,state)){
    //     fprintf(stderr,"Error(BathArray_setBath_i_state): S = %2.1f\n",S);
    //     fprintf(stderr,"Error(BathArray_setBath_i_state): state = %2.1f is out of range\n",state);
    //     exit(EXIT_FAILURE);
    // }

    ba->bath[i]->state = state;

}

void BathArray_setBath_i_detuning(BathArray* ba, const double detuning, int i){
    ba->bath[i]->detuning = detuning;
}

void BathArray_setBath_i_disorder(BathArray* ba, const double disorder, int i){
    ba->bath[i]->disorder = disorder;
}

void BathArray_setBath_i_hypf_j(BathArray* ba, const MatrixXcd hypf, int i, int j){
    ba->bath[i]->hypf[j] = hypf;
} 
 
void BathArray_setBath_i_quad(BathArray* ba, const MatrixXcd quad, int i){

    float S = BathArray_getBath_i_spin(ba,i);
    if (S<1.0){
        if (rank==0){
            printf("Warning(BathArray_setBath_i_quad): BathSpin[%d] S = %2.1f\n",i,S);
            printf("Warning(BathArray_setBath_i_quad): You set the quadrupole but the spin number is larger than 0.5\n");
            printf("Warning(BathArray_setBath_i_quad): I hope you know what you do\n");
        }
    }
    ba->bath[i]->quad = quad;
}

void BathArray_setBath_i_hypf_sub(BathArray* ba, const MatrixXcd hypf_sub, int i){
    ba->bath[i]->hypf_sub = hypf_sub;
}

void BathArray_setBath_i_mainspidx(BathArray* ba, const int mainspidx, int i){
    ba->bath[i]->mainspidx = mainspidx;

}

void BathArray_setProp_nspecies(BathArray* ba, const int nspecies){
    ba->prop_nspecies = nspecies;
}

void BathArray_setProp_names_i(BathArray* ba, const char* name, const int i){
    strcpy(ba->prop_names[i],name);
}

void BathArray_setProp_gyros_i(BathArray* ba, const double gyro, const int i){
    ba->prop_gyros[i] = gyro;
}

void BathArray_setProp_spins_i(BathArray* ba, const float spin, const int i){
    ba->prop_spins[i] = spin;
}

// get
int BathArray_getNspin(BathArray* ba){
    if (ba == NULL){
        return 0;
    }
    return ba->nspin;
}

int BathArray_getProp_nspecies(BathArray* ba){
    return ba->prop_nspecies;
}

char** BathArray_getProp_names(BathArray* ba){
    return ba->prop_names;
}

double* BathArray_getProp_gyros(BathArray* ba){
    return ba->prop_gyros;
}

float* BathArray_getProp_spins(BathArray* ba){
    return ba->prop_spins;
}

BathSpin* BathArray_getBath_i(BathArray* ba, int i){
    return ba->bath[i];
}

char* BathArray_getBath_i_name(BathArray* ba, int i){
    return ba->bath[i]->name;
}

float BathArray_getBath_i_spin(BathArray* ba, int i){
    return ba->bath[i]->spin;
}

double BathArray_getBath_i_gyro(BathArray* ba, int i){
    return ba->bath[i]->gyro;
}

double* BathArray_getBath_i_xyz(BathArray* ba, int i){
    return ba->bath[i]->xyz;
}

float BathArray_getBath_i_state(BathArray* ba, int i){
    return ba->bath[i]->state;
}

double BathArray_getBath_i_detuning(BathArray* ba, int i){
    return ba->bath[i]->detuning;
}

double BathArray_getBath_i_disorder(BathArray* ba, int i){
    return ba->bath[i]->disorder;
}

MatrixXcd BathArray_getBath_i_hypf_j(BathArray* ba, int i, int j){
    return ba->bath[i]->hypf[j];
}

MatrixXcd BathArray_getBath_i_quad(BathArray* ba, int i){
    return ba->bath[i]->quad;
}

MatrixXcd BathArray_getBath_i_hypf_sub(BathArray* ba, int i){
    return ba->bath[i]->hypf_sub;
}

int BathArray_getBath_i_mainspidx(BathArray* ba, int i){
    return ba->bath[i]->mainspidx;
}

// free
void BathArray_freeAll(BathArray* ba){
    BathArray_freeProp_names(ba);
    BathArray_freeProp_gyros(ba);
    BathArray_freeProp_spins(ba);
    BathArray_freeBath(ba);
    freeArray1d((void**)&ba);
}

void BathArray_freeProp_names(BathArray* ba){
    freeChar2d(&(ba->prop_names),ba->prop_nspecies);
}

void BathArray_freeProp_gyros(BathArray* ba){
    freeDouble1d(&(ba->prop_gyros));
}

void BathArray_freeProp_spins(BathArray* ba){
    freeFloat1d(&(ba->prop_spins));
}

void BathArray_freeBath(BathArray* ba){
    for (int i=0; i<ba->nspin; i++){
        BathArray_freeBath_i_hypf(ba,i);
        freeArray1d((void**)&(ba->bath[i]));
    }
}

void BathArray_freeBath_i_hypf(BathArray* ba, int i){
    delete[] ba->bath[i]->hypf;
}

////////////////////////////////////////////////////////////////
// BathSpin

void 
BathSpin_setName(BathSpin* bs, char* name){
    strcpy(bs->name,name);
}

void       
BathSpin_setName_withType(BathSpin* bs, char* name, char* type){
    char name_type[MAX_CHARARRAY_LENGTH];
    snprintf(name_type,MAX_CHARARRAY_LENGTH,"%s_%s",name,type);
    BathSpin_setName(bs,name_type);
}

void 
BathSpin_setSpin(BathSpin* bs, float spin){
    bs->spin = spin;
}

void 
BathSpin_setGyro(BathSpin* bs, double gyro){
    bs->gyro = gyro;
}

void 
BathSpin_setXyz(BathSpin* bs, double* xyz){
    bs->xyz[0] = xyz[0];
    bs->xyz[1] = xyz[1];
    bs->xyz[2] = xyz[2];
}

void
BathSpin_setXyz_fromRxyz(BathSpin* bs, double* xyz0, double* rxyz){
    bs->xyz[0] = xyz0[0] + rxyz[0];
    bs->xyz[1] = xyz0[1] + rxyz[1];
    bs->xyz[2] = xyz0[2] + rxyz[2];
}

void 
BathSpin_setState(BathSpin* bs, float state){

    if (!isSubLevel(bs->spin,state)){
        fprintf(stderr,"Error(BathSpin_setState): S = %2.1f\n",bs->spin);
        fprintf(stderr,"Error(BathSpin_setState): state = %2.1f is out of range\n",state);
        exit(EXIT_FAILURE);
    }
    bs->state = state;
}

void 
BathSpin_setDetuning(BathSpin* bs, double detuning){
    bs->detuning = detuning;
}

void 
BathSpin_setDisorder(BathSpin* bs, double disorder){
    bs->disorder = disorder;
}

void 
BathSpin_setHypf_i(BathSpin* bs, MatrixXcd hypf, int iq){
    bs->hypf[iq] = hypf;
}

void 
BathSpin_setQuad(BathSpin* bs, MatrixXcd quad){
    bs->quad = quad;
}

void 
BathSpin_setQuad_fromEFG(BathSpin* bs, MatrixXcd efg, double eq, float spin){
    // efg : Hartree/Bohr^2
    // eq : 10e-30 m^2

    // plank constant
    double h = 4.135667*pow(10,-15); // eV * sec

    // Unit conversion factors
    double Hatress2eV = 27.211386; // 1Hatree = 27.211386eV
    double BohrRadius2m = 0.5291772*pow(10,-10); // 1Bohr_raidus = 5.29177E-11 m
    double BohrRadiusSq2mSq = pow(BohrRadius2m,2);  // 1Bohr_raidus^2 = (5.29177E-11 m)^2
    double Hz2MHz = pow(10,-6); // 1Hz = 10^-6 MHz
    double UnitConversion = Hatress2eV/BohrRadiusSq2mSq*Hz2MHz;

    // Conversion
    efg = UnitConversion * efg; // Hartree/Bohr^2 -> MHz
    eq = eq * 1.0e-30; // 10e-30 m^2 -> m^2

    // Quadrupole (radkHz)
    bs->quad = MHZ_TO_RADKHZ((eq)/(spin*(2.0*spin-1.0))/h * efg);
}


void 
BathSpin_setHypfSub(BathSpin* bs, MatrixXcd hypf_sub){
    bs->hypf_sub = hypf_sub;
}

void 
BathSpin_setMainspidx(BathSpin* bs, int mainspidx){
    bs->mainspidx = mainspidx;
}


char*      
BathSpin_getName(BathSpin* bs){
    return bs->name;
}

float      
BathSpin_getSpin(BathSpin* bs){
    return bs->spin;
}

double     
BathSpin_getGyro(BathSpin* bs){
    return bs->gyro;
}

double*    
BathSpin_getXyz(BathSpin* bs){
    return bs->xyz;
}

float      
BathSpin_getState(BathSpin* bs){
    return bs->state;
}

double     
BathSpin_getDetuning(BathSpin* bs){
    return bs->detuning;
}

double     
BathSpin_getDisorder(BathSpin* bs){
    return bs->disorder;
}

MatrixXcd  
BathSpin_getHypf_i(BathSpin* bs, int iq){
    return bs->hypf[iq];
}

MatrixXcd  
BathSpin_getQuad(BathSpin* bs){
    return bs->quad;
}

MatrixXcd  
BathSpin_getHypfSub(BathSpin* bs){
    return bs->hypf_sub;
}

int        
BathSpin_getMainspidx(BathSpin* bs){
    return bs->mainspidx;
}
