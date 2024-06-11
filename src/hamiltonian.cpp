#include "../include/hamiltonian.h"
#include "../include/memory.h"
#include <iostream>

/*! @brief point-dipole interaction tensor
 * @param[in] xyz1 position of spin 1
 * @param[in] xyz2 position of spin 2
 * @param[in] gamma1 gyromagnetic ratio of spin 1
 * @param[in] gamma2 gyromagnetic ratio of spin 2
 * @return point-dipole interaction tensor
 * @details The point-dipole interaction tensor is given by
 * Hdd = (1/r^3) * (mu_0 / (4 * pi * hbar) ) * (gamma1 * gamma2 * mu_B^2) * (S1 dot S2 - 3(S1 dot r)(S2 dot r)/r^2
 * unit : (1/r^3) * (mu_0 / (4 * pi * hbar) ) * (gamma1 * gamma2 * mu_B^2)
 * 1/r^3 : A^-3
 * mu_0 : N/A^2
 * hbar : J*s
 * gamma1 : radkHz/G -> radkHz/T
 * gamma2 : radkHz/G -> radkHz/T
 * mu_B : J/T
 * 
 */


// Zeeman Interaction [rad kHz]
MatrixXcd ZeemanInt(double gamma, float* B0, MatrixXcd* Pmat, int m){
    MatrixXcd H(m,m);
    H = MatrixXcd::Zero(m,m);
    H = -1.0 * gamma * B0[0] * Pmat[1] + 
        -1.0 * gamma * B0[1] * Pmat[2] + 
        -1.0 * gamma * B0[2] * Pmat[3];
    // cf. if B0 = (0,0,Bz), H = -gamma*Bz*Z
    return H;
}

// Detuning Interaction [rad kHz]
MatrixXcd DetuningInt(double detun, MatrixXcd* Pmat, int m){
    return detun * Pmat[3];
}

/** 
 * @brief Zeeman Interaction vector[rad kHz]
 * @details The Zeeman Interaction vector is given by
 *          mu_B (B * g) 
 *          g' = eye(3) * -gamma
 *          Zeeman vector = [ Bxgx , Bygy , Bzgz ]
 * @param[in] gamma gyromagnetic ratio [radkHz/G]
 * @param[in] B0 magnetic field vector [G]
 * @return Zeeman Interaction vector [rad kHz]
*/
MatrixXcd calZeemanVector(double gamma, float* B0){
    MatrixXcd vector = MatrixXcd::Zero(1,3);
    for (int i=0; i<3; i++){
        vector(0,i) = -1.0 * gamma * B0[i];
    }
    return vector;
}

/** 
 * @brief Detuning Interaction vector[rad kHz]
 * @details The Detuning Interaction vector is given by
 *          [0,0,detun]
 * @param[in] detun detuning [rad kHz]
 * @return Detuning Interaction vector [rad kHz]
*/
MatrixXcd calDetuningVector(double detun){
    MatrixXcd vector = MatrixXcd::Zero(1,3);
    vector(0,2) = detun;
    return vector;
}

/** 
 * @brief Overhausler Interaction vector[rad kHz]
 * @details The Overhausler Interaction vector is given by
 *          [0,0,overhaus]
 * @param[in] overhaus overhausler effect [rad kHz]
 * @return Overhausler Interaction vector [rad kHz]
*/
MatrixXcd calOverhauserVector(double overhaus){
    MatrixXcd vector = MatrixXcd::Zero(1,3);
    vector(0,2) = overhaus;
    return vector;
}

// Interaction tensor by asumming Point-dipole approximation [rad kHz]
MatrixXcd calPointDipoleTensor(double xyz1[3], double xyz2[3], double gamma1, double gamma2){

    MatrixXcd Adip = MatrixXcd::Zero(3,3);
    //std::cout << "Adip" << std::endl;
    //std::cout << Adip << std::endl;

    double r = dist(xyz1, xyz2);
    double ct = cosTheta(xyz1,xyz2, r);
    double st = sinTheta(xyz1,xyz2, r);
    double cp = cosPhi(xyz1,xyz2);
    double sp = sinPhi(xyz1,xyz2);
    double DS = H_BAR * gamma1 * gamma2 / pow(r,3);
    // (radkHz/G) * (radkHz/G) / (A^3)
    // hbar * hbar * (radkHz/G) * (radkHz/G) / (A^3)

    //parameter print
    // std::cout << "r : " << r << std::endl;
    // std::cout << "ct : " << ct << std::endl;
    // std::cout << "st : " << st << std::endl;
    // std::cout << "cp : " << cp << std::endl;
    // std::cout << "sp : " << sp << std::endl;
    // std::cout << "DS : " << DS << std::endl;

    //Diag term
    Adip(0,0) = 1 - 3 * pow(st,2) * pow(cp,2);
    Adip(1,1) = 1 - 3 * pow(st,2) * pow(sp,2);
    Adip(2,2) = 1 - 3 * pow(ct,2);

    //xy, yx
    Adip(0,1) = -3 * pow(st,2) * cp * sp;
    Adip(1,0) = -3 * pow(st,2) * cp * sp;

    //xz, zx
    Adip(0,2) = -3 * st * ct * cp;
    Adip(2,0) = -3 * st * ct * cp;

    //yz, zy
    Adip(1,2) = -3 * st * ct * sp;
    Adip(2,1) = -3 * st * ct * sp;

    Adip = DS * Adip;

    // std::cout << "Adip" << std::endl;
    // std::cout << Adip << std::endl;
    
    return Adip; // H = S1 * Adip * S2 
}




/**
 * @brief Hamiltonian for interacting two spins
 * @details The Hamiltonian is given by (S1 * Tensor * S2)
 * @param[in] Pmat1 Spin operator vector for first spin
 * @param[in] Tensor 3x3 interaction tensor
 * @param[in] Pmat2 Spin operator vector for second spin
 * @return Hamiltonian interacting two spins 
 *         (Pmat1[0,1,2] * Tensor * Pmat2[0,1,2])
 * @note Pmat format is the following :
 *       Pmat[0] = Pauli matrix X
 *       Pmat[1] = Pauli matrix Y
 *       Pmat[2] = Pauli matrix Z
 * @note The dimension of Returned Hamiltonian is 
 *       expanded based on two spins
*/
MatrixXcd calHamiltonianHeteroInt(MatrixXcd** Pmats, MatrixXcd Tensor, int nSpin, int iSpin, int jSpin){

    // Calculate the dimension of Hamiltonian
    int dimrow = 1;
    int dimcol = 1;
    for (int i=0; i<nSpin; i++){
        // printf("Pmats[%d][0] : (%d, %d)\n",i,Pmats[i][0].rows(),Pmats[i][0].cols());
        dimrow *= Pmats[i][0].rows();
        dimcol *= Pmats[i][0].cols();
    }
    // printf("dimrow : %d\n",dimrow);
    // printf("dimcol : %d\n",dimcol);
    
    // Calculate expanded Pauli Operators
    MatrixXcd SixSjx, SixSjy, SixSjz;
    MatrixXcd SiySjx, SiySjy, SiySjz;
    MatrixXcd SizSjx, SizSjy, SizSjz;

    for (int i=0; i<nSpin; i++){
        
        MatrixXcd SixSjxTmp, SixSjyTmp, SixSjzTmp;
        MatrixXcd SiySjxTmp, SiySjyTmp, SiySjzTmp;
        MatrixXcd SizSjxTmp, SizSjyTmp, SizSjzTmp;

        if (i == iSpin){
            SixSjxTmp = Pmats[iSpin][1] ; SixSjyTmp = Pmats[iSpin][1] ; SixSjzTmp = Pmats[iSpin][1] ;
            SiySjxTmp = Pmats[iSpin][2] ; SiySjyTmp = Pmats[iSpin][2] ; SiySjzTmp = Pmats[iSpin][2] ;
            SizSjxTmp = Pmats[iSpin][3] ; SizSjyTmp = Pmats[iSpin][3] ; SizSjzTmp = Pmats[iSpin][3] ;
        }
        else if (i == jSpin){
            SixSjxTmp = Pmats[jSpin][1] ; SixSjyTmp = Pmats[jSpin][2] ; SixSjzTmp = Pmats[jSpin][3] ;
            SiySjxTmp = Pmats[jSpin][1] ; SiySjyTmp = Pmats[jSpin][2] ; SiySjzTmp = Pmats[jSpin][3] ;
            SizSjxTmp = Pmats[jSpin][1] ; SizSjyTmp = Pmats[jSpin][2] ; SizSjzTmp = Pmats[jSpin][3] ;            
        }
        else{
            SixSjxTmp = Pmats[i][0]     ; SixSjyTmp = Pmats[i][0]     ; SixSjzTmp = Pmats[i][0]     ;
            SiySjxTmp = Pmats[i][0]     ; SiySjyTmp = Pmats[i][0]     ; SiySjzTmp = Pmats[i][0]     ;
            SizSjxTmp = Pmats[i][0]     ; SizSjyTmp = Pmats[i][0]     ; SizSjzTmp = Pmats[i][0]     ;
        }
        
        if (i==0){
            SixSjx = SixSjxTmp; SixSjy = SixSjyTmp; SixSjz = SixSjzTmp;
            SiySjx = SiySjxTmp; SiySjy = SiySjyTmp; SiySjz = SiySjzTmp;
            SizSjx = SizSjxTmp; SizSjy = SizSjyTmp; SizSjz = SizSjzTmp;
        }
        else{
            SixSjx = kron(SixSjx,SixSjxTmp); SixSjy = kron(SixSjy,SixSjyTmp); SixSjz = kron(SixSjz,SixSjzTmp);
            SiySjx = kron(SiySjx,SiySjxTmp); SiySjy = kron(SiySjy,SiySjyTmp); SiySjz = kron(SiySjz,SiySjzTmp);
            SizSjx = kron(SizSjx,SizSjxTmp); SizSjy = kron(SizSjy,SizSjyTmp); SizSjz = kron(SizSjz,SizSjzTmp);
        }
    }

    if (SixSjx.rows() != dimrow || SixSjx.cols() != dimcol){
        perror("Error(calHamiltonianHeteroInt): Dimension of SixSjx is not matched");
        exit(EXIT_FAILURE);
    }

    // Calculate Hamiltonian
    MatrixXcd H = MatrixXcd::Zero(dimrow,dimcol);
    H = Tensor(0,0) * SixSjx + Tensor(0,1) * SixSjy + Tensor(0,2) * SixSjz +
        Tensor(1,0) * SiySjx + Tensor(1,1) * SiySjy + Tensor(1,2) * SiySjz +
        Tensor(2,0) * SizSjx + Tensor(2,1) * SizSjy + Tensor(2,2) * SizSjz;

    return H;
}

/**
 * @brief Matrix expand to the whole Hilbert space
 * @details The matrix is expanded to the whole Hilbert space
 * @param[in] Pmats Spin operator vector for single spin
 * @param[in] Hi Hamiltonian for single spin
 * @param[in] nSpin The number of spins
 * @param[in] iSpin The index of spin
 * @return Expanded Hamiltonian
*/
MatrixXcd expandHamiltonian(MatrixXcd** Pmats, MatrixXcd Hi, int nSpin, int iSpin){

    // Calculate the dimension of Hamiltonian
    int dimrow = 1;
    int dimcol = 1;
    for (int i=0; i<nSpin; i++){
        dimrow *= Pmats[i][0].rows();
        dimcol *= Pmats[i][0].cols();
    }

    // expand Hamiltonian
    MatrixXcd expandedHamiltonian;

    for (int i=0; i<nSpin; i++){
        MatrixXcd tmpMatrix;
        if (i == iSpin){
            tmpMatrix = Hi;
        }
        else{
            tmpMatrix = Pmats[i][0];
        }

        if (i==0){
            expandedHamiltonian = tmpMatrix;
        }
        else{
            expandedHamiltonian = kron(expandedHamiltonian,tmpMatrix);
        }
    }

    if (expandedHamiltonian.rows() != dimrow || expandedHamiltonian.cols() != dimcol){
        perror("Error(expandHamiltonian): Dimension of expandedHamiltonian is not matched");
        exit(EXIT_FAILURE);
    }

    return expandedHamiltonian;
}

/**
 * @brief Hamiltonian for single spin
 * @details The Hamiltonian is given by (Vector * Pmat)
 * @param[in] Vector 1x3 matrix
 * @param[in] Pmat Pauli matrix for single spin
 * @note Pmat format is the following :
 *       Pmat[0] = Pauli matrix I
 *       Pmat[1] = Pauli matrix X
 *       Pmat[2] = Pauli matrix Y
 *       Pmat[3] = Pauli matrix Z
 * @return Hamiltonian for single spin  (Vector * Pmat)
 * @note The dimension of Returned Hamiltonian is 
 *       the same dimension as the single spin
*/
MatrixXcd calHamiltonianSingleInt(MatrixXcd Vector, MatrixXcd* Pmat1){

    if (Vector.rows() == 1 && Vector.cols() == 3){//1x3
        ;
    }else if (Vector.rows() == 3 && Vector.cols() == 1){ //3x1
        Vector = Vector.transpose().eval();
    }else{
        perror("Error(calHamiltonianSingleInt): Dimension of Vector is not matched");
        exit(EXIT_FAILURE);       
    }

    int dimrow = Pmat1[0].rows();
    int dimcol = Pmat1[0].cols();

    MatrixXcd H = MatrixXcd::Zero(dimrow,dimcol);

    for (int i=0; i<3; i++){
        H += Vector(0,i) * Pmat1[i+1]; //Px, Py, Pz
    }

    return H;
}

/**
 * @brief Hamiltonian for self interaction
 * @details The Hamiltonian is given by (Pmat * Tensor * Pmat)
 * @param[in] Pmat Spin operator vector for single spin
 * @param[in] Tensor 3x3 interaction tensor
 * @note Pmat format is the following :
 *       Pmat[0] = Pauli matrix I
 *       Pmat[1] = Pauli matrix X
 *       Pmat[2] = Pauli matrix Y
 *       Pmat[3] = Pauli matrix Z
 * @return Hamiltonian for self interaction (Pmat * Tensor * Pmat)
 * @note The dimension of Returned Hamiltonian is
 *      the same dimension as the single spin
*/
MatrixXcd calHamiltonianSelfInt(MatrixXcd* Pmat, MatrixXcd Tensor){

    int dimrow = Pmat[0].rows();
    int dimcol = Pmat[0].cols();

    MatrixXcd H = MatrixXcd::Zero(dimrow,dimcol);

    for (int i=0; i<3; i++){
        for (int j=0; j<3; j++){
            H += Tensor(i,j) * Pmat[i+1] * Pmat[j+1];
        }
    }

    return H;
}

MatrixXcd Pauli_matrix_I(int mSize){
    MatrixXcd II_h(mSize,mSize);
    for(int a=0; a<mSize; a++){
        for(int b=0; b<mSize; b++){
            if(a==b){
                II_h(a,b)=1;
            }else{
                II_h(a,b)=0;
            }
        }
    }
    return II_h;
}

MatrixXcd Pauli_matrix_X(int mSize){
    MatrixXcd Ix_h(mSize,mSize);
    for(int a=0; a<mSize; a++){
        for(int b=0; b<mSize; b++){
            if((a == (b+1)) || ((a+1)==b)){
                Ix_h(a,b)=sqrt((((double)mSize-1)/2+1)*((a+1)+(1+b)-1)-(a+1)*(b+1))/2;
            }else{
                Ix_h(a,b)=0;
            }
        }
    }
    return Ix_h;
}

MatrixXcd Pauli_matrix_Y(int mSize){
    MatrixXcd Iy_h(mSize,mSize);
    for(int a=0; a<mSize; a++){
        for(int b=0; b<mSize; b++){
            if(a == (b+1)){
                Iy_h(a,b)=doublec(0,sqrt((((double)mSize-1)/2+1)*((a+1)+(1+b)-1)-(a+1)*(b+1))/2);
            }else if((a+1)==b){
                Iy_h(a,b)=doublec(0,-sqrt((((double)mSize-1)/2+1)*((a+1)+(1+b)-1)-(a+1)*(b+1))/2);
            }else{
                Iy_h(a,b)=0;
            }
        }
    }
    return Iy_h;
}

MatrixXcd Pauli_matrix_Z(int mSize){
    MatrixXcd Iz_h(mSize,mSize);
    for(int a=0; a<mSize; a++){
        for(int b=0; b<mSize; b++){
            if(a==b){
                Iz_h(a,b)=(1-(a+1)+((double)mSize-1)/2);
            }else{
                Iz_h(a,b)=0;
            }
        }
    }
    return Iz_h;
}

MatrixXcd* getPauliOperators(int mSize){
    MatrixXcd* spinVector = new MatrixXcd[4];
    spinVector[0] = Pauli_matrix_I(mSize);
    spinVector[1] = Pauli_matrix_X(mSize);
    spinVector[2] = Pauli_matrix_Y(mSize);
    spinVector[3] = Pauli_matrix_Z(mSize);
    return spinVector;
}


MatrixXcd* getGeneralPauliOperators(MatrixXcd alpha, MatrixXcd beta){

    // printMatrixXcd(alpha,"alpha");
    // printMatrixXcd(beta,"beta");

    int dimrow = alpha.rows() * beta.cols();
    int dimcol = alpha.cols() * beta.rows();
    

    MatrixXcd* sigma = new MatrixXcd[4];

    // MatrixXcd gg = kron(alpha,alpha.adjoint());
    MatrixXcd gg = alpha * alpha.adjoint();
    MatrixXcd ge = alpha * beta.adjoint();
    MatrixXcd ee = beta * beta.adjoint();
    MatrixXcd eg = beta * alpha.adjoint();
    
    sigma[0] = MatrixXcd::Identity(dimrow,dimcol);
    sigma[1] = MatrixXcd::Zero(dimrow,dimcol);
    sigma[2] = MatrixXcd::Zero(dimrow,dimcol);
    sigma[3] = MatrixXcd::Zero(dimrow,dimcol);

    sigma[1] = doublec(1.0, 0.0) * ge + doublec(1.0,0.0) * eg; //   |01><10| +  |10><01|
    sigma[2] = doublec(0.0,-1.0) * ge + doublec(0.0,1.0) * eg; // -i|01><10| + i|10><01|
    sigma[3] = doublec(1.0, 0.0) * gg - doublec(1.0,0.0) * ee; //   |00><00| -  |11><11|

    // printMatrixXcd(sigma[0]*sigma[0],"sigma[0]");
    // printMatrixXcd(sigma[1]*sigma[1],"sigma[1]");
    // printMatrixXcd(sigma[2]*sigma[2],"sigma[2]");
    // printMatrixXcd(sigma[3]*sigma[3],"sigma[3]");

    // isInvolutory(sigma[0]);
    // isInvolutory(sigma[1]);
    // isInvolutory(sigma[2]);
    // isInvolutory(sigma[3]);
    
    return sigma;
}

//In mathematics, an involutory matrix is a square matrix that is its own inverse. 
// Mat *  Mat = Identity
bool isInvolutory(MatrixXcd A){

    if ((A.rows() != A.cols())){
        perror("Error(isInvolutory): Matrix is not square");
        exit(EXIT_FAILURE);
    }

    MatrixXcd B = A * A;
    MatrixXcd I = MatrixXcd::Identity(A.rows(),A.cols());

    // printMatrixXcd(B,"B");
    // printMatrixXcd(I,"I");
    if ((isSame(B,I))){
        return true;
    }
    else{
        perror("Error(isInvolutory): Matrix is not involutory");
        exit(EXIT_FAILURE);
    }
}

bool isSame(MatrixXcd A, MatrixXcd B){
    if ((A.rows() == B.rows()) && (A.cols() == B.cols())){
        for (int i=0; i<A.rows(); i++){
            for (int j=0; j<A.cols(); j++){
                if (A(i,j) != B(i,j)){
                    return false;
                }
            }
        }
        return true;
    }
    else{
        return false;
    }
}

// sort in desending order
int* getIndexInOrder(VectorXcd eigenValues){
    int dim = eigenValues.size();

    int* idx = allocInt1d(dim);

    for (int i=0; i<dim; i++){
        idx[i] = i;
    }

    for (int i=0; i<dim; i++){
        for (int j=i+1; j<dim; j++){
            if (eigenValues(idx[i]).real() < eigenValues(idx[j]).real()){
                int tmp = idx[i];
                idx[i] = idx[j];
                idx[j] = tmp;
            }
        }
    }

    return idx;
}

VectorXcd sortEigenValues(VectorXcd eigenValues, int* idx){

    int dim = eigenValues.size();

    VectorXcd sortedEigenValues = VectorXcd::Zero(dim);

    for (int i=0; i<dim; i++){
        sortedEigenValues(i) = eigenValues(idx[i]);
    }

    return sortedEigenValues;
}

MatrixXcd sortEigenVectors(MatrixXcd eigenVectors, int* idx){

    int dimrow = eigenVectors.rows();
    int dimcol = eigenVectors.cols();

    MatrixXcd sortedEigenVectors = MatrixXcd::Zero(dimrow,dimcol);

    for (int i=0; i<dimcol; i++){
        for (int j=0; j<dimrow; j++){
            sortedEigenVectors(j,i) = eigenVectors(j,idx[i]);
        }
    }

    return sortedEigenVectors;
}

