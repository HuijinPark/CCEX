/**
 * @file hamiltonian.cpp
 * @brief Spin operators and the single-spin and two-spin Hamiltonian terms.
 */
#include "../../include/calculation/hamiltonian.h"
#include "../../include/common/memory.h"
#include <iostream>

/**
 * @brief Zeeman Hamiltonian of one spin [rad kHz]
 * @details H = -gamma B0 . S = -gamma (B0x Sx + B0y Sy + B0z Sz)
 * @param[in] gamma gyromagnetic ratio [radkHz/G]
 * @param[in] B0 magnetic field, 3 components [G]
 * @param[in] Pmat spin operators of the spin: Pmat[0] = I, Pmat[1..3] = Sx, Sy, Sz
 * @param[in] m dimension of the spin space (2S+1)
 * @return m x m Zeeman Hamiltonian [rad kHz]
 */
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

/**
 * @brief Point-dipole interaction tensor [rad kHz], H = S1 . Adip . S2
 * @param[in] xyz1 position of spin 1 [angstrom]
 * @param[in] xyz2 position of spin 2 [angstrom]
 * @param[in] gamma1 gyromagnetic ratio of spin 1 [radkHz/G]
 * @param[in] gamma2 gyromagnetic ratio of spin 2 [radkHz/G]
 * @return 3x3 point-dipole tensor [rad kHz]
 * @details Adip = (mu_0 / 4 pi) * hbar * gamma1 * gamma2 / r^3 * (I - 3 n n^T), n = (xyz2 - xyz1)/r,
 *          an angular frequency (the energy divided by hbar). With gamma in radkHz/G and r in angstrom
 *          the prefactor (mu_0 / 4 pi) hbar = 1e-7 * 1.0546e-34 J s, times 1e14 from the two gammas
 *          (radkHz/G = 1e7 rad/(s T)) and 1e30 from 1/r^3, is 1.0546e3 rad/s = 1.0546 rad kHz:
 *          the code's DS = H_BAR * gamma1 * gamma2 / r^3 with H_BAR = 1.054571729.
 */
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
 * @brief A (d_i d_j)-sized operator on spins i and j, identities on the rest
 * @details Row r of the register keeps every other spin's state, so it only reaches the
 *          d_i*d_j columns that differ from r in spins i and j.
 * @param[in] Pmats spin operators of every site; Pmats[k][0] = I fixes the dimension of site k
 * @param[in] nSpin number of sites in the register
 * @param[in] iSpin first site (0-based)
 * @param[in] jSpin second site (0-based)
 * @param[in] small the operator on the space of sites i and j, (d_i d_j) x (d_i d_j)
 * @return the full-register operator
*/
static MatrixXcd embedTwoSites(MatrixXcd** Pmats, int nSpin, int iSpin, int jSpin, const MatrixXcd& small){
    long long dim = 1, si = 1, sj = 1;
    for (int k=0; k<nSpin; k++) dim *= Pmats[k][0].rows();
    for (int k=iSpin+1; k<nSpin; k++) si *= Pmats[k][0].rows();
    for (int k=jSpin+1; k<nSpin; k++) sj *= Pmats[k][0].rows();
    const int di = Pmats[iSpin][0].rows(), dj = Pmats[jSpin][0].rows();
    MatrixXcd H = MatrixXcd::Zero(dim, dim);
    for (long long r=0; r<dim; r++){
        const int x = int((r/si)%di), u = int((r/sj)%dj);
        const long long base = r - x*si - u*sj;
        for (int y=0; y<di; y++) for (int v=0; v<dj; v++)
            H(r, base + y*si + v*sj) = small(x*dj+u, y*dj+v);
    }
    return H;
}

// sum_ab T_ab S_i^a S_j^b as a (d_i d_j)-sized matrix on spins i and j alone.
// The entries are what nine Kronecker chains gave: T_ab * (S_i^a(x,y) * S_j^b(u,v)),
// summed over ab in the same order (the identity factors were exact multiplications by 1).
static MatrixXcd twoSiteOperator(MatrixXcd** Pmats, const MatrixXcd& Tensor, int nSpin, int iSpin, int jSpin){
    if (iSpin<0 || jSpin<0 || iSpin>=nSpin || jSpin>=nSpin || iSpin==jSpin){
        perror("Error(calHamiltonianHeteroInt): spin indices out of range");
        exit(EXIT_FAILURE);
    }
    const int di = Pmats[iSpin][0].rows(), dj = Pmats[jSpin][0].rows();
    MatrixXcd small = MatrixXcd::Zero(di*dj, di*dj);
    for (int x=0; x<di; x++) for (int y=0; y<di; y++)
    for (int u=0; u<dj; u++) for (int v=0; v<dj; v++){
        doublec value = Tensor(0,0) * (Pmats[iSpin][1](x,y) * Pmats[jSpin][1](u,v));
        value += Tensor(0,1) * (Pmats[iSpin][1](x,y) * Pmats[jSpin][2](u,v));
        value += Tensor(0,2) * (Pmats[iSpin][1](x,y) * Pmats[jSpin][3](u,v));
        value += Tensor(1,0) * (Pmats[iSpin][2](x,y) * Pmats[jSpin][1](u,v));
        value += Tensor(1,1) * (Pmats[iSpin][2](x,y) * Pmats[jSpin][2](u,v));
        value += Tensor(1,2) * (Pmats[iSpin][2](x,y) * Pmats[jSpin][3](u,v));
        value += Tensor(2,0) * (Pmats[iSpin][3](x,y) * Pmats[jSpin][1](u,v));
        value += Tensor(2,1) * (Pmats[iSpin][3](x,y) * Pmats[jSpin][2](u,v));
        value += Tensor(2,2) * (Pmats[iSpin][3](x,y) * Pmats[jSpin][3](u,v));
        small(x*dj+u, y*dj+v) = value;
    }
    return small;
}

/**
 * @brief Hamiltonian for interacting two spins
 * @details The Hamiltonian is given by S_i . Tensor . S_j = sum_ab Tensor_ab S_i^a S_j^b
 * @param[in] Pmats spin operators of every site: Pmats[k][0] = I, Pmats[k][1..3] = Sx, Sy, Sz
 * @param[in] Tensor 3x3 interaction tensor
 * @param[in] nSpin number of sites in the register
 * @param[in] iSpin first spin (0-based)
 * @param[in] jSpin second spin (0-based)
 * @return the two-spin Hamiltonian on the whole register (dimension: product of all site dimensions)
*/
MatrixXcd calHamiltonianHeteroInt(MatrixXcd** Pmats, MatrixXcd Tensor, int nSpin, int iSpin, int jSpin){

    // sum_ab T_ab S_i^a S_j^b on the whole register. Every other spin contributes an
    // identity, so the operator is a (d_i d_j)-sized matrix placed at the positions
    // where all the other spins keep their state: built here, then scattered into the
    // full matrix - D * d_i * d_j work instead of forming nine D x D Kronecker chains.
    return embedTwoSites(Pmats, nSpin, iSpin, jSpin, twoSiteOperator(Pmats, Tensor, nSpin, iSpin, jSpin));
}

/**
 * @brief Adds the two-spin Hamiltonian of calHamiltonianHeteroInt to H
 * @details H += calHamiltonianHeteroInt(...), without the full-size temporary: only the
 *          D * d_i * d_j entries the operator has are touched. The sum a caller builds is
 *          the same, bit for bit, as adding the full matrix - the entries outside the
 *          pattern are +0, and adding +0 changes nothing in a sum that started from +0
 *          (such a sum is never -0). For a register of 14 spin-1/2 the full matrix is
 *          4.3 GB, and one was made, zeroed and added for each of the 91 pairs.
*/
void addHamiltonianHeteroInt(MatrixXcd& H, MatrixXcd** Pmats, const MatrixXcd& Tensor, int nSpin, int iSpin, int jSpin){
    const MatrixXcd small = twoSiteOperator(Pmats, Tensor, nSpin, iSpin, jSpin);
    long long dim = 1, si = 1, sj = 1;
    for (int k=0; k<nSpin; k++) dim *= Pmats[k][0].rows();
    for (int k=iSpin+1; k<nSpin; k++) si *= Pmats[k][0].rows();
    for (int k=jSpin+1; k<nSpin; k++) sj *= Pmats[k][0].rows();
    const int di = Pmats[iSpin][0].rows(), dj = Pmats[jSpin][0].rows();
    if (H.rows()!=dim || H.cols()!=dim){
        perror("Error(addHamiltonianHeteroInt): Dimension of H is not matched");
        exit(EXIT_FAILURE);
    }
    for (long long r=0; r<dim; r++){
        const int x = int((r/si)%di), u = int((r/sj)%dj);
        const long long base = r - x*si - u*sj;
        for (int y=0; y<di; y++) for (int v=0; v<dj; v++)
            H(r, base + y*si + v*sj) += small(x*dj+u, y*dj+v);
    }
}

/**
 * @brief Projected Hamiltonian for interacting two spins
 * @details This related to the hyperfine mediated interaction term
 *          This is exactly the same to calHamiltonianHeteroInt but the tensor is the following:
 *          (Sab*Ai)^T * (Sba*Aj) 
*/
MatrixXcd* calHamiltonianHeteroInt_o2(MatrixXcd Sab, MatrixXcd Sba, MatrixXcd** Pmats, MatrixXcd tensor_ib, MatrixXcd tensor_jb, int nSpin, int iSpin, int jSpin){

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
    
    MatrixXcd SabAi_alpha(1,3); 	MatrixXcd SabAj_beta(1,3);	
    MatrixXcd SbaAj_alpha(1,3);		MatrixXcd SbaAi_beta(1,3);		
    MatrixXcd Aij_alpha(3,3);		MatrixXcd Aij_beta(3,3);
    MatrixXcd Ai(3,3);
    MatrixXcd Aj(3,3);
    Ai = tensor_ib; 
    Aj = tensor_jb; 
    
    //Aij_alpha
    // = (Sab * Ai)^T * (Sba * Aj)
    SabAi_alpha = Sab.transpose() * Ai;
    SbaAj_alpha = Sba.transpose() * Aj;
    Aij_alpha   = (SabAi_alpha.transpose()) * SbaAj_alpha;
    
    //printf("Aij_alpha = SabAi_alpha.T * SbaAj_alpha\n");
    //std::cout << Aij_alpha << std::endl;
    
    //Aij_beta
    // = (Sba * Ai)^T * (Sab * Aj)
    SbaAi_beta  = Sba.transpose() * Ai;
    SabAj_beta  = Sab.transpose() * Aj;
    Aij_beta    = (SbaAi_beta.transpose()) * SabAj_beta;
    
    //printf("Aij_beta = SbaAi_beta.T * SabAj_beta\n");
    //std::cout << Aij_beta << std::endl;
    
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
    MatrixXcd* Hij_medi = new MatrixXcd[2];
    MatrixXcd Sab_Ai_Ii_Sba_Aj_Ij = MatrixXcd::Zero(dimrow,dimcol);
    MatrixXcd Sba_Ai_Ii_Sab_Aj_Ij = MatrixXcd::Zero(dimrow,dimcol);

    Sab_Ai_Ii_Sba_Aj_Ij = Aij_alpha(0,0) * SixSjx + Aij_alpha(0,1) * SixSjy + Aij_alpha(0,2) * SixSjz +
                          Aij_alpha(1,0) * SiySjx + Aij_alpha(1,1) * SiySjy + Aij_alpha(1,2) * SiySjz +
                          Aij_alpha(2,0) * SizSjx + Aij_alpha(2,1) * SizSjy + Aij_alpha(2,2) * SizSjz;

    Sba_Ai_Ii_Sab_Aj_Ij = Aij_beta(0,0) * SixSjx + Aij_beta(0,1) * SixSjy + Aij_beta(0,2) * SixSjz +
                          Aij_beta(1,0) * SiySjx + Aij_beta(1,1) * SiySjy + Aij_beta(1,2) * SiySjz +
                          Aij_beta(2,0) * SizSjx + Aij_beta(2,1) * SizSjy + Aij_beta(2,2) * SizSjz;

    Hij_medi[0] = Sab_Ai_Ii_Sba_Aj_Ij;
    Hij_medi[1] = Sba_Ai_Ii_Sab_Aj_Ij;

    return Hij_medi;
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

    // Hi on spin iSpin, identities elsewhere: Hi's entries at the positions where the
    // other spins keep their state (what the Kronecker chain gave, without it).
    if (iSpin<0 || iSpin>=nSpin || Hi.rows()!=Pmats[iSpin][0].rows() || Hi.cols()!=Pmats[iSpin][0].cols()){
        perror("Error(expandHamiltonian): Dimension of expandedHamiltonian is not matched");
        exit(EXIT_FAILURE);
    }
    long long dim = 1, stride = 1;
    for (int k=0; k<nSpin; k++) dim *= Pmats[k][0].rows();
    for (int k=iSpin+1; k<nSpin; k++) stride *= Pmats[k][0].rows();
    const int d = Hi.rows();
    MatrixXcd H = MatrixXcd::Zero(dim, dim);
    for (long long r=0; r<dim; r++){
        const int a = int((r/stride)%d);
        const long long base = r - a*stride;
        for (int b=0; b<d; b++) H(r, base + b*stride) = Hi(a,b);
    }
    return H;
}

/**
 * @brief Adds the expanded single-spin Hamiltonian of expandHamiltonian to H
 * @details H += expandHamiltonian(...), without the full-size temporary (see
 *          addHamiltonianHeteroInt): the same sum, bit for bit.
*/
void addExpandedHamiltonian(MatrixXcd& H, MatrixXcd** Pmats, const MatrixXcd& Hi, int nSpin, int iSpin){
    if (iSpin<0 || iSpin>=nSpin || Hi.rows()!=Pmats[iSpin][0].rows() || Hi.cols()!=Pmats[iSpin][0].cols()){
        perror("Error(addExpandedHamiltonian): Dimension of expandedHamiltonian is not matched");
        exit(EXIT_FAILURE);
    }
    long long dim = 1, stride = 1;
    for (int k=0; k<nSpin; k++) dim *= Pmats[k][0].rows();
    for (int k=iSpin+1; k<nSpin; k++) stride *= Pmats[k][0].rows();
    const int d = Hi.rows();
    if (H.rows()!=dim || H.cols()!=dim){
        perror("Error(addExpandedHamiltonian): Dimension of H is not matched");
        exit(EXIT_FAILURE);
    }
    for (long long r=0; r<dim; r++){
        const int a = int((r/stride)%d);
        const long long base = r - a*stride;
        for (int b=0; b<d; b++) H(r, base + b*stride) += Hi(a,b);
    }
}

/**
 * @brief Hamiltonian for single spin
 * @details The Hamiltonian is given by (Vector * Pmat1)
 * @param[in] Vector 1x3 or 3x1 matrix
 * @param[in] Pmat1 spin operators of the single spin (any spin S, not only Pauli matrices)
 * @note Pmat1 format is the following :
 *       Pmat1[0] = I
 *       Pmat1[1] = Sx
 *       Pmat1[2] = Sy
 *       Pmat1[3] = Sz
 * @return Hamiltonian for single spin  (Vector * Pmat1)
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

// Operator projection
MatrixXcd projOperators(MatrixXcd* sigma, MatrixXcd psia, MatrixXcd psib){

    MatrixXcd proj_sigma = MatrixXcd::Zero(3,1);
    proj_sigma(0,0) = ((psia.transpose()) * sigma[1] * psib)(0,0); //Sab_x
    proj_sigma(1,0) = ((psia.transpose()) * sigma[2] * psib)(0,0); //Sab_y
    proj_sigma(2,0) = ((psia.transpose()) * sigma[3] * psib)(0,0); //Sab_z

    return proj_sigma;
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

