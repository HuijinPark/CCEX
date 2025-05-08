#include "../include/hamiltonian.h"
#include "../include/utilities.h"
#include "../include/memory.h"

bool verbosity = false;

int main(){

    MatrixXcd I = Pauli_matrix_I(2);
    MatrixXcd X = Pauli_matrix_X(2);
    MatrixXcd Y = Pauli_matrix_Y(2);
    MatrixXcd Z = Pauli_matrix_Z(2);
    printInlineMatrixXcd("I",I);
    printInlineMatrixXcd("X",X);
    printInlineMatrixXcd("Y",Y);
    printInlineMatrixXcd("Z",Z);
    

    MatrixXcd alpha = MatrixXcd::Random(2,2);
    MatrixXcd beta = MatrixXcd::Random(2,2);
    MatrixXcd* pauligen = getGeneralPauliOperators(alpha,beta);
    printInlineMatrixXcd("I",pauligen[0]);
    printInlineMatrixXcd("X",pauligen[1]);
    printInlineMatrixXcd("Y",pauligen[2]);
    printInlineMatrixXcd("Z",pauligen[3]);

    MatrixXcd* pauli = getPauliOperators(3);
    printInlineMatrixXcd("I",pauli[0]);
    printInlineMatrixXcd("X",pauli[1]);
    printInlineMatrixXcd("Y",pauli[2]);
    printInlineMatrixXcd("Z",pauli[3]);

    double gamma = GAMMA_ELECTRON;
    float bfield[3] = {30.0, 20.0, 100.0};
    double detuning = 0.0;
    double overhauser = 0.0;
    double xyz1[3] = {1.0, 2.0, 3.0};
    double xyz2[3] = {4.0, 5.0, 6.0};

    MatrixXcd zvec = calZeemanVector(gamma, bfield);
    MatrixXcd detunevec = calDetuningVector(detuning);
    MatrixXcd overvec = calOverhauserVector(overhauser);
    MatrixXcd pdtensor = calPointDipoleTensor(xyz1, xyz2, gamma, gamma);
    printInlineMatrixXcd("Zeeman Vector",zvec);
    printInlineMatrixXcd("Detuning Vector",detunevec);
    printInlineMatrixXcd("Overhauser Vector",overvec);
    printInlineMatrixXcd("Point Dipole Tensor",pdtensor);

    // MatrixXcd calHamiltonianHeteroInt(MatrixXcd** pmats, MatrixXcd tensor, int nspin, int ispin, int jspin);
    // MatrixXcd calHamiltonianSelfInt(MatrixXcd* pmat, MatrixXcd tensor);
    // MatrixXcd calHamiltonianSingleInt(MatrixXcd Vector, MatrixXcd* pmat1);
    // MatrixXcd expandHamiltonian(MatrixXcd** pmats, MatrixXcd Hi, int nspin, int ispin);

    // MatrixXcd ZeemanInt(double gamma, float* B0, MatrixXcd* pmat, int m);
    // MatrixXcd kron(MatrixXcd a, MatrixXcd b);

    // int* getIndexInOrder(VectorXcd eigenValues);
    // VectorXcd sortEigenValues(VectorXcd eigenValues, int* idx);
    // MatrixXcd sortEigenVectors(MatrixXcd eigenVectors, int* idx);

    // bool isInvolutory(MatrixXcd a);
    // bool isSame(MatrixXcd a, MatrixXcd b);

    delete[] pauligen;
    delete[] pauli;

    return 0;
}