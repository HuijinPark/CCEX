#ifndef __CCEX_HAMILTONIAN_H_
#define __CCEX_HAMILTONIAN_H_

#include "utilities.h"

MatrixXcd Pauli_matrix_I(int size);
MatrixXcd Pauli_matrix_X(int size);
MatrixXcd Pauli_matrix_Y(int size);
MatrixXcd Pauli_matrix_Z(int size);
MatrixXcd* getGeneralPauliOperators(MatrixXcd alpha, MatrixXcd beta);
MatrixXcd* getPauliOperators(int size);

MatrixXcd calZeemanVector(double gamma, float* bfield);
MatrixXcd calDetuningVector(double detuning);
MatrixXcd calOverhauserVector(double overhauser);
MatrixXcd calPointDipoleTensor(double xyz1[3], double xyz2[3], double gamma1, double gamma2);

MatrixXcd calHamiltonianHeteroInt(MatrixXcd** pmats, MatrixXcd tensor, int nspin, int ispin, int jspin);
MatrixXcd calHamiltonianSelfInt(MatrixXcd* pmat, MatrixXcd tensor);
MatrixXcd calHamiltonianSingleInt(MatrixXcd Vector, MatrixXcd* pmat1);
MatrixXcd expandHamiltonian(MatrixXcd** pmats, MatrixXcd Hi, int nspin, int ispin);

MatrixXcd ZeemanInt(double gamma, float* B0, MatrixXcd* pmat, int m);

int* getIndexInOrder(VectorXcd eigenValues);
VectorXcd sortEigenValues(VectorXcd eigenValues, int* idx);
MatrixXcd sortEigenVectors(MatrixXcd eigenVectors, int* idx);

bool isInvolutory(MatrixXcd a);
bool isSame(MatrixXcd a, MatrixXcd b);

#endif // __CCEX_HAMILTONIAN_H_