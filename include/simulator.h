#ifndef __CCEX_SIMULATOR_H_
#define __CCEX_SIMULATOR_H_

#include <unsupported/Eigen/MatrixFunctions>
#include "utilities.h"
#include "qubit.h"
#include "bath.h"
#include "defect.h"
#include "output.h"
#include "pulse.h"
#include "cluster.h"
#include "general.h"

void calculate(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cc, Pulse* pulse, Cluster* cluster, Output* output, int*** localClusters);

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, DefectArray* dfa, int nqubit);

MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Output* op);
MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls);

MatrixXcd calPropagatorGcce(QubitArray* qa, MatrixXcd Htot, Pulse* pls, double tfree);

MatrixXcd HamilQubit(QubitArray* qa, BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilBath(BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilQubitBath(QubitArray* qa, BathArray* ba, MatrixXcd** qsigmas, MatrixXcd** bsigmas, Config* cnf);
MatrixXcd* HamilQubitBathSecularApp(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas, Config* cnf);

#endif // __CCEX_SIMULATOR_H_
