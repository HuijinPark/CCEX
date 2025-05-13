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

//MatrixXcd cal_Upulse_withDetuning(QubitArray* qa, Pulse* pulse, double angle, MatrixXcd sigma, MatrixXcd sigma_z);
void calUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse);
void calbUpulses(MatrixXcd* Upulses, BathArray* ba, Pulse* pulse);
void calTDUpulses(MatrixXcd* Upulses, QubitArray* qa, Pulse* pulse, MatrixXcd Hq);
void calculate(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cc, Pulse* pulse, Cluster* cluster, Output* output, int*** localClusters);

BathArray* createBathArray(int* cluster, int nspin, BathArray* ba, DefectArray* dfa, int nqubit);

MatrixXcd* calCoherenceGcce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls, Output* op);
MatrixXcd* calCoherenceCce(QubitArray* qa, BathArray* ba, Config* cnf, Pulse* pls);

MatrixXcd calPropagatorGcce(QubitArray* qa, BathArray* ba, MatrixXcd Htot, Pulse* pls, double tfree, MatrixXcd* Upulses, MatrixXcd* bUpulses);

MatrixXcd HamilQubit(QubitArray* qa, BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilBath(BathArray* ba, MatrixXcd** sigmas, Config* cnf);
MatrixXcd HamilQubitBath(QubitArray* qa, BathArray* ba, MatrixXcd** qsigmas, MatrixXcd** bsigmas, Config* cnf);
MatrixXcd* HamilQubitBathSecularApp(QubitArray* qa, BathArray* ba, MatrixXcd** bsigmas, Config* cnf);
double calQubit_Energy_difference(QubitArray* qa, MatrixXcd Hq);
MatrixXcd cal_TDUpulse(QubitArray* qa, MatrixXcd Hq, Pulse* pulse, MatrixXcd ox, MatrixXcd oy, char axis, double pulse_angle, double Ediff, double detuning_factor, double pulse_time, int Nstep);


MatrixXcd** gatherSigmas(MatrixXcd** qsigmas, MatrixXcd** bsigmas, int nqubit, int nspin);
MatrixXcd** gatherSigmas_singlequbit(MatrixXcd* qsigma, MatrixXcd** bsigmas, int nqubit, int nspin);
int close_state_index(MatrixXcd state, MatrixXcd eigenVectors);

#endif // __CCEX_SIMULATOR_H_
