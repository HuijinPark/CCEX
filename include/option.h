#ifndef __CCEX_OPTION_H_
#define __CCEX_OPTION_H_

#include "general.h"
#include "qubit.h"
#include "cluster.h"
#include "pulse.h"
#include "output.h"
#include "json.h"
#include "utilities.h"

void cJSON_readOptionConfig(Config* cnf, char* fccein); // general.h
void cJSON_readOptionQubitArray(QubitArray* qa, char* fccein); //qubit.h
void cJSON_readOptionCluster(Cluster* clus, char* fccein); // cluster.h
// void cJSON_readOptionBathArray(BathArray* ba, char* fccein); // bath.h
void cJSON_readOptionPulse(Pulse* pulse, char* fccein); // pulse.h
void cJSON_readOptionOutput(Output* output, char* fccein); // output.h



char* cJSON_ReadFccein(char* fccein);
char* cJSON_ReadFilePath(cJSON* root, char* key, bool _default, char* default_value);
char* cJSON_ReadString(cJSON* root, char* key, bool _default, char* default_value);

int cJSON_ReadInt(cJSON* root, char* key, bool _default, int* default_value);
float cJSON_ReadFloat(cJSON* root, char* key, bool _default, float* default_value);
double cJSON_ReadDouble(cJSON* root, char* key, bool _default, double* default_value);
bool cJSON_ReadBool(cJSON* root, char* key, bool _default, bool default_value);

char** cJSON_ReadFilePath1d(int* length, cJSON* root, char* key, bool _default, char** default_value); // memory free is needed
int* cJSON_ReadInt1d(cJSON* root, char* key, bool _default, int* default_value, int size); // memory free is needed
float* cJSON_ReadFloat1d(cJSON* root, char* key, bool _default, float* default_value, int size); // memory free is needed
double* cJSON_ReadDouble1d(cJSON* root, char* key, bool _default, double* default_value, int size); // memory free is needed
double** cJSON_ReadDouble2d(cJSON* root, char* key, bool _default, double** default_value, int row, int col); // memory free is needed

MatrixXcd cJSON_ReadTensor(cJSON* root, char* key, bool _default, MatrixXcd default_value);

// void checkComparable(Config* cnf, QubitArray* qa, BathArray* ba, Cluster* clus, Pulse* pulse, Output* output);
// DefectTotalSpin vs qa->qubit[i]->spin

#endif // __CCEX_OPTION_H_
