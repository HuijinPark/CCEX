#ifndef __CCEX_READER_H_
#define __CCEX_READER_H_

#include "qubit.h"
#include "bath.h"
#include "general.h"


void readQubitfile(QubitArray* qa, Config* cnf);
void readBathfiles(BathArray* ba, QubitArray* qa, Config* cnf);
void readGyrofile(BathArray* ba, Config* cnf);
void readHftensorfile(BathArray* ba, QubitArray* qa, Config* config);
// void readQdtensorfile(BathArray* ba, QubitArray* qa, Config* config);


// Read tensor files
bool READ_BD_vertex(const char* inputfile, double*** vertex, double*** center, double*** normal, char** vertex_condi);
void READ_Tensor(const char* inputfile, double*** Tensor,int numCol);
void READ_Tensor_etc(const char* inputfile, char** names, int nspecies, double*** Tensor,char* condition, int numCol);
void READ_Tensor_const(const char* inputfile, char** names, int nspecies, double** Array,char* condition);
bool READ_Tensor_array(const char* inputfile, double(*Array)[3],char* condition,int numCol);
int READ_Tensor_ver(const char* inputfile,double* SpinFactor, double DefectTotSpin, double* CorrTotSpin);

void printHfInfo_version(int version, double DefectTotSpin, double CorrTotSpin, double SpinFactor);
void printHfInfo_BD(double** vertex, double** center, double** normal, double MinDif[3], double MaxDif[3], bool Usingvertex);
void printHfInfo_etc(double** A_Etc, double* A_Gfactor, char** names, int nspecies, int mode);
void printHfInfo_tensor(double** AtensorArray, int mode);


bool CheckPosition(const double Posi[3],double* refPosi);
bool FIND_AtomPosi(const double refPosi[3], double** TensorValue, int* num);
bool CheckBD_Range(const double Posi[3],const double minRange[3],const double maxRange[3],double err);
bool CheckBD_vertex(const double Posi[3],double** vertex, double** center, double** normal, double err);
bool CheckingUniqeBD(const double Posi[3], double** center, double** normal);
void CreatePlaneInfo(double** vertex, double*** center, double*** normal);
double ReDefinediff(double difXYZ,const  double minRange,const double maxRange,const double Copy_Length);
void CalCenter(double** center, double p1[], double p2[], double p3[], double p4[]);
void CalNormal(double** normal, double u[], double v[]);
void vector_diff(double (*result)[3], double p1[], double p2[]);
double vector_dot(double u[], double v[]);
void vector_cross(double** cross_P, double u[], double v[], bool norm);


#endif // __CCEX_READER_H_