#ifndef __CCEX_CLUSTER_DSJ_H_
#define __CCEX_CLUSTER_DSJ_H_


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "general.h"
#include "bath.h"
#include "utilities.h"
#include "cluster.h"

void clusterizeDsjItb(Cluster* cls, int*** cmap, float*** stmap, BathArray* ba);
//double calCouplingStrength(double* spin1, double* spin2);
void makeDuo13CwithCS(double*** Duo13CwithCS, int nspin, int*** cmap, float*** stmap);
void jointSpinSpin(int** group, int len, int spin1, int spin2);
void jointGroupSpin(int** group, int spin);
void jointGroupGroup(int** group1, int** group2);
void updateGroupChecker(int*** GroupChecker, int* newGroup, int newGroupNumber);
int Delete_Int_2dArr(int*** arr, int row, int col, int delrow);
int addNewGroup(int*** Groups, int* group, int order);
void addInterBathPair(int*** Pairs, int spin1, int spin2);
//void writeClusters(int** arr , int order,const char* fname);
//void writeBathSpins(InputData* ID, Gyro* Gy);
//void createDisjointClusters(CCE_group* CCE, InputData* ID, int order,int clusalgo, int nprocess, int rank);

void sort2DArray(int*** arr,int ROWS, int COLS, int col);
void QuickSort_Double_2d(double ***arr, int col, int left, int right);
//void writeCoherenceEachCluster(double _Complex* arr, int len, const char* fname);
void shrinkClusterO1(Cluster* cls);


#endif // __CCEX_CLUSTER_H_
