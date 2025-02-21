#ifndef __CCEX_OUTPUT_H_
#define __CCEX_OUTPUT_H_

#include "general.h"
#include "bath.h"
#include "utilities.h"

typedef struct {
   
    // Intermediate Save
    // char* InterSaveClusterFile; /**< File name to save clusters' information */
    // char* InterSaveCoherenceFile; /**< File name to save coherence functions of clusters */
    char savemode[20]; /**< Save mode : all |  | averaged */
    // all : coherence function for each cluster and each state : output = output_{cluster}_{state}_wD/nD
    // normal : coherence function for each state (default) : output = output_{state}_wD/nD
    // avg : averaged coherence function for all states : ouput = output_avg_wD/nD
    // info : normal save mode + deep learning information (Azx, Azz, bathfile) : output = output_{state}_wD/nD

    //////////////////////////////////////////////
    // Final Save (Main result)
    // outfile1. actual result 
    // outfile2. the result to compare them to remove the diverged points
    // char* OutputFileHead; /**< File name to save final result */
    // char* OutputFileTail; /**< File name to save final result */
    // char* OutputFileHeadDisjoint; /**< File name to save final result */
    // char* OutputFileTailDisjoint; /**< File name to save final result */
    char* outfile; // output file name with mode (auto generated)
} Output;

Output* Output_init();
void Output_save(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, float deltat, int istate);
void Output_save_all(Output* op, MatrixXcd* result, int nstep, float deltat,int* cluster, int nspin, int istate);
void Output_save_info(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, float deltat, int istate, double Azx, double Azz, char* bathfile);

void Output_setSavemode(Output* op, char* savemode);
char* Output_getSavemode(Output* op);

void Output_allocOutfile(Output* op);
void Output_freeOutfile(Output* op);
void Output_setOutfile(Output* op, char* outfile);
char* Output_getOutfile(Output* op);

void Output_report(Output* op);

void Output_freeAll(Output* op);



#endif // __CCEX_OUTPUT_H_