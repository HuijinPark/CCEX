/**
 * @file output.h
 * @brief The Output structure and the writing of the result files.
 */
#ifndef __CCEX_OUTPUT_H_
#define __CCEX_OUTPUT_H_

#include "input/general.h"
#include "physical/bath.h"
#include "common/utilities.h"

/**
 * @brief Where and what a run writes: the result file name, the output folder and the output.savemode setting.
 */
typedef struct {
   
    // Intermediate Save
    // char* InterSaveClusterFile; /**< File name to save clusters' information */
    // char* InterSaveCoherenceFile; /**< File name to save coherence functions of clusters */
    char savemode[32]; /**< output.savemode (its mode): coherence (default), coherence_tagged, cluster_coherence or cluster_density_matrix */
    // coherence : coherence function for each state (default) : output = output_{state}_wD/nD
    // coherence_tagged : coherence + deep learning information (Azx, Azz, bathfile) : output = output_{state}_wD/nD
    // cluster_coherence : coherence function for each cluster and each state : output = output_{cluster}_{state}_wD/nD
    // cluster_density_matrix : cluster_coherence + each cluster's full density matrix (gCCE)

    //////////////////////////////////////////////
    // Final Save (Main result)
    // outfile1. actual result 
    // outfile2. the result to compare them to remove the diverged points
    // char* OutputFileHead; /**< File name to save final result */
    // char* OutputFileTail; /**< File name to save final result */
    // char* OutputFileHeadDisjoint; /**< File name to save final result */
    // char* OutputFileTailDisjoint; /**< File name to save final result */
    char* outfile; /**< result file name without the _wiDiv/_noDiv ending (made from output.file and the run) */
    bool save_nodiv; /**< write the _noDiv file beside the _wiDiv one (output.savemode.nodiv, default true) */
    char outpath[512]; /**< output.path: the folder of the outputs when -o DIR is not given; "" for none */
} Output;

Output* Output_init();
void Output_save(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, double deltat, int istate);
void Output_save_all(Output* op, MatrixXcd* result, int nstep, double deltat,int* cluster, int nspin, int istate);
/** savemode cluster_density_matrix: op with the outfile "{outfile}_{bath indices}_full[_state{istate}]" written into name (MAX_FILEPATH
 *  characters), the name savemode cluster_coherence gives the cluster with "_full" ("{outfile}_full..." for the zeroth cluster, n = 0).
 *  calCoherenceGcce writes the cluster's full density matrices to the Output it is given. */
Output Output_fullOutput(const Output* op, char* name, const int* cluster, int n, int istate);
void Output_save_info(Output* op, MatrixXcd* result_wD, MatrixXcd* result_nD, int nstep, double deltat, int istate, double Azx, double Azz, char* bathfile);

void Output_setSavemode(Output* op, char* savemode);
char* Output_getSavemode(Output* op);

void Output_allocOutfile(Output* op);
void Output_freeOutfile(Output* op);
void Output_setOutfile(Output* op, char* outfile);
char* Output_getOutfile(Output* op);

void Output_report(Output* op);

void Output_freeAll(Output* op);



#endif // __CCEX_OUTPUT_H_