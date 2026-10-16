/**
 * @file memory_estimate.h
 * @brief The memory check before the calculation: the peak memory per rank against what the node allows.
 */
#ifndef __CCEX_MEMORY_ESTIMATE_H_
#define __CCEX_MEMORY_ESTIMATE_H_
// A67: the memory check before the calculation (no input key). See memory_estimate.cpp.
#include "physical/qubit.h"
#include "physical/bath.h"
#include "physical/defect.h"
#include "physical/pulse.h"
#include "common/output.h"
#include "input/general.h"

// The peak memory of one rank: what it holds now (VmRSS) and what its largest cluster needs on the path the run takes for it
// (memory_estimate.cpp), against the memory one rank of its node may use (Slurm's SLURM_MEM_PER_NODE or SLURM_MEM_PER_CPU,
// else MemAvailable, shared by the ranks of the node). Rank 0 prints one line; a rank whose estimate passes what it may use
// stops the run before any cluster is calculated (ccexRefuse, every rank). Collective over MPI_COMM_WORLD.
void checkClusterMemory(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, Output* op,
                        int order, int*** localclusters);
#endif
