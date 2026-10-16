/**
 * @file mpi_seed.h
 * @brief The random seed, chosen on rank 0 and shared with every rank.
 */
#ifndef CCEX_MPI_SEED_H
#define CCEX_MPI_SEED_H
#include <mpi.h>
#include <ctime>
// The clock is read only on rank zero; explicit seeds follow the same collective.
inline int resolveMpiSeed(int requested, MPI_Comm comm,
                         std::time_t (*clock)(std::time_t*)=std::time) {
    int localRank=0;MPI_Comm_rank(comm,&localRank);
    int seed=requested;
    if(localRank==0 && seed==-1)seed=static_cast<int>(clock(nullptr));
    MPI_Bcast(&seed,1,MPI_INT,0,comm);
    return seed;
}
#endif
