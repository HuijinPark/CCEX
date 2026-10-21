// Simulate ranks crossing a wall-clock second without timing-dependent sleeps.
#include "input/mpi_seed.h"
#include <iostream>
static int localRank=0, calls=0;
static std::time_t rankClock(std::time_t*){++calls;return 1726000000+localRank;}
int main(int argc,char** argv){
    MPI_Init(&argc,&argv);MPI_Comm_rank(MPI_COMM_WORLD,&localRank);
    int failures=0;
    if(resolveMpiSeed(-1,MPI_COMM_WORLD,rankClock)!=1726000000)++failures;
    if(calls!=(localRank==0 ? 1 : 0))++failures;
    const int before=calls;
    if(resolveMpiSeed(123456,MPI_COMM_WORLD,rankClock)!=123456 || before!=calls)++failures;
    int total=0;MPI_Allreduce(&failures,&total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if(localRank==0)std::cout<<(total ? "FAIL" : "PASS")<<" MPI generated/explicit seed and root-only clock\n";
    MPI_Finalize();return total ? 1 : 0;
}
