/**
 * @file memory_estimate.cpp
 * @brief The memory check before the calculation: the peak memory per rank against what the node allows.
 */
// A67 (2026-10-15, user: an automatic memory check before the run, no input key): the peak memory of a rank, estimated
// from its largest cluster on the path the run takes for it, against what one rank of the node may use. The run is stopped
// before any cluster is calculated when the estimate does not fit (multiNV6 filled a node, A60; a large megCCE cluster,
// A61). The counts of dense matrices below were measured against the peak resident memory (VmHWM) of local runs
// ((internal notes) A67) and round up.
#include "../../include/calculation/memory_estimate.h"
#include "../../include/calculation/sparse_propagator.h"
#include "../../include/calculation/gpu/compute_backend.h"
#include "../../include/calculation/gpu/gpu_batch_memory.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <set>
#include <string>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

namespace {

double statusKb(const char* key){            // a line of /proc/self/status, kB; 0 when it cannot be read
#ifdef __APPLE__
    // macOS has no /proc: VmRSS is the task's resident size (until 2026-10-19 it read 0 there, so the estimate left out
    // what the rank held and check_refusals' rj_memory_not_enough was accepted)
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (!strcmp(key,"VmRSS") &&
        task_info(mach_task_self(),MACH_TASK_BASIC_INFO,(task_info_t)&info,&count) == KERN_SUCCESS)
        return double(info.resident_size)/1024.0;
    return 0;
#endif
    double kb = 0;
    if (FILE* fp = fopen("/proc/self/status","r")){
        char line[256];
        const size_t n = strlen(key);
        while (fgets(line,sizeof line,fp)) if (!strncmp(line,key,n) && line[n]==':'){ kb = atof(line+n+1); break; }
        fclose(fp);
    }
    return kb;
}
double memAvailableBytes(){                   // MemAvailable of /proc/meminfo; 0 when it cannot be read
#ifdef __APPLE__
    // macOS: the free and inactive pages, what psutil calls available there
    vm_statistics64_data_t vm;
    mach_msg_type_number_t count = HOST_VM_INFO64_COUNT;
    vm_size_t page = 0;
    if (host_page_size(mach_host_self(),&page) == KERN_SUCCESS &&
        host_statistics64(mach_host_self(),HOST_VM_INFO64,(host_info64_t)&vm,&count) == KERN_SUCCESS)
        return (double(vm.free_count) + double(vm.inactive_count))*double(page);
    return 0;
#endif
    double kb = 0;
    if (FILE* fp = fopen("/proc/meminfo","r")){
        char line[256];
        while (fgets(line,sizeof line,fp)) if (!strncmp(line,"MemAvailable:",13)){ kb = atof(line+13); break; }
        fclose(fp);
    }
    return kb*1024.0;
}
double envNumber(const char* name){           // a positive number in the environment, else 0
    const char* v = getenv(name);
    if (!v || !*v) return 0;
    char* end = nullptr;
    const double x = strtod(v,&end);
    return end != v && x > 0 ? x : 0;
}
// Slurm writes a memory as a number of megabytes, optionally with a unit letter (K, M, G, T)
double slurmMegabytes(const char* name){
    const char* v = getenv(name);
    if (!v || !*v) return 0;
    char* end = nullptr;
    double x = strtod(v,&end);
    if (end == v || !(x > 0)) return 0;
    switch (*end){ case 'K': case 'k': x /= 1024; break; case 'G': case 'g': x *= 1024; break; case 'T': case 't': x *= 1024*1024; break; }
    return x;
}

// What one cluster of bath dimension bdim needs while it is calculated, in bytes, beyond what the rank already holds.
// square: one dense complex matrix of the cluster's space (16 D^2 bytes); liouville: one of its Liouville space (16 D^4).
constexpr long double kLibraries = 8e6;   // what the libraries hold besides the matrices (a few MB at small dimensions)
long double clusterBytes(Config* cnf, Pulse* pls, Output* op, int qdim, double bdim, int nsites){
    const char* kernel = methodKernel(cnf->method);
    const bool gcce = !strcmp(kernel,"gcce"), mecce = !strcmp(kernel,"mecce"), megcce = !strcmp(kernel,"megcce");
    const long double D = (gcce || megcce ? qdim : 1) * (long double)bdim;
    const long double square = 16.0L*D*D;
    const int nstep = Config_getNstep(cnf);
    const int npulse = std::max(0,Pulse_getNpulse(pls));
    const bool manual = !strcasecmp(pls->pulsename,"Manual");
    const bool ensemble = Config_getNstate(cnf) == 0;
    int distinct = 0;                          // the distinct intervals of the sequence (the column blocks, A73 (10))
    for (int i=0; i<=npulse; i++) if (pls->sequence_indices && pls->sequence_indices[i] == i) distinct++;
    if (distinct == 0) distinct = npulse+1;
    if (mecce || megcce){
        // the Lindbladian and, after an odd number of pulses, its flipped form (meCCE) or the pulse superoperators
        // (megCCE), the two half-step propagators, their product and the scaling and squaring of the exponential
        // (measured: meCCE 13.2 at D 32, megCCE 11.6 at D 2 x 16, Hahn)
        const long double liouville = 16.0L*D*D*D*D;
        return (megcce ? 10.0L + npulse + distinct : 14.0L)*liouville + 4.0L*square + kLibraries;
    }
    const bool onGpu = D <= 2147483647.0L && computeMayUseCuda(cnf,int(D));
    const bool allFull = gcce && !strcasecmp(Output_getSavemode(op),"cluster_density_matrix");
    if (onGpu){
        // the host side: the Hamiltonians built and handed to the card, and with cluster_density_matrix every time's
        // matrix (and the read-back of one time batch); the card's own memory is checked by the GPU path
        const int batch = std::max(1,std::min(cnf->gpu_batch_steps,nstep));
        return (gcce ? 5.0L : 3.0L)*square + (allFull ? (long double)(nstep+batch)*square : 0.0L);
    }
    if (D <= 4) return 0;                      // the fixed-size and batched small clusters: nothing a node would notice
    if (sparsePropagatorFor(cnf,long(D))){
        // the Krylov basis (31 vectors), the Hamiltonians in compressed rows (an entry: 16 + 4 bytes; a row holds at most
        // the diagonal, one entry per spin flip and one per pair of flips) and the states
        const long double perRow = 1.0L + 3.0L*nsites + 3.0L*nsites*(nsites-1);
        return 16.0L*D*(31.0L + 8.0L) + 2.0L*20.0L*D*perRow;
    }
    // dense: so many complex D x D matrices at the peak (measured, rounded up; with MKL the eigensolver (zheevd) keeps its
    // vectors and two workspaces of D^2 per solver)
    long double count;
    if (gcce){
        // the Hamiltonian kept, the parts it is built from, the eigenvectors and the solver, rho, U and the pulses expanded
        // to the whole space (npulse), with evolution matrix the free propagators of the intervals (npulse+1)
        // (measured: matrix Hahn 19.7 and 19.0 at D 768 and 1536, matrix CPMG-8 33.7, vector Hahn 11.4 and 10.8)
        const bool matrix = strcasecmp(cnf->evolution,"vector") != 0;
        count = matrix ? 18.0L + npulse + (npulse+1) : 11.0L + npulse;
        if (allFull) count += nstep + 2;       // every time's density matrix, kept to be written (43 at D 384, nstep 20)
    } else if (!ensemble){
        // the two Hamiltonians, the two eigensolvers, M and M^dagger (measured 14.0 and 13.3 at D 1024 and 2048)
        count = 15.0L;
    } else if (manual){
        count = 19.0L + 4.0L*(npulse+1);      // four products per interval (measured 39 and 37, 5 intervals)
    } else {
        // Ramsey, Hahn echo, CPMG (measured 15.3, 17.9 and 24.6 at D 512; 13.8, 15.5 and 22.1 at D 1024)
        count = npulse == 0 ? 16.0L : npulse == 1 ? 18.0L : 25.0L;
    }
    long double bytes = count*square + kLibraries;
    // A62: one state with the time points as matrix columns - per distinct interval two D x columns blocks of phases, and
    // the four blocks of the chain (A73 (10): Manual 17 intervals, D 1024, 64 columns: about 35 MB)
    if (!gcce && !ensemble && cnf->column_evolution && D >= cnf->column_evolution_min_dim){
        const long double columns = std::min(cnf->column_evolution_columns,std::max(1,nstep));
        bytes += 16.0L*D*columns*(2.0L*distinct + 4.0L) + 16.0L*D*4.0L*(npulse+1);
    }
    return bytes;
}

}   // namespace

void checkClusterMemory(QubitArray* qa, BathArray* ba, DefectArray* dfa, Config* cnf, Pulse* pls, Output* op,
                        int order, int*** localclusters){
    // the bath dimension of a site with the defect's sub-spins (as Cluster_costFunction)
    std::vector<double> siteDim(BathArray_getNspin(ba), 1.0);
    for (int i=0; i<BathArray_getNspin(ba); i++){
        siteDim[i] = BathSpin_dim(BathArray_getBath_i(ba,i));
        const int idf = DefectArray_findDefectIndex(dfa,BathArray_getBath_i_name(ba,i));
        for (int isp=0; idf>=0 && isp<DefectArray_getDefect_idf_naddspin(dfa,idf); isp++)
            siteDim[i] *= BathSpin_dim(DefectArray_getSubbath_i_isp(dfa,i,isp));
    }
    const int qdim = QubitArray_dim(qa);
    // the distinct (bath dimension, sites) of this rank's clusters; the path, and so the memory, depends on both
    std::set<std::pair<double,int>> shapes;
    for (int n=1; n<=order; n++)
        for (int ic=1; ic<localclusters[n][0][0]; ic++){
            double d = 1.0;
            for (int j=1; j<=n; j++) d *= siteDim[localclusters[n][ic][j]];
            shapes.insert({d,n});
        }
    const char* kernel = methodKernel(cnf->method);
    const bool qubitInMatrix = !strcmp(kernel,"gcce") || !strcmp(kernel,"megcce");
    long double work = 0;
    double largest = qubitInMatrix ? qdim : 1;   // the zeroth cluster of gCCE and megCCE
    if (qubitInMatrix) work = clusterBytes(cnf,pls,op,qdim,1.0,0);
    for (const auto& s : shapes){
        work = std::max(work, clusterBytes(cnf,pls,op,qdim,s.first,s.second));
        largest = std::max(largest, (qubitInMatrix ? qdim : 1)*s.first);
    }
    // the run's own buffers of a density matrix per time point (quantity dm: the product with and without the division,
    // the zeroth cluster and a cluster's values)
    if (qubitInMatrix && !strcasecmp(Config_getQuantity(cnf),"dm"))
        work += 4.0L*Config_getNstep(cnf)*16.0L*qdim*qdim;
    const double now = statusKb("VmRSS")*1024.0;
    const double estimate = now + double(std::min(work,(long double)1e30));

    // what one rank of this node may use
    MPI_Comm node;
    MPI_Comm_split_type(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&node);
    int nodeRanks = 1, nodeRank = 0;
    MPI_Comm_size(node,&nodeRanks); MPI_Comm_rank(node,&nodeRank);
    double available = 0;
    std::string source;
    const double perNode = slurmMegabytes("SLURM_MEM_PER_NODE"), perCpu = slurmMegabytes("SLURM_MEM_PER_CPU");
    if (perNode > 0){
        available = perNode*1048576.0/nodeRanks;
        source = "SLURM_MEM_PER_NODE";
    } else if (perCpu > 0){
        const double cpus = envNumber("SLURM_CPUS_ON_NODE");
        available = perCpu*1048576.0*(cpus > 0 ? cpus : nodeRanks)/nodeRanks;
        source = "SLURM_MEM_PER_CPU";
    } else {
        double free = nodeRank == 0 ? memAvailableBytes() : 0;     // one reading for the node
        MPI_Bcast(&free,1,MPI_DOUBLE,0,node);
        // MemAvailable leaves out what the ranks hold already, so each rank's share is its part of it and what it holds
        available = free > 0 ? free/nodeRanks + now : 0;
        source = "MemAvailable";
    }
    MPI_Comm_free(&node);

    // the rank whose estimate comes closest to (or passes) what it may use; a rank that cannot tell (0) never refuses
    struct { double ratio; int rank; } mine = {available > 0 ? estimate/available : 0.0, rank}, worst;
    MPI_Allreduce(&mine,&worst,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    double numbers[5] = {estimate, available, largest, double(nodeRanks), now};
    MPI_Bcast(numbers,5,MPI_DOUBLE,worst.rank,MPI_COMM_WORLD);
    double largestAll = largest;
    MPI_Allreduce(&largest,&largestAll,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    const double GB = 1e9;
    if (rank == 0){
        if (numbers[1] > 0)
            printf("      Memory estimate    : %.3g GB per rank (rank %d: %.3g GB held + its largest cluster, dimension %.0f), %.3g GB "
                   "available per rank (%s / %d rank%s on the node)\n", numbers[0]/GB, worst.rank, numbers[4]/GB, numbers[2],
                   numbers[1]/GB, source.c_str(), int(numbers[3]), numbers[3] > 1 ? "s" : "");
        else
            printf("      Memory estimate    : %.3g GB per rank (largest cluster dimension %.0f), available memory unknown (not checked)\n",
                   numbers[0]/GB, largestAll);
        fflush(stdout);
    }
    if (worst.ratio > 1.0){
        const bool dm = !strcasecmp(Config_getQuantity(cnf),"dm");
        const bool full = !strcasecmp(Output_getSavemode(op),"cluster_density_matrix");
        ccexRefusalPrint(stderr,"Error: not enough memory: the largest cluster (dimension %.0f) needs about %.3g GB per rank, and %.3g GB "
                         "is available per rank (%s, %d rank%s on the node). Run fewer ranks per node, a smaller order or smaller "
                         "clusters (rdip, cluster.pcce or cluster.bundle max_dim)%s%s\n", numbers[2], numbers[0]/GB, numbers[1]/GB,
                         source.c_str(), int(numbers[3]), numbers[3] > 1 ? "s" : "",
                         dm ? ", or quantity coherence instead of dm" : "",
                         full ? ", or output.savemode coherence instead of cluster_density_matrix" : "");
        ccexRefuse();
    }
}
