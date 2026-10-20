#!/usr/bin/env bash
# Environment validated on ssh gpu (idm4-gpu); source from the Slurm jobs.
# The jobs run with --export=NONE (idm's login modules, intel/anaconda, are not passed to the GPU node; 2026-10-13, idm4-gpu
# is the gpu partition of idm's Slurm), so the module function comes from modules.sh and the list starts empty.
if ! type module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh
fi
module purge
module load gcc/gcc-11.5.0 intel/oneapi-2023.2 cuda/cuda-13.2 python/python3.13.10
# intel/oneapi-2023.2: what do_compile.sh gpu builds with and the binaries run on (icpx's runtime
# libraries, Intel MPI's mpirun, MKL), so the launcher and the MPI of the binaries are one. Until
# 2026-10-20 the binaries were Open MPI's (gcc/openmpi-4.1.5) and setvars.sh added icpx; the checks,
# which take the Intel MPI of I_MPI_ROOT first (test/README.md), then needed CCEX_MPI.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

