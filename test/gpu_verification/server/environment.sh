#!/usr/bin/env bash
# Environment validated on ssh gpu (idm4-gpu); source from the Slurm jobs.
if ! type module >/dev/null 2>&1; then
    source /etc/profile.d/modules.sh
fi
module load gcc/gcc-11.5.0 gcc/openmpi-4.1.5 cuda/cuda-12.8 python/python3.13.10
# do_compile.sh gpu builds the C++ with icpx when oneAPI is installed, so the binaries need
# its runtime libraries; setvars.sh also puts Intel MPI first on PATH, which is undone here
# (the binaries are built and run with the module's Open MPI).
# Skipped when the environment already has it: setvars.sh then returns 3, which ends a job
# running under `set -e` before it starts.
if [ -r /opt/intel/oneapi/setvars.sh ] && [ -z "${SETVARS_COMPLETED:-}" ]; then
    ompi_bin="$(dirname "$(command -v mpirun)")"
    restore_u=$(set +o | grep nounset); set +u
    . /opt/intel/oneapi/setvars.sh > /dev/null 2>&1
    eval "$restore_u"
    export PATH="$ompi_bin:$PATH"
fi
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

