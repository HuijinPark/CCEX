#!/bin/sh

file="pulse"
./do_compile.sh idm4 ${file}
mpirun -n 1 "${file}.out" | tee proc_${file}
