#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9

mpirun -np 1 ./bin/AATESTPROG --namelist fort.4.t31 --time 1 > AATESTPROG.eo 2>&1


