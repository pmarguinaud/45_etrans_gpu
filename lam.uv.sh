#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9


mpirun -np 1 ./bin/AATESTPROGUV \
  --namelist fort.4.20x20 \
  --time 1 --no-write

exit

IP=86.201.29.85

for f in *.fa
do
  cp $f ~/tmp/
  ssh ecgate scp ~/tmp/$f phi001@$IP:tmp/$f
done

