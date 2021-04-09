#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9


mpirun -np 2 ./bin/AATESTPROG --no-write --alloperm \
  --namelist fort.4.1000x1000 \
  --time 10 # > AATESTPROG.eo 2>&1

exit

for f in *.fa
do
  cp $f ~/tmp/
  ssh ecgate scp ~/tmp/$f phi001@90.76.140.145:tmp/$f
done

