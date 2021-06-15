#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9


if [ 0 -eq 1 ]
then

export DR_HOOK_NOT_MPI=1
export PGI_ACC_DEBUG=1
mpirun -np 1 ./bin/AATESTPROG --no-write \
  --namelist fort.4.100x100 \
  --time 10 --check # > AATESTPROG.eo 2>&1

exit

fi

#--alloperm \
mpirun -np 2 ./bin/AATESTPROG --no-write \
  --namelist fort.4.100x100 \
  --time 10 --check # > AATESTPROG.eo 2>&1

exit

for f in *.fa
do
  cp $f ~/tmp/
  ssh ecgate scp ~/tmp/$f phi001@90.76.140.145:tmp/$f
done

