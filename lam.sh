#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9


if [ 0 -eq 1 ]
then

for nproc in 1 2
do

mpirun -np $nproc ./bin/AATESTPROGDER \
  --namelist fort.4.20x20 \
  --time 1 > AATESTPROG.eo 2>&1

mv AATESTPROG.eo AATESTPROG.eo.$nproc

for f in fort.???
do
  mv $f $f.$nproc
done

done

else

mpirun -np 2 ./bin/AATESTPROGDER \
  --namelist fort.4.20x20 \
  --time 1 > AATESTPROG.eo 2>&1

for f in *.fa
do
  cp $f ~/tmp/
  ssh ecgate scp ~/tmp/$f phi001@90.76.140.145:tmp/$f
done

fi
