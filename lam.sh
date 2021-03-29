#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9

n=000045
n=000001

mpirun -np 1 ./bin/AATESTPROG --namelist fort.4.20x20 --field-file 20x20/AATESTPROG.20x20.gp.$n.dat --time 1  > AATESTPROG.eo 2>&1

rm -f snapshot_*.png

for i in 1 2
do
~/3d/glgrib/glgrib.sh --field[0].path AATESTPROG.fa%SURFFFFF.000$i --field[0].palette.name cold_hot --scene.center.on --colorbar.on 
done

for i in 0 1
do
mv snapshot_000$i.png ~/tmp/.
ssh ecgate scp ~/tmp/snapshot_000$i.png phi001@90.76.140.145:tmp/snapshot_000$i.png
done

