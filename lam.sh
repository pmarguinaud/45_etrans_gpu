#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9

mpirun -np 4 ./bin/AATESTPROG --namelist fort.4.20x20 --field-file 20x20/AATESTPROG.20x20.gp.000055.dat --time 1  > AATESTPROG.eo 2>&1


rm -f snapshot_*.png

~/3d/glgrib/glgrib.sh --field[0].path AATESTPROG.fa%ZZZFFFFF --field[0].palette.name cold_hot --scene.center.on --colorbar.on 
mv snapshot_0000.png ~/tmp/.
ssh ecgate scp ~/tmp/snapshot_0000.png phi001@90.76.140.145:tmp/snapshot_0000.png

