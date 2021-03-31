#!/bin/bash

ulimit -s unlimited
export OMP_NUM_THREADS=1

set -x
set -e

module unload gnu
module load nvhpc/20.9

n=000001
n=000045

mpirun -np 4 ./bin/AATESTPROG --namelist fort.4.20x20 --field-file 20x20/AATESTPROG.20x20.gp.$n.dat --time 10  > AATESTPROG.eo 2>&1



rm -f snapshot_*.png


ff=$(mpirun -np 1 ./bin/lfitools lfilist AATESTPROG.fa 2>/dev/null | grep SURFFF | perl -pe ' $x = eval $_; $_ = "$x->[0] "')

for f in $ff
do
~/3d/glgrib/glgrib.sh --field[0].path AATESTPROG.fa%$f --field[0].palette.name cold_hot --scene.center.on --colorbar.on --render.width 1000
done

for f in snapshot_*.png
do
mv $f ~/tmp/.
ssh ecgate scp ~/tmp/$f phi001@90.76.140.145:tmp/$f
done

