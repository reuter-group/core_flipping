#!/bin/bash
#Analysis part requires python2 so execute conda activate py2

CHARMMEXEC=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

sed -i 's/kmin 40/kmin @fcon1/g' prep/noe_xray.str
sed -i 's/kmax 40/kmax @fcon1/g' prep/noe_xray.str
sed -i 's/kmin 40/kmin @fcon2/g' prep/noe_revr.str
sed -i 's/kmax 40/kmax @fcon2/g' prep/noe_revr.str


DIR=`pwd`

RUNDIR=$DIR/win$i

### Run the simulation
mkdir $RUNDIR
mkdir $RUNDIR/res $RUNDIR/dcd
cp -r variables.inp prep $RUNDIR/
cd $RUNDIR

### timeout -s SIGINT 8h
echo "window$i started"
for j in 0 1 2
do
echo "replica"$j" started"
mpirun -np 1 -x OMP_NUM_THREADS=1 --bind-to none --map-by node $CHARMMEXEC iflat=$i esteps=500000 nsteps=2500000 seed=$RANDOM fmax=$i myrep=$j -i ../msld_flat.inp > output 2> error

done

rm -rf res prep


