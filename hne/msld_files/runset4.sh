#!/bin/bash

#GetSteps.py requires python2.7
#execute this before running subset4.sh = conda activate py2

CHARMMEXEC=/net/orinoco/pga043/CHARMM_47a2/charmm_47a2/build_openmm_blade/charmm 

ii=$run
ibeg=$(( ( $ii - 1 ) * $nitt + 1 ))
iend=$(( $ii * $nitt ))

DIR=`pwd`

name=`cat name`
nnodes=`cat nnodes`
nreps=`cat nreps`

RUNDIR=$DIR/run$i
echo "RUNDIR is run${i}"

# If this is the first iteration, set up the run directory
if [[ $ii -eq 1 ]]
then

mkdir $RUNDIR
mkdir $RUNDIR/res $RUNDIR/dcd $RUNDIR/failed

cp msld_prod.inp $RUNDIR/
cp variables$ini.inp $RUNDIR/variablesprod.inp
cp -r prep $RUNDIR

fi
# Done setting up the run directory

cd $RUNDIR

while [ $ibeg -gt 1 -a `../ALF/GetSteps.py res/${name}_prod$(( $ibeg - 1 )).lmd_0` -ne 50000 ]
do
echo "Error: Run $ibeg incomplete. Going back one step"
ibeg=$(( $ibeg - ( 1 * $nitt ) ))
done

for itt in `seq $ibeg $iend`
do

# Keep trying until simulation completes correctly with 50000 steps in lambda file
while [ `../ALF/GetSteps.py res/${name}_prod${itt}.lmd_0` -ne 50000 ]
do

# If run failed in past, file the failed files before moving on
if [ -f output_${itt} ]; then
  cp output_${itt} output_${itt}_* error_${itt}* failed/
  rm output_${itt} output_${itt}_* error_${itt}*
fi


# Run the simulation
#mpirun -np $(($nreps * $nnodes)) -x OMP_NUM_THREADS=4 --bind-to none --bynode 
mpirun -np $nreps -x OMP_NUM_THREADS=$nnodes --map-by node $CHARMMEXEC nsteps=500000 nsavc=10000 seed=$RANDOM itt=$itt -i msld_prod.inp > output_$itt 2> error_$itt

done

done
