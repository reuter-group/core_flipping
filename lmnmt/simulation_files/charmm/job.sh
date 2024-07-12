#!/bin/bash

export charmm=/net/orinoco/pga043/CHARMM_47a2/charmm_47a2/build_openmm_blade/charmm 
#export charmm=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

export mol=7
export seg=REVR

export dir=5ag5
export run="$dir"/mol"$mol"/flip
export lig="$mol"
export prot=5ag5_prot

#$charmm dir=$dir lig=$lig seg=$seg run=$run -i build_ligand.inp > $run/build_ligand.out
#$charmm dir=$dir run=$run lig=$lig prot=$prot seg=$seg -i complex.inp > $run/complex.out
#$charmm dir=$dir run=$run lig=$lig prot=$prot -i  init-mini.inp > $run/init-mini.out
#$charmm dir=$dir run=$run lig=$lig prot=$prot -i  solvate.inp > $run/solvate.out

export box=`grep 'GREATERVALUE' $run/solvate.out | head -n 1 | awk '{print $4}' | sed 's/^"\(.*\)"$/\1/'`
echo $box
export fft=`python fft.py $box`
echo $fft
export nions=`python get_charge.py $run/solvate.out`

#$charmm dir=$dir run=$run lig=$lig prot=$prot nions=$nions box=$box fft=$fft -i neutralize.inp > $run/neutralize.out
#$charmm dir=$dir run=$run lig=$lig prot=$prot box=$box fft=$fft -i final-mini.inp > $run/final-mini.out

#echo 'set iseed = ' `date +%H%M%S` > $run/seed.stream

#mpirun -np 10 $charmm dir=$dir run=$run lig=$lig prot=$prot box=$box fft=$fft iseed=$RANDOM -i heat.inp > $run/heat.out 

#$charmm dir=$dir run=$run lig=$lig prot=$prot box=$box fft=$fft -i charmm_equil.inp > $run/charmm_equil.out 2> error
#$charmm dir=$dir run=$run lig=$lig prot=$prot box=$box fft=$fft -i charmm_prod.inp > $run/charmm_prod.out 2> error

#---------------------------------------------------
#---------------------------------------------------
#$charmm dir=$dir run=$run lig=$lig prot=$prot box=$box fft=$fft -i orient.inp #> $run/orient.out

