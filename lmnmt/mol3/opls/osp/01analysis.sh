#!/bin/bash
#Analysis part requires python2 so execute conda activate py2

CHARMMEXEC=/net/orinoco/pga043/charmm/49a1_blade/charmm/build_charmm/charmm

for j in 0 1 2
do

$CHARMMEXEC dir=win1 myrep=$j -i rmsf_ligand.inp

done

