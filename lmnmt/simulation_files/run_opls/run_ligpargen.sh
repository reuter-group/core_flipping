#!/bin/bash

## LigParGen run on one of the PC-lab systems
## pga043@129.177.190.61 
## conda activate boss
export BOSSdir=/export/softwares/BOSS/boss

for i in 7
do

echo "$i"

## cgen = CM1A-LBCC working with BOSS5.0

ligpargen -i "$i".pdb -cgen CM1A-LBCC -c 0 -cgenb CM1A-LBCC -cb 0

rm *.z *.q.* *gmx* *tinker* *openmm* *desmond* *xplor* *lammps*

mv "$i".charmm.rtf "$i".rtf
mv "$i".charmm.prm "$i".prm

python opls2charmm.py "$i"

done

## ligpargen: http://zarbi.chem.yale.edu/ligpargen/index.html

## charge model:  1.14*CM1A-LBCC


