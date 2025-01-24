import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis.distances import dist
import sys

# usgae: python get_residues.py noe_xray.str
#	 python get_residues.py noe_revr.str 

u = mda.Universe('minimized.psf', 'minimized.crd')

f = open(f'{sys.argv[1]}', 'r')
out_name = sys.argv[1][-8:-4]
anchors_dict = {}
for lines in f.readlines():
    line = lines.split()
    try:
        if line[0] == 'assign':
            prot_atom = u.select_atoms(f'segid {line[3]} and resid {line[4]} and name {line[5]}')
            lig_atom = u.select_atoms(f'segid {line[9]} and resid {line[10]} and name {line[11]}')
            print(f'{prot_atom.resnames[0]}{prot_atom.resids[0]}-{prot_atom.names[0]}:{lig_atom.names[0]}')
    except IndexError: None 
