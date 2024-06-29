import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import sys

ref = mda.Universe(sys.argv[1], sys.argv[2])  ## XRAY psf and prod crd

mobile = mda.Universe(sys.argv[3], sys.argv[4])  ## REVR psf and prod crd   # we use the first frame

align.alignto(mobile, ref, select="(protein or segid PROT) and name CA", weights="mass")

lig = mobile.select_atoms('segid REVR')

lig.atoms.write('7r.crd')
lig.atoms.write('7r.pdb')

quit()

