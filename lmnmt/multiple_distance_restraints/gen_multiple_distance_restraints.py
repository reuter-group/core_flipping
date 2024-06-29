import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
from numpy.linalg import norm
import matplotlib.pyplot as plt
import os, sys
from math import pi, sin, log
import scipy.integrate
import time

##----------------------
if len(sys.argv) !=3:
   print('Three command line options are required.')
   print('Usage: python gen_multiple_distance_restraints.py psf_file dcd_file')
   print('Exiting...')
   quit()
##-----------------------

# FIXME: Replace with paths to your own trajectory and topology files
u = mda.Universe(sys.argv[1], sys.argv[2], in_memory=False)

# number of frames
print(f'Frames in the trajectory: {u.trajectory.n_frames} ')

plots = False # do you want to generate plots

lig_heavy = u.select_atoms(" not (segid PROT or resname TIP3 or resname COF or resname POT \
                            or resname SOD or resname CLA or resname CL or name H*  or name LP*)")


out_name = lig_heavy[0].resname
#print(lig_heavy)
#quit()

# anchors dict of dict. For each ligand heavy atom there is a dictionary of protein heavy atoms,
# for each of which there is a dictionary of average distance and standard deviation

anchors_dict = {}
for lig_atom in lig_heavy:
    for prot_atom in u.select_atoms(f"(protein or segid PROT) and (around 15 index {lig_atom.index}) and (not around 10 index {lig_atom.index}) \
                                       and (not name H*)"): # protein does not recognise PRT
        anchors_dict[(lig_atom.index,prot_atom.index)]={}
        anchors_dict[(lig_atom.index, prot_atom.index)]["dists"]=[]


for frame in u.trajectory:
    for lig_atom_index, prot_atom_index in anchors_dict.keys():
        distance = dist(mda.AtomGroup([u.atoms[lig_atom_index]]), mda.AtomGroup([u.atoms[prot_atom_index]]), box=frame.dimensions)[2][0]
        anchors_dict[(lig_atom_index,prot_atom_index)]["dists"].append(distance)


# change lists to numpy arrays
for pair in anchors_dict.keys():
    anchors_dict[pair]["dists"] = np.array(anchors_dict[pair]["dists"])

# calculate average and SD
for pair in anchors_dict.keys():
    anchors_dict[pair]["avg_dist"] = anchors_dict[pair]["dists"].mean()
    anchors_dict[pair]["sd_dist"] = anchors_dict[pair]["dists"].std()

# get n pairs with lowest SD
pairs_ordered_sd=[]
for item in sorted(anchors_dict.items(), key=lambda item: item[1]["sd_dist"]):
    pairs_ordered_sd.append(item[0])
    #print(f'Pair: {item[0]}, av distance: {item[1]["avg_dist"]:.2f}, SD: {item[1]["sd_dist"]:.2f}')

# print out indices of top 20 ligand atoms
lig_anchors = []

for i in range(20):
    print(pairs_ordered_sd[i][0], end=' ')
    lig_anchors.append(pairs_ordered_sd[i][0])
    
print(f"\n\nThe number of unique ligand anchors is: {len(set(lig_anchors))}")


## Ensure no Duplicated Atoms in Restraints
unique_pairs_ordered_sd = []
lig_ats = []
recept_ats = []

for pair in pairs_ordered_sd:
    if pair[0] not in lig_ats:
        if pair[1] not in recept_ats:
            unique_pairs_ordered_sd.append(pair)
            lig_ats.append(pair[0])
            recept_ats.append(pair[1])
            if len(unique_pairs_ordered_sd) == 22:
                break



##------------ plots ---------------------
def get_distance(idx1, idx2, u):
    """ Distance in Angstrom"""
    distance = dist(mda.AtomGroup([u.atoms[idx1]]), mda.AtomGroup([u.atoms[idx2]]), box=u.dimensions)[2][0]
    return distance


num_pairs = round(len(unique_pairs_ordered_sd))

if plots == True:

   frames = u.trajectory.n_frames
   fig, ax = plt.subplots(1,1, figsize=(12,6))
   for i, pair in enumerate(unique_pairs_ordered_sd):
       if i > -1:
          ax.plot([x for x in range(frames)], anchors_dict[pair]["dists"],label=f"Pair {pair}")
          ax.set_ylabel("r ($\AA$)")
          if i == num_pairs-1:
              break
       ax.set_xlabel("Frame No")
       ax.legend()
   fig.tight_layout()
   plt.savefig('distance_' + str(out_name) + '.pdf')

# Plot histograms

   rows = round(len(unique_pairs_ordered_sd) / 2)

   fig, axs = plt.subplots(rows, 2, figsize=(18,6*(len(unique_pairs_ordered_sd)%6)))
   axs = axs.flatten()
   for i, pair in enumerate(unique_pairs_ordered_sd):
       axs[i].hist(anchors_dict[pair]["dists"],label=f"Pair {pair}", edgecolor='k')
       axs[i].axvline(anchors_dict[pair]["avg_dist"], color='r', linestyle = "dashed", linewidth=2, label=f'Mean: {anchors_dict[pair]["avg_dist"]:.2f} A\n Std: {anchors_dict[pair]["sd_dist"]:.2f} A')
       axs[i].set_xlabel("r ($\AA$)")
       axs[i].set_ylabel("Count")
       axs[i].legend()
   fig.tight_layout()
   plt.savefig('histogram_' + str(out_name) + '.pdf')

else:
    print('No plots!')
    None
##-----------------------------------------

##---------- Free Energy calculation functions -----------------------------
'''
NOE
       /  0.5*KMIN*(RAVE-RMIN)**2    R<RMIN
      /
     /    0.0                        RMIN<R<RMAX
E(R)=
     \    0.5*KMAX*(RAVE-RMAX)**2    RMAX<RAVE<RLIM
      \
       \  FMAX*(RAVE-(RLIM+RMAX)/2)  RAVE>RLIM


FMAX = 0.0
RSWITch - not using
SEXP - not using
RAVE = R ; no TCON specified in input
NOE inputs:
kmin = kmax = kr
rmin = r0 - r_fb
rmax = r0 + r_fb

r_fb = 0.4 (let's start with this)

original reference:
r0 =  (rmin + rmax ) / 2
r_fb =  rmax - r0 ## The radius of the flat-bottomed region was selected to be as small as possible
kr = kmin or kmax based on r
'''

# Constants
v0 = 1661 # A^3, the standard state volume
T = 298.15 # K
R = 0.0019872041 ## it should be Boltzmann's constant #unit: kcal/(mol*K), wrong => # kcal mol-1, the molar gas constant


def numerical_distance_integrand(r, r0, r_fb, kr):
    """Integrand for harmonic distance restraint. Domain is on [0, infinity], 
    but this will be truncated to [0, 8 RT] for practicality.

    Args:
        r (float): Distance to be integrated, in Angstrom 
        r0 (float): Equilibrium distance, in Angstrom
        r_fb (float): Flat-bottomed radius, in Angstrom
        kr (float): Force constant, in kcal mol-1 A-2

    Returns:
        float: Value of integrand
    """
    r_eff = abs(r - r0) - r_fb
    if r_eff < 0:
        r_eff = 0
    return (r**2)*np.exp(-(kr*r_eff**2)/(2*R*T))

def get_correction(r0, r_fb, kr):
    """Get the free energy of releasing the harmonic distance restraint.
    Domain is on [0, infinity], but this will be truncated to [0, 8 RT] for practicality.
    Args:
        r0 (float): Equilibrium distance, in Angstrom
        r_fb (float): Flat-bottomed radius, in Angstrom
        kr (float): Force constant, in kcal mol-1 A-2
    Returns:
        float: Free energy of releasing the restraint
    """
    dist_at_8RT = 4*np.sqrt((R*T)/kr) + r_fb #fb_radius # Dist. which gives restraint energy = 8 RT
    r_min = max(0, r0-dist_at_8RT)
    r_max = r0 + dist_at_8RT
    integrand = lambda r: numerical_distance_integrand(r, r0, r_fb, kr)
    z_r = scipy.integrate.quad(integrand, r_min, r_max)[0]
    dg = -R*T*log(v0/(4*np.pi*z_r))

    return dg

##---------------------------------------------------------------------

# Select flat-bottomed region to contain 95 % of the probability density from test simulation, then set high force constants
# for the half-harmonic potentials
# Form of restraints: restraints = { (i0, i1): (r01, kl, Dl), (i0,i2): (r02, kl, Dl) }

restraints = {}
## Atom indices starts from 0 in mdAnalysis
DG = []

try:
   os.remove('noe_'+str(out_name).lower() +'.str')
   os.remove('junk_'+str(out_name).lower() +'.log')
except FileNotFoundError:
    print('File does not exists.')




with open('noe_'+str(out_name).lower() +'.str', 'a') as f, open('junk_'+str(out_name).lower() +'.log', 'a') as f1:
     f1.write(f' prot.segid prot.resid at.name - lig.segid lig.resid at.name : avg.distance +/- sd \n')
     f.write(f'NOE \n')
     for pair in unique_pairs_ordered_sd:
         r0 = round(anchors_dict[pair]["avg_dist"],2) # A
         sd = round(anchors_dict[pair]["sd_dist"],2)
         kl = 40 # kcal/(mol * A**2)
         dl = 0.4
         restraints[pair] = (r0, kl, dl)
         #print(pair[0], pair[1], restraints[pair][0])
         #print(f"{u.atoms[pair[0]]} and {u.atoms[pair[1]]}")
         lig = u.atoms[pair[0]]
         prot = u.atoms[pair[1]]
        
         f1.write(f'{prot.segid} {prot.resid} {prot.name} - {lig.segid} {lig.resid} {lig.name} : {r0} +/- {sd} \n')

         f.write(f'   assign sele atom {prot.segid} {prot.resid} {prot.name} end sele atom {lig.segid} {lig.resid} {lig.name} end - \n')
         f.write(f'   kmin {kl} rmin {round(r0-dl,1)} kmax {kl} rmax {round(r0+dl,1)} fmax 2.0\n')
         #f.write(' \n')
         cor = get_correction(r0, dl, kl)
         DG.append(cor)
         #print(f"Free energy of releasing flat-bottomed restraint: {cor:.2f} kcal mol-1")
     f.write('\n')
     f.write(f'   print \n   print anal \n') 
     f.write(f'END \n')
     f.write(f'! Total Free energy of releasing all the restraints: {sum(DG):.2f} kcal/mol \n')
 
print(f'Total Free energy of releasing all the restraints: {sum(DG):.2f} kcal/mol')

quit()
