import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
import sys, os
from MDAnalysis.analysis import rms

path = ['mol6/trial3', 'mol6/flip', 'mol7/trial1', 'mol7/flip', 'mol13/trial1', 'mol13/flip']
lig = ['6', '6',  '7', '7', '13', '13']
resid = ['XRAY', 'FLIP', 'XRAY', 'FLIP', 'XRAY', 'FLIP']
new_name = ['1', '1', '2', '2', '3', '3']

# Create figure and 6 subplots
fig, ax = plt.subplots(nrows=3, ncols=2, figsize=(10, 6))

for i in range(len(path)):
    u = mda.Universe(f'{path[i]}/5ag5_prot-{lig[i]}-neutralized.psf', f'{path[i]}/5ag5_prot-{lig[i]}-prod.dcd', format='DCD', in_memory=False)
    ref = mda.Universe(f'{path[i]}/5ag5_prot-{lig[i]}-neutralized.psf', f'{path[i]}/5ag5_prot-{lig[i]}-finalmini.crd', format='CRD')

    protein = u.select_atoms("protein")

    for segids in range(len(u.segments.segids)):
        segid = u.segments.segids[segids]
        if segid == 'XRAY' or segid == 'REVR':
            #out = segid
            LIG = f'resname {segid}'

    print(LIG)
    
    R = rms.RMSD(u,  # universe to align
                 u,  # reference universe or atomgroup
                 select='backbone',  # group to superimpose and calculate RMSD
                 groupselections=[LIG],  # groups for RMSD
                 ref_frame=0)  # frame index of the reference
    R.run()
    
    df = pd.DataFrame(R.results.rmsd, columns=['Frame', 'Time (ns)', 'Backbone', f'{resid[i].lower()}'])

    if (i%2) == 0:
        k = 0
    else:
        k = 1
    if i > 1:
        j = int( i - i / 2 )
    else:
        j = 0
    
    df.plot(x='Frame', y=['Backbone', f'{resid[i].lower()}'], kind='line', ax=ax[j,k])
    ax[j,k].legend(title=f'compound{new_name[i]}')


#--------------------------------
#----- hide axis tick labels----
ax[0,0].set_xlabel('')
ax[0,1].set_xlabel('')
ax[1,0].set_xlabel('')
ax[1,1].set_xlabel('')
ax[0,0].xaxis.set_ticklabels([])
ax[0,1].xaxis.set_ticklabels([])
ax[1,0].xaxis.set_ticklabels([])
ax[1,1].xaxis.set_ticklabels([])
#-------------------------------
ax[1, 0].set_ylabel(r'RMSD ($\AA$)', fontsize=18)
ax[2, 0].set_xlabel('Simulation Frames', fontsize=18)
ax[2, 1].set_xlabel('Simulation Frames', fontsize=18)

##---------------------------
for p in range(3):
    for q in range(2):
        # Customize the linewidth of the x and y axes
        #ax[p,q] = plt.gca()
        ax[p,q].spines['bottom'].set_linewidth(2)  # x-axis
        ax[p,q].spines['left'].set_linewidth(2)    # y-axis
        ax[p,q].spines['top'].set_linewidth(2)  # x-axis
        ax[p,q].spines['right'].set_linewidth(2)    # y-axis

        # Customize the tick sizes
        ax[p,q].tick_params(axis='both', which='major', labelsize=16, width=2, length=8)  # Major ticks

##--------------------------

#plt.legend(fontsize="10")

# changing the fontsize of ticks
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16)

plt.savefig('lmnmt_opls.png', dpi=600)
plt.savefig('lmnmt_opls.pdf', dpi=600)
plt.savefig('lmnmt_opls.svg', dpi=600)
plt.savefig('lmnmt_opls.tif', dpi=300)

plt.show()

quit()
