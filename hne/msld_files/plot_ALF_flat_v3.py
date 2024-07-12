import numpy as np
import matplotlib.pyplot as plt
import sys

'''
Written by: Parveen Gartan
10 February 2024
updated: 11 June 2024
written for single site systems only
'''

end = int(sys.argv[1])

nsubs=np.loadtxt('nsubs',dtype='int',ndmin=1)
nblocks=np.sum(nsubs)

fix_pop_trans = {}

for block in range(nblocks):
    fix_pop_trans[f'block{block}'] = {}
    fix_pop_trans[f'block{block}']["fix"] = []
    fix_pop_trans[f'block{block}']["pop"] = []
    for var in np.arange(1, end+1, 1):
        try:
            alf_var = np.genfromtxt(f'variables{var}.inp', dtype=str)
            fix_pop_trans[f'block{block}']["fix"].append(float(alf_var[block][-1]))
            
        except FileNotFoundError:
            break
        except IOError:
            break

# get single transitions at a site
fix_pop_trans["trans"] = {}
fix_pop_trans["trans"]["trans"] = []

for var in np.arange(1, end+1, 1):
    dyna = open(f'run{var}/output_0', 'r')
    block = 0
    for lines in dyna:
        if lines.startswith('SINGLE TRANSITIONS>'):
            line = lines.split()
            fix_pop_trans["trans"]['trans'].append(line[-1])
        ## get population of each block
        if lines.startswith('SINGLE POPULATION>'):
            pops = int(lines.split()[-1])
            fix_pop_trans[f'block{block}']["pop"].append(pops) 
            block += 1

# Create figure and 6 subplots
fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(10, 6))

for block in range(nblocks):
        ax[0].plot(fix_pop_trans[f'block{block}']["fix"], label=f'block{block}')
        ax[0].set_title('Fixed Bias', fontsize=18)
        ax[0].legend(frameon=False)
        
        ax[1].plot(fix_pop_trans[f'block{block}']["pop"], label=f'block{block}')
        ax[1].set_title('Single Population', fontsize=18)
        ax[1].legend(frameon=False)

        ax[2].plot(fix_pop_trans["trans"]["trans"], label=f'block{block}')
        ax[2].set_title('Transitions', fontsize=18)
        ax[2].legend(frameon=False)


ax[0].set_xlabel("ALF iteration index", fontsize=18)
ax[1].set_xlabel("ALF iteration index", fontsize=18)
ax[2].set_xlabel("ALF iteration index", fontsize=18)


plt.show()
quit()
##---------------------------
for p in range(1):
    for q in range(3):
        # Customize the linewidth of the x and y axes
        #ax[p,q] = plt.gca()
        ax[p,q].spines['bottom'].set_linewidth(2)  # x-axis
        ax[p,q].spines['left'].set_linewidth(2)    # y-axis
        ax[p,q].spines['top'].set_linewidth(2)  # x-axis
        ax[p,q].spines['right'].set_linewidth(2)    # y-axis

        # Customize the tick sizes
        ax[p,q].tick_params(axis='both', which='major', labelsize=16, width=2, length=8)  # Major ticks

##--------------------------
quit()

#fig.legend(fontsize="10")

# changing the fontsize of ticks
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16)


# Set overall title and tight layout
#fig.suptitle("ALF flattening")
fig.tight_layout()

# Show the plot
#plt.savefig('ALF_iterations.png', dpi=600)
#plt.savefig('ALF_iterations.svg', dpi=600)
#plt.savefig('ALF_iterations.tif', dpi=600)

plt.show()

