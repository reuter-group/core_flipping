import numpy as np
import matplotlib.pyplot as plt
import sys

'''
Written by: Parveen Gartan
10 February 2024
'''
end = int(sys.argv[1])

l   = []
c   = []
s1  = []
s2  = []
x1  = []
x2  = []

for var in np.arange(1, end+1, 1):
    try:
        alf_var = np.genfromtxt(f'../variables{var}.inp', dtype=str)
        lam1 = float(alf_var[0][-1])  ## phi1
        l.append(float(alf_var[1][-1]))  ## phi2
        c.append(float(alf_var[2][-1]))  ## psi1,2
        s1.append(float(alf_var[5][-1]))  ## omega1,2
        s2.append(float(alf_var[6][-1]))  ## omega2,1
        x1.append(float(alf_var[3][-1]))  ## chi1,2
        x2.append(float(alf_var[4][-1]))  ## chi2,1
    except FileNotFoundError:
        break
    except IOError:
        break

pop  = []
trans = []
for i in range(1, end+1):
    try:
        dyna = open(f'../run{i}/output_0', 'r')
        for lines in dyna:
            if lines.startswith('SINGLE TRANSITIONS>'):
                line = lines.split()
                trans.append(line[-1])
            if lines.startswith('SINGLE POPULATION>'):
                #print(lines.split()[-1])
                pop.append(int(lines.split()[-1]))
    except FileNotFoundError:
        break
    except IOError:
        break

xray = pop[::2]
flip = pop[1::2]

# Create figure and 6 subplots
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))

ax[0,0].plot(l, label="Block 3")
ax[0,0].set_title('Fixed Bias', fontsize=18)
ax[0,0].legend(frameon=False)

ax[0,1].plot(c)
ax[0,1].set_title('Quadratic Bias', fontsize=18)
#ax[0,1].legend(frameon=False)

ax[1,0].plot(s1, label='Block2')
ax[1,0].plot(s2, label='Block3')
ax[1,0].set_title('End-point Bias', fontsize=18)
ax[1,0].legend(frameon=False)

ax[1,1].plot(x1, label='Block2')
ax[1,1].plot(x2, label='Block3')
ax[1,1].set_title('Skew Bias', fontsize=18)
ax[1,1].legend(frameon=False)

ax[1,0].set_xlabel("ALF iteration index", fontsize=18)
ax[1,1].set_xlabel("ALF iteration index", fontsize=18)

##---------------------------
for p in range(2):
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

#fig.legend(fontsize="10")

# changing the fontsize of ticks
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16)


# Set overall title and tight layout
#fig.suptitle("ALF flattening")
fig.tight_layout()

# Show the plot
plt.savefig('ALF_iterations.png', dpi=600)
plt.savefig('ALF_iterations.svg', dpi=600)
plt.savefig('ALF_iterations.tif', dpi=600)

plt.show()

