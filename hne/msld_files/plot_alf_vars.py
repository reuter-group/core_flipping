import numpy as np
import matplotlib.pyplot as plt
import os, sys

win = np.arange(0, 1, 0.05, dtype=float)
#print(win)

## G= aλ + bλ(1−λ),  
## a = biasing potential parameters for the fixed term
## b = biasing potential parameters for the quadratic term
## x = lambda

def alf_ld(x, l, c, s1, s2, x1, x2):
    alpha = 0.017
    sigma = 0.18
    return l * x + c * x * (1-x) + (s1 * x * (1-x) / (x + alpha)) + (s2 * x * (1-x) / ((1 - x) + alpha)) + \
           (x1 * (1-x) * (1 - np.exp(-x/sigma))) + (x2 * x * (1 - np.exp(-(1-x)/sigma)))

'''
lam = phi   => fixed
c   = psi   => quadratic
s   = omega => end-point
x   = chi   => skew
'''

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6))

start = int(sys.argv[1])
end = int(sys.argv[2])

for var in np.arange(start, end, 1):
    alf_var = np.genfromtxt(f'variables{var}.inp', dtype=str)
    lam1 = float(alf_var[0][-1])  ## phi1
    lam2 = float(alf_var[1][-1])  ## phi2
    c    = float(alf_var[2][-1])  ## psi1,2
    s1   = float(alf_var[5][-1])  ## omega1,2
    s2   = float(alf_var[6][-1])  ## omega2,1
    x1   = float(alf_var[3][-1])  ## chi1,2
    x2   = float(alf_var[4][-1])  ## chi2,1
    print(f'{lam1}, {lam2}, {c}, {s1}, {s2}, {x1}, {x2}')
    free_energy_curve = alf_ld(win, lam2, c, s1, s2, x1, x2)
    #print(free_energy_curve)
    plt.plot(win, free_energy_curve, label=f'iteration{var}', linestyle='--')

##---------------------------
# changing the fontsize of ticks
plt.xticks(fontsize=16, rotation=0)
plt.yticks(fontsize=16)

# Customize the linewidth of the x and y axes
ax = plt.gca()
ax.spines['bottom'].set_linewidth(2)  # x-axis
ax.spines['left'].set_linewidth(2)    # y-axis
ax.spines['top'].set_linewidth(2)  # x-axis
ax.spines['right'].set_linewidth(2)    # y-axis
##--------------------------

#fig.legend(fontsize="10")
#---------------------------------

#plt.title('ALF free energy curve')
plt.xlabel(r'$\lambda$', fontsize=18)
plt.ylabel('Free Energy (kcal/mol)', fontsize=18)
plt.legend(fontsize="12", frameon=True, ncol=4)
#plt.tight_layout()

plt.savefig('test.svg', dpi=600)
plt.savefig('test.png', dpi=600)
plt.savefig('test.tif', dpi=600)


plt.show()



