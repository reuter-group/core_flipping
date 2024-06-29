import sys, os
import numpy as np

"""
check here for the equation: https://en.wikipedia.org/wiki/Free_energy_perturbation

For going in forward direction,
dE =  E_i+1 - E_i ; both energies calculated from trajectory of window i using lambdas i and i+i

For going in backward direction,
dE =  E_i-1 - E_i ; both energies calculated from trajectory of window i using lambdas i and i-i

"""

kB=0.0019872041 # unit: kcal/(mol*K)
T=298           # unit: K

def perturbation(f1, f0, f2): #f1 = win-1, f0 = win, f2 = win+1
    frwd = f2 - f0
    exp = np.exp(-frwd/(float(kB)*float(T)))
    average = np.average(exp)
    deltaF_frwd = -float(kB)*float(T) * np.log(average)

    bcwd = f1 - f0
    exp = np.exp(-bcwd/(float(kB)*float(T)))
    average = np.average(exp)
    deltaF_bcwd = -float(kB)*float(T) * np.log(average)

    return deltaF_frwd, deltaF_bcwd 

#-----------------------------------------------------------------
#------------------------------------------------------------------
win = 1
#out = 'forxray'
out = 'forflip'

forw = []
back = []
dGs  = []
for myrep in range(3):
    f1 = np.genfromtxt(f'win{win}/{out}/win{win-1}_{myrep}.dat', dtype=None, skip_header=1, usecols=-1, delimiter=' ')
    f0 = np.genfromtxt(f'win{win}/{out}/win{win}_{myrep}.dat', dtype=None, skip_header=1, usecols=-1, delimiter=' ')
    f2 = np.genfromtxt(f'win{win}/{out}/win{win+1}_{myrep}.dat', dtype=None, skip_header=1, usecols=-1, delimiter=' ')
    
    forw.append(perturbation(f1, f0, f2)[0])
    back.append(perturbation(f1, f0, f2)[-1])
    dGs.append(-perturbation(f1, f0, f2)[-1] + perturbation(f1, f0, f2)[0])

final_dG = np.average(dGs)
std      = np.std(dGs)

print(f'Free energy of removing the restraints:')
print(f'{round(final_dG, 2)} +/- {round(std, 2)}')
quit()
