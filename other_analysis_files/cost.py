import numpy as np
import math
import sys, os
import subprocess

#equil= 81
equil = int(sys.argv[1])

print(os.system('head prep/minimized.psf | grep -i natom'))

cpus = [] # time in seconds
cpum = [] # time in minutes
cpuh = [] # time in hours

sim = []

for i in range(1, equil):
    f = open(f'run{i}/output_0', 'r')
    cpu_time = f.readlines()[-1]
    if cpu_time.split()[-1] == 'SECONDS':
        #print('time in seconds')
        time = float(cpu_time.split()[-2])
        cpus.append(time)
    elif cpu_time.split()[-1] == 'MINUTES':
        #print('time in minutes')
        time = float(cpu_time.split()[-2])
        cpum.append(time)
    else:
        print('time in hours')
        time = float(cpu_time.split()[-2])
        cpuh.append(time)

    sim_time = subprocess.check_output(f"grep 'DYNA>' run{i}/output_0 | tail -n 1 ", shell=True) 
    sim.append(float(sim_time.split()[1])*0.002/1000) # nstep*time_step/1000 (ns)

print ('ALF Flattening Statistics:')
print(f'total time in seconds = {np.sum(cpus)}')
print(f'total time in minutes = {np.sum(cpum)}')
print(f'total time in hours = {np.sum(cpuh)}')
print(f'total simulation time in ns = {np.sum(sim)}')
print(f'total time spent in hours = {round((np.sum(cpus)/3600) + (np.sum(cpum)/60) + (np.sum(cpuh)), 2)}')
print(f'============================')


cpus = [] # time in seconds
cpum = [] # time in minutes
cpuh = [] # time in hours

dupl='abcde'
for i in range(len(dupl)):
    for j in range(1, 5+1):
        f = open(f'run{equil}{dupl[i]}/output_{j}_0', 'r')
        cpu_time = f.readlines()[-1]
        if cpu_time.split()[-1] == 'SECONDS':
            #print('time in seconds')
            time = float(cpu_time.split()[-2])
            cpus.append(time)
        elif cpu_time.split()[-1] == 'MINUTES':
            #print('time in minutes')
            time = float(cpu_time.split()[-2])
            cpum.append(time)
        else:
            print('time in hours')
            time = float(cpu_time.split()[-2])
            cpuh.append(time)

print ('Equilibration Statistics:')
print(f'total time in seconds = {np.sum(cpus)}')
print(f'total time in minutes = {np.sum(cpum)}')
print(f'total time in hours = {np.sum(cpuh)}')
print(f'total time spent in hours = {round((np.sum(cpus)/3600) + (np.sum(cpum)/60) + (np.sum(cpuh)), 2)}')
print(f'============================')

if len(sys.argv) <=2:
    prod = f'{equil+1}'
else:
    prod = int(sys.argv[2])

cpus = [] # time in seconds
cpum = [] # time in minutes
cpuh = [] # time in hours

dupl='abcde'
for i in range(len(dupl)):
    for j in range(1, 30+1):
        f = open(f'run{prod}{dupl[i]}/output_{j}_0', 'r')
        cpu_time = f.readlines()[-1]
        if cpu_time.split()[-1] == 'SECONDS':
            #print('time in seconds')
            time = float(cpu_time.split()[-2])
            cpus.append(time)
        elif cpu_time.split()[-1] == 'MINUTES':
            #print('time in minutes')
            time = float(cpu_time.split()[-2])
            cpum.append(time)
        else:
            print('time in hours')
            time = float(cpu_time.split()[-2])
            cpuh.append(time)

print ('Production Statistics:')
print(f'total time in seconds = {np.sum(cpus)}')
print(f'total time in minutes = {np.sum(cpum)}')
print(f'total time in hours = {np.sum(cpuh)}')
print(f'total time spent in hours = {round((np.sum(cpus)/3600) + (np.sum(cpum)/60) + (np.sum(cpuh)), 2)}')
print(f'============================')

quit()
