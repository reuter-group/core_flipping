import sys
import numpy as np

run = sys.argv[1]
prod = 30
ndupl = 5

nblocks=np.loadtxt('nblocks',dtype='int')
nsubs=np.loadtxt('nsubs',dtype='int',ndmin=1)
nreps=np.loadtxt('nreps',dtype='int')
dupl = 'abcdefghijkl'

fpl = np.zeros(shape=(ndupl, prod))
trans = np.zeros(shape=(ndupl, prod))

for i in range(ndupl):
    for j in range(prod):
        with open(f'run{run}{dupl[i]}/output_{j+1}_0', 'r') as f:
            for lines in f.readlines():
                try:
                    if lines.startswith('FRACTION PHYSICAL LIGAND>'):
                        #print(lines)
                        fpl[i,j] = float(lines.split()[-1])
                    if lines.startswith('SINGLE TRANSITIONS>'):
                        print(lines)
                        trans[i,j] = int(lines.split()[-1]) 
                except IndexError:
                    None


print(np.average(fpl, axis=1)) # row wise or average for each replica
print(np.sum(trans, axis=1)) # row wise or average for each replica
print(round(np.average(fpl), 3))
print(np.sum(trans))
#print(fpl)

#grep 'SINGLE TRANSITIONS>' run73a/output_*_0 | awk -F' ' '{sum+=$5;} END{print sum;}'

