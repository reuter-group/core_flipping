import numpy as np
import matplotlib.pyplot as plt
import sys

'''
Written by: Parveen Gartan
14 April 2023
'''

#x = np.arange(1, int(sys.argv[1]), 1)
#print(x)

f = np.genfromtxt('fixed.dat', dtype=None, skip_header=0)
#print(f)
c = np.genfromtxt('c.dat', dtype=None, skip_header=0)
s1 = np.genfromtxt('s1.dat', dtype=None, skip_header=0)
s2 = np.genfromtxt('s2.dat', dtype=None, skip_header=0)
x1 = np.genfromtxt('x1.dat', dtype=None, skip_header=0)
x2 = np.genfromtxt('x2.dat', dtype=None, skip_header=0)

xray = np.genfromtxt('xray.dat', dtype=None, skip_header=0)
flip = np.genfromtxt('revr.dat', dtype=None, skip_header=0)
x = np.arange(0, len(xray), 1)

trans = np.genfromtxt('trans.dat', dtype=None, skip_header=0)

# Create figure and 6 subplots
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10, 6))

ax[0,0].plot(f, label="Block 3")
ax[0, 0].set_title('Fixed Bias')
ax[0,0].legend(frameon=False)

ax[0,1].plot(c)
ax[0, 1].set_title('Quadratic Bias')
#ax[0,1].legend(frameon=False)

ax[0,2].plot(s1, label='Block2')
ax[0,2].plot(s2, label='Block3')
ax[0, 2].set_title('End-point Bias')
ax[0,2].legend(frameon=False)

ax[1,0].plot(x1, label='Block2')
ax[1,0].plot(x2, label='Block3')
ax[1, 0].set_title('Skew Bias')
ax[1,0].legend(frameon=False)

ax[1,1].bar(x, xray, label='Block2', alpha=0.5)
ax[1,1].bar(x, flip, label='Block3', alpha=0.5)
ax[1, 1].set_title('Single Population')
ax[1,1].legend(frameon=False)

ax[1,2].plot(trans)
ax[1, 2].set_title('Single Transitions')
#ax[1,2].legend(frameon=False)

ax[1,0].set_xlabel("#runs")
ax[1,1].set_xlabel("#runs")
ax[1,2].set_xlabel("#runs")

# Set overall title and tight layout
fig.suptitle("ALF flattening")
fig.tight_layout()

# Show the plot
#plt.savefig('ALF_iterations.png', dpi=600)
plt.show()

