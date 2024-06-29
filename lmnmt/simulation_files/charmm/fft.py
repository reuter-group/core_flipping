import sys
import numpy as np

"""
This is extracted from Charlie's pyCHARMM script for protein dynamics
Date: 10 October 2022
"""

# Ensure that FFT grid is product of small primes 2, 3, 5
def is_factor(n):
    if (n % 2 != 0): return False  # favors even number
    while n:
        flag = False
        for x in (2,3,5):
            if n % x == 0:
               n = n / x
               flag = True
               break

        if flag: continue
        break

    if n == 1: return True
    return False

def checkfft(n, margin = 5):
    n = int(n) + margin
    while 1:
        if is_factor(n): break
        else: n += 1
    return n


boxsize = sys.argv[1]

boxhalf = float(boxsize)/2

fft = checkfft(n=np.ceil(boxhalf)*2,margin=0)
print(fft)
