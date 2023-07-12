#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os
import fileinput

z0 = 0
L = 1

def write_to_ascii(filename, tarray, parray):
    with open(filename, 'w') as file:
        nTimes = len(tarray)
        times = ' '.join(str(x) for x in tarray)
        file.write(f"times {nTimes} ({times});\n\n")
        pos = ' '.join(str(x) for x in parray)
        file.write(f"pistonPositions 1 ( {nTimes} ({pos}) );")

#Extracting and loading trajectory data from floaterFoam log file
t, ang = np.loadtxt('Flap.dat', unpack=True, usecols = (0, 1), skiprows=1)
nTimes = len(t)

pistonPosition = L*np.sin(ang*2*math.pi/360)

plt.plot(t, pistonPosition,'-k', linewidth=4)
#plt.plot(t, ang, '-r', linewidth=4)
plt.xlabel('t')
plt.ylabel('Displacement')
plt.title('Piston Position')

plt.show()

write_to_ascii('pistonPosition', t, pistonPosition)
