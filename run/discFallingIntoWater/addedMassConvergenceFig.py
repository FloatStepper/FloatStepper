#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

plt.rcParams.update({
	"text.usetex": True,
	"font.family": "Times New Roman",
'font.size': 14
})

fig, ax = plt.subplots(1)
#fig.suptitle('Added mass evolution for disc falling into water')

#axins = ax.inset_axes([.55, 0.05, .4, .5])
axins = ax.inset_axes([.6, 0.08, .36, .5])

scanName = "meshScan"
os.chdir(scanName)
caseNames = ["nx50ny50", "nx100ny100", "nx200ny200"]
nxVals = ["nx = 50", "nx = 100", "nx = 200"]
for count, nxVal in enumerate(nxLevs):
    caseDir=scanName+'/'+caseNames[count]
    os.system("./extractAddedMass.sh "+caseDir)
    os.chdir(caseDir)
    t = np.loadtxt('tim', unpack=True)
    A22 = np.loadtxt('A22', unpack=True, usecols = (0))
    nEl = min(len(t),len(A22))
    t = t[0:nEl]
    A22 = A22[0:nEl]

    #Plotting position
    line, = plt.plot(t, A22, '-', label=nxVals[count])
    line, = axins.plot(t, A22, '-')
    os.chdir("..")

plt.legend()
axins.set_xlim(0.225, 0.275)
axins.set_ylim(0, 400)
#axins.set_xticklabels([])
#axins.set_yticklabels([])
axins.grid(True)
ax.indicate_inset_zoom(axins, edgecolor="black")
plt.xlabel(r'$t$ [s]')
plt.ylabel(r'$A_{22}$ [kg]')
plt.grid(True)
plt.show()
#plt.savefig('fallingCircleAddedMassConvergence.pdf', bbox_inches='tight')