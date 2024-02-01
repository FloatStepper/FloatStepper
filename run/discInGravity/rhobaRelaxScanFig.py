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

cwd = os.getcwd()
figdir = os.path.join(cwd, "../figures")
#os.mkdir(figdir)
figname="rhob_aRelaxScan.pdf"

def merge(list1, list2):
     
    merged_list = [(list1[i], list2[i]) for i in range(0, len(list1))]
    return merged_list
     
os.chdir("rhobaRelaxScan")

# Check which of the cases have a 3 folder, meaning they did not crash
tend='3'
os.system('ls -d */'+tend+' | cut -d\'y\' -f3 | cut -d\'a\' -f1 > rhob')
os.system('ls -d */'+tend+' | cut -d\"n\" -f4 | cut -d\"/\" -f1 > aRelax')

rhob = np.loadtxt('rhob', unpack=True)
rhobArray = sorted(set(rhob))
aRelax = np.loadtxt('aRelax', unpack=True)
aRelaxArray = sorted(set(aRelax))
rhobAll, aRelaxAll = np.meshgrid(rhobArray,aRelaxArray)

setA = set(merge(rhob, aRelax))
setB = set(merge(rhobAll.flatten(), aRelaxAll.flatten()))
BminusA = list(setB.difference(setA))
rhobBad = [i[0] for i in BminusA]
aRelaxBad = [i[1] for i in BminusA]

rhobs = np.linspace(0.1,1,50)
aLimit = 2/(1+1/rhobs)

#Plotting position
fig, ax = plt.subplots()
#fig.suptitle('Buoyant disk convergence at T = '+tend+'s')
ax.plot(rhob, aRelax,'+r',zorder=2,markersize=10)
ax.plot(rhobBad, aRelaxBad,'xb',zorder=0,markersize=7)
ax.plot(rhobs, aLimit,'-k',zorder=1)
ax.set_xlabel(r'Mass ratio, $\rho_{b}/\rho_{f} [-]$')
ax.set_ylabel(r'Acceleration relaxation, $\gamma$ [-]')
#ax.set_xlabel(r'Density ratio $\frac{\rho_{body}}{\rho_{fluid}} \left(= \frac{1}{r}\right)$')
#ax.set_ylabel(r'Acceleration relaxation, $\gamma$, ($a^{n+1} = \gamma a^{n+1} + (1-\gamma) a^n$)')
ax.legend(['Not diverged', 'Diverged', r'$\gamma = \frac{2}{1+\rho_f/\rho_b}$'],loc='lower right', fancybox=True, framealpha=0.8)
#ax.axis('equal')
#plt.gca().set_adjustable("box")
ax.grid(True)
plt.gcf().set_size_inches(4, 4)
plt.show()
#plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')
