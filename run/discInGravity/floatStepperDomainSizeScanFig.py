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
'font.size': 12
})

cwd = os.getcwd()
figdir = os.path.join(cwd, "../figures")
#os.mkdir(figdir)
figname="floatStepperDomainSizeScan.pdf"

#Import FloatSteppper data
os.chdir('floatStepperDomainSizeScan/R40nx200')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t0 = np.loadtxt('trajectory/t', unpack=True)
ax0, ay0 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e0 = (ay0-ay0[0])/ay0[0]

#Import medium data
os.chdir('../R60nx300')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
ax1, ay1 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e1 = (ay1-ay0[0])/ay0[0]

#Import largest data
os.chdir('../R80nx400')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t2 = np.loadtxt('trajectory/t', unpack=True)
ax2, ay2 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e2 = (ay2-ay0[0])/ay0[0]

#Plotting acceleration
fig, ax = plt.subplots(1)
ax.plot(t0, e0,'o-', zorder=0, linewidth=1)
ax.plot(t1, e1,'x-', zorder=-10, linewidth=1)
ax.plot(t2, e2,'+-', zorder=-15, linewidth=1)
ax.plot(t0, 0*e0, 'k', zorder=0, linewidth=2)

ax.set_xlabel(r'$t$ [s]')
ax.set_ylabel(r'$(a(t) - a_\textrm{Exact})/a_\textrm{Exact}$ [-]')
ax.grid(True)
ax.set_xlim(0, 3)
ax.legend(['Outer rim: 40 m','Outer rim: 60 m','Outer rim: 80 m'], bbox_to_anchor=(1, .8), bbox_transform=ax.transAxes)
plt.gcf().set_size_inches(2.5, 3.5)

plt.show()
os.chdir('../')
#plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')
