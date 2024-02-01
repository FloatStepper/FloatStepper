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
figname="floatStepperVs6DoFRBM.pdf"

#Import FloatSteppper data
os.chdir('bodyDensity0.8FloatStepper')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t0 = np.loadtxt('trajectory/t', unpack=True)
ax0, ay0 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e0 = (ay0-ay0[0])/ay0[0]

#Import stable sixDoFRBM data
os.chdir('../bd0.8ar0.8_nOC1')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
ax1, ay1 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e1 = (ay1-ay0[0])/ay0[0]

#Import unstable sixDoFRBM data
os.chdir('../bd0.8ar0.8_nOC3')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t2 = np.loadtxt('trajectory/t', unpack=True)
ax2, ay2 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e2 = (ay2-ay0[0])/ay0[0]

#Import unstable sixDoFRBM data
os.chdir('../bd0.8ar0.8_nOC5')
#os.chdir('../bodyDensity0.8FloatStepper_largerDomain')
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t3 = np.loadtxt('trajectory/t', unpack=True)
ax3, ay3 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
e3 = (ay3-ay0[0])/ay0[0]

#Plotting acceleration
fig, ax = plt.subplots(1)
ax.plot(t0, e0,'o-', zorder=0, linewidth=1)
ax.plot(t1, e1,'x-', zorder=-20, linewidth=1)
ax.plot(t2, e2,'+-', zorder=-15, linewidth=1)
ax.plot(t3, e3,'d-', zorder=-10, linewidth=1)
ax.plot(t0, 0*e0, 'k', zorder=0, linewidth=2)

ax.set_xlabel(r'$t$ [s]')
ax.set_ylabel(r'$(a(t)-a_{Exact})/a_{Exact}$ [-]')
ax.grid(True)
ax.set_xlim(0, 2.5)
ax.set_ylim(-0.003,0.003)
#ax.set_yticks(np.arange(-.01, .01001, step=0.005)) 
ax.legend(['FloatStepper','sixDoFRigidBodyMotion - 1 OC','sixDoFRigidBodyMotion - 3 OC','sixDoFRigidBodyMotion - 5 OC'])
plt.gcf().set_size_inches(8, 3.5)

#plt.show()
os.chdir('../')
plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')
