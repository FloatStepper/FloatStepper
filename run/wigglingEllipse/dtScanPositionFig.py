#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

cwd = os.getcwd()
figdir = os.path.join(cwd, "../figures")
figname="wigglingEllipseDtScanPositions.pdf"

maxColist=["0.05", "0.1", "0.2", ".5"]
scanName='maxCoScan'
legends = ['0.05','0.1','0.2', '0.5','Exact']
nx='200'
tmax = 3

os.chdir(scanName)

fig, ax = plt.subplots(3)

#Initial x velocity
vx0 = 1

for i in range(len(maxColist)):
    caseDir='maxCo'+maxColist[i]+'nx'+nx
    os.chdir(caseDir)

    #Extracting and loading trajectory data from floaterFoam log file
    os.system("rm -rf trajectory")
    rc = subprocess.call("./extractTrajectory")
    t1 = np.loadtxt('trajectory/t', unpack=True)
    x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
    costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))
    x1 = x1[t1 <= tmax]
    y1 = y1[t1 <= tmax]
    costh = costh[t1 <= tmax]
    msinth = msinth[t1 <= tmax]
    t1 = t1[t1 <= tmax]
    th1 = np.arctan2(-msinth, costh)

    ax[0].plot(t1, x1-vx0*t1, zorder=-i)
    ax[0].set_ylabel(r'$x-v_x(0)t$')
    ax[0].legend(['floaterFoam','KirchhoffKelvin'],loc='upper right')
    ax[0].grid(True)
    ax[1].plot(t1, y1, zorder=-i)
    ax[1].set_ylabel(r'$y$')
    ax[1].grid(True)
    ax[2].plot(t1, th1, zorder=-i)
    ax[2].set_ylabel(r'$\theta$')
    ax[2].set_xlabel(r'$t$')
    ax[2].grid(True)

    os.chdir('../')

#Running KK solver and loading theoretical trajectory
os.chdir('baseCase')
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th2, = np.loadtxt('KKtrajectory', unpack=True, usecols = (0,1,2,3))
x2 = x2[t2 <= tmax]
y2 = y2[t2 <= tmax]
th2 = th2[t2 <= tmax]
t2 = t2[t2 <= tmax]
th2 = th2 + math.pi/2

#Add zero
t2 = np.append([0.0],t2)
x2 = np.append([0],x2)
y2 = np.append([0],y2)
th2 = np.append([0],th2)

ax[0].plot(t2, x2-vx0*t2, '-k',linewidth=3,zorder=-10)
ax[1].plot(t2, y2, '-k',linewidth=3,zorder=-10)
ax[2].plot(t2, th2, '-k',linewidth=3,zorder=-10)

ax[0].legend(legends,loc='upper center',ncol=5)
ax[0].set_ylim(top=0.4)

plt.gcf().set_size_inches(8, 6)
plt.show()
#plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')

