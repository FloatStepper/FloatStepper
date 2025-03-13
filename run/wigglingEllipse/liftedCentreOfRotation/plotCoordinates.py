#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
plt.rcParams.update({'font.size': 12})

os.chdir(sys.argv[1])

#Extracting body trajectory data from log
os.system("rm -rf trajectory")
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))
th1 = np.arctan2(-msinth, costh)

#Import theoretical trajectory data. Cols: Time x y th vx vy omega
#Running KK solver and loading theoretical trajectory
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th2, vx2, vy2, omega2 = np.loadtxt('KKtrajectory', unpack=True)
nx2 = np.cos(th2)
ny2 = np.sin(th2)
th2 = th2 + math.pi/2

#Plotting position
fig, ax = plt.subplots(3)
fig.canvas.set_window_title(os.getcwd())
#fig.suptitle('CFD vs theoretical position of freely moving ellipse')
vx0 = 1
ax[0].plot(t1, x1-vx0*t1)
ax[0].plot(t2, x2-vx0*t2)
ax[0].set_ylabel(r'$x$')
ax[0].legend(['FloaterStepper','Kirchhoff-Kelvin'],loc='upper right')
ax[0].grid(True)
ax[1].plot(t1, y1)
ax[1].plot(t2, y2)
ax[1].set_ylabel(r'$y$')
ax[1].grid(True)
ax[2].plot(t1, th1)
ax[2].plot(t2, th2)
ax[2].set_ylabel(r'$\theta$')
ax[2].set_xlabel(r'$t$')
ax[2].grid(True)

plt.show()
