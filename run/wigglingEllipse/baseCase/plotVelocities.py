#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

#Extracting body trajectory data from log
os.system("rm -rf trajectory")
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
vx1, vy1 = np.loadtxt('trajectory/v', unpack=True, usecols = (0, 1))
#ax, ay = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
#costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))
omega1 = np.loadtxt('trajectory/omega', unpack=True, usecols = (2))
#alpha = np.loadtxt('trajectory/alpha', unpack=True, usecols = (2))

#Import theoretical trajectory data. Cols: Time x y th vx vy omega
#Running KK solver and loading theoretical trajectory
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th2, vx2, vy2, omega2 = np.loadtxt('KKtrajectory', unpack=True)

#Plotting velocities
fig, ax = plt.subplots(3)
fig.suptitle('CFD vs theoretical velocities of freely moving ellipse')
ax[0].plot(t1, vx1,'.r')
ax[0].plot(t2, vx2,'.b')
#ax[0].plot(t1, vx**2 + vy**2,'.k')
#ax[0].plot(t2, vx2**2 + vy2**2,'.k')
ax[0].set_ylabel('v_x')
ax[0].legend(['FloatStepper','KirchhoffKelvin'],loc='upper right')
ax[0].grid(True)
ax[1].plot(t1, vy1,'.r')
ax[1].plot(t2, vy2,'.b')
ax[1].set_ylabel('v_y')
ax[1].grid(True)
ax[2].plot(t1, omega1,'.r')
ax[2].plot(t2, omega2,'.b')
ax[2].set_ylabel('omega')
ax[2].set_xlabel('time')
ax[2].grid(True)
plt.show()

