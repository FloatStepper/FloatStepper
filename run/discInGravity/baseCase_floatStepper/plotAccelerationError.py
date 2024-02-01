#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

#Import CFD data
os.system("rm -rf trajectory")
os.system("rm -rf forceAndTorque")

#Extracting body trajectory data from log
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
#x, y = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
#vx, vy = np.loadtxt('trajectory/v', unpack=True, usecols = (0, 1))
ax1, ay1 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
#costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))
#omega = np.loadtxt('trajectory/omega', unpack=True, usecols = (2))
alpha1 = np.loadtxt('trajectory/alpha', unpack=True, usecols = (2))

#Import theoretical trajectory data. Cols: Time x y th vx vy omega
#Running KK solver and loading theoretical trajectory
#os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
#t2, x2, y2, th2, vx2, vy2, omega2 = np.loadtxt('KKtrajectory', unpack=True)
#ax2 = (vx2[1:] - vx2[:-1])/(t2[1:] - t2[:-1])
#ay2 = (vy2[1:] - vy2[:-1])/(t2[1:] - t2[:-1])
#alpha2 = (omega2[1:] - omega2[:-1])/(t2[1:] - t2[:-1])

#Plotting acceleration
fig, ax = plt.subplots(1)
fig.suptitle('CFD vs theoretical accelerations of freely moving ellipse')
ax.plot(t1, (ay1-ay1[0])/ay1[0],'.-')
ax.plot(t1, ay1*0)
ax.set_ylabel('Relative acceleration error')
ax.legend(['FloatStepper','Exact'],loc='upper right')
ax.grid(True)
#ax[1].plot(t1, ay1,'.r')
#ax[1].plot(t2[:-1], ay2,'.b')
#ax[1].set_ylabel('a_y')
#ax[1].grid(True)
#ax[2].plot(t1, alpha1,'.r')
#ax[2].plot(t2[:-1], alpha2,'.b')
#ax[2].set_ylabel('domegadt')
#ax[2].set_xlabel('time')
#ax[2].grid(True)
plt.show()
