#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

#Extracting and loading trajectory data from floaterFoam log file
os.system("rm -rf trajectory")
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))

#Running KK solver and loading theoretical trajectory
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th, = np.loadtxt('KKtrajectory', unpack=True, usecols = (0,1,2,3))
nx2 = np.cos(th)
ny2 = np.sin(th)

x1tmp = x1
y1tmp = y1
costhtmp = costh
msinthtmp = msinth
x1 = np.interp(t2,t1,x1tmp)
y1 = np.interp(t2,t1,y1tmp)
costh = np.interp(t2,t1,costhtmp)
msinth = np.interp(t2,t1,msinthtmp)

#Distance between orinetation arrows along trajectory
nq = 20 #Approx number of quivers
dn1 = int(np.round(len(x1)/nq))
dn2 = int(np.round(len(x2)/nq))
#fig, ax = plt.subplots()
plt.plot(x2, y2,'-k', linewidth=4)
plt.quiver(x2[0::dn2], y2[0::dn2], nx2[0::dn2], ny2[0::dn2], scale=15, pivot='mid', headlength=0, headwidth = 1, color='black', width=.004)
plt.plot(x1, y1, '-r', linewidth=2)
plt.quiver(x1[0::dn1], y1[0::dn1], -msinth[0::dn1], -costh[0::dn1], scale=15, pivot='mid', headlength=0, headwidth = 1, color='red', width=.002)
plt.xlabel('x')
plt.ylabel('y')
#plt.set_aspect('equal', 'datalim')
plt.axis('equal')
#plt.set(xlim=(-.5, 1), ylim=(0, 3))
plt.title(r'Zero mass ellipse with $v_x(0) = 1, v_y(0) = 0, \omega(0) = 1$')
plt.legend(['Theoretical (Kirchhoff-Kelvin)','FloatStepper'],loc='upper left')

##plt.ylim(-.3,.25)
##plt.text(0, -0.2, r'Start at $t = 0$', horizontalalignment='center',
##    verticalalignment='center')
##plt.text(3.1, -0.2, r'End at $t = 3$', horizontalalignment='center',
##    verticalalignment='center')
#plt.xlim(left = -3, right = 3)
#plt.ylim(ymax = 3, ymin = 0)
plt.show()

