#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

floatStepperCase='baseCase'
analyticCase='analytic'

#Loading CFD results
os.chdir(floatStepperCase)
os.system("rm -rf trajectory")
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
costh1, msinth1 = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))

#Loading analytical results
os.chdir('../'+analyticCase)
os.system("rm -rf trajectory")
rc = subprocess.call("./extractTrajectory")
t2 = np.loadtxt('trajectory/t', unpack=True)
x2, y2 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
costh2, msinth2 = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))

#Interpolating CFD results to analytic time steps
x1tmp = x1
y1tmp = y1
costhtmp = costh1
msinthtmp = msinth1
x1 = np.interp(t2,t1,x1tmp)
y1 = np.interp(t2,t1,y1tmp)
costh1 = np.interp(t2,t1,costhtmp)
msinth1 = np.interp(t2,t1,msinthtmp)

#Plotting trajectories
#Distance between orinetation arrows along trajectory
nq = 20 #Approx number of quivers
dn1 = int(np.round(len(x1)/nq))
dn2 = int(np.round(len(x2)/nq))
plt.plot(x2, y2,'-k', linewidth=4)
plt.quiver(x2[0::dn2], y2[0::dn2], -msinth2[0::dn2], -costh2[0::dn2], scale=15, pivot='mid', headlength=0, headwidth = 1, color='black', width=.002)
plt.plot(x1, y1, '-r', linewidth=2)
plt.quiver(x1[0::dn1], y1[0::dn1], -msinth1[0::dn1], -costh1[0::dn1], scale=15, pivot='mid', headlength=0, headwidth = 1, color='red', width=.002)
plt.xlabel('x')
plt.ylabel('y')
#plt.set_aspect('equal', 'datalim')
plt.axis('equal')
#plt.set(xlim=(-.5, 10), ylim=(-2, 3))
plt.title(r'Zero mass ellipse with $v_x(0) = 1, v_y(0) = 0, \omega(0) = 1$')
plt.legend(['Theoretical (Kelvin-Kirchhoff)','FloatStepper'],loc='upper left')

##plt.ylim(-.3,.25)
##plt.text(0, -0.2, r'Start at $t = 0$', horizontalalignment='center',
##    verticalalignment='center')
##plt.text(3.1, -0.2, r'End at $t = 3$', horizontalalignment='center',
##    verticalalignment='center')
#plt.xlim(left = -10, right = 10)
#plt.ylim(ymax = 5, ymin = -5)
plt.show()
