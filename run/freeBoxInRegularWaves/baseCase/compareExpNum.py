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

#Extracting body trajectory data from log
os.system("rm -rf trajectory")
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 2))
costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 2))
th1 = np.arctan2(-msinth, costh)

tWGs, WG1, WG2, WG3, WG4, WG5 = np.loadtxt('postProcessing/interfaceHeight1/0/height.dat', unpack=True, usecols = (0, 2, 4, 6, 8, 10))

T = 1.2
lmbda = 1.93625
k = 2.0*math.pi/lmbda
d = 0.4

tEta, eta2 = np.loadtxt('../Ren2015data/H010WaveGauge', unpack=True)
tEtaShift = 5/T - (4.175 - 4.583) -.1 #5 - (5.17-4.74)
tEta = tEta + tEtaShift

tSurge, x2 = np.loadtxt('../Ren2015data/H010Surge', unpack=True)
tSurgeShift = tEtaShift
tSurge = tSurge + tSurgeShift

tHeave, y2 = np.loadtxt('../Ren2015data/H010Heave', unpack=True)
tHeaveShift = tSurgeShift
tHeave = tHeave + tHeaveShift

tPitch, th2 = np.loadtxt('../Ren2015data/H010Pitch', unpack=True)
tPitchShift = tSurgeShift
tPitch = tPitch + tPitchShift

x1 = (x1-x1[0])/d
y1 = (y1-d)/d
#th1 = th1*360/2.0/math.pi
th1 = th1/(k*d)
th1 = - th1

t1shift = 0 #1.95-2.44
#t1 = (t1-t1[0])/T + t1shift
t1 = t1/T + t1shift
tWGs = tWGs/T
WG1 = WG1/d
WG2 = WG2/d
WG3 = WG3/d
WG4 = WG4/d
WG5 = WG5/d

#Plotting position
tmin = 6
tmax = 8.5
fig, ax = plt.subplots(4)
fig.set_size_inches(6.5, 10)
ax[0].set_title(r'$H = 0.10$ m, $T = 1.2$ s')
ax[0].plot(tEta, eta2,'^k', label=r"$x = ?$ m")
ax[0].plot(tWGs, WG1,'b', label=r"$x = 0$ m")
#ax[0].plot(tWGs, WG2,'c', label="x = 1")
#ax[0].plot(tWGs, WG3,'g', label="x = 3")
ax[0].plot(tWGs, WG4,'r', label=r"$x = 4$ m")
#ax[0].plot(tWGs, WG5,'m', label="x = 5")
ax[0].set_ylabel(r'$\eta/d$')
ax[0].legend(loc='upper right')
#ax[0].set_xlim([tmin, tmax])
#ax[0].legend(['floaterFoam','Exp'],loc='upper right')
ax[0].grid(True)
ax[0].set_xticklabels([])


ax[1].plot(tSurge, x2,'^k',label=r"Exp")
ax[1].plot(t1, x1,'r', label="CFD")
ax[1].set_ylabel(r'$x/d$')
#ax[1].set_xlim([tmin, tmax])
#ax[1].legend(['floaterFoam','Exp'],loc='upper right')
ax[1].legend()
ax[1].grid(True)
ax[1].set_xticklabels([])

ax[2].plot(tHeave, y2,'^k')
ax[2].plot(t1, y1,'r')
ax[2].set_ylabel(r'$y/d$')
ax[2].grid(True)
#ax[2].set_xlim([tmin, tmax])
ax[2].set_xticklabels([])

ax[3].plot(tPitch, th2,'^k')
ax[3].plot(t1, th1,'r')
ax[3].set_ylabel(r'$\theta/kd$')
ax[3].set_xlabel(r'$t/T$')
ax[3].grid(True)
#ax[3].set_xlim([tmin, tmax])
plt.show()
#plt.savefig('Ren2015H010vsFloatStepper.pdf', bbox_inches='tight')