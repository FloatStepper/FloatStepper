#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

os.system("rm -rf trajectory")
os.system("rm -rf forceAndTorque")

#Extracting body trajectory data from log
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
x, y = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
vx, vy = np.loadtxt('trajectory/v', unpack=True, usecols = (0, 1))
ax, ay = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))
costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))
omega = np.loadtxt('trajectory/omega', unpack=True, usecols = (2))
alpha = np.loadtxt('trajectory/alpha', unpack=True, usecols = (2))

# Swapping costh and sinth to account for -90 degree turn of simulation
tmp = costh
costh = -msinth
msinth = tmp


#Calculating added mass diagonal - note eccentricity hardcoded here - not read
pi = math.pi
R = 1
a = 0.5
rhof = 1
A11 = rhof*pi*R**2*(1-a**2)**2
A22 = rhof*pi*R**2*(1+a**2)**2
A33 = 2*rhof*pi*R**4*a**4
DA = A22-A11

sinth = -msinth
vhatx = costh*vx + sinth*vy
vhaty = -sinth*vx + costh*vy
ahatx = costh*ax + sinth*ay
ahaty = -sinth*ax + costh*ay
Fhatx = DA*vhaty*omega - 0*A11*ahatx
Fhaty = DA*vhatx*omega - 0*A22*ahaty
tau1 = -DA*vhatx*vhaty - 0*A33*alpha
Fx1 = costh*Fhatx - sinth*Fhaty
Fy1 = sinth*Fhatx + costh*Fhaty
#tau1 = tau1 - x*Fy1 + y*Fx1

#Extracting force and torque data from postProcessing directory
rc = subprocess.call("./extractForceAndTorque")
t2 = np.loadtxt('forceAndTorque/t', unpack=True)
Fx2 = np.loadtxt('forceAndTorque/F0x', unpack=True)
Fy2 = np.loadtxt('forceAndTorque/F0y', unpack=True)
tau2 = np.loadtxt('forceAndTorque/tau0', unpack=True)

#Import theoretical trajectory data. Cols: Time x y th vx vy omega
#Running KK solver and loading theoretical trajectory
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t3, x3, y3, th3, vx3, vy3, omega3 = np.loadtxt('KKtrajectory', unpack=True)
ax3 = (vx3[1:] - vx3[:-1])/(t3[1:] - t3[:-1])
ay3 = (vy3[1:] - vy3[:-1])/(t3[1:] - t3[:-1])
alpha3 = (omega3[1:] - omega3[:-1])/(t3[1:] - t3[:-1])
ax3 = np.concatenate(([0.0],ax3))
ay3 = np.concatenate(([0.0],ay3))
alpha3 = np.concatenate(([0.0],alpha3))
costh = np.cos(th3)
sinth = np.sin(th3)
vhatx = costh*vx3 + sinth*vy3
vhaty = -sinth*vx3 + costh*vy3
ahatx = costh*ax3 + sinth*ay3
ahaty = -sinth*ax3 + costh*ay3
Fhatx = DA*vhaty*omega3 - 0*A11*ahatx
Fhaty = DA*vhatx*omega3 - 0*A22*ahaty
tau3 = -DA*vhatx*vhaty - 0*A33*alpha3
Fx3 = costh*Fhatx - sinth*Fhaty
Fy3 = sinth*Fhatx + costh*Fhaty

#Plotting
fig, ax1 = plt.subplots(3)
fig.suptitle('Pressure force vs theoretical for given body velocity')
ax1[0].plot(t1, Fx1,'.r')
ax1[0].plot(t2, Fx2,'.b')
ax1[0].plot(t3, Fx3,'.g')
ax1[0].set_ylabel('Fx')
ax1[0].legend(['Theory for given motion','floaterFoam','Exact orbit'])
#ax1[0].legend(['Theory given motion','floaterFoam'])
ax1[0].grid(True)
ax1[1].plot(t1, Fy1,'.r')
ax1[1].plot(t2, Fy2,'.b')
ax1[1].plot(t3, Fy3,'.g')
ax1[1].set_ylabel('Fy')
ax1[1].grid(True)
ax1[2].plot(t1, tau1,'.r')
ax1[2].plot(t2, tau2,'.b')
ax1[2].plot(t3, tau3,'.g')
ax1[2].set_ylabel('tau')
ax1[2].set_xlabel('time')
ax1[2].grid(True)
plt.show()
