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
'font.size': 16
})

cwd = os.getcwd()
figdir = os.path.join(cwd, "../figures")
#os.mkdir(figdir)
figname="aRelaxConvAndDiv.pdf"

print(figdir)

os.chdir('interIsoFoam_reconPhi/bodyDensity0.8accelerationRelaxation0.8')

#Import CFD data
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t0 = np.loadtxt('trajectory/t', unpack=True)
ax0, ay0 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))

os.chdir('../bodyDensity0.8accelerationRelaxation0.9')

#Import CFD data
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t1 = np.loadtxt('trajectory/t', unpack=True)
ax1, ay1 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))

#Import theoretical trajectory data. Cols: Time x y th vx vy omega
#Running KK solver and loading theoretical trajectory
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th2, vx2, vy2, omega2 = np.loadtxt('KKtrajectory', unpack=True)
ax2 = (vx2[1:] - vx2[:-1])/(t2[1:] - t2[:-1])
ay2 = (vy2[1:] - vy2[:-1])/(t2[1:] - t2[:-1])

#Plotting acceleration
fig, ax = plt.subplots(1)
col = 'r'
ax.plot(t0, ay0,'.-'+col, zorder=0, linewidth=0.5, label=r"$\rho_b/\rho_f = 0.8, \gamma = 0.8$")
ax.plot(t0, ay0,'.'+col, zorder=2)
ax.set_ylabel(r'$\ddot{y}$ [m/s$^2$]')

col = 'b'
ax.plot(t1, ay1,'.-'+col, zorder=-20, linewidth=0.5, label=r"$\rho_b/\rho_f = 0.8, \gamma = 0.9$")
ax.plot(t1, ay1,'.'+col, zorder=-20)

ax.plot(t2[:-1], ay2, 'k', zorder=-10, linewidth=5, label=r"$g_y(\rho_f-\rho_b)/(\rho_f+\rho_b)$")
ax.set_xlabel(r'$t$ [s]')
ax.legend(loc='upper right')
ax.grid(True)
ax.set_xlim(0, 2)
ax.set_ylim(ay2[0]-0.3,ay2[0]+0.3)
plt.gcf().set_size_inches(7, 3.5)

#plt.show()
os.chdir('../../')
plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')
