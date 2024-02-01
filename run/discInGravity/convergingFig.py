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
figname="looseCouplingConv.pdf"

print(figdir)

os.chdir('convergingDivergingCases/rhob1.1aRelax1nOC1')

#Import CFD data
os.system("rm -rf trajectory")

#Extracting body trajectory data
rc = subprocess.call("./extractTrajectory")
t0 = np.loadtxt('trajectory/t', unpack=True)
ax0, ay0 = np.loadtxt('trajectory/a', unpack=True, usecols = (0, 1))

#Import theoretical trajectory data. Cols: Time x y th vx vy omega
#Running KK solver and loading theoretical trajectory
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th2, vx2, vy2, omega2 = np.loadtxt('KKtrajectory', unpack=True)
ax2 = (vx2[1:] - vx2[:-1])/(t2[1:] - t2[:-1])
ay2 = (vy2[1:] - vy2[:-1])/(t2[1:] - t2[:-1])
#alpha2 = (omega2[1:] - omega2[:-1])/(t2[1:] - t2[:-1])

#Plotting acceleration
fig, ax = plt.subplots(1)
col = 'r'
ax.plot(t0, ay0,'.-'+col, zorder=0, linewidth=1, label=r'$\rho_b/\rho_f = 1.1$')
ax.plot(t0, ay0,'.'+col, zorder=2)
ax.set_ylabel(r'$\ddot{y}$ [m/s$^2$]')
#ax.tick_params(axis ='y', labelcolor = col)

ax.plot(t2[:-1], ay2, 'k', zorder=-30, linewidth=4, label=r"$g_y(\rho_f-\rho_b)/(\rho_f+\rho_b)$")
#ax2.plot(t2[:-1], ay2, zorder=0, linewidth=3)
ax.set_xlabel(r'$t$ [s]')
#ax.set_ylabel('acceleration, a')
ax.legend()
ax.grid(True)
ax.set_xlim(0, 2.5)
#ax2.set_ylim(ay2[0]-0.3,ay2[0]+0.3)
#ax.set_ylim(-1e4,1e4)
plt.gcf().set_size_inches(5, 2.5)

plt.show()
os.chdir('../../')
#plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')