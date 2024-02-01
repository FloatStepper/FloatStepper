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
    "font.size":12
})


#Wave height: '04' or '10' cm
H = '10'
T = 1.2
lmbda = 1.93625
k = 2.0*math.pi/lmbda
d = 0.4

cwd = os.getcwd()
figdir = '.' #os.path.join(cwd, "/")
figname = "RenEtAl2015H0"+H+".pdf"


# Load expeimental data
tEta, eta2 = np.loadtxt('Ren2015data/H0'+H+'WaveGauge', unpack=True)
eta2 = d*eta2
#Note: We do not know where Ren2015 measured surface elevation, but assuming it
#is close to the wave make the only difference should be a phase shift.
#This justifies using another time shifting for wave elevation than for body
#coordinates
if H == '04':
    tEtaShift = 4.4
elif H == '10':
    tEtaShift = 3.4
tEta = tEta + tEtaShift
tEta = tEta*T

tSurge, x2 = np.loadtxt('Ren2015data/H0'+H+'Surge', unpack=True)
x2 = d*x2
if H == '04':
    tSurgeShift = 4.37
elif H == '10':
    tSurgeShift = 3.3
tSurge = tSurge + tSurgeShift
tSurge = tSurge*T

tHeave, y2 = np.loadtxt('Ren2015data/H0'+H+'Heave', unpack=True)
y2 = d*y2
tHeaveShift = tSurgeShift
tHeave = tHeave + tHeaveShift
tHeave = tHeave*T

tPitch, th2 = np.loadtxt('Ren2015data/H0'+H+'Pitch', unpack=True)
th2 = th2*k*d*360/(2.0*np.pi)
tPitchShift = tSurgeShift
tPitch = tPitch + tPitchShift
tPitch = tPitch*T

if H == '04':
    tmin = 2*T
    tmax = 9.5*T
elif H == '10':
    tmin = 2*T
    tmax = 8.5*T

expLineColor = '.k'
expLineWidth = 2
# Plot experimental data
fig, ax = plt.subplots(4)
fig.set_size_inches(6.5, 10)
ax[0].set_title(r'$H = 0.'+H+'$ m, $T = 1.2$ s')
ax[0].plot(tEta, eta2, expLineColor, linewidth=expLineWidth, label=r"$x = ?$ m", zorder=-10)
ax[0].set_ylabel(r'$\eta$ [m]')

ax[1].plot(tSurge, x2, expLineColor, linewidth=expLineWidth, label=r"Exp. (Ren et al. 2015)", zorder=-10)
ax[1].set_ylabel(r'$x$ [m]')

ax[2].plot(tHeave, y2, expLineColor, linewidth=expLineWidth, zorder=-10)
ax[2].set_ylabel(r'$y$ [m]')

ax[3].plot(tPitch, th2, expLineColor, linewidth=expLineWidth, zorder=-10)
ax[3].set_ylabel(r'$\theta$ [deg]')
ax[3].set_xlabel(r'$t$ [s]')

caseDirs=['H0'+H+'_x0Eq4m', 'H0'+H+'_x0Eq4m_CFL01', 'H0'+H+'_x0Eq4m_fine']
legends = [r'coarse, CFL $\leq$ 0.5', r'coarse, CFL $\leq$ 0.1', r'fine, CFL $\leq$ 0.5']
lineWidths = ['3', '2', '1']
lineColors = ['r', 'm', 'c']

for i in range(len(caseDirs)):
    caseDir=str(caseDirs[i])
    os.chdir(caseDir)

    # Extracting body trajectory
    os.system("rm -rf trajectory")
    rc = subprocess.call("./extractTrajectory")
    t1 = np.loadtxt('trajectory/t', unpack=True)
    x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 2))
    costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 2))
    th1 = np.arctan2(-msinth, costh)

    tWGs, WG1, WG2, WG3, WG4, WG5 = np.loadtxt('postProcessing/interfaceHeight1/0/height.dat', unpack=True, usecols = (0, 2, 4, 6, 8, 10))


    # Nondimensionalise numerical data
    x1 = (x1-x1[0])
    y1 = (y1-d)
    th1 = th1*360/2.0/math.pi
    #th1 = th1/(k*d)
    th1 = - th1


    t1shift = 0 #1.95-2.44
    #t1 = (t1-t1[0])/T + t1shift
#    t1 = t1/T + t1shift
#    tWGs = tWGs/T
    tWGs = tWGs
    #WG1 = WG1/d
    WG2 = WG2
    #WG3 = WG3/d
    #WG4 = WG4/d
    #WG5 = WG5/d

    # Plotting surface elevation
    #ax[0].plot(tWGs, WG1,'b', label=r"$x = 0$ m")
    ax[0].plot(tWGs, WG2, linewidth=lineWidths[i], label=legends[i], zorder=i)
    #ax[0].plot(tWGs, WG3,'g', label="x = 3")
    #ax[0].plot(tWGs, WG4,'r', label=r"$x = 4$ m")
    #ax[0].plot(tWGs, WG5,'m', label="x = 5")

    # Plotting position
#    ax[1].plot(t1, x1, lineColors[i], linewidth=lineWidths[i], label=legends[i], zorder=-i)
#    ax[2].plot(t1, y1, lineColors[i], linewidth=lineWidths[i], label=legends[i], zorder=-i)
#    ax[3].plot(t1, th1, lineColors[i], linewidth=lineWidths[i], label=legends[i], zorder=-i)
    ax[1].plot(t1, x1, linewidth=lineWidths[i], label=legends[i], zorder=i)
    ax[2].plot(t1, y1, linewidth=lineWidths[i], label=legends[i], zorder=i)
    ax[3].plot(t1, th1, linewidth=lineWidths[i], label=legends[i], zorder=i)
    os.chdir('../')

panels = range(4)
for i in panels:
    ax[i].set_xlim([tmin, tmax])
    ax[i].grid(True)

ax[1].legend()
#ax[1].legend(loc='upper left')
ax[0].set_xticklabels([])
if H == '04':
    ax[1].set_ylim([-.01, .25])
elif H == '10':
    ax[0].set_ylim([-.06, .06])
    ax[1].set_ylim([-.05, 1.2])
ax[1].set_xticklabels([])
ax[2].set_xticklabels([])

plt.show()
#plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')
