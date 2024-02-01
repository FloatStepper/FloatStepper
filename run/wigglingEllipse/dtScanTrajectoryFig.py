#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import math
import subprocess
import os

cwd = os.getcwd()
figdir = os.path.join(cwd, "../figures")
figname="wigglingEllipseDtScanTrajectories.pdf"

scanName='maxCoScan'
maxColist=["0.05", "0.1", "0.2", ".5"]
nx='400'
tmax = 6

os.chdir(scanName)

fig, ax = plt.subplots()

# Major ticks every 20, minor ticks every 5
xmajor_ticks = np.arange(0, 12, .5)
xminor_ticks = np.arange(0, 12, .1)
ymajor_ticks = np.arange(-1, 1, .5)
yminor_ticks = np.arange(-1, 1, .1)

ax.set_xticks(xmajor_ticks)
ax.set_xticks(xminor_ticks, minor=True)
ax.set_yticks(ymajor_ticks)
ax.set_yticks(yminor_ticks, minor=True)

# And a corresponding grid
ax.grid(which='major')
ax.grid(which='minor')

plt.grid(b=True, which='major', color='0.0', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='.5', linestyle=':')

plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
#plt.title(r'Zero mass ellipse with $v_x(0) = 1, v_y(0) = 0, \omega(0) = 1$')

ax.axis('equal')
plt.gca().set_adjustable("box")
#ax.set_xlim(-0.6, 6)
ax.set_ylim(-.4,.6)

dtquiver = 0.4

for i in range(len(maxColist)):
    caseDir='maxCo'+maxColist[i]+'nx'+nx
    os.chdir(caseDir)

    #Extracting and loading trajectory data from floaterFoam log file
    os.system("rm -rf trajectory")
    rc = subprocess.call("./extractTrajectory")
    t1 = np.loadtxt('trajectory/t', unpack=True)
    x1, y1 = np.loadtxt('trajectory/x', unpack=True, usecols = (0, 1))
    costh, msinth = np.loadtxt('trajectory/Q', unpack=True, usecols = (0, 1))
    x1 = x1[t1 <= tmax]
    y1 = y1[t1 <= tmax]
    costh = costh[t1 <= tmax]
    msinth = msinth[t1 <= tmax]
    t1 = t1[t1 <= tmax]
    th1 = np.arctan2(-msinth, costh)
    line, = plt.plot(x1, y1, linewidth=2, zorder=-i,label=str(maxColist[i]))
    col = line.get_color() 

    #Distance between orinetation arrows along trajectory
    tquiv = np.arange(start=t1[0], stop=t1[-1], step=dtquiver)
    xquiv = np.interp(tquiv,t1,x1)
    yquiv = np.interp(tquiv,t1,y1)
    costhquiv = np.interp(tquiv,t1,costh)
    msinthquiv = np.interp(tquiv,t1,msinth)
    plt.quiver(xquiv, yquiv, -msinthquiv, -costhquiv, scale=15, pivot='mid', headlength=0, headwidth = 1, width=.002, zorder=-i, color=col)
    os.chdir('../')

#Running KK solver and loading theoretical trajectory
os.chdir('baseCase')
os.system("KirchhoffKelvinIntegrator2D RKDP45 > log.KirchhoffKelvinIntegrator2D 2>&1")
t2, x2, y2, th2, = np.loadtxt('KKtrajectory', unpack=True, usecols = (0,1,2,3))
x2 = x2[t2 <= tmax]
y2 = y2[t2 <= tmax]
th2 = th2[t2 <= tmax]
t2 = t2[t2 <= tmax]

#Add zeros in front
t2 = np.append([0.0],t2)
x2 = np.append([0],x2)
y2 = np.append([0],y2)
th2 = np.append([math.pi/2],th2)
nx2 = np.cos(th2)
ny2 = np.sin(th2)

x2quiv = np.interp(tquiv,t2,x2)
y2quiv = np.interp(tquiv,t2,y2)
nx2quiv = np.interp(tquiv,t2,nx2)
ny2quiv = np.interp(tquiv,t2,ny2)

plt.plot(x2, y2,'-k', linewidth=4,zorder=-10,label='Exact')
plt.quiver(x2quiv, y2quiv, nx2quiv, ny2quiv, scale=15, pivot='mid', headlength=0, headwidth = 1, color='black', width=.004,zorder=-10)
ax.legend(ncol=5,loc='upper center')

plt.gcf().set_size_inches(10, 3)
#plt.savefig( os.path.join(figdir, figname), bbox_inches='tight')
plt.show()

