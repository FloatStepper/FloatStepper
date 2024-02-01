#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys

inputfilepath=sys.argv[1]

t, x, y, th, vx, vy, omega = np.loadtxt(inputfilepath, unpack=True)
#Fz = np.loadtxt(inputfilepath)

nx = np.cos(th)
ny = np.sin(th)
nq = 5
#fig, ax = plt.subplots()
plt.plot(x,y,'.r')
plt.quiver(x[0::nq], y[0::nq], nx[0::nq], ny[0::nq])
plt.xlabel('x')
plt.ylabel('y')
plt.axis('equal')
plt.title('Ellipse trajectory')
#plt.legend()
plt.show()

