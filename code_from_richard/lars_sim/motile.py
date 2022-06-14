
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import scipy
import matplotlib.animation as animation
#from celluloid import Camera
import matplotlib.tri as mtri
from io import StringIO
import matplotlib.cm as cm

data = np.genfromtxt("./Data/motility/hexagonal.txt")

fig, ax1 = plt.subplots()

left, bottom, width, height = [0.6, 0.3, 0.2, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])


ax1.plot(data[:, 0], data[:, 1], "rs", label="J")
ax1.errorbar(data[::4, 0], data[::4, 1], yerr=np.sqrt(data[::4, 2]), fmt="s", color="red")
ax1.plot(data[:, 0], data[:, 3], "bo", label="M")
ax1.legend()
ax1.set_xlabel(r"$\alpha$")
ax1.set_ylabel(r"$J,M$")
ax1.grid()

#ax2.plot(data[:30, 0], data[:30, 3], "bo")

plt.show()

data = np.genfromtxt("./Data/motility/hexagonal_2.txt")

fig2, ax1 = plt.subplots()

left, bottom, width, height = [0.6, 0.3, 0.2, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])


ax1.plot(data[:, 0], data[:, 1], "rs", label="J")
ax1.errorbar(data[::4, 0], data[::4, 1], yerr=np.sqrt(data[::4, 2]), fmt="s", color="red")
ax1.plot(data[:, 0], data[:, 3], "bo", label="M")
ax1.legend()
ax1.set_xlabel(r"$\alpha$")
ax1.set_ylabel(r"$J,M$")
ax1.grid()

#ax2.plot(data[:30, 0], data[:30, 3], "bo")

plt.show()
