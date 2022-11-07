import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import time
import random as rn
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import matplotlib.patches as mpl_patches
from matplotlib import rcParams
from io import StringIO
import matplotlib.tri as mtri
import matplotlib
matplotlib.use('Agg')

plt.rcParams.update({'font.size': 12})

sq2f = np.loadtxt("square_perc_fhyst_2.txt", delimiter=" ")
sq2b = np.loadtxt("square_perc_bhyst_2.txt", delimiter=" ")

sq3f = np.loadtxt("square_perc_fhyst_3.txt", delimiter=" ")
sq3b = np.loadtxt("square_perc_bhyst_3.txt", delimiter=" ")


fig = plt.figure()
#plt.plot(sq2f[:, 0], sq2f[:, 1], label=r"$J$, forward")
#plt.plot(sq2f[:, 0], sq2f[:, 3], label=r"$M$, forward")
#plt.plot(sq2f[:, 0], sq2f[:, 5], label=r"$w_N$, forward")
plt.plot(sq2f[:, 0], sq2f[:, 7], "m-x", label=r"$c_N$, forward")

#plt.plot(sq2b[:, 0], sq2b[:, 1], label=r"$J$, backward")
#plt.plot(sq2b[:, 0], sq2b[:, 3], label=r"$M$, backward")
#plt.plot(sq2b[:, 0], sq2b[:, 5], label=r"$w_N$, backward")
plt.plot(sq2b[:, 0], sq2b[:, 7], "g-o", label=r"$c_N$, backward")
plt.ylabel(r"$c_N$")
plt.xscale("log")
plt.xlabel(r"$\alpha$")
plt.grid()
plt.legend()
plt.savefig("phase_sq_2_hyst.pdf", dpi=150)


fig2 = plt.figure()
#plt.plot(sq3f[:, 0], sq3f[:, 1], label=r"$J$, forward")
#plt.plot(sq3f[:, 0], sq3f[:, 3], label=r"$M$, forward")
#plt.plot(sq3f[:, 0], sq3f[:, 5], label=r"$w_N$, forward")
plt.plot(sq3f[:, 0], sq3f[:, 7], "m-x", label=r"$c_N$, forward")

#plt.plot(sq3b[:, 0], sq3b[:, 1], label=r"$J$, backward")
#plt.plot(sq3b[:, 0], sq3b[:, 3], label=r"$M$, backward")
#plt.plot(sq3b[:, 0], sq3b[:, 5], label=r"$w_N$, backward")
plt.plot(sq3b[:, 0], sq3b[:, 7], "g-o", label=r"$c_N$, backward")
plt.ylabel(r"$c_N$")
plt.xscale("log")
plt.xlabel(r"$\alpha$")
plt.grid()
plt.legend()
plt.savefig("phase_sq_3_hyst.pdf", dpi=150)
