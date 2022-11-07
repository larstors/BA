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

L = 100
c = np.genfromtxt("tri_coor.txt", delimiter=" ")
f = np.genfromtxt("fourier_sq_3_special.txt", delimiter=" ")
"""
plt.scatter(x=c[1,:], y =c[0,:])
plt.show()
"""

x = np.array([np.ones(100) for i in range(100)])
#y = np.array([np.ones(100)*(i-50) for i in range(100)])

for i in range(100):
    x[:, i] = - np.pi + i * np.pi/50

y = x.T
x = x.flatten()
y = y.flatten()

fig = plt.figure()
f = f/(5500)*1/101
sc = plt.scatter(x=x, y =y, s = 4, vmin = 0, marker="s", c=f, cmap='viridis')
plt.axis([-np.pi, np.pi, -np.pi, np.pi])
plt.colorbar(sc)
plt.xlabel(r"$k_x$")
plt.ylabel(r"$k_y$")
#plt.show()
plt.savefig("squ_3_special.pdf")

fig2, ax = plt.subplots()
sc = ax.scatter(x=x, y =y, s = 4, vmin = -6, marker="s", c=np.log(f), cmap='viridis')
ax.axis([-np.pi, np.pi, -np.pi, np.pi])
cbar = plt.colorbar(sc)
cbar.ax.set_ylabel(r"$\logS(k)$")
ax.set_xlabel(r"$k_x$")
ax.set_ylabel(r"$k_y$")
#plt.show()
plt.savefig("squlog_3_special.pdf")

"""
f1 = np.genfromtxt("fourier_sq_1.txt", delimiter=" ") * 1/(4700*101)
f2 = np.genfromtxt("fourier_sq_2.txt", delimiter=" ") * 1/(5500*101)
f3 = np.genfromtxt("fourier_sq_3.txt", delimiter=" ") * 1/(5800*101)

fig3 = plt.figure()
plt.plot(np.linspace(0, np.pi, len(f1[50*100+50: 50*100+100])), f1[50*100+50: 50*100+100], "-o", label=r"$n_\mathrm{max}=1$")
plt.plot(np.linspace(0, np.pi, len(f2[50*100+50: 50*100+100])), f2[50*100+50: 50*100+100], "-o", label=r"$n_\mathrm{max}=2$")
plt.plot(np.linspace(0, np.pi, len(f3[50*100+50: 50*100+100])), f3[50*100+50: 50*100+100], "-o", label=r"$n_\mathrm{max}=3$")
plt.legend()
plt.yscale("log")
plt.axis([0, np.pi, 0.001, 10**3])
plt.xlabel(r"k")
plt.ylabel(r"$S($k$)$")
plt.savefig("squ_across.pdf")
#plt.show()


Sq = np.asarray(f1[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)

Sq = np.asarray(f2[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)

Sq = np.asarray(f3[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)
"""

fig3 = plt.figure()
plt.plot(np.linspace(0, np.pi, len(f[50*100+50: 50*100+100])), f[50*100+50: 50*100+100], "-o", label=r"$Special case$")
plt.legend()
plt.yscale("log")
plt.axis([0, np.pi, 0.001, 10**3])
plt.xlabel(r"k")
plt.ylabel(r"$S($k$)$")
#plt.savefig("squ_across_special.pdf")
plt.show()


Sq = np.asarray(f[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)
