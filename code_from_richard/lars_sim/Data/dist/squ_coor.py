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
f = np.genfromtxt("fourier_sq_3.txt", delimiter=" ")
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
f = f/(5800)*1/101
sc = plt.scatter(x=x, y =y, s = 4, vmin = 0, marker="s", c=f, cmap='viridis')
plt.colorbar(sc)
#plt.show()
plt.savefig("squ_3.pdf")

fig2 = plt.figure()
sc = plt.scatter(x=x, y =y, s = 4, vmin = -3, marker="s", c=np.log(f), cmap='viridis')
plt.colorbar(sc)
#plt.show()
plt.savefig("squlog_3.pdf")


fig3 = plt.figure()
plt.plot(np.linspace(min(x), max(x), len(f[50*100: 50*100+100])), f[50*100: 50*100+100])
plt.savefig("squ_across_3.pdf")

print(1/f[50*100+99], f[50*100+99])

"""
wtf = np.genfromtxt("tri_wtf.txt", delimiter=" ")

w = wtf[0, :]
fig3 = plt.figure()
plt.hist(w, bins=20)
plt.savefig("xwtf.pdf")

for i in range(len(wtf)):
    print(np.sum(wtf[i, :]))
"""