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
f = np.genfromtxt("fourier_tr_1.txt", delimiter=" ")
"""
plt.scatter(x=c[1,:], y =c[0,:])
plt.show()
"""

"""

fig = plt.figure()
f = f/(4000)*1/101
sc = plt.scatter(x=c[3,:], y =c[2,:], s = 4, vmin = 0, marker="s", c=f, cmap='viridis')
plt.colorbar(sc)
#plt.show()
plt.savefig("tri_1.pdf")

fig2 = plt.figure()
sc = plt.scatter(x=c[3,:], y =c[2,:], s = 4, vmin = -6, marker="s", c=np.log(f), cmap='viridis')
plt.colorbar(sc)
#plt.show()
plt.savefig("trilog_1.pdf")


fig3 = plt.figure()
plt.plot(np.linspace(min(c[3,:]), max(c[3,:]), len(f[50*100: 50*100+100])), f[50*100: 50*100+100])
plt.savefig("tri_1_across.pdf")

print(f[50*100+99])


wtf = np.genfromtxt("tri_wtf.txt", delimiter=" ")

w = wtf[0, :]
fig3 = plt.figure()
plt.hist(w, bins=20)
plt.savefig("xwtf.pdf")

for i in range(len(wtf)):
    print(np.sum(wtf[i, :]))
"""

f1 = np.genfromtxt("fourier_tr_1.txt", delimiter=" ")
f2 = np.genfromtxt("fourier_tr_2.txt", delimiter=" ")
f3 = np.genfromtxt("fourier_tr_3.txt", delimiter=" ")

fig3 = plt.figure()
plt.plot(np.linspace(0, np.pi, len(f1[50*100+50: 50*100+100])), f1[50*100+50: 50*100+100], label="1")
plt.plot(np.linspace(0, np.pi, len(f2[50*100+50: 50*100+100])), f2[50*100+50: 50*100+100], label="2")
plt.plot(np.linspace(0, np.pi, len(f3[50*100+50: 50*100+100])), f3[50*100+50: 50*100+100], label="3")
plt.legend()
#plt.savefig("tri_1_across.pdf")
plt.show()


Sq = np.asarray(f1[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)

Sq = np.asarray(f2[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)

Sq = np.asarray(f3[50*100 + 50: 50*100+100])
i = np.argmax(Sq)
print(i)