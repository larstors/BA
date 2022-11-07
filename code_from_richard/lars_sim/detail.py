#
#   This is for making plots regarding special interest points in the time evolution of the system
#
#
#
#
#


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



f = open("./Data/tri_detail/run_2.txt")
tri_data_n_3 = []
for line in f:
    tri_data_n_3.append(np.loadtxt(StringIO(line), dtype=int))

y_1 = tri_data_n_3[0][::2]
y_2 = tri_data_n_3[1][::2]
y_3 = tri_data_n_3[2][::2]

print(np.sum(y_1))
print(np.sum(y_2))
print(np.sum(y_3))

x_1 = np.argwhere(y_1!=0)
x_2 = np.argwhere(y_2!=0)
x_3 = np.argwhere(y_3!=0)

y_1 = y_1[x_1]#/np.sum(y_1)
y_2 = y_2[x_2]#/np.sum(y_2)
y_3 = y_3[x_3]#/np.sum(y_3)


plt.plot(x_1+1, y_1, "r-o", label=r"$t=0$")
plt.plot(x_2+1, y_2, "g-x", label=r"$t=500$")
plt.plot(x_3+1, y_3, "b-s", label=r"$t=10^4$")
#plt.axis([1, 2500, 1e-3, 1])
plt.legend()
plt.grid()
plt.xlabel(r"$k$")
plt.ylabel(r"$N_k$")
plt.yscale("log")
plt.xscale("log")
#plt.savefig("./plots/tri_detail_n_4.pdf", dpi=150)
plt.show()