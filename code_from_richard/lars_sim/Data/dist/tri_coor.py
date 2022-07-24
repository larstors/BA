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


c = np.genfromtxt("tri_coor.txt", delimiter=" ")
f = np.genfromtxt("fourier_tr_3.txt", delimiter=" ")

"""
plt.scatter(x=c[1,:], y =c[0,:])
plt.show()
"""
f = f/4700*1/101
sc = plt.scatter(x=c[1,:], y =c[0,:], s = 4, vmin = -10, marker="s", c=f, cmap='viridis')
plt.colorbar(sc)
plt.show()

sc = plt.scatter(x=c[1,:], y =c[0,:], s = 4, vmin = -6, marker="s", c=np.log(f), cmap='viridis')
plt.colorbar(sc)
plt.show()