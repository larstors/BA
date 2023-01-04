from wsgiref.headers import tspecials
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
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 15})

L = [100, 100]
n = 3

f = open("./square_snapnumber_alpha0.001_phi1.20_L100_3.txt")
data = []
for line in f:
    data.append(np.loadtxt(StringIO(line), dtype=int))



m = np.genfromtxt("./square_number_alpha0.001_phi1.20_L100_3.txt", delimiter=" ")
print(np.shape(m))

def conv(n):
    """_summary_

    Args:
        n (_type_): _description_

    Returns:
        _type_: _description_  
    """

    x = []
    y = []
    for i in n:
        x.append(i % L[0])
        y.append(int(i/L[0]))
    
    return x, y


fig, ax = plt.subplots()
x1, y1 = conv(data[0][0::2])
ax.scatter(x=x1, y=y1, c=data[0][1::2], cmap="Greens", s=1, vmin=0, vmax=n, marker="s")
plt.savefig("snap1.pdf")

fig1, ax1 = plt.subplots()
x1, y1 = conv(data[1][0::2])
ax1.scatter(x=x1, y=y1, c=data[1][1::2], cmap="Greens", s=1, vmin=0, vmax=n, marker="s")
plt.savefig("snap2.pdf")

fig2, ax2 = plt.subplots()
x1, y1 = conv(data[2][0::2])
ax2.scatter(x=x1, y=y1, c=data[2][1::2], cmap="Greens", s=1, vmin=0, vmax=n, marker="s")
plt.savefig("snap3.pdf")