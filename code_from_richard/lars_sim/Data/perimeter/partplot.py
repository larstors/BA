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
plt.rcParams.update({'font.size': 22})

L = [100, 100]
 
f = open("./square_particles.txt")
data = []
for line in f:
    data.append(np.loadtxt(StringIO(line), dtype=int))
 
def conv(n):
    x = []
    y = []
    for i in n:
        x.append(i % L[0])
        y.append(int(i/L[0]))
    return x, y


plt.scatter(x=conv(data[-1][0::2])[0], y=conv(data[-1][0::2])[1], c=data[-1][1::2], cmap="Greens", s=4, vmin=0, vmax=3, marker="s")
plt.savefig("wooba.pdf")