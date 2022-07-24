"""
small program to produce some lattice examples for latex
"""

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
colormap = cm.inferno

print("L?")
L = int(input())


f = open("triangle.txt")
tr = []
for line in f:
    tr.append(np.loadtxt(StringIO(line), dtype=int))

n_tr = tr[0][1::3]
p_tr = tr[0][2::3]

f = open("hexagonal.txt")
hx = []
for line in f:
    hx.append(np.loadtxt(StringIO(line), dtype=int))
n_hx = hx[0][1::4]
i_hx = hx[0][2::4]
p_hx = hx[0][3::4]

def tri_conv(n):
    
        x = []
        y = []
        for i in n:
            x.append(i % L - 0.5 * int(i / L) - L/4)
            y.append(int(i/L) * np.sqrt(3)/2)
        
        return x, y


# triangle
def trian():
        # indices for entire lattice
        ind = np.array([i for i in range(L*L)])
        # x and y coordinate for each lattice point
        x = tri_conv(ind)[0]
        y = tri_conv(ind)[1]
        
        #for i in range(len(x)):
        #    x[i] += 10

        triangles = []

        for j in range(L-1):
            for i in range(L-1):
                i = i + j * L
                triangles.append([i, i+1, i+L+1])
                triangles.append([i+L, i, i+1+L])
        triang = mtri.Triangulation(x, y, triangles)
        return triang
    
triang = trian()

fig1 = plt.figure()
plt.triplot(triang, 'r-', alpha=1, linewidth=0.5, zorder=1)
plt.scatter(tri_conv(n_tr)[0], tri_conv(n_tr)[1], c=p_tr, cmap="Greens", s=100, vmin=0, vmax=3, marker="s", zorder=2)
#plt.axis("off")
#plt.savefig("tri_ex.pdf", dpi=100)
plt.show()

def hex_cor(n, j):
    """function that converts the index n and site index j into coordinates x and y, with x being horisontal.

        Args:
            n (int): index of particle in lattice

        Returns:
            tuple: x and y coordinates
    """
        
    xf = []
    yf = []

    c30 = np.cos(np.pi/6)
    s30 = np.sin(np.pi/6)

    for k in range(len(n)):
        y = int(n[k]/L)
        x = 2 * (n[k]%L) + j[k] - y
        x += 1
        xf.append(x * c30)
        if (y%2 == 1):
            if (x%2 == 1):
                yf.append(y * (1 + s30))
            else:
                yf.append((y * 3 + 1)/2)
        else:
            if (x%2 == 1):
                yf.append(1.5 * y + 0.5)
            else:
                yf.append(1.5 * y)

    return xf, yf

def hex_cor_int( n, j):
    xf = 0
    yf = 0

    c30 = np.cos(np.pi/6)
    s30 = np.sin(np.pi/6)

        
    y = int(n/L)
    x = 2 * (n%L) + j - y
    x += 1
    xf = (x * c30)

    if (y%2 == 1):
        if (x%2 == 1):
            yf = (y * (1 + s30))
        else:
            yf = ((y * 3 + 1)/2)
    else:
        if (x%2 == 1):
            yf = (1.5 * y + 0.5)
        else:
            yf = (1.5 * y)

    return xf, yf
    
def fill(xl, yl):
    for i in range(L):
        for k in range(2):
            f = 1 - k
            for j in range(L):  
                n = i * L + j

                if (f%2==0):
                    xl.append(hex_cor_int(n, f)[0])
                    yl.append(hex_cor_int(n, f)[1])
                else:
                    xl.append(hex_cor_int(n, f)[0])
                    yl.append(hex_cor_int(n, f)[1])
    return xl, yl



def grid_plot():

    xl = []
    yl = []
    xl, yl = fill(xl, yl)
    for n in range(0, L-1, 1):
        for i in range(L-1):
            coord = [[xl[i + 2*n*L], yl[i+ 2*n*L]], [xl[i + L + 1 + 2*n*L], yl[i+ L + 1 + 2*n*L]], [xl[i+ 2*L + 1+ 2*n*L], yl[i+ 2*L + 1+ 2*n*L]], [xl[i+ 3*L + 1+ 2*n*L], yl[i+ 3*L +1+ 2*n*L]], [xl[i+ L*2+ 2*n*L], yl[i+ L*2+ 2*n*L]], [xl[i+ L+ 2*n*L], yl[i+ L+ 2*n*L]], [xl[i+ 2*n*L], yl[i+ 2*n*L]]]
            xs, ys = zip(*coord)
            plt.plot(xs, ys, 'r-', alpha=1, linewidth=0.5, zorder=1)

    plt.plot([xl[L-1], xl[2*L-1]], [yl[L-1], yl[2*L-1]], 'r-', alpha=1, linewidth=0.5, zorder=1)
    plt.plot([xl[2*L*(L-1)], xl[2*L*L - L]], [yl[2*L*(L-1)], yl[2*L*L - L]], 'r-', alpha=1, linewidth=0.5, zorder=1)


fig2 = plt.figure()
grid_plot()
plt.scatter(hex_cor(n_hx, i_hx)[0], hex_cor(n_hx, i_hx)[1], c=p_hx, cmap="Greens", s=100, vmin=0, vmax=3, marker="s", zorder=2)
plt.axis("off")
plt.savefig("hex_ex.pdf", dpi=100)
#plt.show()

