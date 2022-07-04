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


print("lattice size?")
l = int(input())
L = l
alpha = []
N = 600

print("Occupation number?")
n = int(input())

print("Lattice?")
output = input()

print("Frame?")
f = int(input())

print(f)

def conv(n):
        
        x = []
        y = []
        for i in n:
            x.append(i % L)
            y.append(int(i/L))
        
        return x, y

def tri_conv(n):
    
        x = []
        y = []
        for i in n:
            x.append(i % L - 0.5 * int(i / L))
            y.append(int(i/L) * np.sqrt(3)/2)
        
        return x, y

def hex_conv(n, j):
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

if output[0] == "s":
    # square
    f = open("square.txt")
    square_part = []
    for line in f:
        square_part.append(np.loadtxt(StringIO(line), dtype=int))




elif output[0] == "t":
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

    f = open("triangle.txt")
    square_part = []
    for line in f:
        square_part.append(np.loadtxt(StringIO(line), dtype=int))

elif output[0] == "h":
    #hexagonal
    f = open("hexagonal.txt")
    square_part = []
    for line in f:
        square_part.append(np.loadtxt(StringIO(line), dtype=int))

    fig, ax = plt.subplots(nrows=1, ncols=1)
    def hex_cor_int(n, j):
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
                ax.plot(xs, ys, 'r-', alpha=0.6, linewidth=0.2)

        ax.plot([xl[L-1], xl[2*L-1]], [yl[L-1], yl[2*L-1]], 'r-', alpha=0.6, linewidth=0.2)
        ax.plot([xl[2*L*(L-1)], xl[2*L*L - L]], [yl[2*L*(L-1)], yl[2*L*L - L]], 'r-', alpha=0.6, linewidth=0.2)


    ax.scatter(x=hex_conv(square_part[0][1::4], square_part[0][2::4])[0], y=hex_conv(square_part[0][1::4], square_part[0][2::4])[1], c=square_part[0][3::4], marker="s", s=16, cmap="Greens", vmin=0, vmax=n)
    grid_plot()
    plt.show()



