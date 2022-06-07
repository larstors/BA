"""
This is an addition by Lars. 


"""


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


colormap = cm.inferno


# TODO: better way of doing this...
# parameters
l = int(input())
L = [l, l]
alpha = []
N = 600
n = int(input())






class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50):
        # opening file (note the annoying way to do it because of different number of cols)
        f = open("./Data/stable/hexagonal.txt")
        self.data = []
        for line in f:
            self.data.append(np.loadtxt(StringIO(line), dtype=int))
        
        self.numb = np.genfromtxt("./Data/stable/hexagonal_number.txt", delimiter=" ")
        # Setup the figure and axes
        self.fig, self.ax = plt.subplots(nrows=2, ncols=1, figsize=(20, 20), gridspec_kw={'height_ratios': [15, 2]})
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
                                          init_func=self.setup_plot, blit=True, save_count=len(self.numb)-1)
        
    
    def hex_cor(self, n, j):
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
            y = int(n[k]/L[0])
            x = 2 * (n[k]%L[0]) + j[k] - y
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

    def hex_cor_int(self, n, j):
        xf = 0
        yf = 0

        c30 = np.cos(np.pi/6)
        s30 = np.sin(np.pi/6)

        
        y = int(n/L[0])
        x = 2 * (n%L[0]) + j - y
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
    
    def fill(self, xl, yl):
        for i in range(L[0]):
            for k in range(2):
                f = 1 - k
                for j in range(L[0]):  
                    n = i * L[0] + j

                    if (f%2==0):
                        xl.append(self.hex_cor_int(n, f)[0])
                        yl.append(self.hex_cor_int(n, f)[1])
                    else:
                        xl.append(self.hex_cor_int(n, f)[0])
                        yl.append(self.hex_cor_int(n, f)[1])
        return xl, yl


    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        s = 1.8 * 100
        kos = self.data[0][1::4] # array with n
        jos = self.data[0][2::4] # array with j
        weight = self.data[0][3::4]
        self.scat = self.ax[0].scatter(x=self.hex_cor(kos, jos)[0], y=self.hex_cor(kos, jos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=2, marker="s")
        #self.ax.axis([-0.5*L[0], L[0],-5, L[1]])

        return self.scat,

    def grid_plot(self):

        xl = []
        yl = []
        xl, yl = self.fill(xl, yl)

        for n in range(0, L[0]-1, 1):
            for i in range(L[0]-1):
                coord = [[xl[i + 2*n*L[0]], yl[i+ 2*n*L[0]]], [xl[i + L[0] + 1 + 2*n*L[0]], yl[i+ L[0] + 1 + 2*n*L[0]]], [xl[i+ 2*L[0] + 1+ 2*n*L[0]], yl[i+ 2*L[0] + 1+ 2*n*L[0]]], [xl[i+ 3*L[0] + 1+ 2*n*L[0]], yl[i+ 3*L[0] +1+ 2*n*L[0]]], [xl[i+ L[0]*2+ 2*n*L[0]], yl[i+ L[0]*2+ 2*n*L[0]]], [xl[i+ L[0]+ 2*n*L[0]], yl[i+ L[0]+ 2*n*L[0]]], [xl[i+ 2*n*L[0]], yl[i+ 2*n*L[0]]]]
                xs, ys = zip(*coord)
                self.ax[0].plot(xs, ys, 'r-', alpha=0.6, linewidth=0.2)

        self.ax[0].plot([xl[L[0]-1], xl[2*L[0]-1]], [yl[L[0]-1], yl[2*L[0]-1]], 'r-', alpha=0.6, linewidth=0.2)
        self.ax[0].plot([xl[2*L[0]*(L[0]-1)], xl[2*L[0]*L[0] - L[0]]], [yl[2*L[0]*(L[0]-1)], yl[2*L[0]*L[0] - L[0]]], 'r-', alpha=0.6, linewidth=0.2)
    
    def update(self, i):
        """Update the scatter plot."""
        kos = self.data[i][1::4] # array with n
        jos = self.data[i][2::4] # array with j

        weight = self.data[i][3::4]

        self.ax[0].cla()
        s = 20
                
        self.scat = self.ax[0].scatter(x=self.hex_cor(kos, jos)[0], y=self.hex_cor(kos, jos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=n, marker="s")
        #self.grid_plot() # ! THIS LINE TAKES FOREVER
        self.ax[0].axis([-100, 190, -10, 160])
        if i%2:
            self.ax[1].plot(self.numb[i, 0], self.numb[i, 1], "ko", markersize=1)
        self.ax[1].axis([min(self.numb[:, 0]), max(self.numb[:, 0]), min(self.numb[:, 1]), max(self.numb[:, 1])])
        self.ax[1].set_yscale("log")
        self.ax[1].set_ylabel("Number of Cluster")
        self.ax[1].set_xlabel(r"Time $t$")
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

# square
if True:
    f = open("./Data/stable/square.txt")
    data = []
    for line in f:
        data.append(np.loadtxt(StringIO(line), dtype=int))
        
    numb = np.genfromtxt("./Data/stable/square_number.txt", delimiter=" ")


    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(20, 20), gridspec_kw={'height_ratios': [15, 2]})

    def conv(n):

        x = []
        y = []
        for i in n:
            x.append(i % L[0])
            y.append(int(i/L[0]))
        
        return x, y
        
    def animate_tri(i):
        ax[0].cla()
        ax[0].scatter(x=conv(data[i][1::3])[0], y=conv(data[i][1::3])[1], c=data[i][2::3], cmap="Greens", s=20, vmin=0, vmax=n, marker="s")
        ax[0].axis([-50, 110, -5, 90])
        if i%2:
            ax[1].plot(numb[i, 0], numb[i, 1], "ko", markersize=1)
        ax[1].axis([min(numb[:, 0]), max(numb[:, 0]), min(numb[:, 1]), max(numb[:, 1])])
        ax[1].set_yscale("log")
        ax[1].set_ylabel("Number of Cluster")
        ax[1].set_xlabel(r"Time $t$")
        



    # Then setup FuncAnimation.
    ani_tri = animation.FuncAnimation(fig, animate_tri, interval=1, blit=False, save_count=len(numb)-1)
    ani_tri.save('./gif/Sq_stable_n_2_alpha_0.1.gif', writer='PillowWriter', fps=10)

# triangular
if False:
    f = open("./Data/stable/triangular.txt")
    data = []
    for line in f:
        data.append(np.loadtxt(StringIO(line), dtype=int))
        
    numb = np.genfromtxt("./Data/stable/triangular_number.txt", delimiter=" ")


    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(20, 20), gridspec_kw={'height_ratios': [15, 2]})

    def trian():
        # indices for entire lattice
        ind = np.array([i for i in range(L[0]*L[1])])
        # x and y coordinate for each lattice point
        x = conv(ind)[0]
        y = conv(ind)[1]

        triangles = []

        for j in range(L[1]-1):
            for i in range(L[0]-1):
                i = i + j * L[0]
                triangles.append([i, i+1, i+L[0]+1])
                triangles.append([i+L[0], i, i+1+L[0]])
        triang = mtri.Triangulation(x, y, triangles)
        return triang
    
    triang = trian()


    def conv(n):

        x = []
        y = []
        for i in n:
            x.append(i % L[0] - 0.5 * int(i / L[0]))
            y.append(int(i/L[0]) * np.sqrt(3)/2)
        
        return x, y
        
    def animate_tri(i):
        ax[0].cla()
        ax[0].triplot(triang, 'r-', alpha=0.3, linewidth=0.1)
        ax[0].scatter(x=conv(data[i][1::3])[0], y=conv(data[i][1::3])[1], c=data[i][2::3], cmap="Greens", s=20, vmin=0, vmax=n, marker="s")
        ax[0].axis([-50, 110, -5, 90])
        if i%2:
            ax[1].plot(numb[i, 0], numb[i, 1], "ko", markersize=1)
        ax[1].axis([min(numb[:, 0]), max(numb[:, 0]), min(numb[:, 1]), max(numb[:, 1])])
        ax[1].set_yscale("log")
        ax[1].set_ylabel("Number of Cluster")
        ax[1].set_xlabel(r"Time $t$")
        



    # Then setup FuncAnimation.
    ani_tri = animation.FuncAnimation(fig, animate_tri, interval=1, blit=False, save_count=len(numb)-1)
    ani_tri.save('./gif/Tri_stable_n_2_alpha_0.1.gif', writer='PillowWriter', fps=10)




if __name__ == '__main__':
    a = AnimatedScatter()
    a.ani.save('./gif/Hex_stable_n_3_alpha_0.001_long.gif', fps=10)
    #a.update(0)
    #plt.show()
    #plt.savefig("./lars_sim/triang_test.pdf")
    a = 1


