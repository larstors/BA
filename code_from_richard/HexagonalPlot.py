"""
This is an addition by Lars. 

Here I attempt to read and work with the output of the positional output of the triangular lattice.

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
n_max = int(input())






class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50):
        # opening file (note the annoying way to do it because of different number of cols)
        f = open("./lars_sim/testing/hexagonal.txt")
        self.data = []
        for line in f:
            self.data.append(np.loadtxt(StringIO(line), dtype=int))
        
        # Setup the figure and axes
        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        # Then setup FuncAnimation.
        #self.ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
        #                                  init_func=self.setup_plot, blit=True, save_count=1000)
        
    
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
        self.scat = self.ax.scatter(x=self.hex_cor(kos, jos)[0], y=self.hex_cor(kos, jos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=2, marker="s")
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
                plt.plot(xs, ys, 'r-', alpha=0.6, linewidth=0.2)

        plt.plot([xl[L[0]-1], xl[2*L[0]-1]], [yl[L[0]-1], yl[2*L[0]-1]], 'r-', alpha=0.6, linewidth=0.2)
        plt.plot([xl[2*L[0]*(L[0]-1)], xl[2*L[0]*L[0] - L[0]]], [yl[2*L[0]*(L[0]-1)], yl[2*L[0]*L[0] - L[0]]], 'r-', alpha=0.6, linewidth=0.2)
    
    def ani_particle(self):
        ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
                                          init_func=self.setup_plot, blit=True)
        return ani
    
    def ani_direction(self):
        ani = animation.FuncAnimation(self.fig, self.direction, interval=5, 
                                          init_func=self.setup_direction, blit=True)
        return ani

    def update(self, i):
        """Update the scatter plot."""
        kos = self.data[i][1::4] # array with n
        jos = self.data[i][2::4] # array with j

        weight = self.data[i][3::4]

        self.ax.cla()
        s = 1.8 * 100
                
        self.scat = self.ax.scatter(x=self.hex_cor(kos, jos)[0], y=self.hex_cor(kos, jos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=2, marker="s")
        self.grid_plot()
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,

    def setup_direction(self):
        # working with 1 particle as of now
        self.grid_plot()
        dir = np.loadtxt("hexdir.txt", delimiter="\t", dtype=int)
        n = dir[0,0]
        j = dir[0,1]
        m = dir[0,2]
        x, y = self.hex_cor_int(n, j)
        origin = np.array([[x], [y]])
        dx, dy = self.hex_cor_int(m, (j+1)%2)
        V = np.array([[dx-x, dy-y]])
        plt.quiver(*origin, V[0, 0], V[0, 1])
        scat = plt.scatter(x=x, y=y, cmap="Greens")
        return scat,


    def direction(self, i):
        # working with 1 particle as of now
        plt.clf()
        self.grid_plot()
        dir = np.loadtxt("hexdir.txt", delimiter="\t", dtype=int)
        n = dir[i,0]
        j = dir[i,1]
        m = dir[i,2]
        x, y = self.hex_cor_int(n, j)
        origin = np.array([[x], [y]])
        origin2 = np.array([[x], [y]])
        
        dx, dy = self.hex_cor_int(m, (j+1)%2)
        V = np.array([[dx-x, dy-y]])
        plt.arrow(x, y, V[0, 0], V[0, 1], width=0.1, length_includes_head=True, head_width=0.15, head_length=0.2, color="black")
        #plt.quiver(*origin, V[0, 0], V[0, 1])
        #plt.quiver(*origin2, V[0, 0], V[0, 1], color="red")
        scat = plt.scatter(x=x, y=y, cmap="Greens", s=100)
        return scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    #a.ani.save('./lars_sim/gif/HexStart_n_2_N_750.gif', fps=20)
    #a.direction(0)
    a.update(3)
    plt.show()
    #plt.savefig("hexdir2.pdf", dpi=200)
    #animation = a.ani_direction()
    #animation.save("./lars_sim/HexagonalDirection.gif", fps=10)
    




