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
L = [100, 100]
alpha = []
N = 600





class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50):
        f = open("./triangle.txt")
        self.data = []
        for line in f:
            self.data.append(np.loadtxt(StringIO(line), dtype=int))
        
        self.triang = self.trian()
        # Setup the figure and axes
        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        # Then setup FuncAnimation.
        #self.ani = animation.FuncAnimation(self.fig, self.update, interval=1, 
        #                                  init_func=self.setup_plot, blit=True)
        

    def conv(self, n):
        """function that converts the index n into coordinates x and y, with x being horisontal. Assuming equilateral triangle of sides 1.

        Args:
            n (int): index of particle in lattice

        Returns:
            tuple: x and y coordinates
        """
        x = []
        y = []
        for i in n:
            x.append(i % L[0] - 0.5 * int(i / L[0]))
            y.append(int(i/L[0]) * np.sqrt(3)/2)
        
        return x, y
    
    def trian(self):
        # indices for entire lattice
        ind = np.array([i for i in range(L[0]*L[1])])
        # x and y coordinate for each lattice point
        x = self.conv(ind)[0]
        y = self.conv(ind)[1]

        triangles = []

        for j in range(L[1]-1):
            for i in range(L[0]-1):
                i = i + j * L[0]
                triangles.append([i, i+1, i+L[0]+1])
                triangles.append([i+L[0], i, i+1+L[0]])
        triang = mtri.Triangulation(x, y, triangles)
        return triang

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        s = 1.8 * 10
        kos = self.data[0][1::3]
        weight = self.data[0][2::3]
        self.scat = self.ax.scatter(x=self.conv(kos)[0], y=self.conv(kos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=2, marker="s")
        self.ax.axis([-0.5*L[0], L[0],-5, L[1]])
        self.ax.triplot(self.triang, 'r-', alpha=0.6, linewidth=0.2)
        return self.scat,


    def update(self, i):
        """Update the scatter plot."""
        kos = self.data[i][1::3]
        self.ax.cla()
        s = 1.8 
        weight = self.data[i][2::3]
        self.scat = self.ax.scatter(x=self.conv(kos)[0], y=self.conv(kos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=2, marker="s")
        self.ax.axis([-0.5*L[0], L[0],-5, L[1]])
        self.ax.triplot(self.triang, 'r-', alpha=0.6, linewidth=0.2)
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    #a.ani.save('TriangleScatter.gif', fps=10)
    
    a.update(50)
    plt.show()
    #plt.savefig("./lars_sim/triang_test.pdf")


