"""
This is an addition by Lars. 

Here I attempt to read and work with the output of the positional output of the triangular lattice.

"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import scipy
import matplotlib.animation as animation
#from celluloid import Camera
import matplotlib.tri as mtri

# TODO: better way of doing this...
# parameters
L = [10, 10]
alpha = []
N = 600



data = np.loadtxt("./hexagonal.txt")



class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50):
        self.data = np.loadtxt("./hexagonal.txt")
        # Setup the figure and axes
        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
                                          init_func=self.setup_plot, blit=True)
        
    
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

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        s = 1.8 * 100 / L[0]
        kos = self.data[0][1::3] # array with n
        jos = self.data[0][2::3] # array with j
        self.scat = self.ax.scatter(x=self.hex_cor(kos, jos)[0], y=self.hex_cor(kos, jos)[1], c="k", s=s, vmin=0, vmax=1, marker="s", edgecolor="k")
        self.ax.axis([-0.5*L[0], L[0],-5, L[1]])

        return self.scat,


    def update(self, i):
        """Update the scatter plot."""
        kos = self.data[i][1::3] # array with n
        jos = self.data[i][2::3] # array with j
        self.ax.cla()
        s = 1.8 * 100 / L[0]
        self.scat = self.ax.scatter(x=self.hex_cor(kos, jos)[0], y=self.hex_cor(kos, jos)[1], c="k", s=s, vmin=0, vmax=1, marker="s", edgecolor="k")
        self.ax.axis([-0.5*L[0], L[0],-5, L[1]])

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    a.ani.save('HexagonalScatter.gif', fps=5)


