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



data = np.loadtxt("./triangle.txt")



class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, numpoints=50):
        self.data = np.loadtxt("./triangle.txt")
        self.triang = self.trian()
        # Setup the figure and axes
        self.fig, self.ax = plt.subplots(figsize=(20, 20))
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=5, 
                                          init_func=self.setup_plot, blit=True)
        

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
        s = 1.8 * 100 / L[0]
        kos = self.data[0][1::3]
        weight = self.data[0][2::3]
        self.scat = self.ax.scatter(x=self.conv(kos)[0], y=self.conv(kos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=1, marker="s", edgecolor="k")
        self.ax.axis([-0.5*L[0], L[0],-5, L[1]])
        self.ax.triplot(self.triang, 'r-', alpha=0.6, linewidth=0.2)
        return self.scat,


    def update(self, i):
        """Update the scatter plot."""
        kos = self.data[i][1::3]
        self.ax.cla()
        s = 1.8 * 100 / L[0]
        weight = self.data[i][2::3]
        self.scat = self.ax.scatter(x=self.conv(kos)[0], y=self.conv(kos)[1], c=weight, cmap="Greens", s=s, vmin=0, vmax=1, marker="s", edgecolor="k")
        self.ax.axis([-0.5*L[0], L[0],-5, L[1]])
        self.ax.triplot(self.triang, 'r-', alpha=0.6, linewidth=0.2)
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat,


if __name__ == '__main__':
    a = AnimatedScatter()
    a.ani.save('scatter.gif', fps=5)



"""
# creating figure
fig = plt.figure()
ax = plt.axes(xlim=(-50, 105), ylim=(-5, 100))
line, = ax.plot([], [])


# animation function that will be called sequentially
def animate(i):
"""
"""camera = Camera(plt.figure())
for i in range(len(data)):
    #positions 
    kos = data[i][1::2]
    #plot them
    plt.scatter(x=conv(kos)[0], y=conv(kos)[1], s=1.8, marker="s", c="k")
    #snapshot for movie
    camera.snap()
anim = camera.animate(blit=True)
anim.save('scatter.gif')
"""

"""
k = data[0]


plt.savefig("triangle_plot_test.pdf", dpi=400)"""