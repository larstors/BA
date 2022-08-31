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
L = [100, 100]
alpha = []
N = 600
print("n?")
n = int(input())


colormap = cm.inferno


# triangular
if True:
    f = open("triangle_2.txt")
    data = []
    for line in f:
        data.append(np.loadtxt(StringIO(line), dtype=int))
        


    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 20))
    
    def conv(n):

        x = []
        y = []
        for i in n:
            x.append(i % L[0] - 0.5 * int(i / L[0]))
            y.append(int(i/L[0]) * np.sqrt(3)/2)
        
        return x, y
        


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


    
    def animate_tri(i):
        ax.cla()
        ax.triplot(triang, 'r-', alpha=0.3, linewidth=0.1)
        ax.scatter(x=conv(data[i][1::3])[0], y=conv(data[i][1::3])[1], c=data[i][2::3], cmap="Greens", s=20, vmin=0, vmax=n, marker="s")
        ax.axis([-50, 110, -5, 90])
        ax.set_title(r"$N=8000, \alpha=0.3, n_\mathrm{max}=2$")
    
        



    # Then setup FuncAnimation.
    ani_tri = animation.FuncAnimation(fig, animate_tri, interval=1, blit=False, save_count=399)
    ani_tri.save('tri_2_N_8000.gif', writer='PillowWriter', fps=10)
