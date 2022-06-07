
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

data = np.genfromtxt("./Data/motility/hexagonal.txt")

plt.plot(data[:, 0], data[:, 1], "rs", label="J")
plt.plot(data[:, 0], data[:, 2], "bo", label="M")
plt.show()
print(0.2*10*10*2*3)