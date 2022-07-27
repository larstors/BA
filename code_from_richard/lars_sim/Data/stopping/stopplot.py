import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numba as nb
from matplotlib import rcParams
plt.rcParams.update({'font.size': 35})


print("Lattice? (only first letter)")
lat = input()

# triangle
if lat == "t":
    # n=1
    d1_001 = np.genfromtxt("tri_1_0.001.txt", delimiter=" ")
    d1_01 = np.genfromtxt("tri_1_0.01.txt", delimiter=" ")
    d1_1 = np.genfromtxt("tri_1_0.1.txt", delimiter=" ")

    # n=2
    d2_001 = np.genfromtxt("tri_2_0.001.txt", delimiter=" ")
    d2_01 = np.genfromtxt("tri_2_0.01.txt", delimiter=" ")
    d2_1 = np.genfromtxt("tri_2_0.1.txt", delimiter=" ")

    # n=3
    d3_001 = np.genfromtxt("tri_3_0.001.txt", delimiter=" ")
    d3_01 = np.genfromtxt("tri_3_0.01.txt", delimiter=" ")
    d3_1 = np.genfromtxt("tri_3_0.1.txt", delimiter=" ")


    d1_001 = np.flatten(np.log(d1_001))
    d1_001, h1_001 = np.histogram(d1_001, bins=1000, density=True)
    x1_001 = []
    for i in range(len(d1_001)):
        x1_001.append((h1_001[i+1]+h1_001[i])/2)
    

    plt.plot(np.exp(x1_001), d1_001, "o")
    plt.xscale("log")
    plt.show()

