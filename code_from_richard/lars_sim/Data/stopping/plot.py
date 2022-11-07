import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
plt.rcParams.update({'font.size': 20})

sq_1 = np.genfromtxt("square_1.txt", delimiter=" ")

sq_1 = sq_1.flatten()

#sq_1 = np.log(sq_1)

hist_sq, bin_sq = np.histogram(sq_1, bins=100000, density=True)

x = np.zeros(len(hist_sq))
for i in range(len(x)):
    x[i] = (bin_sq[i] + bin_sq[i+1])/2

plt.plot(x, hist_sq)
plt.xscale("log")

plt.savefig("stop_square.pdf")