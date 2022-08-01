from telnetlib import TM
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from numba import *
from matplotlib import rcParams
from scipy.special import i0e
plt.rcParams.update({'font.size': 35})


tr_data = np.genfromtxt("tri_1_1.000.txt", delimiter=" ")

L = 100
N = 4000
n = 1
rho2 = N**2/(L**4*n**2)
K = 1001

kappa0 = 1/(K*L**2*n**2)*np.sum(tr_data[:K]*tr_data[:K])

Tmax = 10000
corr = np.ones(Tmax)

t = np.arange(0, Tmax)*100

#@njit(parallel=True, nopython=True)
def calc():
    for t in range(Tmax):
        psi = 1/(K*L**2*n**2)*np.sum(tr_data[:K]*tr_data[t:(K+t)]) - rho2
        corr[t] = psi / (kappa0 - rho2)

    
calc()

plt.plot(t, corr)
plt.plot(t, i0e(t))
plt.plot(t, corr/i0e(t))
plt.xscale("log")
plt.yscale("log")
plt.savefig("1.pdf")



# if we do 0.1 step size and have up to 10^6, we get 10^7 steps, but we need to average over a time.
# if b=1 000 000 and u=1 100 000, we have a total of 100 000 extra steps to average over, this should work

