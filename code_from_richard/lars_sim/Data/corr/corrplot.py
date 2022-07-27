from telnetlib import TM
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from numba import *
from matplotlib import rcParams
plt.rcParams.update({'font.size': 35})


sq_data = np.genfromtxt("sq....", delimiter=" ")

L = 100
N = 5800
n = 3
rho2 = N**2/(L**4*n**2)
K = 100000

kappa0 = 1/(K*L**2*n**2)*np.sum(sq_data[:K]*sq_data[:K])

Tmax = 1000000
corr = np.ones(Tmax)

t = np.arange(0, Tmax)*0.1

@njit(parallel=True, nopython=True)
def calc():
    for t in range(Tmax):
        psi = 1/(K*L**2*n**2)*np.sum(sq_data[:K]*sq_data[t:(K+t)]) - rho2
        corr[t] = psi / (kappa0 - rho2)
    



# if we do 0.1 step size and have up to 10^6, we get 10^7 steps, but we need to average over a time.
# if b=1 000 000 and u=1 100 000, we have a total of 100 000 extra steps to average over, this should work

