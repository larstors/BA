from wsgiref.headers import tspecials
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import time
import random as rn
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import matplotlib.patches as mpl_patches
from matplotlib import rcParams
from io import StringIO
import matplotlib.tri as mtri
import matplotlib.patches as mpatches

plt.rcParams.update({'font.size': 15})


def lin_fit(x, a, b):
    return b*x**a


def plot(rho, n):
    low = np.genfromtxt("sq_varkur_rho_"+str(rho)+"_"+str(n)+"_alpha0.001.txt", delimiter=" ")
    med = np.genfromtxt("sq_varkur_rho_"+str(rho)+"_"+str(n)+"_alpha0.010.txt", delimiter=" ")
    hig = np.genfromtxt("sq_varkur_rho_"+str(rho)+"_"+str(n)+"_alpha0.100.txt", delimiter=" ")

    hig_par = opt.curve_fit(lin_fit, hig[:10, 0], hig[:10, 1], p0=(1, 1))[0]

    plt.plot(low[:, 0], low[:, 1], label=r"$\alpha = 10^{-3}$")
    plt.plot(med[:, 0], med[:, 1], label=r"$\alpha = 10^{-2}$")
    plt.plot(hig[:, 0], hig[:, 1], label=r"$\alpha = 10^{-1}$")
    plt.plot(hig[:10, 0], lin_fit(hig[:10, 0], *hig_par), "--", label="Linear fit, a = %g" % hig_par[0])
    plt.legend()
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\langle\Delta x(t)^2\rangle$")
    plt.grid()
    plt.xscale("log")
    plt.yscale("log")
    plt.axis([0.9, 1000, 0.1, 100])
    plt.savefig("msd_comp_sq_n_"+str(n)+"_rho_"+str(rho)+".svg", dpi=200, bbox_inches="tight")




plot(65, 1)


