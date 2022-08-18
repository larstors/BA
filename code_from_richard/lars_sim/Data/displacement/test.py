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

plt.rcParams.update({'font.size': 14})

n = 1

def distribution(array, bin, index):
    x = []
    t, h = np.histogram(array[:, index] - array[:, 0], bins=bin, density=True)
    for k in range(len(t)):
        x.append((h[k+1]+h[k])/2)
    
    return t, x

def gauss_function(x, x0, sigma):
    return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-x0)**2/(2*sigma**2))

def laplace_distr(x, mu, sigma):
    return 1/np.sqrt(2*sigma**2)*np.exp(-np.abs(x-mu)*np.sqrt(2)/sigma)

name = [r"$t=50$", r"$t=150$", r"$t=200$", r"$t=500$", r"$t=1000$", r"$t=3000$"]
distinguisher = ["r-", "g-", "b-"]

# var and kurt
if (n == 1 and True):
    svk_1 = np.genfromtxt("sq_varkur_rho_05_1_test.txt", delimiter=" ")
    #svk_2 = np.genfromtxt("sq_varkur_rho_25_1.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 150]
    thigh = t[t > 300]

    def fit (x, a, b):
        return x ** a * b

    par1 = opt.curve_fit(fit, tlow, svk_1[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_1[-len(thigh):, 1])[0]

    fig = plt.figure()
    plt.plot(tlow, fit(tlow, par1[0], par1[1] + 0.2), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(thigh, fit(thigh, par2[0], par2[1] + 10), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    plt.plot(svk_1[:, 0], svk_1[:, 1], ".", label="var")
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

    fig1 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".", label="kurt")
    plt.show()

    fig2 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 1] - svk_1[:, 3])
    plt.grid()
    plt.show()


    svk_1 = np.genfromtxt("sq_varkur_rho_05_2_long.txt", delimiter=" ")
    #svk_2 = np.genfromtxt("sq_varkur_rho_25_1.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 150]
    thigh = t[t > 300]

    def fit (x, a, b):
        return x ** a * b

    par1 = opt.curve_fit(fit, tlow, svk_1[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_1[-len(thigh):, 1])[0]

    fig = plt.figure()
    plt.plot(tlow, fit(tlow, par1[0], par1[1] + 0.2), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(thigh, fit(thigh, par2[0], par2[1] + 10), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    plt.plot(svk_1[:, 0], svk_1[:, 1], ".", label="var")
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

    fig1 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".", label="kurt")
    plt.show()
