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

plt.rcParams.update({'font.size': 18})

n = 1

def distribution(array, bin, index):
    x = []
    t, h = np.histogram(array[index, :], bins=bin, density=True)
    for k in range(len(t)):
        x.append((h[k+1]+h[k])/2)
    
    return t, x

def gauss_function(x, x0, sigma):
    return 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-x0)**2/(2*sigma**2))

def laplace_distr(x, mu, sigma):
    return 1/np.sqrt(2*sigma**2)*np.exp(-np.abs(x-mu)*np.sqrt(2)/sigma)

name = [r"$t=50$", r"$t=150$", r"$t=200$", r"$t=500$", r"$t=5000$", r"$t=25000$"]
distinguisher = ["r-", "g-", "b-"]

# var and kurt
if (n == 1 and False):
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


    svk_1 = np.genfromtxt("sq_varkur_rho_25_1_test.txt", delimiter=" ")
    #svk_2 = np.genfromtxt("sq_varkur_rho_25_1.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 10]
    thigh = t[t > 80000]

    def fit (x, a, b):
        return x ** a * b

    par1 = opt.curve_fit(fit, tlow, svk_1[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_1[-len(thigh):, 1])[0]

    fig = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 1], ".", label="var")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(thigh, fit(thigh, par2[0], par2[1]), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    #plt.axis([10000, 100000, 1000, 1000000])
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.show()

    fig1 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".", label="kurt")
    plt.show()

if False:
    svk_1 = np.genfromtxt("sq_varkur_rho_25_1_test_test_test.txt", delimiter=" ")
    #svk_2 = np.genfromtxt("sq_varkur_rho_25_1.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 10]
    thigh = t[t > 2000]

    def fit (x, a, b):
        return x ** a * b

    par1 = opt.curve_fit(fit, tlow, svk_1[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_1[-len(thigh):, 1])[0]

    fig = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 1], ".", label="var")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(thigh, fit(thigh, par2[0], par2[1]), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    #plt.axis([10000, 100000, 1000, 1000000])
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    #plt.savefig("L_100_25_3.pdf")
    plt.show()


if False:
    x = np.linspace(1, 1000, 1000)
    def sq(x):
        return x**2
    
    def pos(x):
        return x**2 + x

    plt.plot(x, sq(x), label="square")
    plt.plot(x, pos(x), label="possion")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid()
    plt.show()

if False:
    svk_1 = np.genfromtxt("sq_varkur_rho_25_2_test_test.txt", delimiter=" ")
    #svk_2 = np.genfromtxt("sq_varkur_rho_25_1.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 10]
    thigh = t[t > 20000]
    t3 = t[10000:30000]

    def fit (x, a, b):
        return x ** a * b

    par1 = opt.curve_fit(fit, tlow, svk_1[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_1[-len(thigh):, 1])[0]
    par3 = opt.curve_fit(fit, t3, svk_1[10000:30000, 1])[0]

    fig = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 1], ".", label="var")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(thigh, fit(thigh, par2[0], par2[1]), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    plt.plot(t3, fit(t3, par3[0], par3[1]), "k-", label=r"$\sim t^{%.2f}$" % par3[0])
    #plt.axis([10000, 100000, 1000, 1000000])
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    #plt.savefig("L_100_25_3.pdf")
    plt.show()


if True:
    svk_1 = np.genfromtxt("sq_varkur_rho_25_3_test2.txt", delimiter=" ")
    dist = np.genfromtxt("sq_dist_rho_25_3_test2.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 10]
    thigh = t[-40000:]
    t3 = t[200:1000]

    def fit (x, a, b):
        return x ** a * b

    par1 = opt.curve_fit(fit, tlow, svk_1[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_1[-len(thigh):, 1])[0]
    par3 = opt.curve_fit(fit, t3, svk_1[200:1000, 1])[0]

    fig = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 1], "m.")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(thigh, fit(thigh, par2[0], par2[1]), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    plt.plot(t3, fit(t3, par3[0], par3[1]), "k-", label=r"$\sim t^{%.2f}$" % par3[0])
    #plt.axis([10000, 100000, 1000, 1000000])
    plt.legend()
    plt.grid()
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\Delta x(t)^2$")
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("var_25_3.pdf", bbox_inches="tight")

    fig3 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".")
    plt.grid()
    plt.xlabel(r"$t$")
    plt.ylabel(r"$\kappa(t)$")
    plt.savefig("kur_25_3.pdf", bbox_inches="tight")

    bins = 300

    fig1 = plt.figure()
    for i in range(3):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.xlabel(r"$x(t)-x(0)$")
    plt.ylabel(r"$P(x(t)-x(0))$")
    plt.legend()
    plt.savefig("dist_low_25_3.pdf", dpi=150, bbox_inches="tight")
    

    dis = distribution(dist, bins, 5)
    d1 = np.asarray(dis[0])
    x1 = np.asarray(dis[1])
    
    n = len(x1)
    mean = sum(x1*d1)/n
    sigma = sum(d1*(x1-mean)**2)/n 



    par = opt.curve_fit(gauss_function, xdata=x1, ydata=d1, p0=[mean,sigma])[0]

    print(par)

    fig2 = plt.figure()
    plt.plot(x1, gauss_function(x1, *par), "k-.", label="Gaussian fit")
    for i in range(3, 6):
        dist1 = distribution(dist, bins, i)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.xlabel(r"$x(t)-x(0)$")
    plt.ylabel(r"$P(x(t)-x(0))$")
    plt.axis([-1000, 1000, 10**(-6), 1])
    plt.legend()
    plt.savefig("dist_high_25_3.pdf", dpi=150, bbox_inches="tight")