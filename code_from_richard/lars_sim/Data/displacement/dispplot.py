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

print("n?")
n = int(input())

# var and kurt
if (n == 1 and True):
    svk_1 = np.genfromtxt("sq_varkur_rho_05_1.txt", delimiter=" ")
    svk_2 = np.genfromtxt("sq_varkur_rho_25_1.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 200]
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
    plt.savefig("var_05_1.pdf", dpi=150)

    fig1 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".", label="kurt")
    plt.savefig("kur_05_1.pdf", dpi=150)

    # HIGHER DENSITY
    t = svk_2[:, 0]

    tlow = t[t <= 10]
    tint = t[10:100]
    thigh = t[t > 20000]

    par1 = opt.curve_fit(fit, tlow, svk_2[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_2[-len(thigh):, 1])[0]
    par3 = opt.curve_fit(fit, tint, svk_2[10:100, 1])[0]

    fig = plt.figure()
    plt.plot(svk_2[:, 0], svk_2[:, 1], ".", label="var")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]+.1), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(t[-10000:], fit(t[-10000:], 1, par2[1]+.5), "k-", label=r"$\sim t^{%.2f}$" % 1)
    plt.plot(tint, fit(tint, par3[0], par3[1]+1), linewidth=2, color="black", linestyle="dotted", label=r"$\sim t^{%.2f}$" % par3[0])
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("var_25_1.pdf", dpi=150)

    fig1 = plt.figure()
    plt.plot(np.linspace(min(svk_2[:, 0]), max(svk_2[:, 0]), 2), 3*np.ones(2), "k--")
    plt.plot(svk_2[:, 0], svk_2[:, 2], ".", label="kurt")
    plt.savefig("kur_25_1.pdf", dpi=150)

elif (n == 2 and True):
    svk_1 = np.genfromtxt("sq_varkur_rho_05_2.txt", delimiter=" ")
    svk_2 = np.genfromtxt("sq_varkur_rho_25_2.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 200]
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
    plt.savefig("var_05_2.pdf", dpi=150)

    fig1 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".", label="kurt")
    plt.plot([0, 10000], 3*np.ones(2), "k--")
    plt.savefig("kur_05_2.pdf", dpi=150)

    # HIGHER DENSITY
    t = svk_2[:, 0]

    tlow = t[t <= 10]
    tint = t[60:600]
    thigh = t[t > 7000]

    par1 = opt.curve_fit(fit, tlow, svk_2[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_2[-len(thigh):, 1])[0]
    par3 = opt.curve_fit(fit, tint, svk_2[60:600, 1])[0]

    fig = plt.figure()
    plt.plot(svk_2[:, 0], svk_2[:, 1], ".", label="var")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]+.1), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(t[-5000:], fit(t[-5000:], par2[0], par2[1]+.01), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    plt.plot(tint, fit(tint, par3[0], par3[1]+20), linewidth=2, color="black", linestyle="dotted", label=r"$\sim t^{%.2f}$" % par3[0])
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("var_25_2.pdf", dpi=150)

    fig1 = plt.figure()
    plt.plot(np.linspace(min(svk_2[:, 0]), max(svk_2[:, 0]), 2), 3*np.ones(2), "k--")
    plt.plot(svk_2[:, 0], svk_2[:, 2], ".", label="kurt")
    plt.savefig("kur_25_2.pdf", dpi=150)


elif (n == 3 and True):
    svk_1 = np.genfromtxt("sq_varkur_rho_05_3.txt", delimiter=" ")
    svk_2 = np.genfromtxt("sq_varkur_rho_25_3.txt", delimiter=" ")

    t = svk_1[:, 0]

    tlow = t[t <= 200]
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
    plt.savefig("var_05_3.pdf", dpi=150)

    fig1 = plt.figure()
    plt.plot(svk_1[:, 0], svk_1[:, 2], ".", label="kurt")
    plt.plot([0, 10000], 3*np.ones(2), "k--")
    plt.savefig("kur_05_3.pdf", dpi=150)

    # HIGHER DENSITY
    t = svk_2[:, 0]

    tlow = t[t <= 10]
    tint = t[60:600]
    thigh = t[t > 7000]

    par1 = opt.curve_fit(fit, tlow, svk_2[:len(tlow), 1])[0]
    par2 = opt.curve_fit(fit, thigh, svk_2[-len(thigh):, 1])[0]
    par3 = opt.curve_fit(fit, tint, svk_2[60:600, 1])[0]

    fig = plt.figure()
    plt.plot(svk_2[:, 0], svk_2[:, 1], ".", label="var")
    plt.plot(tlow, fit(tlow, par1[0], par1[1]+.1), "k--", label=r"$\sim t^{%.2f}$" % par1[0])
    plt.plot(t[-5000:], fit(t[-5000:], par2[0], par2[1]+.01), "k-.", label=r"$\sim t^{%.2f}$" % par2[0])
    plt.plot(tint, fit(tint, par3[0], par3[1]+20), linewidth=2, color="black", linestyle="dotted", label=r"$\sim t^{%.2f}$" % par3[0])
    plt.legend()
    plt.xscale("log")
    plt.yscale("log")
    plt.savefig("var_25_3.pdf", dpi=150)

    fig1 = plt.figure()
    plt.plot(np.linspace(min(svk_2[:, 0]), max(svk_2[:, 0]), 2), 3*np.ones(2), "k--")
    plt.plot(svk_2[:, 0], svk_2[:, 2], ".", label="kurt")
    plt.savefig("kur_25_3.pdf", dpi=150)

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
# distribution
if (n == 1 and True):
    dist = np.genfromtxt("sq_dist_rho_05_1.txt", delimiter=" ")
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
    plt.savefig("dist_low_05_1.pdf", dpi=150)

    dis = distribution(dist, bins, 6)
    d1 = np.asarray(dis[0])
    x1 = np.asarray(dis[1])
    
    n = len(x1)
    mean = sum(x1*d1)/n
    sigma = sum(d1*(x1-mean)**2)/n 


    par = opt.curve_fit(laplace_distr, xdata=x1, ydata=d1, p0=[0,sigma])[0]

    print(par)

    fig2 = plt.figure()
    plt.plot(x1, laplace_distr(x1, *par), "k-.", label="Laplace fit")
    for i in range(3, 6):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.xlabel(r"$x(t)-x(0)$")
    plt.ylabel(r"$P(x(t)-x(0))$")
    plt.axis([-3000, 3000, 10**(-6), 1])
    plt.legend()
    plt.savefig("dist_high_05_1.pdf", dpi=150)  

    
    dist = np.genfromtxt("sq_dist_rho_25_1.txt", delimiter=" ")
    bins = 300

    fig1 = plt.figure()
    for i in range(3):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.legend()
    plt.savefig("dist_low_25_1.pdf", dpi=150)

    fig2 = plt.figure()
    for i in range(3, 6):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.legend()
    plt.savefig("dist_high_25_1.pdf", dpi=150)
    

if (n == 2 and True):
    
    dist = np.genfromtxt("sq_dist_rho_05_2.txt", delimiter=" ")
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
    plt.savefig("dist_low_05_2.pdf", dpi=150)

    dis = distribution(dist, bins, 6)
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
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.xlabel(r"$x(t)-x(0)$")
    plt.ylabel(r"$P(x(t)-x(0))$")
    plt.axis([-3000, 3000, 10**(-6), 1])
    plt.legend()
    plt.savefig("dist_high_05_2.pdf", dpi=150)

    # HIGHER DENSITY

    dist = np.genfromtxt("sq_dist_rho_25_2.txt", delimiter=" ")
    bins = 300

    fig1 = plt.figure()
    for i in range(3):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.legend()
    plt.savefig("dist_low_25_2.pdf", dpi=150)

    fig2 = plt.figure()
    for i in range(3, 6):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.legend()
    plt.savefig("dist_high_25_2.pdf", dpi=150)

if (n == 3 and True):
    dist = np.genfromtxt("sq_dist_rho_05_3.txt", delimiter=" ")
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
    plt.savefig("dist_low_05_3.pdf", dpi=150)

    dis = distribution(dist, bins, 6)
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
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.xlabel(r"$x(t)-x(0)$")
    plt.ylabel(r"$P(x(t)-x(0))$")
    plt.axis([-3000, 3000, 10**(-6), 1])
    plt.legend()
    plt.savefig("dist_high_05_3.pdf", dpi=150)  

    # HIGHER DENSITY

    dist = np.genfromtxt("sq_dist_rho_25_3.txt", delimiter=" ")
    bins = 300

    fig1 = plt.figure()
    for i in range(3):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.legend()
    plt.savefig("dist_low_25_3.pdf", dpi=150)

    fig2 = plt.figure()
    for i in range(3, 6):
        dist1 = distribution(dist, bins, i+1)
        d1 = np.asarray(dist1[0])
        x1 = np.asarray(dist1[1])
        plt.plot(x1[d1>0], d1[d1>0], distinguisher[i-3], label=name[i], linewidth=1)
    plt.yscale("log")
    plt.legend()
    plt.savefig("dist_high_25_3.pdf", dpi=150)