import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numba as nb
from matplotlib import rcParams
plt.rcParams.update({'font.size': 30})


print("Lattice? (only first letter)")
lat = input()


bin= 1000

def stopping(x, f):
    """function to calculate the stopping time, i.e. T_stop=sum_t t*s(t)

    Args:
        x (array): time
        f (array): frequency
    """
    #x = np.exp(x*np.log(10))

    return np.dot(x, f)

def time_stop(name, alpha):
    a = []
    dd = []
    xx = []
    for i in range(len(name)):
        d = np.genfromtxt(name[i], delimiter=" ")
        #d = np.log10(d[d>0]).flatten()
        d, h = np.histogram(d, bins=bin, density=True)
        x = []
        for k in range(len(d)):
            x.append((h[k+1]+h[k])/2)
        x = np.array(x)
        dx = x[1]-x[0] # should be same length between each element
        
        a.append(stopping(x, dx*d))
        dd.append(d)
        xx.append(x)
    return dd, xx, a

def fit(x, a, b):
    return x*a + b

def p(a, b):
    return a**b

# triangle
if lat == "t":
    n1 = ["tri_1_0.001.txt", "tri_1_0.005.txt", "tri_1_0.010.txt", "tri_1_0.050.txt", "tri_1_0.100.txt"]
    n2 = ["tri_2_0.001.txt", "tri_2_0.005.txt", "tri_2_0.010.txt", "tri_2_0.050.txt", "tri_2_0.100.txt"]
    n3 = ["tri_3_0.001.txt", "tri_3_0.005.txt", "tri_3_0.010.txt", "tri_3_0.050.txt", "tri_3_0.100.txt"]
    alp = [0.001, 0.005, 0.01, 0.05, 0.1]

    print("Triangular lattice")
    
    Tstop = []

    data1 = time_stop(n1, alp)
    Tstop.append(data1[2])

    fig1 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_tri_1.pdf", dpi=150, bbox_inches="tight")
    #plt.show()

    data1 = time_stop(n2, alp)
    Tstop.append(data1[2])
    fig2 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_tri_2.pdf", dpi=150, bbox_inches="tight")
    #plt.show()

    data1 = time_stop(n3, alp)
    Tstop.append(data1[2])
    fig3 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_tri_3.pdf", dpi=150, bbox_inches="tight")
    #plt.show()


    par1 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[0]))[0]
    par2 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[1]))[0]
    par3 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[2]))[0]
    
    print(par1[0], par2[0], par3[0])

    t = np.linspace(min(alp), max(alp))

    fig4 = plt.figure(figsize=(10, 5))
    plt.plot(t, p(10,fit(np.log10(t), *par1)), "r--", label=r"$\sim \alpha^{%.2f}$" % par1[0])
    plt.plot(t, p(10,fit(np.log10(t), *par2)), "b--", label=r"$\sim \alpha^{%.2f}$" % par2[0])
    plt.plot(t, p(10,fit(np.log10(t), *par3)), "g--", label=r"$\sim \alpha^{%.2f}$" % par3[0])
    plt.plot(alp, p(10, np.array(Tstop[0])), "ro", label=r"$n_\mathrm{max}=1$")
    plt.plot(alp, p(10, np.array(Tstop[1])), "bo", label=r"$n_\mathrm{max}=2$")
    plt.plot(alp, p(10, np.array(Tstop[2])), "go", label=r"$n_\mathrm{max}=3$")
    #plt.legend()
    plt.ylabel(r"${T_\mathrm{stop}^\alpha}$")
    plt.xlabel(r"${\alpha}$")
    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("stop_tri_fit.pdf", dpi=150, bbox_inches="tight")

    alp = np.asarray(alp)

    par1 = opt.curve_fit(fit, np.log10(1/alp), np.array(Tstop[0]))[0]
    par2 = opt.curve_fit(fit, np.log10(1/alp), np.array(Tstop[1]))[0]
    par3 = opt.curve_fit(fit, np.log10(1/alp), np.array(Tstop[2]))[0]
    
    print(par1[0], par2[0], par3[0])

    t = np.linspace(min(1/alp), max(1/alp))
    fig5 = plt.figure(figsize=(10, 5))
    plt.plot(t, p(10,fit(np.log10(t), *par1)), "r--", label=r"$\sim \alpha^{%.2f}$" % par1[0])
    plt.plot(t, p(10,fit(np.log10(t), *par2)), "b--", label=r"$\sim \alpha^{%.2f}$" % par2[0])
    plt.plot(t, p(10,fit(np.log10(t), *par3)), "g--", label=r"$\sim \alpha^{%.2f}$" % par3[0])
    plt.plot(1/alp, p(10, np.array(Tstop[0])), "ro", label=r"$n_\mathrm{max}=1$")
    plt.plot(1/alp, p(10, np.array(Tstop[1])), "bo", label=r"$n_\mathrm{max}=2$")
    plt.plot(1/alp, p(10, np.array(Tstop[2])), "go", label=r"$n_\mathrm{max}=3$")
    #plt.legend()
    plt.ylabel(r"${T_\mathrm{stop}^\alpha}$")
    plt.xlabel(r"${\alpha^{-1}}$")
    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("stop_tri_fit_inv.pdf", dpi=150, bbox_inches="tight")
    

if False:
    dat = np.genfromtxt("square_3_0.600.txt")
    k = dat.flatten()
    y, x = np.histogram(k, bins=40, density=True)
    s = 0
    for i in range(len(y)):
        s+= y[i] * (x[i+1] + x[i])/2

    s = s/np.sum(y)
    print(s)
    l = plt.figure()
    plt.hist(dat, bins=30, density=True)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$T(t)$")
    plt.savefig("hist.pdf", dpi=200)

    dat = np.genfromtxt("square_3_0.001.txt")
    k = dat.flatten()
    y, x = np.histogram(k, bins=40, density=True)
    s = 0
    for i in range(len(y)):
        s+= y[i] * (x[i+1] + x[i])/2

    s = s/np.sum(y)
    print(s)
    f = plt.figure()
    plt.hist(dat, bins=30, density=True)
    plt.xlabel(r"$t$")
    plt.ylabel(r"$T(t)$")
    plt.savefig("hist_low.pdf", dpi=200)

# square
if lat == "s":
    plt.rcParams.update({'font.size': 20})
    n1 = ["square_1_0.001.txt", "square_1_0.002.txt", "square_1_0.003.txt", "square_1_0.004.txt", "square_1_0.005.txt", "square_1_0.006.txt", "square_1_0.007.txt", "square_1_0.008.txt", "square_1_0.009.txt", "square_1_0.010.txt"]
    n2 = ["square_2_0.001.txt", "square_2_0.002.txt", "square_2_0.003.txt", "square_2_0.004.txt", "square_2_0.005.txt", "square_2_0.006.txt", "square_2_0.007.txt", "square_2_0.008.txt", "square_2_0.009.txt", "square_2_0.010.txt"]
    n3 = ["square_3_0.001.txt", "square_3_0.002.txt", "square_3_0.003.txt", "square_3_0.004.txt", "square_3_0.005.txt", "square_3_0.006.txt", "square_3_0.007.txt", "square_3_0.008.txt", "square_3_0.009.txt", "square_3_0.010.txt"]
    alp = [0.001 * i for i in range(1, 11)]

    print("Square lattice")
    
    Tstop = []

    data1 = time_stop(n1, alp)
    Tstop.append(data1[2])
    data1 = time_stop(n2, alp)
    Tstop.append(data1[2])
    data1 = time_stop(n3, alp)
    Tstop.append(data1[2])
    """
    fig1 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_squ_1.pdf", dpi=150, bbox_inches="tight")
    #plt.show()

    data1 = time_stop(n2, alp)
    Tstop.append(data1[2])
    fig2 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_squ_2.pdf", dpi=150, bbox_inches="tight")
    #plt.show()

    data1 = time_stop(n3, alp)
    Tstop.append(data1[2])
    fig3 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_squ_3.pdf", dpi=150, bbox_inches="tight")
    #plt.show()
    """

    """par1 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[0]))[0]
    par2 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[1]))[0]"""
    
    

    t = np.linspace(min(alp), max(alp))

    fig4 = plt.figure(figsize=(10, 5))
    #plt.plot(t, p(10,fit(np.log10(t), *par1)), "r--", label=r"$\sim \alpha^{%.2f}$" % par1[0])
    #plt.plot(t, p(10,fit(np.log10(t), *par2)), "b--", label=r"$\sim \alpha^{%.2f}$" % par2[0])
    #plt.plot(t, p(10,fit(np.log10(t), *par3)), "g--", label=r"$\sim \alpha^{%.2f}$" % par3[0])
    plt.plot(alp, np.array(Tstop[0]), "ro", label=r"$n_\mathrm{max}=1$")
    plt.plot(alp, np.array(Tstop[1]), "bo", label=r"$n_\mathrm{max}=2$")
    plt.plot(alp, np.array(Tstop[2]), "go", label=r"$n_\mathrm{max}=3$")
    #plt.legend()
    plt.ylabel(r"${T_\mathrm{stop}^\alpha}$")
    plt.xlabel(r"${\alpha}$")
    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("stop_squ_fit.pdf", dpi=150, bbox_inches="tight")
    


    """
    n1 = ["square_3_0.001.txt", "square_3_0.002.txt", "square_3_0.003.txt", "square_3_0.004.txt", "square_3_0.005.txt", "square_3_0.006.txt", "square_3_0.007.txt", "square_3_0.008.txt", "square_3_0.009.txt", "square_3_0.010.txt", "square_3_0.020.txt", "square_3_0.030.txt", "square_3_0.040.txt", "square_3_0.050.txt", "square_3_0.060.txt", "square_3_0.070.txt", "square_3_0.080.txt", "square_3_0.090.txt"]#, "square_3_0.100.txt"]
    n2 = ["square_3_0.100.txt", "square_3_0.200.txt", "square_3_0.300.txt", "square_3_0.400.txt", "square_3_0.500.txt", "square_3_0.600.txt"]
    alp2 = [0.100, 0.200, 0.100*3, 0.100*4, 0.100*5, 0.100*6]
    alp = [0.001 * i for i in range(1, 11)]
    for i in range(2, 10):
        alp.append(0.01*i)

    print("Square lattice")
    
    Tstop = []

    data1 = time_stop(n1, alp)
    Tstop.append(data1[2])

    
    Tstop1 = []

    data1 = time_stop(n2, alp2)
    Tstop1.append(data1[2])
    t = np.linspace(min(alp[:7]), max(alp[:7]))

    print(p(10, np.array(Tstop1[0])), np.array(Tstop1[0]))
    par3 = opt.curve_fit(fit, np.log10(alp[:7]), np.log10(np.array(Tstop[0][:7])))[0]
    par4 = opt.curve_fit(fit, np.log10(alp[11:]), np.array(Tstop[0][11:]))[0]
    t1 = np.linspace(min(alp[11:]), max(alp[11:]))

    fig5 = plt.figure(figsize=(10, 5))
    plt.plot(t, p(10,fit(np.log10(t), *par3)), "r--", label=r"$\sim \alpha^{%.2f}$" % par3[0])
    #plt.plot(t1, p(10,fit(np.log10(t1), *par4)), "b--", label=r"$\sim \alpha^{%.2f}$" % par4[0])
    plt.plot(alp, np.array(Tstop[0]), "go", label=r"$n_\mathrm{max}=3$")
    plt.plot(alp2, np.array(Tstop1[0]), "go")
    plt.plot([min(alp), max(alp2)], [1, 1], "-.", label="Random Walk Limit")
    plt.legend()
    plt.ylabel(r"${T_\mathrm{stop}^\alpha}$")
    plt.xlabel(r"${\alpha}$")
    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("stop_squ_fit_3.pdf", dpi=150, bbox_inches="tight")"""


# hexagonal
if lat == "h":
    n1 = ["hex_1_0.001.txt", "hex_1_0.005.txt", "hex_1_0.010.txt", "hex_1_0.050.txt", "hex_1_0.100.txt"]
    n2 = ["hex_2_0.001.txt", "hex_2_0.005.txt", "hex_2_0.010.txt", "hex_2_0.050.txt", "hex_2_0.100.txt"]
    n3 = ["hex_3_0.001.txt", "hex_3_0.005.txt", "hex_3_0.010.txt", "hex_3_0.050.txt", "hex_3_0.100.txt"]
    alp = [0.001, 0.005, 0.01, 0.05, 0.1]

    print("Hexagonal lattice")
    
    Tstop = []

    data1 = time_stop(n1, alp)
    Tstop.append(data1[2])

    fig1 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_hex_1.pdf", dpi=150, bbox_inches="tight")
    #plt.show()

    data1 = time_stop(n2, alp)
    Tstop.append(data1[2])
    fig2 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_hex_2.pdf", dpi=150, bbox_inches="tight")
    #plt.show()

    data1 = time_stop(n3, alp)
    Tstop.append(data1[2])
    fig3 = plt.figure()
    plt.plot(10**data1[2][0]*np.ones(2), np.array([0, max(data1[0][0][ :])]), "r--", label=r"$T_\mathrm{stop}^{0.001}$")
    plt.plot(10**data1[2][2]*np.ones(2), np.array([0, max(data1[0][2][ :])]), "b--", label=r"$T_\mathrm{stop}^{0.01}$")
    plt.plot(10**data1[2][4]*np.ones(2), np.array([0, max(data1[0][4][ :])]), "g--", label=r"$T_\mathrm{stop}^{0.1}$")
    plt.plot(np.power(10,data1[1][0][:]), data1[0][0][:], "r-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,data1[1][2][:]), data1[0][2][:], "b-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,data1[1][4][:]), data1[0][4][:], "g-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.savefig("stop_hex_3.pdf", dpi=150, bbox_inches="tight")
    #plt.show()


    par1 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[0]))[0]
    par2 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[1]))[0]
    par3 = opt.curve_fit(fit, np.log10(alp), np.array(Tstop[2]))[0]
    
    print(par1[0], par2[0], par3[0])

    t = np.linspace(min(alp), max(alp))

    fig4 = plt.figure(figsize=(10, 5))
    plt.plot(t, p(10,fit(np.log10(t), *par1)), "r--", label=r"$\sim \alpha^{%.2f}$" % par1[0])
    plt.plot(t, p(10,fit(np.log10(t), *par2)), "b--", label=r"$\sim \alpha^{%.2f}$" % par2[0])
    plt.plot(t, p(10,fit(np.log10(t), *par3)), "g--", label=r"$\sim \alpha^{%.2f}$" % par3[0])
    plt.plot(alp, p(10, np.array(Tstop[0])), "ro", label=r"$n_\mathrm{max}=1$")
    plt.plot(alp, p(10, np.array(Tstop[1])), "bo", label=r"$n_\mathrm{max}=2$")
    plt.plot(alp, p(10, np.array(Tstop[2])), "go", label=r"$n_\mathrm{max}=3$")
    #plt.legend()
    plt.ylabel(r"${T_\mathrm{stop}^\alpha}$")
    plt.xlabel(r"${\alpha}$")
    plt.grid()
    plt.yscale("log")
    plt.xscale("log")
    plt.savefig("stop_hex_fit.pdf", dpi=150, bbox_inches="tight")