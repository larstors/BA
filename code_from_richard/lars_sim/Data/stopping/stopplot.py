import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import numba as nb
from matplotlib import rcParams
plt.rcParams.update({'font.size': 20})


print("Lattice? (only first letter)")
lat = input()


bin= 1000


# triangle
if lat == "t":
    # n=1
    d1_001 = np.genfromtxt("tri_1_0.001.txt", delimiter=" ")
    d1_01 = np.genfromtxt("tri_1_0.010.txt", delimiter=" ")
    d1_1 = np.genfromtxt("tri_1_0.100.txt", delimiter=" ")

    # n=2
    d2_001 = np.genfromtxt("tri_2_0.001.txt", delimiter=" ")
    d2_01 = np.genfromtxt("tri_2_0.010.txt", delimiter=" ")
    d2_1 = np.genfromtxt("tri_2_0.100.txt", delimiter=" ")

    # n=3
    d3_001 = np.genfromtxt("tri_3_0.001.txt", delimiter=" ")
    d3_01 = np.genfromtxt("tri_3_0.010.txt", delimiter=" ")
    d3_1 = np.genfromtxt("tri_3_0.100.txt", delimiter=" ")


    d1_001 = np.log10(d1_001).flatten()
    d1_001, h1_001 = np.histogram(d1_001.flatten(), bins=bin, density=True)
    x1_001 = []
    d1_01 = np.log10(d1_01).flatten()
    d1_01, h1_01 = np.histogram(d1_01.flatten(), bins=bin, density=True)
    x1_01 = []
    d1_1 = np.log10(d1_1).flatten()
    d1_1, h1_1 = np.histogram(d1_1.flatten(), bins=bin, density=True)
    x1_1 = []
    for i in range(len(d1_001)):
        x1_001.append((h1_001[i+1]+h1_001[i])/2)
        x1_01.append((h1_01[i+1]+h1_01[i])/2)
        x1_1.append((h1_1[i+1]+h1_1[i])/2)


    
    d2_001 = np.log10(d2_001).flatten()
    d2_001, h2_001 = np.histogram(d2_001.flatten(), bins=bin, density=True)
    x2_001 = []
    d2_01 = np.log10(d2_01).flatten()
    d2_01, h2_01 = np.histogram(d2_01.flatten(), bins=bin, density=True)
    x2_01 = []
    d2_1 = np.log10(d2_1).flatten()
    d2_1, h2_1 = np.histogram(d2_1.flatten(), bins=bin, density=True)
    x2_1 = []
    for i in range(len(d2_001)):
        x2_001.append((h2_001[i+1]+h2_001[i])/2)
        x2_01.append((h2_01[i+1]+h2_01[i])/2)
        x2_1.append((h2_1[i+1]+h2_1[i])/2)



    d3_001 = np.log10(d3_001).flatten()
    d3_001, h3_001 = np.histogram(d3_001.flatten(), bins=bin, density=True)
    x3_001 = []
    d3_01 = np.log10(d3_01).flatten()
    d3_01, h3_01 = np.histogram(d3_01.flatten(), bins=bin, density=True)
    x3_01 = []
    d3_1 = np.log10(d3_1).flatten()
    d3_1, h3_1 = np.histogram(d3_1.flatten(), bins=bin, density=True)
    x3_1 = []
    for i in range(len(d3_001)):
        x3_001.append((h3_001[i+1]+h3_001[i])/2)
        x3_01.append((h3_01[i+1]+h3_01[i])/2)
        x3_1.append((h3_1[i+1]+h3_1[i])/2)


    fig1 = plt.figure()
    plt.plot(np.power(10,x1_001), d1_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x1_01), d1_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x1_1), d1_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()

    fig2 = plt.figure()
    plt.plot(np.power(10,x2_001), d2_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x2_01), d2_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x2_1), d2_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()

    fig3 = plt.figure()
    plt.plot(np.power(10,x3_001), d3_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x3_01), d3_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x3_1), d3_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()



# square
if lat == "s":
    # n=1
    d1_001 = np.genfromtxt("square_1_0.001.txt", delimiter=" ")
    d1_01 = np.genfromtxt("square_1_0.010.txt", delimiter=" ")
    d1_1 = np.genfromtxt("square_1_0.100.txt", delimiter=" ")

    # n=2
    d2_001 = np.genfromtxt("square_2_0.001.txt", delimiter=" ")
    d2_01 = np.genfromtxt("square_2_0.010.txt", delimiter=" ")
    d2_1 = np.genfromtxt("square_2_0.100.txt", delimiter=" ")

    # n=3
    d3_001 = np.genfromtxt("square_3_0.001.txt", delimiter=" ")
    d3_01 = np.genfromtxt("square_3_0.010.txt", delimiter=" ")
    d3_1 = np.genfromtxt("square_3_0.100.txt", delimiter=" ")


    d1_001 = np.log10(d1_001).flatten()
    d1_001, h1_001 = np.histogram(d1_001.flatten(), bins=bin, density=True)
    x1_001 = []
    d1_01 = np.log10(d1_01).flatten()
    d1_01, h1_01 = np.histogram(d1_01.flatten(), bins=bin, density=True)
    x1_01 = []
    d1_1 = np.log10(d1_1).flatten()
    d1_1, h1_1 = np.histogram(d1_1.flatten(), bins=bin, density=True)
    x1_1 = []
    for i in range(len(d1_001)):
        x1_001.append((h1_001[i+1]+h1_001[i])/2)
        x1_01.append((h1_01[i+1]+h1_01[i])/2)
        x1_1.append((h1_1[i+1]+h1_1[i])/2)


    
    d2_001 = np.log10(d2_001).flatten()
    d2_001, h2_001 = np.histogram(d2_001.flatten(), bins=bin, density=True)
    x2_001 = []
    d2_01 = np.log10(d2_01).flatten()
    d2_01, h2_01 = np.histogram(d2_01.flatten(), bins=bin, density=True)
    x2_01 = []
    d2_1 = np.log10(d2_1).flatten()
    d2_1, h2_1 = np.histogram(d2_1.flatten(), bins=bin, density=True)
    x2_1 = []
    for i in range(len(d2_001)):
        x2_001.append((h2_001[i+1]+h2_001[i])/2)
        x2_01.append((h2_01[i+1]+h2_01[i])/2)
        x2_1.append((h2_1[i+1]+h2_1[i])/2)



    d3_001 = np.log10(d3_001).flatten()
    d3_001, h3_001 = np.histogram(d3_001.flatten(), bins=bin, density=True)
    x3_001 = []
    d3_01 = np.log10(d3_01).flatten()
    d3_01, h3_01 = np.histogram(d3_01.flatten(), bins=bin, density=True)
    x3_01 = []
    d3_1 = np.log10(d3_1).flatten()
    d3_1, h3_1 = np.histogram(d3_1.flatten(), bins=bin, density=True)
    x3_1 = []
    for i in range(len(d3_001)):
        x3_001.append((h3_001[i+1]+h3_001[i])/2)
        x3_01.append((h3_01[i+1]+h3_01[i])/2)
        x3_1.append((h3_1[i+1]+h3_1[i])/2)


    fig1 = plt.figure()
    plt.plot(np.power(10,x1_001), d1_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x1_01), d1_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x1_1), d1_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()

    fig2 = plt.figure()
    plt.plot(np.power(10,x2_001), d2_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x2_01), d2_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x2_1), d2_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()

    fig3 = plt.figure()
    plt.plot(np.power(10,x3_001), d3_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x3_01), d3_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x3_1), d3_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()


# hexagonal
if lat == "h":
    # n=1
    d1_001 = np.genfromtxt("hex_1_0.001.txt", delimiter=" ")
    d1_01 = np.genfromtxt("hex_1_0.010.txt", delimiter=" ")
    d1_1 = np.genfromtxt("hex_1_0.100.txt", delimiter=" ")

    # n=2
    d2_001 = np.genfromtxt("hex_2_0.001.txt", delimiter=" ")
    d2_01 = np.genfromtxt("hex_2_0.010.txt", delimiter=" ")
    d2_1 = np.genfromtxt("hex_2_0.100.txt", delimiter=" ")

    # n=3
    d3_001 = np.genfromtxt("hex_3_0.001.txt", delimiter=" ")
    d3_01 = np.genfromtxt("hex_3_0.010.txt", delimiter=" ")
    d3_1 = np.genfromtxt("hex_3_0.100.txt", delimiter=" ")


    d1_001 = np.log10(d1_001).flatten()
    d1_001, h1_001 = np.histogram(d1_001.flatten(), bins=bin, density=True)
    x1_001 = []
    d1_01 = np.log10(d1_01).flatten()
    d1_01, h1_01 = np.histogram(d1_01.flatten(), bins=bin, density=True)
    x1_01 = []
    d1_1 = np.log10(d1_1).flatten()
    d1_1, h1_1 = np.histogram(d1_1.flatten(), bins=bin, density=True)
    x1_1 = []
    for i in range(len(d1_001)):
        x1_001.append((h1_001[i+1]+h1_001[i])/2)
        x1_01.append((h1_01[i+1]+h1_01[i])/2)
        x1_1.append((h1_1[i+1]+h1_1[i])/2)


    
    d2_001 = np.log10(d2_001).flatten()
    d2_001, h2_001 = np.histogram(d2_001.flatten(), bins=bin, density=True)
    x2_001 = []
    d2_01 = np.log10(d2_01).flatten()
    d2_01, h2_01 = np.histogram(d2_01.flatten(), bins=bin, density=True)
    x2_01 = []
    d2_1 = np.log10(d2_1).flatten()
    d2_1, h2_1 = np.histogram(d2_1.flatten(), bins=bin, density=True)
    x2_1 = []
    for i in range(len(d2_001)):
        x2_001.append((h2_001[i+1]+h2_001[i])/2)
        x2_01.append((h2_01[i+1]+h2_01[i])/2)
        x2_1.append((h2_1[i+1]+h2_1[i])/2)



    d3_001 = np.log10(d3_001).flatten()
    d3_001, h3_001 = np.histogram(d3_001.flatten(), bins=bin, density=True)
    x3_001 = []
    d3_01 = np.log10(d3_01).flatten()
    d3_01, h3_01 = np.histogram(d3_01.flatten(), bins=bin, density=True)
    x3_01 = []
    d3_1 = np.log10(d3_1).flatten()
    d3_1, h3_1 = np.histogram(d3_1.flatten(), bins=bin, density=True)
    x3_1 = []
    for i in range(len(d3_001)):
        x3_001.append((h3_001[i+1]+h3_001[i])/2)
        x3_01.append((h3_01[i+1]+h3_01[i])/2)
        x3_1.append((h3_1[i+1]+h3_1[i])/2)


    fig1 = plt.figure()
    plt.plot(np.power(10,x1_001), d1_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x1_01), d1_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x1_1), d1_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()

    fig2 = plt.figure()
    plt.plot(np.power(10,x2_001), d2_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x2_01), d2_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x2_1), d2_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()

    fig3 = plt.figure()
    plt.plot(np.power(10,x3_001), d3_001, "-", label=r"$\alpha=0.001$")
    plt.plot(np.power(10,x3_01), d3_01, "-", label=r"$\alpha=0.01$")
    plt.plot(np.power(10,x3_1), d3_1, "-", label=r"$\alpha=0.1$")
    plt.xscale("log")
    plt.legend()
    plt.show()