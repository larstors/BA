import numpy as np
import matplotlib.pyplot as plt
import scipy as sc



if (False):
    """
    We will start with triangular lattice to figure out how well this works/looks
    """

    shiftx, shifty = 0, 0

    tr001 = np.genfromtxt("tri_number_alpha0.001_phi1.20_L100_3.txt", delimiter=" ")
    tr010 = np.genfromtxt("tri_number_alpha0.010_phi1.20_L100_3.txt", delimiter=" ")

    fig, ax = plt.subplots()

    ax.plot(tr001[:, 2]/12000-shiftx, tr001[:, 1]-shifty, label=r"$\alpha=0.001$")
    ax.plot(tr010[:, 2]/12000-shiftx, tr010[:, 1]-shifty, label=r"$\alpha=0.010$")
    ax.plot(tr010[0, 2]/12000-shiftx, tr010[0, 1]-shifty, "ro", label=r"$\alpha=0.010$ Start")
    ax.plot(tr001[0, 2]/12000-shiftx, tr001[0, 1]-shifty, "go", label=r"$\alpha=0.001$ Start")
    
    ax.annotate("S", xy=[0.8, 0.4], xytext=[0.8, 0.4])
    ax.annotate("C", xy=[0.1, 0.6], xytext=[0.1, 0.6])
    ax.annotate("G", xy=[0.1, 0.05], xytext=[0.1, 0.05])
    # ax.spines['left'].set_position('center')
    # ax.spines['bottom'].set_position('center')
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    # ax.xaxis.set_ticks_position('bottom')
    # ax.yaxis.set_ticks_position('left')

    ax.axhline(color='black', lw=0.5, y=.2)
    ax.axvline(color='black', lw=0.5, x=.2)
    
    ax.set_xlabel(r"$c_N$")
    ax.set_ylabel(r"$w_N$")
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    #ax.grid()
    ax.legend()
    ax.set_title("Triangle")
    plt.savefig("Tri_trajectory_3_rho04.pdf", dpi=200)

if (True):
    """
    We will start with square lattice to figure out how well this works/looks
    """

    shiftx, shifty = 0.3, 0.2

    tr001 = np.genfromtxt("square_number_alpha0.001_phi1.20_L100_3.txt", delimiter=" ")
    tr010 = np.genfromtxt("square_number_alpha0.010_phi1.20_L100_3.txt", delimiter=" ")

    fig, ax = plt.subplots()

    ax.plot(tr001[:, 2]-shiftx, tr001[:, 1]-shifty, label=r"$\alpha=0.001$")
    ax.plot(tr010[:, 2]/12000-shiftx, tr010[:, 1]-shifty, label=r"$\alpha=0.010$")
    ax.plot(tr010[0, 2]/12000-shiftx, tr010[0, 1]-shifty, "ro", label=r"$\alpha=0.010$ Start")
    ax.plot(tr001[0, 2]-shiftx, tr001[0, 1]-shifty, "go", label=r"$\alpha=0.001$ Start")
    
    # ax.spines['left'].set_position('center')
    # ax.spines['bottom'].set_position('center')
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    # ax.xaxis.set_ticks_position('bottom')
    # ax.yaxis.set_ticks_position('left')

    ax.axhline(color='black', lw=0.5)
    ax.axvline(color='black', lw=0.5)
    
    ax.set_xlabel(r"$c_N$")
    ax.set_ylabel(r"$w_N$")
    #ax.set_xscale("log")
    #ax.set_yscale("log")
    #ax.grid()
    ax.legend()
    ax.set_title("Square")
    plt.savefig("squ_tra.pdf")

if (False):
    """
    We will start with hexagonal lattice to figure out how well this works/looks
    """

    shiftx, shifty = 0, 0

    tr001 = np.genfromtxt("hex_number_alpha0.001_phi2.40_L100_3.txt", delimiter=" ")
    tr010 = np.genfromtxt("hex_number_alpha0.010_phi2.40_L100_3.txt", delimiter=" ")

    print(min(tr001[:, 2]/24000-shiftx), min(tr010[:, 2]/24000-shiftx))

    fig, ax = plt.subplots()

    ax.plot(tr001[:, 2]/24000-shiftx, tr001[:, 1]-shifty, label=r"$\alpha=0.001$")
    ax.plot(tr010[:, 2]/24000-shiftx, tr010[:, 1]-shifty, label=r"$\alpha=0.010$")
    ax.plot(tr010[0, 2]/24000-shiftx, tr010[0, 1]-shifty, "ro", label=r"$\alpha=0.010$ Start")
    ax.plot(tr001[0, 2]/24000-shiftx, tr001[0, 1]-shifty, "go", label=r"$\alpha=0.001$ Start")
    
    # ax.spines['left'].set_position('center')
    # ax.spines['bottom'].set_position('center')
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    # ax.xaxis.set_ticks_position('bottom')
    # ax.yaxis.set_ticks_position('left')

    #ax.axhline(color='black', lw=0.5)
    #ax.axvline(color='black', lw=0.5)
    
    ax.set_xlabel(r"$c_N$")
    ax.set_ylabel(r"$w_N$")
    ax.set_xscale("log")
    ax.set_yscale("log")
    #ax.grid()
    ax.legend()
    ax.set_title("Hexagonal")
    plt.show()

    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
    plt.tight_layout()

    ax[0].plot(tr001[:, 0], tr001[:, 2]/24000, label=r"$\alpha=0.001, c_N$")
    ax[0].plot(tr001[:, 0], tr001[:, 1], label=r"$\alpha=0.001, w_N$")
    ax[0].plot(tr010[:, 0], tr010[:, 2]/24000, label=r"$\alpha=0.010, c_N$")
    ax[0].plot(tr010[:, 0], tr010[:, 1], label=r"$\alpha=0.010, w_N$")
    ax[0].grid()
    ax[0].set_yscale("log")
    ax[0].legend()
    ax[0].set_ylabel(r"$w_N$")
    ax[0].set_xscale("log")

    ax[1].plot(tr001[:, 0], tr001[:, -1], label=r"$\alpha=0.001, c_N$")
    ax[1].plot(tr010[:, 0], tr010[:, -1], label=r"$\alpha=0.010, c_N$")
    ax[1].grid()
    ax[1].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].legend()
    ax[1].set_ylabel("Mean cluster size")
    ax[1].set_xlabel(r"Time $(t)$")
    plt.show()

if (False):
    tr001 = np.genfromtxt("tri_number_alpha0.001_phi1.20_L100_3.txt", delimiter=" ")
    tr010 = np.genfromtxt("tri_number_alpha0.010_phi1.20_L100_3.txt", delimiter=" ")
    fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
    plt.tight_layout()

    ax[0].plot(tr001[1:-2, 0], tr001[1:-2, 3]*10/11, label=r"$\alpha=0.001, J$")
    ax[0].plot(tr010[1:-2, 0], tr010[1:-2, 3]*10/11, label=r"$\alpha=0.010, J$")
    ax[0].plot(tr001[1:-2, 0], tr001[1:-2, 1]*10/11, label=r"$\alpha=0.001, w_N$")
    ax[0].plot(tr010[1:-2, 0], tr010[1:-2, 1]*10/11, label=r"$\alpha=0.010, w_N$")
    ax[0].plot(tr001[1:-2, 0], tr001[1:-2, 4]/12000*10/11, label=r"$\alpha=0.001, M$")
    ax[0].plot(tr010[1:-2, 0], tr010[1:-2, 4]/12000*10/11, label=r"$\alpha=0.010, M$")
    # ax[0].plot(tr001[1:-2, 0], tr001[1:-2, 2]/12000, label=r"$\alpha=0.001, c_N$")
    # ax[0].plot(tr010[1:-2, 0], tr010[1:-2, 2]/12000, label=r"$\alpha=0.010, c_N$")
    ax[0].grid()
    ax[0].set_yscale("log")
    ax[0].legend(loc="right")
    ax[0].set_ylabel(r"$J$")
    ax[0].set_xscale("log")

    ax[1].plot(tr001[1:-2, 0], tr001[1:-2, -1], label=r"$\alpha=0.001$")
    ax[1].plot(tr010[1:-2, 0], tr010[1:-2, -1], label=r"$\alpha=0.010$")
    ax[1].grid()
    ax[1].set_yscale("log")
    ax[1].set_xscale("log")
    ax[1].legend()
    ax[1].set_ylabel("Mean cluster size")
    ax[1].set_xlabel(r"Time $(t)$")
    plt.show()