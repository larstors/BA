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

#plt.rcParams.update({'font.size': 22})

L = 100

what = int(input())




def func(x):
    return x**2*1/(4*np.pi) - 1/(4*np.pi) + 1

def fit_func(x, a):
    return x**a

# square
if what == 1:
    f = open("square_part.txt")
    square_part = []
    for line in f:
        square_part.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("square_sv.txt")
    square_sv = []
    for line in f:
        square_sv.append(np.loadtxt(StringIO(line), dtype=int))
    
    f = open("square_clust.txt")
    square_cl = []
    for line in f:
        square_cl.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("square_border.txt")
    square_int = []
    for line in f:
        square_int.append(np.loadtxt(StringIO(line), dtype=int))

    def conv(n):
        
        x = []
        y = []
        for i in n:
            x.append(i % L)
            y.append(int(i/L))
        
        return x, y

    """
    # finding out which is the smaller cluster (as this is not done automatically by the c++ algorithm)
    if (len(square_cl[0]) > len(square_cl[1])):
        big = square_cl[0]
        small = square_cl[1]
    else:
        big = square_cl[1]
        small = square_cl[0]
    """

    


    s = []
    v = []
    for i in range(len(square_sv)):
        for j in range(1, len(square_sv[i]), 2):
            s.append(square_sv[i][j])
            v.append(square_sv[i][j-1])

    par = opt.curve_fit(fit_func, s, v)[0]
    print(par)
    """
    x_small = conv(small[2:])[0]
    y_small = conv(small[2:])[1]

    x_big = np.asarray(conv(big[2:])[0])
    y_big = conv(big[2:])[1]
    
    k = np.argwhere(x_big[:] < 0)
    for i in k:
        x_big[i] += 100
    """

    
    fig1, ax1 = plt.subplots()
    """
    left, bottom, width, height = [0.15, 0.4, 0.3, 0.3]
    ax2 = fig1.add_axes([left, bottom, width, height])
    ax2.scatter(x=x_small, y=y_small, marker="s", color="black", s=30)
    ax2.axis([min(x_small) - 5, max(x_small) + 5, min(y_small) - 5, max(y_small) + 5])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_color('green')  
    
    left, bottom, width, height = [0.6, 0.2, 0.3, 0.3]
    ax3 = fig1.add_axes([left, bottom, width, height])
    ax3.scatter(x=x_big, y=y_big, marker="s", color="black", s=1e-1)
    ax3.axis([min(x_big) - 1, max(x_big) + 1, min(y_big) - 5, max(y_big) + 5])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_color('blue')  
    """

    #plt.plot(np.linspace(1, 140), func(np.linspace(1, 140)))
    ax1.plot(np.linspace(1, max(s)), fit_func(np.linspace(1, max(s)), par[0]), "k-", label=r"$x^{%g}$" % par[0])
    for i in range(6):
        ax1.plot(square_sv[i*16][1::2], square_sv[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
    #ax1.plot(small[1], small[0], "gs")
    #ax1.plot(big[1], big[0], "bs")
    ax1.set_xlabel(r"Surface $S$")
    ax1.set_ylabel(r"Volume $V$")
    ax1.grid()
    ax1.legend()
    plt.savefig("sq_n_3.pdf", dpi=150)
    #plt.show()



    fig2, ax1 = plt.subplots()
    """
    left, bottom, width, height = [0.15, 0.4, 0.3, 0.3]
    ax2 = fig2.add_axes([left, bottom, width, height])
    ax2.scatter(x=x_small, y=y_small, marker="s", color="black", s=30)
    ax2.axis([min(x_small) - 5, max(x_small) + 5, min(y_small) - 5, max(y_small) + 5])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_color('green')  
    
    left, bottom, width, height = [0.6, 0.2, 0.3, 0.3]
    ax3 = fig2.add_axes([left, bottom, width, height])
    ax3.scatter(x=x_big, y=y_big, marker="s", color="black", s=1)
    ax3.axis([min(x_big) - 1, max(x_big) + 1, min(y_big) - 5, max(y_big) + 5])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_color('blue')  
    """

    #plt.plot(np.linspace(1, 140), func(np.linspace(1, 140)))
    ax1.plot(np.linspace(1, max(s)), fit_func(np.linspace(1, max(s)), par[0]), "k-", label=r"$x^{%g}$" % par[0])
    for i in range(6):
        ax1.plot(square_sv[i*16][1::2], square_sv[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
    #ax1.plot(small[1], small[0], "gs")
    #ax1.plot(big[1], big[0], "bs")
    ax1.set_xlabel(r"Surface $S$")
    ax1.set_ylabel(r"Volume $V$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid()
    ax1.legend()
    plt.savefig("sq_n_3_log.pdf", dpi=150)
    #plt.show()

    fig3, ax = plt.subplots(nrows=3, ncols=2, figsize=(20, 20))
    plt.tight_layout()
    for i in range(6):
        ax[i//2,i%2].set_title(r"$t=%d\cdot 10^3$" % (16*i))
        ax[i//2,i%2].scatter(x=conv(square_part[16*i][::2])[0], y=conv(square_part[16*i][::2])[1], marker="s", color="black", s=4)
        ax[i//2,i%2].set_xticks([])
        ax[i//2,i%2].set_yticks([])
        ax[i//2,i%2].scatter(x=conv(square_int[16*i][:])[0], y=conv(square_int[16*i][:])[1], marker="s", color="red", s=4)

    #plt.show()
    plt.savefig("snapshot_sq_n_3.pdf", dpi=200)


#tri
if what == 2:
    f = open("tri_part.txt")
    square_part = []
    for line in f:
        square_part.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("tri_sv.txt")
    square_sv = []
    for line in f:
        square_sv.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("tri_clust.txt")
    square_cl = []
    for line in f:
        square_cl.append(np.loadtxt(StringIO(line), dtype=int))
    
    f = open("tri_border.txt")
    square_int = []
    for line in f:
        square_int.append(np.loadtxt(StringIO(line), dtype=int))

    # finding out which is the smaller cluster (as this is not done automatically by the c++ algorithm)
    if (len(square_cl[0]) > len(square_cl[1])):
        big = square_cl[0]
        small = square_cl[1]
    else:
        big = square_cl[1]
        small = square_cl[0]
    

    def tri_conv(n):
    
        x = []
        y = []
        for i in n:
            x.append(i % L - 0.5 * int(i / L))
            y.append(int(i/L) * np.sqrt(3)/2)
        
        return x, y
    
    def trian():
        # indices for entire lattice
        ind = np.array([i for i in range(L*L)])
        # x and y coordinate for each lattice point
        x = tri_conv(ind)[0]
        y = tri_conv(ind)[1]
        
        #for i in range(len(x)):
        #    x[i] += 10

        triangles = []

        for j in range(L-1):
            for i in range(L-1):
                i = i + j * L
                triangles.append([i, i+1, i+L+1])
                triangles.append([i+L, i, i+1+L])
        triang = mtri.Triangulation(x, y, triangles)
        return triang
    
    triang = trian()


    s = []
    v = []
    for i in range(len(square_sv)):
        for j in range(1, len(square_sv[i]), 2):
            s.append(square_sv[i][j])
            v.append(square_sv[i][j-1])

    par = opt.curve_fit(fit_func, s, v)[0]
    print(par)
    
    x_small = tri_conv(small[2:])[0]
    y_small = tri_conv(small[2:])[1]

    x_big = np.asarray(tri_conv(big[2:])[0])
    y_big = tri_conv(big[2:])[1]
    """
    k = np.argwhere(x_big[:] < 0)
    for i in k:
        x_big[i] += 100
    """
    
    fig1, ax1 = plt.subplots()
    """
    left, bottom, width, height = [0.15, 0.4, 0.3, 0.3]
    ax2 = fig1.add_axes([left, bottom, width, height])
    ax2.triplot(triang, 'r-', alpha=0.3, linewidth=0.1)
    ax2.scatter(x=x_small, y=y_small, marker="s", color="black", s=30)
    ax2.axis([min(x_small) - 5, max(x_small) + 5, min(y_small) - 5, max(y_small) + 5])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_color('green')  
    
    left, bottom, width, height = [0.6, 0.2, 0.3, 0.3]
    ax3 = fig1.add_axes([left, bottom, width, height])
    ax3.triplot(triang, 'r-', alpha=0.3, linewidth=0.1)
    ax3.scatter(x=x_big, y=y_big, marker="s", color="black", s=1e-1)
    ax3.axis([min(x_big) - 1, max(x_big) + 1, min(y_big) - 5, max(y_big) + 5])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_color('blue')  
    """

    #plt.plot(np.linspace(1, 140), func(np.linspace(1, 140)))
    ax1.plot(np.linspace(1, max(s)), fit_func(np.linspace(1, max(s)), par[0]), "k-", label=r"$x^{%g}$" % par[0])
    for i in range(6):
        ax1.plot(square_sv[i*16][1::2], square_sv[i*16][::2], "o", label=r"$t=%d\cdot 10^4$" % (16*i))
    #ax1.plot(small[1], small[0], "gs")
    #ax1.plot(big[1], big[0], "bs")
    ax1.set_xlabel(r"Surface $S$")
    ax1.set_ylabel(r"Volume $V$")
    ax1.grid()
    ax1.legend()
    plt.savefig("tr_n_1.pdf", dpi=150)
    #plt.show()



    fig2, ax1 = plt.subplots()
    """
    left, bottom, width, height = [0.15, 0.4, 0.3, 0.3]
    ax2 = fig2.add_axes([left, bottom, width, height])
    ax2.triplot(triang, 'r-', alpha=0.3, linewidth=0.1)
    ax2.scatter(x=x_small, y=y_small, marker="s", color="black", s=30)
    ax2.axis([min(x_small) - 5, max(x_small) + 5, min(y_small) - 5, max(y_small) + 5])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_color('green')  
    
    left, bottom, width, height = [0.6, 0.2, 0.3, 0.3]
    ax3 = fig2.add_axes([left, bottom, width, height])
    ax3.triplot(triang, 'r-', alpha=0.3, linewidth=0.1)
    ax3.scatter(x=x_big, y=y_big, marker="s", color="black", s=1)
    ax3.axis([min(x_big) - 1, max(x_big) + 1, min(y_big) - 5, max(y_big) + 5])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_color('blue')  
    """

    #plt.plot(np.linspace(1, 140), func(np.linspace(1, 140)))
    ax1.plot(np.linspace(1, max(s)), fit_func(np.linspace(1, max(s)), par[0]), "k-", label=r"$x^{%g}$" % par[0])
    for i in range(6):
        ax1.plot(square_sv[i*16][1::2], square_sv[i*16][::2], "o", label=r"$t=%d\cdot 10^4$" % (16*i))
    #ax1.plot(small[1], small[0], "gs")
    #ax1.plot(big[1], big[0], "bs")
    ax1.set_xlabel(r"Surface $S$")
    ax1.set_ylabel(r"Volume $V$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid()
    ax1.legend()
    plt.savefig("tr_n_1_log.pdf", dpi=150)
    #plt.show()



    fig3, ax = plt.subplots(nrows=3, ncols=2, figsize=(20, 20))
    plt.tight_layout()
    for i in range(6):
        ax[i//2,i%2].set_title(r"$t=%d\cdot 10^4$" % (16*i))
        ax[i//2,i%2].scatter(x=tri_conv(square_part[16*i][1::3])[0], y=tri_conv(square_part[16*i][1::3])[1], marker="s", color="black", s=4)
        ax[i//2,i%2].set_xticks([])
        ax[i//2,i%2].set_yticks([])
        ax[i//2,i%2].scatter(x=tri_conv(square_int[16*i][:])[0], y=tri_conv(square_int[16*i][:])[1], marker="s", color="red", s=4)

    plt.savefig("snapshot_tri_n_1.pdf", dpi=200)
    



#hex
if what == 3:
    f = open("hex_part.txt")
    square_part = []
    for line in f:
        square_part.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("hex_sv.txt")
    square_sv = []
    for line in f:
        square_sv.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("hex_clust.txt")
    square_cl = []
    for line in f:
        square_cl.append(np.loadtxt(StringIO(line), dtype=int))

    f = open("hex_border.txt")
    square_int = []
    for line in f:
        square_int.append(np.loadtxt(StringIO(line), dtype=int))


    def conv(n, j):
        """function that converts the index n and site index j into coordinates x and y, with x being horisontal.

        Args:
            n (int): index of particle in lattice

        Returns:
            tuple: x and y coordinates
        """
        
        xf = []
        yf = []

        c30 = np.cos(np.pi/6)
        s30 = np.sin(np.pi/6)

        for k in range(len(n)):
            y = int(n[k]/L)
            x = 2 * (n[k]%L) + j[k] - y
            x += 1
            xf.append(x * c30)
            if (y%2 == 1):
                if (x%2 == 1):
                    yf.append(y * (1 + s30))
                else:
                    yf.append((y * 3 + 1)/2)
            else:
                if (x%2 == 1):
                    yf.append(1.5 * y + 0.5)
                else:
                    yf.append(1.5 * y)

        return xf, yf


    # finding out which is the smaller cluster (as this is not done automatically by the c++ algorithm)
    """if (len(square_cl[0]) > len(square_cl[1])):
        big = square_cl[0]
        small = square_cl[1]
    else:
        big = square_cl[1]
        small = square_cl[0]
    """

    


    s = []
    v = []
    for i in range(len(square_sv)):
        for j in range(1, len(square_sv[i]), 2):
            s.append(square_sv[i][j])
            v.append(square_sv[i][j-1])

    par = opt.curve_fit(fit_func, s, v)[0]
    print(par)
    """
    x_small = conv(small[2::2], small[3::2])[0]
    y_small = conv(small[2::2], small[3::2])[1]

    x_big = np.asarray(conv(big[2::2], big[3::2])[0])
    y_big = conv(big[2::2], big[3::2])[1]
    """
    fig1, ax1 = plt.subplots()
    """
    
    left, bottom, width, height = [0.15, 0.4, 0.3, 0.3]
    ax2 = fig1.add_axes([left, bottom, width, height])
    ax2.scatter(x=x_small, y=y_small, marker="s", color="black", s=30)
    ax2.axis([min(x_small) - 5, max(x_small) + 5, min(y_small) - 5, max(y_small) + 5])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_color('green')  
    
    left, bottom, width, height = [0.6, 0.2, 0.3, 0.3]
    ax3 = fig1.add_axes([left, bottom, width, height])
    ax3.scatter(x=x_big, y=y_big, marker="s", color="black", s=1e-2)
    ax3.axis([min(x_big) - 1, max(x_big) + 1, min(y_big) - 5, max(y_big) + 5])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_color('blue')  
    """

    #plt.plot(np.linspace(1, 140), func(np.linspace(1, 140)))
    ax1.plot(np.linspace(1, max(s)), fit_func(np.linspace(1, max(s)), par[0]), "k-", label=r"$x^{%g}$" % par[0])
    for i in range(6):
        ax1.plot(square_sv[i*16][1::2], square_sv[i*16][::2], "o", label=r"$t=%d\cdot 10^4$" % (16*i))
    #ax1.plot(small[1], small[0], "gs")
    #ax1.plot(big[1], big[0], "bs")
    ax1.set_xlabel(r"Surface $S$")
    ax1.set_ylabel(r"Volume $V$")
    ax1.grid()
    ax1.legend()
    plt.savefig("hx_n_3_a_0.01.pdf", dpi=150)
    #plt.show()


    fig2, ax1 = plt.subplots()
    """
    left, bottom, width, height = [0.15, 0.4, 0.3, 0.3]
    ax2 = fig2.add_axes([left, bottom, width, height])
    ax2.scatter(x=x_small, y=y_small, marker="s", color="black", s=30)
    ax2.axis([min(x_small) - 5, max(x_small) + 5, min(y_small) - 5, max(y_small) + 5])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.spines[:].set_color('green')  
    
    left, bottom, width, height = [0.6, 0.2, 0.3, 0.3]
    ax3 = fig2.add_axes([left, bottom, width, height])
    ax3.scatter(x=x_big, y=y_big, marker="s", color="black", s=1e-2)
    ax3.axis([min(x_big) - 1, max(x_big) + 1, min(y_big) - 5, max(y_big) + 5])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.spines[:].set_color('blue')  

    """
    #plt.plot(np.linspace(1, 140), func(np.linspace(1, 140)))
    ax1.plot(np.linspace(1, max(s)), fit_func(np.linspace(1, max(s)), par[0]), "k-", label=r"$x^{%g}$" % par[0])
    for i in range(6):
        ax1.plot(square_sv[i*16][1::2], square_sv[i*16][::2], "o", label=r"$t=%d\cdot 10^4$" % (16*i))
    #ax1.plot(small[1], small[0], "gs")
    #ax1.plot(big[1], big[0], "bs")
    ax1.set_xlabel(r"Surface $S$")
    ax1.set_ylabel(r"Volume $V$")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.grid()
    ax1.legend()
    plt.savefig("hx_n_3_log_a_0.01.pdf", dpi=150)
    #plt.show()



    fig3, ax = plt.subplots(nrows=3, ncols=2, figsize=(20, 20))
    plt.tight_layout()
    for i in range(6):
        ax[i//2,i%2].set_title(r"$t=%d\cdot 10^4$" % (16*i))
        ax[i//2,i%2].scatter(x=conv(square_part[16*i][1::4], square_part[16*i][2::4])[0], y=conv(square_part[16*i][1::4], square_part[16*i][2::4])[1], marker="s", color="black", s=4)
        ax[i//2,i%2].set_xticks([])
        ax[i//2,i%2].set_yticks([])
        ax[i//2,i%2].scatter(x=conv(square_int[16*i][::2], square_int[16*i][1::2])[0], y=conv(square_int[16*i][::2], square_int[16*i][1::2])[1], marker="s", color="red", s=4)
    
    plt.savefig("snapshot_hx_n_3_a_0.01.pdf", dpi=200)





# plot for 0.2 density ratio
if what == 4:
    # hex 1
    f = open("hex_sv_low1.txt")
    h1 = []
    for line in f:
        h1.append(np.loadtxt(StringIO(line), dtype=int))
    
    # hex 2
    f = open("hex_sv_low2.txt")
    h2 = []
    for line in f:
        h2.append(np.loadtxt(StringIO(line), dtype=int))

    # hex 3
    f = open("hex_sv_low3.txt")
    h3 = []
    for line in f:
        h3.append(np.loadtxt(StringIO(line), dtype=int))

    
    # tri 1
    f = open("tri_sv_low1.txt")
    t1 = []
    for line in f:
        t1.append(np.loadtxt(StringIO(line), dtype=int))

    # tri 2
    f = open("tri_sv_low2.txt")
    t2 = []
    for line in f:
        t2.append(np.loadtxt(StringIO(line), dtype=int))
    
    # tri 3
    f = open("tri_sv_low3.txt")
    t3 = []
    for line in f:
        t3.append(np.loadtxt(StringIO(line), dtype=int))

    # squ 1
    f = open("square_sv_low1.txt")
    s1 = []
    for line in f:
        s1.append(np.loadtxt(StringIO(line), dtype=int))
    
    # squ 2
    f = open("square_sv_low2.txt")
    s2 = []
    for line in f:
        s2.append(np.loadtxt(StringIO(line), dtype=int))
    
    # squ 3
    f = open("square_sv_low3.txt")
    s3 = []
    for line in f:
        s3.append(np.loadtxt(StringIO(line), dtype=int))

    

    def fill(a):
        s = []
        v = []
        for i in range(len(a)):
            for j in range(1, len(a[i]), 2):
                s.append(a[i][j])
                v.append(a[i][j-1])
        return s, v

    print(np.shape(s1))


    # square 1 plot
    if True:
        ps1 = opt.curve_fit(fit_func, fill(s1)[0], fill(s1)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(s1)[0])), fit_func(np.linspace(1, max(fill(s1)[0])), ps1[0]), "k-", label=r"$x^{%g}$" % ps1[0])
        for i in range(6):
            ax1.plot(s1[i*16][1::2], s1[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/sn1.pdf", dpi=150)

    # square 2 plot
    if True:
        ps2 = opt.curve_fit(fit_func, fill(s2)[0], fill(s2)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(s2)[0])), fit_func(np.linspace(1, max(fill(s2)[0])), ps2[0]), "k-", label=r"$x^{%g}$" % ps2[0])
        for i in range(6):
            ax1.plot(s2[i*16][1::2], s2[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/sn2.pdf", dpi=150)
    
    # square 3 plot
    if True:
        ps3 = opt.curve_fit(fit_func, fill(s3)[0], fill(s3)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(s3)[0])), fit_func(np.linspace(1, max(fill(s3)[0])), ps3[0]), "k-", label=r"$x^{%g}$" % ps3[0])
        for i in range(6):
            ax1.plot(s3[i*16][1::2], s3[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/sn3.pdf", dpi=150)
    
    # tri 1 plot
    if True:
        pt1 = opt.curve_fit(fit_func, fill(t1)[0], fill(t1)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(t1)[0])), fit_func(np.linspace(1, max(fill(t1)[0])), pt1[0]), "k-", label=r"$x^{%g}$" % pt1[0])
        for i in range(6):
            ax1.plot(t1[i*16][1::2], t1[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/tn1.pdf", dpi=150)
    
    # tri 2 plot
    if True:
        pt2 = opt.curve_fit(fit_func, fill(t2)[0], fill(t2)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(t2)[0])), fit_func(np.linspace(1, max(fill(t2)[0])), pt2[0]), "k-", label=r"$x^{%g}$" % pt2[0])
        for i in range(6):
            ax1.plot(t2[i*16][1::2], t2[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/tn2.pdf", dpi=150)
    
    # tri 3 plot
    if True:
        pt3 = opt.curve_fit(fit_func, fill(t3)[0], fill(t3)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(t3)[0])), fit_func(np.linspace(1, max(fill(t3)[0])), pt3[0]), "k-", label=r"$x^{%g}$" % pt3[0])
        for i in range(6):
            ax1.plot(t3[i*16][1::2], t3[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/tn3.pdf", dpi=150)
    
    # hex 1 plot
    if True:
        ph1 = opt.curve_fit(fit_func, fill(h1)[0], fill(h1)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(h1)[0])), fit_func(np.linspace(1, max(fill(h1)[0])), ph1[0]), "k-", label=r"$x^{%g}$" % ph1[0])
        for i in range(6):
            ax1.plot(h1[i*16][1::2], h1[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/hn1.pdf", dpi=150)
    
    # hex 2 plot
    if True:
        ph2 = opt.curve_fit(fit_func, fill(h2)[0], fill(h2)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(h2)[0])), fit_func(np.linspace(1, max(fill(h2)[0])), ph2[0]), "k-", label=r"$x^{%g}$" % ph2[0])
        for i in range(6):
            ax1.plot(h2[i*16][1::2], h2[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/hn2.pdf", dpi=150)
    
    # hex 3 plot
    if True:
        ph3 = opt.curve_fit(fit_func, fill(h3)[0], fill(h3)[1])[0]
        fig1, ax1 = plt.subplots()
        ax1.plot(np.linspace(1, max(fill(h3)[0])), fit_func(np.linspace(1, max(fill(h3)[0])), ph3[0]), "k-", label=r"$x^{%g}$" % ph3[0])
        for i in range(6):
            ax1.plot(h3[i*16][1::2], h3[i*16][::2], "o", label=r"$t=%d\cdot 10^3$" % (16*i))
        ax1.set_xlabel(r"Surface $S$")
        ax1.set_ylabel(r"Volume $V$")
        ax1.grid()
        ax1.legend()
        plt.savefig("./plotsurf/hn3.pdf", dpi=150)
    



