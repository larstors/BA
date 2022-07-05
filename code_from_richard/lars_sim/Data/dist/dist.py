import numpy as np
import cmath
import matplotlib.pyplot as plt
import scipy as sc

index = int(input())
"""
fig, ax = plt.subplots()
L = 100
c = ax.pcolormesh(X, Y, np.log(Sf), cmap='viridis')
ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax.axis([X.min(), X.max(), Y.min(), Y.max()])
fig.colorbar(c, ax=ax)

plt.show()
"""

if index == 5:
    L = 100
    X, Y = np.meshgrid(np.linspace(0, L, L), np.linspace(0, L, L))
    m = 101
    f = np.genfromtxt("fourier_test.txt", delimiter=" ")

    rs = np.ones(L*L).reshape(L, L)

    #filling matrix
    for i in range(L):
        for j in range(L):
            rs[i, j] = f[L*i + j]

    rs = rs * 1/(4700**2)

    fig1, ax1 = plt.subplots()

    c = ax1.pcolormesh(X, Y, rs, cmap='viridis', shading="auto")
    ax1.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig1.colorbar(c, ax=ax1)

    #plt.show()
    plt.savefig("test1.pdf")


    f1 = np.genfromtxt("fourier_test1.txt", delimiter=" ")

    rs1 = np.ones(L*L).reshape(L, L)

    #filling matrix
    
    for i in range(L):
        for j in range(L):
            rs1[i, j] = f1[L*i + j]

    rs1 = rs1 * 1/(4700**2)

    fig1, ax1 = plt.subplots()

    c = ax1.pcolormesh(X, Y, rs1, cmap='viridis', shading="auto")
    ax1.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig1.colorbar(c, ax=ax1)

    #plt.show()
    plt.savefig("test2.pdf")

    fig1, ax1 = plt.subplots()
    """
    c = ax1.pcolormesh(X, Y, rs/rs1, cmap='viridis', shading="auto")
    ax1.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig1.colorbar(c, ax=ax1)

    #plt.show()
    plt.savefig("test3.pdf")
    """
    print(rs1.max(), rs.max())



"""
# matrix with positions
x = np.array([np.ones(L) * i for i in range(-50, 50, 1)]).T
y = np.array([np.ones(L) * i for i in range(-50, 50, 1)])

r = []

for i in range(-50, 50, 1):
    for j in range(-50, 50, 1):
        r.append(np.array([i, j]))


r = np.asarray(r)


for k in range(m):
    print(k)
    for i in range(L*L):
        if (i%L == 0):
            print("hey")
        for j in range(L*L):
            
            for l in range(L*L):
                S[k, i] += s_d[k, j]*s_d[k, l]*np.exp(-1j * np.dot(r[i,:], r[j,:] - r[l,:]))


Stot = np.ones(L*L)
for k in range(m):
    Stot[:] += S[k, :]
    
    
import numpy as np
import matplotlib.pyplot as plt
import scipy as sci
"""

# for plotting distributions
L = 100

den_sq_3 = np.genfromtxt("square_dens_3.txt", delimiter=" ")
#den_hx_3 = np.genfromtxt("hex_dens_3.txt", delimiter=" ")
#den_tr_3 = np.genfromtxt("tri_dens_3.txt", delimiter=" ")

den_sq_1 = np.genfromtxt("square_dens_1.txt", delimiter=" ")
#den_hx_1 = np.genfromtxt("hex_dens_1.txt", delimiter=" ")
#den_tr_1 = np.genfromtxt("tri_dens_1.txt", delimiter=" ")


ind = 0


if index == 3:
    n = len(den_sq_3)
    tri_dens = np.ones(L*L*n).reshape(n,L,L)
    


    for k in range(n):
        for i in range(L):
            for j in range(L):
                tri_dens[k, i, j] = den_sq_3[k, L*i + j]

    X, Y = np.meshgrid(np.linspace(0, L, L), np.linspace(0, L, L))

    fig, ax =plt.subplots()

    c = ax.pcolormesh(X, Y, tri_dens[-1, :, :], cmap='viridis', shading="auto")
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig("test_3.pdf", dpi=100)


    ff_tri = np.fft.fft2(tri_dens[ind, :, :])
    fig2, ax2 =plt.subplots()

    c = ax2.pcolormesh(X, Y, np.log(1/4700*np.abs(ff_tri)**2), cmap='viridis', shading="auto")
    ax2.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig2.colorbar(c, ax=ax2)
    plt.savefig("fft_test_3.pdf", dpi=100)


if index == 1:
    n = len(den_sq_1)
    l = L
    tri_dens = np.ones(l*l*n).reshape(n,l,l)
    


    for k in range(n):
        for i in range(l):
            for j in range(l):
                tri_dens[k, i, j] = den_sq_1[k, l*i + j]

    X, Y = np.meshgrid(np.linspace(0, l, l), np.linspace(0, l, l))

    fig, ax =plt.subplots()

    c = ax.pcolormesh(X, Y, tri_dens[0, :, :], cmap='viridis', shading="auto")
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig("test_1.pdf", dpi=100)


    ff_tri = np.fft.fft2(tri_dens[ind, :, :])
    fig2, ax2 =plt.subplots()

    c = ax2.pcolormesh(X, Y, np.log(1/4000*np.abs(ff_tri)**2), cmap='viridis', shading="auto")
    ax2.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig2.colorbar(c, ax=ax2)
    plt.savefig("fft_test_1.pdf", dpi=100)


# distributions
if index == 4:
    s1 = np.genfromtxt("square_1.txt", delimiter=" ")
    s2 = np.genfromtxt("square_2.txt", delimiter=" ")
    s3 = np.genfromtxt("square_3.txt", delimiter=" ")

    h1 = np.genfromtxt("hex_1.txt", delimiter=" ")
    h2 = np.genfromtxt("hex_2.txt", delimiter=" ")
    h3 = np.genfromtxt("hex_3.txt", delimiter=" ")

    t1 = np.genfromtxt("tri_1.txt", delimiter=" ")
    t2 = np.genfromtxt("tri_2.txt", delimiter=" ")
    t3 = np.genfromtxt("tri_3.txt", delimiter=" ")

    s1 = s1/np.sum(s1)
    s2 = s2/np.sum(s2)
    s3 = s3/np.sum(s3)

    t1 = t1/np.sum(t1)
    t2 = t2/np.sum(t2)
    t3 = t3/np.sum(t3)

    h1 = h1/np.sum(h1)
    h2 = h2/np.sum(h2)
    h3 = h3/np.sum(h3)




    X = np.arange(0, 5, 1)
    fig1, ax1 = plt.subplots()
    ax1.bar(X - 0.25, s1, color = 'b', width = 0.25)
    ax1.bar(X + 0.00, s2, color = 'g', width = 0.25)
    ax1.bar(X + 0.25, s3, color = 'r', width = 0.25)
    plt.tick_params(axis='x', width=20)
    plt.savefig("square_dist_0.1.pdf", dpi=100)

    X = np.arange(0, 4, 1)
    fig2, ax2 = plt.subplots()
    ax2.bar(X - 0.25, h1, color = 'b', width = 0.25)
    ax2.bar(X + 0.00, h2, color = 'g', width = 0.25)
    ax2.bar(X + 0.25, h3, color = 'r', width = 0.25)
    ax2.set_xticks([0, 1, 2, 3])
    plt.tick_params(axis='x', width=20)
    plt.savefig("hex_dist_0.1.pdf", dpi=100)

    X = np.arange(0, 7, 1)
    fig3, ax3 = plt.subplots()
    ax3.bar(X - 0.25, t1, color = 'b', width = 0.25)
    ax3.bar(X + 0.00, t2, color = 'g', width = 0.25)
    ax3.bar(X + 0.25, t3, color = 'r', width = 0.25)
    plt.tick_params(axis='x', width=20)
    plt.savefig("tri_dist_0.1.pdf", dpi=100)






