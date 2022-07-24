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
den_sq_3 = np.genfromtxt("square_dens_3.txt", delimiter=" ")
#den_hx_3 = np.genfromtxt("hex_dens_3.txt", delimiter=" ")
den_tr_3 = np.genfromtxt("tri_dens_3.txt", delimiter=" ")

den_sq_1 = np.genfromtxt("square_dens_1.txt", delimiter=" ")
#den_hx_1 = np.genfromtxt("hex_dens_1.txt", delimiter=" ")
#den_tr_1 = np.genfromtxt("tri_dens_1.txt", delimiter=" ")

m_sq = len(den_sq_1)
if index == 5:
    L = 100
    X, Y = np.meshgrid(np.linspace(-np.pi, np.pi, L), np.linspace(-np.pi, np.pi, L))
    m = m_sq
    """
    f = np.genfromtxt("fourier_sq_3.txt", delimiter=" ")
    
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
    """


    f1 = np.genfromtxt("fourier_sq_1.txt", delimiter=" ")

    rs1 = np.ones(L*L).reshape(L, L)

    #filling matrix
    
    for i in range(L):
        for j in range(L):
            rs1[i, j] = f1[L*i + j]

    rs1 = rs1 * 1/(4700*m)

    fig1, ax1 = plt.subplots()

    c = ax1.pcolormesh(X, Y, rs1, cmap='viridis', shading="auto")
    ax1.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig1.colorbar(c, ax=ax1)

    #plt.show()
    plt.savefig("ff_sq_n_1.pdf")

    fig1, ax1 = plt.subplots()

    c = ax1.pcolormesh(X, Y, np.log(rs1), cmap='viridis', shading="auto")
    ax1.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig1.colorbar(c, ax=ax1)
    c.set_clim(-10.0, 0)

    #plt.show()
    plt.savefig("ff_sq_n_1_log.pdf")
    print(rs1.max())



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

    test = 1/5800**2*abs(np.fft.fftshift(np.fft.fft2(tri_dens[54, :, :])))**2
    ff_tri = 1/5800**2*abs(np.fft.fftshift(np.fft.fft2(tri_dens[0, :, :])))**2
    ff_tri_n = 1/5800**2*abs(np.fft.fft2(tri_dens[0, :, :]))**2
    for i in range(1, len(tri_dens)):
        ff_tri += 1/5800**2*abs(np.fft.fftshift(np.fft.fft2(tri_dens[i, :, :])))**2
        ff_tri_n += 1/5800**2*abs(np.fft.fft2(tri_dens[i, :, :]))**2

    fig2, ax2 =plt.subplots()
    X, Y = np.meshgrid(np.linspace(0, L, L-3), np.linspace(0, L, L-3))

    #c = ax2.pcolormesh(X, Y, np.log(1/101*ff_tri[1:,1:]), cmap='viridis', shading="auto")
    c = ax2.pcolormesh(X, Y, (1/n*ff_tri[2:-1, 2:-1]), cmap='viridis', shading="auto")
    
    #ax2.axis([30, 70, 30, 70])
    fig2.colorbar(c, ax=ax2)
    plt.savefig("fft_test_3.pdf", dpi=100)


if index == 1:
    n = len(den_tr_1)
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


    ff_tri = 1/4700**2*abs(np.fft.fftshift(np.fft.fft2(tri_dens[0, :, :])))**2

    for i in range(1, len(tri_dens)):
        ff_tri += 1/4700**2*abs(np.fft.fftshift(np.fft.fft2(tri_dens[i, :, :])))**2

    fig2, ax2 =plt.subplots()

    c = ax2.pcolormesh(X, Y, np.log(1/101*ff_tri[1:, 1:]), cmap='viridis', shading="auto")
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


if index == 6:
    a = np.zeros(100*100)
    k = np.zeros(100*100).reshape(100,100)

    def co(i,j):
        x = i - 50
        y = j - 50
        return x, y

    for i in range(100):
        for j in range(100):
            x = co(i, j)[0]
            y = co(i, j)[1]
            if (np.sqrt((x+1)**2 + (y+1)**2) < 5):
            #if (np.abs(x) < 5 and np.abs(y) > 25):
                a[i*100 + j] = 1
            k[i,j] = a[i*100 + j]
    
    a.reshape(1, 100*100)
    np.savetxt("circ.txt", a, delimiter=" ", newline=" ")
    Y, X = np.meshgrid(np.arange(-50, 50, 1), np.arange(-50, 50, 1))
    fig, ax = plt.subplots()
    c = ax.pcolormesh(X, Y, k, cmap='viridis', shading="auto")
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig("circ.pdf", dpi=100)

    fig, ax = plt.subplots()
    c = ax.pcolormesh(X, Y, 1/np.sum(a)**2*np.abs(np.fft.fftshift(np.fft.fft2(k)))**2, cmap='viridis', shading="auto")
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig("ff_circ1.pdf", dpi=100)

    fig, ax = plt.subplots()
    c = ax.pcolormesh(X, Y, 1/np.sum(a)**2*np.abs(np.fft.fft2(k))**2, cmap='viridis', shading="auto")
    ax.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig.colorbar(c, ax=ax)
    plt.savefig("ff_circ2.pdf", dpi=100)


    f1 = np.genfromtxt("fourier_sq_2.txt", delimiter=" ")

    rs1 = np.ones(L*L).reshape(L, L)

    #filling matrix
    
    for i in range(L):
        for j in range(L):
            rs1[i, j] = f1[L*i + j]

    rs1 = rs1

    fig1, ax1 = plt.subplots()

    c = ax1.pcolormesh(X, Y, 1/101*rs1, cmap='viridis', shading="auto")
    ax1.set_title('pcolormesh')
    # set the limits of the plot to the limits of the data
    ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
    fig1.colorbar(c, ax=ax1)

    #plt.show()
    plt.savefig("circ_f.pdf")

