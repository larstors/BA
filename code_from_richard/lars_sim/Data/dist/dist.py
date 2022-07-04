import numpy as np
import cmath
import matplotlib.pyplot as plt
import scipy as sc

L = 100

s_d = np.genfromtxt("square_dens_1.txt", delimiter=" ")

m = len(s_d)

# 3d matrix with densities
rs = np.ones(m*L*L).reshape(m, L, L)

#filling matrix
for k in range(m):
    for i in range(L):
        for j in range(L):
            rs[k, i, j] = s_d[k, L*i + j]


#matrix for structure factor
S = np.zeros(m*L*L, dtype = 'complex_').reshape(m, L*L)

fs = np.ones(m*L*L, dtype = 'complex_').reshape(m, L, L)

for k in range(m):
    fs[k, :, :] = np.fft.fft2(rs[k, :, :])

Sf = np.zeros(L*L).reshape(L, L)

for k in range(m):
    Sf[:, :] += np.abs(fs[k, :, :])

Sf = 1/m * Sf


Y, X = np.meshgrid(np.arange(0, L, 1), np.arange(0, L, 1))
"""
fig, ax = plt.subplots()

c = ax.pcolormesh(X, Y, np.log(Sf), cmap='viridis')
ax.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax.axis([X.min(), X.max(), Y.min(), Y.max()])
fig.colorbar(c, ax=ax)

plt.show()
"""


f = np.genfromtxt("fourier_test.txt", delimiter=" ")

rs = np.ones(L*L).reshape(L, L)

#filling matrix
for k in range(m):
    for i in range(L):
        for j in range(L):
            rs[i, j] = f[L*i + j]

rs = rs * 1/(4700 * m)

fig1, ax1 = plt.subplots()

c = ax1.pcolormesh(X, Y, rs, cmap='viridis', shading="auto")
ax1.set_title('pcolormesh')
# set the limits of the plot to the limits of the data
ax1.axis([X.min(), X.max(), Y.min(), Y.max()])
fig1.colorbar(c, ax=ax1)

plt.show()

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
    Stot[:] += S[k, :]"""