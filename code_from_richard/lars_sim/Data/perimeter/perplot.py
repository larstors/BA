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

name = [r"$J$", r"$M$", r"$c_N$", r"$w_N$"]

#h2 = np.genfromtxt("hex_2.txt", delimiter=" ")
#h3 = np.genfromtxt("hex_3.txt", delimiter=" ")

s2 = np.genfromtxt("square_2.txt", delimiter=" ")
s3 = np.genfromtxt("square_3.txt", delimiter=" ")

t2 = np.genfromtxt("tri_2.txt", delimiter=" ")
t3 = np.genfromtxt("tri_3.txt", delimiter=" ")


n = len(t2)
m = len(t2[0, ::4])

print(t2[:,2::4].max(), t3[:,2::4].max(), s2[:,2::4].max(), np.argmax(s3))

# triangular

rho = np.array([0.001+i*0.03 for i in range(m)])
alp = np.array([np.log10(0.001*1.5**i) for i in range(n)])

y, x = np.meshgrid(rho, alp)


fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, t2[:,0::4])
c1 = ax[1].pcolormesh(x, y, t2[:,1::4])
c2 = ax[2].pcolormesh(x, y, t2[:,2::4])
c = ax[3].pcolormesh(x, y, t2[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")


plt.savefig("test_t2.pdf", bbox_inches='tight')

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, t3[:,0::4])
c1 = ax[1].pcolormesh(x, y, t3[:,1::4])
c2 = ax[2].pcolormesh(x, y, t3[:,2::4])
c = ax[3].pcolormesh(x, y, t3[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")


plt.savefig("test_t3.pdf", bbox_inches='tight')



# square

n = len(s2)
m = len(s2[0, ::4])

rho = np.array([0.001+i*0.034 for i in range(m)])
alp = np.array([np.log10(0.001*1.5**i) for i in range(n)])

print(rho[-3]*100*100*3, 10**alp[-3])

y, x = np.meshgrid(rho, alp)

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, s2[:,0::4])
c1 = ax[1].pcolormesh(x, y, s2[:,1::4])
c2 = ax[2].pcolormesh(x, y, s2[:,2::4])
c = ax[3].pcolormesh(x, y, s2[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")


plt.savefig("test_s2.pdf", bbox_inches='tight')

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, s3[:,0::4])
c1 = ax[1].pcolormesh(x, y, s3[:,1::4])
c2 = ax[2].pcolormesh(x, y, s3[:,2::4])
c = ax[3].pcolormesh(x, y, s3[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")


plt.savefig("test_s3.pdf", bbox_inches='tight')

