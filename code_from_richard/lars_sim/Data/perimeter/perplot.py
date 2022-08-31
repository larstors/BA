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
import matplotlib
matplotlib.use('Agg')

plt.rcParams.update({'font.size': 17})

name = [r"$J$", r"$M$", r"$c_N$", r"$w_N$"]

h1 = np.genfromtxt("hex_1.txt", delimiter=" ")
h2 = np.genfromtxt("hex_2.txt", delimiter=" ")
h3 = np.genfromtxt("hex_3.txt", delimiter=" ")

s1 = np.genfromtxt("square_1.txt", delimiter=" ")
s2 = np.genfromtxt("square_2.txt", delimiter=" ")
s3 = np.genfromtxt("square_3.txt", delimiter=" ")

t1 = np.genfromtxt("tri_1.txt", delimiter=" ")
t2 = np.genfromtxt("tri_2.txt", delimiter=" ")
t3 = np.genfromtxt("tri_3.txt", delimiter=" ")

def sigm(x, a, b, c):
    return a * (1 / (1+ np.exp(-(x-b)/c)))

def pol(x, m, a, c):
    return a*(x-m)**c


n = len(t1)
m = len(t1[0, ::4])

# triangular
#####################################################################################
rho = np.array([0.001+i*.05 for i in range(m)])
alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])

tfit = np.linspace(min(alp), max(alp), 100)

y, x = np.meshgrid(rho, alp)

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, t1[:,0::4])
c1 = ax[1].pcolormesh(x, y, t1[:,1::4])
c2 = ax[2].pcolormesh(x, y, t1[:,2::4])
c = ax[3].pcolormesh(x, y, t1[:,3::4])

for i in range(4):
    ax[i].plot(tfit, sigm(tfit, 0.475, min(tfit)+.8, .05)+.05, "--", color="red", linewidth=2)

fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_t1.pdf", bbox_inches='tight')

###################################################################################

n = len(t2)
m = len(t2[0, ::4])

rho = np.array([0.001+i*.05 for i in range(m)])
alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])

y, x = np.meshgrid(rho, alp)

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, t2[:,0::4])
c1 = ax[1].pcolormesh(x, y, t2[:,1::4])
c2 = ax[2].pcolormesh(x, y, t2[:,2::4])
c = ax[3].pcolormesh(x, y, t2[:,3::4])

for i in range(4):
    ax[i].plot(np.ones(2)*(-0.75), np.array([0.4, 1]), color="red", linewidth=2)
    ax[i].plot(np.array([-0.4, max(alp)]), np.ones(2)*0.35, "--", color="red", linewidth=2)
    ax[i].plot(np.array([-0.75,-0.4]), np.array([0.4, 0.35]), color="red", linewidth=2)
    ax[i].plot(np.linspace(-3, -0.75), pol(np.linspace(-3, -0.75), -3, 0.012, 5)+0.05, color="red", linewidth=2)

fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")


plt.savefig("perp_t2.pdf", bbox_inches='tight')

##########################################################################################

n = len(t3)
m = len(t3[0, ::4])


rho = np.array([0.001+i*.05 for i in range(m)])
alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])

y, x = np.meshgrid(rho, alp)

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, t3[:,0::4])
c1 = ax[1].pcolormesh(x, y, t3[:,1::4])
c2 = ax[2].pcolormesh(x, y, t3[:,2::4])
c = ax[3].pcolormesh(x, y, t3[:,3::4])
fig.colorbar(c, ax=ax)

for i in range(4):
    ax[i].plot(np.ones(2)*(-0.8), np.array([0.25, 1]), color="red", linewidth=2)
    ax[i].plot(np.array([-0.8, max(alp)]), np.ones(2)*0.25, "--", color="red", linewidth=2)
    ax[i].plot(np.linspace(-3, -0.8), pol(np.linspace(-3, -0.8), -3, 0.012, 5)+0.05, color="red", linewidth=2)

for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")


plt.savefig("perp_t3.pdf", bbox_inches='tight')



# square
##########################################################################################
n = len(s1)
m = len(s1[0, ::4])

rho = np.array([0.001+i*.05 for i in range(m)])
alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])


y, x = np.meshgrid(rho, alp)


fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, s1[:,0::4])
c1 = ax[1].pcolormesh(x, y, s1[:,1::4])
c2 = ax[2].pcolormesh(x, y, s1[:,2::4])
c = ax[3].pcolormesh(x, y, s1[:,3::4])

for i in range(4):
    ax[i].plot(tfit, sigm(tfit, 0.15, min(tfit)+.9, .2)+.45, "--", color="red", linewidth=2)

fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_s1.pdf", bbox_inches='tight')
##########################################################################################
n = len(s2)
m = len(s2[0, ::4])

rho = np.array([0.001+i*.05 for i in range(m)])
alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])


y, x = np.meshgrid(rho, alp)


fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, s2[:,0::4])
c1 = ax[1].pcolormesh(x, y, s2[:,1::4])
c2 = ax[2].pcolormesh(x, y, s2[:,2::4])
c = ax[3].pcolormesh(x, y, s2[:,3::4])
fig.colorbar(c, ax=ax)

for i in range(4):
    ax[i].plot(np.ones(2)*(-0.75), np.array([0.45, 1]), color="red", linewidth=2)
    ax[i].plot(np.array([-0.4, max(alp)]), np.ones(2)*0.4, "--", color="red", linewidth=2)
    ax[i].plot(np.array([-0.75,-0.4]), np.array([0.45, 0.4]), color="red", linewidth=2)
    ax[i].plot(np.linspace(-3, -0.75), pol(np.linspace(-3, -0.75), -3, 0.012, 5)+0.05, color="red", linewidth=2)

for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_s2.pdf", bbox_inches='tight')
##########################################################################################
n = len(s3)
m = len(s3[0, ::4])

rho = np.array([0.001+i*.05 for i in range(m)])
alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])
print(10**alp[-1], rho[7]*100*100*3)

y, x = np.meshgrid(rho, alp)

fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, s3[:,0::4])
c1 = ax[1].pcolormesh(x, y, s3[:,1::4])
c2 = ax[2].pcolormesh(x, y, s3[:,2::4])
c = ax[3].pcolormesh(x, y, s3[:,3::4])

for i in range(4):
    ax[i].plot(np.ones(2)*(-0.75), np.array([0.3, 1]), color="red", linewidth=2)
    ax[i].plot(np.array([-0.75, max(alp)]), np.ones(2)*0.3, "--", color="red", linewidth=2)
    ax[i].plot(np.linspace(-3, -0.75), pol(np.linspace(-3, -0.75), -3, 0.012, 5)+0.05, color="red", linewidth=2)

fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_s3.pdf", bbox_inches='tight')



# hexagonal
##########################################################################################
n = len(h1)
m = len(h1[0, ::4])

rho = np.array([0.001+i*.08 for i in range(m)])
alp = np.array([np.log10(0.001*2.4**i) for i in range(n)])

y, x = np.meshgrid(rho, alp)
fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, h1[:,0::4])
c1 = ax[1].pcolormesh(x, y, h1[:,1::4])
c2 = ax[2].pcolormesh(x, y, h1[:,2::4])
c = ax[3].pcolormesh(x, y, h1[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_h1.pdf", bbox_inches='tight')
##########################################################################################

n = len(h2)
m = len(h2[0, ::4])

rho = np.array([0.001+i*.08 for i in range(m)])
alp = np.array([np.log10(0.001*2.4**i) for i in range(n)])

y, x = np.meshgrid(rho, alp)
fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, h2[:,0::4])
c1 = ax[1].pcolormesh(x, y, h2[:,1::4])
c2 = ax[2].pcolormesh(x, y, h2[:,2::4])
c = ax[3].pcolormesh(x, y, h2[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_h2.pdf", bbox_inches='tight')

##########################################################################################
n = len(h3)
m = len(h3[0, ::4])

rho = np.array([0.001+i*.08 for i in range(m)])
alp = np.array([np.log10(0.001*2.4**i) for i in range(n)])

y, x = np.meshgrid(rho, alp)
fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
plt.tight_layout()
c0 = ax[0].pcolormesh(x, y, h3[:,0::4])
c1 = ax[1].pcolormesh(x, y, h3[:,1::4])
c2 = ax[2].pcolormesh(x, y, h3[:,2::4])
c = ax[3].pcolormesh(x, y, h3[:,3::4])
fig.colorbar(c, ax=ax)
for i in range(4):
    ax[i].axis([x.min(), x.max(), y.min(), y.max()])
    ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
    ax[i].set_title(name[i])
ax[0].set_ylabel(r"$\rho$")
plt.savefig("perp_h3.pdf", bbox_inches='tight')
##########################################################################################



if True:
    n = len(t1)
    m = len(t1[0, ::4])

    # triangular
    #####################################################################################
    rho = np.array([0.001+i*.05 for i in range(m)])
    alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])

    tfit = np.linspace(min(alp), max(alp), 100)

    y, x = np.meshgrid(rho, alp)

    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, t1[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, t1[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, t1[:,2::4], levels=20)
    c = ax[3].contourf(x, y, t1[:,3::4], levels=20)

    for i in range(4):
        ax[i].plot(tfit, sigm(tfit, 0.475, min(tfit)+.8, .05)+.05, "--", color="red", linewidth=2)

    fig.colorbar(c, ax=ax)
    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_t1_contour.pdf", bbox_inches='tight')

    ###################################################################################

    n = len(t2)
    m = len(t2[0, ::4])

    rho = np.array([0.001+i*.05 for i in range(m)])
    alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])

    y, x = np.meshgrid(rho, alp)

    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, t2[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, t2[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, t2[:,2::4], levels=20)
    c = ax[3].contourf(x, y, t2[:,3::4], levels=20)

    for i in range(4):
        ax[i].plot(np.ones(2)*(-0.75), np.array([0.4, 1]), color="red", linewidth=2)
        ax[i].plot(np.array([-0.4, max(alp)]), np.ones(2)*0.35, "--", color="red", linewidth=2)
        ax[i].plot(np.array([-0.75,-0.4]), np.array([0.4, 0.35]), color="red", linewidth=2)
        ax[i].plot(np.linspace(-3, -0.75), pol(np.linspace(-3, -0.75), -3, 0.012, 5)+0.05, color="red", linewidth=2)

    fig.colorbar(c, ax=ax)
    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")


    plt.savefig("perp_t2_contour.pdf", bbox_inches='tight')

    ##########################################################################################

    n = len(t3)
    m = len(t3[0, ::4])


    rho = np.array([0.001+i*.05 for i in range(m)])
    alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])

    y, x = np.meshgrid(rho, alp)

    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, t3[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, t3[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, t3[:,2::4], levels=20)
    c = ax[3].contourf(x, y, t3[:,3::4], levels=20)
    fig.colorbar(c, ax=ax)

    for i in range(4):
        ax[i].plot(np.ones(2)*(-0.8), np.array([0.25, 1]), color="red", linewidth=2)
        ax[i].plot(np.array([-0.8, max(alp)]), np.ones(2)*0.25, "--", color="red", linewidth=2)
        ax[i].plot(np.linspace(-3, -0.8), pol(np.linspace(-3, -0.8), -3, 0.012, 5)+0.05, color="red", linewidth=2)

    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")


    plt.savefig("perp_t3_contour.pdf", bbox_inches='tight')



    # square
    ##########################################################################################
    n = len(s1)
    m = len(s1[0, ::4])

    rho = np.array([0.001+i*.05 for i in range(m)])
    alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])


    y, x = np.meshgrid(rho, alp)


    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, s1[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, s1[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, s1[:,2::4], levels=20)
    c = ax[3].contourf(x, y, s1[:,3::4], levels=20)

    for i in range(4):
        ax[i].plot(tfit, sigm(tfit, 0.15, min(tfit)+.9, .2)+.45, color="red", linewidth=2)

    fig.colorbar(c, ax=ax)
    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_s1_contour.pdf", bbox_inches='tight')
    ##########################################################################################
    n = len(s2)
    m = len(s2[0, ::4])

    rho = np.array([0.001+i*.05 for i in range(m)])
    alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])


    y, x = np.meshgrid(rho, alp)


    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, s2[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, s2[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, s2[:,2::4], levels=20)
    c = ax[3].contourf(x, y, s2[:,3::4], levels=20)
    fig.colorbar(c, ax=ax)

    for i in range(4):
        ax[i].plot(np.ones(2)*(-0.75), np.array([0.45, 1]), color="red", linewidth=2)
        ax[i].plot(np.array([-0.4, max(alp)]), np.ones(2)*0.4, "--", color="red", linewidth=2)
        ax[i].plot(np.array([-0.75,-0.4]), np.array([0.45, 0.4]), color="red", linewidth=2)
        ax[i].plot(np.linspace(-3, -0.75), pol(np.linspace(-3, -0.75), -3, 0.012, 5)+0.05, color="red", linewidth=2)

    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_s2_contour.pdf", bbox_inches='tight')
    ##########################################################################################
    n = len(s3)
    m = len(s3[0, ::4])

    rho = np.array([0.001+i*.05 for i in range(m)])
    alp = np.array([np.log10(0.001*1.78**i) for i in range(n)])
    print(10**alp[-1], rho[7]*100*100*3)

    y, x = np.meshgrid(rho, alp)

    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, s3[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, s3[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, s3[:,2::4], levels=20)
    c = ax[3].contourf(x, y, s3[:,3::4], levels=20)
    fig.colorbar(c, ax=ax)

    for i in range(4):
        ax[i].plot(np.ones(2)*(-0.75), np.array([0.3, 1]), color="red", linewidth=2)
        ax[i].plot(np.array([-0.75, max(alp)]), np.ones(2)*0.3, "--", color="red", linewidth=2)
        ax[i].plot(np.linspace(-3, -0.75), pol(np.linspace(-3, -0.75), -3, 0.012, 5)+0.05, color="red", linewidth=2)

    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_s3_contour.pdf", bbox_inches='tight')



    # hexagonal
    ##########################################################################################
    n = len(h1)
    m = len(h1[0, ::4])

    rho = np.array([0.001+i*.08 for i in range(m)])
    alp = np.array([np.log10(0.001*2.4**i) for i in range(n)])

    y, x = np.meshgrid(rho, alp)
    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, h1[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, h1[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, h1[:,2::4], levels=20)
    c = ax[3].contourf(x, y, h1[:,3::4], levels=20)
    fig.colorbar(c, ax=ax)

    for i in range(4):
        ax[i].plot(tfit, sigm(tfit, 0.15, min(tfit)+.9, .2)+.45, color="red", linewidth=2)

    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_h1_contour.pdf", bbox_inches='tight')
    ##########################################################################################

    n = len(h2)
    m = len(h2[0, ::4])

    rho = np.array([0.001+i*.08 for i in range(m)])
    alp = np.array([np.log10(0.001*2.4**i) for i in range(n)])

    y, x = np.meshgrid(rho, alp)
    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, h2[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, h2[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, h2[:,2::4], levels=20)
    c = ax[3].contourf(x, y, h2[:,3::4], levels=20)
    fig.colorbar(c, ax=ax)
    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_h2_contour.pdf", bbox_inches='tight')

    ##########################################################################################
    n = len(h3)
    m = len(h3[0, ::4])

    rho = np.array([0.001+i*.08 for i in range(m)])
    alp = np.array([np.log10(0.001*2.4**i) for i in range(n)])

    y, x = np.meshgrid(rho, alp)
    fig, ax = plt.subplots(ncols=4, nrows=1, figsize=(15, 5), sharex=True, sharey=True)
    plt.tight_layout()
    c0 = ax[0].contourf(x, y, h3[:,0::4], levels=20)
    c1 = ax[1].contourf(x, y, h3[:,1::4], levels=20)
    c2 = ax[2].contourf(x, y, h3[:,2::4], levels=20)
    c = ax[3].contourf(x, y, h3[:,3::4], levels=20)
    fig.colorbar(c, ax=ax)
    for i in range(4):
        ax[i].axis([x.min(), x.max(), y.min(), y.max()])
        ax[i].set_xlabel(r"$\log_{10}(\alpha)$")
        ax[i].set_title(name[i])
    ax[0].set_ylabel(r"$\rho$")
    plt.savefig("perp_h3_contour.pdf", bbox_inches='tight')
    ##########################################################################################

