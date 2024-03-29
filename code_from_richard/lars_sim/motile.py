
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import sys
import time
import scipy
import matplotlib.animation as animation
#from celluloid import Camera
import matplotlib.tri as mtri
from io import StringIO
import matplotlib.cm as cm
plt.rcParams.update({'font.size': 35})

"""
data_hex_1 = np.genfromtxt("./Data/motility/hexagonal.txt", delimiter=" ")
data_hex_2 = np.genfromtxt("./Data/motility/hexagonal_2.txt", delimiter=" ")
data_hex_3 = np.genfromtxt("./Data/motility/hexagonal_3.txt", delimiter=" ")
data_tri_1 = np.genfromtxt("./Data/motility/triangular.txt", delimiter=" ")
data_tri_2 = np.genfromtxt("./Data/motility/triangular_2.txt", delimiter=" ")
data_tri_3 = np.genfromtxt("./Data/motility/triangular_3.txt", delimiter=" ")
data_squ_1 = np.genfromtxt("./Data/motility/square.txt", delimiter=" ")
data_squ_2 = np.genfromtxt("./Data/motility/square.txt_2.txt", delimiter=" ")
data_squ_3 = np.genfromtxt("./Data/motility/square.txt_3.txt", delimiter=" ")


fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True, sharey=True, figsize=(20,20))
plt.tight_layout()
ax[0, 0].plot(data_squ_1[1::2, 0], data_squ_1[1::2, 1], "rs", label="J")
ax[0, 0].plot(data_squ_1[1::2, 0], data_squ_1[1::2, 3], "bo", label="M")
ax[0, 0].errorbar(data_squ_1[1::2, 0], data_squ_1[1::2, 1], yerr=np.sqrt(data_squ_1[1::2, 2]), fmt="s", color="red")
ax[0, 0].errorbar(data_squ_1[1::2, 0], data_squ_1[1::2, 3], yerr=np.sqrt(data_squ_1[1::2, 4]), fmt="o", color="blue")
ax[0, 0].set_ylabel(r"$J,M$")
ax[0, 0].legend()
#ax[0, 0].set_yscale("log")
ax[0, 0].grid()

ax[0, 1].plot(data_squ_2[1::2, 0], data_squ_2[1::2, 1], "rs", label="J")
ax[0, 1].plot(data_squ_2[1::2, 0], data_squ_2[1::2, 3], "bo", label="M")
ax[0, 1].errorbar(data_squ_2[1::2, 0], data_squ_2[1::2, 1], yerr=np.sqrt(data_squ_2[1::2, 2]), fmt="s", color="red")
ax[0, 1].errorbar(data_squ_2[1::2, 0], data_squ_2[1::2, 3], yerr=np.sqrt(data_squ_2[1::2, 4]), fmt="o", color="blue")
ax[0, 1].legend()
ax[0, 1].grid()

ax[0, 2].plot(data_squ_3[1::2, 0], data_squ_3[1::2, 1], "rs", label="J")
ax[0, 2].plot(data_squ_3[1::2, 0], data_squ_3[1::2, 3], "bo", label="M")
ax[0, 2].errorbar(data_squ_3[1::2, 0], data_squ_3[1::2, 1], yerr=np.sqrt(data_squ_3[1::2, 2]), fmt="s", color="red")
ax[0, 2].errorbar(data_squ_3[1::2, 0], data_squ_3[1::2, 3], yerr=np.sqrt(data_squ_3[1::2, 4]), fmt="o", color="blue")
ax[0, 2].legend()
ax[0, 2].grid()

ax[1, 0].plot(data_hex_1[1::2, 0], data_hex_1[1::2, 1], "rs", label="J")
ax[1, 0].plot(data_hex_1[1::2, 0], data_hex_1[1::2, 3], "bo", label="M")
ax[1, 0].errorbar(data_hex_1[1::2, 0], data_hex_1[1::2, 1], yerr=np.sqrt(data_hex_1[1::2, 2]), fmt="s", color="red")
ax[1, 0].errorbar(data_hex_1[1::2, 0], data_hex_1[1::2, 3], yerr=np.sqrt(data_hex_1[1::2, 4]), fmt="o", color="blue")
ax[1, 0].set_ylabel(r"$J,M$")
ax[1, 0].legend()
ax[1, 0].grid()

ax[1, 1].plot(data_hex_2[1:, 0], data_hex_2[1:, 1], "rs", label="J")
ax[1, 1].plot(data_hex_2[1:, 0], data_hex_2[1:, 3], "bo", label="M")
ax[1, 1].errorbar(data_hex_2[1::2, 0], data_hex_2[1::2, 1], yerr=np.sqrt(data_hex_2[1::2, 2]), fmt="s", color="red")
ax[1, 1].errorbar(data_hex_2[1::2, 0], data_hex_2[1::2, 3], yerr=np.sqrt(data_hex_2[1::2, 4]), fmt="o", color="blue")
ax[1, 1].legend()
ax[1, 1].grid()

ax[1, 2].plot(data_hex_3[1:, 0], data_hex_3[1:, 1], "rs", label="J")
ax[1, 2].plot(data_hex_3[1:, 0], data_hex_3[1:, 3], "bo", label="M")
ax[1, 2].errorbar(data_hex_3[1::2, 0], data_hex_3[1::2, 1], yerr=np.sqrt(data_hex_3[1::2, 2]), fmt="s", color="red")
ax[1, 2].errorbar(data_hex_3[1::2, 0], data_hex_3[1::2, 3], yerr=np.sqrt(data_hex_3[1::2, 4]), fmt="o", color="blue")
ax[1, 2].legend()
ax[1, 2].grid()

ax[2, 0].plot(data_tri_1[1::2, 0], data_tri_1[1::2, 1], "rs", label="J")
ax[2, 0].plot(data_tri_1[1::2, 0], data_tri_1[1::2, 3], "bo", label="M")
ax[2, 0].errorbar(data_tri_1[1::2, 0], data_tri_1[1::2, 1], yerr=np.sqrt(data_tri_1[1::2, 2]), fmt="s", color="red")
ax[2, 0].errorbar(data_tri_1[1::2, 0], data_tri_1[1::2, 3], yerr=np.sqrt(data_tri_1[1::2, 4]), fmt="o", color="blue")
ax[2, 0].set_ylabel(r"$J,M$")
ax[2, 0].set_xlabel(r"$\alpha$")
ax[2, 0].legend()
ax[2, 0].grid()

ax[2, 1].plot(data_tri_2[1::2, 0], data_tri_2[1::2, 1], "rs", label="J")
ax[2, 1].plot(data_tri_2[1::2, 0], data_tri_2[1::2, 3], "bo", label="M")
ax[2, 1].errorbar(data_tri_2[1::2, 0], data_tri_2[1::2, 1], yerr=np.sqrt(data_tri_2[1::2, 2]), fmt="s", color="red")
ax[2, 1].errorbar(data_tri_2[1::2, 0], data_tri_2[1::2, 3], yerr=np.sqrt(data_tri_2[1::2, 4]), fmt="o", color="blue")
ax[2, 1].set_xlabel(r"$\alpha$")
ax[2, 1].legend()
ax[2, 1].grid()

ax[2, 2].plot(data_tri_3[1::2, 0], data_tri_3[1::2, 1], "rs", label="J")
ax[2, 2].plot(data_tri_3[1::2, 0], data_tri_3[1::2, 3], "bo", label="M")
ax[2, 2].errorbar(data_tri_3[1::2, 0], data_tri_3[1::2, 1], yerr=np.sqrt(data_tri_3[1::2, 2]), fmt="s", color="red")
ax[2, 2].errorbar(data_tri_3[1::2, 0], data_tri_3[1::2, 3], yerr=np.sqrt(data_tri_3[1::2, 4]), fmt="o", color="blue")
ax[2, 2].set_xlabel(r"$\alpha$")
ax[2, 2].legend()
ax[2, 2].grid()

plt.savefig("./plots/motility.pdf", dpi=200, bbox_inches='tight')
"""

data_hex_1 = np.genfromtxt("./Data/motility/hexagonal_perc.txt", delimiter=" ")
data_hex_2 = np.genfromtxt(
    "./Data/motility/hexagonal_perc_2.txt", delimiter=" ")
data_hex_3 = np.genfromtxt(
    "./Data/motility/hexagonal_perc_3.txt", delimiter=" ")
data_tri_1 = np.genfromtxt(
    "./Data/motility/triangular_perc.txt", delimiter=" ")
data_tri_2 = np.genfromtxt(
    "./Data/motility/triangular_perc_2.txt", delimiter=" ")
data_tri_3 = np.genfromtxt(
    "./Data/motility/triangular_perc_3.txt", delimiter=" ")
data_squ_1 = np.genfromtxt("./Data/motility/square_perc.txt", delimiter=" ")
data_squ_2 = np.genfromtxt("./Data/motility/square_perc_2.txt", delimiter=" ")
data_squ_3 = np.genfromtxt("./Data/motility/square_perc_3.txt", delimiter=" ")

data_hex_3_det = np.genfromtxt(
    "./Data/motility/hexagonal_perc_details_3.txt", delimiter=" ")
data_tri_3_det = np.genfromtxt(
    "./Data/motility/triangular_perc_details_3.txt", delimiter=" ")
data_squ_3_det = np.genfromtxt(
    "./Data/motility/square_perc_details_3.txt", delimiter=" ")

if False:
    fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                           sharey=True, figsize=(20, 20))
    plt.tight_layout()
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 1], yerr=np.sqrt(
        data_squ_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 3], yerr=np.sqrt(
        data_squ_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 5], yerr=np.sqrt(
        data_squ_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 0].set_ylabel(r"$J,M$")
    #ax[0, 0].set_yscale("log")
    ax[0, 0].legend()
    ax[0, 0].grid()

    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 1], yerr=np.sqrt(
        data_squ_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 3], yerr=np.sqrt(
        data_squ_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 5], yerr=np.sqrt(
        data_squ_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 1].legend()
    ax[0, 1].grid()

    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 1], yerr=np.sqrt(
        data_squ_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 3], yerr=np.sqrt(
        data_squ_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 5], yerr=np.sqrt(
        data_squ_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(
        data_squ_3_det[:, 2]), fmt="s", color="red")
    ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(
        data_squ_3_det[:, 4]), fmt="o", color="blue")
    ax[0, 2].legend()
    ax[0, 2].grid()

    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 1], yerr=np.sqrt(
        data_hex_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 3], yerr=np.sqrt(
        data_hex_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 5], yerr=np.sqrt(
        data_hex_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 0].set_ylabel(r"$J,M$")
    ax[1, 0].legend()
    ax[1, 0].grid()

    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 1], yerr=np.sqrt(
        data_hex_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 3], yerr=np.sqrt(
        data_hex_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 5], yerr=np.sqrt(
        data_hex_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 1].legend()
    ax[1, 1].grid()

    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 1], yerr=np.sqrt(
        data_hex_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 3], yerr=np.sqrt(
        data_hex_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 5], yerr=np.sqrt(
        data_hex_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(
        data_hex_3_det[:, 2]), fmt="s", color="red")
    ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(
        data_hex_3_det[:, 4]), fmt="o", color="blue")
    ax[1, 2].legend()
    ax[1, 2].grid()

    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 1], yerr=np.sqrt(
        data_tri_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 3], yerr=np.sqrt(
        data_tri_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 5], yerr=np.sqrt(
        data_tri_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 0].set_ylabel(r"$J,M$")
    ax[2, 0].set_xlabel(r"$\alpha$")
    ax[2, 0].legend()
    ax[2, 0].grid()

    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 1], yerr=np.sqrt(
        data_tri_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 3], yerr=np.sqrt(
        data_tri_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 5], yerr=np.sqrt(
        data_tri_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 1].set_xlabel(r"$\alpha$")
    ax[2, 1].legend()
    ax[2, 1].grid()

    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 1], yerr=np.sqrt(
        data_tri_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 3], yerr=np.sqrt(
        data_tri_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 5], yerr=np.sqrt(
        data_tri_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(
        data_tri_3_det[:, 2]), fmt="s", color="red")
    ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(
        data_tri_3_det[:, 4]), fmt="o", color="blue")
    ax[2, 2].set_xlabel(r"$\alpha$")
    ax[2, 2].legend()
    ax[2, 2].grid()
    plt.savefig("./plots/motility_perc.pdf", dpi=150, bbox_inches='tight')


"""
ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(data_squ_3_det[:, 2]), fmt="s", color="red")
ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(data_squ_3_det[:, 4]), fmt="o", color="blue")
ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(data_hex_3_det[:, 2]), fmt="s", color="red")
ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(data_hex_3_det[:, 4]), fmt="o", color="blue")
ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(data_tri_3_det[:, 2]), fmt="s", color="red")
ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(data_tri_3_det[:, 4]), fmt="o", color="blue")
"""


data_hex_1_L_50 = np.genfromtxt(
    "./Data/motility/hexagonal_perc_L_50.txt", delimiter=" ")
data_hex_2_L_50 = np.genfromtxt(
    "./Data/motility/hexagonal_perc_L_50_2.txt", delimiter=" ")
data_hex_3_L_50 = np.genfromtxt(
    "./Data/motility/hexagonal_perc_L_50_3.txt", delimiter=" ")
data_tri_1_L_50 = np.genfromtxt(
    "./Data/motility/triangular_perc_L_50.txt", delimiter=" ")
data_tri_2_L_50 = np.genfromtxt(
    "./Data/motility/triangular_perc_L_50_2.txt", delimiter=" ")
data_tri_3_L_50 = np.genfromtxt(
    "./Data/motility/triangular_perc_L_50_3.txt", delimiter=" ")
data_squ_1_L_50 = np.genfromtxt(
    "./Data/motility/square_perc_L_50.txt", delimiter=" ")
data_squ_2_L_50 = np.genfromtxt(
    "./Data/motility/square_perc_L_50_2.txt", delimiter=" ")
data_squ_3_L_50 = np.genfromtxt(
    "./Data/motility/square_perc_L_50_3.txt", delimiter=" ")
"""
data_hex_3_det = np.genfromtxt("./Data/motility/hexagonal_perc_details_3.txt", delimiter=" ")
data_tri_3_det = np.genfromtxt("./Data/motility/triangular_perc_details_3.txt", delimiter=" ")
data_squ_3_det = np.genfromtxt("./Data/motility/square_perc_details_3.txt", delimiter=" ")
"""
if False:
    fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                           sharey=True, figsize=(20, 20))
    plt.tight_layout()
    ax[0, 0].errorbar(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 1], yerr=np.sqrt(
        data_squ_1_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 0].errorbar(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 3], yerr=np.sqrt(
        data_squ_1_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 0].errorbar(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 5], yerr=np.sqrt(
        data_squ_1_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 0].set_ylabel(r"$J,M$")
    #ax[0, 0].set_yscale("log")
    ax[0, 0].legend()
    ax[0, 0].grid()

    ax[0, 1].errorbar(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 1], yerr=np.sqrt(
        data_squ_2_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 1].errorbar(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 3], yerr=np.sqrt(
        data_squ_2_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 1].errorbar(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5], yerr=np.sqrt(
        data_squ_2_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 1].legend()
    ax[0, 1].grid()

    ax[0, 2].errorbar(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 1], yerr=np.sqrt(
        data_squ_3_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 2].errorbar(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 3], yerr=np.sqrt(
        data_squ_3_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 2].errorbar(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5], yerr=np.sqrt(
        data_squ_3_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(data_squ_3_det[:, 2]), fmt="s", color="red")
    #ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(data_squ_3_det[:, 4]), fmt="o", color="blue")
    ax[0, 2].legend()
    ax[0, 2].grid()

    ax[1, 0].errorbar(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 1], yerr=np.sqrt(
        data_hex_1_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 0].errorbar(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 3], yerr=np.sqrt(
        data_hex_1_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 0].errorbar(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5], yerr=np.sqrt(
        data_hex_1_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 0].set_ylabel(r"$J,M$")
    ax[1, 0].legend()
    ax[1, 0].grid()

    ax[1, 1].errorbar(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 1], yerr=np.sqrt(
        data_hex_2_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 1].errorbar(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 3], yerr=np.sqrt(
        data_hex_2_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 1].errorbar(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5], yerr=np.sqrt(
        data_hex_2_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 1].legend()
    ax[1, 1].grid()

    ax[1, 2].errorbar(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 1], yerr=np.sqrt(
        data_hex_3_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 2].errorbar(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 3], yerr=np.sqrt(
        data_hex_3_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 2].errorbar(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5], yerr=np.sqrt(
        data_hex_3_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(data_hex_3_det[:, 2]), fmt="s", color="red")
    #ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(data_hex_3_det[:, 4]), fmt="o", color="blue")
    ax[1, 2].legend()
    ax[1, 2].grid()

    ax[2, 0].errorbar(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 1], yerr=np.sqrt(
        data_tri_1_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 0].errorbar(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 3], yerr=np.sqrt(
        data_tri_1_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 0].errorbar(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5], yerr=np.sqrt(
        data_tri_1_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 0].set_ylabel(r"$J,M$")
    ax[2, 0].set_xlabel(r"$\alpha$")
    ax[2, 0].legend()
    ax[2, 0].grid()

    ax[2, 1].errorbar(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 1], yerr=np.sqrt(
        data_tri_2_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 1].errorbar(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 3], yerr=np.sqrt(
        data_tri_2_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 1].errorbar(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5], yerr=np.sqrt(
        data_tri_2_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 1].set_xlabel(r"$\alpha$")
    ax[2, 1].legend()
    ax[2, 1].grid()

    ax[2, 2].errorbar(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 1], yerr=np.sqrt(
        data_tri_3_L_50[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 2].errorbar(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 3], yerr=np.sqrt(
        data_tri_3_L_50[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 2].errorbar(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5], yerr=np.sqrt(
        data_tri_3_L_50[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(data_tri_3_det[:, 2]), fmt="s", color="red")
    #ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(data_tri_3_det[:, 4]), fmt="o", color="blue")
    ax[2, 2].set_xlabel(r"$\alpha$")
    ax[2, 2].legend()
    ax[2, 2].grid()

    plt.savefig("./plots/motility_perc_L_50.pdf", dpi=150, bbox_inches='tight')


# w comp
if False:

    fig2, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                            figsize=(30, 30), constrained_layout=True)
    # plt.tight_layout()
    ax2 = ax[0, 0].twinx()
    ax2.errorbar(data_squ_1[:, 0], data_squ_1[:, 5]*4700, yerr=np.sqrt(
        data_squ_1[:, 6])*4700*0, fmt="s", color="red", label=r"$w, L=100$")
    ax2.errorbar(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 5]*1200, yerr=np.sqrt(
        data_squ_1_L_50[:, 6])*1200*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[0, 0].errorbar(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 5], yerr=np.sqrt(
        data_squ_1_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 5], yerr=np.sqrt(
        data_squ_1[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[0, 0].set_ylabel(r"$w_N$")
    ax2.set_ylabel(r"$w$")
    ax2.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[0, 0].legend(loc=1)
    ax[0, 0].grid()

    ax3 = ax[0, 1].twinx()
    ax3.errorbar(data_squ_2[:, 0], data_squ_2[:, 5]*5500, yerr=np.sqrt(
        data_squ_2[:, 6])*5500*0, fmt="s", color="red", label=r"$w, L=100$")
    ax3.errorbar(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5]*1400, yerr=np.sqrt(
        data_squ_2_L_50[:, 6])*1400*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[0, 1].errorbar(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5], yerr=np.sqrt(
        data_squ_2_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 5], yerr=np.sqrt(
        data_squ_2[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[0, 1].set_ylabel(r"$w_N$")
    ax3.set_ylabel(r"$w$")
    ax3.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[0, 1].legend(loc=1)
    ax[0, 1].grid()

    ax4 = ax[0, 2].twinx()
    ax4.errorbar(data_squ_3[:, 0], data_squ_3[:, 5]*5800, yerr=np.sqrt(
        data_squ_3[:, 6])*5800*0, fmt="s", color="red", label=r"$w, L=100$")
    ax4.errorbar(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5]*1500, yerr=np.sqrt(
        data_squ_3_L_50[:, 6])*1500*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[0, 2].errorbar(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5], yerr=np.sqrt(
        data_squ_3_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 5], yerr=np.sqrt(
        data_squ_3[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[0, 2].set_ylabel(r"$w_N$")
    ax4.set_ylabel(r"$w$")
    ax4.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[0, 2].legend(loc=1)
    ax[0, 2].grid()

    ax5 = ax[1, 0].twinx()
    ax5.errorbar(data_hex_1[:, 0], data_hex_1[:, 5]*11000, yerr=np.sqrt(
        data_hex_1[:, 6])*11000*0, fmt="s", color="red", label=r"$w, L=100$")
    ax5.errorbar(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5]*2800, yerr=np.sqrt(
        data_hex_1_L_50[:, 6])*2800*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[1, 0].errorbar(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5], yerr=np.sqrt(
        data_hex_1_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 5], yerr=np.sqrt(
        data_hex_1[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[1, 0].set_ylabel(r"$w_N$")
    ax5.set_ylabel(r"$w$")
    ax5.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[1, 0].legend(loc=1)
    ax[1, 0].grid()

    ax6 = ax[1, 1].twinx()
    ax6.errorbar(data_hex_2[:, 0], data_hex_2[:, 5]*13200, yerr=np.sqrt(
        data_hex_2[:, 6])*13200*0, fmt="s", color="red", label=r"$w, L=100$")
    ax6.errorbar(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5]*3300, yerr=np.sqrt(
        data_hex_2_L_50[:, 6])*3300*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[1, 1].errorbar(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5], yerr=np.sqrt(
        data_hex_2_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 5], yerr=np.sqrt(
        data_hex_2[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[1, 1].set_ylabel(r"$w_N$")
    ax6.set_ylabel(r"$w$")
    ax6.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[1, 1].legend(loc=1)
    ax[1, 1].grid()

    ax7 = ax[1, 2].twinx()
    ax7.errorbar(data_hex_3[:, 0], data_hex_3[:, 5]*14100, yerr=np.sqrt(
        data_hex_3[:, 6])*14100*0, fmt="s", color="red", label=r"$w, L=100$")
    ax7.errorbar(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5]*3500, yerr=np.sqrt(
        data_hex_3_L_50[:, 6])*3500*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[1, 2].errorbar(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5], yerr=np.sqrt(
        data_hex_3_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 5], yerr=np.sqrt(
        data_hex_3[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[1, 2].set_ylabel(r"$w_N$")
    ax7.set_ylabel(r"$w$")
    ax7.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[1, 2].legend(loc=1)
    ax[1, 2].grid()

    ax8 = ax[2, 0].twinx()
    ax8.errorbar(data_tri_1[:, 0], data_tri_1[:, 5]*4000, yerr=np.sqrt(
        data_tri_1[:, 6])*4000*0, fmt="s", color="red", label=r"$w, L=100$")
    ax8.errorbar(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5]*1000, yerr=np.sqrt(
        data_tri_1_L_50[:, 6])*1000*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[2, 0].errorbar(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5], yerr=np.sqrt(
        data_tri_1_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 5], yerr=np.sqrt(
        data_tri_1[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[2, 0].set_ylabel(r"$w_N$")
    ax8.set_ylabel(r"$w$")
    ax8.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[2, 0].legend(loc=1)
    ax[2, 0].grid()

    ax9 = ax[2, 1].twinx()
    ax9.errorbar(data_tri_2[:, 0], data_tri_2[:, 5]*4500, yerr=np.sqrt(
        data_tri_2[:, 6])*4500*0, fmt="s", color="red", label=r"$w, L=100$")
    ax9.errorbar(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5]*1100, yerr=np.sqrt(
        data_tri_2_L_50[:, 6])*1100*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[2, 1].errorbar(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5], yerr=np.sqrt(
        data_tri_2_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 5], yerr=np.sqrt(
        data_tri_2[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[2, 1].set_ylabel(r"$w_N$")
    ax9.set_ylabel(r"$w$")
    ax9.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[2, 1].legend(loc=1)
    ax[2, 1].grid()

    ax21 = ax[2, 2].twinx()
    ax21.errorbar(data_tri_3[:, 0], data_tri_3[:, 5]*4700, yerr=np.sqrt(
        data_tri_3[:, 6])*4700*0, fmt="s", color="red", label=r"$w, L=100$")
    ax21.errorbar(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5]*1200, yerr=np.sqrt(
        data_tri_3_L_50[:, 6])*1200*0, fmt="x", color="red", label=r"$w, L=50$")
    ax[2, 2].errorbar(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5], yerr=np.sqrt(
        data_tri_3_L_50[:, 6])*0, fmt="x", color="black", label=r"$w_N, L=50$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 5], yerr=np.sqrt(
        data_tri_3[:, 6])*0, fmt="s", color="black", label=r"$w_N, L=100$")
    ax[2, 2].set_ylabel(r"$w_N$")
    ax21.set_ylabel(r"$w$")
    ax21.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    ax[2, 2].legend(loc=1)
    ax[2, 2].grid()

    plt.savefig("./plots/w_comp.pdf", dpi=150, bbox_inches='tight')

# M comp
if False:
    fig2, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                            figsize=(30, 30), constrained_layout=True)
    # plt.tight_layout()
    ax2 = ax[0, 0].twinx()
    ax2.plot(data_squ_1[:, 0], data_squ_1[:, 5-2]*4700,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax2.plot(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 5-2]
             * 1200, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[0, 0].plot(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[0, 0].plot(data_squ_1[:, 0], data_squ_1[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    ax[0, 0].set_ylabel(r"$M$")
    # ax2.set_ylabel(r"$NM$")
    #ax2.legend(loc=7, bbox_to_anchor=(1.05, 1))
    #ax[0, 0].set_yscale("log")
    #ax[0, 0].legend(loc=1, bbox_to_anchor=(1.05, 1))
    ax[0, 0].grid()

    ax3 = ax[0, 1].twinx()
    ax3.plot(data_squ_2[:, 0], data_squ_2[:, 5-2]*5500,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax3.plot(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5-2]
             * 1400, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[0, 1].plot(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[0, 1].plot(data_squ_2[:, 0], data_squ_2[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    #ax[0, 1].set_ylabel(r"$M$")
    # ax3.set_ylabel(r"$NM$")
    # ax3.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[0, 1].legend(loc=1)
    ax[0, 1].grid()

    ax4 = ax[0, 2].twinx()
    ax4.plot(data_squ_3[:, 0], data_squ_3[:, 5-2]*5800,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax4.plot(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5-2]
             * 1500, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[0, 2].plot(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[0, 2].plot(data_squ_3[:, 0], data_squ_3[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    #ax[0, 2].set_ylabel(r"$M$")
    ax4.set_ylabel(r"$NM$")
    # ax4.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[0, 2].legend(loc=1)
    ax[0, 2].grid()

    ax5 = ax[1, 0].twinx()
    ax5.plot(data_hex_1[:, 0], data_hex_1[:, 5-2]*11000,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax5.plot(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5-2]
             * 2800, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[1, 0].plot(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[1, 0].plot(data_hex_1[:, 0], data_hex_1[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    ax[1, 0].set_ylabel(r"$M$")
    # ax5.set_ylabel(r"$NM$")
    # ax5.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[1, 0].legend(loc=1)
    ax[1, 0].grid()

    ax6 = ax[1, 1].twinx()
    ax6.plot(data_hex_2[:, 0], data_hex_2[:, 5-2]*13200,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax6.plot(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5-2]
             * 3300, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[1, 1].plot(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[1, 1].plot(data_hex_2[:, 0], data_hex_2[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    ##ax[1, 1].set_ylabel(r"$M$")
    # ax6.set_ylabel(r"$NM$")
    # ax6.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[1, 1].legend(loc=1)
    ax[1, 1].grid()

    ax7 = ax[1, 2].twinx()
    ax7.plot(data_hex_3[:, 0], data_hex_3[:, 5-2]*14100,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax7.plot(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5-2]
             * 3500, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[1, 2].plot(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[1, 2].plot(data_hex_3[:, 0], data_hex_3[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    #ax[1, 2].set_ylabel(r"$M$")
    ax7.set_ylabel(r"$NM$")
    # ax7.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[1, 2].legend(loc=1)
    ax[1, 2].grid()

    ax8 = ax[2, 0].twinx()
    ax8.plot(data_tri_1[:, 0], data_tri_1[:, 5-2]*4000,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax8.plot(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5-2]
             * 1000, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[2, 0].plot(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[2, 0].plot(data_tri_1[:, 0], data_tri_1[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    ax[2, 0].set_ylabel(r"$M$")
    ax[2, 0].set_xlabel(r"$\alpha$")
    # ax8.set_ylabel(r"$NM$")
    # ax8.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[2, 0].legend(loc=1)
    ax[2, 0].grid()

    ax9 = ax[2, 1].twinx()
    ax9.plot(data_tri_2[:, 0], data_tri_2[:, 5-2]*4500,
             "d", color="red", label=r"$NM, L=100$", ms=10)
    ax9.plot(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5-2]
             * 1100, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[2, 1].plot(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[2, 1].plot(data_tri_2[:, 0], data_tri_2[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    #ax[2, 1].set_ylabel(r"$M$")
    # ax9.set_ylabel(r"$NM$")
    ax[2, 1].set_xlabel(r"$\alpha$")
    # ax9.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[2, 1].legend(loc=1)
    ax[2, 1].grid()

    ax21 = ax[2, 2].twinx()
    ax21.plot(data_tri_3[:, 0], data_tri_3[:, 5-2]*4700,
              "d", color="red", label=r"$NM, L=100$", ms=10)
    ax21.plot(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5-2]
              * 1200, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[2, 2].plot(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5-2],
                  "x", color="black", label=r"$M, L=50$", ms=10)
    ax[2, 2].plot(data_tri_3[:, 0], data_tri_3[:, 5-2], "s",
                  color="black", label=r"$M, L=100$", ms=10)
    #ax[2, 2].set_ylabel(r"$M$")
    ax[2, 2].set_xlabel(r"$\alpha$")
    ax21.set_ylabel(r"$NM$")
    # ax21.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[2, 2].legend(loc=1)
    ax[2, 2].grid()

    plt.savefig("./plots/M_comp.pdf", dpi=300, bbox_inches='tight')

# M comp 1
if False:
    fig2, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                            figsize=(30, 30), constrained_layout=True)
    # plt.tight_layout()
    #ax2 = ax[0,0].twinx()
    #ax2.plot(data_squ_1[:, 0], data_squ_1[:, 5-2]*4700, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax2.plot(data_squ_1_L_50[:, 0], data_squ_1_L_50[:, 5-2]*1200, "x", color="red", label=r"$NM, L=50$", ms=10)
    ax[0, 0].plot(data_squ_1_L_50[:, 0], data_squ_1[:, 5-2] -
                  data_squ_1_L_50[:, 5-2], "s-", color="black", label=r"$M$", ms=10)
    ax[0, 0].plot(data_squ_1_L_50[:, 0], (data_squ_1[:, 5-2]*4700-data_squ_1_L_50[:,
                                                                                  5-2]*1200)/(4700-1200), "s-", color="red", label=r"$NM$", ms=10)
    #ax[0, 0].plot(data_squ_1[:, 0], data_squ_1[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[0, 0].set_ylabel(r"$\Delta M, \Delta NM$")
    # ax2.set_ylabel(r"$NM$")
    #ax2.legend(loc=7, bbox_to_anchor=(1.05, 1))
    #ax[0, 0].set_yscale("log")
    ax[0, 0].legend()
    ax[0, 0].grid()

    #ax3 = ax[0,1].twinx()
    ax[0, 1].plot(data_squ_2_L_50[:, 0], data_squ_2[:, 5-2]-data_squ_2_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[0, 1].plot(data_squ_2_L_50[:, 0], (data_squ_2[:, 5-2]*5500-data_squ_2_L_50[:, 5-2]*1400) /
                  (5500-1400), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    #ax3.plot(data_squ_2[:, 0], data_squ_2[:, 5-2]*5500, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax3.plot(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5-2]*1400, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[0, 1].plot(data_squ_2_L_50[:, 0], data_squ_2_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[0, 1].plot(data_squ_2[:, 0], data_squ_2[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    #ax[0, 1].set_ylabel(r"$M$")
    # ax3.set_ylabel(r"$NM$")
    # ax3.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[0, 1].legend(loc=1)
    ax[0, 1].grid()

    #ax4 = ax[0,2].twinx()
    #ax4.plot(data_squ_3[:, 0], data_squ_3[:, 5-2]*5800, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax4.plot(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5-2]*1500, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[0, 2].plot(data_squ_3_L_50[:, 0], data_squ_3_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[0, 2].plot(data_squ_3[:, 0], data_squ_3[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[0, 2].plot(data_squ_3_L_50[:, 0], data_squ_3[:, 5-2]-data_squ_3_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[0, 2].plot(data_squ_3_L_50[:, 0], (data_squ_3[:, 5-2]*5800-data_squ_3_L_50[:, 5-2]*1500) /
                  (5800-1500), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    #ax[0, 2].set_ylabel(r"$M$")
    # ax4.set_ylabel(r"$NM$")
    # ax4.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[0, 2].legend(loc=1)
    ax[0, 2].grid()

    #ax5 = ax[1,0].twinx()
    #ax5.plot(data_hex_1[:, 0], data_hex_1[:, 5-2]*11000, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax5.plot(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5-2]*2800, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[1, 0].plot(data_hex_1_L_50[:, 0], data_hex_1_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[1, 0].plot(data_hex_1[:, 0], data_hex_1[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[1, 0].plot(data_hex_1_L_50[:, 0], data_hex_1[:, 5-2]-data_hex_1_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[1, 0].plot(data_hex_1_L_50[:, 0], (data_hex_1[:, 5-2]*11000-data_hex_1_L_50[:, 5-2]*2800) /
                  (11000-2800), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    ax[1, 0].set_ylabel(r"$\Delta M, \Delta NM$")
    # ax5.set_ylabel(r"$NM$")
    # ax5.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[1, 0].legend(loc=1)
    ax[1, 0].grid()

    #ax6 = ax[1,1].twinx()
    #ax6.plot(data_hex_2[:, 0], data_hex_2[:, 5-2]*13200, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax6.plot(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5-2]*3300, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[1, 1].plot(data_hex_2_L_50[:, 0], data_hex_2_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[1, 1].plot(data_hex_2[:, 0], data_hex_2[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[1, 1].plot(data_hex_2_L_50[:, 0], data_hex_2[:, 5-2]-data_hex_2_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[1, 1].plot(data_hex_2_L_50[:, 0], (data_hex_2[:, 5-2]*13200-data_hex_2_L_50[:, 5-2]*3300) /
                  (13200-3300), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    ##ax[1, 1].set_ylabel(r"$M$")
    # ax6.set_ylabel(r"$NM$")
    # ax6.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[1, 1].legend(loc=1)
    ax[1, 1].grid()

    #ax7 = ax[1,2].twinx()
    #ax7.plot(data_hex_3[:, 0], data_hex_3[:, 5-2]*14100, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax7.plot(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5-2]*3500, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[1, 2].plot(data_hex_3_L_50[:, 0], data_hex_3_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[1, 2].plot(data_hex_3[:, 0], data_hex_3[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[1, 2].plot(data_hex_3_L_50[:, 0], data_hex_3[:, 5-2]-data_hex_3_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[1, 2].plot(data_hex_3_L_50[:, 0], (data_hex_3[:, 5-2]*14100-data_hex_3_L_50[:, 5-2]*3500) /
                  (14100-3500), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    #ax[1, 2].set_ylabel(r"$M$")
    # ax7.set_ylabel(r"$NM$")
    # ax7.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[1, 2].legend(loc=1)
    ax[1, 2].grid()

    #ax8 = ax[2,0].twinx()
    #ax8.plot(data_tri_1[:, 0], data_tri_1[:, 5-2]*4000, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax8.plot(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5-2]*1000, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[2, 0].plot(data_tri_1_L_50[:, 0], data_tri_1_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[2, 0].plot(data_tri_1[:, 0], data_tri_1[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[2, 0].plot(data_tri_1_L_50[:, 0], data_tri_1[:, 5-2]-data_tri_1_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[2, 0].plot(data_tri_1_L_50[:, 0], (data_tri_1[:, 5-2]*4000-data_tri_1_L_50[:, 5-2]*1000) /
                  (4000-1000), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    ax[2, 0].set_ylabel(r"$\Delta M, \Delta NM$")
    ax[2, 0].set_xlabel(r"$\alpha$")
    # ax8.set_ylabel(r"$NM$")
    # ax8.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[2, 0].legend(loc=1)
    ax[2, 0].grid()

    #ax9 = ax[2,1].twinx()
    #ax9.plot(data_tri_2[:, 0], data_tri_2[:, 5-2]*4500, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax9.plot(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5-2]*1100, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[2, 1].plot(data_tri_2_L_50[:, 0], data_tri_2_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[2, 1].plot(data_tri_2[:, 0], data_tri_2[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    ax[2, 1].plot(data_tri_2_L_50[:, 0], data_tri_2[:, 5-2]-data_tri_2_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[2, 1].plot(data_tri_2_L_50[:, 0], (data_tri_2[:, 5-2]*4500-data_tri_2_L_50[:, 5-2]*1100) /
                  (4500-1100), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    #ax[2, 1].set_ylabel(r"$M$")
    # ax9.set_ylabel(r"$NM$")
    ax[2, 1].set_xlabel(r"$\alpha$")
    # ax9.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[2, 1].legend(loc=1)
    ax[2, 1].grid()

    #ax21 = ax[2,2].twinx()
    #ax21.plot(data_tri_3[:, 0], data_tri_3[:, 5-2]*4700, "d", color="red", label=r"$NM, L=100$", ms=10)
    #ax21.plot(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5-2]*1200, "x", color="red", label=r"$NM, L=50$", ms=10)
    #ax[2, 2].plot(data_tri_3_L_50[:, 0], data_tri_3_L_50[:, 5-2], "x", color="black", label=r"$M, L=50$", ms=10)
    #ax[2, 2].plot(data_tri_3[:, 0], data_tri_3[:, 5-2], "s", color="black", label=r"$M, L=100$", ms=10)
    #ax[2, 2].set_ylabel(r"$M$")
    ax[2, 2].plot(data_tri_3_L_50[:, 0], data_tri_3[:, 5-2]-data_tri_3_L_50[:,
                                                                            5-2], "s-", color="black", label=r"$M_{100}-M_{50}$", ms=10)
    ax[2, 2].plot(data_tri_3_L_50[:, 0], (data_tri_3[:, 5-2]*4700-data_tri_3_L_50[:, 5-2]*1200) /
                  (4700-1200), "s-", color="red", label=r"$(M_{100}N_{100}-M_{50}N_{50})/(N_{100}-N_{50})$", ms=10)
    ax[2, 2].set_xlabel(r"$\alpha$")
    # ax21.set_ylabel(r"$NM$")
    # ax21.legend(loc=7)
    #ax[0, 0].set_yscale("log")
    #ax[2, 2].legend(loc=1)
    ax[2, 2].grid()

    plt.savefig("./plots/M_comp_diff.pdf", dpi=300, bbox_inches='tight')

# Hysteresis
if True:
    plt.rcParams.update({'font.size': 40})
    hex_for_L_50 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_fhyst_L_50_3.txt", delimiter=" ")
    hex_back_L_50 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_bhyst_L_50_3.txt", delimiter=" ")
    hex_for_L_100 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_fhyst_3.txt", delimiter=" ")
    hex_back_L_100 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_bhyst_3.txt", delimiter=" ")
    tri_for_L_50 = np.genfromtxt(
        "./Data/motility/triangular_perc_fhyst_L_50_3.txt", delimiter=" ")
    tri_back_L_50 = np.genfromtxt(
        "./Data/motility/triangular_perc_bhyst_L_50_3.txt", delimiter=" ")
    tri_for_L_100 = np.genfromtxt(
        "./Data/motility/triangular_perc_fhyst_3.txt", delimiter=" ")
    tri_back_L_100 = np.genfromtxt(
        "./Data/motility/triangular_perc_bhyst_3.txt", delimiter=" ")
    squ_for_L_50 = np.genfromtxt(
        "./Data/motility/square_perc_fhyst_L_50_3.txt", delimiter=" ")
    squ_back_L_50 = np.genfromtxt(
        "./Data/motility/square_perc_bhyst_L_50_3.txt", delimiter=" ")
    squ_for_L_100 = np.genfromtxt(
        "./Data/motility/square_perc_fhyst_3.txt", delimiter=" ")
    squ_back_L_100 = np.genfromtxt(
        "./Data/motility/square_perc_bhyst_3.txt", delimiter=" ")

    marks = 15

    fig3, ax = plt.subplots(nrows=2, ncols=3, figsize=(30, 20))
    plt.tight_layout()
    ax[0, 0].set_title("Hexagonal")
    ax[0, 0].plot(hex_for_L_50[:, 0], hex_for_L_50[:, 1], ">--",
                  color="red", label=r"$J$, forward", markersize=marks)
    #ax[0,0].plot(hex_for_L_50[:, 0], hex_for_L_50[:, 3], "s", color="green", label=r"$M$, forward", markersize=marks)
    ax[0, 0].plot(hex_back_L_50[:, 0], hex_back_L_50[:, 1], "*--",
                  color="black", label=r"$J$, backward", markersize=marks)
    #ax[0,0].plot(hex_back_L_50[:, 0], hex_back_L_50[:, 3],  "o", color="red", label=r"$M$, backward", markersize=marks)
    ax[0, 0].grid()
    ax[0, 0].set_ylabel(r"$J,M$")
    ax[0, 0].legend()

    ax[1, 0].plot(hex_for_L_100[:, 0], hex_for_L_100[:, 1],  ">--",
                  color="red", label=r"$J$, forward", markersize=marks)
    #ax[1,0].plot(hex_for_L_100[:, 0], hex_for_L_100[:, 3],  "s", color="green", label=r"$M$, forward", markersize=marks)
    ax[1, 0].plot(hex_back_L_100[:, 0], hex_back_L_100[:, 1],  "*--",
                  color="black", label=r"$J$, backward", markersize=marks)
    #ax[1,0].plot(hex_back_L_100[:, 0], hex_back_L_100[:, 3],  "o", color="red", label=r"$M$, backward", markersize=marks)
    ax[1, 0].set_ylabel(r"$J,M$")
    ax[1, 0].grid()
    ax[1, 0].set_xlabel(r"$\alpha$")

    ax[0, 1].set_title("Triangular")
    ax[0, 1].plot(tri_for_L_50[:, 0], tri_for_L_50[:, 1],  ">--",
                  color="red", label=r"$J$, forward", markersize=marks)
    #ax[0,1].plot(tri_for_L_50[:, 0], tri_for_L_50[:, 3],  "s", color="green", label=r"$M$, forward", markersize=marks)
    ax[0, 1].plot(tri_back_L_50[:, 0], tri_back_L_50[:, 1],  "*--",
                  color="black", label=r"$J$, backward", markersize=marks)
    #ax[0,1].plot(tri_back_L_50[:, 0], tri_back_L_50[:, 3],  "o", color="red", label=r"$M$, backward", markersize=marks)
    ax[0, 1].grid()

    ax[1, 1].plot(tri_for_L_100[:, 0], tri_for_L_100[:, 1],  ">--",
                  color="red", label=r"$J$, forward", markersize=marks)
    #ax[1,1].plot(tri_for_L_100[:, 0], tri_for_L_100[:, 3],  "s", color="green", label=r"$M$, forward", markersize=marks)
    ax[1, 1].plot(tri_back_L_100[:, 0], tri_back_L_100[:, 1],  "*--",
                  color="black", label=r"$J$, backward", markersize=marks)
    #ax[1,1].plot(tri_back_L_100[:, 0], tri_back_L_100[:, 3],  "o", color="red", label=r"$M$, backward", markersize=marks)
    ax[1, 1].set_xlabel(r"$\alpha$")
    ax[1, 1].grid()

    ax[0, 2].set_title("Square")
    ax[0, 2].plot(squ_for_L_50[:, 0], squ_for_L_50[:, 1],  ">--",
                  color="red", label=r"$J$, forward", markersize=marks)
    #ax[0,2].plot(squ_for_L_50[:, 0], squ_for_L_50[:, 3],  "s", color="green", label=r"$M$, forward", markersize=marks)
    ax[0, 2].plot(squ_back_L_50[:, 0], squ_back_L_50[:, 1],  "*--",
                  color="black", label=r"$J$, backward", markersize=marks)
    #ax[0,2].plot(squ_back_L_50[:, 0], squ_back_L_50[:, 3],  "o", color="red", label=r"$M$, backward", markersize=marks)
    ax[0, 2].grid()

    ax[1, 2].plot(squ_for_L_100[:, 0], squ_for_L_100[:, 1],  ">--",
                  color="red", label=r"$J$, forward", markersize=marks)
    #ax[1,2].plot(squ_for_L_100[:, 0], squ_for_L_100[:, 3],  "s", color="green", label=r"$M$, forward", markersize=marks)
    ax[1, 2].plot(squ_back_L_100[:, 0], squ_back_L_100[:, 1],  "*--",
                  color="black", label=r"$J$, backward", markersize=marks)
    #ax[1,2].plot(squ_back_L_100[:, 0], squ_back_L_100[:, 3],  "o", color="red", label=r"$M$, backward", markersize=marks)
    ax[1, 2].set_xlabel(r"$\alpha$")
    ax[1, 2].grid()

    for a, color in zip([ax[0, 0], ax[1, 0], ax[0, 1], ax[1, 1], ax[0, 2], ax[1, 2]], ['black', 'black', 'blue', 'blue', 'blue', 'blue']):
        for ticks in a.xaxis.get_ticklines() + a.yaxis.get_ticklines():
            ticks.set_color(color)
        for pos in ['top', 'bottom', 'right', 'left']:
            a.spines[pos].set_edgecolor(color)

    plt.savefig("./plots/hyst_1.svg", dpi=500, bbox_inches='tight')


"""
The following is for 0.2 percolation threshold!

"""


if False:
    data_hex_1 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_low.txt", delimiter=" ")
    data_hex_2 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_low_2.txt", delimiter=" ")
    data_hex_3 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_low_3.txt", delimiter=" ")
    data_tri_1 = np.genfromtxt(
        "./Data/motility/triangular_perc_low.txt", delimiter=" ")
    data_tri_2 = np.genfromtxt(
        "./Data/motility/triangular_perc_low_2.txt", delimiter=" ")
    data_tri_3 = np.genfromtxt(
        "./Data/motility/triangular_perc_low_3.txt", delimiter=" ")
    data_squ_1 = np.genfromtxt(
        "./Data/motility/square_perc_low.txt", delimiter=" ")
    data_squ_2 = np.genfromtxt(
        "./Data/motility/square_perc_low_2.txt", delimiter=" ")
    data_squ_3 = np.genfromtxt(
        "./Data/motility/square_perc_low_3.txt", delimiter=" ")
    fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                           sharey=True, figsize=(20, 20))
    plt.tight_layout()
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 1], yerr=np.sqrt(
        data_squ_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 3], yerr=np.sqrt(
        data_squ_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 5], yerr=np.sqrt(
        data_squ_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 0].set_ylabel(r"$J,M$")
    #ax[0, 0].set_yscale("log")
    ax[0, 0].legend()
    ax[0, 0].grid()
    ax[0, 0].axis([-0.005, 0.055, -0.005, 1.005])

    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 1], yerr=np.sqrt(
        data_squ_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 3], yerr=np.sqrt(
        data_squ_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 5], yerr=np.sqrt(
        data_squ_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 1].legend()
    ax[0, 1].grid()

    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 1], yerr=np.sqrt(
        data_squ_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 3], yerr=np.sqrt(
        data_squ_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 5], yerr=np.sqrt(
        data_squ_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(data_squ_3_det[:, 2]), fmt="s", color="red")
    #ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(data_squ_3_det[:, 4]), fmt="o", color="blue")
    ax[0, 2].legend()
    ax[0, 2].grid()

    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 1], yerr=np.sqrt(
        data_hex_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 3], yerr=np.sqrt(
        data_hex_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 5], yerr=np.sqrt(
        data_hex_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 0].set_ylabel(r"$J,M$")
    ax[1, 0].legend()
    ax[1, 0].grid()

    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 1], yerr=np.sqrt(
        data_hex_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 3], yerr=np.sqrt(
        data_hex_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 5], yerr=np.sqrt(
        data_hex_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 1].legend()
    ax[1, 1].grid()

    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 1], yerr=np.sqrt(
        data_hex_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 3], yerr=np.sqrt(
        data_hex_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 5], yerr=np.sqrt(
        data_hex_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(data_hex_3_det[:, 2]), fmt="s", color="red")
    #ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(data_hex_3_det[:, 4]), fmt="o", color="blue")
    ax[1, 2].legend()
    ax[1, 2].grid()

    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 1], yerr=np.sqrt(
        data_tri_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 3], yerr=np.sqrt(
        data_tri_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 5], yerr=np.sqrt(
        data_tri_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 0].set_ylabel(r"$J,M$")
    ax[2, 0].set_xlabel(r"$\alpha$")
    ax[2, 0].legend()
    ax[2, 0].grid()

    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 1], yerr=np.sqrt(
        data_tri_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 3], yerr=np.sqrt(
        data_tri_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 5], yerr=np.sqrt(
        data_tri_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 1].set_xlabel(r"$\alpha$")
    ax[2, 1].legend()
    ax[2, 1].grid()

    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 1], yerr=np.sqrt(
        data_tri_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 3], yerr=np.sqrt(
        data_tri_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 5], yerr=np.sqrt(
        data_tri_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(data_tri_3_det[:, 2]), fmt="s", color="red")
    #ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(data_tri_3_det[:, 4]), fmt="o", color="blue")
    ax[2, 2].set_xlabel(r"$\alpha$")
    ax[2, 2].legend()
    ax[2, 2].grid()
    plt.savefig("./plots/motility_perc_low.pdf", dpi=150, bbox_inches='tight')


if False:
    data_hex_1 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_low_L_50.txt", delimiter=" ")
    data_hex_2 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_low_L_50_2.txt", delimiter=" ")
    data_hex_3 = np.genfromtxt(
        "./Data/motility/hexagonal_perc_low_L_50_3.txt", delimiter=" ")
    data_tri_1 = np.genfromtxt(
        "./Data/motility/triangular_perc_low_L_50.txt", delimiter=" ")
    data_tri_2 = np.genfromtxt(
        "./Data/motility/triangular_perc_low_L_50_2.txt", delimiter=" ")
    data_tri_3 = np.genfromtxt(
        "./Data/motility/triangular_perc_low_L_50_3.txt", delimiter=" ")
    data_squ_1 = np.genfromtxt(
        "./Data/motility/square_perc_low_L_50.txt", delimiter=" ")
    data_squ_2 = np.genfromtxt(
        "./Data/motility/square_perc_low_L_50_2.txt", delimiter=" ")
    data_squ_3 = np.genfromtxt(
        "./Data/motility/square_perc_low_L_50_3.txt", delimiter=" ")
    fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True,
                           sharey=True, figsize=(20, 20))
    plt.tight_layout()
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 1], yerr=np.sqrt(
        data_squ_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 3], yerr=np.sqrt(
        data_squ_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 5], yerr=np.sqrt(
        data_squ_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 0].set_ylabel(r"$J,M$")
    #ax[0, 0].set_yscale("log")
    ax[0, 0].legend()
    ax[0, 0].grid()
    ax[0, 0].axis([-0.005, 0.055, -0.005, 1.005])

    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 1], yerr=np.sqrt(
        data_squ_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 3], yerr=np.sqrt(
        data_squ_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 5], yerr=np.sqrt(
        data_squ_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[0, 1].legend()
    ax[0, 1].grid()

    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 1], yerr=np.sqrt(
        data_squ_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 3], yerr=np.sqrt(
        data_squ_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 5], yerr=np.sqrt(
        data_squ_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(data_squ_3_det[:, 2]), fmt="s", color="red")
    #ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(data_squ_3_det[:, 4]), fmt="o", color="blue")
    ax[0, 2].legend()
    ax[0, 2].grid()

    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 1], yerr=np.sqrt(
        data_hex_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 3], yerr=np.sqrt(
        data_hex_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 5], yerr=np.sqrt(
        data_hex_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 0].set_ylabel(r"$J,M$")
    ax[1, 0].legend()
    ax[1, 0].grid()

    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 1], yerr=np.sqrt(
        data_hex_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 3], yerr=np.sqrt(
        data_hex_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 5], yerr=np.sqrt(
        data_hex_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[1, 1].legend()
    ax[1, 1].grid()

    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 1], yerr=np.sqrt(
        data_hex_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 3], yerr=np.sqrt(
        data_hex_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 5], yerr=np.sqrt(
        data_hex_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(data_hex_3_det[:, 2]), fmt="s", color="red")
    #ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(data_hex_3_det[:, 4]), fmt="o", color="blue")
    ax[1, 2].legend()
    ax[1, 2].grid()

    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 1], yerr=np.sqrt(
        data_tri_1[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 3], yerr=np.sqrt(
        data_tri_1[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 5], yerr=np.sqrt(
        data_tri_1[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 0].set_ylabel(r"$J,M$")
    ax[2, 0].set_xlabel(r"$\alpha$")
    ax[2, 0].legend()
    ax[2, 0].grid()

    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 1], yerr=np.sqrt(
        data_tri_2[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 3], yerr=np.sqrt(
        data_tri_2[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 5], yerr=np.sqrt(
        data_tri_2[:, 6]), fmt="x", color="black", label=r"$w_N$")
    ax[2, 1].set_xlabel(r"$\alpha$")
    ax[2, 1].legend()
    ax[2, 1].grid()

    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 1], yerr=np.sqrt(
        data_tri_3[:, 2]), fmt="s", color="red", label=r"$J$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 3], yerr=np.sqrt(
        data_tri_3[:, 4]), fmt="o", color="blue", label=r"$M$")
    ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 5], yerr=np.sqrt(
        data_tri_3[:, 6]), fmt="x", color="black", label=r"$w_N$")
    #ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(data_tri_3_det[:, 2]), fmt="s", color="red")
    #ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(data_tri_3_det[:, 4]), fmt="o", color="blue")
    ax[2, 2].set_xlabel(r"$\alpha$")
    ax[2, 2].legend()
    ax[2, 2].grid()
    plt.savefig("./plots/motility_perc_low_L_50.pdf",
                dpi=150, bbox_inches='tight')
