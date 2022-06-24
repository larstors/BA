
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
plt.rcParams.update({'font.size': 22})

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
data_hex_2 = np.genfromtxt("./Data/motility/hexagonal_perc_2.txt", delimiter=" ")
data_hex_3 = np.genfromtxt("./Data/motility/hexagonal_perc_3.txt", delimiter=" ")
data_tri_1 = np.genfromtxt("./Data/motility/triangular_perc.txt", delimiter=" ")
data_tri_2 = np.genfromtxt("./Data/motility/triangular_perc_2.txt", delimiter=" ")
data_tri_3 = np.genfromtxt("./Data/motility/triangular_perc_3.txt", delimiter=" ")
data_squ_1 = np.genfromtxt("./Data/motility/square_perc.txt", delimiter=" ")
data_squ_2 = np.genfromtxt("./Data/motility/square_perc_2.txt", delimiter=" ")
data_squ_3 = np.genfromtxt("./Data/motility/square_perc_3.txt", delimiter=" ")

data_hex_3_det = np.genfromtxt("./Data/motility/hexagonal_perc_details_3.txt", delimiter=" ")
data_tri_3_det = np.genfromtxt("./Data/motility/triangular_perc_details_3.txt", delimiter=" ")
data_squ_3_det = np.genfromtxt("./Data/motility/square_perc_details_3.txt", delimiter=" ")

fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True, sharey=True, figsize=(20,20))
plt.tight_layout()
ax[0, 0].plot(data_squ_1[:, 0], data_squ_1[:, 1], "rs", label="J")
ax[0, 0].plot(data_squ_1[:, 0], data_squ_1[:, 3], "bo", label="M")
ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 1], yerr=np.sqrt(data_squ_1[:, 2]), fmt="s", color="red")
ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 3], yerr=np.sqrt(data_squ_1[:, 4]), fmt="o", color="blue")
ax[0, 0].set_ylabel(r"$J,M$")
#ax[0, 0].set_yscale("log")
ax[0, 0].legend()
ax[0, 0].grid()

ax[0, 1].plot(data_squ_2[:, 0], data_squ_2[:, 1], "rs", label="J")
ax[0, 1].plot(data_squ_2[:, 0], data_squ_2[:, 3], "bo", label="M")
ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 1], yerr=np.sqrt(data_squ_2[:, 2]), fmt="s", color="red")
ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 3], yerr=np.sqrt(data_squ_2[:, 4]), fmt="o", color="blue")
ax[0, 1].legend()
ax[0, 1].grid()

ax[0, 2].plot(data_squ_3[:, 0], data_squ_3[:, 1], "rs", label="J")
ax[0, 2].plot(data_squ_3[:, 0], data_squ_3[:, 3], "bo", label="M")
ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 1], yerr=np.sqrt(data_squ_3[:, 2]), fmt="s", color="red")
ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 3], yerr=np.sqrt(data_squ_3[:, 4]), fmt="o", color="blue")
ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(data_squ_3_det[:, 2]), fmt="s", color="red")
ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(data_squ_3_det[:, 4]), fmt="o", color="blue")
ax[0, 2].legend()
ax[0, 2].grid()

ax[1, 0].plot(data_hex_1[:, 0], data_hex_1[:, 1], "rs", label="J")
ax[1, 0].plot(data_hex_1[:, 0], data_hex_1[:, 3], "bo", label="M")
ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 1], yerr=np.sqrt(data_hex_1[:, 2]), fmt="s", color="red")
ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 3], yerr=np.sqrt(data_hex_1[:, 4]), fmt="o", color="blue")
ax[1, 0].set_ylabel(r"$J,M$")
ax[1, 0].legend()
ax[1, 0].grid()

ax[1, 1].plot(data_hex_2[:, 0], data_hex_2[:, 1], "rs", label="J")
ax[1, 1].plot(data_hex_2[:, 0], data_hex_2[:, 3], "bo", label="M")
ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 1], yerr=np.sqrt(data_hex_2[:, 2]), fmt="s", color="red")
ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 3], yerr=np.sqrt(data_hex_2[:, 4]), fmt="o", color="blue")
ax[1, 1].legend()
ax[1, 1].grid()

ax[1, 2].plot(data_hex_3[:, 0], data_hex_3[:, 1], "rs", label="J")
ax[1, 2].plot(data_hex_3[:, 0], data_hex_3[:, 3], "bo", label="M")
ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 1], yerr=np.sqrt(data_hex_3[:, 2]), fmt="s", color="red")
ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 3], yerr=np.sqrt(data_hex_3[:, 4]), fmt="o", color="blue")
ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(data_hex_3_det[:, 2]), fmt="s", color="red")
ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(data_hex_3_det[:, 4]), fmt="o", color="blue")
ax[1, 2].legend()
ax[1, 2].grid()

ax[2, 0].plot(data_tri_1[:, 0], data_tri_1[:, 1], "rs", label="J")
ax[2, 0].plot(data_tri_1[:, 0], data_tri_1[:, 3], "bo", label="M")
ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 1], yerr=np.sqrt(data_tri_1[:, 2]), fmt="s", color="red")
ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 3], yerr=np.sqrt(data_tri_1[:, 4]), fmt="o", color="blue")
ax[2, 0].set_ylabel(r"$J,M$")
ax[2, 0].set_xlabel(r"$\alpha$")
ax[2, 0].legend()
ax[2, 0].grid()

ax[2, 1].plot(data_tri_2[:, 0], data_tri_2[:, 1], "rs", label="J")
ax[2, 1].plot(data_tri_2[:, 0], data_tri_2[:, 3], "bo", label="M")
ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 1], yerr=np.sqrt(data_tri_2[:, 2]), fmt="s", color="red")
ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 3], yerr=np.sqrt(data_tri_2[:, 4]), fmt="o", color="blue")
ax[2, 1].set_xlabel(r"$\alpha$")
ax[2, 1].legend()
ax[2, 1].grid()

ax[2, 2].plot(data_tri_3[:, 0], data_tri_3[:, 1], "rs", label="J")
ax[2, 2].plot(data_tri_3[:, 0], data_tri_3[:, 3], "bo", label="M")
ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 1], yerr=np.sqrt(data_tri_3[:, 2]), fmt="s", color="red")
ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 3], yerr=np.sqrt(data_tri_3[:, 4]), fmt="o", color="blue")
ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(data_tri_3_det[:, 2]), fmt="s", color="red")
ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(data_tri_3_det[:, 4]), fmt="o", color="blue")
ax[2, 2].set_xlabel(r"$\alpha$")
ax[2, 2].legend()
ax[2, 2].grid()

plt.savefig("./plots/motility_perc.pdf", dpi=150, bbox_inches='tight')


data_hex_1 = np.genfromtxt("./Data/motility/hexagonal_perc_L_50.txt", delimiter=" ")
data_hex_2 = np.genfromtxt("./Data/motility/hexagonal_perc_L_50_2.txt", delimiter=" ")
data_hex_3 = np.genfromtxt("./Data/motility/hexagonal_perc_L_50_3.txt", delimiter=" ")
data_tri_1 = np.genfromtxt("./Data/motility/triangular_perc_L_50.txt", delimiter=" ")
data_tri_2 = np.genfromtxt("./Data/motility/triangular_perc_L_50_2.txt", delimiter=" ")
data_tri_3 = np.genfromtxt("./Data/motility/triangular_perc_L_50_3.txt", delimiter=" ")
data_squ_1 = np.genfromtxt("./Data/motility/square_perc_L_50.txt", delimiter=" ")
data_squ_2 = np.genfromtxt("./Data/motility/square_perc_L_50_2.txt", delimiter=" ")
data_squ_3 = np.genfromtxt("./Data/motility/square_perc_L_50_3.txt", delimiter=" ")
"""
data_hex_3_det = np.genfromtxt("./Data/motility/hexagonal_perc_details_3.txt", delimiter=" ")
data_tri_3_det = np.genfromtxt("./Data/motility/triangular_perc_details_3.txt", delimiter=" ")
data_squ_3_det = np.genfromtxt("./Data/motility/square_perc_details_3.txt", delimiter=" ")
"""
fig, ax = plt.subplots(nrows=3, ncols=3,  sharex=True, sharey=True, figsize=(20,20))
plt.tight_layout()
ax[0, 0].plot(data_squ_1[:, 0], data_squ_1[:, 1], "rs", label="J")
ax[0, 0].plot(data_squ_1[:, 0], data_squ_1[:, 3], "bo", label="M")
ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 1], yerr=np.sqrt(data_squ_1[:, 2]), fmt="s", color="red")
ax[0, 0].errorbar(data_squ_1[:, 0], data_squ_1[:, 3], yerr=np.sqrt(data_squ_1[:, 4]), fmt="o", color="blue")
ax[0, 0].set_ylabel(r"$J,M$")
#ax[0, 0].set_yscale("log")
ax[0, 0].legend()
ax[0, 0].grid()

ax[0, 1].plot(data_squ_2[:, 0], data_squ_2[:, 1], "rs", label="J")
ax[0, 1].plot(data_squ_2[:, 0], data_squ_2[:, 3], "bo", label="M")
ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 1], yerr=np.sqrt(data_squ_2[:, 2]), fmt="s", color="red")
ax[0, 1].errorbar(data_squ_2[:, 0], data_squ_2[:, 3], yerr=np.sqrt(data_squ_2[:, 4]), fmt="o", color="blue")
ax[0, 1].legend()
ax[0, 1].grid()

ax[0, 2].plot(data_squ_3[:, 0], data_squ_3[:, 1], "rs", label="J")
ax[0, 2].plot(data_squ_3[:, 0], data_squ_3[:, 3], "bo", label="M")
ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 1], yerr=np.sqrt(data_squ_3[:, 2]), fmt="s", color="red")
ax[0, 2].errorbar(data_squ_3[:, 0], data_squ_3[:, 3], yerr=np.sqrt(data_squ_3[:, 4]), fmt="o", color="blue")
#ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 1], yerr=np.sqrt(data_squ_3_det[:, 2]), fmt="s", color="red")
#ax[0, 2].errorbar(data_squ_3_det[:, 0], data_squ_3_det[:, 3], yerr=np.sqrt(data_squ_3_det[:, 4]), fmt="o", color="blue")
ax[0, 2].legend()
ax[0, 2].grid()

ax[1, 0].plot(data_hex_1[:, 0], data_hex_1[:, 1], "rs", label="J")
ax[1, 0].plot(data_hex_1[:, 0], data_hex_1[:, 3], "bo", label="M")
ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 1], yerr=np.sqrt(data_hex_1[:, 2]), fmt="s", color="red")
ax[1, 0].errorbar(data_hex_1[:, 0], data_hex_1[:, 3], yerr=np.sqrt(data_hex_1[:, 4]), fmt="o", color="blue")
ax[1, 0].set_ylabel(r"$J,M$")
ax[1, 0].legend()
ax[1, 0].grid()

ax[1, 1].plot(data_hex_2[:, 0], data_hex_2[:, 1], "rs", label="J")
ax[1, 1].plot(data_hex_2[:, 0], data_hex_2[:, 3], "bo", label="M")
ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 1], yerr=np.sqrt(data_hex_2[:, 2]), fmt="s", color="red")
ax[1, 1].errorbar(data_hex_2[:, 0], data_hex_2[:, 3], yerr=np.sqrt(data_hex_2[:, 4]), fmt="o", color="blue")
ax[1, 1].legend()
ax[1, 1].grid()

ax[1, 2].plot(data_hex_3[:, 0], data_hex_3[:, 1], "rs", label="J")
ax[1, 2].plot(data_hex_3[:, 0], data_hex_3[:, 3], "bo", label="M")
ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 1], yerr=np.sqrt(data_hex_3[:, 2]), fmt="s", color="red")
ax[1, 2].errorbar(data_hex_3[:, 0], data_hex_3[:, 3], yerr=np.sqrt(data_hex_3[:, 4]), fmt="o", color="blue")
#ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 1], yerr=np.sqrt(data_hex_3_det[:, 2]), fmt="s", color="red")
#ax[1, 2].errorbar(data_hex_3_det[:, 0], data_hex_3_det[:, 3], yerr=np.sqrt(data_hex_3_det[:, 4]), fmt="o", color="blue")
ax[1, 2].legend()
ax[1, 2].grid()

ax[2, 0].plot(data_tri_1[:, 0], data_tri_1[:, 1], "rs", label="J")
ax[2, 0].plot(data_tri_1[:, 0], data_tri_1[:, 3], "bo", label="M")
ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 1], yerr=np.sqrt(data_tri_1[:, 2]), fmt="s", color="red")
ax[2, 0].errorbar(data_tri_1[:, 0], data_tri_1[:, 3], yerr=np.sqrt(data_tri_1[:, 4]), fmt="o", color="blue")
ax[2, 0].set_ylabel(r"$J,M$")
ax[2, 0].set_xlabel(r"$\alpha$")
ax[2, 0].legend()
ax[2, 0].grid()

ax[2, 1].plot(data_tri_2[:, 0], data_tri_2[:, 1], "rs", label="J")
ax[2, 1].plot(data_tri_2[:, 0], data_tri_2[:, 3], "bo", label="M")
ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 1], yerr=np.sqrt(data_tri_2[:, 2]), fmt="s", color="red")
ax[2, 1].errorbar(data_tri_2[:, 0], data_tri_2[:, 3], yerr=np.sqrt(data_tri_2[:, 4]), fmt="o", color="blue")
ax[2, 1].set_xlabel(r"$\alpha$")
ax[2, 1].legend()
ax[2, 1].grid()

ax[2, 2].plot(data_tri_3[:, 0], data_tri_3[:, 1], "rs", label="J")
ax[2, 2].plot(data_tri_3[:, 0], data_tri_3[:, 3], "bo", label="M")
ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 1], yerr=np.sqrt(data_tri_3[:, 2]), fmt="s", color="red")
ax[2, 2].errorbar(data_tri_3[:, 0], data_tri_3[:, 3], yerr=np.sqrt(data_tri_3[:, 4]), fmt="o", color="blue")
#ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 1], yerr=np.sqrt(data_tri_3_det[:, 2]), fmt="s", color="red")
#ax[2, 2].errorbar(data_tri_3_det[:, 0], data_tri_3_det[:, 3], yerr=np.sqrt(data_tri_3_det[:, 4]), fmt="o", color="blue")
ax[2, 2].set_xlabel(r"$\alpha$")
ax[2, 2].legend()
ax[2, 2].grid()

plt.savefig("./plots/motility_perc_L_50.pdf", dpi=150, bbox_inches='tight')