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

class distribution:
    def __init__(self):
        #self.data_tri = np.genfromtxt("tridist.txt", delimiter=" ")
        #self.data_tri_nr = np.genfromtxt("tridist_nr.txt", delimiter=" ")
        #self.data_hex = np.genfromtxt("hexdist.txt", delimiter=" ")
        #self.data_hex_nr = np.genfromtxt("hexdist_nr.txt", delimiter=" ")
        self.work = True
    
    def fit(self, x, b, c):
        N = -b * np.log(x) - c*x
        return N - np.log(np.sum(np.exp(N)))

    def plot_particles_tri(self):
        # for normalization
        self.fig, (self.ax1, self.ax2) =  plt.subplots(nrows=1, ncols=2, figsize=(20, 20))
        y = self.data_tri[0::2]
        y_nr = self.data_tri_nr[0::2]
        sum_hist = np.sum(y)
        sum_hist_nr = np.sum(y_nr)
        x = np.argwhere(y != 0)
        x_nr = np.argwhere(y_nr != 0)

        y = y[x] / sum_hist
        y_nr = y_nr[x_nr] / sum_hist_nr

        n = 100
        n_nr = 100

        p = opt.curve_fit(self.fit, x[:n, 0]+1, np.log(y[:n, 0]), [2.0, 5.0/max(x)], bounds=(0, np.inf))[0]
        p_nr = opt.curve_fit(self.fit, x_nr[:n_nr, 0]+1, np.log(y_nr[:n_nr, 0]), [2.0, 5.0/max(x_nr)], bounds=(0, np.inf))[0]

        self.ax1.loglog(x+1, y, "r-", label="hist")
        self.ax1.loglog(x[:n, 0]+1, np.exp(self.fit(x[:n, 0]+1, *p)), "k--", label=r"$x^{-%.2g}e^{-%.2g x}$" % (p[0], p[1]))
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Size by area")
        self.ax1.legend()

        self.ax2.loglog(x_nr+1, y_nr, label="hist")
        self.ax2.loglog(x_nr[:n_nr, 0]+1, np.exp(self.fit(x_nr[:n_nr, 0]+1, *p_nr)), "k--", label=r"$x^{-%.2g}e^{-%.2g x}$" % (p_nr[0], p_nr[1]))
        self.ax2.set_xlabel("Cluster size")
        self.ax2.set_ylabel("Frequency")
        self.ax2.set_title("Size by number")
        self.ax2.legend()

    def plot_particles_hex(self):
        # for normalization
        self.fig, (self.ax1, self.ax2) =  plt.subplots(nrows=1, ncols=2, figsize=(20, 20))
        y = self.data_hex[0::2]
        y_nr = self.data_hex_nr[0::2]
        sum_hist = np.sum(y)
        sum_hist_nr = np.sum(y_nr)
        x = np.argwhere(y != 0)
        x_nr = np.argwhere(y_nr != 0)

        y = y[x] / sum_hist
        y_nr = y_nr[x_nr] / sum_hist_nr

        x_int = x + 1

        p = opt.curve_fit(self.fit, x[:20, 0]+1, np.log(y[:20, 0]), [2.0, 5.0/max(x)], bounds=(0, np.inf))[0]
        p_nr = opt.curve_fit(self.fit, x_nr[:32, 0]+1, np.log(y_nr[:32, 0]), [2.0, 5.0/max(x_nr)], bounds=(0, np.inf))[0]

        self.ax1.loglog(x+1, y, "r-", label="hist")
        self.ax1.loglog(x[:20, 0]+1, np.exp(self.fit(x[:20, 0]+1, *p)), "k--", label=r"$x^{-%.2g}e^{-%.2g x}$" % (p[0], p[1]))
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Size by area")
        self.ax1.legend()

        self.ax2.loglog(x_nr+1, y_nr, label="hist")
        self.ax2.loglog(x_nr[:32, 0]+1, np.exp(self.fit(x_nr[:32, 0]+1, *p_nr)), "k--", label=r"$x^{-%.2g}e^{-%.2g x}$" % (p_nr[0], p_nr[1]))
        self.ax2.set_xlabel("Cluster size")
        self.ax2.set_ylabel("Frequency")
        self.ax2.set_title("Size by number")
        self.ax2.legend()
    

    def plot_combination(self):
        # hexagonal
        self.fig, (self.ax1, self.ax2) =  plt.subplots(nrows=1, ncols=2, figsize=(20, 20))
        y_hex = self.data_hex[0::2]
        y_hex_nr = self.data_hex_nr[0::2]
        sum_hist_hex = np.sum(y_hex)
        sum_hist_hex_nr = np.sum(y_hex_nr)
        x_hex = np.argwhere(y_hex != 0)
        x_hex_nr = np.argwhere(y_hex_nr != 0)

        y_hex = y_hex[x_hex] / sum_hist_hex
        y_hex_nr = y_hex_nr[x_hex_nr] / sum_hist_hex_nr

        #triangular
        y_tri = self.data_tri[0::2]
        y_tri_nr = self.data_tri_nr[0::2]
        sum_hist_tri = np.sum(y_tri)
        sum_hist_tri_nr = np.sum(y_tri_nr)
        x_tri = np.argwhere(y_tri != 0)
        x_tri_nr = np.argwhere(y_tri_nr != 0)

        y_tri = y_tri[x_tri] / sum_hist_tri
        y_tri_nr = y_tri_nr[x_tri_nr] / sum_hist_tri_nr

        self.ax1.loglog(x_hex+1, y_hex, "g", label="hex")
        self.ax1.loglog(x_tri+1, y_tri, "r", label="tri")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Size by area")
        self.ax1.legend()

        self.ax2.loglog(x_hex_nr+1, y_hex_nr, "g", label="hex")
        self.ax2.loglog(x_tri_nr+1, y_tri_nr, "r", label="tri")
        self.ax2.set_xlabel("Cluster size")
        self.ax2.set_ylabel("Frequency")
        self.ax2.set_title("Size by number")
        self.ax2.legend()

    def plot_comparison_hex(self):
        self.fig, self.ax1 =  plt.subplots(nrows=1, ncols=1, figsize=(20, 20))
        y = self.data_hex[0::2]
        y_nr = self.data_hex_nr[0::2]
        sum_hist = np.sum(y)
        sum_hist_nr = np.sum(y_nr)
        x = np.argwhere(y != 0)
        x_nr = np.argwhere(y_nr != 0)
        n = np.intersect1d(x, x_nr)
        diff = np.abs(y[n]/sum_hist - y_nr[n]/sum_hist_nr)
        
        y = y[x] / sum_hist
        y_nr = y_nr[x_nr] / sum_hist_nr
        
        self.ax1.loglog(x+1, y, "g", label="# by area")
        self.ax1.loglog(x_nr+1, y_nr, "m", label="# by number")
        self.ax1.loglog(n, diff, "--b", label="difference")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Comparison of Hex. Lattice")
        self.ax1.legend()

    def plot_comparison_tri(self):
        self.fig, self.ax1 =  plt.subplots(nrows=1, ncols=1, figsize=(20, 20))
        y = self.data_tri[0::2]
        y_nr = self.data_tri_nr[0::2]
        sum_hist = np.sum(y)
        sum_hist_nr = np.sum(y_nr)
        x = np.argwhere(y != 0)
        x_nr = np.argwhere(y_nr != 0)
        n = np.intersect1d(x, x_nr)
        diff = np.abs(y[n]/sum_hist - y_nr[n]/sum_hist_nr)
        
        y = y[x] / sum_hist
        y_nr = y_nr[x_nr] / sum_hist_nr
        
        self.ax1.loglog(x+1, y, "g", label="# by area")
        self.ax1.loglog(x_nr+1, y_nr, "m", label="# by number")
        self.ax1.loglog(n+1, diff, "--b", label="difference")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Comparison of Tri. Lattice")
        self.ax1.legend()

     
    def plot_time_evolution(self):
        plt.rcParams.update({'font.size': 18})
        # ############################### SQUARE DATA #################################
        data_1_square = np.genfromtxt("./number/square_number_alpha0.100_phi0.47_L100.txt", delimiter=" ")
        data_2_square = np.genfromtxt("./number/square_number_alpha0.010_phi0.47_L100.txt", delimiter=" ")
        data_7_square = np.genfromtxt("./number/square_number_alpha0.001_phi0.47_L100.txt", delimiter=" ")
        data_10_square = np.genfromtxt("./number/square_number_alpha0.000_phi0.47_L100.txt", delimiter=" ")

        data_3_square = np.genfromtxt("./number/square_number_alpha0.100_phi0.55_L100_2.txt", delimiter=" ")
        data_4_square = np.genfromtxt("./number/square_number_alpha0.010_phi0.55_L100_2.txt", delimiter=" ")
        data_8_square = np.genfromtxt("./number/square_number_alpha0.001_phi0.55_L100_2.txt", delimiter=" ")
        data_11_square = np.genfromtxt("./number/square_number_alpha0.000_phi0.55_L100_2.txt", delimiter=" ")

        data_5_square = np.genfromtxt("./number/square_number_alpha0.100_phi0.58_L100_3.txt", delimiter=" ")
        data_6_square = np.genfromtxt("./number/square_number_alpha0.010_phi0.58_L100_3.txt", delimiter=" ")
        data_9_square = np.genfromtxt("./number/square_number_alpha0.001_phi0.58_L100_3.txt", delimiter=" ")
        data_12_square = np.genfromtxt("./number/square_number_alpha0.000_phi0.58_L100_3.txt", delimiter=" ")
        # ############################### TRI DATA #################################
        data_1_tri = np.genfromtxt("./number/tri_number_alpha0.100_phi0.40_L100.txt", delimiter=" ")
        data_2_tri = np.genfromtxt("./number/tri_number_alpha0.010_phi0.40_L100.txt", delimiter=" ")
        data_7_tri = np.genfromtxt("./number/tri_number_alpha0.001_phi0.40_L100.txt", delimiter=" ")
        data_10_tri = np.genfromtxt("./number/tri_number_alpha0.000_phi0.40_L100.txt", delimiter=" ")

        data_3_tri = np.genfromtxt("./number/tri_number_alpha0.100_phi0.45_L100_2.txt", delimiter=" ")
        data_4_tri = np.genfromtxt("./number/tri_number_alpha0.010_phi0.45_L100_2.txt", delimiter=" ")
        data_8_tri = np.genfromtxt("./number/tri_number_alpha0.001_phi0.45_L100_2.txt", delimiter=" ")
        data_11_tri = np.genfromtxt("./number/tri_number_alpha0.000_phi0.45_L100_2.txt", delimiter=" ")

        data_5_tri = np.genfromtxt("./number/tri_number_alpha0.100_phi0.47_L100_3.txt", delimiter=" ")
        data_6_tri = np.genfromtxt("./number/tri_number_alpha0.010_phi0.47_L100_3.txt", delimiter=" ")
        data_9_tri = np.genfromtxt("./number/tri_number_alpha0.001_phi0.47_L100_3.txt", delimiter=" ")
        data_12_tri = np.genfromtxt("./number/tri_number_alpha0.000_phi0.47_L100_3.txt", delimiter=" ")
        # ############################### HEX DATA #################################
        data_1_hex = np.genfromtxt("./number/hex_number_alpha0.100_phi1.10_L100.txt", delimiter=" ")
        data_2_hex = np.genfromtxt("./number/hex_number_alpha0.010_phi1.10_L100.txt", delimiter=" ")
        data_7_hex = np.genfromtxt("./number/hex_number_alpha0.001_phi1.10_L100.txt", delimiter=" ")
        data_10_hex = np.genfromtxt("./number/hex_number_alpha0.000_phi1.10_L100.txt", delimiter=" ")

        data_3_hex = np.genfromtxt("./number/hex_number_alpha0.100_phi1.32_L100_2.txt", delimiter=" ")
        data_4_hex = np.genfromtxt("./number/hex_number_alpha0.010_phi1.32_L100_2.txt", delimiter=" ")
        data_8_hex = np.genfromtxt("./number/hex_number_alpha0.001_phi1.32_L100_2.txt", delimiter=" ")
        data_11_hex = np.genfromtxt("./number/hex_number_alpha0.000_phi1.32_L100_2.txt", delimiter=" ")

        data_5_hex = np.genfromtxt("./number/hex_number_alpha0.100_phi1.41_L100_3.txt", delimiter=" ")
        data_6_hex = np.genfromtxt("./number/hex_number_alpha0.010_phi1.41_L100_3.txt", delimiter=" ")
        data_9_hex = np.genfromtxt("./number/hex_number_alpha0.001_phi1.41_L100_3.txt", delimiter=" ")
        data_12_hex = np.genfromtxt("./number/hex_number_alpha0.000_phi1.41_L100_3.txt", delimiter=" ")

        # Square
        if True:
            # ###################################### FIG 1 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_1_square[1:, 0], data_1_square[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_2_square[1:, 0], data_2_square[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_7_square[1:, 0], data_7_square[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_10_square[1:, 0], data_10_square[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_1_square[1:, 0], data_1_square[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_2_square[1:, 0], data_2_square[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_7_square[1:, 0], data_7_square[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_10_square[1:, 0], data_10_square[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/square_1.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 2 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_3_square[1:, 0], data_3_square[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_4_square[1:, 0], data_4_square[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_8_square[1:, 0], data_8_square[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_11_square[1:, 0], data_11_square[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_3_square[1:, 0], data_3_square[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_4_square[1:, 0], data_4_square[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_8_square[1:, 0], data_8_square[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_11_square[1:, 0], data_11_square[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/square_2.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 3 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_5_square[1:, 0], data_5_square[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_6_square[1:, 0], data_6_square[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_9_square[1:, 0], data_9_square[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_12_square[1:, 0], data_12_square[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_5_square[1:, 0], data_5_square[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_6_square[1:, 0], data_6_square[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_9_square[1:, 0], data_9_square[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_12_square[1:, 0], data_12_square[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/square_3.pdf", dpi=200, bbox_inches='tight')

            # comparison of same alpha different n
            if True:
                # ###################################### FIG 1 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_1_square[1:, 0], data_1_square[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_3_square[1:, 0], data_3_square[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_5_square[1:, 0], data_5_square[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_1_square[1:, 0], data_1_square[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_3_square[1:, 0], data_3_square[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_5_square[1:, 0], data_5_square[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/square_comp_n_alpha_0.1.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 2 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_2_square[1:, 0], data_2_square[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_4_square[1:, 0], data_4_square[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_6_square[1:, 0], data_6_square[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_2_square[1:, 0], data_2_square[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_4_square[1:, 0], data_4_square[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_6_square[1:, 0], data_6_square[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/square_comp_n_alpha_0.01.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 3 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_7_square[1:, 0], data_7_square[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_8_square[1:, 0], data_8_square[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_9_square[1:, 0], data_9_square[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_7_square[1:, 0], data_7_square[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_8_square[1:, 0], data_8_square[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_9_square[1:, 0], data_9_square[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/square_comp_n_alpha_0.001.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 4 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_10_square[1:, 0], data_10_square[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_11_square[1:, 0], data_11_square[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_12_square[1:, 0], data_12_square[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_10_square[1:, 0], data_10_square[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_11_square[1:, 0], data_11_square[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_12_square[1:, 0], data_12_square[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/square_comp_n_alpha_0.pdf", dpi=200, bbox_inches='tight')

        # Tri
        if False:
            # ###################################### FIG 1 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_1_tri[1:, 0], data_1_tri[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_2_tri[1:, 0], data_2_tri[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_7_tri[1:, 0], data_7_tri[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_10_tri[1:, 0], data_10_tri[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_1_tri[1:, 0], data_1_tri[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_2_tri[1:, 0], data_2_tri[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_7_tri[1:, 0], data_7_tri[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_10_tri[1:, 0], data_10_tri[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/tri_1.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 2 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_3_tri[1:, 0], data_3_tri[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_4_tri[1:, 0], data_4_tri[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_8_tri[1:, 0], data_8_tri[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_11_tri[1:, 0], data_11_tri[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_3_tri[1:, 0], data_3_tri[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_4_tri[1:, 0], data_4_tri[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_8_tri[1:, 0], data_8_tri[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_11_tri[1:, 0], data_11_tri[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/tri_2.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 3 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_5_tri[1:, 0], data_5_tri[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_6_tri[1:, 0], data_6_tri[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_9_tri[1:, 0], data_9_tri[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_12_tri[1:, 0], data_12_tri[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_5_tri[1:, 0], data_5_tri[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_6_tri[1:, 0], data_6_tri[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_9_tri[1:, 0], data_9_tri[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_12_tri[1:, 0], data_12_tri[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/tri_3.pdf", dpi=200, bbox_inches='tight')

            # comparison of same alpha different n
            if True:
                # ###################################### FIG 1 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_1_tri[1:, 0], data_1_tri[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_3_tri[1:, 0], data_3_tri[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_5_tri[1:, 0], data_5_tri[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_1_tri[1:, 0], data_1_tri[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_3_tri[1:, 0], data_3_tri[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_5_tri[1:, 0], data_5_tri[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/tri_comp_n_alpha_0.1.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 2 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_2_tri[1:, 0], data_2_tri[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_4_tri[1:, 0], data_4_tri[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_6_tri[1:, 0], data_6_tri[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_2_tri[1:, 0], data_2_tri[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_4_tri[1:, 0], data_4_tri[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_6_tri[1:, 0], data_6_tri[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/tri_comp_n_alpha_0.01.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 3 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_7_tri[1:, 0], data_7_tri[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_8_tri[1:, 0], data_8_tri[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_9_tri[1:, 0], data_9_tri[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_7_tri[1:, 0], data_7_tri[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_8_tri[1:, 0], data_8_tri[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_9_tri[1:, 0], data_9_tri[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/tri_comp_n_alpha_0.001.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 4 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_10_tri[1:, 0], data_10_tri[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_11_tri[1:, 0], data_11_tri[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_12_tri[1:, 0], data_12_tri[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_10_tri[1:, 0], data_10_tri[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_11_tri[1:, 0], data_11_tri[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_12_tri[1:, 0], data_12_tri[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/tri_comp_n_alpha_0.pdf", dpi=200, bbox_inches='tight')

        # Hex
        if False:
            # ###################################### FIG 1 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_1_hex[1:, 0], data_1_hex[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_2_hex[1:, 0], data_2_hex[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_7_hex[1:, 0], data_7_hex[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_10_hex[1:, 0], data_10_hex[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_1_hex[1:, 0], data_1_hex[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_2_hex[1:, 0], data_2_hex[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_7_hex[1:, 0], data_7_hex[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_10_hex[1:, 0], data_10_hex[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/hex_1.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 2 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_3_hex[1:, 0], data_3_hex[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_4_hex[1:, 0], data_4_hex[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_8_hex[1:, 0], data_8_hex[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_11_hex[1:, 0], data_11_hex[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")


            ax[1].plot(data_3_hex[1:, 0], data_3_hex[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_4_hex[1:, 0], data_4_hex[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_8_hex[1:, 0], data_8_hex[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_11_hex[1:, 0], data_11_hex[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/hex_2.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 3 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_5_hex[1:, 0], data_5_hex[1:, 1], label=r"$\alpha=0.1$")
            ax[0].plot(data_6_hex[1:, 0], data_6_hex[1:, 1], label=r"$\alpha=0.01$")
            ax[0].plot(data_9_hex[1:, 0], data_9_hex[1:, 1], label=r"$\alpha=0.001$")
            ax[0].plot(data_12_hex[1:, 0], data_12_hex[1:, 1], label=r"$\alpha=0$")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_5_hex[1:, 0], data_5_hex[1:, 2], label=r"$\alpha=0.1$")
            ax[1].plot(data_6_hex[1:, 0], data_6_hex[1:, 2], label=r"$\alpha=0.01$")
            ax[1].plot(data_9_hex[1:, 0], data_9_hex[1:, 2], label=r"$\alpha=0.001$")
            ax[1].plot(data_12_hex[1:, 0], data_12_hex[1:, 2], label=r"$\alpha=0$")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/hex_3.pdf", dpi=200, bbox_inches='tight')

            # n comparison
            if True:
                # ###################################### FIG 1 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_1_hex[1:, 0], data_1_hex[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_3_hex[1:, 0], data_3_hex[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_5_hex[1:, 0], data_5_hex[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_1_hex[1:, 0], data_1_hex[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_3_hex[1:, 0], data_3_hex[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_5_hex[1:, 0], data_5_hex[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/hex_comp_n_alpha_0.1.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 2 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_2_hex[1:, 0], data_2_hex[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_4_hex[1:, 0], data_4_hex[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_6_hex[1:, 0], data_6_hex[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_2_hex[1:, 0], data_2_hex[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_4_hex[1:, 0], data_4_hex[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_6_hex[1:, 0], data_6_hex[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/hex_comp_n_alpha_0.01.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 3 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_7_hex[1:, 0], data_7_hex[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_8_hex[1:, 0], data_8_hex[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_9_hex[1:, 0], data_9_hex[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_7_hex[1:, 0], data_7_hex[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_8_hex[1:, 0], data_8_hex[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_9_hex[1:, 0], data_9_hex[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/hex_comp_n_alpha_0.001.pdf", dpi=200, bbox_inches='tight')

                # ###################################### FIG 4 ################################################
                fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
                plt.tight_layout()

                ax[0].plot(data_10_hex[1:, 0], data_10_hex[1:, 1], label=r"$n_\mathrm{max}=1$")
                ax[0].plot(data_11_hex[1:, 0], data_11_hex[1:, 1], label=r"$n_\mathrm{max}=2$")
                ax[0].plot(data_12_hex[1:, 0], data_12_hex[1:, 1], label=r"$n_\mathrm{max}=3$")
                ax[0].grid()
                ax[0].set_yscale("log")
                ax[0].legend()
                ax[0].set_ylabel(r"$w_N$")
                ax[0].set_xscale("log")

                ax[1].plot(data_10_hex[1:, 0], data_10_hex[1:, 2], label=r"$n_\mathrm{max}=1$")
                ax[1].plot(data_11_hex[1:, 0], data_11_hex[1:, 2], label=r"$n_\mathrm{max}=2$")
                ax[1].plot(data_12_hex[1:, 0], data_12_hex[1:, 2], label=r"$n_\mathrm{max}=3$")
                ax[1].grid()
                ax[1].set_yscale("log")
                ax[1].set_xscale("log")
                ax[1].legend()
                ax[1].set_ylabel("Mean cluster size")
                ax[1].set_xlabel(r"Time $t$")

                plt.savefig("./plots/hex_comp_n_alpha_0.pdf", dpi=200, bbox_inches='tight')

        # comparison of lattices given n and alpha
        if True:
            # ###################################### FIG 1 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_1_square[1:, 0], data_1_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_1_tri[1:, 0], data_1_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_1_hex[1:, 0], data_1_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_1_square[1:, 0], data_1_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_1_tri[1:, 0], data_1_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_1_hex[1:, 0], data_1_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_1_alpha_0.1.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 10 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_10_square[1:, 0], data_10_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_10_tri[1:, 0], data_10_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_10_hex[1:, 0], data_10_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_10_square[1:, 0], data_10_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_10_tri[1:, 0], data_10_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_10_hex[1:, 0], data_10_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_1_alpha_0.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 11 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_11_square[1:, 0], data_11_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_11_tri[1:, 0], data_11_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_11_hex[1:, 0], data_11_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_11_square[1:, 0], data_11_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_11_tri[1:, 0], data_11_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_11_hex[1:, 0], data_11_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_2_alpha_0.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 12 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_12_square[1:, 0], data_12_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_12_tri[1:, 0], data_12_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_12_hex[1:, 0], data_12_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_12_square[1:, 0], data_12_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_12_tri[1:, 0], data_12_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_12_hex[1:, 0], data_12_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_3_alpha_0.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 2 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_2_square[1:, 0], data_2_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_2_tri[1:, 0], data_2_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_2_hex[1:, 0], data_2_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_2_square[1:, 0], data_2_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_2_tri[1:, 0], data_2_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_2_hex[1:, 0], data_2_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_1_alpha_0.01.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 3 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_7_square[1:, 0], data_7_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_7_tri[1:, 0], data_7_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_7_hex[1:, 0], data_7_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_7_square[1:, 0], data_7_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_7_tri[1:, 0], data_7_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_7_hex[1:, 0], data_7_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_1_alpha_0.001.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 4 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_3_square[1:, 0], data_3_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_3_tri[1:, 0], data_3_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_3_hex[1:, 0], data_3_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_3_square[1:, 0], data_3_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_3_tri[1:, 0], data_3_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_3_hex[1:, 0], data_3_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_2_alpha_0.1.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 5 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_4_square[1:, 0], data_4_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_4_tri[1:, 0], data_4_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_4_hex[1:, 0], data_4_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_4_square[1:, 0], data_4_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_4_tri[1:, 0], data_4_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_4_hex[1:, 0], data_4_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_2_alpha_0.01.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 6 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_8_square[1:, 0], data_8_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_8_tri[1:, 0], data_8_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_8_hex[1:, 0], data_8_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_8_square[1:, 0], data_8_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_8_tri[1:, 0], data_8_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_8_hex[1:, 0], data_8_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_2_alpha_0.001.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 7 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_5_square[1:, 0], data_5_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_5_tri[1:, 0], data_5_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_5_hex[1:, 0], data_5_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_5_square[1:, 0], data_5_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_5_tri[1:, 0], data_5_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_5_hex[1:, 0], data_5_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_3_alpha_0.1.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 8 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_6_square[1:, 0], data_6_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_6_tri[1:, 0], data_6_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_6_hex[1:, 0], data_6_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_6_square[1:, 0], data_6_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_6_tri[1:, 0], data_6_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_6_hex[1:, 0], data_6_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_3_alpha_0.01.pdf", dpi=200, bbox_inches='tight')

            # ###################################### FIG 1 ################################################

            fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(12, 8), sharex=True)
            plt.tight_layout()

            ax[0].plot(data_9_square[1:, 0], data_9_square[1:, 1], label=r"Sq.Lat.")
            ax[0].plot(data_9_tri[1:, 0], data_9_tri[1:, 1], label=r"Tr.Lat.")
            ax[0].plot(data_9_hex[1:, 0], data_9_hex[1:, 1], label=r"Hx.Lat.")
            ax[0].grid()
            ax[0].set_yscale("log")
            ax[0].legend()
            ax[0].set_ylabel(r"$w_N$")
            ax[0].set_xscale("log")

            ax[1].plot(data_9_square[1:, 0], data_9_square[1:, 2], label=r"Sq.Lat.")
            ax[1].plot(data_9_tri[1:, 0], data_9_tri[1:, 2], label=r"Tr.Lat.")
            ax[1].plot(data_9_hex[1:, 0], data_9_hex[1:, 2], label=r"Hx.Lat.")
            ax[1].grid()
            ax[1].set_yscale("log")
            ax[1].set_xscale("log")
            ax[1].legend()
            ax[1].set_ylabel("Mean cluster size")
            ax[1].set_xlabel(r"Time $t$")

            plt.savefig("./plots/lat_comp_n_3_alpha_0.001.pdf", dpi=200, bbox_inches='tight')

 
if __name__ == "__main__":
    dist = distribution()

    #dist.plot_lattice_type()
    dist.plot_time_evolution()

    #dist.plot_heat()
    #dist.plot_alpha()
    #dist.plot_cum_hex()
    #plt.show()
    #dist.plot_square_L()
    #dist.plot_hex_L()
    #dist.plot_tri_L()
    #dist.plot_cum_tri()
    #plt.show()

    #dist.plot_particles_tri()
    #plt.show()
    #dist.plot_particles_hex()
    #plt.show()
    #dist.plot_combination()
    #plt.show()
    #dist.plot_comparison_hex()
    #plt.show()
    #dist.plot_comparison_tri()
    #plt.show()

        
