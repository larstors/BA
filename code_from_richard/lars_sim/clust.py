import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import time
import random as rn



class distribution:
    def __init__(self):
        self.data_tri = np.genfromtxt("tridist.txt", delimiter=" ")
        self.data_tri_nr = np.genfromtxt("tridist_nr.txt", delimiter=" ")
        self.data_hex = np.genfromtxt("hexdist.txt", delimiter=" ")
        self.data_hex_nr = np.genfromtxt("hexdist_nr.txt", delimiter=" ")

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

        self.ax1.loglog(x, y)
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Size by area")

        self.ax2.loglog(x_nr, y_nr)
        self.ax2.set_xlabel("Cluster size")
        self.ax2.set_ylabel("Frequency")
        self.ax2.set_title("Size by number")

        

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

        self.ax1.loglog(x, y)
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Size by area")

        self.ax2.loglog(x_nr, y_nr)
        self.ax2.set_xlabel("Cluster size")
        self.ax2.set_ylabel("Frequency")
        self.ax2.set_title("Size by number")
    

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

        self.ax1.loglog(x_hex, y_hex, "g", label="hex")
        self.ax1.loglog(x_tri, y_tri, "r", label="tri")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Size by area")
        self.ax1.legend()

        self.ax2.loglog(x_hex_nr, y_hex_nr, "g", label="hex")
        self.ax2.loglog(x_tri_nr, y_tri_nr, "r", label="tri")
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
        
        self.ax1.loglog(x, y, "g", label="# by area")
        self.ax1.loglog(x_nr, y_nr, "m", label="# by number")
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
        
        self.ax1.loglog(x, y, "g", label="# by area")
        self.ax1.loglog(x_nr, y_nr, "m", label="# by number")
        self.ax1.loglog(n, diff, "--b", label="difference")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_title("Comparison of Tri. Lattice")
        self.ax1.legend()


if __name__ == "__main__":
    dist = distribution()
    dist.plot_combination()
    plt.show()
    dist.plot_comparison_hex()
    plt.show()
    dist.plot_comparison_tri()
    plt.show()

        
