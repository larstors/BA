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


if __name__ == "__main__":
    dist = distribution()
    dist.plot_particles_tri()
    plt.show()
    dist.plot_particles_hex()
    plt.show()
    dist.plot_combination()
    plt.show()
    dist.plot_comparison_hex()
    plt.show()
    dist.plot_comparison_tri()
    plt.show()

        
