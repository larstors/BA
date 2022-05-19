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

    def plot_cum_hex(self):
        self.fig, (self.ax1, self.ax2) =  plt.subplots(nrows=1, ncols=2, figsize=(20, 20))
        y = self.data_hex[0::2]
        y_nr = self.data_hex_nr[0::2]
        sum_hist = np.sum(y)
        sum_hist_nr = np.sum(y_nr)
        x = np.argwhere(y != 0)
        x_nr = np.argwhere(y_nr != 0)
        
        
        y = y[x] / sum_hist
        y_nr = y_nr[x_nr] / sum_hist_nr


        y_c = np.cumsum(y)
        y_c_nr = np.cumsum(y_nr)
        

        self.ax1.plot(x+1, y_c, "g", label="# by area")
        self.ax1.plot(x_nr+1, y_c_nr, "m", label="# by number")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_xscale("log")
        self.ax1.set_title("CDF of Hex. Lattice cluster dist")
        self.ax1.legend()

        self.ax2.plot(x+1, y, "g", label="# by area")
        self.ax2.plot(x_nr+1, y_nr, "m", label="# by number")
        self.ax2.set_xlabel("Cluster size")
        self.ax2.set_ylabel("Frequency")
        self.ax2.set_xscale("log")
        self.ax2.set_title("PDF of Hex. Lattice cluster dist")
        self.ax2.legend()
    
    def plot_cum_tri(self):
        self.fig, self.ax1 =  plt.subplots(nrows=1, ncols=1, figsize=(20, 20))
        y = self.data_tri[0::2]
        y_nr = self.data_tri_nr[0::2]
        sum_hist = np.sum(y)
        sum_hist_nr = np.sum(y_nr)
        x = np.argwhere(y != 0)
        x_nr = np.argwhere(y_nr != 0)
        
        y = y[x] / sum_hist
        y_nr = y_nr[x_nr] / sum_hist_nr

        y_c = np.cumsum(y)
        y_c_nr = np.cumsum(y_nr)
        
        self.ax1.plot(x+1, y_c, "g", label="# by area")
        self.ax1.plot(x_nr+1, y_c_nr, "m", label="# by number")
        self.ax1.set_xlabel("Cluster size")
        self.ax1.set_ylabel("Frequency")
        self.ax1.set_xscale("log")
        self.ax1.set_title("CDF of Tri. Lattice cluster dist")
        self.ax1.legend()
    
    def plot_square_L(self):
        """Plots cumulative distribution for differing lattice sizes
        """
        data_L_10 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.10_L10.txt", delimiter=" ")[0::2]
        data_L_20 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.10_L20.txt", delimiter=" ")[0::2]
        data_L_50 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.10_L50.txt", delimiter=" ")[0::2]
        data_L_100 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.10_L100.txt", delimiter=" ")[0::2]
        data_L_200 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.10_L200.txt", delimiter=" ")[0::2]
        data_L_500 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.10_L500.txt", delimiter=" ")[0::2]
        
        """
        data_L_nr_10 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.10_L10.txt", delimiter=" ")[0::2]
        data_L_nr_20 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.10_L20.txt", delimiter=" ")[0::2]
        data_L_nr_50 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.10_L50.txt", delimiter=" ")[0::2]
        data_L_nr_100 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.10_L100.txt", delimiter=" ")[0::2]
        data_L_nr_200 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.10_L200.txt", delimiter=" ")[0::2]
        data_L_nr_500 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.10_L500.txt", delimiter=" ")[0::2]
        """


        x_10 = np.argwhere(data_L_10 != 0)
        x_20 = np.argwhere(data_L_20 != 0)
        x_50 = np.argwhere(data_L_50 != 0)
        x_100 = np.argwhere(data_L_100 != 0)
        x_200 = np.argwhere(data_L_200 != 0)
        x_500 = np.argwhere(data_L_500 != 0)
        """
        x_nr_10 = np.argwhere(data_L_nr_10 != 0)
        x_nr_20 = np.argwhere(data_L_nr_20 != 0)
        x_nr_50 = np.argwhere(data_L_nr_50 != 0)
        x_nr_100 = np.argwhere(data_L_nr_100 != 0)
        x_nr_200 = np.argwhere(data_L_nr_200 != 0)
        x_nr_500 = np.argwhere(data_L_nr_500 != 0)
        """
        data_L_10 = data_L_10[x_10]/np.sum(data_L_10[x_10])         
        data_L_20 = data_L_20[x_20]/np.sum(data_L_20[x_20])
        data_L_50 = data_L_50[x_50]/np.sum(data_L_50[x_50])
        data_L_100 = data_L_100[x_100]/np.sum(data_L_100[x_100])
        data_L_200 = data_L_200[x_200]/np.sum(data_L_200[x_200])
        data_L_500 = data_L_500[x_500]/np.sum(data_L_500[x_500])
        """
        data_L_nr_10 = data_L_nr_10[x_nr_10]/np.sum(data_L_nr_10[x_nr_10])         
        data_L_nr_20 = data_L_nr_20[x_nr_20]/np.sum(data_L_nr_20[x_nr_20])
        data_L_nr_50 = data_L_nr_50[x_nr_50]/np.sum(data_L_nr_50[x_nr_50])
        data_L_nr_100 = data_L_nr_100[x_nr_100]/np.sum(data_L_nr_100[x_nr_100])
        data_L_nr_200 = data_L_nr_200[x_nr_200]/np.sum(data_L_nr_200[x_nr_200])
        data_L_nr_500 = data_L_nr_500[x_nr_500]/np.sum(data_L_nr_500[x_nr_500])
        """

        fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
        
        plt.rcParams.update({'font.size': 18})
        plt.title("Square Lattice")
        plt.loglog(x_10 + 1, np.cumsum(data_L_10), "r-", label="L=10")
        plt.loglog(x_20 + 1, np.cumsum(data_L_20), "g-", label="L=20")
        plt.loglog(x_50 + 1, np.cumsum(data_L_50), "b-", label="L=50")
        plt.loglog(x_100[::5] + 1, np.cumsum(data_L_100)[::5], "m-x", label="L=100")
        plt.loglog(x_200[::5] + 1, np.cumsum(data_L_200)[::5], "k-o", label="L=200")
        plt.loglog(x_500[::5] + 1, np.cumsum(data_L_500)[::5], "y-s", label="L=500")
        plt.xlabel(r"Cluster size $k$")
        plt.ylabel(r"CDF $C(k)$")
        plt.legend()
        plt.grid()
        plt.savefig("./plots/square_L.pdf", dpi=200)
        
        # ###############################################################################################
        # the following is for higher densities
        data_L_10 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.50_L10.txt", delimiter=" ")[0::2]
        data_L_20 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.50_L20.txt", delimiter=" ")[0::2]
        data_L_50 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.50_L50.txt", delimiter=" ")[0::2]
        data_L_100 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.50_L100.txt", delimiter=" ")[0::2]
        data_L_200 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.50_L200.txt", delimiter=" ")[0::2]
        data_L_500 = np.genfromtxt("./latticesize/square_alpha0.01_phi0.50_L500.txt", delimiter=" ")[0::2]
        
        """
        data_L_nr_10 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.50_L10.txt", delimiter=" ")[0::2]
        data_L_nr_20 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.50_L20.txt", delimiter=" ")[0::2]
        data_L_nr_50 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.50_L50.txt", delimiter=" ")[0::2]
        data_L_nr_100 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.50_L100.txt", delimiter=" ")[0::2]
        data_L_nr_200 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.50_L200.txt", delimiter=" ")[0::2]
        data_L_nr_500 = np.genfromtxt("./latticesize/square_nr_alpha0.01_phi0.50_L500.txt", delimiter=" ")[0::2]
        """


        x_10 = np.argwhere(data_L_10 != 0)
        x_20 = np.argwhere(data_L_20 != 0)
        x_50 = np.argwhere(data_L_50 != 0)
        x_100 = np.argwhere(data_L_100 != 0)
        x_200 = np.argwhere(data_L_200 != 0)
        x_500 = np.argwhere(data_L_500 != 0)
        """
        x_nr_10 = np.argwhere(data_L_nr_10 != 0)
        x_nr_20 = np.argwhere(data_L_nr_20 != 0)
        x_nr_50 = np.argwhere(data_L_nr_50 != 0)
        x_nr_100 = np.argwhere(data_L_nr_100 != 0)
        x_nr_200 = np.argwhere(data_L_nr_200 != 0)
        x_nr_500 = np.argwhere(data_L_nr_500 != 0)
        """
        data_L_10 = data_L_10[x_10]/np.sum(data_L_10[x_10])         
        data_L_20 = data_L_20[x_20]/np.sum(data_L_20[x_20])
        data_L_50 = data_L_50[x_50]/np.sum(data_L_50[x_50])
        data_L_100 = data_L_100[x_100]/np.sum(data_L_100[x_100])
        data_L_200 = data_L_200[x_200]/np.sum(data_L_200[x_200])
        data_L_500 = data_L_500[x_500]/np.sum(data_L_500[x_500])
        """
        data_L_nr_10 = data_L_nr_10[x_nr_10]/np.sum(data_L_nr_10[x_nr_10])         
        data_L_nr_20 = data_L_nr_20[x_nr_20]/np.sum(data_L_nr_20[x_nr_20])
        data_L_nr_50 = data_L_nr_50[x_nr_50]/np.sum(data_L_nr_50[x_nr_50])
        data_L_nr_100 = data_L_nr_100[x_nr_100]/np.sum(data_L_nr_100[x_nr_100])
        data_L_nr_200 = data_L_nr_200[x_nr_200]/np.sum(data_L_nr_200[x_nr_200])
        data_L_nr_500 = data_L_nr_500[x_nr_500]/np.sum(data_L_nr_500[x_nr_500])
        """

        fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
        
        plt.rcParams.update({'font.size': 18})
        plt.title("Square Lattice")
        plt.loglog(x_10 + 1, np.cumsum(data_L_10), "r-", label="L=10")
        plt.loglog(x_20 + 1, np.cumsum(data_L_20), "g-", label="L=20")
        plt.loglog(x_50 + 1, np.cumsum(data_L_50), "b-", label="L=50")
        plt.loglog(x_100[::5] + 1, np.cumsum(data_L_100)[::5], "m-x", label="L=100")
        plt.loglog(x_200[::5] + 1, np.cumsum(data_L_200)[::5], "k-o", label="L=200")
        plt.loglog(x_500[::5] + 1, np.cumsum(data_L_500)[::5], "y-s", label="L=500")
        plt.xlabel(r"Cluster size $k$")
        plt.ylabel(r"CDF $C(k)$")
        plt.legend()
        plt.grid()
        plt.savefig("./plots/square_L_phi_0.5.pdf", dpi=200)
        

    def plot_hex_L(self):
        """Plots cumulative distribution for differing lattice sizes
        """
        data_L_10 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi0.20_L10.txt", delimiter=" ")[0::2]
        data_L_20 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi0.20_L20.txt", delimiter=" ")[0::2]
        data_L_50 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi0.20_L50.txt", delimiter=" ")[0::2]
        data_L_100 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi0.20_L100.txt", delimiter=" ")[0::2]
        data_L_200 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi0.20_L200.txt", delimiter=" ")[0::2]
        data_L_500 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi0.20_L500.txt", delimiter=" ")[0::2]
        
        """
        data_L_nr_10 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi0.20_L10.txt", delimiter=" ")[0::2]
        data_L_nr_20 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi0.20_L20.txt", delimiter=" ")[0::2]
        data_L_nr_50 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi0.20_L50.txt", delimiter=" ")[0::2]
        data_L_nr_100 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi0.20_L100.txt", delimiter=" ")[0::2]
        data_L_nr_200 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi0.20_L200.txt", delimiter=" ")[0::2]
        data_L_nr_500 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi0.20_L500.txt", delimiter=" ")[0::2]
        """

        x_10 = np.argwhere(data_L_10 != 0)
        x_20 = np.argwhere(data_L_20 != 0)
        x_50 = np.argwhere(data_L_50 != 0)
        x_100 = np.argwhere(data_L_100 != 0)
        x_200 = np.argwhere(data_L_200 != 0)
        x_500 = np.argwhere(data_L_500 != 0)

        """
        x_nr_10 = np.argwhere(data_L_nr_10 != 0)
        x_nr_20 = np.argwhere(data_L_nr_20 != 0)
        x_nr_50 = np.argwhere(data_L_nr_50 != 0)
        x_nr_100 = np.argwhere(data_L_nr_100 != 0)
        x_nr_200 = np.argwhere(data_L_nr_200 != 0)
        x_nr_500 = np.argwhere(data_L_nr_500 != 0)
        """

        data_L_10 = data_L_10[x_10]/np.sum(data_L_10[x_10])         
        data_L_20 = data_L_20[x_20]/np.sum(data_L_20[x_20])
        data_L_50 = data_L_50[x_50]/np.sum(data_L_50[x_50])
        data_L_100 = data_L_100[x_100]/np.sum(data_L_100[x_100])
        data_L_200 = data_L_200[x_200]/np.sum(data_L_200[x_200])
        data_L_500 = data_L_500[x_500]/np.sum(data_L_500[x_500])

        """
        data_L_nr_10 = data_L_nr_10[x_nr_10]/np.sum(data_L_nr_10[x_nr_10])         
        data_L_nr_20 = data_L_nr_20[x_nr_20]/np.sum(data_L_nr_20[x_nr_20])
        data_L_nr_50 = data_L_nr_50[x_nr_50]/np.sum(data_L_nr_50[x_nr_50])
        data_L_nr_100 = data_L_nr_100[x_nr_100]/np.sum(data_L_nr_100[x_nr_100])
        data_L_nr_200 = data_L_nr_200[x_nr_200]/np.sum(data_L_nr_200[x_nr_200])
        data_L_nr_500 = data_L_nr_500[x_nr_500]/np.sum(data_L_nr_500[x_nr_500])
        """

        fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
        plt.rcParams.update({'font.size': 18})
        plt.title("Hexagonal Lattice")
        plt.loglog(x_10 + 1, np.cumsum(data_L_10), "r-", label="L=10")
        plt.loglog(x_20 + 1, np.cumsum(data_L_20), "g-", label="L=20")
        plt.loglog(x_50 + 1, np.cumsum(data_L_50), "b-", label="L=50")
        plt.loglog(x_100[::5] + 1, np.cumsum(data_L_100)[::5], "m-x", label="L=100")
        plt.loglog(x_200[::5] + 1, np.cumsum(data_L_200)[::5], "k-o", label="L=200")
        plt.loglog(x_500[::5] + 1, np.cumsum(data_L_500)[::5], "y-s", label="L=500")
        plt.xlabel(r"Cluster size $k$")
        plt.ylabel(r"CDF $C(k)$")
        plt.legend()
        plt.grid()
        plt.savefig("./plots/hexagonal_L.pdf", dpi=200)

        # ###################################################################################
        # The following is for higher densities
        data_L_10 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi1.00_L10.txt", delimiter=" ")[0::2]
        data_L_20 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi1.00_L20.txt", delimiter=" ")[0::2]
        data_L_50 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi1.00_L50.txt", delimiter=" ")[0::2]
        data_L_100 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi1.00_L100.txt", delimiter=" ")[0::2]
        data_L_200 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi1.00_L200.txt", delimiter=" ")[0::2]
        data_L_500 = np.genfromtxt("./latticesize/hexagonal_alpha0.01_phi1.00_L500.txt", delimiter=" ")[0::2]
        
        """
        data_L_nr_10 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi1.00_L10.txt", delimiter=" ")[0::2]
        data_L_nr_20 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi1.00_L20.txt", delimiter=" ")[0::2]
        data_L_nr_50 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi1.00_L50.txt", delimiter=" ")[0::2]
        data_L_nr_100 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi1.00_L100.txt", delimiter=" ")[0::2]
        data_L_nr_200 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi1.00_L200.txt", delimiter=" ")[0::2]
        data_L_nr_500 = np.genfromtxt("./latticesize/hexagonal_nr_alpha0.01_phi1.00_L500.txt", delimiter=" ")[0::2]
        """

        x_10 = np.argwhere(data_L_10 != 0)
        x_20 = np.argwhere(data_L_20 != 0)
        x_50 = np.argwhere(data_L_50 != 0)
        x_100 = np.argwhere(data_L_100 != 0)
        x_200 = np.argwhere(data_L_200 != 0)
        x_500 = np.argwhere(data_L_500 != 0)

        """
        x_nr_10 = np.argwhere(data_L_nr_10 != 0)
        x_nr_20 = np.argwhere(data_L_nr_20 != 0)
        x_nr_50 = np.argwhere(data_L_nr_50 != 0)
        x_nr_100 = np.argwhere(data_L_nr_100 != 0)
        x_nr_200 = np.argwhere(data_L_nr_200 != 0)
        x_nr_500 = np.argwhere(data_L_nr_500 != 0)
        """

        data_L_10 = data_L_10[x_10]/np.sum(data_L_10[x_10])         
        data_L_20 = data_L_20[x_20]/np.sum(data_L_20[x_20])
        data_L_50 = data_L_50[x_50]/np.sum(data_L_50[x_50])
        data_L_100 = data_L_100[x_100]/np.sum(data_L_100[x_100])
        data_L_200 = data_L_200[x_200]/np.sum(data_L_200[x_200])
        data_L_500 = data_L_500[x_500]/np.sum(data_L_500[x_500])

        """
        data_L_nr_10 = data_L_nr_10[x_nr_10]/np.sum(data_L_nr_10[x_nr_10])         
        data_L_nr_20 = data_L_nr_20[x_nr_20]/np.sum(data_L_nr_20[x_nr_20])
        data_L_nr_50 = data_L_nr_50[x_nr_50]/np.sum(data_L_nr_50[x_nr_50])
        data_L_nr_100 = data_L_nr_100[x_nr_100]/np.sum(data_L_nr_100[x_nr_100])
        data_L_nr_200 = data_L_nr_200[x_nr_200]/np.sum(data_L_nr_200[x_nr_200])
        data_L_nr_500 = data_L_nr_500[x_nr_500]/np.sum(data_L_nr_500[x_nr_500])
        """

        fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
        plt.rcParams.update({'font.size': 18})
        plt.title("Hexagonal Lattice")
        plt.loglog(x_10 + 1, np.cumsum(data_L_10), "r-", label="L=10")
        plt.loglog(x_20 + 1, np.cumsum(data_L_20), "g-", label="L=20")
        plt.loglog(x_50 + 1, np.cumsum(data_L_50), "b-", label="L=50")
        plt.loglog(x_100[::5] + 1, np.cumsum(data_L_100)[::5], "m-x", label="L=100")
        plt.loglog(x_200[::5] + 1, np.cumsum(data_L_200)[::5], "k-o", label="L=200")
        plt.loglog(x_500[::5] + 1, np.cumsum(data_L_500)[::5], "y-s", label="L=500")
        plt.xlabel(r"Cluster size $k$")
        plt.ylabel(r"CDF $C(k)$")
        plt.legend()
        plt.grid()
        plt.savefig("./plots/hexagonal_L_phi_0.5.pdf", dpi=200)


    def plot_tri_L(self):
        """Plots cumulative distribution for differing lattice sizes
        """
        data_L_10 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.10_L10.txt", delimiter=" ")[0::2]
        data_L_20 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.10_L20.txt", delimiter=" ")[0::2]
        data_L_50 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.10_L50.txt", delimiter=" ")[0::2]
        data_L_100 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.10_L100.txt", delimiter=" ")[0::2]
        data_L_200 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.10_L200.txt", delimiter=" ")[0::2]
        data_L_500 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.10_L500.txt", delimiter=" ")[0::2]
        
        """
        data_L_nr_10 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.10_L10.txt", delimiter=" ")[0::2]
        data_L_nr_20 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.10_L20.txt", delimiter=" ")[0::2]
        data_L_nr_50 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.10_L50.txt", delimiter=" ")[0::2]
        data_L_nr_100 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.10_L100.txt", delimiter=" ")[0::2]
        data_L_nr_200 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.10_L200.txt", delimiter=" ")[0::2]
        data_L_nr_500 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.10_L500.txt", delimiter=" ")[0::2]
        """
        x_10 = np.argwhere(data_L_10 != 0)
        x_20 = np.argwhere(data_L_20 != 0)
        x_50 = np.argwhere(data_L_50 != 0)
        x_100 = np.argwhere(data_L_100 != 0)
        x_200 = np.argwhere(data_L_200 != 0)
        x_500 = np.argwhere(data_L_500 != 0)

        """
        x_nr_10 = np.argwhere(data_L_nr_10 != 0)
        x_nr_20 = np.argwhere(data_L_nr_20 != 0)
        x_nr_50 = np.argwhere(data_L_nr_50 != 0)
        x_nr_100 = np.argwhere(data_L_nr_100 != 0)
        x_nr_200 = np.argwhere(data_L_nr_200 != 0)
        x_nr_500 = np.argwhere(data_L_nr_500 != 0)
        """

        data_L_10 = data_L_10[x_10]/np.sum(data_L_10[x_10])         
        data_L_20 = data_L_20[x_20]/np.sum(data_L_20[x_20])
        data_L_50 = data_L_50[x_50]/np.sum(data_L_50[x_50])
        data_L_100 = data_L_100[x_100]/np.sum(data_L_100[x_100])
        data_L_200 = data_L_200[x_200]/np.sum(data_L_200[x_200])
        data_L_500 = data_L_500[x_500]/np.sum(data_L_500[x_500])

        """
        data_L_nr_10 = data_L_nr_10[x_nr_10]/np.sum(data_L_nr_10[x_nr_10])         
        data_L_nr_20 = data_L_nr_20[x_nr_20]/np.sum(data_L_nr_20[x_nr_20])
        data_L_nr_50 = data_L_nr_50[x_nr_50]/np.sum(data_L_nr_50[x_nr_50])
        data_L_nr_100 = data_L_nr_100[x_nr_100]/np.sum(data_L_nr_100[x_nr_100])
        data_L_nr_200 = data_L_nr_200[x_nr_200]/np.sum(data_L_nr_200[x_nr_200])
        data_L_nr_500 = data_L_nr_500[x_nr_500]/np.sum(data_L_nr_500[x_nr_500])
        """

        fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
        plt.rcParams.update({'font.size': 18})
        plt.title("Triangular Lattice")
        plt.loglog(x_10 + 1, np.cumsum(data_L_10), "r-", label="L=10")
        plt.loglog(x_20 + 1, np.cumsum(data_L_20), "g-", label="L=20")
        plt.loglog(x_50 + 1, np.cumsum(data_L_50), "b-", label="L=50")
        plt.loglog(x_100[::5] + 1, np.cumsum(data_L_100)[::5], "m-x", label="L=100")
        plt.loglog(x_200[::5] + 1, np.cumsum(data_L_200)[::5], "k-o", label="L=200")
        plt.loglog(x_500[::5] + 1, np.cumsum(data_L_500)[::5], "y-s", label="L=500")
        plt.xlabel(r"Cluster size $k$")
        plt.ylabel(r"CDF $C(k)$")
        plt.legend()
        plt.grid()
        plt.savefig("./plots/triangular_L.pdf", dpi=200)

        
        # ###############################################################################################
        # The following is for higher density
        data_L_10 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.50_L10.txt", delimiter=" ")[0::2]
        data_L_20 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.50_L20.txt", delimiter=" ")[0::2]
        data_L_50 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.50_L50.txt", delimiter=" ")[0::2]
        data_L_100 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.50_L100.txt", delimiter=" ")[0::2]
        data_L_200 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.50_L200.txt", delimiter=" ")[0::2]
        data_L_500 = np.genfromtxt("./latticesize/triangular_alpha0.01_phi0.50_L500.txt", delimiter=" ")[0::2]
        
        """
        data_L_nr_10 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.50_L10.txt", delimiter=" ")[0::2]
        data_L_nr_20 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.50_L20.txt", delimiter=" ")[0::2]
        data_L_nr_50 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.50_L50.txt", delimiter=" ")[0::2]
        data_L_nr_100 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.50_L100.txt", delimiter=" ")[0::2]
        data_L_nr_200 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.50_L200.txt", delimiter=" ")[0::2]
        data_L_nr_500 = np.genfromtxt("./latticesize/triangular_nr_alpha0.01_phi0.50_L500.txt", delimiter=" ")[0::2]
        """
        x_10 = np.argwhere(data_L_10 != 0)
        x_20 = np.argwhere(data_L_20 != 0)
        x_50 = np.argwhere(data_L_50 != 0)
        x_100 = np.argwhere(data_L_100 != 0)
        x_200 = np.argwhere(data_L_200 != 0)
        x_500 = np.argwhere(data_L_500 != 0)

        """
        x_nr_10 = np.argwhere(data_L_nr_10 != 0)
        x_nr_20 = np.argwhere(data_L_nr_20 != 0)
        x_nr_50 = np.argwhere(data_L_nr_50 != 0)
        x_nr_100 = np.argwhere(data_L_nr_100 != 0)
        x_nr_200 = np.argwhere(data_L_nr_200 != 0)
        x_nr_500 = np.argwhere(data_L_nr_500 != 0)
        """

        data_L_10 = data_L_10[x_10]/np.sum(data_L_10[x_10])         
        data_L_20 = data_L_20[x_20]/np.sum(data_L_20[x_20])
        data_L_50 = data_L_50[x_50]/np.sum(data_L_50[x_50])
        data_L_100 = data_L_100[x_100]/np.sum(data_L_100[x_100])
        data_L_200 = data_L_200[x_200]/np.sum(data_L_200[x_200])
        data_L_500 = data_L_500[x_500]/np.sum(data_L_500[x_500])

        """
        data_L_nr_10 = data_L_nr_10[x_nr_10]/np.sum(data_L_nr_10[x_nr_10])         
        data_L_nr_20 = data_L_nr_20[x_nr_20]/np.sum(data_L_nr_20[x_nr_20])
        data_L_nr_50 = data_L_nr_50[x_nr_50]/np.sum(data_L_nr_50[x_nr_50])
        data_L_nr_100 = data_L_nr_100[x_nr_100]/np.sum(data_L_nr_100[x_nr_100])
        data_L_nr_200 = data_L_nr_200[x_nr_200]/np.sum(data_L_nr_200[x_nr_200])
        data_L_nr_500 = data_L_nr_500[x_nr_500]/np.sum(data_L_nr_500[x_nr_500])
        """

        fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
        plt.rcParams.update({'font.size': 18})
        plt.title("Triangular Lattice")
        plt.loglog(x_10 + 1, np.cumsum(data_L_10), "r-", label="L=10")
        plt.loglog(x_20 + 1, np.cumsum(data_L_20), "g-", label="L=20")
        plt.loglog(x_50 + 1, np.cumsum(data_L_50), "b-", label="L=50")
        plt.loglog(x_100 + 1, np.cumsum(data_L_100), "m-x", label="L=100")
        plt.loglog(x_200 + 1, np.cumsum(data_L_200), "k-o", label="L=200")
        plt.loglog(x_500 + 1, np.cumsum(data_L_500), "y-s", label="L=500")
        plt.xlabel(r"Cluster size $k$")
        plt.ylabel(r"CDF $C(k)$")
        plt.legend()
        plt.grid()
        plt.savefig("./plots/triangular_L_phi_0.5.pdf", dpi=200)


    def plot_alpha(self):
        # square
        if False:
            data_alpha_00001 = np.genfromtxt("./tumblerate/square_alpha0.0001000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_0001 = np.genfromtxt("./tumblerate/square_alpha0.0010000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_001 = np.genfromtxt("./tumblerate/square_alpha0.0100000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_01 = np.genfromtxt("./tumblerate/square_alpha0.1000000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_05 = np.genfromtxt("./tumblerate/square_alpha0.5000000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_1 = np.genfromtxt("./tumblerate/square_alpha1.0000000_phi0.50_L100.txt", delimiter=" ")[0::2]

            """
            data_nr_alpha_0.0001 = np.genfromtxt()
            data_nr_alpha_0.001 = np.genfromtxt()
            data_nr_alpha_0.01 = np.genfromtxt()
            data_nr_alpha_0.1 = np.genfromtxt()
            data_nr_alpha_0.5 = np.genfromtxt()
            data_nr_alpha_1 = np.genfromtxt()
            """
            x_00001 = np.argwhere(data_alpha_00001 != 0)
            x_0001 = np.argwhere(data_alpha_0001 != 0)
            x_001 = np.argwhere(data_alpha_001 != 0)
            x_01 = np.argwhere(data_alpha_01 != 0)
            x_05 = np.argwhere(data_alpha_05 != 0)
            x_1 = np.argwhere(data_alpha_1 != 0)

            """
            x_nr_00001 = np.argwhere(data_nr_alpha_00001 != 0)
            x_nr_0001 = np.argwhere(data_nr_alpha_0001 != 0)
            x_nr_001 = np.argwhere(data_nr_alpha_001 != 0)
            x_nr_01 = np.argwhere(data_nr_alpha_01 != 0)
            x_nr_05 = np.argwhere(data_nr_alpha_05 != 0)
            x_nr_1 = np.argwhere(data_nr_alpha_1 != 0)
            """
            data_alpha_00001 = data_alpha_00001[x_00001]/np.sum(data_alpha_00001[x_00001])
            data_alpha_0001 = data_alpha_0001[x_0001]/np.sum(data_alpha_0001[x_0001])
            data_alpha_001 = data_alpha_001[x_001]/np.sum(data_alpha_001[x_001])
            data_alpha_01 = data_alpha_01[x_01]/np.sum(data_alpha_01[x_01])
            data_alpha_05 = data_alpha_05[x_05]/np.sum(data_alpha_05[x_05])
            data_alpha_1 = data_alpha_1[x_1]/np.sum(data_alpha_1[x_1])

            """
            data_nr_alpha_00001 = data_nr_alpha_00001[x_nr_00001]/np.sum(data_nr_alpha_00001[x_nr_00001])
            data_nr_alpha_0001 = data_nr_alpha_0001[x_nr_0001]/np.sum(data_nr_alpha_0001[x_nr_0001])
            data_nr_alpha_001 = data_nr_alpha_001[x_nr_001]/np.sum(data_nr_alpha_001[x_nr_001])
            data_nr_alpha_01 = data_nr_alpha_01[x_nr_01]/np.sum(data_nr_alpha_01[x_nr_01])
            data_nr_alpha_05 = data_nr_alpha_05[x_nr_05]/np.sum(data_nr_alpha_05[x_nr_05])
            data_nr_alpha_1 = data_nr_alpha_1[x_nr_1]/np.sum(data_nr_alpha_1[x_nr_1])
            """

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, np.cumsum(data_alpha_00001), "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, np.cumsum(data_alpha_0001), "g-", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha.pdf", dpi=200)
            
            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_high.pdf", dpi=200)
            
            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, data_alpha_00001, "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, data_alpha_0001, "g--", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_pdf.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_pdf_high.pdf", dpi=200)

            # ################################################################################
            # The following is for lower denisties

            data_alpha_00001 = np.genfromtxt("./tumblerate/square_alpha0.0001000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_0001 = np.genfromtxt("./tumblerate/square_alpha0.0010000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_001 = np.genfromtxt("./tumblerate/square_alpha0.0100000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_01 = np.genfromtxt("./tumblerate/square_alpha0.1000000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_05 = np.genfromtxt("./tumblerate/square_alpha0.5000000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_1 = np.genfromtxt("./tumblerate/square_alpha1.0000000_phi0.10_L100.txt", delimiter=" ")[0::2]

            """
            data_nr_alpha_0.0001 = np.genfromtxt()
            data_nr_alpha_0.001 = np.genfromtxt()
            data_nr_alpha_0.01 = np.genfromtxt()
            data_nr_alpha_0.1 = np.genfromtxt()
            data_nr_alpha_0.5 = np.genfromtxt()
            data_nr_alpha_1 = np.genfromtxt()
            """
            x_00001 = np.argwhere(data_alpha_00001 != 0)
            x_0001 = np.argwhere(data_alpha_0001 != 0)
            x_001 = np.argwhere(data_alpha_001 != 0)
            x_01 = np.argwhere(data_alpha_01 != 0)
            x_05 = np.argwhere(data_alpha_05 != 0)
            x_1 = np.argwhere(data_alpha_1 != 0)

            """
            x_nr_00001 = np.argwhere(data_nr_alpha_00001 != 0)
            x_nr_0001 = np.argwhere(data_nr_alpha_0001 != 0)
            x_nr_001 = np.argwhere(data_nr_alpha_001 != 0)
            x_nr_01 = np.argwhere(data_nr_alpha_01 != 0)
            x_nr_05 = np.argwhere(data_nr_alpha_05 != 0)
            x_nr_1 = np.argwhere(data_nr_alpha_1 != 0)
            """
            data_alpha_00001 = data_alpha_00001[x_00001]/np.sum(data_alpha_00001[x_00001])
            data_alpha_0001 = data_alpha_0001[x_0001]/np.sum(data_alpha_0001[x_0001])
            data_alpha_001 = data_alpha_001[x_001]/np.sum(data_alpha_001[x_001])
            data_alpha_01 = data_alpha_01[x_01]/np.sum(data_alpha_01[x_01])
            data_alpha_05 = data_alpha_05[x_05]/np.sum(data_alpha_05[x_05])
            data_alpha_1 = data_alpha_1[x_1]/np.sum(data_alpha_1[x_1])

            """
            data_nr_alpha_00001 = data_nr_alpha_00001[x_nr_00001]/np.sum(data_nr_alpha_00001[x_nr_00001])
            data_nr_alpha_0001 = data_nr_alpha_0001[x_nr_0001]/np.sum(data_nr_alpha_0001[x_nr_0001])
            data_nr_alpha_001 = data_nr_alpha_001[x_nr_001]/np.sum(data_nr_alpha_001[x_nr_001])
            data_nr_alpha_01 = data_nr_alpha_01[x_nr_01]/np.sum(data_nr_alpha_01[x_nr_01])
            data_nr_alpha_05 = data_nr_alpha_05[x_nr_05]/np.sum(data_nr_alpha_05[x_nr_05])
            data_nr_alpha_1 = data_nr_alpha_1[x_nr_1]/np.sum(data_nr_alpha_1[x_nr_1])
            """

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, np.cumsum(data_alpha_00001), "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, np.cumsum(data_alpha_0001), "g-", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_phi_low.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_high_phi_low.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, data_alpha_00001, "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, data_alpha_0001, "g--", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_phi_low_pdf.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/square_alpha_high_phi_low_pdf.pdf", dpi=200)

        # hexagonal
        if False:
            data_alpha_00001 = np.genfromtxt("./tumblerate/hexagonal_alpha0.0001000_phi1.00_L100.txt", delimiter=" ")[0::2]
            data_alpha_0001 = np.genfromtxt("./tumblerate/hexagonal_alpha0.0010000_phi1.00_L100.txt", delimiter=" ")[0::2]
            data_alpha_001 = np.genfromtxt("./tumblerate/hexagonal_alpha0.0100000_phi1.00_L100.txt", delimiter=" ")[0::2]
            data_alpha_01 = np.genfromtxt("./tumblerate/hexagonal_alpha0.1000000_phi1.00_L100.txt", delimiter=" ")[0::2]
            data_alpha_05 = np.genfromtxt("./tumblerate/hexagonal_alpha0.5000000_phi1.00_L100.txt", delimiter=" ")[0::2]
            data_alpha_1 = np.genfromtxt("./tumblerate/hexagonal_alpha1.0000000_phi1.00_L100.txt", delimiter=" ")[0::2]

            """
            data_nr_alpha_0.0001 = np.genfromtxt()
            data_nr_alpha_0.001 = np.genfromtxt()
            data_nr_alpha_0.01 = np.genfromtxt()
            data_nr_alpha_0.1 = np.genfromtxt()
            data_nr_alpha_0.5 = np.genfromtxt()
            data_nr_alpha_1 = np.genfromtxt()
            """
            x_00001 = np.argwhere(data_alpha_00001 != 0)
            x_0001 = np.argwhere(data_alpha_0001 != 0)
            x_001 = np.argwhere(data_alpha_001 != 0)
            x_01 = np.argwhere(data_alpha_01 != 0)
            x_05 = np.argwhere(data_alpha_05 != 0)
            x_1 = np.argwhere(data_alpha_1 != 0)

            """
            x_nr_00001 = np.argwhere(data_nr_alpha_00001 != 0)
            x_nr_0001 = np.argwhere(data_nr_alpha_0001 != 0)
            x_nr_001 = np.argwhere(data_nr_alpha_001 != 0)
            x_nr_01 = np.argwhere(data_nr_alpha_01 != 0)
            x_nr_05 = np.argwhere(data_nr_alpha_05 != 0)
            x_nr_1 = np.argwhere(data_nr_alpha_1 != 0)
            """
            data_alpha_00001 = data_alpha_00001[x_00001]/np.sum(data_alpha_00001[x_00001])
            data_alpha_0001 = data_alpha_0001[x_0001]/np.sum(data_alpha_0001[x_0001])
            data_alpha_001 = data_alpha_001[x_001]/np.sum(data_alpha_001[x_001])
            data_alpha_01 = data_alpha_01[x_01]/np.sum(data_alpha_01[x_01])
            data_alpha_05 = data_alpha_05[x_05]/np.sum(data_alpha_05[x_05])
            data_alpha_1 = data_alpha_1[x_1]/np.sum(data_alpha_1[x_1])

            """
            data_nr_alpha_00001 = data_nr_alpha_00001[x_nr_00001]/np.sum(data_nr_alpha_00001[x_nr_00001])
            data_nr_alpha_0001 = data_nr_alpha_0001[x_nr_0001]/np.sum(data_nr_alpha_0001[x_nr_0001])
            data_nr_alpha_001 = data_nr_alpha_001[x_nr_001]/np.sum(data_nr_alpha_001[x_nr_001])
            data_nr_alpha_01 = data_nr_alpha_01[x_nr_01]/np.sum(data_nr_alpha_01[x_nr_01])
            data_nr_alpha_05 = data_nr_alpha_05[x_nr_05]/np.sum(data_nr_alpha_05[x_nr_05])
            data_nr_alpha_1 = data_nr_alpha_1[x_nr_1]/np.sum(data_nr_alpha_1[x_nr_1])
            """

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, np.cumsum(data_alpha_00001), "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, np.cumsum(data_alpha_0001), "g-", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha.pdf", dpi=200)
            
            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_high.pdf", dpi=200)
            
            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, data_alpha_00001, "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, data_alpha_0001, "g--", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_pdf.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_pdf_high.pdf", dpi=200)

            # ################################################################################
            # The following is for lower denisties

            data_alpha_00001 = np.genfromtxt("./tumblerate/hexagonal_alpha0.0001000_phi0.20_L100.txt", delimiter=" ")[0::2]
            data_alpha_0001 = np.genfromtxt("./tumblerate/hexagonal_alpha0.0010000_phi0.20_L100.txt", delimiter=" ")[0::2]
            data_alpha_001 = np.genfromtxt("./tumblerate/hexagonal_alpha0.0100000_phi0.20_L100.txt", delimiter=" ")[0::2]
            data_alpha_01 = np.genfromtxt("./tumblerate/hexagonal_alpha0.1000000_phi0.20_L100.txt", delimiter=" ")[0::2]
            data_alpha_05 = np.genfromtxt("./tumblerate/hexagonal_alpha0.5000000_phi0.20_L100.txt", delimiter=" ")[0::2]
            data_alpha_1 = np.genfromtxt("./tumblerate/hexagonal_alpha1.0000000_phi0.20_L100.txt", delimiter=" ")[0::2]

            """
            data_nr_alpha_0.0001 = np.genfromtxt()
            data_nr_alpha_0.001 = np.genfromtxt()
            data_nr_alpha_0.01 = np.genfromtxt()
            data_nr_alpha_0.1 = np.genfromtxt()
            data_nr_alpha_0.5 = np.genfromtxt()
            data_nr_alpha_1 = np.genfromtxt()
            """
            x_00001 = np.argwhere(data_alpha_00001 != 0)
            x_0001 = np.argwhere(data_alpha_0001 != 0)
            x_001 = np.argwhere(data_alpha_001 != 0)
            x_01 = np.argwhere(data_alpha_01 != 0)
            x_05 = np.argwhere(data_alpha_05 != 0)
            x_1 = np.argwhere(data_alpha_1 != 0)

            """
            x_nr_00001 = np.argwhere(data_nr_alpha_00001 != 0)
            x_nr_0001 = np.argwhere(data_nr_alpha_0001 != 0)
            x_nr_001 = np.argwhere(data_nr_alpha_001 != 0)
            x_nr_01 = np.argwhere(data_nr_alpha_01 != 0)
            x_nr_05 = np.argwhere(data_nr_alpha_05 != 0)
            x_nr_1 = np.argwhere(data_nr_alpha_1 != 0)
            """
            data_alpha_00001 = data_alpha_00001[x_00001]/np.sum(data_alpha_00001[x_00001])
            data_alpha_0001 = data_alpha_0001[x_0001]/np.sum(data_alpha_0001[x_0001])
            data_alpha_001 = data_alpha_001[x_001]/np.sum(data_alpha_001[x_001])
            data_alpha_01 = data_alpha_01[x_01]/np.sum(data_alpha_01[x_01])
            data_alpha_05 = data_alpha_05[x_05]/np.sum(data_alpha_05[x_05])
            data_alpha_1 = data_alpha_1[x_1]/np.sum(data_alpha_1[x_1])

            """
            data_nr_alpha_00001 = data_nr_alpha_00001[x_nr_00001]/np.sum(data_nr_alpha_00001[x_nr_00001])
            data_nr_alpha_0001 = data_nr_alpha_0001[x_nr_0001]/np.sum(data_nr_alpha_0001[x_nr_0001])
            data_nr_alpha_001 = data_nr_alpha_001[x_nr_001]/np.sum(data_nr_alpha_001[x_nr_001])
            data_nr_alpha_01 = data_nr_alpha_01[x_nr_01]/np.sum(data_nr_alpha_01[x_nr_01])
            data_nr_alpha_05 = data_nr_alpha_05[x_nr_05]/np.sum(data_nr_alpha_05[x_nr_05])
            data_nr_alpha_1 = data_nr_alpha_1[x_nr_1]/np.sum(data_nr_alpha_1[x_nr_1])
            """

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, np.cumsum(data_alpha_00001), "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, np.cumsum(data_alpha_0001), "g-", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_phi_low.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_high_phi_low.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, data_alpha_00001, "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, data_alpha_0001, "g--", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_phi_low_pdf.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/hex_alpha_high_phi_low_pdf.pdf", dpi=200)

        # triangular
        if True:
            data_alpha_00001 = np.genfromtxt("./tumblerate/triangular_alpha0.0001000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_0001 = np.genfromtxt("./tumblerate/triangular_alpha0.0010000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_001 = np.genfromtxt("./tumblerate/triangular_alpha0.0100000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_01 = np.genfromtxt("./tumblerate/triangular_alpha0.1000000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_05 = np.genfromtxt("./tumblerate/triangular_alpha0.5000000_phi0.50_L100.txt", delimiter=" ")[0::2]
            data_alpha_1 = np.genfromtxt("./tumblerate/triangular_alpha1.0000000_phi0.50_L100.txt", delimiter=" ")[0::2]

            """
            data_nr_alpha_0.0001 = np.genfromtxt()
            data_nr_alpha_0.001 = np.genfromtxt()
            data_nr_alpha_0.01 = np.genfromtxt()
            data_nr_alpha_0.1 = np.genfromtxt()
            data_nr_alpha_0.5 = np.genfromtxt()
            data_nr_alpha_1 = np.genfromtxt()
            """
            x_00001 = np.argwhere(data_alpha_00001 != 0)
            x_0001 = np.argwhere(data_alpha_0001 != 0)
            x_001 = np.argwhere(data_alpha_001 != 0)
            x_01 = np.argwhere(data_alpha_01 != 0)
            x_05 = np.argwhere(data_alpha_05 != 0)
            x_1 = np.argwhere(data_alpha_1 != 0)

            """
            x_nr_00001 = np.argwhere(data_nr_alpha_00001 != 0)
            x_nr_0001 = np.argwhere(data_nr_alpha_0001 != 0)
            x_nr_001 = np.argwhere(data_nr_alpha_001 != 0)
            x_nr_01 = np.argwhere(data_nr_alpha_01 != 0)
            x_nr_05 = np.argwhere(data_nr_alpha_05 != 0)
            x_nr_1 = np.argwhere(data_nr_alpha_1 != 0)
            """
            data_alpha_00001 = data_alpha_00001[x_00001]/np.sum(data_alpha_00001[x_00001])
            data_alpha_0001 = data_alpha_0001[x_0001]/np.sum(data_alpha_0001[x_0001])
            data_alpha_001 = data_alpha_001[x_001]/np.sum(data_alpha_001[x_001])
            data_alpha_01 = data_alpha_01[x_01]/np.sum(data_alpha_01[x_01])
            data_alpha_05 = data_alpha_05[x_05]/np.sum(data_alpha_05[x_05])
            data_alpha_1 = data_alpha_1[x_1]/np.sum(data_alpha_1[x_1])

            """
            data_nr_alpha_00001 = data_nr_alpha_00001[x_nr_00001]/np.sum(data_nr_alpha_00001[x_nr_00001])
            data_nr_alpha_0001 = data_nr_alpha_0001[x_nr_0001]/np.sum(data_nr_alpha_0001[x_nr_0001])
            data_nr_alpha_001 = data_nr_alpha_001[x_nr_001]/np.sum(data_nr_alpha_001[x_nr_001])
            data_nr_alpha_01 = data_nr_alpha_01[x_nr_01]/np.sum(data_nr_alpha_01[x_nr_01])
            data_nr_alpha_05 = data_nr_alpha_05[x_nr_05]/np.sum(data_nr_alpha_05[x_nr_05])
            data_nr_alpha_1 = data_nr_alpha_1[x_nr_1]/np.sum(data_nr_alpha_1[x_nr_1])
            """

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, np.cumsum(data_alpha_00001), "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, np.cumsum(data_alpha_0001), "g-", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha.pdf", dpi=200)
            
            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_high.pdf", dpi=200)
            
            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, data_alpha_00001, "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, data_alpha_0001, "g--", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_pdf.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_pdf_high.pdf", dpi=200)

            # ################################################################################
            # The following is for lower denisties

            data_alpha_00001 = np.genfromtxt("./tumblerate/triangular_alpha0.0001000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_0001 = np.genfromtxt("./tumblerate/triangular_alpha0.0010000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_001 = np.genfromtxt("./tumblerate/triangular_alpha0.0100000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_01 = np.genfromtxt("./tumblerate/triangular_alpha0.1000000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_05 = np.genfromtxt("./tumblerate/triangular_alpha0.5000000_phi0.10_L100.txt", delimiter=" ")[0::2]
            data_alpha_1 = np.genfromtxt("./tumblerate/triangular_alpha1.0000000_phi0.10_L100.txt", delimiter=" ")[0::2]

            """
            data_nr_alpha_0.0001 = np.genfromtxt()
            data_nr_alpha_0.001 = np.genfromtxt()
            data_nr_alpha_0.01 = np.genfromtxt()
            data_nr_alpha_0.1 = np.genfromtxt()
            data_nr_alpha_0.5 = np.genfromtxt()
            data_nr_alpha_1 = np.genfromtxt()
            """
            x_00001 = np.argwhere(data_alpha_00001 != 0)
            x_0001 = np.argwhere(data_alpha_0001 != 0)
            x_001 = np.argwhere(data_alpha_001 != 0)
            x_01 = np.argwhere(data_alpha_01 != 0)
            x_05 = np.argwhere(data_alpha_05 != 0)
            x_1 = np.argwhere(data_alpha_1 != 0)

            """
            x_nr_00001 = np.argwhere(data_nr_alpha_00001 != 0)
            x_nr_0001 = np.argwhere(data_nr_alpha_0001 != 0)
            x_nr_001 = np.argwhere(data_nr_alpha_001 != 0)
            x_nr_01 = np.argwhere(data_nr_alpha_01 != 0)
            x_nr_05 = np.argwhere(data_nr_alpha_05 != 0)
            x_nr_1 = np.argwhere(data_nr_alpha_1 != 0)
            """
            data_alpha_00001 = data_alpha_00001[x_00001]/np.sum(data_alpha_00001[x_00001])
            data_alpha_0001 = data_alpha_0001[x_0001]/np.sum(data_alpha_0001[x_0001])
            data_alpha_001 = data_alpha_001[x_001]/np.sum(data_alpha_001[x_001])
            data_alpha_01 = data_alpha_01[x_01]/np.sum(data_alpha_01[x_01])
            data_alpha_05 = data_alpha_05[x_05]/np.sum(data_alpha_05[x_05])
            data_alpha_1 = data_alpha_1[x_1]/np.sum(data_alpha_1[x_1])

            """
            data_nr_alpha_00001 = data_nr_alpha_00001[x_nr_00001]/np.sum(data_nr_alpha_00001[x_nr_00001])
            data_nr_alpha_0001 = data_nr_alpha_0001[x_nr_0001]/np.sum(data_nr_alpha_0001[x_nr_0001])
            data_nr_alpha_001 = data_nr_alpha_001[x_nr_001]/np.sum(data_nr_alpha_001[x_nr_001])
            data_nr_alpha_01 = data_nr_alpha_01[x_nr_01]/np.sum(data_nr_alpha_01[x_nr_01])
            data_nr_alpha_05 = data_nr_alpha_05[x_nr_05]/np.sum(data_nr_alpha_05[x_nr_05])
            data_nr_alpha_1 = data_nr_alpha_1[x_nr_1]/np.sum(data_nr_alpha_1[x_nr_1])
            """

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, np.cumsum(data_alpha_00001), "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, np.cumsum(data_alpha_0001), "g-", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_phi_low.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, np.cumsum(data_alpha_001), "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, np.cumsum(data_alpha_01), "m-x", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, np.cumsum(data_alpha_05), "k-o", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, np.cumsum(data_alpha_1), "y-s", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"CDF $C(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_high_phi_low.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_00001 + 1, data_alpha_00001, "r-", label=r"$\alpha=0.0001$")
            plt.loglog(x_0001 + 1, data_alpha_0001, "g--", label=r"$\alpha=0.001$")
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_phi_low_pdf.pdf", dpi=200)

            fig = plt.figure(figsize=(9.6, 7.2), tight_layout=True)
            plt.rcParams.update({'font.size': 18})
            plt.loglog(x_001 + 1, data_alpha_001, "b-", label=r"$\alpha=0.01$")
            plt.loglog(x_01 + 1, data_alpha_01, "m--", label=r"$\alpha=0.1$")
            plt.loglog(x_05 + 1, data_alpha_05, "k-", label=r"$\alpha=0.5$")
            plt.loglog(x_1 + 1, data_alpha_1, "y--", label=r"$\alpha=1$")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"PDF $p(k)$")
            plt.legend()
            plt.grid()
            plt.savefig("./plots/tri_alpha_high_phi_low_pdf.pdf", dpi=200)

if __name__ == "__main__":
    dist = distribution()

    dist.plot_alpha()
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

        
