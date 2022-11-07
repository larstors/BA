import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

l = 1

sq = np.genfromtxt("square_alpha0.010_phi0.47_L100_1.txt" , delimiter=" ")[2::2]
tr = np.genfromtxt("triangular_alpha0.010_phi0.40_L100_1.txt" , delimiter=" ")[2::2]
hx = np.genfromtxt("hexagonal_alpha0.010_phi1.10_L100_1.txt", delimiter=" ")[2::2]


ks = np.arange(2,len(sq)+2)[sq>l]
kt = np.arange(2,len(tr)+2)[tr>l]
kh = np.arange(2,len(hx)+2)[hx>l]

def tpow(x, b, c):
        N = -b * np.log(x)-c*x
        return N - np.log(np.sum(np.exp(N)))

sq = sq[sq>l]
tr = tr[tr>l]
hx = hx[hx>l]


Zs = np.sum(sq)
Zt = np.sum(tr)
Zh = np.sum(hx)




ps = [2.0,1]
ps = opt.curve_fit(tpow, ks, np.log(sq/Zs), ps, bounds=([0, 1/max(ks)],[np.inf, 1]))[0]

pt = [2.0,1]
pt = opt.curve_fit(tpow, kt, np.log(tr/Zt), pt, bounds=([0, 1/max(kt)],[np.inf, 1]))[0]

ph = [2.0,1]
ph = opt.curve_fit(tpow, kh, np.log(hx/Zh), ph, bounds=([0, 1/max(kh)],[np.inf, 1]))[0]

twodp = lambda x: f'{float(f"{x:.2g}"):g}'

fig = plt.figure()
plt.loglog(ks, np.exp(tpow(ks,*ps)), '-', label=f"$A^{{-{twodp(ps[0])}}} e^{{-A/{twodp(1/ps[1])}}}$")
plt.loglog(ks, sq/Zs, 'o', label="Sim.")
plt.legend()
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.grid()
plt.savefig("sq_cl_dist_n1.pdf", dpi=100)


fig1 = plt.figure()
plt.loglog(kh, np.exp(tpow(kh,*ph)), '-', label=f"$A^{{-{twodp(ph[0])}}} e^{{-A/{twodp(1/ph[1])}}}$")
plt.loglog(kh, hx/Zh, 'o', label="Sim.")
plt.legend()
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.grid()
plt.savefig("hx_cl_dist_n1.pdf", dpi=100)


fig2 = plt.figure()
plt.loglog(kt, np.exp(tpow(kt,*pt)), '-', label=f"$A^{{-{twodp(pt[0])}}} e^{{-A/{twodp(1/pt[1])}}}$")
plt.loglog(kt, tr/Zt, 'o', label="Sim.")
plt.legend()
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.grid()
plt.savefig("tr_cl_dist_n1.pdf", dpi=100)





sq = np.genfromtxt("square_alpha0.010_phi0.12_L100_1.txt" , delimiter=" ")[2::2]
tr = np.genfromtxt("triangular_alpha0.010_phi0.10_L100_1.txt" , delimiter=" ")[2::2]
hx = np.genfromtxt("hexagonal_alpha0.010_phi0.28_L100_1.txt", delimiter=" ")[2::2]


ks = np.arange(2,len(sq)+2)[sq>0]
kt = np.arange(2,len(tr)+2)[tr>0]
kh = np.arange(2,len(hx)+2)[hx>0]

def tpow(x, b, c):
        N = -b * np.log(x)-c*x
        return N - np.log(np.sum(np.exp(N)))

sq = sq[sq>0]
tr = tr[tr>0]
hx = hx[hx>0]


Zs = np.sum(sq)
Zt = np.sum(tr)
Zh = np.sum(hx)


print(1/max(kh))

ps = [2.0,1]
ps = opt.curve_fit(tpow, ks, np.log(sq/Zs), ps, bounds=([0, 1/max(ks)],[np.inf, 1]))[0]

pt = [2.0,1]
pt = opt.curve_fit(tpow, kt, np.log(tr/Zt), pt, bounds=([0, 1/max(kt)],[np.inf, 1]))[0]

ph = [2,1]
ph = opt.curve_fit(tpow, kh, np.log(hx/Zh), ph, bounds=([0, 1/max(kh)],[np.inf, 1]))[0]

twodp = lambda x: f'{float(f"{x:.2g}"):g}'



fig = plt.figure()
plt.loglog(ks, np.exp(tpow(ks,*ps)), '-', label=f"$A^{{-{twodp(ps[0])}}} e^{{-A/{twodp(1/ps[1])}}}$")
plt.loglog(ks, sq/Zs, 'o', label="Sim.")
plt.legend()
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.grid()
plt.savefig("sq_cl_dist_n1_low.pdf", dpi=100)


fig1 = plt.figure()
plt.loglog(kh, np.exp(tpow(kh,*ph)), '-', label=f"$A^{{-{twodp(ph[0])}}} e^{{-A/{twodp(1/ph[1])}}}$")
plt.loglog(kh, hx/Zh, 'o', label="Sim.")
plt.legend()
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.grid()
plt.savefig("hx_cl_dist_n1_low.pdf", dpi=100)


fig2 = plt.figure()
plt.loglog(kt, np.exp(tpow(kt,*pt)), '-', label=f"$A^{{-{twodp(pt[0])}}} e^{{-A/{twodp(1/pt[1])}}}$")
plt.loglog(kt, tr/Zt, 'o', label="Sim.")
plt.legend()
plt.xlabel(r"$A$")
plt.ylabel(r"$p(A)$")
plt.grid()
plt.savefig("tr_cl_dist_n1_low.pdf", dpi=100)


