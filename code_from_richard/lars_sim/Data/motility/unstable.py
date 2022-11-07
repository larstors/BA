import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt
plt.rcParams.update({'font.size': 20})

tr = np.genfromtxt("tri_long.txt", delimiter=" ")

t = tr[:, 0]
val = tr[:, 1]
"""
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), sharey=True, gridspec_kw={'width_ratios': [9, 4]})
#plt.tight_layout()
ax[0].plot(t, val, "o")
ax[0].set_xlabel(r"$t$")
ax[0].set_ylabel(r"$J$")
ax[0].set_xticklabels([r"0", r'$2\cdot 10^6$', r'$4\cdot 10^6$', r'$6\cdot 10^6$', r'$8\cdot 10^6$', r'$10\cdot 10^6$'])
ax[0].grid()

ax[1].hist(val, bins=30, orientation="horizontal")
ax[1].set_xscale("log")



plt.savefig("tri_unstable_long.pdf", dpi=200, bbox_inches='tight')

"""

tr_long = np.genfromtxt("tri_longlong.txt", delimiter=" ")

dur = []

start = 0
finish = 0
k1 = 0
k2 = 0

for i in range(len(tr_long)):

    if (i == 0 and tr_long[i, 1] < 0.05):
        start = tr_long[i, 0]
        k1 = 1
    elif k1 == 0 and tr_long[i, 1] < 0.05 and tr_long[i-1, 1] > 0.05:
        start = tr_long[i, 0]
        k1 = 1
    
    if i != len(tr_long) - 1 and tr_long[i, 1] < 0.05 and tr_long[i+1, 1] > 0.05:
        finish = tr_long[i, 0]
        k2 = 1
    

    if (k1 != 0 and k2 != 0):
        k1 = 0
        k2 = 0
        dur.append(finish - start)


avg = np.mean(dur)
var = np.var(dur)
print(avg, np.sqrt(var))

durh = []

start = 0
finish = 0
k1 = 0
k2 = 0

for i in range(len(tr_long)):

    if (i == 0 and tr_long[i, 1] > 0.05):
        start = tr_long[i, 0]
        k1 = 1
    elif k1 == 0 and tr_long[i, 1] > 0.05 and tr_long[i-1, 1] < 0.05:
        start = tr_long[i, 0]
        k1 = 1
    
    if i != len(tr_long) - 1 and tr_long[i, 1] > 0.05 and tr_long[i+1, 1] < 0.05:
        finish = tr_long[i, 0]
        k2 = 1
    

    if (k1 != 0 and k2 != 0):
        k1 = 0
        k2 = 0
        durh.append(finish - start)


avg = np.mean(durh)
var = np.var(durh)
print(avg, np.sqrt(var))


n_bins = int(len(dur)/2)

hist, bin_edges = np.histogram(dur, bins=n_bins)
x = []
tx = []
for i in range(len(hist)):
    tx.append((bin_edges[i+1] + bin_edges[i])/2)
    x.append((i+1))

histh, bin_edgesh = np.histogram(durh, bins=n_bins)
xh = []
txh = []
for i in range(len(hist)):
    txh.append((bin_edgesh[i+1] + bin_edgesh[i])/2)
    xh.append((i+1))

x = np.array(x)
xh = np.array(xh)

xx = x[hist>0]
xxh = xh[histh>0]

xhh = hist[hist>0]
xhhh = histh[histh>0]

def fit(x, a):
    return a*np.exp(-a*x) 

p = opt.curve_fit(fit, xx, xhh/sum(hist))[0]
ph = opt.curve_fit(fit, xxh, xhhh/sum(histh))[0]

print(p, ph)
print(histh[histh>0])

txx = np.linspace(min(tx), max(tx), 100)
txxh = np.linspace(min(txh), max(txh), 100)
xx = np.linspace(min(x), max(x), 100)
xxh = np.linspace(min(xh), max(xh), 100)



fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(20, 10))

colors = [r'$J<0.05$', r'$J>0.05$']
ax[0].plot(txx, fit(xx, *p), "-r", label="exp. fit")
ax[0].plot(tx, hist/sum(hist), "s")
ax[0].bar(tx, hist/sum(hist), label=colors[0], width=350000, align="center", color="green")
ax[0].legend(prop={'size': 10})
ax[0].set_xlabel(r"$t$")


ax[1].plot(txxh, fit(xxh, *ph), "-r", label="exp. fit")
ax[1].plot(txh, histh/sum(histh), "s")
ax[1].bar(txh, histh/sum(histh), label=colors[1], width=40000, align="center", color="green")
ax[1].legend(prop={'size': 10})
ax[1].set_xlabel(r"$t$")


plt.savefig("hist_long_tri.pdf", dpi=200)




