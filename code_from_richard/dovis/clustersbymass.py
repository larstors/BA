"""
Comparison of cluster distribution time series by mass
"""

# Rather than averaging across multiple time series, this compares multiple
# time series, locating the time point that corresponds to a particular mass
# (mass = mean cluster size)

from dovis.core import *
import sys, numpy as np, matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


def cdists_bymass(data, dist, masses=None):
    """
    Reads the sequence cluster time series from data.
    The corresponding element of the returned sequence is a tuple (parameters, timeseries)
    where (the inappropriately-named) timeseries is a dict with the mass as the key and
    value the tuple (time, histogram). Herem the time is the first time at which a given
    mass is exceeded, and histogram the corresponding normalised 1-based histogram.

    If the masses argument is provided, it should be a list of target masses at
    which histograms are captured. dist should be one of pdist, vdist
    """

    retval = []
    timeseries = None
    masses = None if masses is None else sorted(masses)

    P = None
    for S in data:
        if isinstance(S, Parameters):
            # This starts a new timeseries
            if P is not None:
                retval += [(P,timeseries)]
            P = S
            timeseries = {}
            mseq = None if masses is None else masses[:]
        else:
            hist = np.array(getattr(S, dist)) # Copy the histogram
            hist = hist / np.sum(hist) # Normalise
            mass = np.sum( hist * (1+np.arange(len(hist))) )
            if mseq is None:
                if mass not in timeseries:
                    timeseries[mass] = (S.time, hist)
            elif mseq and mass >= mseq[0]:
                timeseries[mass] = (S.time, hist)
                # Leapfrog any requested masses that are smaller than that just achieved
                while mseq and mass >= mseq[0]:
                    mseq.pop(0)

    if timeseries is not None:
        retval += [(P,timeseries)]

    return retval

def cdist_compar(args, dist):
    """
    Compare multiple distribution functions with the same mass
    """

    data = process_files(args.file)
    masses = None if args.means is None else map(float, args.means.split(','))
    twosf = lambda x: f'{float(f"{x:.2g}"):g}'

    tss = cdists_bymass(data, dist, masses)
    max_k = 0
    for P, timeseries in tss:
        for mass, (time, hist) in timeseries.items():
            plt.plot(1+np.arange(len(hist)),hist,'o' if P.beta is None else ':' if P.beta == 0 else '-',label=f"$\\bar{{k}}={twosf(mass)}$ $t={time:g}$")
            if len(hist) >= max_k:
                max_k = len(hist)+1

    # Plot fixed asymptotes
    ks = np.arange(0,max_k)
    for alpha, rho in set( (P.alpha, P.density) for P, timeseries in tss ):
        b = 0.5*alpha**2
        r = (np.sqrt(b**2+4*b*(alpha-b)*rho*(1-rho))-b) / 2/(alpha-b)
        g = b/alpha * (1-rho-r)/r
        plt.plot(ks, g*(1+g)**(-ks), '--')

    plt.xlabel(r'Cluster size, $k$')
    plt.ylabel(r'Frequency, $P(k)$')
    plt.xlim(adjust_range(plt.xlim(), args.xrange))
    plt.ylim(adjust_range(plt.ylim(), args.yrange))
    plt.legend()

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

pcomp = lambda args: cdist_compar(args, 'pdist')
vcomp = lambda args: cdist_compar(args, 'vdist')

def cdist_movie(args, dist):
    """
    Animation showing mass-indexed distribution functions over time
    """
    tss = cdists_bymass(process_files(args.file), dist)

    # Establish a reasonable range of masses over which comparison is possible
    m0, m1 = None, None
    for P, timeseries in tss:
        minm, maxm = min(timeseries.keys()), max(timeseries.keys())
        if m0 is None or minm > m0:
            m0 = minm
        if m1 is None or maxm < m1:
            m1 = maxm

    if m1 < m0:
        print("Masses did not overlap, no movie produced")
        return

    masses = np.arange(m0, m1+args.mstep, args.mstep)
    histseqs = []
    paramseq = []

    # Extract the histograms that lie closest to the desired masses
    for P, timeseries in tss:
        paramseq += [ P ]
        histseqs += [ np.array(list(timeseries.values()),dtype=object)[np.argmin(abs(masses-np.array(list(timeseries.keys()))[:,np.newaxis]),axis=0)] ]
    del(tss) # This will release any histograms we won't be plotting

    # Get maximum k and histogram heights
    max_k = max( len(hist) for histseq in histseqs for (time, hist) in histseq )
    max_h = max( max(hist) for histseq in histseqs for (time, hist) in histseq )

    # First set up the figure, the axis, and the plot element we want to animate
    fig, ax = plt.subplots(1,1)
    plots = [ ax.plot([], [], 'o' if P.beta is None else ':' if P.beta == 0 else '-')[0] for P in paramseq ]

    txt = ax.text(0.5,0.95,'',ha='center',transform=ax.transAxes)
    ax.set_xlim(adjust_range((0,max_k), args.xrange))
    ax.set_ylim(adjust_range((0,max_h), args.yrange))
    plt.xlabel(r'Cluster size, $k$')
    plt.ylabel(r'Frequency, $P(k)$')

    # Plot fixed asymptotes
    ks = np.arange(0,max_k)
    for alpha, rho in set( (P.alpha, P.density) for P in paramseq ):
        b = 0.5*alpha**2
        r = (np.sqrt(b**2+4*b*(alpha-b)*rho*(1-rho))-b) / 2/(alpha-b)
        g = b/alpha * (1-rho-r)/r
        ax.plot(ks, g*(1+g)**(-ks), '--')


    threesf = lambda x: f'{float(f"{x:.3g}"):g}'
    def update(n):
        for histseq, plot in zip(histseqs, plots):
            time, hist = histseq[n]
            ks = np.arange(1,len(hist)+1)
            plot.set_data(ks, hist)
        txt.set_text(f"$\\bar{{k}}\\approx{threesf(masses[n])}$")
        return plots + [txt]

    ani = FuncAnimation(fig, update, frames=len(masses), interval=50, repeat=True, blit=True, save_count=sys.maxsize)

    if args.fig is None:
        plt.show()
    else:
        ani.save(args.fig)

pmmov = lambda args: cdist_movie(args, 'pdist')
vmmov = lambda args: cdist_movie(args, 'vdist')

def time_of_mass(args, dist):
    """
    Compare the times at which particlar masses are reached
    """

    data = process_files(args.file)
    masses = None if args.means is None else map(float, args.means.split(','))

    for P, timeseries in cdists_bymass(data, dist, masses):
        x,y = zip(*[(mass, time) for mass, (time,hist) in timeseries.items()])
        plt.semilogy(x,y)

    plt.xlabel(r'mean cluster size, $\bar{k}$')
    plt.ylabel(r'time, $t$')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

ptom = lambda args: time_of_mass(args, 'pdist')
vtom = lambda args: time_of_mass(args, 'vdist')

def time_of_time(args, dist):
    """
    Plot of relative time that different series reach the same mass
    """
    tss = cdists_bymass(process_files(args.file), dist)

    # Establish a reasonable range of masses over which comparison is possible
    m0, m1 = None, None
    for P, timeseries in tss:
        minm, maxm = min(timeseries.keys()), max(timeseries.keys())
        if m0 is None or minm > m0:
            m0 = minm
        if m1 is None or maxm < m1:
            m1 = maxm

    if m1 < m0:
        print("Masses did not overlap, no movie produced")
        return

    masses = np.arange(m0, m1+args.mstep, args.mstep)
    timeseqs = []
    paramseq = []

    # Extract the times that lie closest to the desired masses
    for P, timeseries in tss:
        paramseq += [P]
        timeseqs += [ np.array([time for (time,hist) in timeseries.values()])[np.argmin(abs(masses-np.array(list(timeseries.keys()))[:,np.newaxis]),axis=0)] ]
    del(tss) # This will release redundant data

    # Plot ratio of each time
    baseline = timeseqs.pop(0)
    paramseq.pop(0)
    mask = baseline>0
    for P, times in zip(paramseq,timeseqs):
        plt.plot(masses[mask], times[mask]/baseline[mask], 'o' if P.beta is None else ':' if P.beta == 0 else '-')

    plt.xlabel(r'Mean cluster size, $\bar{k}$')
    plt.ylabel(r'Relative time')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

ptot = lambda args: time_of_time(args, 'pdist')
vtot = lambda args: time_of_time(args, 'vdist')
