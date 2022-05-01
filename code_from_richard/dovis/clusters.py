"""
Cluster size distributions
"""

from dovis.core import *
import sys, numpy as np, matplotlib.pyplot as plt, itertools
from matplotlib.animation import FuncAnimation
from scipy.optimize import curve_fit
from scipy.special import gamma

def cluster_dist(args, *cycle):
    """
    Plot full distribution of cluster sizes
    """
    cycle = itertools.cycle(cycle)

    data = process_files(args.file)

    plots = []
    params = []

    def tpow(x, b, c):
        N = -b * np.log(x)-c*x
        return N - np.log(np.sum(np.exp(N)))

    # Accumulate the histogram, and normalise at the end, I guess
    hist = np.array([0])
    for S in data:
        if isinstance(S, Parameters):
            P = S # Current set of parameters
        elif isinstance(S, ClusterDistribution):
            dname = next(cycle)
            hist = getattr(S,dname)
            ks = np.arange(1,len(hist)+1)[hist>0] # Cluster size
            hist = hist[hist>0] # Corresponding number
            Z = np.sum(hist) # Total number
            plots += plt.loglog(ks, hist/Z, 'o' if dname == 'pdist' else 'x')
            params += [(P.alpha, P.density if dname == 'pdist' else P.vdensity, P.strshape)]

    if plots:
        # Quick and dirty fit to the last data series
        twodp = lambda x: f'{float(f"{x:.2g}"):g}'
        p = [2.0,5.0/max(ks)]
        p = curve_fit(tpow, ks, np.log(hist/Z), p, bounds=(0,np.inf))[0]
        plt.loglog(ks, np.exp(tpow(ks,*p)), '-', label=f"$k^{{-{twodp(p[0])}}} e^{{-k/{twodp(1/p[1])}}}$")
        plt.legend()


    # Identify any parameters that vary between series and put these in a legend
    if params:
        pnames = [r'\alpha',r'\rho',r'L']
        for n in range(len(pnames)):
            if all( x[n] == params[0][n] for x in params ):
                pnames[n] = None
        if any(m is not None for m in pnames):
            for plot, pcomb in zip(plots, params):
                plot.set_label(' '.join( f'${x}={y}$' for x,y in zip(pnames,[pcomb[0],float(pcomb[1]),re.sub('x',r'\\times',pcomb[2])]) if x is not None ))
            plt.legend()

    plt.xlabel(r'Cluster size, $k$')
    plt.ylabel(r'Frequency, $P(k)$')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

# Subtypes according to the order in which we expect clusters
pdist = lambda args: cluster_dist(args, 'pdist')
vdist = lambda args: cluster_dist(args, 'vdist')
pvdist = lambda args: cluster_dist(args, 'pdist', 'vdist')
vpdist = lambda args: cluster_dist(args, 'vdist', 'pdist')

def cluster_mean(args, ctype):
    """
    Plot mean cluster sizes as function of system volume.
    ctype should be 'particle' or 'vacancy'
    """
    data = process_files(args.file)

    means = {}

    for S in data:
        if isinstance(S, Parameters):
            P = S # Current set of parameters
            rho = P.density if ctype == 'particle' else 1-P.density
            if (P.alpha, rho) not in means:
                means[P.alpha, rho] = []
        elif isinstance(S, ClusterDistribution):
            means[P.alpha,rho] += [ (P.volume, S.pmean if ctype == 'particle' else S.vmean)]

    if means:
        # Identify any parameters that vary between series for legend construction purposes
        pnames = [r'\alpha',r'\rho']
        base = next(iter(means))
        for n in range(len(pnames)):
            if all( x[n] == base[n] for x in means ):
                pnames[n] = None
        for pcomb, dataset in means.items():
            x, y = zip(*sorted(dataset,key = lambda vm: vm[0]))
            plt.plot(x,y,'o',label=' '.join( f'${u}={v}$' for u,v in zip(pnames,[pcomb[0],float(pcomb[1])]) if u is not None ))
        if any( x is not None for x in pnames ):
            plt.legend()

    plt.xlabel(r'System volume, $V$')
    plt.ylabel(r'Mean cluster size, $\langle k \rangle$')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

pmean = lambda args: cluster_mean(args, 'particle')
vmean = lambda args: cluster_mean(args, 'vacancy')

def cdists_timeseries(data, dist, normed=True):
    """
    Transform each sequence of cluster distributions into a time-ordered map from
    time to distribution.

    Returns [(params, dists)] where dists[t] is a histogram, normalised unless
    normed=False.
    NB: t can be None, in which case this is a stationary distribution.
    """

    P = None
    retval = []

    for S in data:
        if isinstance(S, Parameters):
            if P is not None and dists: # Don't bother with empty series
                retval += [(P, dists)]
            P = S # current set of parameters
            dists = {} # current sequence of distributions
        elif isinstance(S, ClusterDistribution):
            # Insert a normalised histogram into dists
            hist = np.array(getattr(S, dist))
            if normed:
                dists[S.time] = hist / np.sum(hist)
            else:
                dists[S.time] = hist

    if P is not None and dists:
        retval += [(P, dists)]
    return retval

def csds_of_t(args, dist):
    """
    Static plot of cluster size distributions at different times
    dist should be one of pdist or vdist
    """

    kmax = 0
    labels = 0
    cds = cdists_timeseries(process_files(args.file), dist, normed=not args.raw)
    for P, hists in cds:
        # Determine the times at which we will plot the distribution
        times = [ t for t in hists if t is not None ]
        if not times: continue
        if args.times is not None:
            tset = set()
            for time in ( float(x) for x in args.times.split(',') ):
                # Find the time that's closest to args.time
                tset |= {times[np.argmin(np.abs(np.array(times)-time))]}
            times = sorted(tset)
        else:
            times = sorted(times)

        for t in times:
            hist = hists[t]
            kmax = max(kmax, len(hist))
            ks = np.arange(1,len(hist)+1)
            if args.raw:
                # Rescale to continuum form
                delta = np.sqrt(P.alpha/P.beta/P.density/(1.0-P.density)) * np.sum(hist)
                ks = np.sqrt(P.beta*(1.0-P.density)/P.alpha/P.density) * ks
                hist = P.alpha / P.beta / (1.0-P.density) / delta * hist

                # One-parameter qss
                lm = 1.0/(1.0+delta)
                def qss(x,m):
                    A  = lm / (2*lm - m - 1)
                    pm = np.sqrt(complex(lm**2 + 4*(lm-m)))
                    kp = (lm+pm)/2
                    km = (lm-pm)/2
                    Ap = (kp*(-A-m)-(m-lm)*km*(1-A))/lm/(kp-km)
                    Am = (km*(-A-m)-(m-lm)*kp*(1-A))/lm/(km-kp)
                    return A*np.exp(-x) + (Ap*np.exp(-kp*x)+Am*np.exp(-km*x)).real

                p = curve_fit(qss, ks, hist, (1.5,))[0]
                plt.plot(ks, qss(ks, *p), '--', label=rf"$\Delta={delta:.2g}$ $m={p[0]:.2f}$")
                labels += 1

            # TODO: adjust labels only according to what changes between data sets
            plt.plot(ks, hist, 'o' if P.beta is None else ':' if P.beta == 0 else '-', label=f"$t={t}$")
            labels += 1

    if args.raw:
        plt.axhline(c='k',lw=0.5)
    else:
        for alpha, rho in set( (P.alpha, P.density) for P, hists in cds ):
            # Search for a stationary distribution in the input
            for P, hists in cds:
                if P.alpha == alpha and P.density == rho and None in hists:
                    hist = hists[None]
                    kmax = max(kmax, len(hist))
                    # TODO: identfy in some way?
                    plt.plot(ks, hist, '--')
                    break
            else:
                # Default exponential histogram from single-particle theory
                # NB: we could deprecate this and require an explicit stationary
                # distribution to be provided?
                b = 0.5*alpha**2
                r = (np.sqrt(b**2+4*b*(alpha-b)*rho*(1-rho))-b) / 2/(alpha-b)
                g = b/alpha * (1-rho-r)/r
                ks = np.arange(1,len(hist)+1)
                plt.plot(ks, g*(1+g)**(-ks), '--')

    plt.xlim(adjust_range(plt.xlim(), args.xrange))
    plt.ylim(adjust_range(plt.ylim(), args.yrange))
    if labels <= 6:
        plt.legend() # For now, quick and dirty legend-size restriction

    if args.raw:
        plt.xlabel(r'Scaled cluster size, $\sqrt{\beta(1-\phi)/\alpha\phi} k$')
        plt.ylabel(r'Scaled deficit, $\alpha \epsilon_k / \beta (1-\phi) \Delta$')
    else:
        plt.xlabel(r'Cluster size, $k$')
        plt.ylabel(r'Frequency, $P(k)$')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

poft = lambda args: csds_of_t(args, 'pdist')
voft = lambda args: csds_of_t(args, 'vdist')

def linfit_moft(t, m, asymp=np.log):
    """
    Perform linear fit with assumed asymptotic form to m(t) data.

    Returns a tuple containing:
        - coefficients in the geometric formula
        - fitted t values
        - fitted m values
        - sum of residuals
    """

    x = np.array(t)
    y = m[x>1]
    x = x[x>1]
    coeff = np.vstack( (1/np.sqrt(x), np.ones_like(x), 1/asymp(x)) ).T
    res = np.linalg.lstsq(coeff, 1/y, rcond=None)[0]

    return (
        [1/res[0], 1/res[1], res[2]/res[1]**2],
        x,
        1/(coeff@res),
        np.linalg.norm(y-1/(coeff@res))
    )

def geofit_moft(t, m, asymp=np.log, p0=[1,0,0], positive=True):
    """
    Perform nonlinear fit with assumed asymptotic form to m(t) data.

    p0 should be the coefficients (e.g., obtained from linfit_moft)
    positive determines if only positive coefficients should be tried

    Returns a tuple containing:
        - fitted coefficients in the geometric formula
        - fitted t values
        - fitted m values
        - sum of residuals
    """

    x = np.array(t)
    y = m[x>1]
    x = x[x>1]

    if positive: # Force parameters to be positive
        p0 = [ max(0,_) for _ in p0 ]
        bounds = (0,np.inf)
    else:
        bounds = None

    def geo(x, a,m8,b):
        return 1 / ( 1/(a*np.sqrt(x)) + 1/(m8-b/asymp(x)) )

    p,_ = curve_fit(geo, x, y, p0=p0, bounds=bounds, ftol=None)
    yf = geo(x, *p)

    return ( p, x, yf, np.linalg.norm(y-yf) )

def claw(args):
    """
    Variation of density with time (should be constant, this serves as a check)
    """

    cds = cdists_timeseries(process_files(args.file), 'pdist', normed=False)
    for P, hists in cds:
        times = np.array(sorted(list(hists.keys())))

        rho = np.zeros(len(times))
        for n,t in enumerate(times):
            hist = hists[t]
            rho[n] = np.sum(np.arange(1,len(hist)+1) * hist)

        times, rho = prune(times, rho)

        # TODO: add distinguishing labels as required
        plt.plot(times, rho, 'o' if P.beta is None else ':' if P.beta == 0 else '-')

    plt.xlabel(r'time, $t$')
    plt.ylabel(r'density, $\rho(t)$')
    plt.legend()
    plt.xlim(adjust_range(plt.xlim(), args.xrange))
    plt.ylim(adjust_range(plt.ylim(), args.yrange))


    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

def mass_oft(args, dist):
    """
    Variation of mass with time
    """

    # Assumed asymptote
    asymp = np.log if args.asymp==0 else lambda u: u**args.asymp

    cds = cdists_timeseries(process_files(args.file), dist)
    for P, hists in cds:
        times = np.array(sorted(list(hists.keys())))

        cmeans = np.zeros(len(times))
        for n,t in enumerate(times):
            hist = hists[t]
            cmeans[n] = np.sum(np.arange(1,len(hist)+1) * hist)

        times, cmeans = prune(times, cmeans)

        # TODO: add distinguishing labels as required
        plt.plot(times, cmeans, 'o' if P.beta is None else ':' if P.beta == 0 else '-')

        if P.beta is None:
            # Attempt to fit the full timeseries with a specified asymptote

            # Linearised first
            p0, x, y, resid = linfit_moft(times, cmeans, asymp)
            plt.plot(x, y, ':', lw=2, label='linear')
            print(p0, resid)

            # Now nonlinear
            p, x, y, resid = geofit_moft(times, cmeans, asymp, p0)
            plt.plot(x, y, '--', lw=2, c='k', label='geometric')
            print(p, resid)

    plt.xlabel(r'time, $t$')
    plt.ylabel(r'mean cluster size, $\bar{k}(t)$')
    plt.legend()
    plt.xlim(adjust_range(plt.xlim(), args.xrange))
    plt.ylim(adjust_range(plt.ylim(), args.yrange))


    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

mpoft = lambda args: mass_oft(args, 'pdist')
mvoft = lambda args: mass_oft(args, 'vdist')

def rhotot_oft(args, dist):

    cds = cdists_timeseries(process_files(args.file), dist, normed=False)

    for P, hists in cds:
        times = np.array(sorted(list(hists.keys())))

        rhotots = np.zeros(len(times))
        negs = np.zeros(len(times))
        for n,t in enumerate(times):
            hist = hists[t]
            rhotots[n] = np.sum(hist)
            negs[n] = np.argmax(hist<0)

        t, r = prune(times, rhotots)
        # plt.semilogy(t,r)
        # coeff = np.vstack( (np.ones_like(t), np.log(t) ) ).T
        # print(np.linalg.lstsq(coeff, np.log(r), rcond=None))
        t, n = prune(times, negs)
        plt.semilogx(t,n)
        # print(np.linalg.lstsq(coeff, np.log(n), rcond=None))
        # t, z = prune(times, negs**2 * rhotots)
        # plt.plot(t,z)


    plt.xlabel(r'time, $t$')
    plt.ylabel(r'$\rho_{tot}(t)$')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

rtpoft = lambda args: rhotot_oft(args, 'pdist')
rtvoft = lambda args: rhotot_oft(args, 'vdist')

def dmass_oft(args, dist):
    """
    Attempt to find derivative of mass wrt time with poor-man's loess fit
    """

    cds = cdists_timeseries(process_files(args.file), dist)
    for P, hists in cds:
        times = np.array(sorted(list(hists.keys())))

        cmeans = np.zeros(len(times))
        for n,t in enumerate(times):
            hist = hists[t]
            cmeans[n] = np.sum(np.arange(1,len(hist)+1) * hist)

        # times, cmeans = prune(times, cmeans)
        #
        # # TODO: add distinguishing labels as required
        # plt.plot(times, cmeans, 'o' if P.beta is None else ':' if P.beta == 0 else '-')

        # Fit local polynomial through a moving window; plot derivative
        win = 9 # Should be about right for a quadratic
        fxs = []
        fys = []
        for n in range(len(cmeans)-win):
            xs, ys = times[n:n+win], cmeans[n:n+win]
            xm, ym = xs[win//2], ys[win//2]
            coeff = np.vstack(( np.ones_like(xs), xs-xm, (xs-xm)**2 )).T
            res = np.linalg.lstsq(coeff, ys-ym, rcond=None)[0]
            fxs += [xm]
            fys += [res[1]]

        fxs, fys = np.array(fxs), np.array(fys)
        plt.loglog(fxs, fys*(fxs)*(np.log(fxs)**2), 'o' if P.beta is None else ':' if P.beta == 0 else '-')

    plt.xlabel(r'time, $t$')
    plt.ylabel(r'cluster growth rate, $\dot{k}(t)$')
    plt.legend()
    plt.xlim(adjust_range(plt.xlim(), args.xrange))
    plt.ylim(adjust_range(plt.ylim(), args.yrange))


    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

dmpoft = lambda args: dmass_oft(args, 'pdist')
dmvoft = lambda args: dmass_oft(args, 'vdist')


def residuals_moft(args, dist):
    """
    Residuals of the linear/nonlinear fits as a function of asymptote
    """

    cds = cdists_timeseries(process_files(args.file), dist)
    for P, hists in cds:
        times = np.array(sorted(list(hists.keys())))

        cmeans = np.zeros(len(times))
        for n,t in enumerate(times):
            hist = hists[t]
            cmeans[n] = np.sum(np.arange(1,len(hist)+1) * hist) / np.sum(hist)

        times, cmeans = prune(times, cmeans)

        x = np.array(times)
        y = cmeans[x>1]
        x = x[x>1]

        rates = np.linspace(0,0.5)
        lresids = np.zeros_like(rates)
        nresids = np.zeros_like(rates)

        for n,rate in enumerate(rates):

            # Get residuals of the linear fit
            asymp = np.log if rate==0 else lambda u:u**rate
            p0, x, y, lresids[n] = linfit_moft(times, cmeans, asymp)

            # Get residuals of the nonlinear fit
            p, x, y, nresids[n] = geofit_moft(times, cmeans, asymp, p0)

        plt.axhline(min(lresids), linestyle='--', color='0.7')
        plt.axvline(rates[np.argmin(lresids)], linestyle='--', color='0.7')
        plt.plot(rates, lresids, 'x', label='linear')
        plt.axhline(min(nresids), linestyle=':', color='0.7')
        plt.axvline(rates[np.argmin(nresids)], linestyle=':', color='0.7')
        plt.plot(rates, nresids, '+', label='geometric')
        plt.legend()
        plt.xlabel(r'Approach to stationarity, $\nu$')
        plt.ylabel(r'Goodness of fit, $|\vec{k}-\vec{k}_{fit}|$')

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)


residp = lambda args: residuals_moft(args, 'pdist')
residv = lambda args: residuals_moft(args, 'vdist')


def cluster_movie(args, dist):
    """
    Animation of cluster size distributions
    """

    normed = not args.raw
    cds = cdists_timeseries(process_files(args.file), dist, normed=normed)

    # Establish a reasonable set of times at which all series are defined; we assume
    # regularly spaced measurements
    t0,t1,interval = None,None,None # t0 is going to be the latest start, t1 the earliest finish; interval the longest sampling interval
    stationary = [] # Place to keep stationary distributions
    for P, hists in cds:
        # Extract stationary distributions
        if None in hists:
            stationary += [(P, hists[None])]
            del hists[None]
        if not hists: continue
        mint, maxt = min(hists.keys()), max(hists.keys())
        if t0 is None or mint > t0:
            t0 = mint
        if t1 is None or maxt < t1:
            t1 = maxt
        dt = (maxt-mint)/(len(hists.keys())-1)
        if interval is None or dt > interval:
            interval = dt

    if t1 < t0:
        print("Time series did not overlap, no movie produced")
        return

    # Get the set of times that we want histograms at
    times = np.arange(t0, t1+interval, interval)
    histseqs = []
    paramseq = []

    # Extract the histograms that lie closest to the desired times
    for P, hists in cds:
        if not hists: continue
        paramseq += [P]
        histseqs += [ np.array(list(hists.values()),dtype=object)[np.argmin(abs(times-np.array(list(hists.keys()))[:,np.newaxis]),axis=0)] ]

    del(cds) # This will release any histograms we won't be plotting

    # Get maximum k and histogram heights
    max_k = max( len(hist) for histseq in histseqs for hist in histseq )
    max_h = max( max(hist) for histseq in histseqs for hist in histseq )
    min_h = min( 0, min(min(hist) for histseq in histseqs for hist in histseq) )
    # First set up the figure, the axis, and the plot element we want to animate
    fig, ax = plt.subplots(1,1)
    plots = [ ax.plot([], [], 'o' if P.beta is None else ':' if P.beta == 0 else '-')[0] for P in paramseq ]

    # Plot fixed asymptotes
    if args.raw:
        plt.axhline(ls='--',c='C1')
    else:
        for alpha, rho in set( (P.alpha, P.density) for P in paramseq ):
            for P, hist in stationary:
                if P.alpha == alpha and P.density == rho:
                    ks = np.arange(1,len(hist)+1)
                    pks = ks[hist>0] # Cluster size
                    nhist = hist[hist>0] # Corresponding number
                    Z = np.sum(nhist) # Total number
                    if normed:
                        nhist = nhist/Z
                    # TODO: identfy in some way?
                    plt.plot(pks, nhist, '--')
                    break
            else:
                # Default exponential histogram from single-particle theory
                # NB: we could deprecate this and require an explicit stationary
                # distribution to be provided?
                b = 0.5*alpha**2
                r = (np.sqrt(b**2+4*b*(alpha-b)*rho*(1-rho))-b) / 2/(alpha-b)
                g = b/alpha * (1-rho-r)/r
                pre = g if normed else g*g*rho
                ks = np.arange(1,max_k+1)
                plt.plot(ks, pre*(1+g)**(1-ks), '--')

    txt = ax.text(0.5,0.95,'',ha='center',transform=ax.transAxes)
    ax.set_xlim(adjust_range((0,max_k), args.xrange))
    ax.set_ylim(adjust_range((min_h,max_h), args.yrange))
    plt.xlabel(r'Cluster size, $k$')
    if args.raw:
        plt.ylabel(r'Correction, $u(k)$')
    else:
        plt.ylabel(r'Frequency, $P(k)$')

    def update(n):
        for histseq, plot in zip(histseqs, plots):
            hist = histseq[n]
            ks = np.arange(1,len(hist)+1)
            if args.raw:
                d = np.sum(hist)
                hist = hist/d

            plot.set_data(ks, hist)
        txt.set_text(f"t={times[n]}")
        return plots + [txt]

    ani = FuncAnimation(fig, update, frames=len(times), interval=50, repeat=True, blit=True, save_count=sys.maxsize)

    if args.fig is None:
        plt.show()
    else:
        ani.save(args.fig)

pmov = lambda args: cluster_movie(args, 'pdist')
vmov = lambda args: cluster_movie(args, 'vdist')
