"""
Track mean square displacement
"""

from dovis.core import *
import numpy as np, matplotlib.pyplot as plt

def msd(args):
    # Track the mean square displacement of particles/vacancies as a function of time
    # Can average over multiple runs with the same parameter values

    data = process_files(args.file)
    P = next(data)
    R = Ruler(P.shape)
    t0 = None # Initial time; None indicates that we are waiting for a new timeseries
    nsd = {} # Number of square displacement recorded at the time point given by the key
    ssde = {} # Sum of square end-to-end position on periodic lattice at time points matching those in nsd
    ssdi = {} # Sum of square displacements accumulated between snapshots (can grow indefinitely) at point points matching those in nsd

    valid_snapshots = (VacancySnapshot, ParticleSnapshot)

    for S in data:
        if isinstance(S, Parameters):
            # Check that any subsequent datasets are run at the same parameter values
            # (We could subdivide into different series; and reset valid_snapshots too)
            if P.shape != S.shape or P.N != S.N or P.alpha != S.alpha:
                raise RuntimeError("Can't (yet) combine data from runs at different parameter values")
            # Await new timeseries
            t0 = None
        elif isinstance(S, valid_snapshots):
            valid_snapshots = type(S) # From now on, we only allow snapshots of the same type
            if t0 is None:
                x0, x, t0 = S.coordinates, S.coordinates, S.time
                dxi = np.zeros_like(x0, dtype=float) # Total displacement, ignoring periodicity
            else:
                dt = S.time-t0 # Total time elapsed since initial condition
                dx = R(x0, S.coordinates) # End-to-end distance respecting peridicity
                dxi += R(x, S.coordinates) # Accumulated total displacement between snapshots
                x = S.coordinates # Keep track of last snapshot
                if dt in nsd:
                    nsd[dt]  += len(dx)
                    ssde[dt] += np.sum(dx**2)
                    ssdi[dt] += np.sum(dxi**2)
                else:
                    nsd[dt]  = len(dx)
                    ssde[dt] = np.sum(dx**2)
                    ssdi[dt] = np.sum(dxi**2)
        else:
            raise RuntimeError("Require compact snapshots of the same type for MSD analysis")

    # Now extract the values at the different time points; we won't aggregate nearby times just yet
    ts = np.array(sorted(nsd.keys()))
    mdse = np.array([ssde[t]/nsd[t] for t in ts])
    mdsi = np.array([ssdi[t]/nsd[t] for t in ts])

    plt.loglog(ts, mdse)
    plt.loglog(ts, mdsi)

    sat = np.argmax(mdsi > 1.1*mdse) # Estimate point where data saturate

    # Do early and late time fits, subject to enough data
    if sat > 10:
        # Early-time fit
        a,b = np.linalg.lstsq( np.vstack( (np.ones_like(ts[:sat]), np.log(ts[:sat]) )).T, np.log(mdsi[:sat]), rcond=None)[0]
        plt.loglog(ts, np.exp(a)*ts**b, 'k--', label=rf'$t^{{{b:.2f}}}$')
    if len(ts)-sat > 10:
        # Late-time fit
        a,b = np.linalg.lstsq( np.vstack( (np.ones_like(ts[sat:]), np.log(ts[sat:]) )).T, np.log(mdsi[sat:]), rcond=None)[0]
        plt.loglog(ts, np.exp(a)*ts**b, 'k:', label=rf'$t^{{{b:.2f}}}$')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$\langle \delta x(t)^2 \rangle$')
    plt.legend()

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)
