"""
Distance to nearest neighbours
"""

from dovis.core import *
import numpy as np, matplotlib.pyplot as plt
from collections import Counter

def nndist(args):
    # Distribution of distances between nearest neighbours (particles or vacancies but not both)
    data = process_files(args.file)
    P = next(data)
    V = P.volume # Total number of sites on the lattice
    R = Ruler(P.shape)

    # Smallest circle of integer radius centred on x that the point y falls inside (including boundary)
    cdist = lambda x,y: int(np.ceil( np.sqrt(np.sum(R(x,y)**2)) ))

    # Count the number of times a given cdist is represented on the lattice
    O = np.array(np.unravel_index(0, P.shape)) # Location of the origin
    count = Counter( cdist(O, np.array(np.unravel_index(n, P.shape))) for n in range(1,V) )

    maxdist = max(count.keys())+1 # upper bound on the nn dist

    # Set up an empirical distribution
    empirical = np.zeros(maxdist, dtype=int)

    valid_snapshots = (ParticleSnapshot, VacancySnapshot)
    N = None # Number of particles / vacancies we are dealing with

    ## Build the empirical distribution
    for S in data:
        if isinstance(S, Parameters):
            if sum(empirical) > 0:
                plt.plot(np.arange(1,maxdist),empirical[1:]/np.sum(empirical), 'o', label=rf'$\alpha={P.alpha}$')
            if S.shape != P.shape or S.N != P.N:
                print("Warning: system size or density has changed; bailing out")
                break
            P = S
            empirical[:] = 0

        elif isinstance(S, valid_snapshots):

            coords = S.coordinates

            if N is None:
                N = len(coords)
                valid_snapshots = type(S)

            # We shuffle the order of the coordinates to ensure a different focal particle each time
            np.random.shuffle(coords)
            # The smallest cdist is the one we want
            empirical[min( cdist(coords[0], coords[n]) for n in range(1,N) )] += 1
        else:
            raise RuntimeError("Require compact snapshots of the same type for nndist analysis")

    if sum(empirical) > 0:
        plt.plot(np.arange(1,maxdist),empirical[1:]/np.sum(empirical), 'o', label=rf'$\alpha={P.alpha}$')

    ## Obtain the baseline distribution

    # Convert to a cumulative distribution vs[n] = sum(vs[r] for r<=n)
    # Also get the probability that all these sites are empty with a given number of particles / vacancies
    v = 0
    empty = np.ones(maxdist)
    for r in range(1, maxdist):
        v += count[r]
        empty[r] = np.prod( np.arange(V-N-v+1,V-N+1) / np.arange(V-v,V) )

    # Convert to distribution over cdists
    gas = np.roll(empty,1) - empty
    plt.plot(np.arange(1,maxdist),gas[1:], ':', label='gas')

    plt.xlabel(r'$k$')
    if valid_snapshots is VacancySnapshot:
        plt.ylabel(r'$P_v(k)$')
    else:
        plt.ylabel(r'$P_p(k)$')
    plt.legend()

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)
