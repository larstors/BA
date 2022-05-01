#!/usr/bin/env python3

"""
Exact solution of the equation d rho_n / d t = alpha [ rho_{n-1} - 2 rho_n + rho_{n+1} ]
"""

import argparse
import numpy as np
from scipy.special import ive

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--length",   "-L", type=int,   default=1000,  help="Lattice pseudo-length")
parser.add_argument("--density",  "-r", type=float, default=0.3,   help="Particle density")
parser.add_argument("--tumble",   "-t", type=float, default=0.001, help="Tumble rate")
parser.add_argument("--truncate", "-n", type=int,   default=100,   help="Cluster size truncation")
parser.add_argument("--bessel",   "-k", type=int,   default=100,   help="Bessel function truncation")
parser.add_argument("--until",    "-u", type=float, default=1000,  help="Time to run for once measurements started")
parser.add_argument("--every",    "-e", type=float, default=50,    help="Measurement interval")

args = parser.parse_args()


def rho(n,t,phi,kmax):
    """
    Cluster density at sites n and times t starting from an initial condition
    with particle volume fraction phi.

    Returns an array where element [i,j] corresponds to site n[i] at time t[j] / 2 alpha.
    """

    # Initial condition
    ks = np.arange(kmax)[:,np.newaxis]
    rho0 = (1-phi)**2 * phi**(1+ks)

    # Density at site n
    return np.sum(rho0[...,np.newaxis] * (ive((n-1-ks)[...,np.newaxis], t) - ive((n+1+ks)[...,np.newaxis], t)), axis=0)

sol = rho(np.arange(1,args.truncate), 2*args.tumble*np.arange(0,args.until+args.every,args.every), args.density, args.bessel)

# Produce output compatible with C++ code
print(f"""\
# L = [ {args.length} ]
# N = {int(args.length*args.density)}
# alpha = {args.tumble}
# beta = 0
# output = mfiter
# initial = 0
# interval = {args.every}\
""")
for dist in sol.T:
    print(' '.join(map(str,dist)))
