#!/usr/bin/env python3

"""
Python implementation of mean-field equations for cluster distribution
"""

import numpy as np, argparse
from scipy.integrate import solve_ivp
from dovis.core import *

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--length",    "-L", type=int,   default=1000,  help="Lattice pseudo-length")
parser.add_argument("--density",   "-r", type=float, default=0.3,   help="Particle density")
parser.add_argument("--tumble",    "-t", type=float, default=0.001, help="Tumble rate")
parser.add_argument('--initial',   "-i", choices=['homogeneous','clustered','file'], default='homogeneous', help="Initial condition: homogeneous|clustered")
parser.add_argument("--source",    "-s", choices=['monomer','effective', 'polymer'], default='monomer', help="Source term: monomer|effective|polymer")
parser.add_argument("--delta",     "-d", action='store_true', help="Solve for distance to stationary state")
parser.add_argument("--kernel",    "-k", choices=['stationary','single'], default='stationary', help="Polymeric chipping kernel: stationary|single")
parser.add_argument("--pairkernel", "-K", choices=[None,'stationary','single'], default=None, help="Pair chipping kernel (if different): stationary|single")
parser.add_argument("--truncate",  "-n", type=int,   default=100,   help="Cluster size truncation")
parser.add_argument("--until",     "-u", type=float, default=1000,  help="Time to run for once measurements started")
parser.add_argument("--every",     "-e", type=float, default=50,    help="Measurement interval")
parser.add_argument('file', default=['-'], nargs='*', help='File to read data from, where required [stdin]')

args = parser.parse_args()

# For convenience
alpha, phi, M = args.tumble, args.density, args.truncate
beta = 0.5*alpha**2

if args.source == 'monomer':
    # These are Peter's original equations which have a monomeric sourcing term
    # Set up a delta function on the first site of an M-1 site lattice (site 1 of the original lattice)
    site1 = np.zeros(M-1)
    site1[0] = 1

    if args.delta:

        # The system to solve, epsilon = rho - rho_s
        # y[0] is epsilon_tot, y[1:] is epsilon/epsilon_tot
        rhost = (np.sqrt(beta**2 + 4*beta*(alpha-beta)*phi*(1-phi))-beta)/2/(alpha-beta)
        src = (1-phi) * (rhost/phi)**2 * (1-rhost/phi)**np.arange(M-1)

        def odesys(t, y):
            ytot = y[0]
            y = np.hstack(([0], y[1:]))
            # Decay rate of ytot
            dytot = (
                - (alpha + beta*((1-phi)/(rhost+ytot)-1)) * y[1]
                + beta * ((1-phi)/phi*rhost/(rhost+ytot) - 1)
            )
            # Create finite derivatives valid for n=1..M-1
            dyplus = (np.roll(y,-1) - y)[1:] # y[n+1]-y[n]
            dyminus = (np.roll(y,1) - y)[1:] # y[n-1]-y[n]

            return np.hstack((
                [dytot * ytot],
                -dytot * y[1:] +
                alpha*(dyplus+dyminus) +
                beta*((1-phi)/(rhost+ytot)-1)*dyplus +
                beta*(src/(rhost+ytot) - site1)
            ))

    else:
        # The system to solve, y is rho
        def odesys(t, y):
            ytot = np.sum(y)
            # Create finite derivatives valid for n=1..M-1
            dyplus = (np.roll(y,-1) - y)[1:] # y[n+1]-y[n]
            dyminus = (np.roll(y,1) - y)[1:] # y[n-1]-y[n]
            return np.hstack(([0], alpha*(dyplus+dyminus) + beta*(1-phi-ytot)*(site1 + dyplus/ytot)))

elif args.source == 'effective':
    # These are Richard's effective equations which are obtained from the stationary kernel
    # and throwing away all the terms we don't like the look of
    decay = (1/2)**np.arange(2,M) * np.arange(1,M-1)

    # The system to solve, y is rho
    def odesys(t, y):
        ytot = np.sum(y)
        # Create finite derivatives valid for n=1..M-1
        dyplus = (np.roll(y,-1) - y)[2:] # y[n+1]-y[n]
        dyminus = (np.roll(y,1) - y)[2:] # y[n-1]-y[n]
        return np.hstack(([0,0], 6*alpha*(dyplus+dyminus) + beta*(1-phi-ytot)*(decay + 4*dyplus/ytot)))


elif args.source == 'polymer':
    # These are Richard's modified equations, with have a polymeric sourcing term.
    #
    # Available chipping kernels:
    # kernel[t](x,y) is the probability that y particles chip off from a domain of size x
    kernel = {
        'stationary': lambda x,y: (1/2)**y * ((y>0)*(y<x-1) + 4*(x==y)*(y>1)),
        'single':     lambda x,y: 1.0*((y==1)*(x>2) + (y==2)*(x==2))
    }

    # Select a kernel combination
    kern = kernel[args.kernel]
    pkern = kern if args.pairkernel is None else kernel[args.pairkernel]

    n,m = np.indices((M,M))
    Amat = (kern(m,m-n) - 2*(n==m))*(n>1)
    Cmat = (2*(pkern(m,m-n) - (n==m)))*(n>1)

    n,m,r = np.indices((M,M,M))
    Bmat = kern(r,n-m)*(n>1)

    n,m,r,k = np.indices((M,M,M,M)) # NB: Potentially expensive for big M
    Dmat = np.sum( (pkern(m,k) * pkern(r,n-k))*(n>1), axis=3) # K axis
    del(n,m,r,k) # Save space

    def odesys(t, y):
        ytot = np.sum(y)
        return alpha * (Amat @ y + Bmat/ytot @ y @ y) + beta*(1-phi-ytot)/ytot * (Cmat @ y + Dmat/ytot @ y @ y)
else:
    print(f"Unrecognised source type: {args.source}")
    exit(1)

# Set up the initial condition
if args.initial == 'homogeneous':
    y = (1-phi)**2 * phi**(np.arange(M)) # Initial condition
elif args.initial == 'clustered':
    y = phi * (np.arange(M)-1)*0.5**(np.arange(M)+2)
elif args.initial == 'file':
    # Read the input until we find a cluster size distribution
    y = np.zeros(M)
    for S in process_files(args.file):
        if isinstance(S, ClusterDistribution):
            k = min(M-1,len(S.pdist))
            y[1:k+1] = S.pdist[:k]
            break
    else:
        print(f"Could not locate an initial condition in input data")
        exit(1)
else:
    print(f"Unrecognised initial condition: {args.initial}")
    exit(1)

y[0] = 0 # Boundary condition
if args.source == 'polymer':
    y[1] = 0 # Remove monomers with polymeric sourcing

if args.delta:
    y[1:] -= phi * (rhost/phi)**2 * (1-rhost/phi)**np.arange(M-1)

    y[0] = np.sum(y)
    y[1:] /= y[0]


sol = solve_ivp(odesys, (0,args.until), y, t_eval=np.arange(0,args.until+args.every,args.every))

# Produce output compatible with C++ code
print(f"""\
# L = [ {args.length} ]
# N = {int(args.length*phi)}
# alpha = {alpha}
# beta = {beta}
# output = mfiter
# initial = 0
# interval = {args.every}\
""")
if args.delta:
    for dist in sol.y.T:
        print(' '.join(map(str,dist[0]*dist[1:])))
else:
    for dist in sol.y.T:
        print(' '.join(map(str,dist[1:])))
