#!/usr/bin/env python3
"""
Attempt to find the stationary state
"""

import numpy as np, argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--length",    "-L", type=int,   default=1000,  help="Lattice pseudo-length")
parser.add_argument("--density",   "-r", type=float, default=0.3,   help="Particle density")
parser.add_argument("--tumble",    "-t", type=float, default=0.001, help="Tumble rate")
#parser.add_argument('--initial',   "-i", choices=['homogeneous','clustered'], default='homogeneous', help="Initial condition: homogeneous|clustered")
parser.add_argument("--kernel",    "-k", choices=['stationary','single'], default='stationary', help="Polymeric chipping kernel: stationary|single")
parser.add_argument("--pairkernel", "-K", choices=[None,'stationary','single'], default=None, help="Pair chipping kernel (if different): stationary|single")
parser.add_argument("--truncate",  "-n", type=int,   default=100,   help="Cluster size truncation")
parser.add_argument("--tol",       "-T", type=float, default=1e-6,  help='Convergence tolerance')

args = parser.parse_args()


alpha, phi, M = args.tumble, args.density, args.truncate
beta = 0.5*alpha**2

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

# Initial guess
y = phi * (np.arange(M)-1)*0.5**(np.arange(M)+2)
y[:2] = 0
ytot = np.sum(y)

# Iterations
while True:

    # Linearisation
    L = alpha * (Amat + Bmat/ytot @ y) + beta*(1-phi-ytot)/ytot * (Cmat + Dmat/ytot @ y)
    # Remove everything with n<2
    L = L[2:,2:]
    # Replace the bottom row with mass conservation constraint
    L[-1,:] = np.arange(2,M)

    rhs = np.zeros(M-2)
    rhs[-1] = args.density

    # Find the solution of the linearised equations
    ynew = np.hstack(([0,0],np.linalg.solve(L,rhs)))
    ytotnew = np.sum(ynew)

    # Determine how much we have moved
    dy = max(abs(ynew-y))/max(ynew)
    dytot = abs(ytotnew-ytot)/ytotnew

    # Advance to the next iteration
    y, ytot = ynew, ytotnew
    if dy < args.tol and dytot < args.tol:
        break

# Output in the common format
print(f"""\
# L = [ {args.length} ]
# N = {int(args.length*phi)}
# alpha = {alpha}
# beta = {beta}
# output = stationary\
""")
print(' '.join(map(str,y[1:])))
