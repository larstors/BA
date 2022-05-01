#!/usr/bin/env python3

"""
Determine the dense effective initial condition
"""

from numpy.random import binomial
import numpy as np, matplotlib.pyplot as plt
import random, re, argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--length",    "-L", type=int,   default=1000,  help="Lattice pseudo-length")
parser.add_argument("--density",   "-r", type=float, default=0.3,   help="Particle density")
parser.add_argument("--samples",   "-s", type=int,   default=1000,  help="Number of samples for histogram")
parser.add_argument("--jamtime",   "-j", type=float, default=50.0,  help="Jamming onset time (simulation data only)")
parser.add_argument("--xrange",    "-x",                            help="X-axis range")
parser.add_argument('file', nargs='*', help='File to read data from, where required [use - for stdin]')
args = parser.parse_args()

phi = args.density # quicker to type

def domainsizes(L, N, all=True):
    """
    Return the domain sizes for a set of N particles on an L site lattice
    """

    Nr = binomial(N, 0.5)

    # Actual initial condition in string representation so we can do re replaces
    init = ['+']*Nr + ['-']*(N-Nr) + ['.']*(L-N)
    random.shuffle(init)
    init = ''.join(init)

    if Nr == 0:
        final = '-'*N
    elif Nr == N:
        final = '+'*N
    else:
        # Rotate the configuration so that the first particle is + and the last is -
        # (This will always be possible since 0<Nr<N)
        # Also strip any vacancies before the first + and after the last -
        if re.match(r'^\.*\+', init) is None or re.search(r'-\.*$', init) is None:
            # We need to rotate
            mat = re.search(r'-\.*\+', init)
            final = init[mat.end()-1:] + init[:mat.start()+1]
        else:
            final = init

        # Now close up all + - pairs
        final = re.sub(r'\+\.{2,}-', '.+-.', final)
        final = re.sub(r'\+\.-', lambda m: '.+-' if random.random() < 0.5 else '+-.', final)

        # Now close up ++ and -- pairs, leaving a space behind the particle that catches up
        while re.search(r'\+\.+\+', final):
            final = re.sub(r'\+\.+\+', '.++', final, count=1)

        while re.search(r'-\.+-', final):
            final = re.sub(r'-\.+-', '--.', final, count=1)

        final = re.sub(r'^\.+', '', final)
        final = re.sub(r'\.+$', '', final)

    return [len(x) for x in re.split(r'\.+', final) if all or re.fullmatch(r'\++-+',x)]

# Build up a cluster size distribution
hist = np.array([],dtype=int)
for n in range(args.samples):
    dhist = np.bincount(domainsizes(args.length, int(args.length*phi), False))
    if len(hist) >= len(dhist):
        hist[:len(dhist)] += dhist
    else:
        dhist[:len(hist)] += hist
        hist = dhist

xs = np.arange(1,len(hist))
yes = phi * hist[1:] / np.sum(xs*hist[1:])

yts = phi * 0.5**(xs+2) * (xs-1)  # Low density only

print(yes[1]) # dimers
print(phi*(1-phi)/16) # * (1 + phi*(1-phi)*(2+phi)**2/(2-phi)**2) )

# plt.plot(xs, yes, 'x')
# plt.plot(xs, yts, ':')
# plt.plot(xs, yrs)

# Add plots from any files that are given. We will grab the second
# cluster distribution from each file as an estimate of the effective distribution
if args.file:
    from dovis.core import *
    find_jam = False
    for S in process_files(args.file):
        if isinstance(S, Parameters):
            find_jam = True
        elif find_jam and isinstance(S, ClusterDistribution):
            if S.time >= args.jamtime:
                xss = 1+np.arange(len(S.pdist))
                yss = phi * S.pdist / np.sum(xss*S.pdist)
                print(yss[1])
                # plt.plot(xss, yss, 'o')
                find_jam = False
    # plt.xlim(adjust_range(plt.xlim(), args.xrange))

# plt.show()
