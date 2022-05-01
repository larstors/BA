#!/usr/bin/env python3

"""
List files matching specific filename criteria
"""

import argparse, os, re, numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--datadir', '-D', default='data', help='Data directory')
parser.add_argument('--density', '-r', type=float, help='Particle/vacancy density')
parser.add_argument('--vacancies', '-v', action='store_true', help='Select by vacancy properties')
parser.add_argument('--tumble', '-t', help='Tumble rate')
parser.add_argument('--sort', '-s', choices=['V'], help='Sort key' )
parser.add_argument('output', help='Output type')
args = parser.parse_args()

matching = []

for f in os.listdir(args.datadir):
    if not re.search(f'-{args.output}\.dat$', f): continue # Wrong data type
    # Determine the system size from the filename
    mat = re.match(r'L([\dx]+)', f)
    if mat is None: continue # Can't grok system size spec
    L = [ int(x) for x in mat.group(1).split('x') ]
    V = np.prod(L)
    # TODO: apply system size / dimensionality filters here
    if args.density is not None:
        N = int(V * args.density)
        if args.vacancies:
            N = V-N
        if not re.search(f'-N{N}-', f): continue # Wrong number of particles
    if args.tumble is not None:
        if not re.search(f'-t{args.tumble}-', f): continue # Wrong tumble rate

    if args.sort is not None:
        if args.sort == 'V':
            matching += [(f,V)]
    else:
        matching += [f]

if args.sort is not None:
    matching = [ x[0] for x in sorted(matching, key=lambda y: y[1]) ]

for f in matching:
    print(f"{args.datadir}/{f}")
