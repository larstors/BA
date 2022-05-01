#!/usr/bin/env python3

"""
Visualise the crumble output
"""

import importlib, argparse

# Available visualisation types, and where to find them
# Map from type to module; will call module.type(args) to generate output
visualisations = {
    'inefficiency': 'performance', 'gain': 'performance',
    'arrows': 'movies', 'blocks': 'movies',
    'nndist': 'nearest', 'msd': 'displacement',
    'pdist': 'clusters', 'vdist': 'clusters', 'pvdist': 'clusters', 'vpdist': 'clusters',
    'pmean': 'clusters', 'vmean': 'clusters',
    'poft': 'clusters', 'voft': 'clusters', 'mpoft': 'clusters', 'mvoft': 'clusters',
    'rtpoft': 'clusters', 'rtvoft': 'clusters', 'dmpoft': 'clusters', 'dmvoft': 'clusters', 'claw': 'clusters',
    'residp':'clusters', 'residv':'clusters', 'pmov': 'clusters', 'vmov': 'clusters',
    'pcomp': 'clustersbymass', 'vcomp': 'clustersbymass', 'ptom': 'clustersbymass', 'vtom': 'clustersbymass',
    'pmmov': 'clustersbymass', 'vmmov': 'clustersbymass', 'ptot': 'clustersbymass', 'vtot': 'clustersbymass'
}

# TODO: could move to subcommands rather than have the quasi-mandatory -t option
# This would also allow different options for the different subcommands

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--active', '-a', action='store_true', help='(movies only) distinguish active particles')
parser.add_argument('--times', '-T', help='(poft/voft only) timepoint(s)')
parser.add_argument('--means', '-M', help='(pofm/vofm only) mean cluster size(s)')
parser.add_argument('--mstep', '-S', default=1.0, type=float, help='(pmmov/vmmov only) mean cluster size increase between frames')
parser.add_argument('--xrange', '-x', help='(where supported) x axis range')
parser.add_argument('--yrange', '-y', help='(where supported) y axis range')
parser.add_argument('--raw', '-R', action='store_true', help='(where supported) no renormalisation of cluster distribution')
parser.add_argument('--asymp', '-A', default=0, type=float, help='(mpoft/mvoft) assumed exponent of approach to stationarity (0=log)')
parser.add_argument('--type', '-t', choices=visualisations.keys(), default='inefficiency', help='Plot type: ' + ('|'.join(visualisations.keys())) + ' [inefficiency]')
parser.add_argument('--fig', '-F', help='Figure filename [screen]')
parser.add_argument('file', default=['-'], nargs='*', help='File to read data from [stdin]')
args = parser.parse_args()

# Import the relevant module and call the function with the same name as type
getattr(importlib.import_module('dovis.'+visualisations[args.type]), args.type)(args)
