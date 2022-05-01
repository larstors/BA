#!/usr/bin/env python3

"""
Accumulate cluster distribution time series from different runs
"""

import argparse, numpy as np
import dovis.core as core

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('file', default=['-'], nargs='*', help='File to read data from [stdin]')
args = parser.parse_args()

sets = {
    'pdist':{}, 'vdist':{}
}

P = None

for S in core.process_files(args.file):
    if isinstance(S, core.Parameters):
        if P is not None:
            # Check commensurate
            if S.shape != P.shape or S.N != P.N or P.alpha != S.alpha:
                print("Parameters don't match, skipping")
                continue
        P = S
    elif isinstance(S, core.ClusterDistribution):
        for dist, hists in sets.items():
            hist = np.array(getattr(S, dist)) # Copy the histogram
            if S.time in hists:
                # Make sure the longer of the histograms is stored
                if len(hist) > len(hists[S.time]):
                    hist, hists[S.time] = hists[S.time], hist
                # Add the shorter histogram
                hists[S.time][:len(hist)] += hist
            else:
                hists[S.time] = hist

times = sorted(set(sets['pdist']) | set(sets['vdist']))

# Replicate the file header - nb, these may not be super-accurate
print(f"""\
# L = [ {' '.join(str(x) for x in P.shape)} ]
# N = { P.N }
# alpha = { P.alpha }
# output = clusters
# initial = { times[0] }
# interval = { (times[-1]-times[0])/(len(times)-1) }
# localaverage = { P.localaverage }
# localinterval = { P.localinterval }\
""")

# Now zip up the histograms
for time in times:
    pdist = sets['pdist'][time] if time in sets['pdist'] else []
    vdist = sets['vdist'][time] if time in sets['vdist'] else []
    cdist = np.zeros(2*max(len(pdist), len(vdist)), dtype=int)
    cdist[0:2*len(pdist):2] = pdist
    cdist[1:2*len(vdist)+1:2] = vdist
    print(' '.join([str(x) for x in cdist]))
