"""
Plot of how well the simulation algorithm performs relative to the naive
implementation
"""

from dovis.core import *
import matplotlib.pyplot as plt

def performance_plot(args):

    # Collate data from file or files
    statistics = {}
    snapshots, svalue = 0, 0.0

    for obj in process_files(args.file):
        if isinstance(obj, Parameters):
            if snapshots:
                statistics[(P.density, P.alpha)] = svalue/snapshots
            P = obj # currently-applicable parameters
            snapshots, svalue = 0, 0.0
        elif isinstance(obj, FullSnapshot):
            snapshots += 1
            svalue += (obj.inefficiency if args.type == 'inefficiency' else obj.gain)

    if snapshots:
        statistics[(P.density, P.alpha)] = svalue/snapshots

    densities = sorted(set( x[0] for x in statistics.keys() ))
    alphas = sorted(set( x[1] for x in statistics.keys() ))

    for alpha in alphas:
        xs = []
        ys = []
        for density in densities:
            if (density, alpha) not in statistics: continue
            xs += [density]
            ys += [statistics[(density,alpha)]]
        plt.plot(xs, ys, label=rf'$\alpha={alpha}$')

    plt.legend()

    if args.fig is None:
        plt.show()
    else:
        plt.savefig(args.fig)

# Both the inefficiency and gain types map to performance_plot
inefficiency = gain = performance_plot
