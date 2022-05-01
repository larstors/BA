"""
Core routines for processing the output from the C++ code
"""

import re, numpy as np, sys, fractions

# These classes provide a python view of different parts of the simulation output

class Parameters:
    "Simulation run parameters"

    def __init__(self, shape, N, alpha, beta, initial, interval, localaverage=None, localinterval=None):
        self.shape = tuple(shape)
        self.N = N
        self.alpha = alpha
        self.beta = beta
        self.initial = initial
        self.interval = interval
        self.localaverage = localaverage
        self.localinterval = localinterval

        # self.dirs[d] holds the vector that an agent with direction d = 1,2,..., 2 len(shape)
        # moves in
        self.dirs = np.zeros((2*len(shape)+1,len(shape)))
        for d in range(1,2*len(shape)+1):
            self.dirs[d][len(shape)-(d+1)//2] = 1 if d%2 else -1

    @property
    def strshape(self):
        "System size as a string (e.g., 100x100)"
        return 'x'.join( [str(x) for x in self.shape] )

    @property
    def volume(self):
        "Total number of sites on the lattice"
        return np.prod(self.shape)

    @property
    def density(self):
        "Particle density"
        return fractions.Fraction(int(self.N),int(self.volume))

    @property
    def vdensity(self):
        "Vacancy density"
        return 1-self.density

    def __str__(self):
        return f"{self.shape=} {self.N=} {self.alpha=}"


class FullSnapshot:
    "Full information about particles, directions and active status"

    def __init__(self, time, directions, active):
        self.time = time
        self.directions = directions

        self.active = active.astype(bool) # Active sites
        self.occupied = self.directions > 0 # Occupied sites

        # Count neigbours
        self.mobile = np.zeros_like(self.occupied, dtype=int)
        for axis in range(len(self.occupied.shape)):
            self.mobile += np.roll(self.occupied, 1, axis=axis)
            self.mobile += np.roll(self.occupied, -1, axis=axis)
        # Sites where a particle can move (at least one vacant nb)
        self.mobile = (self.mobile < 2*len(self.occupied.shape)) & self.occupied

    @property
    def inefficiency(self):
        return np.sum(self.active) / np.sum(self.mobile)

    @property
    def gain(self):
        return np.sum(self.occupied) / np.sum(self.active)


class VacancySnapshot:
    "Limited information about vacancy (and hence particle) location"

    def __init__(self, time, shape, ids, locations):
        self.time = time
        self.ids = ids
        self.locations = locations

        # Construct occupied array as we would have for full snapshot
        self.occupied = np.ones(np.prod(shape), dtype=bool)
        self.occupied[locations] = False
        self.occupied = self.occupied.reshape(shape)

    @property
    def coordinates(self):
        "Element [n,d] of the returned array contains the dth coordinate of vacancy n"
        return np.array(np.unravel_index(self.locations[np.argsort(self.ids)], self.occupied.shape)).T


class ParticleSnapshot:
    "Limited information about particle (and hence vacancy) location"

    def __init__(self, time, shape, ids, locations):
        self.time = time
        self.ids = ids
        self.locations = locations

        # Construct occupied array as we would have for full snapshot
        self.occupied = np.zeros(np.prod(shape), dtype=bool)
        self.occupied[locations] = True
        self.occupied = self.occupied.reshape(shape)

    @property
    def coordinates(self):
        "Element [n,d] of the returned array contains the dth coordinate of particle n"
        return np.array(np.unravel_index(self.locations[np.argsort(self.ids)], self.occupied.shape)).T


class ClusterDistribution:
    "Distribution of particle and vacancy cluster sizes"

    # pdist[n] is number of particle clusters of size n+1
    # vdist[n] is number of vacancy clusters of size n+1
    # Both arrays are of the same length at a given timepoint, but may vary between timepoints

    def __init__(self, time, pdist, vdist):
        self.time = time
        self.pdist = pdist
        self.vdist = vdist

    @property
    def pmean(self):
        "Mean particle cluster size"
        return np.sum(self.pdist*(1+np.arange(len(self.pdist)))/np.sum(self.pdist))

    @property
    def vmean(self):
        "Mean vacancy cluster size"
        return np.sum(self.vdist*(1+np.arange(len(self.vdist)))/np.sum(self.vdist))

def process(lines):
    """
    Iterate over a sequence of lines and spit out Parameters and Snapshots
    as we have them. Any given Snapshot corresponds to the most recent Parameters.
    """
    pupdate = False # True when Parameters have been updated
    for line in lines:
        if line[0] == '#':
            if pupdate == False:
                shape = N = alpha = beta = output = initial = interval = localaverage = localinterval = None
            # In header section; collate parameters
            if mat := re.search(r'\bL = \[ (((\d+) )+)\]', line):
                shape = [int(x) for x in reversed(mat.group(1).split())]
                pupdate = True
            if mat := re.search(r'\bN = (\d+)', line):
                N = int(mat.group(1))
                pupdate = True
            if mat := re.search(r'\balpha = ([\d.]+)', line):
                # Isotropic tumbling
                alpha = float(mat.group(1))
                pupdate = True
            if mat := re.search(r'\balpha = \[ ((([\d.]+) )+)\]', line):
                # Anisotropic tumbling: for now, just average the tumble rate
                alpha = np.mean( [float(x) for x in mat.group(1).split()] )
                pupdate = True
            if mat := re.search(r'\bbeta = ([\d.e+-]+)', line):
                # Effective dimer creation rate (theory only)
                beta = float(mat.group(1))
                pupdate = True
            if mat := re.search(r'\boutput = (\w+)', line):
                output = mat.group(1)
                pupdate = True
            if mat := re.search(r'\binitial = ([\d.]+)', line):
                initial = float(mat.group(1))
                pupdate = True
            if mat := re.search(r'\binterval = ([\d.]+)', line):
                interval = float(mat.group(1))
                pupdate = True
            if mat := re.search(r'\blocalaverage = (\d+)', line):
                localaverage = int(mat.group(1))
                pupdate = True
            if mat := re.search(r'\blocalinterval = ([\d.]+)', line):
                localinterval = float(mat.group(1))
                pupdate = True
        else:
            # In snapshot section; yield parameters if any are pending
            if pupdate:
                yield Parameters(shape, N, alpha, beta, initial, interval, localaverage, localinterval)
                pupdate = False
                sid = 0 # Where we are in a sequence of snapshots
            if output == 'snapshots':
                # Get combined directions, active pairs
                das = np.array([int(x) for x in line.split()])
                yield FullSnapshot(initial + sid*interval, das[0::2].reshape(shape), das[1::2].reshape(shape))
                sid += 1
            elif output == 'vacancies':
                # Get vacancy id, site pairs
                vss = np.array([int(x) for x in line.split()])
                yield VacancySnapshot(initial + sid*interval, shape, vss[0::2], vss[1::2])
                sid += 1
            elif output == 'particles':
                # Get particle id, site pairs
                pss = np.array([int(x) for x in line.split()])
                yield ParticleSnapshot(initial + sid*interval, shape, pss[0::2], pss[1::2])
                sid += 1
            elif output == 'clusters':
                # Get particle, vacancy distribution
                pvd = np.array([int(x) for x in line.split()])
                yield ClusterDistribution(initial + sid*interval, pvd[0::2], pvd[1::2])
                sid += 1
            elif output == 'mfiter':
                # Get particle cluster size distribution (fractional)
                # This output type does not have a vacancy cluster size distribution
                pd = np.array([float(x) for x in line.split()])
                yield ClusterDistribution(initial + sid*interval, pd, None)
                sid += 1
            elif output == 'stationary':
                # Same as mfiter, but with no time associated to the histogram
                pd = np.array([float(x) for x in line.split()])
                yield ClusterDistribution(None, pd, None)

    # Return any dangling parameters
    if pupdate:
        yield Parameters(shape, N, alpha, beta, initial, interval, localaverage, localinterval)

def process_files(filenames):
    """
    Process all files in the sequence of filenames provided
    """
    for file in filenames:
        with (sys.stdin if file == '-' else open(file)) as fh:
            yield from process(fh)

def process_block(filenames, skip=0):
    """
    Return a single block from the sequence of filenames. A block is
    defined as a Parameters object followed by a sequence of non-Parameters
    objects. The specified number of blocks is skipped.
    """
    block = -1
    for file in filenames:
        with (sys.stdin if file == '-' else open(file)) as fh:
            for obj in process(fh):
                if isinstance(obj, Parameters):
                    block += 1
                    if block > skip: return
                if block == skip:
                    yield obj


class Ruler:
    "Measure the shortest distance between two points on the lattice, respecting periodicity"

    def __init__(self, shape):
        self.L = np.array(shape)

    def __call__(self, r0, r1):
        """
        Returns the vector r1-r0 corresponding to the shortest distance between them
        NB: r0 and r1 can be arrays as long as the final dimension matches that of the
        lattice.
        """
        return self.L/2 - (self.L/2-(r1-r0))%self.L


def adjust_range(r, spec=None):
    """
    Take a tuple r=(a,b) and a specification string "[c:]d"
    and return (max(a,c), min(b,d)).
    If spec is None, r is returned. This makes it easy to use with command line
    arguments
    """
    if spec is not None:
        s = ([None] + spec.split(':'))[-2:]
        if s[0]: r = (float(s[0]), r[1])
        if s[1]: r = (r[0], float(s[1]))
    return r


def prune(*xs, points=100):
    """
    Reduce each array in xs to the desired number of points by averaging

    Returns an iterable over the pruned arrays
    """
    for x in xs:
        n = len(x)
        if n>points:
            k,r = divmod(n, points)
            yield np.mean(x[:n-r].reshape(points,k), axis=1)
        else:
            yield x
