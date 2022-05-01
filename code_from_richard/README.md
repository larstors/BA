# Condensed Run-and-tUMBLE (crumble) simulation / visualisation code

Original author: Richard Blythe <R.A.Blythe@ed.ac.uk>

Developments and improvements are encouraged under the MIT license.

Note that the command-line parsing code CLI11.hpp is distributed under its
own license.


## Prerequisites

- a C++17 compiler to build the simulation, e.g. [clang++](https://clang.llvm.org) or [g++](https://gcc.gnu.org)
- [Python 3.8](https://www.python.org) or later for the visualisation code
- [matplotlib 3.3.4](https://matplotlib.org) + dependencies (earlier versions may work)
- [ffmpeg 101442](https://evermeet.cx/ffmpeg/) (earlier versions may work)

## Model

The code simulates particles on a $`d`$-dimensional hypercubic lattice of configurable
extent in each direction ($`L_1, L_2, \ldots, L_d`$). At most one particle occupies
a site, and each particle has a notion of its current direction. At unit rate,
a particle attempts to hop into the site in front, with the move rejected if the
receiving site is occupied. At rate $`\frac{\alpha_k}{d}`$, a particle changes its
direction so as to align with lattice axis $`k`$ (with equal probability in the
positive and negative directions). The total reorientation rate is $`\bar{\alpha}`$,
the arithmetic mean of the rates $`\alpha_k`$. We include the factor $`d`$ so that
isotropic reorientation with rate $`\alpha`$ is obtained if all $`\alpha_k`$ are
set to $`\alpha`$. Note that a possible outcome of the reorientation is that the
particle has the same direction as it did before the reorientation event occurred.

## Simulation algorithm

The simulation uses a waiting-time algorithm, wherein waiting times to the next
event for a given particle are generated from an exponential distribution with
the appropriate rate.

As described in [these slides](https://owncloud.gwdg.de/index.php/s/8fqKIueXrShkjJZ),
there is a single event type per particle, referred to as a _hop attempt_. For
the hop attempt at the head of the event queue (i.e., with the earliest system
time) we execute the following steps:
- if the particle cannot move in any direction (all 2d neigbouring sites are
occupied) the hop is not attempted, the particle becomes inactive and no
future event is queued.
- otherwise, the hop direction that applied at the previous hop attempt
is retained with probability $`\exp(-\bar{\alpha}\theta)`$ where $`\theta`$ is the time
since the last hop attempt. (Note: the case where a particle cannot move and
becomes inactive does not constitute a hop attempt for this purpose). Otherwise
the movement axis $`k`$ is chosen with probability $`\frac{\alpha_k}{d\bar{\alpha}}`$,
and the direction of motion set to be in the positive or negative direction along
that axis with equal probability.
- then, if the site in front of the particle (in its direction of motion) is
vacant, the particle moves into that site. Any inactive particles are adjacent
to the departure site are reactivated by assigning a waiting time to the next
hop attempt from an exponential distribution with unit mean.
- the particle that has just moved is always active and also needs a waiting
time assigned.

The initial condition is one in which $`N`$ of the $`L_1 L_2 \ldots L_d`$ sites
are occupied, these sampled from a uniform distribution. The initial directions
are also assigned uniformly to each particle. Active particles are then identified
and a waiting time to its first hop attempt assigned from the exponential distribution
with unit mean.

## Implementation

The event queue is managed by a C++ standard library `priority_queue` object.
The code that looks after this is in `Scheduler.hpp`. An event is a combination of
a system time at which the event is scheduled to occur, and a function that is
called when the event occurs. This allows a lot of flexibility: the C++ lambda
function syntax allows this function to perform arbitrary actions, including
scheduling future events. (This is how the event queue is continually replenished).
The standard priority queue orders with the 'largest' item at the head of the queue:
therefore the comparison function that the queue utilises has the direction of comparison
reversed so that earlier events happen first. It is mainly the definition of this function
that uses C++17 specific features; it can be backported to C++11 (albeit with uglier syntax)
if people are saddled with elderly compilers.

The simulation has been designed with the expectation that a particle hop is the most
expensive event (and that in a dense system, these will be rare). For each site, we
maintain a record of how many neighbours are occupied: this makes the test as to whether
to deactivate a particle very simple. When a particle hops, we do need to update the
neighbours of both the source and target sites to reflect the particle movement. In doing
so we can also check for particle activation: since it is possible for inactive particles
to be in the event queue, we need an additional flag to record this (the information
cannot be gleaned simply by looking at whether all the adjacent sites are occupied).

We also track both particles and vacancies with a numeric id. In the initial condition
each occupied site receives an id in the range $`0`$ to $`N-1`$ which uniquely identifies
the particle. Each vacant site receives an id in the range $`N`$ to $`(L_1 L_2 \ldots L_d) - 1`$
which uniquely identifies the vacancy. When a particle moves into a vacancy, the two
ids are exchanged.

The code includes a large number of consistency checks which slow the simulation
down but can be disabled by defining the NDEBUG macro.

## Building the simulation

A shell script `build.sh` is provided to facilitate building the simulation.
Only the two provided `.hpp` files, and standard library routines, are referenced
so it should build on any compliant C++17 compiler. The default is `clang++`
but this can be changed to e.g. `g++` by setting the `CXX` environment variable.
Two build modes are supported. The default is for all optimisations to be enabled,
including the suppression of consistency checks. While developing the code, it
is advised to build in debug mode, which can be achieved by adding `--debug` as
an option to the build script.

To summarise:

```
./build.sh
```

performs a default build with `clang++` (optimised for speed; consistency checks
disabled).

```
CXX=g++ ./build.sh
```

performs a default build with an alternate compiler (here `g++`).

```
./build.sh --debug
```

performs a debug build (no speed optimisations; consistency checks enabled).

In all cases the executable is called `crumble`.

In compiling, both the `-Wall` and `-Werror` switches are enabled, so that
all warnings are reported and halt compilation.

## Running the simulation

By default, a system of $`L=100`$ sites and $`N=50`$ particles with a reorientation
rate $`\alpha=1`$ is simulated. (Recall that the hop rate is always unity). This
can be configured on the command line with the following options:
- `-L` sets the system size; to simulate a multi-dimensional system, separate
the size of each dimension with a space. E.g. `-L 50 100` simulates a 50x100
system.
- `-N` sets the number of particles. No check is made to ensure that this does
not exceed the number of sites on the lattice. E.g. `-N100` puts 100 particles
on the lattice.
- `-t` sets the reorientation rate(s). E.g., `-t0.01` sets the total reorientation
rate to 1/100th of the hop rate, and the particles reorient uniformly over the
$`2d`$ possible different directions. In a two dimensional simulation, `-t 1.5 0.5`
would have a total reorientation rate equal to the hop rate (because the mean
of 1.5 and 0.5 is 1), and with reorientation into the first dimension taking
place three times as often as into the second dimension.

The output comprises a header that records the simulation parameters and output
options. Each header line begins with a `#` character so it can be skipped easily.
The default output is a snapshot, one per line, that indicates the state of each site
with two numbers ($`D`$, $`A`$) where $`D`$ is the direction of motion and $`A`$ is
a binary variable that indicates if the site is active or not. This output is
designed mostly to assess how well the code is functioning, rather than for
doing serious data analysis.

Sites are output in `Fortran` order: that is, if $`x_1, x_2, \ldots`$ are positions
along each axis with $`0 \le x_i < L_i`$, then $`x_1`$ changes fastest, $`x_2`$ next
fastest and so on. (This is opposite to the standard C order, and also opposite
to the order that the `numpy` routines used by the Python visualisation code
expects; the dimensions are reversed in the analysis code, so in 2d the first
dimension specifies the height and the second the width, which is opposite to
what might be expected).

The direction $`D`$ is interpreted as follows. The integral part of $`(D+1)/2`$ determines
which axis (in the Fortran ordering) the particle is aligned on. That is, $`D=1,2`$
have the particle moving along the axis of length $`L_1`$; $`D=3,4`$ have the particle
moving along the axis of length $`L_2`$ and so on. If $`D`$ is odd, the particle is
moving in the positive direction along that axis, otherwise it is moving in the
negative direction.

Snapshots are generated at regular intervals during a measurement phase that begins
after an initial burnin period has expired. Recall that the unit of time is set
by the mean time between particle hop attempts. The default burnin period is 1000
units; the default measurement period is 5000 units and the measurement interval
is 2.5 units. These can be changed with the following command-line options:
- `-b` sets the burnin time (e.g., `-b2000` increases it to 2000 units)
- `-u` sets the measurement ('until') time (e.g., `-b1000` reduces it to 1000 units)
- `-e` sets the measurement interval ('every') (e.g., `-e10` increases it to 10 units)

Other output formats can be specified as a single positional argument to the
command line.

The arguments `particles` and `vacancies` generate a compact
output format is also available, which is designed for analysis of
the particle or vacancy dynamics. It comprises a pair of numbers ($`K`$,$`n`$)
where $`K`$ is the id of the particle or vacancy and $`n`$ is the site number
where the particle or vacancy is located. The particle/vacancy id $`K`$ is an
integer that runs from $`0`$ to $`N-1`$ (for particles) and $`N`$ to $`V-1`$
(for vacancies, where $`V`$ is the total number of sites on the lattice). In
both cases the id follows the particle or vacancy as it moves around the lattice.
The site numbering follows the same Fortran ordering as described above. The value
$`n \mod L_1`$ gives the position along the first dimension,
$`(n \div L_1) \mod L_2`$ along the second dimension, and so on. The
(mutually-exclusive) command-line switches `-P` and `-V` select the compact
output that tracks particles and vacancies, respectively.

The argument `clusters` generates the cluster size distribution of particles
and vacancies. Each line of the output comprises a sequence of pairs of count
$`p_1, v_1`$, $`p_2, v_2`$, $`p_3, v_3`$, $`\ldots`$, where $`p_k`$ is the number
of particle clusters of size $`k`$, and $`v_k`$ is the number of vacancy clusters
of size $`k`$. The number of pairs output can vary from one time point to the next,
depending on the size of the largest cluster. A particle cluster is defined as a set
of neighbouring sites all containing particles bounded by vacancies. A vacancy
cluster comprises a set of vacancies bounded by particles.

By default a single measurement is taken at each measurement interval and summed
over and output as a single histogram at the end. This saves generating large
data files and summing over them (slowly) in the Python analysis code. It
is possible to obtain a locally-averaged histogram at each measurement interval
to track the time dependence of the cluster size distribution. The options:
- `-a` sets the number of timepoints to average over at the start of each
measurement interval
- `-i` sets the interval between these local averaging timepoints.
Note that the number of timepoints multiplied by the interval must not exceed
the measurement interval set by `-e`.

## Visualising the output

The Python script `visualise.py` can visualise the output from `crumble` in a
variety of cases. Unless specific information about site activity or particle
direction is required, it can utilise either the full or the compact output
format. By default the script reads from standard input but data can also
be read from files. Some visualisations combine multiple simulation outputs;
others (most notably movies) make use of the first data set presented.

Because the script makes heavy use of pure Python loops, generating a visualisation
can take many orders of magnitude longer than running the original simulation...

The type of visualisation is specified with the `-t` option to the script.
We briefly describe each below.

### Animation of particles showing direction of motion

This is selected with the `-t arrows` option.

The script should be presented
with a sequence of full snapshots from a single simulation run. (Additional
sequences are ignored). It can visualise one and two-dimensional systems. One
frame is generated per measurement. Usually it is most convenient to pipe in the output
from `crumble` directly to `visualise.py`:

```
./crumble -L 50 50 | ./visualise.py -t arrows
```

The movie can be exported in mp4 format, as long as the `ffmpeg` binary is installed
somewhere on your path. Use the `-F` option for this, e.g.

```
./crumble -L 50 50 | ./visualise.py -t arrows -F movie.mp4
```

In my version of matplotlib, there appears to be a bug that prevents previous
frames being cleared before the next one is drawn. (Either that or I have misunderstood
how animations are supposed to work, which is possible, as the documentation is a little
terse). The alternative 'blocks' style of animation does not suffer from this problem.

By default, all particles have the same colour. Add the `-a` switch to dim inactive
particles: this gives a sense of how well the event queue adapts to the particles
that can actually move.


### Animation of particle locations only

This is selected with the `-t blocks` option.

This is a simplification of the `arrows` visualisation that simply displays
particles as shaded squares in a grid. It can work with both full snapshots and
the compact output. Again, active particles can be differentiated
with the `-a` option, and the movie export via `ffmpeg` seems more reliable.


### Algorithm performance graphs

These are selected with the `-t inefficiency` and `-t gain` options.

Each takes multiple sequences of full simulation outputs (e.g., as separate files
or from the concatenated output of several runs) and generates a graph of inefficiency
or gain as a function of particle density, with a different series for different
values of $`\alpha`$.

Inefficiency is defined as the ratio between the length of the event queue
and its shortest possible length given the state of the system, averaged over
timepoints. (The shortest possible length is the number of particles that are
adjacent to at least one vacancy, and therefore could move if pointing in the
right direction). This quantifies how well the algorithm adapts to the presence
of immobile particles.

Gain is defined as the ratio between the maximum length of the
event queue (which is equal to the number of particles) and the actual length
of the event queue, averaged over timepoints. This quantifies how much is gained
relative to a naive algorithm that does not have any optimisations.

The graphs can be exported as PDF with the `-F` option, e.g.,

```
./visualise.py *.dat -t inefficiency -F alldata-inefficiency.pdf
```


### Inter-particle and inter-vacancy gap distributions

These are selected with the `-t nndist` option. This analysis requires the
condensed output (in principle it could be written also to use the full output
but the condensed output is much easier to work with). It computes a histogram
of distance between nearest neighbours (rounded up to the nearest integer in
more than 1d): at each timepoint a focal particle or vacancy is chosen at random
and the distance to its nearest neighbour calculated. The plot compares against
the distribution for a gas (same number of particles or vacancies distributed
uniformly over the lattice).

Separate data series are accepted; however they can only differ by the tumble
rate. You also can't mix and match particle and vacancy distributions, even if
their number is the same.

These histograms can again be exported as PDF with the `-F` option.


### Mean-square displacement plots

These are selected with the `-t msd` option. This analysis requires the condensed
output, as particle or vacancy ids are needed to determine how far they have moved
between snapshots. At each time point in the output, the distance moved relative
to the initial timepoint in the output is calculated. (This is the shortest distance
across all periodic images of the system). The square of this distance
is then averaged over all particles/vacancies and multiple runs. In this way we
measure the mean-square particle/vacancy displacement as a function of time.

For example, the following command generates a plot of the mean-square vacancy
displacement on a system with 100 sites and 5 vacancies, averaged over 100 runs:

```
(repeat 100 do; ./crumble vacancies -L 100 -N 95 -t0.01; done) | ./visualise.py -t msd
```


### Cluster size distributions

These are selected with the `-t pdist` (for particles) and `-t vdist` (vacancies)
options. This analysis requires the clusters output.

This visualisation is very new and the details are likely to change.
