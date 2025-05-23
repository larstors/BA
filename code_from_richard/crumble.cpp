/**
 * crumble: condensed run-and-tumble simulation
 */

/*
General notes on the code from Lars:

  - I did not touch the decoder, so they probably do not work

  - The main change I did to the code is:

    (a) I added two lattices. It is probably not the best implementation, but I basically copied
        the lattice class by Richard and made the appropriate adjustments
    (b) I included a variable upper bound for allowed particles per site. This is not 100% perfect, as I didnt
        find a way to globally introduce this while also defining the size in the classes. I have commented
        the appropriate places in the code where this takes place (or rather, does not take place)

  - Because I re-use some of the output multiple times to save computational effort, I changed the output of the
    script to give me txt files, that I then read in in python for further consideration/analysis. While not the
    most efficient it allows for some flexible output of different data.

  - The actual output has become quite a mess. Especially since I started looking at many different things.
    Hence, one has to navigate the appropriate "if" statements to find the proper output. I am well aware that
    this is not he most efficient way to do this. Thus, I would for now not recommend looking into the actual
    plotting scripts, because the consist of a spaghetti of code that I would not wish upon my worst enemy.

  - As you might have guessed, I am not the most experienced coder, especially not in C++ (this is my first time
    using it...), so there are bound to be sections that are wildly inefficient. As for the python scripts I went
    for a brute force approach. Most parts are just repetitive plotting with some interludes of transforming the
    data into suitable formats.

*/

/* TODO LIST

  - add displacement to Tri and Hex
  - fix the function mess:
    - dont need to recalculate the cluster sizes each time -> also less code
    - optimise them -> feel like there are some unnecessary loops there somewhere
  - overall optimisation
    - code is very slow for some systems/outputs


*/

#include <vector>
#include <list>
#include <map>
#include <valarray>
#include <numeric>
#include <functional>
#include <random>
#include <iostream>
#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>

#include "CLI11.hpp"
#include "Scheduler.hpp"

using namespace std;
using direction_t = unsigned char;
using hist_t = std::valarray<unsigned>;
using vec = std::vector<unsigned>;
using vec_d = std::vector<double>;

/**
 * Model parameters. We put these in a struct so we can pass as a single
 * argument to Lattice and also so all the defaults are in one place.
 */
struct Parameters
{
  std::vector<unsigned> L{100};   // Lengths of each lattice dimension
  unsigned N = 50;                // Number of particles
  std::vector<double> alpha{1.0}; // Scaled tumble rate in each dimension (actual rate is alpha/d)
  // ! Addition by Lars
  unsigned n_max = 1;  // Max occupation number on each lattice site
  bool tagged = false; // boolean for tagged particle, by default false
};

// We need to forward declare these so the compiler can cope with the friend declarations
template <typename Engine>
class SnapshotWriter;
template <typename Engine>
class VacancyWriter;
template <typename Engine>
class ParticleWriter;
template <typename Engine>
class ParticleWriterWID;
template <typename Engine>
class ClusterWriter;

template <typename Engine>
class ParticleWriterDetails;

// new writers for the two new lattices
template <typename Engine>
class TriangleParticleWriter;
template <typename Engine>
class HexagonalParticleWriter;
template <typename Engine>
class HexDirectionWriter;

// ! There will be some changes is this class, mostly the implementation of occupation number

template <typename Engine>
class Lattice
{

  // For output in different formats
  friend class SnapshotWriter<Engine>;
  friend class VacancyWriter<Engine>;
  friend class ParticleWriter<Engine>;
  friend class ClusterWriter<Engine>;
  friend class ParticleWriterWID<Engine>;
  friend class ParticleWriterDetails<Engine>;

  static constexpr unsigned n_max = 10; // ! Nicer way of doing this?

  // Data associated with each site; by default all of these are set to zero
  // which represents a vacant site
  // now we use vectors here, as we allow up to n_max particles per site
  struct Site
  {
    std::vector<unsigned> id = std::vector<unsigned>(n_max);              // Particle / vacancy id
    std::vector<bool> occupied = std::vector<bool>(n_max);                // There is a particle here
    std::vector<bool> active = std::vector<bool>(n_max);                  // A move event is scheduled
    direction_t neighbours;                                               // Number of neighbours that are occupied
    std::vector<direction_t> direction = std::vector<direction_t>(n_max); // Direction of last hop attempt
    std::vector<double> hoptime = std::vector<double>(n_max);             // Time of last hop attempt
    unsigned present = 0;                                                 // Number of particles present at site. Has to be <= n_max
    std::vector<double> last_jump = std::vector<double>(n_max);           // time of last jump made
    std::vector<int> x_coor = std::vector<int>(n_max);                    // x coordinate (that actually just counts steps in each direction)
  };

  Parameters P;                                         // A local copy of the model parameters
  std::vector<Site> sites;                              // Representation of the sites
  Engine &rng;                                          // Source of noise: this is a reference as there should only be one of these!
  std::discrete_distribution<unsigned> anyway, initial; // Distribution over tumble directions
  std::exponential_distribution<double> run;            // Distribution of times between run and tumble events
  Scheduler S;                                          // Keeps track of the event queue

  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours. We have periodic boundary conditions
  auto neighbours(unsigned n) const
  {
    std::vector<unsigned> nbs(2 * P.L.size());
    unsigned below = 1;
    for (unsigned d = 0; d < P.L.size(); ++d)
    {
      unsigned L = P.L[d];
      unsigned above = below * L;
      // x is the position along dimension d
      // y is the contribution to the site index from all other dimensions
      unsigned x = (n / below) % L;
      unsigned y = (n % below) + (n / above) * above;
      // Neighbours in the increasing and decreasing directions
      nbs[2 * d] = y + ((x + 1) % L) * below;
      nbs[2 * d + 1] = y + ((x + L - 1) % L) * below;
      below = above;
    }
    return nbs;
  }

  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours in the forward direction along each axis.
  // We still have periodic boundary conditions
  auto forward_neighbours(unsigned n) const
  {
    std::vector<unsigned> nbs(P.L.size());
    unsigned below = 1;
    for (unsigned d = 0; d < P.L.size(); ++d)
    {
      unsigned L = P.L[d];
      unsigned above = below * L;
      // x is the position along dimension d
      // y is the contribution to the site index from all other dimensions
      unsigned x = (n / below) % L;
      unsigned y = (n % below) + (n / above) * above;
      // Neighbours in the increasing and decreasing directions
      nbs[d] = y + ((x + 1) % L) * below;
      below = above;
    }
    return nbs;
  }

  // Place a particle with given direction and hop time at site n;
  // neighbouring sites will be accordingly adjusted
  void place(unsigned n, unsigned id, direction_t d, double t, unsigned index, double x)
  {
    sites[n].id[index] = id;
    sites[n].direction[index] = d;
    sites[n].hoptime[index] = t;
    sites[n].last_jump[index] = t;
    if (!sites[n].occupied[index])
    {
      sites[n].occupied[index] = true;
      for (const auto &m : neighbours(n))
        ++sites[m].neighbours;
    }
    sites[n].present++;
    sites[n].x_coor[index] = x;
  }

  // Schedule a hop event for a particle at site n
  void schedule(unsigned n, unsigned index)
  {
    assert(sites[n].occupied[index]);
    S.schedule(run(rng), [this, n, index]()
               {
      assert(sites[n].occupied[index]);
      assert(sites[n].active[index]);
      // If there are no local vacancies, mark this particle as inactive and exit
      if(sites[n].neighbours == 2*P.L.size() * P.n_max) {
        // std::cout << "Can't move from "; decode(n); std::cout << " deactivating" << std::endl;
        sites[n].active[index] = false;
      } else {
        // because of some c++ stuff for alpha=0 I did this, probably not the best solution (we also never went below this
        // limit in the thesis, so it didn't cause any trouble)
        if (P.alpha[0] - 1e-5 > 0){
            if (tumble(rng) < S.time() - sites[n].hoptime[index]) {
              sites[n].direction[index] = anyway(rng);
            }
        }
        sites[n].hoptime[index] = S.time();
        // Get the sites adjacent to the departure site
        auto dnbs = neighbours(n);
        if (sites[dnbs[sites[n].direction[index]]].present < P.n_max) {
          // get the lowest index at the target site where we can put the particle
          auto itr = std::find(sites[dnbs[sites[n].direction[index]]].occupied.begin(), sites[dnbs[sites[n].direction[index]]].occupied.end(), false);
          unsigned k = std::distance(sites[dnbs[sites[n].direction[index]]].occupied.begin(), itr);
          // perform jump
          if (k < (P.n_max) && !sites[dnbs[sites[n].direction[index]]].occupied[k]){

                    // Get the id of the vacancy that is being displaced
                    unsigned vid = sites[dnbs[sites[n].direction[index]]].id[k];

                    // Deactive the departure site; also mark it empty
                    sites[n].occupied[index] = sites[n].active[index] = false;
                    
                    // as it moves from n, we have one less present at n
                    sites[n].present--;
                    assert(sites[n].present >= 0);


                    // if it looks to the right add 1
                    if (sites[n].direction[index] == 0){
                      sites[n].x_coor[index]++;
                    }
                    // if it looks to the left remove 1
                    else if (sites[n].direction[index] == 1){
                      sites[n].x_coor[index]--;
                    }

                    // Place a particle on the target site; it has the same direction, hoptime etc. as the departing particle
                    place(dnbs[sites[n].direction[index]], sites[n].id[index], sites[n].direction[index], sites[n].hoptime[index], k, sites[n].x_coor[index]);
                    
                    // Move the vacancy id onto the departure site
                    sites[n].id[index] = vid;
                    

                    // Now go through the neighbours of the departure site, update neighbour count and activate any
                    // that can now move. Note the particle this is at the target site is included in this list
                    // and will be activated accordingly
                    for (const auto& m : dnbs) {
                        --sites[m].neighbours;
                        for (unsigned i = 0; i < P.n_max; i++){
                          if (sites[m].occupied[i] && !sites[m].active[i]) schedule(m, i);
                        }
                    }
                    
                  }
                
                }
                else {
                    // std::cout << "Didn't move from "; decode(n); std::cout << std::endl;
                    // This site is still active, so schedule another hop
                    schedule(n, index);
                }
            } });
    sites[n].active[index] = true;
    // std::cout << "Scheduled "; decode(n); std::cout << std::endl;
  }

  // For testing
  void decode(unsigned n)
  {
    std::cout << "[ ";
    for (const auto &L : P.L)
    {
      std::cout << (n % L) << " ";
      n /= L;
    }
    std::cout << "]";
  }

  // Check the lattice state is consistent
  bool consistent()
  {
    unsigned active = 0, occupied = 0;
    std::set<unsigned> ids;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned k = 0; k < P.n_max; k++)
      {
        // Check each site has a unique id
        if (ids.count(sites[n].id[k]))
          return false;
        ids.insert(sites[n].id[k]);
        // Check that empty sites are also inactive
        if (!sites[n].occupied[k])
        {
          if (sites[n].active[k])
            return false;
          // Nothing left to do if empty
          continue;
        }
        ++occupied;
        // Check that the neighbour count is correct
        unsigned nbs = 0;
        for (const auto &m : neighbours(n))
        {
          for (unsigned i = 0; i < P.n_max; i++)
          {
            if (sites[m].occupied[i])
              ++nbs;
          }
        }
        if (nbs != sites[n].neighbours)
          return false;
        // Check that mobile particles are active
        if (nbs < 4 * P.n_max && !sites[n].active[k])
          return false;
        if (sites[n].active[k])
          ++active;
      }
    }
    // Check we've not lost any particles
    return occupied == P.N && active == S.pending();
  }

public:
  // have to declare exp distribution here to properly access it later on for the hysteresis
  std::exponential_distribution<double> tumble;
  Lattice(const Parameters &P, Engine &rng) : P(P),                                                                           // NB: this takes a copy
                                              sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
                                              rng(rng),                                                                       // NB: this takes a reference
                                              run(1),                                                                         // run distribution
                                              tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size())   // Tumble time generator: set to the average of the given tumble rates
  {
    // Set up the tumble direction distribution
    std::vector<double> drates(2 * P.L.size());
    // and one for initial distribution during initialisation (this is needed as alpha=0 doesn't properly work here, at least not for me...)
    std::vector<double> drates_initial(4);
    for (unsigned d = 0; d < P.L.size(); ++d)
    {
      drates[2 * d] = drates[2 * d + 1] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
      // because I always use isotropic alpha, I just define the initial direction distribution as this
      drates_initial[2 * d] = drates_initial[2 * d + 1] = 1.0;
    }
    anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());
    initial = std::discrete_distribution<unsigned>(drates_initial.begin(), drates_initial.end());

    // vector with all site indices (so n and i combined)
    std::vector<unsigned> position;
    for (unsigned i = 0; i < P.n_max; i++)
    {
      for (unsigned n = 0; n < sites.size(); n++)
      {
        position.push_back(i * sites.size() + n);
      }
    }

    // Now placing particles
    unsigned possibilities = position.size();
    unsigned unplaced = P.N;
    unsigned id = 0;
    while (unplaced > 0)
    {
      // take an index out of the remaining sites
      unsigned index = std::uniform_int_distribution<unsigned>(0, possibilities - 1)(rng);
      unsigned l = position[index];
      // position of lattice
      unsigned n = l % sites.size();
      unsigned i = l / sites.size();

      // placing particle and removing site from list
      place(n, id, initial(rng), 0.0, i, 0);
      position.erase(position.begin() + index);
      position.push_back(l);

      // increase id so no duplicate
      id++;

      // take away possibility to take already used indices
      possibilities--;
      unplaced--;
    }

    // add id to vacancies
    for (unsigned index = 0; index < P.n_max; index++)
    {
      for (unsigned n = 0; n < sites.size(); ++n)
      {
        if (sites[n].occupied[index] == false)
        {
          sites[n].id[index] = id;
          id++;
        }
      }
    }

    // Activate particles that can move, and schedule a hop accordingly
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned k = 0; k < P.n_max; ++k)
      {
        if (sites[n].occupied[k] && sites[n].neighbours < 4 * P.n_max)
          schedule(n, k);
      }
    }
    assert(consistent());
  }

  // Iterates the simulation until the given time; returns the actual time run to
  double run_until(double time)
  {
    while (S.advance(time))
      ;
    assert(consistent());
    return S.time();
  }

  // function that allows us to adjust alpha mid run (needed for hysteresis)
  void set_new_lambda(std::exponential_distribution<double> *exp_dis, double val)
  {
    typename std::exponential_distribution<double>::param_type new_lambda(val);
    exp_dis->param(new_lambda);
  }

  // Sample all particle directions from the distribution that applies at the current instant
  // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
  void realise_directions()
  {
    for (auto &site : sites)
    {
      for (unsigned k = 0; k < P.n_max; k++)
      {
        if (!site.occupied[k])
          continue;
        if (tumble(rng) < S.time() - site.hoptime[k])
        {
          // At least one tumble event happened since we last attempted a hop; so randomise the direction
          site.direction[k] = anyway(rng);
        }
        // Advance the hop time so we don't generate temporal paradox
        site.hoptime[k] = S.time();
      }
    }
  }

  // calculates the fraction of the system that is jammed
  double motility_fraction()
  {
    double count = 0;
    // loop over entire lattice
    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned k = 0; k < P.n_max; k++)
      {
        // only consider particles
        if (!sites[n].occupied[k])
          continue;
        if (tumble(rng) < S.time() - sites[n].hoptime[k])
        {
          // At least one tumble event happened since we last attempted a hop; so randomise the direction
          sites[n].direction[k] = anyway(rng);
        }
        // if particle is jammed add one to count
        auto dnbs = neighbours(n);
        if (sites[dnbs[sites[n].direction[k]]].present == P.n_max)
          count++;
      }
    }
    // return jammed fraction
    return count / double(P.N);
  }

  // Return cluster size distributions.
  // Element 2n   contains the number of clusters of particles of size n+1
  // Element 2n+1 contains the number of clusters of vacancies of size n+1
  hist_t cluster_distributions() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    hist_t dists(2 * maxsize);
    for (const auto &kv : clusters)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present == 0)
        dists[2 * kv.second.size() - 1]++;
      else
        dists[2 * kv.second.size() - 2]++;
    }
    return dists;
  }

  // I did this far too late, but here is a function to generate the map of clusters.
  // TODO add this instead of tedious other way of doing it
  map<unsigned, std::list<unsigned>> clusters() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    return clusters_nr;
  }

  // while the cluster distribution above only consideres lattice sites, this one also consideres the particles present
  hist_t cluster_distribution_particle_number() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    hist_t dists_nr(2 * maxsize_nr);
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present == 0)
        dists_nr[2 * kv.second.size() - 1]++;
      else
        dists_nr[2 * kv.second.size() - 2]++;
    }
    return dists_nr;
  }

  // Function to determine the size of largest cluster. Note that we will only regard particle clusters here
  // (at least so far)
  size_t max_cluster_size()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    std::size_t max_s = 1;
    for (const auto &kv : clusters)
    {

      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        // std::cout << "-------------------" << endl;
        // std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
        max_s = std::max(max_s, kv.second.size());
        // std::cout << max_s << endl;
      }
    }

    // std::cout << max_s << endl;

    return max_s;
  }

  // function for max cluster size by counting particles
  size_t max_cluster_size_nr()
  {
    std::map<unsigned, std::list<unsigned>> clusters_nr = clusters();

    // maxsize above also counts vacancies
    std::size_t max_s = 1;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
        max_s = std::max(max_s, kv.second.size());
    }
    return max_s;
  }

  // Function to determine the mean cluster size by simple taking the sum off all clusters (not vacant ones) and dividing by number of clusters. Note that we will only regard particle clusters here
  // (at least so far)
  double avg_cluster_size_nr()
  {
    std::map<unsigned, std::list<unsigned>> clusters_nr = clusters();

    // variable for mean and the count of clusters
    double mean = 0;
    double count = 0;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        mean += kv.second.size();
        count += 1;
      }
    }

    mean = mean / count;

    return mean;
  }

  // function to return the number of particle clusters
  unsigned number_cluster()
  {
    std::map<unsigned, std::list<unsigned>> clusters_nr = clusters();

    // variable for mean and the count of clusters
    unsigned count = 0;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        count += 1;
      }
    }

    return count;
  }

  // determine surface and volume of clusters by lattice sites, not particle number
  vec surface_volume()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // vector for output of surface and volume of clusters
    // the index 2i is for volume and 2i+1 is for surface
    vec output;

    for (const auto &kv : clusters)
    {
      // variable for surface and check
      unsigned surf = 0;
      unsigned check = 0;
      // check whether particle cluster
      if (sites[kv.first].present > 0)
      {
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours(n))
          {
            if (sites[m].present == 0)
            {
              check++;
              continue;
            }
          }
          if (check > 0)
            surf++;
        }
        // std::cout << surf << " " << kv.second.size() << endl;
        output.push_back(kv.second.size());
        output.push_back(surf);
      }
    }
    return output;
  }

  // determine surface and volume of clusters by particle number
  vec surface_volume_nr()
  {
    std::map<unsigned, std::list<unsigned>> cluster = clusters();

    // vector for output of surface and volume of clusters
    // the index 2i is for volume and 2i+1 is for surface
    vec output;
    vec index;
    for (const auto &kv : cluster)
    {
      // variable for surface and check
      unsigned surf = 0;
      unsigned check = 0;
      // check whether particle cluster
      if (sites[kv.first].present > 0)
      {
        for (const auto &n : kv.second)
        {
          // for duplicate index
          unsigned dup = 0;
          for (const auto &m : index)
          {
            if (m == n)
              dup++;
          }
          index.push_back(n);
          check = 0;
          if (dup == 0)
          {
            for (const auto &m : neighbours(n))
            {
              if (sites[m].present == 0)
              {
                check++;
              }
            }
          }
          if (check > 0)
            surf += sites[n].present;
        }
        output.push_back(kv.second.size());
        output.push_back(surf);
      }
    }
    return output;
  }

  // vector of position of all sites on surface, only used for visual representation of clusters so far
  vec cluster_surface()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // output vector for position of all surface sites
    vec output;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      unsigned check = 0;
      if (sites[n].present > 0)
      {
        for (const auto &m : neighbours(n))
        {
          if (sites[m].present < 1)
          {
            check++;
            continue;
          }
        }
        if (check > 0)
          output.push_back(n);
      }
    }
    return output;
  }

  // function to give out position of a cluster in the lattice that has size between lower and upper bound
  // ! take out? I am not using this anymore, I think (unsure so letting it stay for now)
  vec single_cluster(unsigned clust_min, unsigned clust_max)
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // vector for output of surface and volume of clusters
    vec output;

    for (const auto &kv : clusters)
    {
      // check whether particle cluster
      if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max)
      {
        // to get surface and volume
        unsigned surf = 0;
        unsigned check = 0;
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours(n))
          {
            if (sites[m].present == 0)
            {
              check++;
              continue;
            }
          }
          if (check > 0)
            surf++;
        }
        // add surface and volume to beginning of output
        output.push_back(kv.second.size());
        output.push_back(surf);
        // add position of particles
        for (const auto &n : kv.second)
        {
          output.push_back(n);
        }
        break;
      }
    }

    return output;
  }

  // checks whether any clusters of size between bounds is present (measured by lattice sites, not particles)
  unsigned clust_size(unsigned clust_min, unsigned clust_max)
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    unsigned output = 0;

    for (const auto &kv : clusters)
    {
      // check whether particle cluster
      if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max)
      {
        output = 1;
      }
    }

    return output;
  }

  // histogram of how many neighbouring sites are occupied
  hist_t particle_neighbour_dist()
  {
    hist_t dist(5);
    unsigned count = 0;
    for (unsigned n = 0; n < sites.size(); n++)
    {
      count = 0;
      if (sites[n].present > 0)
      {
        for (const auto &m : neighbours(n))
        {
          if (sites[m].present > 0)
            count++;
        }
        dist[count]++;
      }
    }
    return dist;
  }

  // local density (smoothed out a bit, i.e. weighted by 1/2 for site and 1/8 for neighbours)
  // the part regarding structure factor didnt end up in the thesis, but this is the starting point.
  vec_d density()
  {
    vec_d den;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      double local = 0.5 * double(sites[n].present);
      for (const auto &m : neighbours(n))
      {
        local += 1.0 / 8.0 * double(sites[m].present);
      }
      den.push_back(local);
    }

    return den;
  }

  // stopping time -> array of time since last jump occured for each particle
  vec_d stopping(double t)
  {
    vec_d stopping_time;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < P.n_max; i++)
      {
        if (sites[n].occupied[i] == true)
          stopping_time.push_back(t - sites[n].last_jump[i]);
      }
    }
    return stopping_time;
  }

  // function for returning perimeter order parameter (see thesis)
  double perimeter()
  {
    // get distribution and surface-volume
    hist_t dist = cluster_distribution_particle_number();
    vec sv = surface_volume_nr();

    double second_moment = 0;
    double first_moment = 0;

    // calculate appropriate stuff
    for (unsigned i = 0; i < dist.size(); i += 2)
    {
      double count = 0;
      for (unsigned k = 0; k < sv.size(); k += 2)
      {
        if (sv[k] == ((i + 2) / 2))
        {
          second_moment += ((i + 2) / 2) * sv[k + 1];
          // std::cout << sv[k] << " " << sv[k+1] << " " << (i+2)/2 << endl;
          count++;
        }
      }
      first_moment += dist[i] * ((i + 2) / 2);
    }
    return second_moment / (first_moment * double(P.N));
  }

  // function for returning w(_N) order parameter (see thesis)
  double w_N()
  {
    // get distribution and surface-volume
    hist_t dist = cluster_distribution_particle_number();

    double second_moment = 0;
    double first_moment = 0;

    // calculate appropriate stuff
    for (unsigned i = 0; i < dist.size(); i += 2)
    {
      second_moment += dist[i] * ((i + 2) / 2) * ((i + 2) / 2);
      first_moment += dist[i] * ((i + 2) / 2);
    }
    return second_moment / first_moment;
  }

  // array with occupation numbers
  vec occ_array()
  {
    vec output;
    for (unsigned i = 0; i < sites.size(); i++)
    {
      output.push_back(sites[i].present);
    }
    return output;
  }

  // function to get displacement of particles
  vec_d x_coor_3()
  {
    vec_d x = vec_d(P.N);
    for (int n = 0; n < sites.size(); n++)
    {
      for (int i = 0; i < P.n_max; i++)
      {
        if (sites[n].occupied[i])
          x[sites[n].id[i]] = sites[n].x_coor[i];
      }
    }
    return x;
  }

  // histogram of distribution of n on sites
  hist_t occupation_prob()
  {
    hist_t output(P.n_max + 1);
    for (unsigned n = 0; n < sites.size(); n++)
    {
      output[sites[n].present]++;
    }
    return output;
  }
};

// ! Following classes are additions by Lars. Note that in the triangular case it is just the one as above with some very slight adjustments
// ! while the one for the hexagonal lattice has more drastic changes.

// TODO add diffusion and displacement for Tri and Hex!

template <typename Engine>
class Triangle_lattice
{

  // For output in different formats
  friend class SnapshotWriter<Engine>;
  friend class VacancyWriter<Engine>;
  friend class TriangleParticleWriter<Engine>;
  friend class ClusterWriter<Engine>;

  Parameters P; // A local copy of the model parameters

  static constexpr unsigned n_max = 10; // ! Nicer way of doing this?

  // Data associated with each site; by default all of these are set to zero
  // which represents a vacant site
  struct Site
  {
    std::vector<unsigned> id = std::vector<unsigned>(n_max);              // Particle / vacancy id
    std::vector<bool> occupied = std::vector<bool>(n_max);                // There is a particle here
    std::vector<bool> active = std::vector<bool>(n_max);                  // A move event is scheduled
    direction_t neighbours;                                               // Number of neighbours that are occupied
    std::vector<direction_t> direction = std::vector<direction_t>(n_max); // Direction of last hop attempt
    std::vector<double> hoptime = std::vector<double>(n_max);             // Time of last hop attempt
    unsigned present = 0;                                                 // Number of particles present at site. Has to be <= n_max
    std::vector<double> last_jump = std::vector<double>(n_max);           // time of last jump made
    std::vector<bool> tagged = std::vector<bool>(n_max);                  // for tagged particle
    std::vector<int> per_hop = std::vector<int>(n_max);                   // for how many times particle jumped periodic boundary.
                                                                          //+1 for pos, -1 for neg direction
  };

  std::vector<Site> sites;                              // Representation of the sites
  Engine &rng;                                          // Source of noise: this is a reference as there should only be one of these!
  std::discrete_distribution<unsigned> anyway, initial; // Distribution over tumble directions
  std::exponential_distribution<double> run;            // Distribution of times between run and tumble events
  Scheduler S;                                          // Keeps track of the event queue

  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours. We have periodic boundary conditions
  auto neighbours(unsigned n) const
  {

    // this is only valid for 2d
    int L = P.L[0];
    std::vector<int> nbs(2 * P.L.size() + 2);

    int x = n % L;
    int y = (n / L);

    // for diagonals. There surely is a better way than this
    int xk, yk;
    int xm = 0;
    int ym;
    int k = n + L + 1;
    int m = int(n) - L - 1;

    xk = k % L;
    if (m % L >= 0)
      xm = m % L;
    else if (m % L < 0)
    {
      xm = L + m % L;
    }

    if (((k - 1) / L - L) >= 0)
    {
      yk = 0;
    }
    else
    {
      yk = ((k - 1) / L);
    }
    if ((m + 1) < 0)
    {

      ym = L - 1;
    }
    else
    {
      int l = std::abs(int(n - L));
      // if (l < 0) l = -l;
      ym = (l / L);
    }
    nbs[0] = (n + 1) % L + y * L;
    if ((int(n) - 1) % L < 0)
      nbs[1] = L + (int(n) - 1) % L + y * L;
    else if ((int(n) - 1) % L + y * L >= 0)
      nbs[1] = (int(n) - 1) % L + y * L;

    nbs[2] = x + ((y + 1) % L) * L;
    nbs[3] = x + ((y - 1 + L) % L) * L;

    nbs[4] = xk + yk * L;
    nbs[5] = xm + ym * L;
    return nbs;
  }

  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours in the forward direction along each axis.
  // We still have periodic boundary conditions
  auto forward_neighbours(unsigned n) const
  {
    // this is only valid for 2d
    int L = P.L[0];
    std::vector<int> nbs(P.L.size() + 1);

    int x = n % L;
    int y = (n / L);

    // for diagonals. There surely is a better way than this
    int xk, yk;
    int k = n + L + 1;

    xk = k % L;

    if (((k - 1) / L - L) >= 0)
    {
      yk = 0;
    }
    else
    {
      yk = ((k - 1) / L);
    }
    nbs[0] = (n + 1) % L + y * L;

    nbs[1] = x + ((y + 1) % L) * L;

    nbs[2] = xk + yk * L;

    return nbs;
  }

  // Place a particle with given direction and hop time at site n;
  // neighbouring sites will be accordingly adjusted
  void place(unsigned n, unsigned id, direction_t d, double t, unsigned index)
  {
    sites[n].id[index] = id;
    sites[n].direction[index] = d;
    sites[n].hoptime[index] = t;
    sites[n].last_jump[index] = t;
    if (!sites[n].occupied[index])
    {
      sites[n].occupied[index] = true;
      for (const auto &m : neighbours(n))
        ++sites[m].neighbours;
    }

    sites[n].present++; // increase the count of particles on lattice site
  }

  // Schedule a hop event for a particle at site n
  void schedule(unsigned n, unsigned index)
  {
    assert(sites[n].occupied[index]);
    // here we added the index to be given as well, so that we can distinguish single particles for higher occupation number
    S.schedule(run(rng), [this, n, index]()
               {
            assert(sites[n].occupied[index]);
            assert(sites[n].active[index]);
            // If there are no local vacancies, mark this particle as inactive and exit
            if (sites[n].neighbours == 6 * P.n_max) {
                sites[n].active[index] = false;

            }
            else {                
                // Make sure that for alpha=0 we dont get some weird changes.
                // ! Note that for lower alpha this threshold has to be adjusted
                if (P.alpha[0] - 1e-5 > 0){
                  if (tumble(rng) < S.time() - sites[n].hoptime[index]) {
                    sites[n].direction[index] = anyway(rng);
                  }
                }
                sites[n].hoptime[index] = S.time();
                // Get the sites adjacent to the departure site
                auto dnbs = neighbours(n);
                if (sites[dnbs[sites[n].direction[index]]].present < P.n_max) {
                  // Find first spot that is empty in target site and then calculate actual index.
                  auto itr = std::find(sites[dnbs[sites[n].direction[index]]].occupied.begin(), sites[dnbs[sites[n].direction[index]]].occupied.end(), false);
                  unsigned k = std::distance(sites[dnbs[sites[n].direction[index]]].occupied.begin(), itr);
                  if (k < (P.n_max) && !sites[dnbs[sites[n].direction[index]]].occupied[k]){
                    
                    // Get the id of the vacancy that is being displaced
                    unsigned vid = sites[dnbs[sites[n].direction[index]]].id[k];
                    // std::cout << "Moving from "; decode(n); std::cout << " deactivating" << std::endl;
                    // Deactive the departure site; also mark it empty
                    sites[n].occupied[index] = sites[n].active[index] = false;
                    
                    // as it moves from n, we have one less present at n
                    sites[n].present--;
                    assert(sites[n].present >= 0);
                    
                    
                    // Place a particle on the target site; it has the same direction and hoptime as the departing particle
                    // std::cout << "Moving to "; decode(dnbs[sites[n].direction]); std::cout << " placing" << std::endl;
                    place(dnbs[sites[n].direction[index]], sites[n].id[index], sites[n].direction[index], sites[n].hoptime[index], k);
                    
                    // Move the vacancy id onto the departure site
                    sites[n].id[index] = vid;
                    
                    
                    // Now go through the neighbours of the departure site, update neighbour count and activate any
                    // that can now move. Note the particle this is at the target site is included in this list
                    // and will be activated accordingly
                    for (const auto& m : dnbs) {
                        --sites[m].neighbours;
                        for (unsigned i = 0; i < P.n_max; i++){
                          if (sites[m].occupied[i] && !sites[m].active[i]) schedule(m, i);
                        }
                    }
                  }
                }
                else {
                    // std::cout << "Didn't move from "; decode(n); std::cout << std::endl;
                    // This site is still active, so schedule another hop
                    schedule(n, index);
                }
            } });
    sites[n].active[index] = true;
    // std::cout << "Scheduled "; decode(n); std::cout << std::endl;
  }

  // For testing
  void decode(unsigned n)
  {
    std::cout << "[ ";
    for (const auto &L : P.L)
    {
      std::cout << (n % L) << " ";
      n /= L;
    }
    std::cout << "]";
  }

  // Check the lattice state is consistent
  bool consistent()
  {
    unsigned active = 0, occupied = 0;
    std::set<unsigned> ids;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned k = 0; k < P.n_max; k++)
      {
        // Check each site has a unique id

        if (ids.count(sites[n].id[k]))
          return false;
        ids.insert(sites[n].id[k]);
        // Check that empty sites are also inactive
        if (!sites[n].occupied[k])
        {
          if (sites[n].active[k])
            return false;
          // Nothing left to do if empty
          continue;
        }

        // Check that the neighbour count is correct
        ++occupied;
        // ! Can we move this out one loop?
        unsigned nbs = 0;
        for (const auto &m : neighbours(n))
        {
          for (unsigned i = 0; i < P.n_max; i++)
          {
            if (sites[m].occupied[i])
              ++nbs;
          }
        }
        if (nbs != sites[n].neighbours)
        {
          return false;
        }
        // Check that mobile particles are active
        if (nbs < 6 * P.n_max && !sites[n].active[k])
          return false;
        if (sites[n].active[k])
          ++active;
      }
    }
    // Check we've not lost any particles
    return occupied == P.N && active == S.pending();
  }

public:
  std::exponential_distribution<double> tumble;
  Triangle_lattice(const Parameters &P, Engine &rng) : P(P),                                                                           // NB: this takes a copy
                                                       sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
                                                       rng(rng),                                                                       // NB: this takes a reference
                                                       run(1),                                                                         // On average, it runs once per time unit
                                                       tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size())   // Tumble time generator: set to the average of the given tumble rates
  {
    // Set up the tumble direction distribution
    std::vector<double> drates(2 * P.L.size() + 2);
    std::vector<double> drates_initial(6);
    for (unsigned d = 0; d < 3; ++d)
    {
      drates[2 * d] = drates[2 * d + 1] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
      drates_initial[2 * d] = drates_initial[2 * d + 1] = 1.0;
    }
    anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());
    initial = std::discrete_distribution<unsigned>(drates_initial.begin(), drates_initial.end());

    std::vector<unsigned> position;
    for (unsigned i = 0; i < P.n_max; i++)
    {
      for (unsigned n = 0; n < sites.size(); n++)
      {
        position.push_back(i * sites.size() + n);
      }
    }

    unsigned possibilities = position.size();
    unsigned unplaced = P.N;
    unsigned id = 0;
    while (unplaced > 0)
    {
      unsigned index = std::uniform_int_distribution<unsigned>(0, possibilities - 1)(rng);

      unsigned l = position[index];

      unsigned n = l % sites.size();
      unsigned i = l / sites.size();

      place(n, id, initial(rng), 0.0, i);

      position.erase(position.begin() + index);
      position.push_back(l);

      id++;
      possibilities--;
      unplaced--;
    }

    for (unsigned index = 0; index < P.n_max; index++)
    {
      for (unsigned n = 0; n < sites.size(); ++n)
      {
        if (sites[n].occupied[index] == false)
        {
          sites[n].id[index] = id;
          id++;
        }
      }
    }

    // Activate particles that can move, and schedule a hop accordingly
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned k = 0; k < P.n_max; ++k)
      {
        if (sites[n].occupied[k] && sites[n].neighbours < 6 * P.n_max)
          schedule(n, k);
      }
    }
    assert(consistent());
  }

  void set_new_lambda(std::exponential_distribution<double> *exp_dis, double val)
  {
    typename std::exponential_distribution<double>::param_type new_lambda(val);
    exp_dis->param(new_lambda);
  }
  // Iterates the simulation until the given time; returns the actual time run to
  double run_until(double time)
  {
    while (S.advance(time))
      ;
    assert(consistent());
    return S.time();
  }

  map<unsigned, std::list<unsigned>> clusters() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    return clusters_nr;
  }
  // Sample all particle directions from the distribution that applies at the current instant
  // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
  void realise_directions()
  {
    for (auto &site : sites)
    {
      for (unsigned k = 0; k < P.n_max; k++)
      {
        if (!site.occupied[k])
          continue;
        if (tumble(rng) < S.time() - site.hoptime[k])
        {
          // At least one tumble event happened since we last attempted a hop; so randomise the direction
          site.direction[k] = anyway(rng);
        }
        // Advance the hop time so we don't generate temporal paradox
        site.hoptime[k] = S.time();
      }
    }
  }

  double motility_fraction()
  {
    double count = 0;
    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned k = 0; k < P.n_max; k++)
      {
        if (!sites[n].occupied[k])
          continue;
        if (tumble(rng) < S.time() - sites[n].hoptime[k])
        {
          // At least one tumble event happened since we last attempted a hop; so randomise the direction
          sites[n].direction[k] = anyway(rng);
        }
        auto dnbs = neighbours(n);
        if (sites[dnbs[sites[n].direction[k]]].present == P.n_max)
        {
          count++;
        }
      }
    }
    return count / double(P.N);
  }

  // Return cluster size distributions.
  // Element 2n   contains the number of clusters of particles of size n+1
  // Element 2n+1 contains the number of clusters of vacancies of size n+1
  hist_t cluster_distributions() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    hist_t dists(2 * maxsize);
    for (const auto &kv : clusters)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present == 0)
        dists[2 * kv.second.size() - 1]++;
      else
        dists[2 * kv.second.size() - 2]++;
    }
    return dists;
  }

  hist_t cluster_distributions_particle_numbers() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    hist_t dists_nr(2 * maxsize_nr);
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present == 0)
        dists_nr[2 * kv.second.size() - 1]++;
      else
        dists_nr[2 * kv.second.size() - 2]++;
    }
    return dists_nr;
  }

  // Function to determine the size of largest cluster. Note that we will only regard particle clusters here
  // (at least so far)
  size_t max_cluster_size()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    std::size_t max_s = 1;
    for (const auto &kv : clusters)
    {

      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        // std::cout << "-------------------" << endl;
        // std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
        max_s = std::max(max_s, kv.second.size());
        // std::cout << max_s << endl;
      }
    }

    // std::cout << max_s << endl;

    return max_s;
  }

  size_t max_cluster_size_nr()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {

          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    std::size_t max_s = 1;
    for (const auto &kv : clusters_nr)
    {

      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        // std::cout << "-------------------" << endl;
        // std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
        max_s = std::max(max_s, kv.second.size());
        // std::cout << max_s << endl;
      }
    }

    // std::cout << max_s << endl;

    return max_s;
  }

  // Function to determine the mean cluster size by simple taking the sum off all clusters (not vacant ones) and dividing by number of clusters. Note that we will only regard particle clusters here
  // (at least so far)
  double avg_cluster_size()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // variable for mean and the count of clusters
    double mean = 0;
    double count = 0;
    for (const auto &kv : clusters)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        mean += kv.second.size();
        count += 1;
      }
    }

    mean = mean / count;

    return mean;
  }

  double avg_cluster_size_nr()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {

          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    // variable for mean and the count of clusters
    double mean = 0;
    double count = 0;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        mean += kv.second.size();
        count += 1;
      }
    }

    mean = mean / count;

    return mean;
  }

  // function to return the number of particle clusters
  unsigned number_cluster()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (sites[n].present == 0)
        clusters_nr[n] = std::list<unsigned>(1, n);
      else
        clusters_nr[n] = std::list<unsigned>(sites[n].present, n);
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof_nr[n], small = memberof_nr[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {

          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters_nr[large].size() < clusters_nr[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters_nr[small])
            memberof_nr[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
          // Remove the smaller cluster from the map
          clusters_nr.erase(small);
          // Keep track of the largest cluster
          maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
        }
      }
    }

    unsigned count = 0;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[kv.first].present > 0)
      {
        count += 1;
      }
    }

    return count;
  }

  vec surface_volume()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // vector for output of surface and volume of clusters
    vec output;

    for (const auto &kv : clusters)
    {
      // variable for surface and check
      unsigned surf = 0;
      unsigned check = 0;
      // check whether particle cluster
      if (sites[kv.first].present > 0)
      {
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours(n))
          {
            if (sites[m].present == 0)
            {
              check++;
              continue;
            }
          }
          if (check > 0)
            surf++;
        }

        output.push_back(kv.second.size());
        output.push_back(surf);
      }
    }

    return output;
  }

  vec cluster_surface()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // vector for output of surface and volume of clusters
    vec output;

    for (const auto &kv : clusters)
    {
      // variable for surface and check
      unsigned check = 0;
      // check whether particle cluster
      if (sites[kv.first].present > 0)
      {
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours(n))
          {
            if (sites[m].present == 0)
            {
              check++;
              continue;
            }
          }
          if (check > 0)
            output.push_back(n);
        }
      }
    }

    return output;
  }
  // function to give out position of a cluster in the lattice that has size between lower and upper bound
  vec single_cluster(unsigned clust_min, unsigned clust_max)
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // vector for output of surface and volume of clusters
    vec output;

    for (const auto &kv : clusters)
    {
      // check whether particle cluster
      if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max)
      {
        // to get surface and volume
        unsigned surf = 0;
        unsigned check = 0;
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours(n))
          {
            if (sites[m].present == 0)
            {
              check++;
              continue;
            }
          }
          if (check > 0)
            surf++;
        }
        // add surface and volume to beginning of output
        output.push_back(kv.second.size());
        output.push_back(surf);
        // add position of particles
        for (const auto &n : kv.second)
        {
          output.push_back(n);
        }
        break;
      }
    }

    return output;
  }

  unsigned clust_size(unsigned clust_min, unsigned clust_max)
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
      clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (const auto &m : forward_neighbours(n))
      {
        unsigned large = memberof[n], small = memberof[m];
        // continue on if they are part of the same cluster
        if (small == large)
          continue;
        // continue on if they are vacant - not vacant and vise versa
        else if (sites[n].present == 0 && sites[m].present != 0)
          continue;
        else if (sites[n].present != 0 && sites[m].present == 0)
          continue;
        else
        {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if (clusters[large].size() < clusters[small].size())
            std::swap(large, small);
          // Update the cluster number for all sites in the smaller one
          for (const auto &site : clusters[small])
            memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    // vector for output of surface and volume of clusters
    unsigned output = 0;

    for (const auto &kv : clusters)
    {
      // check whether particle cluster
      if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max)
      {
        output = 1;
      }
    }

    return output;
  }

  hist_t particle_neighbour_dist()
  {
    hist_t dist(7);
    unsigned count = 0;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      count = 0;
      if (sites[n].present > 0)
      {
        for (const auto &m : neighbours(n))
        {
          if (sites[m].present > 0)
            count++;
        }
        dist[count]++;
      }
    }

    return dist;
  }

  vec_d density()
  {
    vec_d den;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      double local = 0.5 * double(sites[n].present);
      for (const auto &m : neighbours(n))
      {
        local += 1.0 / 12.0 * double(sites[m].present);
      }
      den.push_back(local);
    }

    return den;
  }

  vec_d stopping(double t)
  {
    vec_d stopping_time;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < P.n_max; i++)
      {
        if (sites[n].occupied[i] == true)
          stopping_time.push_back(t - sites[n].last_jump[i]);
      }
    }
    return stopping_time;
  }

  vec surface_volume_nr()
  {
    std::map<unsigned, std::list<unsigned>> cluster = clusters();

    // vector for output of surface and volume of clusters
    // the index 2i is for volume and 2i+1 is for surface
    vec output;
    vec index;
    for (const auto &kv : cluster)
    {
      // variable for surface and check
      unsigned surf = 0;
      unsigned check = 0;
      // check whether particle cluster
      if (sites[kv.first].present > 0)
      {
        for (const auto &n : kv.second)
        {
          // for duplicate index
          unsigned dup = 0;
          for (const auto &m : index)
          {
            if (m == n)
              dup++;
          }
          index.push_back(n);
          check = 0;
          if (dup == 0)
          {
            for (const auto &m : neighbours(n))
            {
              if (sites[m].present == 0)
              {
                check++;
              }
            }
          }
          if (check > 0)
            surf += sites[n].present;
        }
        output.push_back(kv.second.size());
        output.push_back(surf);
      }
    }
    return output;
  }
  // function for returning perimeter order parameter (see thesis)
  double perimeter()
  {
    // get distribution and surface-volume
    hist_t dist = cluster_distributions_particle_numbers();
    vec sv = surface_volume_nr();

    double second_moment = 0;
    double first_moment = 0;

    // calculate appropriate stuff
    for (unsigned i = 0; i < dist.size(); i += 2)
    {
      for (unsigned k = 0; k < sv.size(); k += 2)
      {
        if (sv[k] == ((i + 2) / 2))
        {
          second_moment += ((i + 2) / 2) * sv[k + 1];
        }
      }
      first_moment += dist[i] * ((i + 2) / 2);
    }
    return second_moment / (first_moment * double(P.N));
  }

  // function for returning w(_N) order parameter (see thesis)
  double w_N()
  {
    // get distribution and surface-volume
    hist_t dist = cluster_distributions_particle_numbers();

    double second_moment = 0;
    double first_moment = 0;

    // calculate appropriate stuff
    for (unsigned i = 0; i < dist.size(); i += 2)
    {
      second_moment += dist[i] * ((i + 2) / 2) * ((i + 2) / 2);
      first_moment += dist[i] * ((i + 2) / 2);
    }
    return second_moment / first_moment;
  }

  vec occ_array()
  {
    vec output;
    for (unsigned i = 0; i < sites.size(); i++)
    {
      output.push_back(sites[i].present);
    }
    return output;
  }

  // x coordinate of tagged particle. Double so I can copy paste to Tri and Hex ^^
  double x_coor()
  {
    double x = 0;
    for (int n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < P.n_max; i++)
      {
        if (sites[n].occupied[i] && sites[n].tagged[i])
          x = double(n % P.L[0]) + sites[n].per_hop[i] * int(P.L[0]);
      }
    }
    return x;
  }

  hist_t perimeter_dist()
  { 
    hist_t clust_size = cluster_distributions_particle_numbers();
    vec sv = surface_volume_nr();
    hist_t dists_nr(clust_size.size());

    for (unsigned k = 0; k < sv.size(); k += 2)
    {
        dists_nr[sv[k+1]]++;
    }

    return dists_nr;
  }

};

template <typename Engine>
class Hexagonal_lattice
{

  // For output in different formats
  friend class SnapshotWriter<Engine>;
  friend class VacancyWriter<Engine>;
  friend class HexagonalParticleWriter<Engine>;
  friend class ClusterWriter<Engine>;
  friend class HexDirectionWriter<Engine>;

  Parameters P; // A local copy of the model parameters
  // TODO Really need a better way of doing this
  static constexpr unsigned n_max = 10;

  // Data associated with each site; by default all of these are set to zero
  // which represents a vacant site

  // We will use the following way to identify particles and their respective positions
  // If index%2 == 0 it is the upper position, if it is index%2 == 1 it is the lower position
  struct Site
  {
    std::vector<unsigned> id = std::vector<unsigned>(2 * n_max);              // Particle / vacancy id
    std::vector<bool> occupied = std::vector<bool>(2 * n_max);                // There is a particle here
    std::vector<bool> active = std::vector<bool>(2 * n_max);                  // A move event is scheduled
    std::vector<int> neighbours = std::vector<int>(2);                        // Number of neighbours that are occupied
    std::vector<direction_t> direction = std::vector<direction_t>(2 * n_max); // Preferred direction of movement
    std::vector<double> hoptime = std::vector<double>(2 * n_max);             // Time of last hop attempt
    std::vector<int> present = std::vector<int>(2);                           // Number of particles present at each site in unit cell
    std::vector<unsigned> current_dir = std::vector<unsigned>(2 * n_max);     // Contains the lattice site it points at.
    std::vector<double> last_jump = std::vector<double>(2 * n_max);           // time of last jump made
    std::vector<double> x_coor = std::vector<double>(n_max);                  // x coordinate (unlike for square lattice we need to consider actual distance here, hence double)
  };

  /*           0        5         2
                x       x       x
                  x  1  x  0  x
                    x   x   x
                      x x x
                2     x x x    5
                      x x x
                    x   x   x
                  x  3  x  4  x
                x       x       x
              3         4         1
  */
  // explanation of the above: To ensure that particles actually move in a, somewhat, constant direction before a tumble event happens
  // we divide the directions into 6 regions (as shown above). These regions are numbered, and the directions are chosen accordingly by
  // seeing at which index of the lattice site we are at. Note that we have to places at each lattice sites, i.e. even and uneven indices.
  // the function preference_direction below choses what direction to take based on this.

  std::vector<Site> sites;                              // Representation of the sites
  Engine &rng;                                          // Source of noise: this is a reference as there should only be one of these!
  std::discrete_distribution<unsigned> anyway, initial; // Distribution over tumble directions
  std::exponential_distribution<double> run;            // Distribution of times between run and tumble events
  Scheduler S;                                          // Keeps track of the event queue
  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours. We have periodic boundary conditions

  // function to choose direction based on what site and which directional orientation it has
  auto preference_direction(unsigned n, unsigned index, direction_t direction)
  {
    int dir = 0;

    if (direction == 0)
    {
      if (index % 2 == 0)
        dir = 5;
      else
        dir = 2;
    }
    else if (direction == 1)
    {
      if (index % 2 == 0)
        dir = 5;
      else
        dir = 0;
    }
    else if (direction == 2)
    {
      if (index % 2 == 0)
        dir = 3;
      else
        dir = 0;
    }
    else if (direction == 3)
    {
      if (index % 2 == 0)
        dir = 3;
      else
        dir = 4;
    }
    else if (direction == 4)
    {
      if (index % 2 == 0)
        dir = 1;
      else
        dir = 4;
    }
    else if (direction == 5)
    {
      if (index % 2 == 0)
        dir = 1;
      else
        dir = 2;
    }
    return dir;
  }

  auto neighbours_dir(unsigned n) const
  {

    // this is only valid for 2d

    std::vector<unsigned> nbs(6);
    int L1 = P.L[0];
    int L2 = P.L[1];

    unsigned x = n % L1;
    unsigned y = (n / L1);

    int left = 0;

    if ((int(n) - 1) % L1 < 0)
      left = L1 + (int(n) - 1) % L1 + y * L1; // left
    else if ((int(n) - 1) % L1 + y * L1 >= 0)
      left = (int(n) - 1) % L1 + y * L1;

    nbs[0] = n;                            // same site
    nbs[1] = n;                            // same site
    nbs[2] = (n + 1) % L1 + y * L1;        // right
    nbs[3] = left;                         // left
    nbs[4] = x + ((y - 1 + L1) % L2) * L1; // down
    nbs[5] = x + ((y + 1) % L2) * L1;      // up
    return nbs;
  }

  auto neighbours(unsigned n, unsigned index) const
  {

    // this is only valid for 2d

    std::vector<int> nbs(3);
    int L1 = P.L[0];
    int L2 = P.L[1];

    int x = n % L1;
    int y = (n / L1);

    if (index)
    {
      nbs[0] = n;                            // same site
      nbs[1] = (n + 1) % L1 + y * L1;        // right
      nbs[2] = x + ((y - 1 + L1) % L2) * L1; // down
    }
    else
    {
      nbs[0] = n;
      if ((int(n) - 1) % L1 < 0)
        nbs[1] = L1 + (int(n) - 1) % L1 + y * L1; // left
      else if ((int(n) - 1) % L1 + y * L1 >= 0)
        nbs[1] = (int(n) - 1) % L1 + y * L1; // left                            // same site
      // nbs[1] = (n - 1) % L1 + y * L1;
      nbs[2] = x + ((y + 1) % L2) * L1; // up
    }
    return nbs;
  }

  auto forward_neighbours(unsigned n, unsigned index) const
  {

    // this is only valid for 2d
    std::vector<unsigned> nbs(2);

    int L1 = P.L[0];
    int L2 = P.L[1];

    unsigned x = n % L1;
    unsigned y = (n / L1);

    if (index)
    {
      nbs[0] = n; // same site
      nbs[0] = n; // TODO this is not ideal, better way?
    }
    else
    {
      nbs[1] = (n - 1) % L1 + y * L1;   // left
      nbs[2] = x + ((y + 1) % L2) * L1; // up
    }

    return nbs;
  }

  // function to calculate x-and-y coordinate
  vec_d coordinates(unsigned n, unsigned index){
    n = int(n); // don't want trouble in case negative
    index = int(index);

    // final x and y coordinates
    double xf, yf;

    // define specific  cos and sin
    double c30 = cos(M_PI / 6.0);
    double s30 = sin(M_PI / 6.0);

    // intermediate x and y
    int y = n / P.L[0];
    int x = 2 * (n%P.L[0]) + index - y + 1;

    xf = x * c30;

    if (y%2 == 1){
      if (x%2 == 1) yf = y * (1 + s30);
      else yf = (3 * y + 1)/2.0;
    }
    else {
      if (x%2 == 1) yf = (3 * y + 1)/2.0;
      else yf = 1.5 * y;
    }
    
    vec_d coord = {xf, yf};
    return coord;
  }

  // Place a particle with given direction and hop time at site n;
  // neighbouring sites will be accordingly adjusted
  void place(unsigned n, unsigned id, direction_t d, double t, unsigned index, unsigned dir, double x)
  {
    sites[n].id[index] = id;
    sites[n].direction[index] = d;
    sites[n].hoptime[index] = t;
    sites[n].last_jump[index] = t;
    if (!sites[n].occupied[index])
    {
      sites[n].occupied[index] = true;
      for (const auto &m : neighbours(n, index % 2))
      {
        ++sites[m].neighbours[(index + 1) % 2];
      }
    }
    sites[n].present[index % 2]++;
    sites[n].current_dir[index] = dir;
    sites[n].x_coor[index] = x;
  }

  // Schedule a hop event for a particle at site n
  void schedule(unsigned n, unsigned index)
  {
    assert(sites[n].occupied[index]);
    planned_moves++;
    S.schedule(run(rng), [this, n, index]()
              {
            assert(sites[n].occupied[index]);
            assert(sites[n].active[index]);
            // If there are no local vacancies, mark this particle as inactive and exit
            if (sites[n].neighbours[index%2] == 3 * P.n_max) {
                // std::cout << "Can't move from "; decode(n); std::cout << " deactivating" << std::endl;
                sites[n].active[index] = false;
            }


            else {
                // if(std::uniform_real_distribution<double>()(rng)>=std::exp(-P.alpha*(S.time()-sites[n].hoptime))) {
                // Make sure that for alpha=0 we dont get some weird changes.
                // ! Note that for lower alpha this threshold has to be adjusted
                if (P.alpha[0] - 1e-5 > 0){
                  if (tumble(rng) < S.time() - sites[n].hoptime[index]) {
                    sites[n].direction[index] = anyway(rng);
                  }
                }
                sites[n].hoptime[index] = S.time();
                // Get the sites adjacent to the departure site
                int dir = preference_direction(n, index, sites[n].direction[index]);

                auto dnbs = neighbours_dir(n);
                
                if (sites[dnbs[dir]].present[(index + 1)%2] < P.n_max){
                  actual_moves++;
                  unsigned ind = 0; // index for target space in direction site
                  for (unsigned k = 1 - index%2; k < 2 * P.n_max; k+=2){
                    if (!sites[dnbs[dir]].occupied[k]){
                      ind = k;
                      break;
                    }
                  }
                  assert(!sites[dnbs[dir]].active[ind]);
                  // Get the id of the vacancy that is being displaced
                  unsigned vid = sites[dnbs[dir]].id[ind];
                  // Deactive the departure site; also mark it empty
                  sites[n].occupied[index] = sites[n].active[index] = false;
                  // Place a particle on the target site; it has the same direction and hoptime as the departing particle
                  int dir_next = preference_direction(dnbs[dir], ind, sites[n].direction[index]);
                  auto dnbs_next = neighbours_dir(dnbs[dir]);

                  // adjust position of particle
                  vec_d xold = coordinates(n, index%2);
                  vec_d xnew = coordinates(dnbs[dir], (index+1)%2);

                  sites[n].x_coor[n] += xnew[0] - xold[0];

                  place(dnbs[dir], sites[n].id[index], sites[n].direction[index], sites[n].hoptime[index], ind, dnbs_next[dir_next], sites[n].x_coor[n]);
                  // Move the vacancy id onto the departure site
                  sites[n].id[index] = vid;
                  sites[n].present[index%2] -= 1;
                  // Now go through the neighbours of the departure site, update neighbour count and activate any
                  // that can now move. Note the particle this is at the target site is included in this list
                  // and will be activated accordingly
                  for (const auto& m : neighbours(n, index%2)) {
                      --sites[m].neighbours[(index+1)%2];
                      for (unsigned k = 1 - index%2; k < 2 * P.n_max; k+=2){
                        if (sites[m].occupied[k] && !sites[m].active[k]) schedule(m, k);
                      }
                  }
                  
                }
                else {
                    // std::cout << "Didn't move from "; decode(n); std::cout << std::endl;
                    // This site is still active, so schedule another hop
                    schedule(n, index);
                }
            } });
    sites[n].active[index] = true;
    // std::cout << "Scheduled "; decode(n); std::cout << std::endl;
  }

  // For testing
  void decode(unsigned n)
  {
    std::cout << "[ ";
    for (const auto &L : P.L)
    {
      std::cout << (n % L) << " ";
      n /= L;
    }
    std::cout << "]";
  }

  // Check the lattice state is consistent
  bool consistent()
  {
    unsigned active = 0, occupied = 0;
    std::set<unsigned> ids;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Check each site has a unique id
      for (unsigned i = 0; i < 2 * P.n_max; i++)
      {
        if (ids.count(sites[n].id[i]))
        {
          return false;
        }
        ids.insert(sites[n].id[i]);
        // Check that empty sites are also inactive
        if (!sites[n].occupied[i])
        {
          if (sites[n].active[i])
          {
            return false;
          }
          // Nothing left to do if empty
          continue;
        }

        // Check that the neighbour count is correct
        occupied++;
        unsigned nbs = 0;
        for (const auto &m : neighbours(n, i % 2))
        {
          for (unsigned k = 1 - i % 2; k < 2 * P.n_max; k += 2)
          {
            if (sites[m].occupied[k])
              ++nbs;
          }
        }
        if (nbs != sites[n].neighbours[i % 2])
        {
          return false;
        }
        // Check that mobile particles are active
        if (nbs < 3 * P.n_max && !sites[n].active[i])
        {
          return false;
        }
        if (sites[n].active[i])
          ++active;
      }
    }
    // Check we've not lost any particles
    return occupied == P.N && active == S.pending();
  }

public:
  unsigned planned_moves;
  unsigned actual_moves;
  std::exponential_distribution<double> tumble;

  Hexagonal_lattice(const Parameters &P, Engine &rng) : P(P),                                                                           // NB: this takes a copy
                                                        sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
                                                        rng(rng),                                                                       // NB: this takes a reference
                                                        run(1),
                                                        tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
  {
    // Set up the tumble direction distribution
    std::vector<double> drates(6);
    std::vector<double> drates_initial(6);
    for (unsigned d = 0; d < 3; ++d)
    {
      drates[2 * d] = drates[2 * d + 1] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
      drates_initial[2 * d] = drates_initial[2 * d + 1] = 1.0;
    }
    anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());
    initial = std::discrete_distribution<unsigned>(drates_initial.begin(), drates_initial.end());

    std::vector<unsigned> position;
    for (unsigned i = 0; i < 2 * P.n_max; i++)
    {
      for (unsigned n = 0; n < sites.size(); n++)
      {
        position.push_back(i * sites.size() + n);
      }
    }

    planned_moves = 0;
    actual_moves = 0;

    unsigned possibilities = position.size();
    unsigned unplaced = P.N;
    unsigned id = 0;
    while (unplaced > 0)
    {
      unsigned index = std::uniform_int_distribution<unsigned>(0, possibilities - 1)(rng);

      unsigned l = position[index];

      unsigned n = l % sites.size();
      unsigned i = l / sites.size();

      direction_t direction = initial(rng);
      int dir = preference_direction(n, i, direction);
      auto dnbs = neighbours_dir(n);
      vec_d coord = coordinates(n, i%2);
      place(n, id, direction, 0.0, i, dnbs[dir], coord[0]);

      position.erase(position.begin() + index);
      position.push_back(l);

      id++;
      possibilities--;
      unplaced--;
    }

    for (unsigned index = 0; index < 2 * P.n_max; index++)
    {
      for (unsigned n = 0; n < sites.size(); ++n)
      {
        if (sites[n].occupied[index] == false)
        {
          sites[n].id[index] = id;
          id++;
        }
      }
    }

    // Activate particles that can move, and schedule a hop accordingly
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2 * P.n_max; i++)
      {
        if (sites[n].occupied[i] && sites[n].neighbours[i % 2] < 3 * P.n_max)
        {
          schedule(n, i);
          sites[n].active[i] = true;
        }
      }
    }

    assert(consistent());
  }

  // Iterates the simulation until the given time; returns the actual time run to
  double run_until(double time)
  {
    while (S.advance(time))
      ;
    assert(consistent());
    return S.time();
  }

  void set_new_lambda(std::exponential_distribution<double> *exp_dis, double val)
  {
    typename std::exponential_distribution<double>::param_type new_lambda(val);
    exp_dis->param(new_lambda);
  }
  // Sample all particle directions from the distribution that applies at the current instant
  // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
  void realise_directions()
  {
    for (auto &site : sites)
    {
      for (unsigned i = 0; i < 2 * P.n_max; i++)
      {
        if (!site.occupied[i])
          continue;
        if (tumble(rng) < S.time() - site.hoptime[i])
        {
          // At least one tumble event happened since we last attempted a hop; so randomise the direction
          site.direction[i] = anyway(rng);
        }
        // Advance the hop time so we don't generate temporal paradox
        site.hoptime[i] = S.time();
      }
    }
  }

  // function that calculates the fraction of immobile particles
  double motility_fraction()
  {
    double count = 0;
    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < 2 * P.n_max; i++)
      {
        if (!sites[n].occupied[i])
          continue;
        if (tumble(rng) < S.time() - sites[n].hoptime[i])
        {
          sites[n].direction[i] = anyway(rng);
        }
        sites[n].hoptime[i] = S.time();
        // Get the sites adjacent to the departure site
        int dir = preference_direction(n, i, sites[n].direction[i]);
        auto dnbs = neighbours_dir(n);
        if (sites[dnbs[dir]].present[(i + 1) % 2] == P.n_max)
        {
          count++;
        }
      }
    }
    return count / double(P.N);
  }

  // Return cluster size distributions.
  // This will only calculate cluster size by area, not by particles contained
  // TODO use if present == 0 to distinguish between area sites
  hist_t cluster_distributions() const
  {

    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }

    hist_t dists(2 * maxsize);
    for (const auto &kv : clusters)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] == 0)
        dists[2 * kv.second.size() - 1]++;
      else
        dists[2 * kv.second.size() - 2]++;
    }
    return dists;
  }

  // Maybe find a better name, this is very long....
  // Important note here: I am considering a site with at least one particle present as a
  // particle cluster. This means that the overall number of id's is not conserved.
  // For this, consider n particles and n_max=2. If all n are on different site, we lose n vacant id's
  // If, on the other hand, the n are on n/2 sites we lose 0 vacant id's
  hist_t cluster_distribution_particle_number() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        if (sites[n].present[i] == 0)
          clusters_nr[2 * n + i] = std::list<unsigned>(1, 2 * n + i);
        else
        {
          clusters_nr[2 * n + i] = std::list<unsigned>(sites[n].present[i], 2 * n + i);
          /*
          for (unsigned k = 0; k < sites[n].present[i]; k++){
            pres[k] = 2*n + i;
          }
          clusters_nr[2*n + i] = pres;
          */
        }
      }
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof_nr[2 * n + i], small = memberof_nr[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters_nr[large].size() < clusters_nr[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters_nr[small])
              memberof_nr[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
            // Remove the smaller cluster from the map
            clusters_nr.erase(small);
            // Keep track of the largest cluster
            maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
          }
        }
      }
    }

    hist_t dists_nr(2 * maxsize_nr);
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] == 0)
        dists_nr[2 * kv.second.size() - 1]++;
      else
        dists_nr[2 * kv.second.size() - 2]++;
    }
    return dists_nr;
  }

  map<unsigned, std::list<unsigned>> clusters() const
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        if (sites[n].present[i] == 0)
          clusters_nr[2 * n + i] = std::list<unsigned>(1, 2 * n + i);
        else
        {
          clusters_nr[2 * n + i] = std::list<unsigned>(sites[n].present[i], 2 * n + i);
          /*
          for (unsigned k = 0; k < sites[n].present[i]; k++){
            pres[k] = 2*n + i;
          }
          clusters_nr[2*n + i] = pres;
          */
        }
      }
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof_nr[2 * n + i], small = memberof_nr[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters_nr[large].size() < clusters_nr[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters_nr[small])
              memberof_nr[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
            // Remove the smaller cluster from the map
            clusters_nr.erase(small);
            // Keep track of the largest cluster
            maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
          }
        }
      }
    }

    return clusters_nr;
  }

  size_t max_cluster_size()
  {

    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }
    std::size_t max_s = 1;
    for (const auto &kv : clusters)
    {

      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        // std::cout << "-------------------" << endl;
        // std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
        max_s = std::max(max_s, kv.second.size());
        // std::cout << max_s << endl;
      }
    }
    // std::cout << max_s << endl;
    return max_s;
  }

  size_t max_cluster_size_nr()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        if (sites[n].present[i] == 0)
          clusters_nr[2 * n + i] = std::list<unsigned>(1, 2 * n + i);
        else
        {
          clusters_nr[2 * n + i] = std::list<unsigned>(sites[n].present[i], 2 * n + i);
          /*
          for (unsigned k = 0; k < sites[n].present[i]; k++){
            pres[k] = 2*n + i;
          }
          clusters_nr[2*n + i] = pres;
          */
        }
      }
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof_nr[2 * n + i], small = memberof_nr[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters_nr[large].size() < clusters_nr[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters_nr[small])
              memberof_nr[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
            // Remove the smaller cluster from the map
            clusters_nr.erase(small);
            // Keep track of the largest cluster
            maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
          }
        }
      }
    }
    std::size_t max_s = 1;
    for (const auto &kv : clusters_nr)
    {

      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        // std::cout << "-------------------" << endl;
        // std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
        max_s = std::max(max_s, kv.second.size());
        // std::cout << max_s << endl;
      }
    }

    // std::cout << max_s << endl;

    return max_s;
  }

  // Function to determine the mean cluster size by simple taking the sum off all clusters (not vacant ones) and dividing by number of clusters. Note that we will only regard particle clusters here
  // (at least so far)
  double avg_cluster_size()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }

    // variable for mean and the count of clusters
    double mean = 0;
    double count = 0;
    for (const auto &kv : clusters)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        mean += kv.second.size();
        count += 1;
      }
    }

    mean = mean / count;

    return mean;
  }

  double avg_cluster_size_nr()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        if (sites[n].present[i] == 0)
          clusters_nr[2 * n + i] = std::list<unsigned>(1, 2 * n + i);
        else
        {
          clusters_nr[2 * n + i] = std::list<unsigned>(sites[n].present[i], 2 * n + i);
          /*
          for (unsigned k = 0; k < sites[n].present[i]; k++){
            pres[k] = 2*n + i;
          }
          clusters_nr[2*n + i] = pres;
          */
        }
      }
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof_nr[2 * n + i], small = memberof_nr[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters_nr[large].size() < clusters_nr[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters_nr[small])
              memberof_nr[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
            // Remove the smaller cluster from the map
            clusters_nr.erase(small);
            // Keep track of the largest cluster
            maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
          }
        }
      }
    }

    // variable for mean and the count of clusters
    double mean = 0;
    double count = 0;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        mean += kv.second.size();
        count += 1;
      }
    }

    mean = mean / count;
    return mean;
  }

  // function to return the number of particle clusters
  unsigned number_cluster()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof_nr(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters_nr;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        if (sites[n].present[i] == 0)
          clusters_nr[2 * n + i] = std::list<unsigned>(1, 2 * n + i);
        else
        {
          clusters_nr[2 * n + i] = std::list<unsigned>(sites[n].present[i], 2 * n + i);
          /*
          for (unsigned k = 0; k < sites[n].present[i]; k++){
            pres[k] = 2*n + i;
          }
          clusters_nr[2*n + i] = pres;
          */
        }
      }
    }
    // Keep track of the size of the largest cluster
    std::size_t maxsize_nr = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof_nr[2 * n + i], small = memberof_nr[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters_nr[large].size() < clusters_nr[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters_nr[small])
              memberof_nr[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters_nr[large].splice(clusters_nr[large].end(), clusters_nr[small]);
            // Remove the smaller cluster from the map
            clusters_nr.erase(small);
            // Keep track of the largest cluster
            maxsize_nr = std::max(maxsize_nr, clusters_nr[large].size());
          }
        }
      }
    }

    // variable for mean and the count of clusters
    unsigned count = 0;
    for (const auto &kv : clusters_nr)
    {
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        count += 1;
      }
    }

    return count;
  }

  vec surface_volume()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }

    // variable for mean and the count of clusters
    vec output;
    for (const auto &kv : clusters)
    {
      unsigned surf = 0;
      unsigned check = 0;
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours((n - n % 2) / 2, n % 2))
          {
            if (sites[m].present[(n + 1) % 2] == 0)
              check++;
          }
          if (check > 0)
            surf++;
        }
        output.push_back(kv.second.size());
        output.push_back(surf);
      }
    }

    return output;
  }
  vec surface_volume_nr()
  {
    // vector for output of surface and volume of clusters
    // the index 2i is for volume and 2i+1 is for surface
    std::map<unsigned, std::list<unsigned>> cluster = clusters();
    vec output;
    vec index;
    for (const auto &kv : cluster)
    {
      unsigned surf = 0;
      unsigned check = 0;
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        for (const auto &n : kv.second)
        {
          unsigned dup = 0;
          for (const auto &m : index)
          {
            if (m == (n - n % 2) / 2)
              dup++;
          }
          index.push_back((n - n % 2) / 2);
          check = 0;
          if (dup == 0)
          {
            for (const auto &m : neighbours((n - n % 2) / 2, n % 2))
            {
              if (sites[m].present[(n + 1) % 2] == 0)
                check++;
            }
          }
          if (check > 0)
            surf += sites[(n - n % 2) / 2].present[n % 2];
        }
        output.push_back(kv.second.size());
        output.push_back(surf);
      }
    }

    return output;
  }
  vec cluster_surface()
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }

    // variable for mean and the count of clusters
    vec output;
    for (const auto &kv : clusters)
    {
      unsigned check = 0;
      // instead of occupied we check whether there are any particles present
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0)
      {
        for (const auto &n : kv.second)
        {
          check = 0;
          for (const auto &m : neighbours((n - n % 2) / 2, n % 2))
          {
            if (sites[m].present[(n + 1) % 2] == 0)
              check++;
          }
          if (check > 0)
          {
            output.push_back((n - n % 2) / 2);
            output.push_back(n % 2);
          }
        }
      }
    }

    return output;
  }
  // function to give out position of a cluster in the lattice that has size between lower and upper bound
  vec single_cluster(unsigned clust_min, unsigned clust_max)
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }

    // vector for output of surface and volume of clusters
    vec output;
    std::cout << "test" << endl;

    for (const auto &kv : clusters)
    {
      // check whether particle cluster
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max)
      {
        // to get surface and volume
        unsigned surf = 0;
        unsigned check = 0;
        for (const auto &n : kv.second)
        {
          if (sites[(n - n % 2) / 2].present[n % 2] > 0)
          {
            check = 0;
            for (const auto &m : neighbours((n - n % 2) / 2, n % 2))
            {
              if (sites[m].present[(n + 1) % 2] == 0)
              {
                check++;
                continue;
              }
            }
            if (check > 0)
              surf++;
          }
        }
        std::cout << "break" << endl;
        // add surface and volume to beginning of output
        output.push_back(kv.second.size());
        output.push_back(surf);
        // add position of particles
        for (const auto &n : kv.second)
        {
          if (sites[(n - n % 2) / 2].present[n % 2] > 0)
          {
            output.push_back((n - n % 2) / 2);
            output.push_back(n % 2);
          }
        }
        break;
      }
    }

    return output;
  }

  unsigned clust_size(unsigned clust_min, unsigned clust_max)
  {
    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(2 * sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        clusters[2 * n + i] = std::list<unsigned>(1, 2 * n + i); // Single-element list comprising the lattice site
      }
    }

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // Loop over neigbours m in one direction only so we visit each bond once
      for (unsigned i = 0; i < 2; i++)
      {
        // TODO figure out why forward neighbour makes algorithm not work
        for (const auto &m : neighbours(n, i))
        {
          unsigned large = memberof[2 * n + i], small = memberof[2 * m + (i + 1) % 2];
          // If they are in the same cluster we can move on
          if (small == large)
            continue;
          // If one of them is empty but the other isn't we move on
          else if (sites[n].present[i] == 0 && sites[m].present[(i + 1) % 2] != 0)
            continue;
          else if (sites[n].present[i] != 0 && sites[m].present[(i + 1) % 2] == 0)
            continue;
          // merge clusters
          else
          {
            // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
            if (clusters[large].size() < clusters[small].size())
              std::swap(large, small);
            // Update the cluster number for all sites in the smaller one
            for (const auto &site : clusters[small])
              memberof[site] = large;
            // Add the members of the smaller cluster onto the end of the larger one
            clusters[large].splice(clusters[large].end(), clusters[small]);
            // Remove the smaller cluster from the map
            clusters.erase(small);
            // Keep track of the largest cluster
            maxsize = std::max(maxsize, clusters[large].size());
          }
        }
      }
    }

    // vector for output of surface and volume of clusters
    unsigned output = 0;

    for (const auto &kv : clusters)
    {
      // check whether particle cluster
      if (sites[(kv.first - kv.first % 2) / 2].present[kv.first % 2] > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max)
      {
        output = 1;
      }
    }

    return output;
  }

  hist_t particle_neighbour_dist()
  {
    hist_t dist(4);
    unsigned count = 0;
    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        count = 0;
        if (sites[n].present[i] > 0)
        {
          for (const auto &m : neighbours(n, i))
          {
            if (sites[m].present[(i + 1) % 2] > 0)
              count++;
          }
          dist[count]++;
        }
      }
    }
    return dist;
  }

  vec_d density()
  {
    vec_d den;
    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        double local = 0.5 * double(sites[n].present[i]);
        for (const auto &m : neighbours(n, i))
        {
          local += 1.0 / 6.0 * double(sites[m].present[(i + 1) % 2]);
        }
        den.push_back(local);
      }
    }
    return den;
  }

  vec_d stopping(double t)
  {
    vec_d stopping_time;

    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < 2 * P.n_max; i++)
      {
        if (sites[n].occupied[i] == true)
          stopping_time.push_back(t - sites[n].last_jump[i]);
      }
    }
    return stopping_time;
  }

  // function for returning perimeter order parameter (see thesis)
  double perimeter()
  {
    // get distribution and surface-volume
    hist_t dist = cluster_distribution_particle_number();

    vec sv = surface_volume_nr();

    double second_moment = 0;
    double first_moment = 0;

    // calculate appropriate stuff
    for (unsigned i = 0; i < dist.size(); i += 2)
    {
      for (unsigned k = 0; k < sv.size(); k += 2)
      {
        if (sv[k] == ((i + 2) / 2))
        {

          second_moment += ((i + 2) / 2) * sv[k + 1];
        }
      }
      first_moment += dist[i] * ((i + 2) / 2);
    }

    return second_moment / (first_moment * double(P.N));
  }

  vec perimeter_dist()
  {
    vec sv = surface_volume_nr();
    vec dist;
    for (unsigned k = 0; k < sv.size(); k += 2)
    {
      if (sv[k] > 0) dist.push_back(sv[k+1]);
    }

    return dist;
  }

  // function for returning w(_N) order parameter (see thesis)
  double w_N()
  {
    // get distribution and surface-volume
    hist_t dist = cluster_distribution_particle_number();

    double second_moment = 0;
    double first_moment = 0;

    // calculate appropriate stuff
    for (unsigned i = 0; i < dist.size(); i += 2)
    {
      second_moment += dist[i] * ((i + 2) / 2) * ((i + 2) / 2);
      first_moment += dist[i] * ((i + 2) / 2);
    }
    return second_moment / first_moment;
  }
  vec occ_array()
  {
    vec output;
    for (unsigned n = 0; n < sites.size(); n++)
    {
      for (unsigned i = 0; i < 2; i++)
      {
        output.push_back(sites[n].present[i]);
      }
    }
    return output;
  }
};

// TODO not only 0 component
template <typename Engine>
class SnapshotWriter
{
  const Lattice<Engine> &L;

public:
  SnapshotWriter(const Lattice<Engine> &L) : L(L) {}

  // Output pairs (d, a) indicating the direction d (numbered from 1) or 0 if site is vacant
  // and a indicating if the site is active

  friend std::ostream &operator<<(std::ostream &out, const SnapshotWriter &SW)
  {

    for (const auto &site : SW.L.sites)
    {
      out << (site.occupied[0] ? unsigned(site.direction[0] + 1) : 0);
      out << " " << site.active[0] << " ";
    }

    return out;
  }
};

// just a small output to see if the directions are properly implemented
template <typename Engine>
class HexDirectionWriter
{
  const Hexagonal_lattice<Engine> &L;

public:
  HexDirectionWriter(const Hexagonal_lattice<Engine> &L, ofstream &outfile) : L(L) {}

  // Output tuple (n, j, d), with n being lattice site, j being index in lattice site and
  // d the lattice site where particle is pointing

  friend std::ostream &operator<<(std::ostream &out, const HexDirectionWriter &HSW)
  {

    const auto &sites = HSW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {

      for (unsigned j = 0; j < 2 * 2; j++)
      {
        if (sites[n].occupied[j])
          out << n << "\t" << j << "\t" << sites[n].current_dir[j] << "\t" << int(sites[n].direction[j]);
      }
    }

    return out;
  }
};

template <typename Engine>
class VacancyWriter
{
  const Lattice<Engine> &L;

public:
  VacancyWriter(const Lattice<Engine> &L) : L(L) {}

  // Output pairs (v, n) indicating the vacancy id and the site on which it is found

  friend std::ostream &operator<<(std::ostream &out, const VacancyWriter &VW)
  {

    const auto &sites = VW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      if (!sites[n].occupied[0])
        out << sites[n].id[0] << " " << n << " ";
    }

    return out;
  }
};

template <typename Engine>
class ParticleWriter
{
  const Lattice<Engine> &L;

public:
  ParticleWriter(const Lattice<Engine> &L, ofstream &outfile) : L(L) {}

  // Output pairs (n, p) indicating site and number of particles present

  friend std::ostream &operator<<(std::ostream &out, const ParticleWriter &PW)
  {

    const auto &sites = PW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {

      if (sites[n].present > 0)
        out << n << " " << sites[n].present << " ";
    }

    return out;
  }
};

template <typename Engine>
class ParticleWriterDetails
{
  const Lattice<Engine> &L;

public:
  ParticleWriterDetails(const Lattice<Engine> &L, ofstream &outfile) : L(L) {}

  // Output pairs (n, p, c) indicating site, number of particles present, and if it is on circumference of cluster


  friend std::ostream &operator<<(std::ostream &out, const ParticleWriterDetails &PW)
  {

    const auto &sites = PW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      int c = 0;
      if (sites[n].present > 0){
        std::vector<unsigned> nbs(4);
        unsigned below = 1;
      for (unsigned d = 0; d < 4; ++d)
      {
        unsigned L = PW.L.P.L[0];
        unsigned above = below * L;
        // x is the position along dimension d
        // y is the contribution to the site index from all other dimensions
        unsigned x = (n / below) % L;
        unsigned y = (n % below) + (n / above) * above;
        // Neighbours in the increasing and decreasing directions
        nbs[2 * d] = y + ((x + 1) % L) * below;
        nbs[2 * d + 1] = y + ((x + L - 1) % L) * below;
        below = above;
      }
        for (const auto &m : nbs){
          if (sites[m].present == 0) c = 1;
        }

        out << n << " " << sites[n].present << " " << c << " ";
      }
    }

    return out;
  }
};

template <typename Engine>
class TriangleParticleWriter
{
  const Triangle_lattice<Engine> &L;

public:
  TriangleParticleWriter(const Triangle_lattice<Engine> &L, ofstream &outfile) : L(L) {}

  // Output pairs (v, n, p) indicating the particle id and the site on which it is found. p is the number of particles
  // present on site n
  friend std::ostream &operator<<(std::ostream &out, const TriangleParticleWriter &PW)
  {

    const auto &sites = PW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      // do I really need id here?
      if (sites[n].present > 0)
        out << sites[n].id[0] << " " << n << " " << sites[n].present << " ";
    }

    return out;
  }
};


template <typename Engine>
class HexagonalParticleWriter
{
  const Hexagonal_lattice<Engine> &L;
  // Parameters P;
public:
  HexagonalParticleWriter(const Hexagonal_lattice<Engine> &L, ofstream &outfile) : L(L) {}

  // Output is (v, n, j, p) as above with j being the binary site argument and p being number of
  // particles present
  friend std::ostream &operator<<(std::ostream &out, const HexagonalParticleWriter &HPW)
  {

    const auto &sites = HPW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned j = 0; j < 2; j++)
      {
        if (sites[n].present[j] > 0)
          out << sites[n].id[j] << "\t" << n << "\t" << j << "\t" << sites[n].present[j] << "\t";
      }
    }

    return out;
  }
  /*
    const auto& sites = HPW.L.sites; // Save typing

    unsigned L = P.L[0];

    double c30 = std::cos(pi/6.0);
    double s30 = std::sin(pi/6.0);
    double xf, yf;

    // we have to loop over sites and indices
    for(unsigned n=0; n<sites.size(); ++n) {
      for (unsigned i = 0; i < 2; i++){
        if(sites[n].occupied[i]){
        // calculating position for plotting
        // Using algorithm/code provided by Filipe Cunha Thewes
        unsigned y = n/P.L[0]; //assuming quadratic lattice, i.e. L[0] = L[1]
        unsigned x = 2 * (n%P.L[0]) + i - n/P.L[0];

        x++;
        xf = x * c30;

        if (y%2){
          if (x%2) yf = y * (1 + s30);
          else yf = (y * 3.0 + 1.0);
        }
        else {
          if (x%2) yf = 1.5 * y + 0.5;
          else yf = 1.5 * y;
        }
        // output is: id x y
        out << sites[n].id[i] << " " << xf << " " << yf << "";
        }
      }
    }

    return out;
  }
  */
};

template <typename Engine>
class ClusterWriter

{
  const Lattice<Engine> &L;

public:
  ClusterWriter(const Lattice<Engine> &L) : L(L) {}

  friend std::ostream &operator<<(std::ostream &out, const ClusterWriter &SW)
  {
    for (const auto &v : SW.L.cluster_distributions())
      out << v << " ";
    return out;
  }
};


template <typename Engine>
class ParticleWriterWID
{
  const Lattice<Engine> &L;

public:
  ParticleWriterWID(const Lattice<Engine> &L, ofstream &outfile) : L(L) {}

  // Output pairs (id, n, p) indicating site and number of particles present

  friend std::ostream &operator<<(std::ostream &out, const ParticleWriterWID &PW)
  {

    const auto &sites = PW.L.sites; // Save typing

    for (unsigned n = 0; n < sites.size(); ++n)
    {
      for (unsigned i = 0; i < PW.L.n_max; i++)
      {
        if (sites[n].occupied[i])
          out << sites[n].id[i] << " " << n << " " << sites[n].present << " ";
      }
    }

    return out;
  }
};

int main(int argc, char *argv[])
{

  // Load up default parameters
  Parameters P;

  // Set up command-line overrides
  CLI::App app{"Condensed Run-and-Tumble model (crumble)"};

  app.add_option("-L,--lengths", P.L, "Lattice lengths");
  app.add_option("-N,--particles", P.N, "Number of particles");
  app.add_option("-t,--tumble", P.alpha, "Tumble rate");
  app.add_option("-n,--occupation", P.n_max, "Max occupation of a site");
  app.add_option("-z,--tagged", P.tagged, "tagged particle");

  // Output parameters
  std::string output = "";
  std::string lattice_type = "";
  double burnin = 1000, until = 5000, every = 2.5;
  unsigned localaverage = 0;
  double localinterval = 10.0;

  unsigned details = 0; // for detailed
  unsigned region = 0; // for region of alpha in analysis, so can run parallel

  app.add_option("-o, --output", output, "Output type: snapshots|particles|vacancies|clusters");
  app.add_option("-l, --lattice_type", lattice_type, "Lattice type: square|triangular|hexagonal");

  app.add_option("-b,--burnin", burnin, "Time to run before starting measurements");
  app.add_option("-u,--until", until, "Time to run for once measurements started");
  app.add_option("-e,--every", every, "Measurement interval");
  app.add_option("-a,--localaverage", localaverage, "Number of local averages (clusters only; 0=stationary state)");
  app.add_option("-i,--localinterval", localinterval, "Interval between local averages");
  app.add_option("-d,--details", details, "for detailed output do 1, else 0 (default)");
  app.add_option("-r,--reg", region, "any integer >= 0 to name file");

  CLI11_PARSE(app, argc, argv);

  if (output[0] == 'p')
    output = "particles";
  else if (output[0] == 'v')
    output = "vacancies";
  else if (output[0] == 'c')
    output = "clusters";
  else if (output[0] == 'z')
    output = "tagged"; // for tagged particle
  else if (output[0] == 't')
    output = "stopping time"; // output for stopping time distribution
  else if (output[0] == 'd')
    output = "distribution"; // for distribution of neighbours (maybe also density?)
  else if (output[0] == 'l')
    output = "lagging"; // this is for output of hysteresis (l and lagging for greek roots of word)
  else if (output[0] == 'a')
    output = "area"; // for area/surface analysis
  else if (output[0] == 'w')
    output = "weighted"; // for weighted distribution
  else if (output[0] == 'm')
    output = "motility"; // for outputting the motility of the system
  else if (output[0] == 's')
    output = "stable"; // this is for making gifs that show how a system stabilizes...
  else if (output[0] == 'n')
    output = "number"; // this is for output of time evolution of cluster number and mean cluster size
  else if (output[0] == 'f')
    output = "function"; // this output is for testing different features. It is not static, as of now, so one should not uses this  without proper inspection of what it does
  else if (output[0] == 'h')
    output = "heatmap"; // this is for heatmap of clustersizes // ! thinking of removing this as not relevant anymore
  else if (output[0] == 'x')
    output = "perimeter"; // (first?) output for perimeter order parameter
  else if (output[0] == 'y')
    output = "correlation"; // data for overlap function, see site-site correlation in thesis // ! this needs a lot of work
  else if (output[0] == 'k')
    output = "detailed balance";
  else if (output[0] == 'q')
    output = "surface scaling";
  else
    output = "snapshots";

  if (lattice_type[0] == 's')
    lattice_type = "square";
  else if (lattice_type[0] == 't')
    lattice_type = "triangular";
  else if (lattice_type[0] == 'h')
    lattice_type = "hexagonal";

  if (localaverage > 0)
  {
    if (output != "clusters")
    {
      std::cerr << "Can only perform local averages in clusters mode" << std::endl;
      return 1;
    }
    if ((localaverage - 1) * localinterval >= every)
    {
      std::cerr << "Local averaging window larger than measurement internal" << std::endl;
      return 1;
    }
  }

  // Ouput name
  // alpha
  std::ostringstream alpha;
  alpha << std::fixed;
  alpha << std::setprecision(3);
  alpha << P.alpha[0];
  std::string alpha_p = alpha.str();
  // phi
  std::ostringstream phi;
  phi << std::fixed;
  phi << std::setprecision(2);
  phi << double(P.N) / double(P.L[0] * P.L[1]);
  std::string phi_p = phi.str();

  // occ
  std::ostringstream occ;
  occ << std::fixed;
  occ << std::setprecision(0);
  occ << P.n_max;
  std::string occ_p = occ.str();


  string txt = ".txt";
  string tumb = "alpha" + alpha_p;
  string dens = "phi" + phi_p;
  string size = "L" + std::to_string(P.L[0]);
  string path = "./lars_sim/Data/clustdist/";
  string txtoutput = path + lattice_type + "_" + tumb + "_" + dens + "_" + size + "_" + occ_p + txt;
  string txtoutput_nr = path + lattice_type + "_nr" + "_" + tumb + "_" + dens + "_" + size + "_" + occ_p + txt;

  std::cout << "# lattice = " << lattice_type << std::endl;

  // if (P.n_max > 1){
  //  Depending on what lattice, the output may be different
  if (lattice_type == "square")
  {
    // Initialise a random number generator and set up the model
    std::mt19937 rng((std::random_device())());
    Lattice L(P, rng);

    // Snapshot sequence has a header that sets out the simulation parameters
    std::cout << "# L = [ ";
    for (const auto &L : P.L)
      std::cout << L << " ";
    std::cout << "]" << std::endl;
    std::cout << "# N = " << P.N << std::endl;
    std::cout << "# alpha = [ ";
    for (const auto &alpha : P.alpha)
      std::cout << alpha << " ";
    std::cout << "]" << std::endl;
    std::cout << "# output = " << output << std::endl;
    std::cout << "# initial = " << burnin << std::endl;
    std::cout << "# interval = " << every << std::endl;
    if (localaverage > 0)
    {
      std::cout << "# localaverage = " << localaverage << std::endl;
      std::cout << "# localinterval = " << localinterval << std::endl;
    }
    std::cout << "# occupation number = " << P.n_max << std::endl;

    double t = 0;

    if (output == "clusters")
    {
      int sample_count = 0;
      double w_n = 0;
      double J = 0;
      double c_n = 0;


      if (localaverage == 0)
      {
        // hist_t sumhist;
        hist_t sumhist_nr;

        // for better statistics, we sample from multiple different
        // iterations of the independent realitsations of the lattice
        for (unsigned it = 0; it < 30; it++){
          cout << "Start \n" << endl;
          // set time to 0
          double t = 0.0;
          // initiate new lattice
          Lattice L(P, rng);
          // We sum the histograms over all measurements
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = L.run_until(burnin + n * every);
            // hist_t hist = L.cluster_distributions();
            hist_t hist_nr = L.cluster_distribution_particle_number();

            // if (hist.size() > sumhist.size())
            // {
            //   hist[std::slice(0, sumhist.size(), 1)] += sumhist;
            //   sumhist = std::move(hist);
            // }
            // else
            // {
            //   sumhist[std::slice(0, hist.size(), 1)] += hist;
            // }
            // with occ number as well

            if (hist_nr.size() > sumhist_nr.size())
            {
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            }
            else
            {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }

            J += L.motility_fraction();
            c_n += L.perimeter();
            w_n += L.w_N();
            sample_count++;
          }

          // make example snapshots
          if (it == 28){
            ofstream snapshot;
            string SnapshotName = "./lars_sim/Data/clustdist/snapshot_L"+std::to_string(P.L[0])+"_"+alpha_p+"_"+size+"_"+dens+"_"+occ_p+txt;
            snapshot.open(SnapshotName);
            snapshot << ParticleWriterDetails(L, snapshot) << endl;
          }
          cout << "Iteration " << it << " complete \n" << endl;
        }
        // output into file
        ofstream outfile_nr, outfile_order;
        // outfile.open(txtoutput);
        // for (const auto &k : sumhist)
        //   outfile << k << " ";
        // outfile << endl;
        outfile_nr.open(txtoutput_nr);
        for (const auto &k : sumhist_nr)
          outfile_nr << k << " ";
        outfile_nr << endl;

        w_n = w_n / double(P.N);

        J = J / sample_count;
        w_n = w_n / sample_count;
        c_n = c_n / sample_count;


        
        string OrderPara = "./lars_sim/Data/clustdist/order_L"+std::to_string(P.L[0])+"_"+alpha_p+"_"+size+"_"+dens+"_"+occ_p+txt;
        outfile_order.open(OrderPara);

        outfile_order << J << " " << " " << c_n << " " << w_n << " " << endl;
      }
      else
      {
        // We perform local averages of the histograms at the given measurement interval
        for (unsigned n = 0; t < burnin + until; ++n)
        {
          t = L.run_until(burnin + n * every);
          hist_t sumhist = L.cluster_distributions();
          hist_t sumhist_nr = L.cluster_distribution_particle_number();
          for (unsigned m = 1; m < localaverage; ++m)
          {
            t = L.run_until(burnin + n * every + m * localinterval);
            hist_t hist = L.cluster_distributions();
            hist_t hist_nr = L.cluster_distribution_particle_number();
            // Add hist to sumhist taking into account they may have different lengths
            if (hist.size() > sumhist.size())
            {
              hist[std::slice(0, sumhist.size(), 1)] += sumhist;
              sumhist = std::move(hist);
            }
            else
            {
              sumhist[std::slice(0, hist.size(), 1)] += hist;
            }
            if (hist_nr.size() > sumhist_nr.size())
            {
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            }
            else
            {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }
          }
          ofstream outfile, outfile_nr;
          outfile.open(txtoutput);
          for (const auto &k : sumhist)
            outfile << k << " ";
          outfile << endl;
          outfile_nr.open(txtoutput_nr);
          for (const auto &k : sumhist_nr)
            outfile_nr << k << " ";
          outfile_nr << endl;
        }
      }
    }
    else if (output == "heatmap")
    { // ! OUTDATED
      ofstream outfile_nr, outfile_avg_nr;
      outfile_nr.open("./lars_sim/heatmap/square_nr_alpha_N_n_3.txt");
      outfile_avg_nr.open("./lars_sim/heatmap/square_nr_alpha_N_n_3_avg.txt");
      for (double alp = 1e-6; alp <= 1.0; alp *= 1.3)
      {
        for (unsigned N = 100 * P.n_max; N <= P.L[0] * P.L[0] * P.n_max; N += 200 * P.n_max)
        {
          Parameters P_h;
          P_h.N = N;
          P_h.alpha[0] = P.alpha[1] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          std::size_t maxsize_nr = 1;
          Lattice LB(P_h, rng);
          t = 0;
          // also taking mean over all used timesteps
          double mean_nr = 0.0;
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            maxsize_nr = std::max(maxsize_nr, LB.max_cluster_size_nr());
            mean_nr += LB.avg_cluster_size_nr();
            count++;
          }

          mean_nr = mean_nr / count;

          outfile_nr << maxsize_nr << " ";
          outfile_avg_nr << mean_nr << " ";
        }
        outfile_nr << endl;
        outfile_avg_nr << endl;
      }
    }
    else if (output == "number")
    { // TODO substitute the old methods with the functions that you have created
      ofstream outfile, snap;
      string name2 = "./lars_sim/Data/trajectory/square_snap";
      string name = "./lars_sim/Data/trajectory/square_";
      string name_number2 = name2 + "number" + "_" + tumb + "_" + dens + "_" + size + "_" + occ_p + txt;
      string name_number = name + "number" + "_" + tumb + "_" + dens + "_" + size + "_" + occ_p + txt;
      outfile.open(name_number);
      snap.open(name_number2);
      snap << ParticleWriter(L, snap) << endl;
      vec_d res, clnr, cn, j, m;
      for (unsigned n = 0; burnin + pow(1.01, n) * every < burnin + until; n++)
      {
        res.push_back(0);
        clnr.push_back(0);
        cn.push_back(0);
        j.push_back(0);
        m.push_back(0);
      }

      unsigned iterations = 10;
      unsigned wow = 0;
      // average the curve over "iterations" many individual runs
      for (int q = 0; q <= iterations; q++)
      {
        // have to start an new, independent run for each iteration
        Lattice L(P, rng);
        t = L.run_until(burnin);
        hist_t hist = L.cluster_distributions();
        // calculate the stuff at t=0
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i += 2)
        {
          second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
          first_moment += hist[i] * ((i + 2) / 2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        res[0] += weighted;
        clnr[0] += L.avg_cluster_size_nr();
        cn[0] += L.perimeter();
        j[0] += L.motility_fraction();
        m[0] += L.max_cluster_size_nr();
        // now let time run and do the appropriate calculations
        for (unsigned n = 1; t < burnin + until; n++)
        {
          t = L.run_until(burnin + pow(1.01, n) * every);
          hist_t hist = L.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i += 2)
          {
            second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment += hist[i] * ((i + 2) / 2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          res[n] += weighted;
          clnr[n] += L.avg_cluster_size_nr();
          cn[n] += L.perimeter();
          j[n] += L.motility_fraction();
          m[n] += L.max_cluster_size_nr();
          if (wow < 2 && n == 500)
          {
            snap << ParticleWriter(L, snap) << endl;
            wow++;
          }
          else if (wow < 2 && n == res.size() - 1)
          {
            snap << ParticleWriter(L, snap) << endl;
            wow++;
          }
        }
        cout << "We are at " << q << endl;
      }
      for (unsigned i = 0; i <= res.size(); i++)
      {
        outfile << burnin + pow(1.01, i) * every << " " << res[i] / double(iterations + 1) << " " << cn[i] / double(iterations + 1) << " " << j[i] / double(iterations + 1) << " " << m[i] / double(iterations + 1) << " " << double(clnr[i]) / double(iterations + 1) << endl;
      }
    }
    else if (output == "lagging")
    {
      ofstream outfile, backward;
      string name = "./lars_sim/Data/phase/square_perc_fhyst";
      string outputname = name + "_" + occ_p + ".txt";
      outfile.open(outputname);
      name = "./lars_sim/Data/phase/square_perc_bhyst";
      outputname = name + "_" + occ_p + ".txt";
      backward.open(outputname);
      // foward hysteresis, i.e. start below critical point and move up
      Lattice LB(P, rng);
      double tmax = burnin + until;
      unsigned c = 0;

      for (double al = 0.03; al < 0.5; al *= 1.15)
      {
        // introducing new alpha
        P.alpha[0] = P.alpha[1] = al;
        LB.set_new_lambda(&LB.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
        std::vector<double> values_mot;
        std::vector<double> values_mas;
        std::vector<double> values_wei;
        std::vector<double> values_per;
        double mean = 0;
        double rel_mass = 0;
        double count = 0;
        double weighted = 0;
        double cov_w = 0;
        double s = 0;
        for (unsigned n = 0; t < (c + 1) * tmax; ++n)
        {
          t = LB.run_until(burnin + n * every + c * tmax);

          values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
          values_mot.push_back(LB.motility_fraction());
          rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
          mean += LB.motility_fraction();
          values_per.push_back(LB.perimeter() / double(P.N));
          s += LB.perimeter() / double(P.N);
          count++;
          // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

          hist_t hist = LB.cluster_distribution_particle_number();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i += 2)
          {
            second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment += hist[i] * ((i + 2) / 2);
          }
          // std::cout << second_moment << " " << first_moment << endl;
          values_wei.push_back(second_moment / first_moment * 1.0 / double(P.N));
          weighted += second_moment / first_moment * 1.0 / double(P.N);
        }
        mean = mean / count;
        rel_mass = rel_mass / count;
        weighted = weighted / count;
        s = s / count;

        double cov_mot = 0;
        for (auto &val : values_mot)
        {
          cov_mot += pow(val - mean, 2);
        }
        cov_mot = cov_mot / (values_mot.size() - 1);

        double cov_mas = 0;
        for (auto &val : values_mas)
        {
          cov_mas += pow(val - rel_mass, 2);
        }
        cov_mas = cov_mas / (values_mas.size() - 1);

        for (auto &val : values_wei)
        {
          cov_w += pow(val - weighted, 2);
        }
        cov_w = cov_w / (values_wei.size() - 1);

        double cov_s = 0;
        for (auto &val : values_per)
        {
          cov_s += pow(val - s, 2);
        }
        cov_s = cov_s / (values_per.size() - 1);

        outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << " " << s << " " << cov_s << endl;
        c++;
      }

      Lattice LT(P, rng);
      c = 0;
      t = 0;
      for (double al = 0.5; al > 0.03; al /= 1.15)
      {
        // introducing new alpha
        P.alpha[0] = P.alpha[1] = al;
        LT.set_new_lambda(&LT.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
        std::vector<double> values_mot_b;
        std::vector<double> values_mas_b;
        std::vector<double> values_wei_b;
        double mean_b = 0;
        double rel_mass_b = 0;
        double count_b = 0;
        double weighted_b = 0;
        double cov_w_b = 0;
        double s_b = 0;
        std::vector<double> values_per_b;
        for (unsigned n = 0; t < (c + 1) * tmax; ++n)
        {
          t = LT.run_until(burnin + n * every + c * tmax);

          values_mas_b.push_back(double(LT.max_cluster_size_nr()) / double(P.N));
          values_mot_b.push_back(LT.motility_fraction());
          rel_mass_b += double(LT.max_cluster_size_nr()) / double(P.N);
          mean_b += LT.motility_fraction();
          values_per_b.push_back(LT.perimeter() / double(P.N));
          s_b += LT.perimeter() / double(P.N);
          count_b++;
          // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

          hist_t hist_b = LT.cluster_distribution_particle_number();
          double second_moment_b = 0;
          double first_moment_b = 0;
          for (unsigned i = 0; i < hist_b.size(); i += 2)
          {
            second_moment_b += hist_b[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment_b += hist_b[i] * ((i + 2) / 2);
          }
          // std::cout << second_moment << " " << first_moment << endl;
          values_wei_b.push_back(second_moment_b / first_moment_b * 1.0 / double(P.N));
          weighted_b += second_moment_b / first_moment_b * 1.0 / double(P.N);
        }
        mean_b = mean_b / count_b;
        rel_mass_b = rel_mass_b / count_b;
        weighted_b = weighted_b / count_b;
        s_b = s_b / count_b;

        double cov_mot_b = 0;
        for (auto &val : values_mot_b)
        {
          cov_mot_b += pow(val - mean_b, 2);
        }
        cov_mot_b = cov_mot_b / (values_mot_b.size() - 1);

        double cov_mas_b = 0;
        for (auto &val : values_mas_b)
        {
          cov_mas_b += pow(val - rel_mass_b, 2);
        }
        cov_mas_b = cov_mas_b / (values_mas_b.size() - 1);

        for (auto &val : values_wei_b)
        {
          cov_w_b += pow(val - weighted_b, 2);
        }
        cov_w_b = cov_w_b / (values_wei_b.size() - 1);

        double cov_s_b = 0;
        for (auto &val : values_per_b)
        {
          cov_s_b += pow(val - s_b, 2);
        }
        cov_s_b = cov_s_b / (values_per_b.size() - 1);

        backward << al << " " << mean_b << " " << cov_mot_b << " " << rel_mass_b << " " << cov_mas_b << " " << weighted_b << " " << cov_w_b << " " << s_b << " " << cov_s_b << endl;
        c++;
      }
    }
    else if (output == "motility")
    {
      if (details == 0)
      {
        ofstream outfile;
        string name = "./lars_sim/Data/motility/square_perc_low_L_50";
        string outputname = name + "_" + occ_p + ".txt";
        outfile.open(outputname);
        for (double al = 0.0; al < 0.05; al += 0.001)
        {
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = al;
          Lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          std::vector<double> values_wei;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          double weighted = 0;
          double cov_w = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
            mean += LB.motility_fraction();
            count++;
            // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distribution_particle_number();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i += 2)
            {
              second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
              first_moment += hist[i] * ((i + 2) / 2);
            }
            // std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0 / double(P.N));
            weighted += second_moment / first_moment * 1.0 / double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;

          double cov_mot = 0;
          for (auto &val : values_mot)
          {
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot / (values_mot.size() - 1);

          double cov_mas = 0;
          for (auto &val : values_mas)
          {
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas / (values_mas.size() - 1);

          for (auto &val : values_wei)
          {
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w / (values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
      }
      else if (details == 1)
      {
        ofstream outfile;
        string name = "./lars_sim/Data/motility/square_perc_details";
        string outputname = name + "_" + occ_p + ".txt";
        outfile.open(outputname);
        for (double al = 0.085; al < 0.095; al += 0.0005)
        {
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = al;
          Lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
            mean += LB.motility_fraction();
            count++;
            // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;
          }
          mean = mean / count;
          rel_mass = rel_mass / count;

          double cov_mot = 0;
          for (auto &val : values_mot)
          {
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot / (values_mot.size() - 1);

          double cov_mas = 0;
          for (auto &val : values_mas)
          {
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas / (values_mas.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << endl;
        }
      }
    }
    else if (output == "stable")
    {
      ofstream part, numb;
      part.open("./lars_sim/Data/stable/square.txt");
      numb.open("./lars_sim/Data/stable/square_number.txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = L.run_until(burnin + n * every);
        // only doing a positional output here
        numb << t << " " << L.number_cluster() << endl;
        part << ParticleWriter(L, part) << endl;
      }
    }
    else if (output == "weighted")
    {
      ofstream output, part;
      part.open("./lars_sim/Data/weighted/square_part.txt");
      output.open("./lars_sim/Data/weighted/square.txt");

      for (double n = 1; t < burnin + until; n++)
      {
        t = L.run_until(burnin + n * every);
        part << ParticleWriter(L, part) << endl;
        hist_t hist = L.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i += 2)
        {
          second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
          first_moment += hist[i] * ((i + 2) / 2);
        }
        // std::cout << second_moment << " " << first_moment << endl;
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        output << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
      }
    }
    else if (output == "distribution")
    {
      ofstream outfile, outfile2;
      outfile.open("./lars_sim/Data/dist/square_" + occ_p + "_special.txt");
      outfile2.open("./lars_sim/Data/dist/square_dens_" + occ_p + "_special.txt");
      hist_t dist(5);
      for (double n = 0; t < burnin + until; n++)
      {
        t = L.run_until(burnin + n * every);
        vec_d dens = L.density();
        for (const auto &m : dens)
          outfile2 << m << " ";
        outfile2 << endl;
        hist_t dr = L.particle_neighbour_dist();
        dist[std::slice(0, dr.size(), 1)] += dr;
      }
      for (const auto &m : dist)
        outfile << m << " ";
      outfile << endl;
    }
    else if (output == "area")
    {
      ofstream surf, part, clust, border;
      surf.open("./lars_sim/Data/surf/square_sv" + occ_p + ".txt");
      part.open("./lars_sim/Data/surf/square_part" + occ_p + ".txt");
      // clust.open("./lars_sim/Data/surf/square_clust"+occ_p+".txt");
      border.open("./lars_sim/Data/surf/square_border" + occ_p + ".txt");

      // unsigned so only one cluster each gets printed
      unsigned check_1 = 0;
      unsigned check_2 = 0;
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = L.run_until(burnin + n * every);
        vec a = L.surface_volume();
        for (const auto &k : a)
          surf << k << " ";
        surf << endl;
        part << ParticleWriter(L, part) << endl;
        for (const auto &n : L.cluster_surface())
          border << n << " ";
        border << endl;

        if (L.clust_size(8, 10) == 1 && check_1 == 0)
        {
          vec n = L.single_cluster(8, 10);
          for (const auto &m : n)
            clust << m << " ";
          clust << endl;
          check_1++;
        }

        if (L.clust_size(3700, 4000) == 1 && check_2 == 0)
        {
          vec n = L.single_cluster(3700, 4000);
          for (const auto &m : n)
            clust << m << " ";
          clust << endl;
          check_2++;
        }
      }
    }
    else if (output == "stopping time")
    {
      if (details == 0){
        ofstream outfile;

        outfile.open("./lars_sim/Data/stopping/square_" + occ_p + "_" + alpha_p + ".txt");
        for (unsigned n = 0; t < burnin + until; ++n)
        {
          t = L.run_until(burnin + n * every);
          vec_d st = L.stopping(t);
          for (const auto &m : st)
            outfile << m << " ";
          outfile << endl;
        }
      }
      else if (details == 1){
        ofstream outfile;
        outfile.open("./lars_sim/Data/stopping/square_" + occ_p + "_stoppingtimes.txt");
        vec_d alp;
        for (int i = 0; i < 3; i++){
          for (int j = 1; j < 10; j++){
            alp.push_back(0.001 * j * pow(10, i));
          }
        }
        // looping over all alpha
        for (const auto & a: alp){
          double t = 0;
          P.alpha[0] = P.alpha[1] = a;
          cout << " --------------------------------------------- " << endl; 
          cout << "Alpha is " << P.alpha[0] << endl;
          // new lattice
          Lattice L(P, rng);

          
          int n_bins = 1000;
          double stopping_time = 0.0;
          vec_d Tstop;
          vec_d hist(n_bins);
          // loop over time
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = L.run_until(burnin + n * every);
            vec_d st = L.stopping(t);
            Tstop.insert(end(Tstop), begin(st), end(st));
          }
          // make histogram
          // bin width
          double width = (*max_element(Tstop.begin(), Tstop.end()) - *min_element(Tstop.begin(), Tstop.end())) / double(n_bins + 1);
          cout << "Width of histogram bins " << width << endl;

          for (const auto &m : Tstop){
            for (int k = 0; k < n_bins; k++){
              if (m - (k + 1) * width < 0){
                hist[k]++;
                break;
              }
            }
          }
          int sum_hist = accumulate(hist.begin(), hist.end(), 0);
          // calculate mean
          for (int k = 0; k < hist.size(); k++){
            stopping_time += width * (k + 0.5) * hist[k] / double(sum_hist);
          }
          outfile << P.alpha[0] << " " << stopping_time << endl;
        }
      }
    }
    else if (output == "perimeter")
    {
      ofstream outfile, pars, probs;
      outfile.open("./lars_sim/Data/perimeter/square_" + occ_p + ".txt");
      pars.open("./lars_sim/Data/perimeter/square_pars_" + occ_p + ".txt");
      probs.open("./lars_sim/Data/perimeter/square_probs_" + occ_p + ".txt");
      for (double alp = 1e-3; alp <= 1e2; alp *= 1.45)
      {
        for (double dens = 0.001; dens < 1.0; dens += .02)
        {
          Parameters P_h;
          P_h.N = unsigned(P.L[0] * P.L[0] * P.n_max * dens);
          P_h.alpha[0] = P.alpha[1] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          Lattice LB(P_h, rng);
          t = 0;

          // output for parameters (curious to see the values)
          pars << P_h.N << endl;

          // jamming, max cl size, perimeter, weighted mean cl size
          double j = 0;
          double m = 0;
          double s = 0;
          double w = 0;

          vec_d prob(P.n_max + 1);

          // count for mean vals
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            j += LB.motility_fraction();
            m += LB.max_cluster_size_nr();
            s += LB.perimeter();
            w += LB.w_N();
            hist_t p = LB.occupation_prob();
            for (unsigned i = 0; i < P.n_max + 1; i++)
            {
              prob[i] += double(p[i]);
            }
            count++;
          }

          j = j / count;
          m = m / count;
          s = s / count;
          w = w / count;

          m = m / double(P_h.N);
          w = w / double(P_h.N);

          for (unsigned i = 0; i < P.n_max + 1; i++)
          {
            prob[i] /= double(P.L[0] * P.L[1] * count);
            probs << prob[i] << " ";
          }
          probs << endl;

          outfile << j << " " << m << " " << s << " " << w << " ";
        }
        outfile << endl;
      }
    }
    else if (output == "correlation")
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/corr/square_" + occ_p + "_" + alpha_p + ".txt");
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = L.run_until(burnin + n * every);
        vec output = L.occ_array();
        for (const auto &m : output)
          outfile << m << " ";
        outfile << endl;
      }
    }
    else if (output == "function")
    {
      if (details == 0)
      {
        ofstream outfile, dist;
        outfile.open("./lars_sim/Data/displacement/sq_varkur_rho_25_" + occ_p + "_test.txt");
        dist.open("./lars_sim/Data/displacement/sq_dist_rho_05_" + occ_p + "_test.txt");
        vector<long double> variance;
        vector<long double> kurtosis;
        // number of configurations
        unsigned nr_configuration = 10;
        map<unsigned, vec_d> x;
        vec dist_time = {50, 150, 200, 500, 5000, unsigned(until - 1)};
        for (unsigned k = 0; t < until; ++k)
        {
          // add part to map for distribution.
          for (const auto &m : dist_time)
          {
            if (m == k)
              x.insert(pair<unsigned, vec_d>(k, {}));
          }
          // fill the arrays with 0, probably not the best way to do it.
          t = L.run_until(k * every);
          if (k > 0)
          {
            variance.push_back(0);
            kurtosis.push_back(0);
          }
        }

        for (unsigned i = 0; i < nr_configuration; i++)
        {
          if (i % (nr_configuration / 10) == 0)
          {
            for (int o = 0; o < i / (nr_configuration / 10); o++)
            {
              cout << "#";
            }
            cout << endl;
          }
          Lattice L(P, rng);
          t = 0;
          for (unsigned n = 1; t < until; ++n)
          {
            t = L.run_until(n * every);
            vec_d pos = L.x_coor_3();
            for (unsigned k = 0; k < P.N; k++)
            {
              variance[n - 1] += pow(pos[k], 2);
              kurtosis[n - 1] += pow(pos[k], 4);
            }
            for (unsigned m = 0; m < dist_time.size(); m++)
            {
              if (dist_time[m] == n)
              {

                x[m].insert(x[m].end(), pos.begin(), pos.end());
                // cout << m << " " << x[m].size() << " " << pos.size() << endl;
              }
            }
          }
        }

        // calculate variance and kurtosis
        for (unsigned n = 0; n < variance.size(); n++)
        {
          variance[n] /= double(P.N * nr_configuration);
          kurtosis[n] /= (double(P.N * nr_configuration) * pow(variance[n], 2));
          outfile << (n + 1) << " " << variance[n] << " " << kurtosis[n] << endl;
        }

        // fill distribution file
        for (unsigned n = 0; n < dist_time.size(); n++)
        {
          for (unsigned k = 0; k < nr_configuration * P.N; k++)
          {
            // cout << x[n].size() << endl;
            dist << x[n][k] << " ";
          }
          dist << endl;
        }
      }
      // ############################ TESTING SOME MORE STUFF #################################
      else if (details == 1)
      {
        ofstream outfile;
        outfile.open("./lars_sim/Data/displacement/sq_varkur_rho_25_" + occ_p + "_test_test_test.txt");
        vector<long double> variance;
        vector<long double> kurtosis;
        // number of configurations
        unsigned nr_configuration = 1000;
        // map<unsigned, vec_d> x;
        for (unsigned k = 0; t < until; ++k)
        {
          t = L.run_until(k * every);
          if (k > 0)
          {
            variance.push_back(0);
            kurtosis.push_back(0);
          }
          // x.insert(pair<unsigned, vec_d> (k, {}));
        }

        for (unsigned i = 0; i < nr_configuration; i++)
        {
          if (i % (nr_configuration / 10) == 0)
          {
            for (int o = 0; o < i / (nr_configuration / 10); o++)
            {
              cout << "#";
            }
            cout << endl;
          }
          Lattice L(P, rng);
          t = 0;
          for (unsigned n = 1; t < until; ++n)
          {
            t = L.run_until(n * every);
            vec_d pos = L.x_coor_3();
            for (unsigned k = 0; k < P.N; k++)
            {
              variance[n - 1] += pow(pos[k], 2);
              kurtosis[n - 1] += pow(pos[k], 4);
            }
          }
        }

        for (unsigned n = 0; n < variance.size(); n++)
        {
          variance[n] /= double(P.N * nr_configuration);
          kurtosis[n] /= (double(P.N * nr_configuration) * pow(variance[n], 2));
          outfile << (n + 1) << " " << variance[n] << " " << kurtosis[n] << endl;
        }
      }
    }

    else if (output == "tagged")
    {
      ofstream outfile, dist;
      outfile.open("./lars_sim/Data/displacement/detailed_displacement/sq_varkur_rho_60_" + occ_p + "_" + tumb + ".txt");
      dist.open("./lars_sim/Data/displacement/detailed_displacement/sq_dist_rho_60_" + occ_p + "_" + tumb + ".txt");
      vector<long double> variance;
      vector<long double> kurtosis;
      // number of configurations
      unsigned nr_configuration = 1;
      map<unsigned, vec_d> x;
      vec dist_time = {50, 150, 200, 500, 5000, unsigned(until - 1)};
      for (unsigned k = 0; k * every < until; ++k)
      {
        // add part to map for distribution.
        for (const auto &m : dist_time)
        {
          if (m == k)
            x.insert(pair<unsigned, vec_d>(k, {}));
        }
        // fill the arrays with 0, probably not the best way to do it.
        if (k > 0)
        {
          variance.push_back(0);
          kurtosis.push_back(0);
        }
      }
      // need to iterate over several lattice to achieve good statistics
      for (unsigned i = 0; i < nr_configuration; i++)
      {
        // just a little progress bar so I am not kept in the dark for too long
        if (nr_configuration >= 10)
        {
          if (i % (nr_configuration / 10) == 0)
          {
            for (int o = 0; o < i / (nr_configuration / 10); o++)
            {
              cout << "#";
            }
            cout << endl;
          }
        }
        else
        {
          for (int o = 0; o <= i; o++)
          {
            cout << "#";
          }
          cout << endl;
        }
        Lattice L(P, rng);
        t = 0;
        for (unsigned n = 1; t < until; ++n)
        {
          t = L.run_until(n * every);
          vec_d pos = L.x_coor_3();
          for (unsigned k = 0; k < P.N; k++)
          {
            variance[n - 1] += pow(pos[k], 2);
            kurtosis[n - 1] += pow(pos[k], 4);
          }
          for (unsigned m = 0; m < dist_time.size(); m++)
          {
            if (dist_time[m] == n)
            {

              x[m].insert(x[m].end(), pos.begin(), pos.end());
              // cout << m << " " << x[m].size() << " " << pos.size() << endl;
            }
          }
        }
      }

      // calculate variance and kurtosis
      for (unsigned n = 0; n < variance.size(); n++)
      {
        variance[n] /= double(P.N * nr_configuration);
        kurtosis[n] /= (double(P.N * nr_configuration) * pow(variance[n], 2));
        outfile << (n + 1) << " " << variance[n] << " " << kurtosis[n] << endl;
      }

      // fill distribution file
      for (unsigned n = 0; n < dist_time.size(); n++)
      {
        for (unsigned k = 0; k < nr_configuration * P.N; k++)
        {
          // cout << x[n].size() << endl;
          dist << x[n][k] << " ";
        }
        dist << endl;
      }
    }
    else if (output == "particles")
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/for_latex/square_" + occ_p + ".txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = L.run_until(burnin + n * every);
        // only doing a positional output here

        outfile << ParticleWriter(L, outfile) << endl;
      }
    }
    else if (output == "detailed balance"){
      ofstream probs;
      probs.open("./lars_sim/Data/perimeter/square_probs_" + occ_p + ".txt");
        for (double dens = 0.001; dens < 1.0; dens += .02)
        {
          Parameters P_h;
          P_h.N = unsigned(P.L[0] * P.L[0] * P.n_max * dens);
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          Lattice LB(P_h, rng);
          t = 0;

          vec_d prob(P.n_max + 1);

          // count for mean vals
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            hist_t p = LB.occupation_prob();
            for (unsigned i = 0; i < P.n_max + 1; i++)
            {
              prob[i] += double(p[i]);
            }
            count++;
          }

          probs << dens << " ";
          for (unsigned i = 0; i < P.n_max + 1; i++)
          {
            prob[i] /= double(P.L[0] * P.L[1] * count);
            probs << prob[i] << " ";
          }
          probs << endl;
        }
    }
    else
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/for_latex/square_" + occ_p + ".txt");
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = L.run_until(burnin + n * every);
        if (output == "particles")
          outfile << ParticleWriter(L, outfile) << std::endl;
        else if (output == "vacancies")
          outfile << VacancyWriter(L) << std::endl;
        else
        {
          // Ensure that all the particle directions are at the simulation time
          L.realise_directions();
          outfile << SnapshotWriter(L) << std::endl;
        }
      }
    }
  }
  else if (lattice_type == "triangular")
  {
    // Initialise a random number generator and set up the model
    std::mt19937 rng((std::random_device())());
    Triangle_lattice TL(P, rng);

    // Snapshot sequence has a header that sets out the simulation parameters
    std::cout << "# L = [ ";
    for (const auto &L : P.L)
      std::cout << L << " ";
    std::cout << "]" << std::endl;
    std::cout << "# N = " << P.N << std::endl;
    std::cout << "# alpha = [ ";
    for (const auto &alpha : P.alpha)
      std::cout << alpha << " ";
    std::cout << "]" << std::endl;
    std::cout << "# output = " << output << std::endl;
    std::cout << "# initial = " << burnin << std::endl;
    std::cout << "# interval = " << every << std::endl;
    if (localaverage > 0)
    {
      std::cout << "# localaverage = " << localaverage << std::endl;
      std::cout << "# localinterval = " << localinterval << std::endl;
    }
    std::cout << "# occupation number = " << P.n_max << std::endl;

    double t = 0;

    if (output == "clusters")
    {
      if (localaverage == 0)
      {
        // We sum the histograms over all measurements
        hist_t sumhist;
        hist_t sumhist_nr;
        for (unsigned n = 0; t < burnin + until; ++n)
        {
          t = TL.run_until(burnin + n * every);
          hist_t hist = TL.cluster_distributions();
          hist_t hist_nr = TL.cluster_distributions_particle_numbers();
          // Only area
          if (hist.size() > sumhist.size())
          {
            hist[std::slice(0, sumhist.size(), 1)] += sumhist;
            sumhist = std::move(hist);
          }
          else
          {
            sumhist[std::slice(0, hist.size(), 1)] += hist;
          }

          // with occ number as well
          if (hist_nr.size() > sumhist_nr.size())
          {
            hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
            sumhist_nr = std::move(hist_nr);
          }
          else
          {
            sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
          }
        }
        // output for each of the distributions
        ofstream outfile, outfile_nr;
        outfile.open(txtoutput);
        for (const auto &k : sumhist)
          outfile << k << " ";
        outfile << endl;
        outfile_nr.open(txtoutput_nr);
        for (const auto &k : sumhist_nr)
          outfile_nr << k << " ";
        outfile_nr << endl;
      }
      else
      {
        // We perform local averages of the histograms at the given measurement interval
        for (unsigned n = 0; t < burnin + until; ++n)
        {
          t = TL.run_until(burnin + n * every);
          hist_t sumhist = TL.cluster_distributions();
          hist_t sumhist_nr = TL.cluster_distributions_particle_numbers();
          for (unsigned m = 1; m < localaverage; ++m)
          {
            t = TL.run_until(burnin + n * every + m * localinterval);
            hist_t hist = TL.cluster_distributions();
            hist_t hist_nr = TL.cluster_distributions_particle_numbers();
            // Add hist to sumhist taking into account they may have different lengths
            if (hist.size() > sumhist.size())
            {
              hist[std::slice(0, sumhist.size(), 1)] += sumhist;
              sumhist = std::move(hist);
            }
            else
            {
              sumhist[std::slice(0, hist.size(), 1)] += hist;
            }
            if (hist_nr.size() > sumhist_nr.size())
            {
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            }
            else
            {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }
          }
          ofstream outfile, outfile_nr;
          outfile.open(txtoutput);
          for (const auto &k : sumhist)
            outfile << k << " ";
          outfile << endl;
          outfile_nr.open(txtoutput_nr);
          for (const auto &k : sumhist_nr)
            outfile_nr << k << " ";
          outfile_nr << endl;
        }
      }
    }
    else if (output == "heatmap")
    {
      ofstream outfile, outfile_avg, outfile_nr, outfile_avg_nr;
      outfile.open("./lars_sim/heatmap/tri_alpha_N_n_2.txt");
      outfile_avg.open("./lars_sim/heatmap/tri_alpha_N_n_2_avg.txt");
      outfile_nr.open("./lars_sim/heatmap/tri_nr_alpha_N_n_2.txt");
      outfile_avg_nr.open("./lars_sim/heatmap/tri_nr_alpha_N_n_2_avg.txt");
      for (double alp = 1e-6; alp <= 1.0; alp *= 1.3)
      {
        for (unsigned N = 100 * P.n_max; N <= P.L[0] * P.L[0] * P.n_max; N += 200 * P.n_max)
        {
          Parameters P_h;
          P_h.N = N;
          P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          std::size_t maxsize = 1;
          std::size_t maxsize_nr = 1;
          Triangle_lattice LB(P_h, rng);
          t = 0;

          // also taking mean over all used timesteps
          double mean = 0.0;
          double mean_nr = 0.0;
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            maxsize = std::max(maxsize, LB.max_cluster_size());
            maxsize_nr = std::max(maxsize_nr, LB.max_cluster_size_nr());
            mean += LB.avg_cluster_size();
            mean_nr += LB.avg_cluster_size_nr();
            count += 1;
          }

          mean = mean / count;
          mean_nr = mean_nr / count;

          outfile << maxsize << " ";
          outfile_avg << mean << " ";
          outfile_nr << maxsize_nr << " ";
          outfile_avg_nr << mean_nr << " ";
        }
        outfile << endl;
        outfile_avg << endl;
        outfile_nr << endl;
        outfile_avg_nr << endl;
      }
    }
    else if (output == "function")
    {
      ofstream outfile;
      outfile.open("./lars_sim/ActivePerc/CStransition/NewScaling_"+std::to_string(region) +"_"+ lattice_type + "_" + dens + "_" + size + "_" + occ_p + txt);
      // outfile.open("./NewScaling_"+ lattice_type + "_" + dens + "_" + size + "_" + occ_p + txt);

      // pars.open("./lars_sim/Data/perimeter/tri_pars_"+occ_p+"test.txt");
      for (double alp = 0.01 * pow(1.16, 4.0*region); alp <= 1 * pow(1.16, 4.0*(region+1)); alp *= 1.16)
      {
        // do the evaluation in python for mean and variance (errorbars) for now
        for (int n_iteration = 0; n_iteration <= 3; n_iteration++){
          Parameters P_h;
          P_h.N = P.N;
          P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          Triangle_lattice LB(P_h, rng);
          t = 0;

          std::cout << P_h.N << n_iteration << endl;

          // jamming, max cl size, perimeter, weighted mean cl size
          double j = 0;
          double m = 0;
          double s = 0;
          double w = 0;

          // count for mean vals
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            j += LB.motility_fraction();
            m += LB.max_cluster_size_nr();
            s += LB.perimeter();
            w += LB.w_N();
            count++;
          }

          j = j / count;
          m = m / count;
          s = s / count;
          w = w / count;

          m = m / double(P_h.N);
          // s = s;
          w = w / double(P_h.N);

          outfile << alp << " " << j << " " << m << " " << s << " " << w << " ";
          outfile << endl;
        }
      }
    }
    else if (output == "surface scaling")
    {
      ofstream outfile, Distoutfile;
      outfile.open("./lars_sim/ActivePerc/CStransition/SurfaceScaling_"+ lattice_type + "_" + dens + "_" + size + "_" + occ_p + txt); 
      Distoutfile.open("./lars_sim/ActivePerc/CStransition/SurfaceDist_"+ lattice_type + "_" + dens + "_" + size + "_" + occ_p + txt);

      for (double rho = 1e-3; rho <= 1; rho *= 1.5)
      {
        int N = rho * pow(P.L[0], 2) * P.n_max;

        double t = 0;
        
        // We sum the histograms over all measurements
        // vec dist;
        // vec dist2;
        // for (int length = 0; length <= 2 * N; length++){
        //   dist.push_back(0);
        //   dist2.push_back(0);
        // }
        hist_t part;
        // do the evaluation in python for mean and variance (errorbars) for now
        for (int n_iteration = 0; n_iteration <= 3; n_iteration++){
          Parameters P_h;
          P_h.N = N;
          P_h.alpha[0] = P.alpha[1] = P.alpha[2] = P.alpha[0];
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          Triangle_lattice LB(P_h, rng);
          t = 0;

          std::cout << P_h.N << " " << n_iteration << endl;

          // jamming, max cl size, perimeter, weighted mean cl size
          double j = 0;
          double m = 0;
          double s = 0;
          double w = 0;

          // count for mean vals
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
              // vec hist_nr = TL.surface_volume_nr();
              hist_t part_nr = LB.cluster_distributions();
              // for (unsigned k = 0; k < hist_nr.size(); k += 2)
              //   {
              //       dist[hist_nr[k]]++;
              //   }
              // cout << "start";
              // for (const auto &k: part_nr){
              //     cout << k << " ";
              //   }
              if (part_nr.size() > part.size())
              {
                part_nr[std::slice(0, part.size(), 1)] += part;
                part = std::move(part_nr);
              }
              else
              {
                part[std::slice(0, part_nr.size(), 1)] += part_nr;
              }
            
            j += LB.motility_fraction();
            m += LB.max_cluster_size_nr();
            s += LB.perimeter();
            w += LB.w_N();
            count++;
          }

          j = j / count;
          m = m / count;
          s = s / count;
          w = w / count;

          m = m / double(P_h.N);
          // s = s;
          w = w / double(P_h.N);

          outfile << rho << " " << j << " " << m << " " << s << " " << w << " " << N << " ";
          outfile << endl;
        
        }
        // for (const auto &k : dist)
        //   Distoutfile << k << " ";
        // Distoutfile << endl;
        for (const auto &k: part){
          Distoutfile << k << " ";
          // if (k != 0) cout << k << " ";
        }
        Distoutfile << endl;
      }
    }
    else if (output == "stopping time")
    {
      ofstream outfile;

      outfile.open("./lars_sim/Data/stopping/tri_" + occ_p + "_" + alpha_p + ".txt");
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = TL.run_until(burnin + n * every);
        vec_d st = TL.stopping(t);
        for (const auto &m : st)
          outfile << m << " ";
        outfile << endl;
      }
    }
    else if (output == "number")
    {
      ofstream outfile;
      string name = "./lars_sim/Data/testing/tri_";
      string name_number = name + "number" + "_" + tumb + "_" + dens + "_" + size + "_" + occ_p + txt;
      outfile.open(name_number);

      vec_d res, clnr, cn, j, m;
      for (unsigned n = 0; burnin + pow(1.01, n) * every < burnin + until; n++)
      {
        res.push_back(0);
        clnr.push_back(0);
        cn.push_back(0);
        j.push_back(0);
        m.push_back(0);
      }

      unsigned iterations = 10;

      for (int q = 0; q <= iterations; q++)
      {
        Triangle_lattice TL(P, rng);
        t = TL.run_until(burnin);
        hist_t hist = TL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i += 2)
        {
          second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
          first_moment += hist[i] * ((i + 2) / 2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        res[0] += weighted;
        clnr[0] += TL.avg_cluster_size_nr();
        cn[0] += TL.perimeter();
        j[0] += TL.motility_fraction();
        m[0] += TL.max_cluster_size_nr();
        for (unsigned n = 1; t < burnin + until; n++)
        {
          t = TL.run_until(burnin + pow(1.01, n) * every);
          hist_t hist = TL.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i += 2)
          {
            second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment += hist[i] * ((i + 2) / 2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          res[n] += weighted;
          clnr[n] += TL.avg_cluster_size_nr();
          cn[n] += TL.perimeter();
          j[n] += TL.motility_fraction();
          m[n] += TL.max_cluster_size_nr();
        }
        cout << "We are at " << q << endl;
      }
      for (unsigned i = 0; i <= res.size(); i++)
      {
        outfile << burnin + pow(1.01, i) * every << " " << res[i] / double(iterations + 1) << " " << cn[i] / double(iterations + 1) << " " << j[i] / double(iterations + 1) << " " << m[i] / double(iterations + 1) << " " << double(clnr[i]) / double(iterations + 1) << endl;
      }
    }
    else if (output == "weighted")
    {
      ofstream output, part;
      part.open("./lars_sim/Data/weighted/tri_part.txt");
      output.open("./lars_sim/Data/weighted/tri.txt");

      for (double n = 1; t < burnin + until; n++)
      {
        t = TL.run_until(burnin + n * every);
        part << TriangleParticleWriter(TL, part) << endl;
        hist_t hist = TL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i += 2)
        {
          second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
          first_moment += hist[i] * ((i + 2) / 2);
        }
        // std::cout << second_moment << " " << first_moment << endl;
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        output << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
      }
    }
    else if (output == "stable")
    {
      ofstream part, numb;
      part.open("./lars_sim/Data/stable/triangular.txt");
      numb.open("./lars_sim/Data/stable/triangular_number.txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = TL.run_until(burnin + n * every);
        // only doing a positional output here
        numb << t << " " << TL.number_cluster() << endl;
        part << TriangleParticleWriter(TL, part) << endl;
      }
    }
    else if (output == "motility")
    {
      if (details == 0)
      {
        ofstream outfile;
        string name = "./lars_sim/Data/motility/triangular_perc";
        string outputname = name + "_" + occ_p + ".txt";
        outfile.open(outputname);
        for (double phi = 0.01; phi < 1 - pow(1 - 0.5, 1 / double(P.n_max)); phi += 0.03)
        {
          P.N = int(phi * P.L[0] * P.L[1] * P.n_max);
          for (double al = 0.005; al < 0.2; al += 0.005)
          {
            // defining lattice for new alpha
            P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
            Triangle_lattice LB(P, rng);
            t = 0;
            std::vector<double> values_mot;
            std::vector<double> values_mas;
            std::vector<double> values_wei;
            double mean = 0;
            double rel_mass = 0;
            double count = 0;
            double weighted = 0;
            double cov_w = 0;
            for (unsigned n = 0; t < burnin + until; ++n)
            {
              t = LB.run_until(burnin + n * every);
              values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
              values_mot.push_back(LB.motility_fraction());
              rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
              mean += LB.motility_fraction();
              count++;
              // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

              hist_t hist = LB.cluster_distributions_particle_numbers();
              double second_moment = 0;
              double first_moment = 0;
              for (unsigned i = 0; i < hist.size(); i += 2)
              {
                second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
                first_moment += hist[i] * ((i + 2) / 2);
              }
              // std::cout << second_moment << " " << first_moment << endl;
              values_wei.push_back(second_moment / first_moment * 1.0 / double(P.N));
              weighted += second_moment / first_moment * 1.0 / double(P.N);
            }
            mean = mean / count;
            rel_mass = rel_mass / count;
            weighted = weighted / count;

            double cov_mot = 0;
            for (auto &val : values_mot)
            {
              cov_mot += pow(val - mean, 2);
            }
            cov_mot = cov_mot / (values_mot.size() - 1);

            double cov_mas = 0;
            for (auto &val : values_mas)
            {
              cov_mas += pow(val - rel_mass, 2);
            }
            cov_mas = cov_mas / (values_mas.size() - 1);

            for (auto &val : values_wei)
            {
              cov_w += pow(val - weighted, 2);
            }
            cov_w = cov_w / (values_wei.size() - 1);

            outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << " " << phi << endl;
          }
        }
      }
      else if (details == 1)
      {
        ofstream outfile;
        string name = "./lars_sim/Data/motility/triangular_perc_details";
        string outputname = name + "_" + occ_p + ".txt";
        outfile.open(outputname);
        for (double al = 0.06; al < 0.07; al += 0.0005)
        {
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
          Triangle_lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
            mean += LB.motility_fraction();
            count++;
            // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;
          }
          mean = mean / count;
          rel_mass = rel_mass / count;

          double cov_mot = 0;
          for (auto &val : values_mot)
          {
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot / (values_mot.size() - 1);

          double cov_mas = 0;
          for (auto &val : values_mas)
          {
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas / (values_mas.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << endl;
        }
      }
    }
    else if (output == "area")
    {
      ofstream surf, part, clust, border;
      surf.open("./lars_sim/Data/surf/tri_sv" + occ_p + ".txt");
      part.open("./lars_sim/Data/surf/tri_part" + occ_p + ".txt");
      clust.open("./lars_sim/Data/surf/tri_clust" + occ_p + ".txt");
      border.open("./lars_sim/Data/surf/tri_border" + occ_p + ".txt");
      // unsigned so only one cluster each gets printed
      // unsigned check_1 = 0;
      // unsigned check_2 = 0;
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = TL.run_until(burnin + n * every);
        vec a = TL.surface_volume();
        for (const auto &k : a)
          surf << k << " ";
        surf << endl;
        part << TriangleParticleWriter(TL, part) << endl;
        for (const auto &n : TL.cluster_surface())
          border << n << " ";
        border << endl;
      }
    }
    else if (output == "lagging")
    {
      ofstream outfile, backward;
      string outputname = "./lars_sim/ActivePerc/CStransition/HystF"+ lattice_type + "_" + dens + "_" + size + "_" + occ_p + txt;
      outfile.open(outputname);
      outputname = "./lars_sim/ActivePerc/CStransition/HystB"+ lattice_type + "_" + dens + "_" + size + "_" + occ_p + txt;
      backward.open(outputname);
      // foward hysteresis, i.e. start below critical point and move up
      Triangle_lattice LB(P, rng);
      double tmax = burnin + until;
      unsigned c = 0;
      for (double al = 0.02; al <= 0.03; al += 0.0005)
      {
        cout << al << endl;
        // introducing new alpha
        P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
        LB.set_new_lambda(&LB.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
        std::vector<double> values_mot;
        std::vector<double> values_mas;
        std::vector<double> values_wei;
        std::vector<double> values_per;
        double mean = 0;
        double rel_mass = 0;
        double count = 0;
        double weighted = 0;
        double cov_w = 0;
        double s = 0;
        for (unsigned n = 0; t < (c + 1) * tmax; ++n)
        {
          t = LB.run_until(burnin + n * every + c * tmax);

          values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
          values_mot.push_back(LB.motility_fraction());
          rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
          mean += LB.motility_fraction();
          values_per.push_back(LB.perimeter() / double(P.N));
          s += LB.perimeter();
          count++;
          // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

          hist_t hist = LB.cluster_distributions_particle_numbers();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i += 2)
          {
            second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment += hist[i] * ((i + 2) / 2);
          }
          // std::cout << second_moment << " " << first_moment << endl;
          values_wei.push_back(second_moment / first_moment * 1.0 / double(P.N));
          weighted += second_moment / first_moment * 1.0 / double(P.N);
        }
        mean = mean / count;
        rel_mass = rel_mass / count;
        weighted = weighted / count;
        s = s / count;

        double cov_mot = 0;
        for (auto &val : values_mot)
        {
          cov_mot += pow(val - mean, 2);
        }
        cov_mot = cov_mot / (values_mot.size() - 1);

        double cov_mas = 0;
        for (auto &val : values_mas)
        {
          cov_mas += pow(val - rel_mass, 2);
        }
        cov_mas = cov_mas / (values_mas.size() - 1);

        for (auto &val : values_wei)
        {
          cov_w += pow(val - weighted, 2);
        }
        cov_w = cov_w / (values_wei.size() - 1);

        double cov_s = 0;
        for (auto &val : values_per)
        {
          cov_s += pow(val - s, 2);
        }
        cov_s = cov_s / (values_per.size() - 1);

        outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << " " << s << " " << cov_s << endl;
        c++;
      }

      Triangle_lattice LT(P, rng);
      c = 0;
      t = 0;
      for (double al = 0.03; al >= 0.02; al -= 0.0005)
      {
        cout << al << endl;
        // introducing new alpha
        P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
        LT.set_new_lambda(&LT.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
        std::vector<double> values_mot_b;
        std::vector<double> values_mas_b;
        std::vector<double> values_wei_b;
        double mean_b = 0;
        double rel_mass_b = 0;
        double count_b = 0;
        double weighted_b = 0;
        double cov_w_b = 0;
        double s_b = 0;
        std::vector<double> values_per_b;
        for (unsigned n = 0; t < (c + 1) * tmax; ++n)
        {
          t = LT.run_until(burnin + n * every + c * tmax);

          values_mas_b.push_back(double(LT.max_cluster_size_nr()) / double(P.N));
          values_mot_b.push_back(LT.motility_fraction());
          rel_mass_b += double(LT.max_cluster_size_nr()) / double(P.N);
          mean_b += LT.motility_fraction();
          values_per_b.push_back(LT.perimeter() / double(P.N));
          s_b += LT.perimeter();
          count_b++;
          // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

          hist_t hist_b = LT.cluster_distributions_particle_numbers();
          double second_moment_b = 0;
          double first_moment_b = 0;
          for (unsigned i = 0; i < hist_b.size(); i += 2)
          {
            second_moment_b += hist_b[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment_b += hist_b[i] * ((i + 2) / 2);
          }
          // std::cout << second_moment << " " << first_moment << endl;
          values_wei_b.push_back(second_moment_b / first_moment_b * 1.0 / double(P.N));
          weighted_b += second_moment_b / first_moment_b * 1.0 / double(P.N);
        }
        mean_b = mean_b / count_b;
        rel_mass_b = rel_mass_b / count_b;
        weighted_b = weighted_b / count_b;
        s_b = s_b / count_b;

        double cov_mot_b = 0;
        for (auto &val : values_mot_b)
        {
          cov_mot_b += pow(val - mean_b, 2);
        }
        cov_mot_b = cov_mot_b / (values_mot_b.size() - 1);

        double cov_mas_b = 0;
        for (auto &val : values_mas_b)
        {
          cov_mas_b += pow(val - rel_mass_b, 2);
        }
        cov_mas_b = cov_mas_b / (values_mas_b.size() - 1);

        for (auto &val : values_wei_b)
        {
          cov_w_b += pow(val - weighted_b, 2);
        }
        cov_w_b = cov_w_b / (values_wei_b.size() - 1);

        double cov_s_b = 0;
        for (auto &val : values_per_b)
        {
          cov_s_b += pow(val - s_b, 2);
        }
        cov_s_b = cov_s_b / (values_per_b.size() - 1);

        backward << al << " " << mean_b << " " << cov_mot_b << " " << rel_mass_b << " " << cov_mas_b << " " << weighted_b << " " << cov_w_b << " " << s_b << " " << cov_s_b << endl;
        c++;
      }
    }
    else if (output == "distribution")
    {
      ofstream outfile, outfile2;
      outfile.open("./lars_sim/Data/dist/tri_" + occ_p + ".txt");
      outfile2.open("./lars_sim/Data/dist/tri_dens_" + occ_p + ".txt");
      hist_t dist(7);
      for (double n = 0; t < burnin + until; n++)
      {

        t = TL.run_until(burnin + n * every);
        vec_d dens = TL.density();
        for (const auto &y : dens)
          outfile2 << y << " ";
        outfile2 << endl;

        hist_t dr = TL.particle_neighbour_dist();
        dist[std::slice(0, dr.size(), 1)] += dr;
      }

      for (const auto &m : dist)
        outfile << m << " ";
      outfile << endl;
    }
    else if (output == "perimeter")
    {
      ofstream outfile, pars;
      outfile.open("./lars_sim/Data/perimeter/tri_" + occ_p + "_t.txt");
      pars.open("./lars_sim/Data/perimeter/tri_pars_" + occ_p + "_t.txt");
      unsigned check = 0;
      for (double alp = 1e-3; alp <= 100.0; alp *= 1.78)
      {
        for (double dens = 0.001; dens < 1.0; dens += .05)
        {
          Parameters P_h;
          P_h.N = unsigned(P.L[0] * P.L[0] * P.n_max * dens);
          P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          Triangle_lattice LB(P_h, rng);
          t = 0;

          // output for parameters (curious to see the values)
          if (check < 7)
          {
            pars << P_h.N << endl;
            check++;
          }

          // jamming, max cl size, perimeter, weighted mean cl size
          double j = 0;
          double m = 0;
          double s = 0;
          double w = 0;

          // count for mean vals
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            j += LB.motility_fraction();
            m += LB.max_cluster_size_nr();
            s += LB.perimeter();
            w += LB.w_N();
            count++;
          }

          j = j / count;
          m = m / count;
          s = s / count;
          w = w / count;

          m = m / double(P_h.N);
          s = s / double(P_h.N);
          w = w / double(P_h.N);

          outfile << j << " " << m << " " << s << " " << w << " ";
        }
        outfile << endl;
      }
    }
    else if (output == "correlation")
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/corr/tri_" + occ_p + "_" + alpha_p + ".txt");
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = TL.run_until(burnin + n * every);
        vec output = TL.occ_array();
        for (const auto &m : output)
          outfile << m << " ";
        outfile << endl;
      }
    }
    else if (output == "particles")
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/for_latex/triangle_" + occ_p + ".txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = TL.run_until(burnin + n * every);
        // only doing a positional output here

        outfile << TriangleParticleWriter(TL, outfile) << endl;
      }
    }
    else
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/for_latex/triangle_" + occ_p + ".txt");
      // outfile << "# L = [ ";
      // for(const auto& L: P.L) outfile << L << " ";
      // outfile << "]" << endl;
      // outfile << "# N = " << P.N << endl;
      // outfile << "# alpha = [ ";
      // for(const auto& alpha: P.alpha) outfile << alpha << " ";
      // outfile << "]" << endl;
      // outfile << "# output = " << output << endl;
      // outfile << "# initial = " << burnin << endl;
      // outfile << "# interval = " << every << endl;

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = TL.run_until(burnin + n * every);
        // only doing a positional output here
        outfile << TriangleParticleWriter(TL, outfile) << endl;
      }
    }
  }
  else if (lattice_type == "hexagonal")
  {

    // Initialise a random number generator and set up the model
    std::mt19937 rng((std::random_device())());
    Hexagonal_lattice HL(P, rng);

    // Snapshot sequence has a header that sets out the simulation parameters
    std::cout << "# L = [ ";
    for (const auto &L : P.L)
      std::cout << L << " ";
    std::cout << "]" << std::endl;
    std::cout << "# N = " << P.N << std::endl;
    std::cout << "# alpha = [ ";
    for (const auto &alpha : P.alpha)
      std::cout << alpha << " ";
    std::cout << "]" << std::endl;
    std::cout << "# output = " << output << std::endl;
    std::cout << "# initial = " << burnin << std::endl;
    std::cout << "# interval = " << every << std::endl;
    if (localaverage > 0)
    {
      std::cout << "# localaverage = " << localaverage << std::endl;
      std::cout << "# localinterval = " << localinterval << std::endl;
    }
    std::cout << "# occupation number = " << P.n_max << std::endl;

    double t = 0;

    if (output == "clusters")
    {
      if (localaverage == 0)
      {
        hist_t sumhist;
        hist_t sumhist_nr;
        // We sum the histograms over all measurements
        ofstream outfile_part;
        outfile_part.open("./lars_sim/gif/hexagonal.txt");

        for (unsigned n = 0; t < burnin + until; ++n)
        {
          t = HL.run_until(burnin + n * every);
          hist_t hist = HL.cluster_distributions();
          hist_t hist_nr = HL.cluster_distribution_particle_number();
          // particle output for comparison
          outfile_part << HexagonalParticleWriter(HL, outfile_part) << endl;
          // Only area
          if (hist.size() > sumhist.size())
          {
            hist[std::slice(0, sumhist.size(), 1)] += sumhist;
            sumhist = std::move(hist);
          }
          else
          {
            sumhist[std::slice(0, hist.size(), 1)] += hist;
          }

          // with occ number as well

          if (hist_nr.size() > sumhist_nr.size())
          {
            hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
            sumhist_nr = std::move(hist_nr);
          }
          else
          {
            sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
          }
        }
        ofstream outfile, outfile_nr;
        outfile.open(txtoutput);
        for (const auto &k : sumhist)
          outfile << k << " ";
        outfile << endl;
        outfile_nr.open(txtoutput_nr);
        for (const auto &k : sumhist_nr)
          outfile_nr << k << " ";
        outfile_nr << endl;
      }
      else
      {
        // We perform local averages of the histograms at the given measurement interval
        for (unsigned n = 0; t < burnin + until; ++n)
        {
          t = HL.run_until(burnin + n * every);
          hist_t sumhist = HL.cluster_distributions();
          hist_t sumhist_nr = HL.cluster_distribution_particle_number();
          for (unsigned m = 1; m < localaverage; ++m)
          {
            t = HL.run_until(burnin + n * every + m * localinterval);
            hist_t hist = HL.cluster_distributions();
            hist_t hist_nr = HL.cluster_distribution_particle_number();
            // Add hist to sumhist taking into account they may have different lengths
            if (hist.size() > sumhist.size())
            {
              hist[std::slice(0, sumhist.size(), 1)] += sumhist;
              sumhist = std::move(hist);
            }
            else
            {
              sumhist[std::slice(0, hist.size(), 1)] += hist;
            }
            if (hist_nr.size() > sumhist_nr.size())
            {
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            }
            else
            {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }
          }
          ofstream outfile, outfile_nr;
          outfile.open(txtoutput);
          for (const auto &k : sumhist)
            outfile << k << " ";
          outfile << endl;
          outfile_nr.open(txtoutput_nr);
          for (const auto &k : sumhist_nr)
            outfile_nr << k << " ";
          outfile_nr << endl;
        }
      }
    }
    else if (output == "particles")
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/for_latex/hexagonal_" + occ_p + ".txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        // only doing a positional output here

        outfile << HexagonalParticleWriter(HL, outfile) << endl;
      }
    }
    else if (output == "lagging")
    {
      ofstream outfile, backward;
      string name = "./lars_sim/Data/motility/hexagonal_perc_fhyst";
      string outputname = name + "_" + occ_p + ".txt";
      outfile.open(outputname);
      name = "./lars_sim/Data/motility/hexagonal_perc_bhyst";
      outputname = name + "_" + occ_p + ".txt";
      backward.open(outputname);
      // foward hysteresis, i.e. start below critical point and move up
      Hexagonal_lattice LB(P, rng);
      double tmax = burnin + until;
      unsigned c = 0;
      for (double al = 0.10; al < 0.126; al += 0.002)
      {
        // introducing new alpha
        P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
        LB.set_new_lambda(&LB.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
        std::vector<double> values_mot;
        std::vector<double> values_mas;
        std::vector<double> values_wei;
        double mean = 0;
        double rel_mass = 0;
        double count = 0;
        double weighted = 0;
        double cov_w = 0;
        for (unsigned n = 0; t < (c + 1) * tmax; ++n)
        {
          t = LB.run_until(burnin + n * every + c * tmax);

          values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
          values_mot.push_back(LB.motility_fraction());
          rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
          mean += LB.motility_fraction();
          count++;
          // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

          hist_t hist = LB.cluster_distribution_particle_number();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i += 2)
          {
            second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment += hist[i] * ((i + 2) / 2);
          }
          // std::cout << second_moment << " " << first_moment << endl;
          values_wei.push_back(second_moment / first_moment * 1.0 / double(P.N));
          weighted += second_moment / first_moment * 1.0 / double(P.N);
        }
        mean = mean / count;
        rel_mass = rel_mass / count;
        weighted = weighted / count;

        double cov_mot = 0;
        for (auto &val : values_mot)
        {
          cov_mot += pow(val - mean, 2);
        }
        cov_mot = cov_mot / (values_mot.size() - 1);

        double cov_mas = 0;
        for (auto &val : values_mas)
        {
          cov_mas += pow(val - rel_mass, 2);
        }
        cov_mas = cov_mas / (values_mas.size() - 1);

        for (auto &val : values_wei)
        {
          cov_w += pow(val - weighted, 2);
        }
        cov_w = cov_w / (values_wei.size() - 1);

        outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        c++;
      }

      Hexagonal_lattice LT(P, rng);
      c = 0;
      t = 0;
      for (double al = 0.126; al > 0.10; al -= 0.002)
      {
        // introducing new alpha
        P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
        LT.set_new_lambda(&LT.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
        std::vector<double> values_mot_b;
        std::vector<double> values_mas_b;
        std::vector<double> values_wei_b;
        double mean_b = 0;
        double rel_mass_b = 0;
        double count_b = 0;
        double weighted_b = 0;
        double cov_w_b = 0;
        for (unsigned n = 0; t < (c + 1) * tmax; ++n)
        {
          t = LT.run_until(burnin + n * every + c * tmax);

          values_mas_b.push_back(double(LT.max_cluster_size_nr()) / double(P.N));
          values_mot_b.push_back(LT.motility_fraction());
          rel_mass_b += double(LT.max_cluster_size_nr()) / double(P.N);
          mean_b += LT.motility_fraction();
          count_b++;
          // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

          hist_t hist_b = LT.cluster_distribution_particle_number();
          double second_moment_b = 0;
          double first_moment_b = 0;
          for (unsigned i = 0; i < hist_b.size(); i += 2)
          {
            second_moment_b += hist_b[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment_b += hist_b[i] * ((i + 2) / 2);
          }
          // std::cout << second_moment << " " << first_moment << endl;
          values_wei_b.push_back(second_moment_b / first_moment_b * 1.0 / double(P.N));
          weighted_b += second_moment_b / first_moment_b * 1.0 / double(P.N);
        }
        mean_b = mean_b / count_b;
        rel_mass_b = rel_mass_b / count_b;
        weighted_b = weighted_b / count_b;

        double cov_mot_b = 0;
        for (auto &val : values_mot_b)
        {
          cov_mot_b += pow(val - mean_b, 2);
        }
        cov_mot_b = cov_mot_b / (values_mot_b.size() - 1);

        double cov_mas_b = 0;
        for (auto &val : values_mas_b)
        {
          cov_mas_b += pow(val - rel_mass_b, 2);
        }
        cov_mas_b = cov_mas_b / (values_mas_b.size() - 1);

        for (auto &val : values_wei_b)
        {
          cov_w_b += pow(val - weighted_b, 2);
        }
        cov_w_b = cov_w_b / (values_wei_b.size() - 1);

        backward << al << " " << mean_b << " " << cov_mot_b << " " << rel_mass_b << " " << cov_mas_b << " " << weighted_b << " " << cov_w_b << endl;
        c++;
      }
    }
    else if (output == "perimeter")
    {
      ofstream outfile, pars;
      outfile.open("./lars_sim/Data/perimeter/hex_" + occ_p + "_t.txt");
      pars.open("./lars_sim/Data/perimeter/hex_pars_" + occ_p + "_t.txt");
      unsigned check = 0;
      for (double alp = 1e-3 * pow(2.4, 10); alp <= 100.0; alp *= 2.4)
      {
        for (double dens = 0.001; dens < 1.0; dens += .08)
        {
          Parameters P_h;
          P_h.N = unsigned(P.L[0] * P.L[0] * P.n_max * 2 * dens);
          P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          Hexagonal_lattice LB(P_h, rng);
          t = 0;

          // output for parameters (curious to see the values)
          if (check < 7)
          {
            pars << P_h.N << endl;
            check++;
          }

          // jamming, max cl size, perimeter, weighted mean cl size
          double j = 0;
          double m = 0;
          double s = 0;
          double w = 0;
          // count for mean vals
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            j += LB.motility_fraction();
            m += LB.max_cluster_size_nr();

            s += LB.perimeter();
            w += LB.w_N();

            count++;
          }

          j = j / count;
          m = m / count;
          s = s / count;
          w = w / count;

          m = m / double(P_h.N);
          s = s / double(P_h.N);
          w = w / double(P_h.N);

          outfile << j << " " << m << " " << s << " " << w << " ";
        }
        outfile << endl;
      }
    }
    else if (output == "motility")
    {
      if (details == 0)
      {
        ofstream outfile;
        string name = "./lars_sim/Data/motility/hexagonal_perc_testing";
        string outputname = name + "_" + occ_p + ".txt";
        outfile.open(outputname);
        for (double al = 0.0; al < 0.2; al += 0.005)
        {
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
          Hexagonal_lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          std::vector<double> values_wei;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          double weighted = 0;
          double cov_w = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
            mean += LB.motility_fraction();
            count++;
            // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distribution_particle_number();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i += 2)
            {
              second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
              first_moment += hist[i] * ((i + 2) / 2);
            }
            // std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0 / double(P.N));
            weighted += second_moment / first_moment * 1.0 / double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;

          double cov_mot = 0;
          for (auto &val : values_mot)
          {
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot / (values_mot.size() - 1);

          double cov_mas = 0;
          for (auto &val : values_mas)
          {
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas / (values_mas.size() - 1);

          for (auto &val : values_wei)
          {
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w / (values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
      }
      else if (details == 1)
      {
        ofstream outfile;
        string name = "./lars_sim/Data/motility/hexagonal_perc_details";
        string outputname = name + "_" + occ_p + ".txt";
        outfile.open(outputname);
        for (double al = 0.105; al < 0.12; al += 0.0005)
        {
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
          Hexagonal_lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr()) / double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr()) / double(P.N);
            mean += LB.motility_fraction();
            count++;
            // outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;
          }
          mean = mean / count;
          rel_mass = rel_mass / count;

          double cov_mot = 0;
          for (auto &val : values_mot)
          {
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot / (values_mot.size() - 1);

          double cov_mas = 0;
          for (auto &val : values_mas)
          {
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas / (values_mas.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << endl;
        }
      }
    }
    else if (output == "stable")
    {
      ofstream part, numb;
      part.open("./lars_sim/Data/stable/hexagonal.txt");
      numb.open("./lars_sim/Data/stable/hexagonal_number.txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        // only doing a positional output here
        numb << t << " " << HL.number_cluster() << endl;
        part << HexagonalParticleWriter(HL, part) << endl;
      }
    }
    else if (output == "weighted")
    {
      ofstream output, part;
      part.open("./lars_sim/Data/weighted/hex_part.txt");
      output.open("./lars_sim/Data/weighted/hex.txt");

      for (double n = 1; t < burnin + until; n++)
      {
        t = HL.run_until(burnin + n * every);
        part << HexagonalParticleWriter(HL, part) << endl;
        hist_t hist = HL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i += 2)
        {
          second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
          first_moment += hist[i] * ((i + 2) / 2);
        }
        // std::cout << second_moment << " " << first_moment << endl;
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        output << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
      }
    }
    else if (output == "area")
    {
      ofstream surf, part, clust, border;
      surf.open("./lars_sim/Data/surf/hex_sv" + occ_p + ".txt");
      part.open("./lars_sim/Data/surf/hex_part" + occ_p + ".txt");
      // clust.open("./lars_sim/Data/surf/hex_clust_low"+occ_p+".txt");
      border.open("./lars_sim/Data/surf/hex_border" + occ_p + ".txt");
      // unsigned so only one cluster each gets printed
      // unsigned check_1 = 0;
      // unsigned check_2 = 0;
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        vec a = HL.surface_volume();
        for (const auto &k : a)
          surf << k << " ";
        surf << endl;
        part << HexagonalParticleWriter(HL, part) << endl;
        for (const auto &n : HL.cluster_surface())
          border << n << " ";
        border << endl;

        /*
        if (HL.clust_size(8, 10) == 1 && check_1 == 0){
          vec n =HL.single_cluster(8, 10);
          for (const auto& m : n) clust << m << " ";
          clust << endl;
          check_1++;
        }

        if (HL.clust_size(9000, 10000) == 1 && check_2 == 0){
          vec n = HL.single_cluster(9000, 10000);
          for (const auto& m : n) clust << m << " ";
          clust << endl;
          check_2++;
        }
        */
      }
    }
    else if (output == "heatmap")
    {
      ofstream outfile, outfile_avg, outfile_nr, outfile_avg_nr;
      outfile.open("./lars_sim/heatmap/hex_alpha_N_n_2.txt");
      outfile_avg.open("./lars_sim/heatmap/hex_alpha_N_n_2_avg.txt");
      outfile_nr.open("./lars_sim/heatmap/hex_nr_alpha_N_n_2.txt");
      outfile_avg_nr.open("./lars_sim/heatmap/hex_nr_alpha_N_n_2_avg.txt");
      for (double alp = 1e-6; alp <= 1.0; alp *= 1.3)
      {
        for (unsigned N = 200 * P.n_max; N <= 2 * P.L[0] * P.L[0] * P.n_max; N += 400 * P.n_max)
        {
          Parameters P_h;
          P_h.N = N;
          P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
          P_h.L = P.L;
          P_h.n_max = P.n_max;
          std::size_t maxsize = 1;
          std::size_t maxsize_nr = 1;
          Hexagonal_lattice LB(P_h, rng);
          t = 0;

          // also taking mean over all used timesteps
          double mean = 0.0;
          double mean_nr = 0.0;
          double count = 0;
          for (unsigned n = 0; t < burnin + until; ++n)
          {
            t = LB.run_until(burnin + n * every);
            maxsize = std::max(maxsize, LB.max_cluster_size());
            maxsize_nr = std::max(maxsize_nr, LB.max_cluster_size_nr());
            mean += LB.avg_cluster_size();
            mean_nr += LB.avg_cluster_size_nr();
            count += 1;
          }

          mean = mean / count;
          mean_nr = mean_nr / count;

          outfile << maxsize << " ";
          outfile_avg << mean << " ";
          outfile_nr << maxsize_nr << " ";
          outfile_avg_nr << mean_nr << " ";
        }
        outfile << endl;
        outfile_avg << endl;
        outfile_nr << endl;
        outfile_avg_nr << endl;
      }
    }
    else if (output == "number")
    {
      ofstream outfile;
      string name = "./lars_sim/Data/trajectory/hex_";
      string name_number = name + "number" + "_" + tumb + "_" + dens + "_" + size + "_" + occ_p + txt;
      outfile.open(name_number);

      vec_d res, clnr, cn, j, m;
      for (unsigned n = 0; burnin + pow(1.01, n) * every < burnin + until; n++)
      {
        res.push_back(0);
        clnr.push_back(0);
        cn.push_back(0);
        j.push_back(0);
        m.push_back(0);
      }

      unsigned iterations = 10;

      for (int q = 0; q <= iterations; q++)
      {

        t = HL.run_until(burnin);
        hist_t hist = HL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i += 2)
        {
          second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
          first_moment += hist[i] * ((i + 2) / 2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        res[0] += weighted;
        clnr[0] += HL.avg_cluster_size_nr();
        cn[0] += HL.perimeter();
        j[0] += HL.motility_fraction();
        m[0] += HL.max_cluster_size_nr();
        for (unsigned n = 1; t < burnin + until; n++)
        {
          t = HL.run_until(burnin + pow(1.01, n) * every);
          hist_t hist = HL.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i += 2)
          {
            second_moment += hist[i] * ((i + 2) / 2) * ((i + 2) / 2);
            first_moment += hist[i] * ((i + 2) / 2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          res[n] += weighted;
          clnr[n] += HL.avg_cluster_size_nr();
          cn[n] += HL.perimeter();
          j[n] += HL.motility_fraction();
          m[n] += HL.max_cluster_size_nr();
        }
        cout << "We are at " << q << endl;
      }
      for (unsigned i = 0; i <= res.size(); i++)
      {
        outfile << burnin + pow(1.01, i) * every << " " << res[i] / double(iterations + 1) << " " << cn[i] / double(iterations + 1) << " " << j[i] / double(iterations + 1) << " " << m[i] / double(iterations + 1) << " " << double(clnr[i]) / double(iterations + 1) << endl;
      }
    }
    else if (output == "snapshots")
    {
      ofstream outfile;
      outfile.open("./lars_sim/gif/hexdir.txt");

      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        // only doing a positional output here
        HL.realise_directions();
        outfile << HexDirectionWriter(HL, outfile) << endl;
      }
    }
    else if (output == "stopping time")
    {
      ofstream outfile;

      outfile.open("./lars_sim/Data/stopping/hex_" + occ_p + "_" + alpha_p + ".txt");
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        vec_d st = HL.stopping(t);
        for (const auto &m : st)
          outfile << m << " ";
        outfile << endl;
      }
    }
    else if (output == "function")
    {
      ofstream outfile;
      outfile.open("./lars_sim/testing/hexagonal.txt");


      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        // only doing a positional output here
        outfile << HexagonalParticleWriter(HL, outfile) << endl;
      }
    }
    else if (output == "distribution")
    {
      ofstream outfile, outfile2;
      outfile.open("./lars_sim/Data/dist/hex_" + occ_p + ".txt");
      outfile2.open("./lars_sim/Data/dist/hex_dens_" + occ_p + ".txt");
      hist_t dist(4);
      for (double n = 0; t < burnin + until; n++)
      {
        t = HL.run_until(burnin + n * every);
        vec_d dens = HL.density();
        for (const auto &m : dens)
          outfile2 << m << " ";
        outfile2 << endl;
        hist_t dr = HL.particle_neighbour_dist();
        dist[std::slice(0, dr.size(), 1)] += dr;
      }
      for (const auto &m : dist)
        outfile << m << " ";
      outfile << endl;
    }
    else if (output == "correlation")
    {
      ofstream outfile;
      outfile.open("./lars_sim/Data/corr/hex_" + occ_p + "_" + alpha_p + ".txt");
      for (unsigned n = 0; t < burnin + until; ++n)
      {
        t = HL.run_until(burnin + n * every);
        vec output = HL.occ_array();
        for (const auto &m : output)
          outfile << m << " ";
        outfile << endl;
      }
    }
  }
  return 0;
}
