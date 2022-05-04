/**
  * crumble: condensed run-and-tumble simulation
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

#include "CLI11.hpp"
#include "Scheduler.hpp"

using namespace std;
using direction_t = unsigned char;
using hist_t = std::valarray<unsigned>;

// needed for calculating position of particles on hexagonal lattice
// TODO see if needed somewhere else, otherwise this can be moved to HexagonalParticleWriter
const double pi = std::acos(-1);

/**
  * Model parameters. We put these in a struct so we can pass as a single
  * argument to Lattice and also so all the defaults are in one place.
  */
struct Parameters {
  std::vector<unsigned>   L {100};  // Lengths of each lattice dimension
  unsigned                N=50;     // Number of particles
  std::vector<double> alpha {1.0};  // Scaled tumble rate in each dimension (actual rate is alpha/d)
  // ! Addition by Lars
  unsigned                n_max=1;  // Max occupation number on each lattice site
};

// We need to forward declare these so the compiler can cope with the friend declarations
template<typename Engine> class SnapshotWriter;
template<typename Engine> class VacancyWriter;
template<typename Engine> class ParticleWriter;
template<typename Engine> class ClusterWriter;

template<typename Engine> class TriangleParticleWriter;
template<typename Engine> class HexagonalParticleWriter;

template<typename Engine>
class Lattice {

  // For output in different formats
  friend class SnapshotWriter<Engine>;
  friend class VacancyWriter<Engine>;
  friend class ParticleWriter<Engine>;
  friend class ClusterWriter<Engine>;

  // Data associated with each site; by default all of these are set to zero
  // which represents a vacant site
  struct Site {
    unsigned id; // Particle / vacancy id
    bool occupied; // There is a particle here
    bool active; // A move event is scheduled
    direction_t neighbours; // Number of neighbours that are occupied
    direction_t direction; // Direction of last hop attempt
    double hoptime; // Time of last hop attempt
  };

  Parameters P; // A local copy of the model parameters
  std::vector<Site> sites; // Representation of the sites
  Engine& rng; // Source of noise: this is a reference as there should only be one of these!
  std::discrete_distribution<unsigned> anyway; // Distribution over tumble directions
  std::exponential_distribution<double> run, tumble; // Distribution of times between run and tumble events
  Scheduler S; // Keeps track of the event queue

  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours. We have periodic boundary conditions
  auto neighbours(unsigned n) const {
    std::vector<unsigned> nbs(2*P.L.size());
    unsigned below = 1;
    for(unsigned d=0; d<P.L.size(); ++d) {
      unsigned L = P.L[d];
      unsigned above = below * L;
      // x is the position along dimension d
      // y is the contribution to the site index from all other dimensions
      unsigned x = (n/below) % L;
      unsigned y = (n%below) + (n/above) * above;
      // Neighbours in the increasing and decreasing directions
      nbs[2*d] = y + ((x+1)%L)*below;
      nbs[2*d+1] = y + ((x+L-1)%L)*below;
      below = above;
    }
    return nbs;
  }

  // Given an index into sites, return a sequence of indices corresponding to
  // its neighbours in the forward direction along each axis.
  // We still have periodic boundary conditions
  auto forward_neighbours(unsigned n) const {
    std::vector<unsigned> nbs(P.L.size());
    unsigned below = 1;
    for(unsigned d=0; d<P.L.size(); ++d) {
      unsigned L = P.L[d];
      unsigned above = below * L;
      // x is the position along dimension d
      // y is the contribution to the site index from all other dimensions
      unsigned x = (n/below) % L;
      unsigned y = (n%below) + (n/above) * above;
      // Neighbours in the increasing and decreasing directions
      nbs[d] = y + ((x+1)%L)*below;
      below = above;
    }
    return nbs;
  }

  // Place a particle with given direction and hop time at site n;
  // neighbouring sites will be accordingly adjusted
  void place(unsigned n, unsigned id, direction_t d, double t) {
    sites[n].id = id;
    sites[n].direction = d;
    sites[n].hoptime = t;
    if(!sites[n].occupied) {
      sites[n].occupied = true;
      for(const auto&m: neighbours(n)) ++sites[m].neighbours;
    }
  }

  // Schedule a hop event for a particle at site n
  void schedule(unsigned n) {
    assert(sites[n].occupied);
    S.schedule(run(rng), [this,n]() {
      assert(sites[n].occupied);
      assert(sites[n].active);
      // If there are no local vacancies, mark this particle as inactive and exit
      if(sites[n].neighbours == 2*P.L.size()) {
        // std::cout << "Can't move from "; decode(n); std::cout << " deactivating" << std::endl;
        sites[n].active = false;
      } else {
        // if(std::uniform_real_distribution<double>()(rng)>=std::exp(-P.alpha*(S.time()-sites[n].hoptime))) {
        if(tumble(rng) < S.time() - sites[n].hoptime) {
          // At least one tumble event happened since we last attempted a hop; so randomise the direction
          sites[n].direction = anyway(rng);
          //std::cout << "Now moving in direction " << int(sites[n].direction) << " from "; decode(n); std::cout << std::endl;
        }
        sites[n].hoptime = S.time();
        // Get the sites adjacent to the departure site
        auto dnbs = neighbours(n);
        if(!sites[dnbs[sites[n].direction]].occupied) {
          assert(!sites[dnbs[sites[n].direction]].active);
          // Get the id of the vacancy that is being displaced
          unsigned vid = sites[dnbs[sites[n].direction]].id;
          // std::cout << "Moving from "; decode(n); std::cout << " deactivating" << std::endl;
          // Deactive the departure site; also mark it empty
          sites[n].occupied = sites[n].active = false;
          // Place a particle on the target site; it has the same direction and hoptime as the departing particle
          // std::cout << "Moving to "; decode(dnbs[sites[n].direction]); std::cout << " placing" << std::endl;
          place(dnbs[sites[n].direction], sites[n].id, sites[n].direction, sites[n].hoptime);
          // Move the vacancy id onto the departure site
          sites[n].id = vid;
          // Now go through the neighbours of the departure site, update neighbour count and activate any
          // that can now move. Note the particle this is at the target site is included in this list
          // and will be activated accordingly
          for(const auto& m: dnbs) {
            --sites[m].neighbours;
            if(sites[m].occupied && !sites[m].active) schedule(m);
          }
        } else {
          // std::cout << "Didn't move from "; decode(n); std::cout << std::endl;
          // This site is still active, so schedule another hop
          schedule(n);
        }
      }
    });
    sites[n].active = true;
    // std::cout << "Scheduled "; decode(n); std::cout << std::endl;
  }


  // For testing
  void decode(unsigned n) {
    std::cout << "[ ";
    for(const auto& L: P.L) {
      std::cout << (n % L) << " ";
      n /= L;
    }
    std::cout << "]";
  }

  // Check the lattice state is consistent
  bool consistent() {
    unsigned active = 0, occupied = 0;
    std::set<unsigned> ids;
    for(unsigned n=0; n<sites.size(); ++n) {
      // Check each site has a unique id
      if(ids.count(sites[n].id)) return false;
      ids.insert(sites[n].id);
      // Check that empty sites are also inactive
      if(!sites[n].occupied) {
        if(sites[n].active) return false;
        // Nothing left to do if empty
        continue;
      }
      // Check that the neighbour count is correct
      ++occupied;
      unsigned nbs = 0;
      for(const auto& m: neighbours(n)) {
        if(sites[m].occupied) ++nbs;
      }
      if(nbs != sites[n].neighbours) return false;
      // Check that mobile particles are active
      if(nbs < 2*P.L.size() && !sites[n].active) return false;
      if(sites[n].active) ++active;
    }
    // Check we've not lost any particles
    return occupied == P.N && active == S.pending();
  }


public:

  Lattice(const Parameters& P, Engine& rng) :
    P(P), // NB: this takes a copy
    sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
    rng(rng), // NB: this takes a reference
    tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0)/P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
  {
    // Set up the tumble direction distribution
    std::vector<double> drates(2*P.L.size());
    for(unsigned d=0;d < P.L.size(); ++d) {
       drates[2*d] = drates[2*d+1] = d<P.alpha.size() ? P.alpha[d]/tumble.lambda() : 1.0;
    }
    anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());

    // Place particles on the lattice
    unsigned unplaced = P.N; // Current number of particles remaining to be placed
    for(unsigned n=0; n<sites.size(); ++n) {
      // Number of remaining sites where partcles could be placed is sites.size()-n, unplaced of which need to be filled
      if(std::uniform_int_distribution<unsigned>(1,sites.size()-n)(rng) <= unplaced) {
        // Particles get numbered from 0
        place(n, P.N-unplaced, anyway(rng), 0.0);
        --unplaced;
      } else {
        // Vacancies get numbered from N
        sites[n].id = unplaced + n;
      }
    }
    assert(unplaced==0);

    // Activate particles that can move, and schedule a hop accordingly
    for(unsigned n=0; n<sites.size(); ++n) {
      if(sites[n].occupied && sites[n].neighbours < 2*P.L.size()) schedule(n);
    }
    assert(consistent());
  }

  // Iterates the simulation until the given time; returns the actual time run to
  double run_until(double time) {
    while(S.advance(time));
    assert(consistent());
    return S.time();
  }

  // Sample all particle directions from the distribution that applies at the current instant
  // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
  void realise_directions() {
    for(auto& site: sites) {
      if(!site.occupied) continue;
      if(tumble(rng) < S.time() - site.hoptime) {
        // At least one tumble event happened since we last attempted a hop; so randomise the direction
        site.direction = anyway(rng);
      }
      // Advance the hop time so we don't generate temporal paradox
      site.hoptime = S.time();
    }
  }

  // Return cluster size distributions.
  // Element 2n   contains the number of clusters of particles of size n+1
  // Element 2n+1 contains the number of clusters of vacancies of size n+1
  hist_t cluster_distributions() const {

    // Lookup table of cluster membership by lattice site
    std::vector<unsigned> memberof(sites.size());
    // Initially, this is just the site id as each site is its own cluster
    std::iota(memberof.begin(), memberof.end(), 0);
    // Create also a map of clusters each containing a list of its members
    std::map<unsigned, std::list<unsigned>> clusters;
    for(unsigned n=0; n<sites.size(); ++n) clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

    // Keep track of the size of the largest cluster
    std::size_t maxsize = 1;

    for(unsigned n=0; n<sites.size(); ++n) {
      // Loop over neigbours m in one direction only so we visit each bond once
      for(const auto& m: forward_neighbours(n)) {
        unsigned large = memberof[n], small = memberof[m];
        // Merge clusters at n and m if they both have the same occupation, and are not already part of the same cluster
        if(sites[n].occupied == sites[m].occupied && small != large) {
          // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
          if(clusters[large].size() < clusters[small].size()) std::swap(large,small);
          // Update the cluster number for all sites in the smaller one
          for(const auto& site: clusters[small]) memberof[site] = large;
          // Add the members of the smaller cluster onto the end of the larger one
          clusters[large].splice(clusters[large].end(), clusters[small]);
          // Remove the smaller cluster from the map
          clusters.erase(small);
          // Keep track of the largest cluster
          maxsize = std::max(maxsize, clusters[large].size());
        }
      }
    }

    hist_t dists(2*maxsize);
    for(const auto& kv: clusters) {
      dists[2*kv.second.size()-1-sites[kv.first].occupied]++;
    }
    return dists;
  }

};

// ! Following classes are additions by Lars. Note that in the triangular case it is just the one as above with some very slight adjustments
// ! while the one for the hexagonal lattice has some more drastic changes. 


// TODO FIGURE OUT WHERE THE SEGFAULT IS?!?!??!?!?
template<typename Engine>
class Triangle_lattice {

    // For output in different formats
    friend class SnapshotWriter<Engine>;
    friend class VacancyWriter<Engine>;
    friend class TriangleParticleWriter<Engine>;
    friend class ClusterWriter<Engine>;

    Parameters P; // A local copy of the model parameters

    static constexpr unsigned n_max = 2; // ! Nicer way of doing this?
    
    

    // Data associated with each site; by default all of these are set to zero
    // which represents a vacant site
    struct Site {
        std::vector<unsigned> id = std::vector<unsigned>(n_max); // Particle / vacancy id
        std::vector<bool> occupied = std::vector<bool>(n_max); // There is a particle here
        std::vector<bool> active = std::vector<bool>(n_max); // A move event is scheduled
        direction_t neighbours; // Number of neighbours that are occupied
        std::vector<direction_t> direction = std::vector<direction_t>(n_max); // Direction of last hop attempt
        std::vector<double> hoptime = std::vector<double>(n_max); // Time of last hop attempt
        unsigned present; // Number of particles present at site. Has to be <= n_max
    };

    std::vector<Site> sites; // Representation of the sites
    Engine& rng; // Source of noise: this is a reference as there should only be one of these!
    std::discrete_distribution<unsigned> anyway; // Distribution over tumble directions
    std::exponential_distribution<double> run, tumble; // Distribution of times between run and tumble events
    Scheduler S; // Keeps track of the event queue

    // Given an index into sites, return a sequence of indices corresponding to
    // its neighbours. We have periodic boundary conditions
    auto neighbours(unsigned n) const {

        // this is only valid for 2d
        int L = P.L[0];
        std::vector<int> nbs(2 * P.L.size() + 2);

        int x = n % L;
        int y = (n/L);

        // for diagonals. There surely is a better way than this
        int xk, yk, xm, ym;
        int k = n + L + 1;
        int m = n - L - 1;

        xk = k%L;
        xm = m%L;

        if (((k-1)/L - L) >= 0){ 
          yk = 0;
        }else{
          yk = ((k-1)/L);
        }
        if ((m+1) < 0){ 
          ym = L - 1;
        }else{
          int l = n - L;
          if (l < 0) l = -l;
          ym = (l/L);
        }
        nbs[0] = (n + 1) % L + y * L;
        nbs[1] = (n - 1) % L + y * L;

        nbs[2] = x + ((y + 1) % L) * L;
        nbs[3] = x + ((y - 1 + L) % L) * L;

        nbs[4] = xk + yk * L;
        nbs[5] = xm + ym * L;

        return nbs;
    }

    // Given an index into sites, return a sequence of indices corresponding to
    // its neighbours in the forward direction along each axis.
    // We still have periodic boundary conditions
    auto forward_neighbours(unsigned n) const {
        // this is only valid for 2d
        int L = P.L[0];
        std::vector<int> nbs(P.L.size() + 1);



        int x = n % L;
        int y = (n/L);

        // for diagonals. There surely is a better way than this
        int xk, yk;
        int k = n + L + 1;

        xk = k%L;

        if (((k-1)/L - L) >= 0){ 
          yk = 0;
        }else{
          yk = ((k-1)/L);
        }
        nbs[0] = (n + 1) % L + y * L;

        nbs[1] = x + ((y + 1) % L) * L;

        nbs[2] = xk + yk * L;

        return nbs;
    }

    // Place a particle with given direction and hop time at site n;
    // neighbouring sites will be accordingly adjusted
    void place(unsigned n, unsigned id, direction_t d, double t, unsigned index) {
        sites[n].id[index] = id;
        sites[n].direction[index] = d;
        sites[n].hoptime[index] = t;
        
        if (!sites[n].occupied[index]) {
            sites[n].occupied[index] = true;
            for (const auto& m : neighbours(n)){
             ++sites[m].neighbours;
            }
        }
        
        sites[n].present++;
    }

    // Schedule a hop event for a particle at site n
    void schedule(unsigned n, unsigned index) {
        assert(sites[n].occupied[index]);
        S.schedule(run(rng), [this, n, index]() {
            assert(sites[n].occupied[index]);
            assert(sites[n].active[index]);
            // If there are no local vacancies, mark this particle as inactive and exit
            if (sites[n].neighbours == 6 * P.n_max) {
                // std::cout << "Can't move from "; decode(n); std::cout << " deactivating" << std::endl;
                sites[n].active[index] = false;
            }
            else {
                  // if(std::uniform_real_distribution<double>()(rng)>=std::exp(-P.alpha*(S.time()-sites[n].hoptime))) {
                if (tumble(rng) < S.time() - sites[n].hoptime[index]) {
                    // At least one tumble event happened since we last attempted a hop; so randomise the direction
                    sites[n].direction[index] = anyway(rng);
                    //std::cout << "Now moving in direction " << int(sites[n].direction) << " from "; decode(n); std::cout << std::endl;
                }
                sites[n].hoptime[index] = S.time();
                // Get the sites adjacent to the departure site
                
                auto dnbs = neighbours(n);
                if (sites[dnbs[sites[n].direction[index]]].present < P.n_max) {
                  auto itr = std::lower_bound(sites[n].occupied.begin(), sites[n].occupied.end(), false);
                  unsigned k = std::distance(sites[n].occupied.begin(), itr);
                  if (k <= (P.n_max-1) && !sites[dnbs[sites[n].direction[index]]].occupied[k]){
                     
                    // Get the id of the vacancy that is being displaced
                    unsigned vid = sites[dnbs[sites[n].direction[index]]].id[k];
                    // std::cout << "Moving from "; decode(n); std::cout << " deactivating" << std::endl;
                    // Deactive the departure site; also mark it empty
                    sites[n].occupied[index] = sites[n].active[index] = false;
                    
                    // as it moves from n, we have one less present at n
                    sites[n].present--;
                    
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
            }
            });
        sites[n].active[index] = true;
        // std::cout << "Scheduled "; decode(n); std::cout << std::endl;
    }


    // For testing
    void decode(unsigned n) {
        std::cout << "[ ";
        for (const auto& L : P.L) {
            std::cout << (n % L) << " ";
            n /= L;
        }
        std::cout << "]";
    }

    // Check the lattice state is consistent
    bool consistent() {
        unsigned active = 0, occupied = 0;
        std::set<unsigned> ids;
        for (unsigned n = 0; n < sites.size(); ++n) {
          for (unsigned k = 0; k < P.n_max; k++){
            // Check each site has a unique id
            if (ids.count(sites[n].id[k])) return false;
            ids.insert(sites[n].id[k]);
            // Check that empty sites are also inactive
            if (!sites[n].occupied[k]) {
                if (sites[n].active[k]) return false;
                // Nothing left to do if empty
                continue;
            }
            // Check that the neighbour count is correct
            ++occupied;
            // ! Can we move this out one loop?
            unsigned nbs = 0;
            for (const auto& m : neighbours(n)) {
              for (unsigned i = 0; i < P.n_max; i++){
                if (sites[m].occupied[i]) ++nbs;
              }
            }
            if (nbs != sites[n].neighbours) return false;
            // Check that mobile particles are active
            if (nbs < 6*P.n_max && !sites[n].active[k]) return false;
            if (sites[n].active[k]) ++active;
          }
        }
        // Check we've not lost any particles
        return occupied == P.N && active == S.pending();
    }


public:

    Triangle_lattice(const Parameters& P, Engine& rng) :
        P(P), // NB: this takes a copy
        sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
        rng(rng), // NB: this takes a reference
        tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
    {
        // Set up the tumble direction distribution
        std::vector<double> drates(2 * P.L.size() + 2);
        for (unsigned d = 0; d < 3; ++d) {
            drates[2 * d] = drates[2 * d + 1] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
        }
        anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());
        

        // Place particles on the lattice
        unsigned unplaced = P.N; // Current number of particles remaining to be placed
        for (unsigned n = 0; n < sites.size(); ++n) {
            // Number of remaining sites where partcles could be placed is sites.size()-n, unplaced of which need to be filled
            if (std::uniform_int_distribution<unsigned>(1, sites.size() - n)(rng) <= unplaced) {
                // For ease we only place one particle per site in the initial configuration
                place(n, P.N - unplaced, anyway(rng), 0.0, 0);
                --unplaced;
                // TODO better way of doing this? We need (n_max - 1)*L**2 more id's for vacancies
                for (unsigned k = 1; k < P.n_max; k++){
                  sites[n].id[k] = n + k * P.L[0] * P.L[1];
                }
            }
            // vacancies
            else {
              for (unsigned k = 0; k < P.n_max; k++){
                if (k == 0) sites[n].id[k] = unplaced + n;
                else sites[n].id[k] = n + k * P.L[0] * P.L[1];
              }
            }
        }
        assert(unplaced == 0);
        
        // Activate particles that can move, and schedule a hop accordingly
        for (unsigned n = 0; n < sites.size(); ++n) {
          for (unsigned k = 0; k < P.n_max; ++k){
            if (sites[n].occupied[k] && sites[n].neighbours < 3 * P.n_max) schedule(n, k);
          }
        }
        assert(consistent());
    }

    // Iterates the simulation until the given time; returns the actual time run to
    double run_until(double time) {
        while (S.advance(time));
        assert(consistent());
        return S.time();
    }

    // Sample all particle directions from the distribution that applies at the current instant
    // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
    void realise_directions() {
        for (auto& site : sites) {
          for (unsigned k = 0; k < P.n_max; k++){
            if (!site.occupied[k]) continue;
            if (tumble(rng) < S.time() - site.hoptime[k]) {
                // At least one tumble event happened since we last attempted a hop; so randomise the direction
                site.direction[k] = anyway(rng);
            }
            // Advance the hop time so we don't generate temporal paradox
            site.hoptime[k] = S.time();
          }
        }
    }

    // Return cluster size distributions.
    // Element 2n   contains the number of clusters of particles of size n+1
    // Element 2n+1 contains the number of clusters of vacancies of size n+1
    hist_t cluster_distributions() const {
        // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof.begin(), memberof.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters;
        for (unsigned n = 0; n < sites.size(); ++n) clusters[n] = std::list<unsigned>(1, n); // Single-element list comprising the lattice site

        // Keep track of the size of the largest cluster
        std::size_t maxsize = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof[n], small = memberof[m];
                // Merge clusters at n and m if they both have the same occupation, and are not already part of the same cluster
                // TODO check whether the [0] argument actually works. 
                if (sites[n].occupied[0] == sites[m].occupied[0] && small != large) {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters[large].size() < clusters[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters[small]) memberof[site] = large;
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
        for (const auto& kv : clusters) {
            for (unsigned i = 0; i < P.n_max; i++){
              dists[2 * kv.second.size() - 1 - sites[kv.first].occupied[i]]++;
            }
        }
        return dists;
    }

};

template<typename Engine>
class Hexagonal_lattice {

    // For output in different formats
    friend class SnapshotWriter<Engine>;
    friend class VacancyWriter<Engine>;
    friend class HexagonalParticleWriter<Engine>;
    friend class ClusterWriter<Engine>;

    Parameters P; // A local copy of the model parameters
    

    // Data associated with each site; by default all of these are set to zero
    // which represents a vacant site
    struct Site {
        std::vector<unsigned> id  = std::vector<unsigned>(2); // Particle / vacancy id
        std::vector<bool> occupied = std::vector<bool>(2); // There is a particle here
        std::vector<bool> active = std::vector<bool>(2); // A move event is scheduled
        std::vector<direction_t> neighbours = std::vector<direction_t>(2); // Number of neighbours that are occupied
        std::vector<direction_t> direction = std::vector<direction_t>(2); // Direction of last hop attempt
        std::vector<double> hoptime = std::vector<double>(2); // Time of last hop attempt
    };


    std::vector<Site> sites; // Representation of the sites
    Engine& rng; // Source of noise: this is a reference as there should only be one of these!
    std::discrete_distribution<unsigned> anyway; // Distribution over tumble directions
    std::exponential_distribution<double> run, tumble; // Distribution of times between run and tumble events
    Scheduler S; // Keeps track of the event queue
    // Given an index into sites, return a sequence of indices corresponding to
    // its neighbours. We have periodic boundary conditions
    auto neighbours(unsigned n, unsigned index) const {

        // this is only valid for 2d

        std::vector<unsigned> nbs(3);
        unsigned L1 = P.L[0];
        unsigned L2 = P.L[1];

        unsigned x = n % L1;
        unsigned y = (n/L1);

        if (index){

          nbs[0] = n;                         // same site
          nbs[1] = (n + 1) % L1 + y * L1;     // right
          nbs[2] = x + ((y - 1) % L2) * L1;   // down

        } else{
          
          nbs[0] = n;                         // same site
          nbs[1] = (n - 1) % L1 + y * L1;     // left
          nbs[2] = x + ((y + 1) % L2) * L1;   // up
          

        }

        return nbs;
    }

    // Place a particle with given direction and hop time at site n;
    // neighbouring sites will be accordingly adjusted
    void place(unsigned n, unsigned id, direction_t d, double t, unsigned index) {
        sites[n].id[index] = id;
        sites[n].direction[index] = d;
        sites[n].hoptime[index] = t;
        if (!sites[n].occupied[index]) {
            sites[n].occupied[index] = true;
            for (const auto& m : neighbours(n, index)) ++sites[m].neighbours[(index+1)%2];
        }
    }

    // Schedule a hop event for a particle at site n
    void schedule(unsigned n, unsigned index) {
        assert(sites[n].occupied[index]);
        S.schedule(run(rng), [this, n, index]() {
            assert(sites[n].occupied[index]);
            assert(sites[n].active[index]);
            // If there are no local vacancies, mark this particle as inactive and exit
            if (sites[n].neighbours[index] == 3) {
                // std::cout << "Can't move from "; decode(n); std::cout << " deactivating" << std::endl;
                sites[n].active[index] = false;
            }
            else {
                // if(std::uniform_real_distribution<double>()(rng)>=std::exp(-P.alpha*(S.time()-sites[n].hoptime))) {
                if (tumble(rng) < S.time() - sites[n].hoptime[index]) {
                    sites[n].direction[index] = anyway(rng);
                }
                sites[n].hoptime[index] = S.time();
                // Get the sites adjacent to the departure site
                auto dnbs = neighbours(n, index);
                if (!sites[dnbs[sites[n].direction[index]]].occupied[(index+1)%2]) {
                    assert(!sites[dnbs[sites[n].direction[index]]].active[(index+1)%2]);
                    // Get the id of the vacancy that is being displaced
                    unsigned vid = sites[dnbs[sites[n].direction[index]]].id[(index+1)%2];
                    // std::cout << "Moving from "; decode(n); std::cout << " deactivating" << std::endl;
                    // Deactive the departure site; also mark it empty
                    sites[n].occupied[index] = sites[n].active[index] = false;
                    // Place a particle on the target site; it has the same direction and hoptime as the departing particle
                    // std::cout << "Moving to "; decode(dnbs[sites[n].direction]); std::cout << " placing" << std::endl;
                    place(dnbs[sites[n].direction[index]], sites[n].id[index], sites[n].direction[index], sites[n].hoptime[index], (index+1)%2);
                    // Move the vacancy id onto the departure site
                    sites[n].id[index] = vid;
                    // Now go through the neighbours of the departure site, update neighbour count and activate any
                    // that can now move. Note the particle this is at the target site is included in this list
                    // and will be activated accordingly
                    for (const auto& m : dnbs) {
                        --sites[m].neighbours[(index+1)%2];
                        if (sites[m].occupied[(index+1)%2] && !sites[m].active[(index+1)%2]) schedule(m, (index+1)%2);
                    }
                }
                else {
                    // std::cout << "Didn't move from "; decode(n); std::cout << std::endl;
                    // This site is still active, so schedule another hop
                    schedule(n, index);
                }
            }
            });
        sites[n].active[index] = true;
        // std::cout << "Scheduled "; decode(n); std::cout << std::endl;
    }


    // For testing
    void decode(unsigned n) {
        std::cout << "[ ";
        for (const auto& L : P.L) {
            std::cout << (n % L) << " ";
            n /= L;
        }
        std::cout << "]";
    }

    // Check the lattice state is consistent
    bool consistent() {
        unsigned active = 0, occupied = 0;
        std::set<unsigned> ids;
        for (unsigned n = 0; n < sites.size(); ++n) {
            // Check each site has a unique id
            for (unsigned i = 0; n < 2; i++){
              if (ids.count(sites[n].id[i])) return false;
              ids.insert(sites[n].id[i]);
              // Check that empty sites are also inactive
              if (!sites[n].occupied[i]) {
                  if (sites[n].active[i]) return false;
                  // Nothing left to do if empty
                  continue;
              }
              // Check that the neighbour count is correct
              ++occupied;
              unsigned nbs = 0;
              for (const auto& m : neighbours(n, i)) {
                  if (sites[m].occupied[(i+1)%2]) ++nbs;
              }
              if (nbs != sites[n].neighbours[i]) return false;
              // Check that mobile particles are active
              if (nbs < 3 && !sites[n].active[i]) return false;
              if (sites[n].active[i]) ++active;
            }
        }
        // Check we've not lost any particles
        return occupied == P.N && active == S.pending();
    }


public:

    Hexagonal_lattice(const Parameters& P, Engine& rng) :
        P(P), // NB: this takes a copy
        sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
        rng(rng), // NB: this takes a reference
        tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
    {
        // Set up the tumble direction distribution
        std::vector<double> drates(3);
        for (unsigned d = 0; d < 3; ++d) {
            drates[d] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
        }
        anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end()); 

        // Place particles on the lattice
        unsigned unplaced = P.N; // Current number of particles remaining to be placed
        for (unsigned n = 0; n < sites.size(); ++n) {
            // Number of remaining sites where partcles could be placed is sites.size()-n, unplaced of which need to be filled
            if (std::uniform_int_distribution<unsigned>(1, sites.size() - n)(rng) <= unplaced) {
                // Check whether to place it on 0 or on 1
                if (std::uniform_int_distribution<int>(1, 2)(rng) < 2){
                  // place particle at 0
                  // Particles get numbered from 0
                  place(n, P.N - unplaced, anyway(rng), 0.0, 0);
                  --unplaced;
                  // TODO better way of doing this? We need L**2 more id's for vacancies
                  sites[n].id[1] = n + P.L[0] * P.L[1];
                } else{
                  // place particle at 1
                  // Particles get numbered from 0
                  place(n, P.N - unplaced, anyway(rng), 0.0, 1);
                  --unplaced;
                  sites[n].id[0] = n + P.L[0] * P.L[1];
                }
                
            }
            else {
                // don't need to randomize these, since they are vacancies
                sites[n].id[0] = unplaced + n;
                sites[n].id[1] = n + P.L[0] * P.L[1];
            }
        }
        assert(unplaced == 0);

        // Activate particles that can move, and schedule a hop accordingly
        for (unsigned n = 0; n < sites.size(); ++n) {
          for (unsigned i = 0; i < 2; i++){
            if (sites[n].occupied[i] && sites[n].neighbours[i] < 3) schedule(n, i);
          }
        }
        assert(consistent());
    }

    // Iterates the simulation until the given time; returns the actual time run to
    double run_until(double time) {
        while (S.advance(time));
        assert(consistent());
        return S.time();
    }

    // Sample all particle directions from the distribution that applies at the current instant
    // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
    void realise_directions() {
        for (auto& site : sites) {
          for (unsigned i = 0; i < 2; i++){
            if (!site.occupied[i]) continue;
            if (tumble(rng) < S.time() - site.hoptime[i]) {
                // At least one tumble event happened since we last attempted a hop; so randomise the direction
                site.direction[i] = anyway(rng);
            }
            // Advance the hop time so we don't generate temporal paradox
            site.hoptime[i] = S.time();
          }
        }
    }

    // Return cluster size distributions.
    // Element 2n   contains the number of clusters of particles of size n+1
    // Element 2n+1 contains the number of clusters of vacancies of size n+1
    hist_t cluster_distributions() const {

        // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof(2 * sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof.begin(), memberof.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters;
        for (unsigned n = 0; n < sites.size(); ++n){
          for (unsigned i = 0; i < 2; i++){
            clusters[2*n + i] = std::list<unsigned>(1, 2*n + i); // Single-element list comprising the lattice site
          }
        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
          for (unsigned i = 0; i < 2; i++){
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // Merge clusters at n and m if they both have the same occupation, and are not already part of the same cluster
                if (sites[n].occupied[i] == sites[m].occupied[(i+1)%2] && small != large) {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters[large].size() < clusters[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters[small]) memberof[site] = large;
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
        for (const auto& kv : clusters) {
          for (unsigned i = 0; i < 2; i++){
            dists[2 * kv.second.size() - 1 - sites[kv.first].occupied[i]]++;
          }
        }
        return dists;
    }

};



template<typename Engine>
class SnapshotWriter {
  const Lattice<Engine>& L;
public:
  SnapshotWriter(const Lattice<Engine>& L) : L(L) { }

  // Output pairs (d, a) indicating the direction d (numbered from 1) or 0 if site is vacant
  // and a indicating if the site is active

  friend std::ostream& operator << (std::ostream& out, const SnapshotWriter& SW) {

    for(const auto& site: SW.L.sites) {
      out << (site.occupied ? unsigned(site.direction+1) : 0);
      out << " " << site.active << " ";
    }

    return out;
  }

};

template<typename Engine>
class VacancyWriter {
  const Lattice<Engine>& L;
public:
  VacancyWriter(const Lattice<Engine>& L) : L(L) { }

  // Output pairs (v, n) indicating the vacancy id and the site on which it is found

  friend std::ostream& operator << (std::ostream& out, const VacancyWriter& VW) {

    const auto& sites = VW.L.sites; // Save typing

    for(unsigned n=0; n<sites.size(); ++n) {
      if(!sites[n].occupied) out << sites[n].id << " " << n << " ";
    }

    return out;
  }

};

template<typename Engine>
class ParticleWriter {
  const Lattice<Engine>& L;
public:
  ParticleWriter(const Lattice<Engine>& L) : L(L) { }

  // Output pairs (v, n) indicating the particle id and the site on which it is found

  friend std::ostream& operator << (std::ostream& out, const ParticleWriter& PW) {

    const auto& sites = PW.L.sites; // Save typing

    for(unsigned n=0; n<sites.size(); ++n) {
      if(sites[n].occupied) out << sites[n].id << " " << n << " ";
    }

    return out;
  }

};


template<typename Engine>
class TriangleParticleWriter {
  const Triangle_lattice<Engine>& L;
public:
  TriangleParticleWriter(const Triangle_lattice<Engine>& L, ofstream& outfile) : L(L) { }

  // Output pairs (v, n, p) indicating the particle id and the site on which it is found. p is the number of particles
  // present on site n
friend std::ostream& operator << (std::ostream& out, const TriangleParticleWriter& PW) {

    const auto& sites = PW.L.sites; // Save typing

    for(unsigned n=0; n<sites.size(); ++n) {
      for (unsigned i = 0; i < sites[n].occupied.size(); i++){
        if(sites[n].occupied[i]) out << sites[n].id[i] << " " << n << " " << sites[n].present;
      }
    }

    return out;
  }

};

template<typename Engine>
class HexagonalParticleWriter {
  const Hexagonal_lattice<Engine>& L;
  //Parameters P;
public:
  HexagonalParticleWriter(const Hexagonal_lattice<Engine>& L, ofstream& outfile) : L(L) { }

  // Output is (v, n, j) as above with j being the binary site argument
friend std::ostream& operator << (std::ostream& out, const HexagonalParticleWriter& HPW) {

    const auto& sites = HPW.L.sites; // Save typing

    for(unsigned n=0; n<sites.size(); ++n) {
      for (unsigned j = 0; j < 2; j++){
        if(sites[n].occupied[j]) out << sites[n].id[j] << " " << n << " " << j << " ";    
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

template<typename Engine>
class ClusterWriter {
  const Lattice<Engine>& L;
public:
  ClusterWriter(const Lattice<Engine>& L) : L(L) { }

  friend std::ostream& operator << (std::ostream& out, const ClusterWriter& SW) {
    for(const auto& v: SW.L.cluster_distributions()) out << v << " ";
    return out;
  }
};




int main(int argc, char* argv[]) {

  // Load up default parameters
  Parameters P;

  // Set up command-line overrides
  CLI::App app{"Condensed Run-and-Tumble model (crumble)"};

  app.add_option("-L,--lengths",    P.L,      "Lattice lengths");
  app.add_option("-N,--particles",  P.N,      "Number of particles");
  app.add_option("-t,--tumble",     P.alpha,  "Tumble rate");
  app.add_option("-n,--occupation",     P.n_max,  "Max occupation of a site");

  // Output parameters
  std::string output = "";
  std::string lattice_type = ""; 
  double burnin = 1000, until = 5000, every = 2.5;
  unsigned localaverage = 0;
  double localinterval = 10.0;

  app.add_option("-o, --output", output, "Output type: snapshots|particles|vacancies|clusters");
  app.add_option("-l, --lattice_type", lattice_type, "Lattice type: square|triangular(being worked on)|hexagonal(being worked on)");

  app.add_option("-b,--burnin",        burnin,        "Time to run before starting measurements");
  app.add_option("-u,--until",         until,         "Time to run for once measurements started");
  app.add_option("-e,--every",         every,         "Measurement interval");
  app.add_option("-a,--localaverage",  localaverage,  "Number of local averages (clusters only; 0=stationary state)");
  app.add_option("-i,--localinterval", localinterval, "Interval between local averages");

  CLI11_PARSE(app, argc, argv);

  if(output[0] == 'p') output = "particles";
  else if(output[0] == 'v') output = "vacancies";
  else if(output[0] == 'c') output = "clusters";
  else output = "snapshots";

  if(lattice_type[0] == 's') lattice_type = "square";
  else if(lattice_type[0] == 't') lattice_type = "triangular";
  else if(lattice_type[0] == 'h') lattice_type = "hexagonal";

  if(localaverage > 0) {
    if(output != "clusters") {
      std::cerr << "Can only perform local averages in clusters mode" << std::endl;
      return 1;
    }
    if((localaverage-1)*localinterval >= every) {
      std::cerr << "Local averaging window larger than measurement internal" << std::endl;
      return 1;
    }
  }


  // Depending on what lattice, the output may be different
  if (lattice_type == "square"){
    // Initialise a random number generator and set up the model
    std::mt19937 rng((std::random_device())());
    Lattice L(P, rng);

    // Snapshot sequence has a header that sets out the simulation parameters
    std::cout << "# L = [ ";
    for(const auto& L: P.L) std::cout << L << " ";
    std::cout << "]" << std::endl;
    std::cout << "# N = " << P.N << std::endl;
    std::cout << "# alpha = [ ";
    for(const auto& alpha: P.alpha) std::cout << alpha << " ";
    std::cout << "]" << std::endl;
    std::cout << "# output = " << output << std::endl;
    std::cout << "# initial = " << burnin << std::endl;
    std::cout << "# interval = " << every << std::endl;
    if(localaverage > 0) {
      std::cout << "# localaverage = " << localaverage << std::endl;
      std::cout << "# localinterval = " << localinterval << std::endl;
    }

    double t = 0;

    if(output == "clusters") {
      if(localaverage == 0) {
        // We sum the histograms over all measurements
        hist_t sumhist;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          hist_t hist = L.cluster_distributions();
          // Add hist to sumhist taking into account they may have different lengths
          if(hist.size() > sumhist.size()) {
            hist[std::slice(0,sumhist.size(),1)] += sumhist;
            sumhist = std::move(hist);
          } else {
            sumhist[std::slice(0,hist.size(),1)] += hist;
          }
        }
        for(const auto& k: sumhist) std::cout << k << " ";
        std::cout << std::endl;
      } else {
        // We perform local averages of the histograms at the given measurement interval
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          hist_t sumhist = L.cluster_distributions();
          for(unsigned m=1; m<localaverage; ++m) {
            t = L.run_until(burnin + n*every + m*localinterval);
            hist_t hist = L.cluster_distributions();
            // Add hist to sumhist taking into account they may have different lengths
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }
          }
          for(const auto& k: sumhist) std::cout << k << " ";
          std::cout << std::endl;
        }
      }
    } else {
      for(unsigned n=0; t < burnin + until; ++n) {
        t = L.run_until(burnin + n * every);
        if (output == "particles") std::cout << ParticleWriter(L) << std::endl;
        else if(output == "vacancies") std::cout << VacancyWriter(L) << std::endl;
        else {
          // Ensure that all the particle directions are at the simulation time
          L.realise_directions();
          std::cout << SnapshotWriter(L) << std::endl;
        }
      }

    }
  } else if(lattice_type == "triangular"){
    // Initialise a random number generator and set up the model
    std::mt19937 rng((std::random_device())());
    Triangle_lattice TL(P, rng);

    // Snapshot sequence has a header that sets out the simulation parameters
    std::cout << "# L = [ ";
    for(const auto& L: P.L) std::cout << L << " ";
    std::cout << "]" << std::endl;
    std::cout << "# N = " << P.N << std::endl;
    std::cout << "# alpha = [ ";
    for(const auto& alpha: P.alpha) std::cout << alpha << " ";
    std::cout << "]" << std::endl;
    std::cout << "# output = " << output << std::endl;
    std::cout << "# initial = " << burnin << std::endl;
    std::cout << "# interval = " << every << std::endl;
    if(localaverage > 0) {
      std::cout << "# localaverage = " << localaverage << std::endl;
      std::cout << "# localinterval = " << localinterval << std::endl;
    }

    double t = 0;

    if(output == "clusters") {
      if(localaverage == 0) {
        // We sum the histograms over all measurements
        hist_t sumhist;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          hist_t hist = TL.cluster_distributions();
          // Add hist to sumhist taking into account they may have different lengths
          if(hist.size() > sumhist.size()) {
            hist[std::slice(0,sumhist.size(),1)] += sumhist;
            sumhist = std::move(hist);
          } else {
            sumhist[std::slice(0,hist.size(),1)] += hist;
          }
        }
        for(const auto& k: sumhist) std::cout << k << " ";
        std::cout << std::endl;
      } else {
        // We perform local averages of the histograms at the given measurement interval
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          hist_t sumhist = TL.cluster_distributions();
          for(unsigned m=1; m<localaverage; ++m) {
            t = TL.run_until(burnin + n*every + m*localinterval);
            hist_t hist = TL.cluster_distributions();
            // Add hist to sumhist taking into account they may have different lengths
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }
          }
          for(const auto& k: sumhist) std::cout << k << " ";
          std::cout << std::endl;
        }
      }
    } else {
      ofstream outfile;
      outfile.open("triangle.txt");
      //outfile << "# L = [ ";
      //for(const auto& L: P.L) outfile << L << " ";
      //outfile << "]" << endl;
      //outfile << "# N = " << P.N << endl;
      //outfile << "# alpha = [ ";
      //for(const auto& alpha: P.alpha) outfile << alpha << " ";
      //outfile << "]" << endl;
      //outfile << "# output = " << output << endl;
      //outfile << "# initial = " << burnin << endl;
      //outfile << "# interval = " << every << endl;

      for(unsigned n=0; t < burnin + until; ++n) {
        t = TL.run_until(burnin + n * every);
        // only doing a positional output here
        
        
        outfile << TriangleParticleWriter(TL, outfile) << endl;

      }

    }


  } else if(lattice_type == "hexagonal"){
    // Initialise a random number generator and set up the model
    std::mt19937 rng((std::random_device())());
    Hexagonal_lattice HL(P, rng);

    // Snapshot sequence has a header that sets out the simulation parameters
    std::cout << "# L = [ ";
    for(const auto& L: P.L) std::cout << L << " ";
    std::cout << "]" << std::endl;
    std::cout << "# N = " << P.N << std::endl;
    std::cout << "# alpha = [ ";
    for(const auto& alpha: P.alpha) std::cout << alpha << " ";
    std::cout << "]" << std::endl;
    std::cout << "# output = " << output << std::endl;
    std::cout << "# initial = " << burnin << std::endl;
    std::cout << "# interval = " << every << std::endl;
    if(localaverage > 0) {
      std::cout << "# localaverage = " << localaverage << std::endl;
      std::cout << "# localinterval = " << localinterval << std::endl;
    }

    double t = 0;

    if(output == "clusters") {
      if(localaverage == 0) {
        // We sum the histograms over all measurements
        hist_t sumhist;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          hist_t hist = HL.cluster_distributions();
          // Add hist to sumhist taking into account they may have different lengths
          if(hist.size() > sumhist.size()) {
            hist[std::slice(0,sumhist.size(),1)] += sumhist;
            sumhist = std::move(hist);
          } else {
            sumhist[std::slice(0,hist.size(),1)] += hist;
          }
        }
        for(const auto& k: sumhist) std::cout << k << " ";
        std::cout << std::endl;
      } else {
        // We perform local averages of the histograms at the given measurement interval
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          hist_t sumhist = HL.cluster_distributions();
          for(unsigned m=1; m<localaverage; ++m) {
            t = HL.run_until(burnin + n*every + m*localinterval);
            hist_t hist = HL.cluster_distributions();
            // Add hist to sumhist taking into account they may have different lengths
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }
          }
          for(const auto& k: sumhist) std::cout << k << " ";
          std::cout << std::endl;
        }
      }
    } else {
      ofstream outfile;
      outfile.open("hexagonal.txt");
      //outfile << "# L = [ ";
      //for(const auto& L: P.L) outfile << L << " ";
      //outfile << "]" << endl;
      //outfile << "# N = " << P.N << endl;
      //outfile << "# alpha = [ ";
      //for(const auto& alpha: P.alpha) outfile << alpha << " ";
      //outfile << "]" << endl;
      //outfile << "# output = " << output << endl;
      //outfile << "# initial = " << burnin << endl;
      //outfile << "# interval = " << every << endl;

      for(unsigned n=0; t < burnin + until; ++n) {
        t = HL.run_until(burnin + n * every);
        // only doing a positional output here
        
        outfile << HexagonalParticleWriter(HL, outfile) << endl;

      }

    }



  }
  // maybe do other lattices as well ?
  return 0;
}
