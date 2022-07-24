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
template<typename Engine> class HexDirectionWriter;

// ! There will be some changes is this class, mostly the implementation of occupation number

template<typename Engine>
class Lattice {

  // For output in different formats
  friend class SnapshotWriter<Engine>;
  friend class VacancyWriter<Engine>;
  friend class ParticleWriter<Engine>;
  friend class ClusterWriter<Engine>;

  static constexpr unsigned n_max = 3; // ! Nicer way of doing this?
    
    

    // Data associated with each site; by default all of these are set to zero
    // which represents a vacant site
    struct Site {
        std::vector<unsigned> id = std::vector<unsigned>(n_max); // Particle / vacancy id
        std::vector<bool> occupied = std::vector<bool>(n_max); // There is a particle here
        std::vector<bool> active = std::vector<bool>(n_max); // A move event is scheduled
        direction_t neighbours; // Number of neighbours that are occupied
        std::vector<direction_t> direction = std::vector<direction_t>(n_max); // Direction of last hop attempt
        std::vector<double> hoptime = std::vector<double>(n_max); // Time of last hop attempt
        unsigned present = 0; // Number of particles present at site. Has to be <= n_max
        std::vector<double> last_jump = std::vector<double>(n_max); // time of last jump made
    };

  Parameters P; // A local copy of the model parameters
  std::vector<Site> sites; // Representation of the sites
  Engine& rng; // Source of noise: this is a reference as there should only be one of these!
  std::discrete_distribution<unsigned> anyway, initial; // Distribution over tumble directions
  std::exponential_distribution<double> run; // Distribution of times between run and tumble events
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
  void place(unsigned n, unsigned id, direction_t d, double t, unsigned index) {
    sites[n].id[index] = id;
    sites[n].direction[index] = d;
    sites[n].hoptime[index] = t;
    sites[n].last_jump[index] = t;
    if(!sites[n].occupied[index]) {
      sites[n].occupied[index] = true;
      for(const auto&m: neighbours(n)) ++sites[m].neighbours;
    }
    sites[n].present++;
  }

  // Schedule a hop event for a particle at site n
  void schedule(unsigned n, unsigned index) {
    assert(sites[n].occupied[index]);
    S.schedule(run(rng), [this,n, index]() {
      assert(sites[n].occupied[index]);
      assert(sites[n].active[index]);
      // If there are no local vacancies, mark this particle as inactive and exit
      if(sites[n].neighbours == 2*P.L.size() * P.n_max) {
        // std::cout << "Can't move from "; decode(n); std::cout << " deactivating" << std::endl;
        sites[n].active[index] = false;
      } else {
        // if(std::uniform_real_distribution<double>()(rng)>=std::exp(-P.alpha*(S.time()-sites[n].hoptime))) {
        if (P.alpha[0] - 1e-5 > 0){
            if (tumble(rng) < S.time() - sites[n].hoptime[index]) {
              sites[n].direction[index] = anyway(rng);
            }
        }
        sites[n].hoptime[index] = S.time();
        // Get the sites adjacent to the departure site
        auto dnbs = neighbours(n);
        if (sites[dnbs[sites[n].direction[index]]].present < P.n_max) {
          auto itr = std::find(sites[dnbs[sites[n].direction[index]]].occupied.begin(), sites[dnbs[sites[n].direction[index]]].occupied.end(), false);
          unsigned k = std::distance(sites[dnbs[sites[n].direction[index]]].occupied.begin(), itr);
          if (k < (P.n_max) && !sites[dnbs[sites[n].direction[index]]].occupied[k]){
                    //std::cout << "here" << endl;
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
            }
            });
        sites[n].active[index] = true;
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
  std::exponential_distribution<double> tumble;
  Lattice(const Parameters& P, Engine& rng) :
    P(P), // NB: this takes a copy
    sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
    rng(rng), // NB: this takes a reference
    run(1),
    tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0)/P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
    {
      // Set up the tumble direction distribution
      std::vector<double> drates(2*P.L.size());
      std::vector<double> drates_initial(4);
      for(unsigned d=0;d < P.L.size(); ++d) {
         drates[2*d] = drates[2*d+1] = d<P.alpha.size() ? P.alpha[d]/tumble.lambda() : 1.0;
         drates_initial[2*d] = drates_initial[2*d+1] = 1.0;
      }
      anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());
      initial = std::discrete_distribution<unsigned>(drates_initial.begin(), drates_initial.end());

      std::vector<unsigned>position;
        for(unsigned i = 0; i < P.n_max; i++){
          for(unsigned n = 0; n < sites.size(); n++){
            position.push_back(i*sites.size() + n);
          }
        }


        unsigned possibilities = position.size();
        unsigned unplaced = P.N;
        unsigned id = 0;
        while (unplaced > 0){
          unsigned index = std::uniform_int_distribution<unsigned>(0, possibilities)(rng);
          
          unsigned l = position[index];


          unsigned n = l%sites.size();
          unsigned i = l/sites.size();

          place(n, id, initial(rng), 0.0, i);

          position.erase(position.begin()+index);
          position.push_back(l);
          
          id++;
          possibilities--;
          unplaced--;
        }

        for (unsigned index = 0; index < P.n_max; index++){
          for (unsigned n = 0; n < sites.size(); ++n) {
            if (sites[n].occupied[index] == false){
              sites[n].id[index] = id;
              id++;
            }
          }
        }
      /*
      // Place particles on the lattice
      unsigned unplaced = P.N; // Current number of particles remaining to be placed
      unsigned id_vac = P.N;
      for (unsigned index = 0; index < P.n_max; index++){
        for (unsigned n = 0; n < sites.size(); ++n) {
            // Number of remaining sites where partcles could be placed is sites.size()-n, unplaced of which need to be filled
            if (std::uniform_int_distribution<unsigned>(1, sites.size() - n)(rng) <= unplaced) {
                // For ease we only place one particle per site in the initial configuration
                place(n, P.N - unplaced, anyway(rng), 0.0, index);
                --unplaced;
            }
            // vacancies
            else {
              sites[n].id[index] = id_vac;
              id_vac++;
            }
        }
      }
      assert(unplaced == 0);
      */
        // Activate particles that can move, and schedule a hop accordingly
      for (unsigned n = 0; n < sites.size(); ++n) {
        for (unsigned k = 0; k < P.n_max; ++k){
        if (sites[n].occupied[k] && sites[n].neighbours < 4 * P.n_max) schedule(n, k);
        }
      }
      assert(consistent());


    }

  // Iterates the simulation until the given time; returns the actual time run to
  double run_until(double time) {
    while(S.advance(time));
    assert(consistent());
    return S.time();
  }

  void set_new_lambda(std::exponential_distribution<double> *exp_dis, double val){
      typename std::exponential_distribution<double>::param_type new_lambda(val);
      exp_dis->param(new_lambda);
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

    double motility_fraction(){
      double count = 0;
      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned k = 0; k < P.n_max; k++){
            if (!sites[n].occupied[k]) continue;
            if (tumble(rng) < S.time() - sites[n].hoptime[k]) {
                // At least one tumble event happened since we last attempted a hop; so randomise the direction
                sites[n].direction[k] = anyway(rng);
            }
            auto dnbs = neighbours(n);
            if (sites[dnbs[sites[n].direction[k]]].present == P.n_max){
              count++;
            }
        }
      }
      return count/double(P.N);
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present == 0) dists[2 * kv.second.size() - 1]++;
            else dists[2 * kv.second.size() - 2]++;
        }
        return dists;
    }

    map<unsigned, std::list<unsigned>> cluster() const {
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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

    hist_t cluster_distribution_particle_number() const {
        // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present == 0) dists_nr[2 * kv.second.size() - 1]++;
            else dists_nr[2 * kv.second.size() - 2]++;
        }
        return dists_nr;
    }

    
    // Function to determine the size of largest cluster. Note that we will only regard particle clusters here 
    // (at least so far)
    size_t max_cluster_size(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        std::size_t max_s = 1;
        for (const auto& kv : clusters) {

            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) {
              //std::cout << "-------------------" << endl;
              //std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
              max_s = std::max(max_s, kv.second.size()); 
              //std::cout << max_s << endl;
            
            }
            
        }

        //std::cout << max_s << endl;

        return max_s;
    }

    size_t max_cluster_size_nr(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {

            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) {
              //std::cout << "-------------------" << endl;
              //std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
              max_s = std::max(max_s, kv.second.size()); 
              //std::cout << max_s << endl;
            
            }
            
        }

        //std::cout << max_s << endl;

        return max_s;
    }

    // Function to determine the mean cluster size by simple taking the sum off all clusters (not vacant ones) and dividing by number of clusters. Note that we will only regard particle clusters here 
    // (at least so far)
    double avg_cluster_size(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // variable for mean and the count of clusters
        double mean = 0; 
        double count = 0; 
        for (const auto& kv : clusters) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) { 
              mean += kv.second.size();
              count += 1;
            }
        }

        mean = mean / count;
        return mean;
    }

    double avg_cluster_size_nr(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) { 
              mean += kv.second.size();
              count += 1;
            }
        }

        mean = mean / count;

        return mean;
    }

    // function to return the number of particle clusters
    unsigned number_cluster(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // variable for mean and the count of clusters
        unsigned count = 0; 
        for (const auto& kv : clusters) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) { 
              count += 1;
            }
        }

        return count;
    }


    vec surface_volume(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        // the index 2i is for volume and 2i+1 is for surface
        vec output;
        
        for (const auto& kv : clusters) {
            // variable for surface and check 
            unsigned surf = 0;
            unsigned check = 0;
            // check whether particle cluster
            if (sites[kv.first].present > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present == 0){
                    check++;
                    continue;
                  }
                }
                if (check > 0) surf++;
              }
              //std::cout << surf << " " << kv.second.size() << endl;
              output.push_back(kv.second.size());
              output.push_back(surf);
            }
        }

        return output;
    }

    vec surface_volume_nr(){
        std::map<unsigned, std::list<unsigned>> clusters = cluster();

        // vector for output of surface and volume of clusters
        // the index 2i is for volume and 2i+1 is for surface
        vec output;
        
        for (const auto& kv : clusters) {
            // variable for surface and check 
            unsigned surf = 0;
            unsigned check = 0;
            // check whether particle cluster
            if (sites[kv.first].present > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present == 0){
                    check++;
                    continue;
                  }
                }
                if (check > 0) surf += sites[n].present;
              }
              //std::cout << surf << " " << kv.second.size() << endl;
              output.push_back(kv.second.size());
              output.push_back(surf);
            }
        }

        return output;
    }
    vec cluster_surface(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        vec output;
        
        for (unsigned n = 0; n < sites.size(); n++){
          unsigned check =0;
          if (sites[n].present > 0){
            for (const auto& m: neighbours(n)){
              if (sites[m].present < 1){
                    check++;
                    continue;
                  }
            }
            if (check > 0) output.push_back(n);
          }
        }

        /*
        for (const auto& kv : clusters) {
            // variable for surface and check 
            unsigned check = 0;
            // check whether particle cluster
            if (sites[kv.first].present > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present < 1){
                    check++;
                    continue;
                  }
                }
                if (check > 0) output.push_back(n);
                
              }
            }
        }
        */

        return output;
    }
    // function to give out position of a cluster in the lattice that has size between lower and upper bound
    vec single_cluster(unsigned clust_min, unsigned clust_max){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        vec output;
        
        for (const auto& kv : clusters) {
            // check whether particle cluster
            if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max) { 
              // to get surface and volume 
              unsigned surf = 0;
              unsigned check = 0;
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present == 0){
                    check++;
                    continue;
                  }
                }
                if (check > 0) surf++;
              }
              // add surface and volume to beginning of output
              output.push_back(kv.second.size());
              output.push_back(surf);
              // add position of particles
              for (const auto& n : kv.second){
                output.push_back(n);
              }
              break;
            }
        }

        return output;
    }

    unsigned clust_size(unsigned clust_min, unsigned clust_max){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        unsigned output = 0;
        
        for (const auto& kv : clusters) {
            // check whether particle cluster
            if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max) { 
              output = 1;
            }
        }

        return output;
    }

    hist_t particle_neighbour_dist(){
      hist_t dist(5);
      unsigned count = 0;
      for (unsigned n = 0; n < sites.size(); n++){
        count = 0;
        if (sites[n].present > 0){
        for (const auto& m : neighbours(n)){
          if (sites[m].present > 0) count++;
        }
        dist[count]++;
        }
        
      }
      return dist;
    }

    vec_d density(){
      vec_d den;

      for (unsigned n = 0; n < sites.size(); n++){
          double local = 0.5 * double(sites[n].present);
          for (const auto& m : neighbours(n)){
            local += 1.0/8.0 * double(sites[m].present);
          }
          den.push_back(local);
        }
      
      return den;
    }
    
    vec_d stopping(double t){ 
      vec_d stopping_time;

      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned i = 0; i < P.n_max; i++){
          if (sites[n].occupied[i] == true) stopping_time.push_back(t - sites[n].last_jump[i]);
        }
        
      }
      return stopping_time;
    }

    // function for returning perimeter order parameter (see thesis)
    double perimeter(){
      // get distribution and surface-volume 
      hist_t dist = cluster_distribution_particle_number();
      vec sv = surface_volume_nr();


      double second_moment = 0;
      double first_moment = 0;

      // calculate appropriate stuff
      for (unsigned i = 0; i < dist.size(); i+=2){
        for (unsigned k = 0; k < sv.size(); k+=2){
          if (sv[k] == ((i+2)/2)){
            second_moment += ((i+2)/2) * sv[k+1];
          }
        }
        first_moment += dist[i] * ((i+2)/2);
        
      }
      return second_moment/first_moment;
    }
};


// ! Following classes are additions by Lars. Note that in the triangular case it is just the one as above with some very slight adjustments
// ! while the one for the hexagonal lattice has more drastic changes. 


template<typename Engine>
class Triangle_lattice {

    // For output in different formats
    friend class SnapshotWriter<Engine>;
    friend class VacancyWriter<Engine>;
    friend class TriangleParticleWriter<Engine>;
    friend class ClusterWriter<Engine>;

    Parameters P; // A local copy of the model parameters

    static constexpr unsigned n_max = 5; // ! Nicer way of doing this?
    
    

    // Data associated with each site; by default all of these are set to zero
    // which represents a vacant site
    struct Site {
        std::vector<unsigned> id = std::vector<unsigned>(n_max); // Particle / vacancy id
        std::vector<bool> occupied = std::vector<bool>(n_max); // There is a particle here
        std::vector<bool> active = std::vector<bool>(n_max); // A move event is scheduled
        direction_t neighbours; // Number of neighbours that are occupied
        std::vector<direction_t> direction = std::vector<direction_t>(n_max); // Direction of last hop attempt
        std::vector<double> hoptime = std::vector<double>(n_max); // Time of last hop attempt
        unsigned present = 0; // Number of particles present at site. Has to be <= n_max
        std::vector<double> last_jump = std::vector<double>(n_max); // time of last jump made
    };

    std::vector<Site> sites; // Representation of the sites
    Engine& rng; // Source of noise: this is a reference as there should only be one of these!
    std::discrete_distribution<unsigned> anyway, initial; // Distribution over tumble directions
    std::exponential_distribution<double> run; // Distribution of times between run and tumble events
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
        xm = (m+L)%L;

        if (((k-1)/L - L) >= 0){ 
          yk = 0;
        }else{
          yk = ((k-1)/L);
        }
        if ((m+1) < 0){ 
          ym = L - 1;
        }else{
          int l = std::abs(int(n - L));
          //if (l < 0) l = -l;
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
        sites[n].last_jump[index] = t;
        if (!sites[n].occupied[index]) {
            
            sites[n].occupied[index] = true;
            //std::cout << "n " << n << endl;
            for (const auto& m : neighbours(n)){
              //std::cout << "m " << m << endl;
             ++sites[m].neighbours;
            }
            //std::cout << "---------------------------------" << endl;
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
                
                auto dnbs = neighbours(n);
                if (sites[dnbs[sites[n].direction[index]]].present < P.n_max) {
                  auto itr = std::find(sites[dnbs[sites[n].direction[index]]].occupied.begin(), sites[dnbs[sites[n].direction[index]]].occupied.end(), false);
                  unsigned k = std::distance(sites[dnbs[sites[n].direction[index]]].occupied.begin(), itr);
                  if (k < (P.n_max) && !sites[dnbs[sites[n].direction[index]]].occupied[k]){
                    
                    //std::cout << "here" << endl;
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
    std::exponential_distribution<double> tumble;
    Triangle_lattice(const Parameters& P, Engine& rng) :
        P(P), // NB: this takes a copy
        sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
        rng(rng), // NB: this takes a reference
        run(1), // On average, it runs once per time unit
        tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
    {
        // Set up the tumble direction distribution
        std::vector<double> drates(2 * P.L.size() + 2);
        std::vector<double> drates_initial(6);
        for (unsigned d = 0; d < 3; ++d) {
            drates[2 * d] = drates[2 * d + 1] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
            drates_initial[2*d] = drates_initial[2*d+1] = 1.0;
        }
        anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end());
        initial = std::discrete_distribution<unsigned>(drates_initial.begin(), drates_initial.end());
        
        std::vector<unsigned>position;
        for(unsigned i = 0; i < P.n_max; i++){
          for(unsigned n = 0; n < sites.size(); n++){
            position.push_back(i*sites.size() + n);
          }
        }


        unsigned possibilities = position.size();
        unsigned unplaced = P.N;
        unsigned id = 0;
        while (unplaced > 0){
          unsigned index = std::uniform_int_distribution<unsigned>(0, possibilities)(rng);
          
          unsigned l = position[index];


          unsigned n = l%sites.size();
          unsigned i = l/sites.size();

          place(n, id, initial(rng), 0.0, i);

          position.erase(position.begin()+index);
          position.push_back(l);
          
          id++;
          possibilities--;
          unplaced--;
        }

        for (unsigned index = 0; index < P.n_max; index++){
          for (unsigned n = 0; n < sites.size(); ++n) {
            if (sites[n].occupied[index] == false){
              sites[n].id[index] = id;
              id++;
            }
          }
        }

        /*
        // Place particles on the lattice
        unsigned unplaced = P.N; // Current number of particles remaining to be placed
        unsigned id_vac = P.N;
        for (unsigned index = 0; index < P.n_max; index++){
          for (unsigned n = 0; n < sites.size(); ++n) {
            // Number of remaining sites where partcles could be placed is sites.size()-n, unplaced of which need to be filled
            if (std::uniform_int_distribution<unsigned>(1, sites.size() - n)(rng) <= unplaced) {
                // For ease we only place one particle per site in the initial configuration
                place(n, P.N - unplaced, anyway(rng), 0.0, index);
                --unplaced;
            }
            // vacancies
            else {
              sites[n].id[index] = id_vac;
              id_vac++;
            }
          }
        }
        assert(unplaced == 0);
        */
        // Activate particles that can move, and schedule a hop accordingly
        for (unsigned n = 0; n < sites.size(); ++n) {
          for (unsigned k = 0; k < P.n_max; ++k){
            if (sites[n].occupied[k] && sites[n].neighbours < 6 * P.n_max) schedule(n, k);
          }
        }
        assert(consistent());

    }

    void set_new_lambda(std::exponential_distribution<double> *exp_dis, double val){
      typename std::exponential_distribution<double>::param_type new_lambda(val);
      exp_dis->param(new_lambda);
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

    double motility_fraction(){
      double count = 0;
      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned k = 0; k < P.n_max; k++){
            if (!sites[n].occupied[k]) continue;
            if (tumble(rng) < S.time() - sites[n].hoptime[k]) {
                // At least one tumble event happened since we last attempted a hop; so randomise the direction
                sites[n].direction[k] = anyway(rng);
            }
            auto dnbs = neighbours(n);
            if (sites[dnbs[sites[n].direction[k]]].present == P.n_max){
              count++;
            }
        }
      }
      return count/double(P.N);
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present == 0) dists[2 * kv.second.size() - 1]++;
            else dists[2 * kv.second.size() - 2]++;
        }
        return dists;
    }


    hist_t cluster_distributions_particle_numbers() const {
        // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present == 0) dists_nr[2 * kv.second.size() - 1]++;
            else dists_nr[2 * kv.second.size() - 2]++;
        }
        return dists_nr;
    }

    // Function to determine the size of largest cluster. Note that we will only regard particle clusters here 
    // (at least so far)
    size_t max_cluster_size(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        std::size_t max_s = 1;
        for (const auto& kv : clusters) {

            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) {
              //std::cout << "-------------------" << endl;
              //std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
              max_s = std::max(max_s, kv.second.size()); 
              //std::cout << max_s << endl;
            
            }
            
        }

        //std::cout << max_s << endl;

        return max_s;
    }

    size_t max_cluster_size_nr(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {

            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) {
              //std::cout << "-------------------" << endl;
              //std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
              max_s = std::max(max_s, kv.second.size()); 
              //std::cout << max_s << endl;
            
            }
            
        }

        //std::cout << max_s << endl;

        return max_s;
    }

    // Function to determine the mean cluster size by simple taking the sum off all clusters (not vacant ones) and dividing by number of clusters. Note that we will only regard particle clusters here 
    // (at least so far)
    double avg_cluster_size(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // variable for mean and the count of clusters
        double mean = 0; 
        double count = 0; 
        for (const auto& kv : clusters) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) { 
              mean += kv.second.size();
              count += 1;
            }
        }

        mean = mean / count;

        return mean;
    }

    double avg_cluster_size_nr(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) { 
              mean += kv.second.size();
              count += 1;
            }
        }

        mean = mean / count;

        return mean;
    }

    // function to return the number of particle clusters
    unsigned number_cluster(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          if (sites[n].present == 0) clusters_nr[n] = std::list<unsigned>(1, n);
          else clusters_nr[n] = std::list<unsigned>(sites[n].present, n);

        }
        // Keep track of the size of the largest cluster
        std::size_t maxsize_nr = 1;

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
            for (const auto& m : forward_neighbours(n)) {
                unsigned large = memberof_nr[n], small = memberof_nr[m];
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
                    
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[kv.first].present > 0) { 
              count += 1;
            }
        }

        return count;
    }

    vec surface_volume(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        vec output;
        
        for (const auto& kv : clusters) {
            // variable for surface and check 
            unsigned surf = 0;
            unsigned check = 0;
            // check whether particle cluster
            if (sites[kv.first].present > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present == 0){
                    check++;
                    continue;
                  }
                }
                if (check > 0) surf++;
              }

              output.push_back(kv.second.size());
              output.push_back(surf);
            }
        }

        return output;
    }

    vec cluster_surface(){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        vec output;
        
        for (const auto& kv : clusters) {
            // variable for surface and check 
            unsigned check = 0;
            // check whether particle cluster
            if (sites[kv.first].present > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present == 0){
                    check++;
                    continue;
                  }
                }
                if (check > 0) output.push_back(n);
              }
            }
        }

        return output;
    }
    // function to give out position of a cluster in the lattice that has size between lower and upper bound
    vec single_cluster(unsigned clust_min, unsigned clust_max){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        vec output;
        
        for (const auto& kv : clusters) {
            // check whether particle cluster
            if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max) { 
              // to get surface and volume 
              unsigned surf = 0;
              unsigned check = 0;
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m : neighbours(n)){
                  if (sites[m].present == 0){
                    check++;
                    continue;
                  }
                }
                if (check > 0) surf++;
              }
              // add surface and volume to beginning of output
              output.push_back(kv.second.size());
              output.push_back(surf);
              // add position of particles
              for (const auto& n : kv.second){
                output.push_back(n);
              }
              break;
            }
        }

        return output;
    }

    unsigned clust_size(unsigned clust_min, unsigned clust_max){
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
                // continue on if they are part of the same cluster
                if (small == large) continue;
                // continue on if they are vacant - not vacant and vise versa
                else if (sites[n].present == 0 && sites[m].present != 0) continue;
                else if (sites[n].present != 0 && sites[m].present == 0) continue;
                else {
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

        // vector for output of surface and volume of clusters
        unsigned output = 0;
        
        for (const auto& kv : clusters) {
            // check whether particle cluster
            if (sites[kv.first].present > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max) { 
              output = 1;
            }
        }

        return output;
    }

    hist_t particle_neighbour_dist(){
      hist_t dist(7);
      unsigned count = 0;
      
      for (unsigned n = 0; n < sites.size(); n++){
        count = 0;
        if (sites[n].present > 0){
        for (const auto& m : neighbours(n)){
          if (sites[m].present > 0) count++;
        }
        dist[count]++;
        }
        
      }

      return dist;
    }

    vec_d density(){
      vec_d den;
      

      for (unsigned n = 0; n < sites.size(); n++){
          double local = 0.5 * double(sites[n].present);
          for (const auto& m : neighbours(n)){
            local += 1.0/12.0 * double(sites[m].present);
          }
          den.push_back(local);
        }
      
      
      return den;
    }

    vec_d stopping(double t){ 
      vec_d stopping_time;

      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned i = 0; i < P.n_max; i++){
          if (sites[n].occupied[i] == true) stopping_time.push_back(t - sites[n].last_jump[i]);
        }
        
      }
      return stopping_time;
    }
};


template<typename Engine>
class Hexagonal_lattice {

    // For output in different formats
    friend class SnapshotWriter<Engine>;
    friend class VacancyWriter<Engine>;
    friend class HexagonalParticleWriter<Engine>;
    friend class ClusterWriter<Engine>;
    friend class HexDirectionWriter<Engine>;

    Parameters P; // A local copy of the model parameters
    // TODO Really need a better way of doing this
    static constexpr unsigned n_max = 3;

    // Data associated with each site; by default all of these are set to zero
    // which represents a vacant site

    // We will use the following way to identify particles and their respective positions
    // If index%2 == 0 it is the upper position, if it is index%2 == 1 it is the lower position
    struct Site {
        std::vector<unsigned> id  = std::vector<unsigned>(2 * n_max); // Particle / vacancy id
        std::vector<bool> occupied = std::vector<bool>(2 * n_max); // There is a particle here
        std::vector<bool> active = std::vector<bool>(2 * n_max); // A move event is scheduled
        std::vector<direction_t> neighbours = std::vector<direction_t>(2); // Number of neighbours that are occupied
        std::vector<direction_t> direction = std::vector<direction_t>(2 * n_max); // Preferred direction of movement
        std::vector<double> hoptime = std::vector<double>(2 * n_max); // Time of last hop attempt
        std::vector<int> present = std::vector<int>(2); // Number of particles present at each site in unit cell
        std::vector<unsigned> current_dir = std::vector<unsigned>(2*n_max); // Contains the lattice site it points at.
        std::vector<double> last_jump = std::vector<double>(2*n_max); // time of last jump made
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

    std::vector<Site> sites; // Representation of the sites
    Engine& rng; // Source of noise: this is a reference as there should only be one of these!
    std::discrete_distribution<unsigned> anyway, initial; // Distribution over tumble directions
    std::exponential_distribution<double> run; // Distribution of times between run and tumble events
    Scheduler S; // Keeps track of the event queue
    // Given an index into sites, return a sequence of indices corresponding to
    // its neighbours. We have periodic boundary conditions


    // function to choose direction based on what site and which directional orientation it has
    auto preference_direction(unsigned n, unsigned index, direction_t direction){
      int dir = 0;
      
      if (direction == 0){
          if (index%2==0) dir = 5;
          else dir = 2;
      }
      else if (direction == 1){
          if (index%2==0) dir = 5;
          else dir = 0;
      }
      else if (direction == 2){
          if (index%2==0) dir = 3;
          else dir = 0;
      }
      else if (direction == 3){
          if (index%2==0) dir = 3;
          else dir = 4;
      }
      else if (direction == 4){
          if (index%2==0) dir = 1;
          else dir = 4;
      }
      else if (direction == 5){
          if (index%2==0) dir = 1;
          else dir = 2;
      }
      return dir;
    }


    auto neighbours_dir(unsigned n) const {

        // this is only valid for 2d

        std::vector<unsigned> nbs(6);
        int L1 = P.L[0];
        int L2 = P.L[1];

        unsigned x = n % L1;
        unsigned y = (n/L1);


        nbs[0] = n;                               // same site
        nbs[1] = n;                               // same site
        nbs[2] = (n + 1) % L1 + y * L1;           // right
        nbs[3] = (n - 1) % L1 + y * L1;           // left
        nbs[4] = x + ((y - 1 + L1) % L2) * L1;    // down
        nbs[5] = x + ((y + 1) % L2) * L1;         // up
        return nbs;
    }

    auto neighbours(unsigned n, unsigned index) const {

        // this is only valid for 2d

        std::vector<unsigned> nbs(3);
        int L1 = P.L[0];
        int L2 = P.L[1];

        unsigned x = n % L1;
        unsigned y = (n/L1);

        if(index){
          nbs[0] = n;                               // same site
          nbs[1] = (n + 1) % L1 + y * L1;           // right
          nbs[2] = x + ((y - 1 + L1) % L2) * L1;    // down
        } else {
          nbs[0] = n;                               // same site
          nbs[1] = (n - 1) % L1 + y * L1;           // left
          nbs[2] = x + ((y + 1) % L2) * L1;         // up
        }
        return nbs;
    }


    auto forward_neighbours(unsigned n, unsigned index) const {

        // this is only valid for 2d
        std::vector<unsigned> nbs(2);
        
        int L1 = P.L[0];
        int L2 = P.L[1];

        unsigned x = n % L1;
        unsigned y = (n/L1);

        if (index){
          nbs[0] = n;                         // same site
          nbs[0] = n;                         // TODO this is not ideal, better way?                       

        } else{
          nbs[1] = (n - 1) % L1 + y * L1;     // left
          nbs[2] = x + ((y + 1) % L2) * L1;   // up
          

        }

        return nbs;
    }



    

    // Place a particle with given direction and hop time at site n;
    // neighbouring sites will be accordingly adjusted
    void place(unsigned n, unsigned id, direction_t d, double t, unsigned index, unsigned dir) {
        sites[n].id[index] = id;
        sites[n].direction[index] = d;
        sites[n].hoptime[index] = t;
        sites[n].last_jump[index] = t;
        if (!sites[n].occupied[index]) {
            sites[n].occupied[index] = true;
            for (const auto& m : neighbours(n, index%2)) ++sites[m].neighbours[(index+1)%2];
        }
        sites[n].present[index%2]++;
        sites[n].current_dir[index] = dir;
    }


    // Schedule a hop event for a particle at site n
    void schedule(unsigned n, unsigned index) {
        assert(sites[n].occupied[index]);
        planned_moves++;
        S.schedule(run(rng), [this, n, index]() {
            assert(sites[n].occupied[index]);
            assert(sites[n].active[index]);
            // If there are no local vacancies, mark this particle as inactive and exit
            if (sites[n].neighbours[index] == 3 * P.n_max) {
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
                  //std::cout << n << " " << int(sites[n].direction[index]) << endl;
                  assert(!sites[dnbs[dir]].active[ind]);
                  // Get the id of the vacancy that is being displaced
                  unsigned vid = sites[dnbs[dir]].id[ind];
                  // std::cout << "Moving from "; decode(n); std::cout << " deactivating" << std::endl;
                  // Deactive the departure site; also mark it empty
                  sites[n].occupied[index] = sites[n].active[index] = false;
                  //std::cout << "before:\t n\t" << sites[n].present[index] << "\t m\t" << sites[dnbs[sites[n].direction[index]]].present[ind] << endl;
                  // Place a particle on the target site; it has the same direction and hoptime as the departing particle
                  // std::cout << "Moving to "; decode(dnbs[sites[n].direction]); std::cout << " placing" << std::endl;
                  int dir_next = preference_direction(dnbs[dir], ind, sites[n].direction[index]);
                  auto dnbs_next = neighbours_dir(dnbs[dir]);
                  place(dnbs[dir], sites[n].id[index], sites[n].direction[index], sites[n].hoptime[index], ind, dnbs_next[dir_next]);
                  // Move the vacancy id onto the departure site
                  sites[n].id[index] = vid;
                  sites[n].present[index%2] -= 1;
                  //if (sites[n].present[index%2] == -1) std::cout << "n: " << n << " j: " << index << endl;
                  //std::cout << "after:\t n\t" << sites[n].present[index] << "\t m\t" << sites[dnbs[sites[n].direction[index]]].present[ind] << endl;

                  // Now go through the neighbours of the departure site, update neighbour count and activate any
                  // that can now move. Note the particle this is at the target site is included in this list
                  // and will be activated accordingly
                  for (const auto& m : dnbs) {
                      --sites[m].neighbours[(index+1)%2];
                      for (unsigned k = 1 - index%2; k < 2 * P.n_max; k+=2)
                      if (sites[m].occupied[k] && !sites[m].active[k]) schedule(m, k);
                  }

                }
                /*
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
                */
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
              occupied += sites[n].present[i];
              unsigned nbs = 0;
              for (const auto& m : neighbours(n, i)) {
                  for (unsigned k = 1 - i%2; k < 2 * P.n_max; k+=2){
                    if (sites[m].occupied[(i+1)%2]) ++nbs;
                  }
              }
              if (nbs != sites[n].neighbours[i]) return false;
              // Check that mobile particles are active
              if (nbs < 3 * P.n_max && !sites[n].active[i]) return false;
              for (unsigned k = i; k < 2 * P.n_max; k+=2){
                if (sites[n].active[k]) ++active;
              }
            }
        }
        // Check we've not lost any particles
        return occupied == P.N && active == S.pending();
    }


public:
    unsigned planned_moves;
    unsigned actual_moves;
    std::exponential_distribution<double> tumble;

    Hexagonal_lattice(const Parameters& P, Engine& rng) :
        P(P), // NB: this takes a copy
        sites(std::accumulate(P.L.begin(), P.L.end(), 1, std::multiplies<unsigned>())), // Initialise lattice with empty sites
        rng(rng), // NB: this takes a reference
        run(1),
        tumble(std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size()) // Tumble time generator: set to the average of the given tumble rates
    {
        // Set up the tumble direction distribution
        std::vector<double> drates(6);
        std::vector<double> drates_initial(6);
        for (unsigned d = 0; d < 3; ++d) {
            drates[2*d] = drates[2*d + 1] = d < P.alpha.size() ? P.alpha[d] / tumble.lambda() : 1.0;
            drates_initial[2*d] = drates_initial[2*d + 1] = 1.0;
        }
        anyway = std::discrete_distribution<unsigned>(drates.begin(), drates.end()); 
        initial = std::discrete_distribution<unsigned>(drates_initial.begin(), drates_initial.end());
        

        std::vector<unsigned>position;
        for(unsigned i = 0; i < 2*P.n_max; i++){
          for(unsigned n = 0; n < sites.size(); n++){
            position.push_back(i*sites.size() + n);
          }
        }

        
        planned_moves = 0;
        actual_moves = 0;

        

        unsigned possibilities = position.size();
        unsigned unplaced = P.N;
        unsigned id = 0;
        while (unplaced > 0){
          unsigned index = std::uniform_int_distribution<unsigned>(0, possibilities)(rng);
          
          unsigned l = position[index];


          unsigned n = l%sites.size();
          unsigned i = l/sites.size();

          direction_t direction = initial(rng);
          int dir = preference_direction(n, i, direction);
          auto dnbs = neighbours_dir(n);

          place(n, id, direction, 0.0, i, dnbs[dir]);

          position.erase(position.begin()+index);
          position.push_back(l);
          
          id++;
          possibilities--;
          unplaced--;
        }

        for (unsigned index = 0; index < P.n_max; index++){
          for (unsigned n = 0; n < sites.size(); ++n) {
            if (sites[n].occupied[index] == false){
              sites[n].id[index] = id;
              id++;
            }
          }
        }
        /*
        // Place particles on the lattice
        unsigned unplaced = P.N; // Current number of particles remaining to be placed
        unsigned id_vac = P.N;
        // the outer loop is to allow for overcrowding of the cells. We have to sites per lattice site, and on each
        // of these we allow n_max particles. We therefore have to loop over all these as we can place particles 
        // on all of them
        for (unsigned index = 0; index < 2 * P.n_max; index++){
          for (unsigned n = 0; n < sites.size(); ++n) {
            // Number of remaining sites where partcles could be placed is sites.size()-n, unplaced of which need to be filled
            if (std::uniform_int_distribution<unsigned>(1, sites.size() - n)(rng) <= unplaced) {
                direction_t direction = anyway(rng);
                int dir = preference_direction(n, index, direction);
                auto dnbs = neighbours_dir(n);
                place(n, P.N - unplaced, direction, 0.0, index, dnbs[dir]);
                --unplaced;        
            }
            else {
                // don't need to randomize these, since they are vacancies
                sites[n].id[index] = id_vac;
                id_vac++;
            }
          }
        }
        assert(unplaced == 0);
        */
        // Activate particles that can move, and schedule a hop accordingly
        for (unsigned n = 0; n < sites.size(); ++n) {
          for (unsigned i = 0; i < 2 * P.n_max; i++){
            if (sites[n].occupied[i] && sites[n].neighbours[i%2] < 3 * P.n_max) schedule(n, i);
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

    void set_new_lambda(std::exponential_distribution<double> *exp_dis, double val){
      typename std::exponential_distribution<double>::param_type new_lambda(val);
      exp_dis->param(new_lambda);
    }
    // Sample all particle directions from the distribution that applies at the current instant
    // (This is needed if you want "realistic" snapshots, rather than the state at a mixture of times)
    void realise_directions() {
        for (auto& site : sites) {
          for (unsigned i = 0; i < 2 * P.n_max; i++){
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

    // function that calculates the fraction of immobile particles
    double motility_fraction(){
      double count = 0;
      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned i = 0; i < 2 * P.n_max; i++){
          if (!sites[n].occupied[i]) continue;
          if (tumble(rng) < S.time() - sites[n].hoptime[i]){
            sites[n].direction[i] = anyway(rng);
          }
          sites[n].hoptime[i] = S.time();
          // Get the sites adjacent to the departure site
          int dir = preference_direction(n, i, sites[n].direction[i]);
          auto dnbs = neighbours_dir(n);
          if (sites[dnbs[dir]].present[(i + 1)%2] == P.n_max){
            count++;
          }
        }
      }
      return count/double(P.N);
    }

    // Return cluster size distributions.
    // This will only calculate cluster size by area, not by particles contained
    // TODO use if present == 0 to distinguish between area sites
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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
            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] == 0) dists[2 * kv.second.size() - 1]++;
            else dists[2 * kv.second.size() - 2]++;
        }
        return dists;
    }

    // Maybe find a better name, this is very long....
    // Important note here: I am considering a site with at least one particle present as a
    // particle cluster. This means that the overall number of id's is not conserved.
    // For this, consider n particles and n_max=2. If all n are on different site, we lose n vacant id's
    // If, on the other hand, the n are on n/2 sites we lose 0 vacant id's
    hist_t cluster_distribution_particle_number() const {
        // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(2 * sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          for (unsigned i = 0; i < 2; i++){
            if (sites[n].present[i] == 0) clusters_nr[2*n + i] = std::list<unsigned>(1, 2*n + i);
            else {
              clusters_nr[2*n + i] = std::list<unsigned>(sites[n].present[i], 2*n+i);
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

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
          for (unsigned i = 0; i < 2; i++){
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof_nr[2*n + i], small = memberof_nr[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
          // instead of occupied we check whether there are any particles present
          if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] == 0) dists_nr[2 * kv.second.size() - 1]++;
          else dists_nr[2 * kv.second.size() - 2]++;
        }
        return dists_nr;
    }


    
    size_t max_cluster_size(){
        
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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
        std::size_t max_s = 1;
        for (const auto& kv : clusters) {

            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) {
              //std::cout << "-------------------" << endl;
              //std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
              max_s = std::max(max_s, kv.second.size()); 
              //std::cout << max_s << endl;
            
            }
            
        }
        //std::cout << max_s << endl;
        return max_s;
    }

    size_t max_cluster_size_nr(){
            // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(2 * sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          for (unsigned i = 0; i < 2; i++){
            if (sites[n].present[i] == 0) clusters_nr[2*n + i] = std::list<unsigned>(1, 2*n + i);
            else {
              clusters_nr[2*n + i] = std::list<unsigned>(sites[n].present[i], 2*n+i);
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

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
          for (unsigned i = 0; i < 2; i++){
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof_nr[2*n + i], small = memberof_nr[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {

            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) {
              //std::cout << "-------------------" << endl;
              //std::cout << kv.second.size() << " " << sites[kv.first].present << endl;
              max_s = std::max(max_s, kv.second.size()); 
              //std::cout << max_s << endl;
            
            }
            
        }

        //std::cout << max_s << endl;

        return max_s;
    }

    // Function to determine the mean cluster size by simple taking the sum off all clusters (not vacant ones) and dividing by number of clusters. Note that we will only regard particle clusters here 
    // (at least so far)
    double avg_cluster_size(){
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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

        // variable for mean and the count of clusters
        double mean = 0; 
        double count = 0; 
        for (const auto& kv : clusters) {
            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) { 
              mean += kv.second.size();
              count += 1;
            }
        }

        mean = mean / count;

        return mean;
    }

    double avg_cluster_size_nr(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(2 * sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          for (unsigned i = 0; i < 2; i++){
            if (sites[n].present[i] == 0) clusters_nr[2*n + i] = std::list<unsigned>(1, 2*n + i);
            else {
              clusters_nr[2*n + i] = std::list<unsigned>(sites[n].present[i], 2*n+i);
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

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
          for (unsigned i = 0; i < 2; i++){
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof_nr[2*n + i], small = memberof_nr[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) { 
              mean += kv.second.size();
              count += 1;
            }
        }

        mean = mean / count;
        return mean;
    }

    // function to return the number of particle clusters
    unsigned number_cluster(){
      // Lookup table of cluster membership by lattice site
        std::vector<unsigned> memberof_nr(2 * sites.size());
        // Initially, this is just the site id as each site is its own cluster
        std::iota(memberof_nr.begin(), memberof_nr.end(), 0);
        // Create also a map of clusters each containing a list of its members
        std::map<unsigned, std::list<unsigned>> clusters_nr;
        for (unsigned n = 0; n < sites.size(); ++n){
          for (unsigned i = 0; i < 2; i++){
            if (sites[n].present[i] == 0) clusters_nr[2*n + i] = std::list<unsigned>(1, 2*n + i);
            else {
              clusters_nr[2*n + i] = std::list<unsigned>(sites[n].present[i], 2*n+i);
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

        for (unsigned n = 0; n < sites.size(); ++n) {
            // Loop over neigbours m in one direction only so we visit each bond once
          for (unsigned i = 0; i < 2; i++){
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof_nr[2*n + i], small = memberof_nr[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
                    // Ensure we have large and small the right way round (this makes the algorithm slightly more efficient)
                    if (clusters_nr[large].size() < clusters_nr[small].size()) std::swap(large, small);
                    // Update the cluster number for all sites in the smaller one
                    for (const auto& site : clusters_nr[small]) memberof_nr[site] = large;
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
        for (const auto& kv : clusters_nr) {
            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) { 
              count += 1;
            }
        }

        return count;
    }

    vec surface_volume(){
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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

        // variable for mean and the count of clusters
        vec output; 
        for (const auto& kv : clusters) {
            unsigned surf = 0;
            unsigned check = 0;
            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m: neighbours((n-n%2)/2, n%2)){
                  if (sites[m].present[(n+1)%2] == 0) check++;
                }
                if (check > 0) surf++;
              }
              output.push_back(kv.second.size());
              output.push_back(surf);
            }
        }

        return output;
    }

    vec cluster_surface(){
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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

        // variable for mean and the count of clusters
        vec output; 
        for (const auto& kv : clusters) {
            unsigned check = 0;
            // instead of occupied we check whether there are any particles present
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0) { 
              for (const auto& n : kv.second){
                check = 0;
                for (const auto& m: neighbours((n-n%2)/2, n%2)){
                  if (sites[m].present[(n+1)%2] == 0) check++;
                }
                if (check > 0){
                  output.push_back((n-n%2)/2);
                  output.push_back(n%2);
                }
              }
            }
        }

        return output;
    }
    // function to give out position of a cluster in the lattice that has size between lower and upper bound
    vec single_cluster(unsigned clust_min, unsigned clust_max){
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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


        // vector for output of surface and volume of clusters
        vec output;
        std::cout << "test" << endl;
        
        for (const auto& kv : clusters) {
            // check whether particle cluster
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max) { 
              // to get surface and volume 
              unsigned surf = 0;
              unsigned check = 0;
              for (const auto& n : kv.second){
                  if (sites[(n-n%2)/2].present[n%2] > 0){
                    check = 0;
                    for (const auto& m : neighbours((n-n%2)/2, n%2)){
                      if (sites[m].present[(n+1)%2] == 0){
                        check++;
                        continue;
                      }
                    }
                    if (check > 0) surf++;
                  }
              }
              std::cout << "break" << endl;
              // add surface and volume to beginning of output
              output.push_back(kv.second.size());
              output.push_back(surf);
              // add position of particles
              for (const auto& n : kv.second){
                  if (sites[(n-n%2)/2].present[n%2] > 0){ 
                    output.push_back((n-n%2)/2);
                    output.push_back(n%2);
                  }
              }
              break;
            }
        }
        
        return output;
    }

    unsigned clust_size(unsigned clust_min, unsigned clust_max){
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
            // TODO figure out why forward neighbour makes algorithm not work
            for (const auto& m : neighbours(n, i)) {
                unsigned large = memberof[2*n + i], small = memberof[2*m + (i+1)%2];
                // If they are in the same cluster we can move on
                if (small == large) continue;
                // If one of them is empty but the other isn't we move on
                else if (sites[n].present[i] == 0 && sites[m].present[(i+1)%2] != 0) continue;
                else if (sites[n].present[i] != 0 && sites[m].present[(i+1)%2] == 0) continue;
                // merge clusters
                else {
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


        // vector for output of surface and volume of clusters
        unsigned output = 0;
        
        for (const auto& kv : clusters) {
            // check whether particle cluster
            if (sites[(kv.first - kv.first%2)/2].present[kv.first%2] > 0 && kv.second.size() >= clust_min && kv.second.size() <= clust_max) { 
              output = 1;
            }
        }

        return output;
    }

    hist_t particle_neighbour_dist(){
      hist_t dist(4);
      unsigned count = 0;
      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned i = 0; i < 2; i++){
          count = 0;
          if (sites[n].present[i] > 0){
          for (const auto& m : neighbours(n, i)){
            if (sites[m].present[(i+1)%2] > 0) count++;
          }
          dist[count]++;
          }
        }
      }
      return dist;
    }

    vec_d density(){
      vec_d den;

      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned i = 0; i < 2; i++){
          double local = 0.5 * double(sites[n].present[i]);
          for (const auto& m : neighbours(n, i)){
            local += 1.0/6.0 * double(sites[m].present[(i+1)%2]);
          }
          den.push_back(local);
        }
      }
      return den;
    }

    vec_d stopping(double t){ 
      vec_d stopping_time;

      for (unsigned n = 0; n < sites.size(); n++){
        for (unsigned i = 0; i < 2*P.n_max; i++){
          if (sites[n].occupied[i] == true) stopping_time.push_back(t - sites[n].last_jump[i]);
        }
        
      }
      return stopping_time;
    }
};


// TODO not only 0 component
template<typename Engine>
class SnapshotWriter {
  const Lattice<Engine>& L;
public:
  SnapshotWriter(const Lattice<Engine>& L) : L(L) { }

  // Output pairs (d, a) indicating the direction d (numbered from 1) or 0 if site is vacant
  // and a indicating if the site is active

  friend std::ostream& operator << (std::ostream& out, const SnapshotWriter& SW) {

    for(const auto& site: SW.L.sites) {
      out << (site.occupied[0] ? unsigned(site.direction[0]+1) : 0);
      out << " " << site.active[0] << " ";
    }

    return out;
  }

};


// just a small output to see if the directions are properly implemented
template<typename Engine>
class HexDirectionWriter {
  const Hexagonal_lattice<Engine>& L;
public:
  HexDirectionWriter(const Hexagonal_lattice<Engine>& L,  ofstream& outfile) : L(L) { }

  // Output tuple (n, j, d), with n being lattice site, j being index in lattice site and 
  // d the lattice site where particle is pointing

  friend std::ostream& operator << (std::ostream& out, const HexDirectionWriter& HSW) {

    const auto& sites = HSW.L.sites; // Save typing
    

    for(unsigned n=0; n<sites.size(); ++n) {
      
      for (unsigned j = 0; j < 2 * 2; j++){
        if(sites[n].occupied[j]) out <<  n << "\t" << j << "\t" << sites[n].current_dir[j] << "\t" << int(sites[n].direction[j]);    
      }
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
      if(!sites[n].occupied[0]) out << sites[n].id[0] << " " << n << " ";
    }

    return out;
  }

};

template<typename Engine>
class ParticleWriter {
  const Lattice<Engine>& L;
public:
  ParticleWriter(const Lattice<Engine>& L, ofstream& outfile) : L(L) { }

  // Output pairs (n, p) indicating site and number of particles present

  friend std::ostream& operator << (std::ostream& out, const ParticleWriter& PW) {

    const auto& sites = PW.L.sites; // Save typing

    for(unsigned n=0; n<sites.size(); ++n) {
    
      if(sites[n].present > 0) out << n << " " << sites[n].present << " ";
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
      // do I really need id here?
      if(sites[n].present > 0) out << sites[n].id[0] << " " << n << " " << sites[n].present << " ";
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

  // Output is (v, n, j, p) as above with j being the binary site argument and p being number of
  // particles present
friend std::ostream& operator << (std::ostream& out, const HexagonalParticleWriter& HPW) {

    const auto& sites = HPW.L.sites; // Save typing

    for(unsigned n=0; n<sites.size(); ++n) {
      for (unsigned j = 0; j < 2; j++){
        if(sites[n].present[j] > 0) out << sites[n].id[j] << "\t" << n << "\t" << j << "\t" << sites[n].present[j] << "\t";    
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

  unsigned details = 0; // for detailed 

  app.add_option("-o, --output", output, "Output type: snapshots|particles|vacancies|clusters");
  app.add_option("-l, --lattice_type", lattice_type, "Lattice type: square|triangular|hexagonal");

  app.add_option("-b,--burnin",        burnin,        "Time to run before starting measurements");
  app.add_option("-u,--until",         until,         "Time to run for once measurements started");
  app.add_option("-e,--every",         every,         "Measurement interval");
  app.add_option("-a,--localaverage",  localaverage,  "Number of local averages (clusters only; 0=stationary state)");
  app.add_option("-i,--localinterval", localinterval, "Interval between local averages");
  app.add_option("-d,--details", details, "for detailed output do 1, else 0 (default)");

  CLI11_PARSE(app, argc, argv);

  if(output[0] == 'p') output = "particles";
  else if(output[0] == 'v') output = "vacancies";
  else if(output[0] == 'c') output = "clusters";
  else if(output[0] == 't') output = "stopping time"; // output for stopping time distribution
  else if(output[0] == 'd') output = "distribution"; // for distribution of neighbours (maybe also density?)
  else if(output[0] == 'l') output = "lagging"; // this is for output of hysteresis (l and lagging for greek roots of word)
  else if(output[0] == 'a') output = "area"; // for area/surface analysis
  else if(output[0] == 'w') output = "weighted"; // for weighted distribution
  else if(output[0] == 'm') output = "motility"; // for outputting the motility of the system
  else if(output[0] == 's') output = "stable"; // this is for making gifs that show how a system stabilizes...
  else if(output[0] == 'n') output = "number"; // this is for output of time evolution of cluster number and mean cluster size
  else if(output[0] == 'f') output = "function"; // this output is for testing different features. It is not static, as of now, so one should not uses this  without proper inspection of what it does
  else if(output[0] == 'h') output = "heatmap"; // this is for heatmap of clustersizes
  else if(output[0] == 'x') output = "perimeter"; // (first?) output for perimeter order parameter
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
  phi << double(P.N)/double(P.L[0] * P.L[1]);
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
  string txtoutput = path+lattice_type+"_"+tumb+"_"+ dens+"_"+size+"_"+occ_p+txt;
  string txtoutput_nr = path+lattice_type+"_nr"+"_"+tumb+"_"+ dens+"_"+size+"_"+occ_p+txt;

  if (P.n_max > 1){
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
      std::cout << "# occupation number = " << P.n_max << std::endl;

      double t = 0;


      if(output == "clusters") {
        if(localaverage == 0) {
          hist_t sumhist;
          hist_t sumhist_nr;
          // We sum the histograms over all measurements
          for(unsigned n=0; t < burnin + until; ++n) {
            t = L.run_until(burnin + n * every);
            hist_t hist = L.cluster_distributions();
            hist_t hist_nr = L.cluster_distribution_particle_number();

            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }

            // with occ number as well

            if(hist_nr.size() > sumhist_nr.size()){
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            } else {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }

          }
          ofstream outfile, outfile_nr;
          outfile.open(txtoutput);
          for(const auto& k: sumhist) outfile << k << " ";
          outfile << endl;
          outfile_nr.open(txtoutput_nr);
          for(const auto& k: sumhist_nr) outfile_nr << k << " ";
          outfile_nr << endl;
        } else {
          // We perform local averages of the histograms at the given measurement interval
          for(unsigned n=0; t < burnin + until; ++n) {
            t = L.run_until(burnin + n * every);
            hist_t sumhist = L.cluster_distributions();
            hist_t sumhist_nr = L.cluster_distribution_particle_number();
            for(unsigned m=1; m<localaverage; ++m) {
              t = L.run_until(burnin + n*every + m*localinterval);
              hist_t hist = L.cluster_distributions();
              hist_t hist_nr = L.cluster_distribution_particle_number();
              // Add hist to sumhist taking into account they may have different lengths
              if(hist.size() > sumhist.size()) {
                hist[std::slice(0,sumhist.size(),1)] += sumhist;
                sumhist = std::move(hist);
              } else {
                sumhist[std::slice(0,hist.size(),1)] += hist;
              }
              if(hist_nr.size() > sumhist_nr.size()){
                hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
                sumhist_nr = std::move(hist_nr);
              } else {
                sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
              }
            }
            ofstream outfile, outfile_nr;
            outfile.open(txtoutput);
            for(const auto& k: sumhist) outfile << k << " ";
            outfile << endl;
            outfile_nr.open(txtoutput_nr);
            for(const auto& k: sumhist_nr) outfile_nr << k << " ";
            outfile_nr << endl;
          }
        } 
      } 
      else if (output == "heatmap"){
        ofstream outfile, outfile_avg, outfile_nr, outfile_avg_nr;
        outfile.open("./lars_sim/heatmap/square_alpha_N_n_3.txt");
        outfile_avg.open("./lars_sim/heatmap/square_alpha_N_n_3_avg.txt");
        outfile_nr.open("./lars_sim/heatmap/square_nr_alpha_N_n_3.txt");
        outfile_avg_nr.open("./lars_sim/heatmap/square_nr_alpha_N_n_3_avg.txt");
        for (double alp = 1e-6; alp <= 1.0; alp*=1.3){
          for (unsigned N = 100 * P.n_max; N <= P.L[0]*P.L[0]*P.n_max; N += 200*P.n_max){
            Parameters P_h;
            P_h.N = N;
            P_h.alpha[0] = P.alpha[1] = alp;
            P_h.L = P.L;
            P_h.n_max = P.n_max;
            std::size_t maxsize = 1;
            std::size_t maxsize_nr = 1;
            Lattice LB(P_h, rng);
            t = 0;
            // also taking mean over all used timesteps
            double mean = 0.0;
            double mean_nr = 0.0;
            double count = 0;
            for(unsigned n=0; t < burnin + until; ++n){
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
      else if (output == "number"){
        ofstream outfile;
        string name = "./lars_sim/number/square_";
        string name_number = name+"number"+"_"+tumb+"_"+dens+"_"+size+"_"+occ_p+txt;
        outfile.open(name_number);


        t = L.run_until(burnin);
        hist_t hist = L.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i+=2){
          second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
          first_moment += hist[i] * ((i+2)/2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        outfile << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
        for(double n=1; t < burnin + until; n*=1.01) {
          t = L.run_until(burnin + n * every);
          hist_t hist = L.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i+=2){
            second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
            first_moment += hist[i] * ((i+2)/2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          outfile << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
        }
      
      } 
      else if (output == "lagging"){
        ofstream outfile, backward;
        string name = "./lars_sim/Data/motility/square_perc_fhyst";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        name = "./lars_sim/Data/motility/square_perc_bhyst";
        outputname = name+"_"+occ_p+".txt";
        backward.open(outputname);
        // foward hysteresis, i.e. start below critical point and move up
        Lattice LB(P, rng);
        double tmax = burnin + until;
        unsigned c = 0;
        for (double al = 0.0875; al < 0.094 ; al+=0.0001625){
          // introducing new alpha
          P.alpha[0] = P.alpha[1] = al;
          LB.set_new_lambda(&LB.tumble, std::accumulate(P.alpha.begin(), P.alpha.end(), 0.0) / P.alpha.size());
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          std::vector<double> values_wei;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          double weighted = 0;
          double cov_w = 0;
          for(unsigned n=0; t < (c+1)*tmax; ++n) {
            t = LB.run_until(burnin + n * every + c*tmax);

            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distribution_particle_number();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          
          }
           mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);
          

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
          c++;
        }

        Lattice LT(P, rng);
        c = 0;
        t = 0;
        for (double al = 0.094; al > 0.0875 ; al-=0.0001625){
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
          for(unsigned n=0; t < (c+1)*tmax; ++n) {
            t = LT.run_until(burnin + n * every + c*tmax);

            values_mas_b.push_back(double(LT.max_cluster_size_nr())/double(P.N));
            values_mot_b.push_back(LT.motility_fraction());
            rel_mass_b += double(LT.max_cluster_size_nr())/double(P.N);
            mean_b += LT.motility_fraction();
            count_b++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist_b = LT.cluster_distribution_particle_number();
            double second_moment_b = 0;
            double first_moment_b = 0;
            for (unsigned i = 0; i < hist_b.size(); i+=2){
              second_moment_b += hist_b[i] * ((i+2)/2) * ((i+2)/2);
              first_moment_b += hist_b[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei_b.push_back(second_moment_b / first_moment_b * 1.0/double(P.N));
            weighted_b += second_moment_b / first_moment_b * 1.0/double(P.N);
          
          }
           mean_b = mean_b / count_b;
          rel_mass_b = rel_mass_b / count_b;
          weighted_b= weighted_b / count_b;


          double cov_mot_b = 0;
          for (auto& val : values_mot_b){
            cov_mot_b += pow(val - mean_b, 2);
          }
          cov_mot_b = cov_mot_b/(values_mot_b.size() - 1);

          double cov_mas_b = 0;
          for (auto& val : values_mas_b){
            cov_mas_b += pow(val - rel_mass_b, 2);
          }
          cov_mas_b = cov_mas_b/(values_mas_b.size() - 1);

          for (auto& val : values_wei_b){
            cov_w_b += pow(val - weighted_b, 2);
          }
          cov_w_b = cov_w_b/(values_wei_b.size() - 1);

          backward << al << " " << mean_b << " " << cov_mot_b << " " << rel_mass_b << " " << cov_mas_b << " " << weighted_b << " " << cov_w_b << endl;
          c++;
        }
      }
      else if (output == "motility"){
        if (details==0){
        ofstream outfile;
        string name = "./lars_sim/Data/motility/square_perc_low";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        for (double al = 0.0; al < 0.05 ; al+=0.001){
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1]  = al;
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
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distribution_particle_number();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
      
      } else if (details==1){
        ofstream outfile;
        string name = "./lars_sim/Data/motility/square_perc_details";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        for (double al = 0.085; al < 0.095 ; al+=0.0005){
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1]  = al;
          Lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;
          }
          mean = mean / count;
          rel_mass = rel_mass / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << endl;
        }

      }
      }
      else if (output == "stable"){
        ofstream part, numb;
        part.open("./lars_sim/Data/stable/square.txt");
        numb.open("./lars_sim/Data/stable/square_number.txt");
        
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          // only doing a positional output here
          numb << t << " " << L.number_cluster() << endl;
          part << ParticleWriter(L, part) << endl;
        }
      }
      else if (output == "weighted"){
         ofstream output, part;
         part.open("./lars_sim/Data/weighted/square_part.txt");
         output.open("./lars_sim/Data/weighted/square.txt");
         
         for(double n=1; t < burnin + until; n++) {
            t = L.run_until(burnin + n * every);
            part << ParticleWriter(L, part) << endl;
            hist_t hist = L.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            double weighted = second_moment / first_moment;
            weighted = weighted / double(P.N);
            output << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
         }
      }
      else if (output == "distribution"){
        ofstream outfile, outfile2;
        outfile.open("./lars_sim/Data/dist/square_"+occ_p+".txt");
        outfile2.open("./lars_sim/Data/dist/square_dens_"+occ_p+".txt");
        hist_t dist(5);
        for(double n = 0; t < burnin + until; n++) {
          t = L.run_until(burnin + n * every);
          vec_d dens = L.density();
          for (const auto& m : dens) outfile2 << m << " ";
          outfile2 << endl;
          hist_t dr = L.particle_neighbour_dist();
          dist[std::slice(0,dr.size(),1)] += dr;
        }
        for (const auto& m : dist) outfile << m << " ";
        outfile << endl;
      }
      else if (output=="area"){
        ofstream surf, part, clust, border;
        surf.open("./lars_sim/Data/surf/square_sv_low"+occ_p+".txt");
        part.open("./lars_sim/Data/surf/square_part_low"+occ_p+".txt");
        clust.open("./lars_sim/Data/surf/square_clust_low"+occ_p+".txt");
        border.open("./lars_sim/Data/surf/square_border_low"+occ_p+".txt");

        // unsigned so only one cluster each gets printed
        unsigned check_1 = 0;
        unsigned check_2 = 0;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          vec a = L.surface_volume();
          for (const auto& k : a) surf << k << " ";
          surf << endl;
          part << ParticleWriter(L, part) << endl;
          for (const auto& n : L.cluster_surface()) border << n << " ";
          border << endl;

          if (L.clust_size(8, 10) == 1 && check_1 == 0){
            vec n =L.single_cluster(8, 10);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_1++;
          }

          if (L.clust_size(3700, 4000) == 1 && check_2 == 0){
            vec n = L.single_cluster(3700, 4000);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_2++;
          }

        }
      }
      else if (output == "stopping time"){
        ofstream outfile;
        
        outfile.open("./lars_sim/Data/stopping/square_"+occ_p+".txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          vec_d st = L.stopping(t);
          for (const auto& m : st) outfile << m << " ";
          outfile << endl;
        }
      }
      else {
        ofstream outfile;
        outfile.open("./lars_sim/gif/square.txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          if (output == "particles") outfile << ParticleWriter(L, outfile) << std::endl;
          else if(output == "vacancies") outfile << VacancyWriter(L) << std::endl;
          else {
            // Ensure that all the particle directions are at the simulation time
            L.realise_directions();
            outfile << SnapshotWriter(L) << std::endl;
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
      std::cout << "# occupation number = " << P.n_max << std::endl;

      double t = 0;

      if(output == "clusters") {
        if(localaverage == 0) {
          // We sum the histograms over all measurements
          hist_t sumhist;
          hist_t sumhist_nr;
          for(unsigned n=0; t < burnin + until; ++n) {
            t = TL.run_until(burnin + n * every);
            hist_t hist = TL.cluster_distributions();
            hist_t hist_nr = TL.cluster_distributions_particle_numbers();
            // Only area
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }

            // with occ number as well

            if(hist_nr.size() > sumhist_nr.size()){
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            } else {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }

          }
          // output for each of the distributions 
          ofstream outfile, outfile_nr;
          outfile.open(txtoutput);
          for(const auto& k: sumhist) outfile << k << " ";
          outfile << endl;
          outfile_nr.open(txtoutput_nr);
          for(const auto& k: sumhist_nr) outfile_nr << k << " ";
          outfile_nr << endl;
        } else {
          // We perform local averages of the histograms at the given measurement interval
          for(unsigned n=0; t < burnin + until; ++n) {
            t = TL.run_until(burnin + n * every);
            hist_t sumhist = TL.cluster_distributions();
            hist_t sumhist_nr = TL.cluster_distributions_particle_numbers();
            for(unsigned m=1; m<localaverage; ++m) {
              t = TL.run_until(burnin + n*every + m*localinterval);
              hist_t hist = TL.cluster_distributions();
              hist_t hist_nr = TL.cluster_distributions_particle_numbers();
              // Add hist to sumhist taking into account they may have different lengths
              if(hist.size() > sumhist.size()) {
                hist[std::slice(0,sumhist.size(),1)] += sumhist;
                sumhist = std::move(hist);
              } else {
                sumhist[std::slice(0,hist.size(),1)] += hist;
              }
              if(hist_nr.size() > sumhist_nr.size()){
                hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
                sumhist_nr = std::move(hist_nr);
              } else {
                sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
              }
            }
            ofstream outfile, outfile_nr;
            outfile.open(txtoutput);
            for(const auto& k: sumhist) outfile << k << " ";
            outfile << endl;
            outfile_nr.open(txtoutput_nr);
            for(const auto& k: sumhist_nr) outfile_nr << k << " ";
            outfile_nr << endl;
          }
        }
      } 
      else if (output == "heatmap"){
        ofstream outfile, outfile_avg, outfile_nr, outfile_avg_nr;
        outfile.open("./lars_sim/heatmap/tri_alpha_N_n_2.txt");
        outfile_avg.open("./lars_sim/heatmap/tri_alpha_N_n_2_avg.txt");
        outfile_nr.open("./lars_sim/heatmap/tri_nr_alpha_N_n_2.txt");
        outfile_avg_nr.open("./lars_sim/heatmap/tri_nr_alpha_N_n_2_avg.txt");
        for (double alp = 1e-6; alp <= 1.0; alp*=1.3){
          for (unsigned N = 100*P.n_max; N <= P.L[0]*P.L[0]*P.n_max; N += 200*P.n_max){
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
            for(unsigned n=0; t < burnin + until; ++n){
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
      else if (output == "function"){
        ofstream outfile;
        outfile.open("./lars_sim/Data/motility/tri_long.txt");
        for(unsigned n = 0; t < burnin + until; n++) {
            t = TL.run_until(burnin + n * every);
            outfile << t << " " << TL.motility_fraction() << " " << double(TL.max_cluster_size_nr())/double(P.N) << endl;
        }
      }
      else if (output == "stopping time"){
        ofstream outfile;
        
        outfile.open("./lars_sim/Data/stopping/tri_"+occ_p+".txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          vec_d st = TL.stopping(t);
          for (const auto& m : st) outfile << m << " ";
          outfile << endl;
        }
      }
      else if (output == "number"){
        ofstream outfile;
        string name = "./lars_sim/number/tri_";
        string name_number = name+"number"+"_"+tumb+"_"+dens+"_"+size+"_"+occ_p+txt;
        outfile.open(name_number);


        t = TL.run_until(burnin);
        hist_t hist = TL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i+=2){
          second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
          first_moment += hist[i] * ((i+2)/2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        outfile << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
        for(double n=1; t < burnin + until; n*=1.01) {
          t = TL.run_until(burnin + n * every);
          hist_t hist = TL.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i+=2){
            second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
            first_moment += hist[i] * ((i+2)/2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          outfile << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
        }
      
      }
      else if (output == "weighted"){
         ofstream output, part;
         part.open("./lars_sim/Data/weighted/tri_part.txt");
         output.open("./lars_sim/Data/weighted/tri.txt");
         
         for(double n=1; t < burnin + until; n++) {
            t = TL.run_until(burnin + n * every);
            part << TriangleParticleWriter(TL, part) << endl;
            hist_t hist = TL.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            double weighted = second_moment / first_moment;
            weighted = weighted / double(P.N);
            output << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
         }
      }
      else if (output == "stable"){
        ofstream part, numb;
        part.open("./lars_sim/Data/stable/triangular.txt");
        numb.open("./lars_sim/Data/stable/triangular_number.txt");
        
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          // only doing a positional output here
          numb << t << " " << TL.number_cluster() << endl;
          part << TriangleParticleWriter(TL, part) << endl;
        }
      }else if (output == "motility"){
        if (details==0){
        ofstream outfile;
        string name = "./lars_sim/Data/motility/triangular_perc_low";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        for (double al = 0.0; al < 0.05 ; al+=0.001){
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
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distributions_particle_numbers();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
        } else if (details==1){
          ofstream outfile;
        string name = "./lars_sim/Data/motility/triangular_perc_details";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        for (double al = 0.06; al < 0.07 ; al+=0.0005){
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
          Triangle_lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;
          }
          mean = mean / count;
          rel_mass = rel_mass / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << endl;
        }

        }
      }
      else if (output=="area"){
        ofstream surf, part, clust, border;
        surf.open("./lars_sim/Data/surf/tri_sv_low"+occ_p+".txt");
        part.open("./lars_sim/Data/surf/tri_part_low"+occ_p+".txt");
        clust.open("./lars_sim/Data/surf/tri_clust_low"+occ_p+".txt");
        border.open("./lars_sim/Data/surf/tri_border_low"+occ_p+".txt");
        // unsigned so only one cluster each gets printed
        //unsigned check_1 = 0;
        //unsigned check_2 = 0;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          vec a = TL.surface_volume();
          for (const auto& k : a) surf << k << " ";
          surf << endl;
          part << TriangleParticleWriter(TL, part) << endl;
          for (const auto& n : TL.cluster_surface()) border << n << " ";
          border << endl;
          /*
          if (TL.clust_size(8, 10) == 1 && check_1 == 0){
            vec n = TL.single_cluster(8, 10);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_1++;
          }

          if (TL.clust_size(1470, 1520) == 1 && check_2 == 0){
            vec n = TL.single_cluster(1470, 1520);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_2++;
          }
          */
        }
      }
      else if (output == "lagging"){
        ofstream outfile, backward;
        string name = "./lars_sim/Data/motility/triangular_perc_fhyst";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        name = "./lars_sim/Data/motility/triangular_perc_bhyst";
        outputname = name+"_"+occ_p+".txt";
        backward.open(outputname);
        // foward hysteresis, i.e. start below critical point and move up
        Triangle_lattice LB(P, rng);
        double tmax = burnin + until;
        unsigned c = 0;
        for (double al = 0.0625; al < 0.07 ; al+=0.0001875){
          // introducing new alpha and updating tumble distribution.

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
          for(unsigned n=0; t < (c+1)*tmax; ++n) {
            t = LB.run_until(burnin + n * every + c*tmax);

            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distributions_particle_numbers();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          
          }
           mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);
          

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
          c++;
        }

        Triangle_lattice LT(P, rng);
        c = 0;
        t = 0;
        for (double al = 0.07; al > 0.0625 ; al-=0.0001875){
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
          for(unsigned n=0; t < (c+1)*tmax; ++n) {
            t = LT.run_until(burnin + n * every + c*tmax);

            values_mas_b.push_back(double(LT.max_cluster_size_nr())/double(P.N));
            values_mot_b.push_back(LT.motility_fraction());
            rel_mass_b += double(LT.max_cluster_size_nr())/double(P.N);
            mean_b += LT.motility_fraction();
            count_b++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist_b = LT.cluster_distributions_particle_numbers();
            double second_moment_b = 0;
            double first_moment_b = 0;
            for (unsigned i = 0; i < hist_b.size(); i+=2){
              second_moment_b += hist_b[i] * ((i+2)/2) * ((i+2)/2);
              first_moment_b += hist_b[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei_b.push_back(second_moment_b / first_moment_b * 1.0/double(P.N));
            weighted_b += second_moment_b / first_moment_b * 1.0/double(P.N);
          
          }
           mean_b = mean_b / count_b;
          rel_mass_b = rel_mass_b / count_b;
          weighted_b= weighted_b / count_b;


          double cov_mot_b = 0;
          for (auto& val : values_mot_b){
            cov_mot_b += pow(val - mean_b, 2);
          }
          cov_mot_b = cov_mot_b/(values_mot_b.size() - 1);

          double cov_mas_b = 0;
          for (auto& val : values_mas_b){
            cov_mas_b += pow(val - rel_mass_b, 2);
          }
          cov_mas_b = cov_mas_b/(values_mas_b.size() - 1);

          for (auto& val : values_wei_b){
            cov_w_b += pow(val - weighted_b, 2);
          }
          cov_w_b = cov_w_b/(values_wei_b.size() - 1);

          backward << al << " " << mean_b << " " << cov_mot_b << " " << rel_mass_b << " " << cov_mas_b << " " << weighted_b << " " << cov_w_b << endl;
          c++;
        }
      }
      else if (output == "distribution"){
        ofstream outfile, outfile2;
        outfile.open("./lars_sim/Data/dist/tri_"+occ_p+".txt");
        outfile2.open("./lars_sim/Data/dist/tri_dens_"+occ_p+".txt");
        hist_t dist(7);
        for(double n = 0; t < burnin + until; n++) {
          
          t = TL.run_until(burnin + n * every);
          vec_d dens = TL.density();
          for (const auto& y : dens) outfile2 << y << " ";
          outfile2 << endl;
          
          hist_t dr = TL.particle_neighbour_dist();
          dist[std::slice(0,dr.size(),1)] += dr;
          
        }
        
        for (const auto& m : dist) outfile << m << " ";
        outfile << endl;
      }
      else {
        ofstream outfile;
        outfile.open("./lars_sim/gif/triangle.txt");
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
      std::cout << "# occupation number = " << P.n_max << std::endl;

      double t = 0;

      if(output == "clusters") {
        if(localaverage == 0) {
          hist_t sumhist;
          hist_t sumhist_nr;
          // We sum the histograms over all measurements
          ofstream outfile_part;
          outfile_part.open("./lars_sim/gif/hexagonal.txt");

          for(unsigned n=0; t < burnin + until; ++n) {
            t = HL.run_until(burnin + n * every);
            hist_t hist = HL.cluster_distributions();
            hist_t hist_nr = HL.cluster_distribution_particle_number();
            // particle output for comparison
            outfile_part << HexagonalParticleWriter(HL, outfile_part) << endl;
            // Only area
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }

            // with occ number as well

            if(hist_nr.size() > sumhist_nr.size()){
              hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
              sumhist_nr = std::move(hist_nr);
            } else {
              sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
            }

          }
          ofstream outfile, outfile_nr;
          outfile.open(txtoutput);
          for(const auto& k: sumhist) outfile << k << " ";
          outfile << endl;
          outfile_nr.open(txtoutput_nr);
          for(const auto& k: sumhist_nr) outfile_nr << k << " ";
          outfile_nr << endl;
        } else {
          // We perform local averages of the histograms at the given measurement interval
          for(unsigned n=0; t < burnin + until; ++n) {
            t = HL.run_until(burnin + n * every);
            hist_t sumhist = HL.cluster_distributions();
            hist_t sumhist_nr = HL.cluster_distribution_particle_number();
            for(unsigned m=1; m<localaverage; ++m) {
              t = HL.run_until(burnin + n*every + m*localinterval);
              hist_t hist = HL.cluster_distributions();
              hist_t hist_nr = HL.cluster_distribution_particle_number();
              // Add hist to sumhist taking into account they may have different lengths
              if(hist.size() > sumhist.size()) {
                hist[std::slice(0,sumhist.size(),1)] += sumhist;
                sumhist = std::move(hist);
              } else {
                sumhist[std::slice(0,hist.size(),1)] += hist;
              }
              if(hist_nr.size() > sumhist_nr.size()){
                hist_nr[std::slice(0, sumhist_nr.size(), 1)] += sumhist_nr;
                sumhist_nr = std::move(hist_nr);
              } else {
                sumhist_nr[std::slice(0, hist_nr.size(), 1)] += hist_nr;
              }
            }
            ofstream outfile, outfile_nr;
            outfile.open(txtoutput);
            for(const auto& k: sumhist) outfile << k << " ";
            outfile << endl;
            outfile_nr.open(txtoutput_nr);
            for(const auto& k: sumhist_nr) outfile_nr << k << " ";
            outfile_nr << endl;
          }
        }
      } else if (output == "particles") {
        ofstream outfile;
        outfile.open("./lars_sim/gif/hexagonal.txt");

        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          // only doing a positional output here
          
          outfile << HexagonalParticleWriter(HL, outfile) << endl;

        }

      }
      else if (output == "lagging"){
        ofstream outfile, backward;
        string name = "./lars_sim/Data/motility/hexagonal_perc_fhyst";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        name = "./lars_sim/Data/motility/hexagonal_perc_bhyst";
        outputname = name+"_"+occ_p+".txt";
        backward.open(outputname);
        // foward hysteresis, i.e. start below critical point and move up
        Hexagonal_lattice LB(P, rng);
        double tmax = burnin + until;
        unsigned c = 0;
        for (double al = 0.10; al < 0.126 ; al+=0.002){
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
          for(unsigned n=0; t < (c+1)*tmax; ++n) {
            t = LB.run_until(burnin + n * every + c*tmax);

            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distribution_particle_number();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          
          }
           mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);
          

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
          c++;
        }

        Hexagonal_lattice LT(P, rng);
        c = 0;
        t = 0;
        for (double al = 0.126; al > 0.10 ; al-=0.002){
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
          for(unsigned n=0; t < (c+1)*tmax; ++n) {
            t = LT.run_until(burnin + n * every + c*tmax);

            values_mas_b.push_back(double(LT.max_cluster_size_nr())/double(P.N));
            values_mot_b.push_back(LT.motility_fraction());
            rel_mass_b += double(LT.max_cluster_size_nr())/double(P.N);
            mean_b += LT.motility_fraction();
            count_b++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist_b = LT.cluster_distribution_particle_number();
            double second_moment_b = 0;
            double first_moment_b = 0;
            for (unsigned i = 0; i < hist_b.size(); i+=2){
              second_moment_b += hist_b[i] * ((i+2)/2) * ((i+2)/2);
              first_moment_b += hist_b[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei_b.push_back(second_moment_b / first_moment_b * 1.0/double(P.N));
            weighted_b += second_moment_b / first_moment_b * 1.0/double(P.N);
          
          }
           mean_b = mean_b / count_b;
          rel_mass_b = rel_mass_b / count_b;
          weighted_b= weighted_b / count_b;


          double cov_mot_b = 0;
          for (auto& val : values_mot_b){
            cov_mot_b += pow(val - mean_b, 2);
          }
          cov_mot_b = cov_mot_b/(values_mot_b.size() - 1);

          double cov_mas_b = 0;
          for (auto& val : values_mas_b){
            cov_mas_b += pow(val - rel_mass_b, 2);
          }
          cov_mas_b = cov_mas_b/(values_mas_b.size() - 1);

          for (auto& val : values_wei_b){
            cov_w_b += pow(val - weighted_b, 2);
          }
          cov_w_b = cov_w_b/(values_wei_b.size() - 1);

          backward << al << " " << mean_b << " " << cov_mot_b << " " << rel_mass_b << " " << cov_mas_b << " " << weighted_b << " " << cov_w_b << endl;
          c++;
        }
      }

      else if (output == "motility"){
        if (details==0){
        ofstream outfile;
        string name = "./lars_sim/Data/motility/hexagonal_perc_low";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        for (double al = 0.0; al < 0.05 ; al+=0.001){
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
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distribution_particle_number();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        
        }
      
      } else if (details==1){
        ofstream outfile;
        string name = "./lars_sim/Data/motility/hexagonal_perc_details";
        string outputname = name+"_"+occ_p+".txt";
        outfile.open(outputname);
        for (double al = 0.105; al < 0.12 ; al+=0.0005){
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1] = P.alpha[2] = al;
          Hexagonal_lattice LB(P, rng);
          t = 0;
          std::vector<double> values_mot;
          std::vector<double> values_mas;
          double mean = 0;
          double rel_mass = 0;
          double count = 0;
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;
          }
          mean = mean / count;
          rel_mass = rel_mass / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << endl;
        }
      }
      }
      else if (output == "stable"){
        ofstream part, numb;
        part.open("./lars_sim/Data/stable/hexagonal.txt");
        numb.open("./lars_sim/Data/stable/hexagonal_number.txt");
        
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          // only doing a positional output here
          numb << t << " " << HL.number_cluster() << endl;
          part << HexagonalParticleWriter(HL, part) << endl;
        }
      }
      else if (output == "weighted"){
         ofstream output, part;
         part.open("./lars_sim/Data/weighted/hex_part.txt");
         output.open("./lars_sim/Data/weighted/hex.txt");
         
         for(double n=1; t < burnin + until; n++) {
            t = HL.run_until(burnin + n * every);
            part << HexagonalParticleWriter(HL, part) << endl;
            hist_t hist = HL.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            double weighted = second_moment / first_moment;
            weighted = weighted / double(P.N);
            output << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
         }
      }
      else if (output=="area"){
        ofstream surf, part, clust, border;
        surf.open("./lars_sim/Data/surf/hex_sv_low"+occ_p+".txt");
        part.open("./lars_sim/Data/surf/hex_part_low"+occ_p+".txt");
        //clust.open("./lars_sim/Data/surf/hex_clust_low"+occ_p+".txt");
        border.open("./lars_sim/Data/surf/hex_border_low"+occ_p+".txt");
        // unsigned so only one cluster each gets printed
        //unsigned check_1 = 0;
        //unsigned check_2 = 0;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          vec a = HL.surface_volume();
          for (const auto& k : a) surf << k << " ";
          surf << endl;
          part << HexagonalParticleWriter(HL, part) << endl;
          for (const auto& n : HL.cluster_surface()) border << n << " ";
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
      else if (output == "heatmap"){
        ofstream outfile, outfile_avg, outfile_nr, outfile_avg_nr;
        outfile.open("./lars_sim/heatmap/hex_alpha_N_n_2.txt");
        outfile_avg.open("./lars_sim/heatmap/hex_alpha_N_n_2_avg.txt");
        outfile_nr.open("./lars_sim/heatmap/hex_nr_alpha_N_n_2.txt");
        outfile_avg_nr.open("./lars_sim/heatmap/hex_nr_alpha_N_n_2_avg.txt");
        for (double alp = 1e-6; alp <= 1.0; alp*=1.3){
          for (unsigned N = 200 * P.n_max; N <= 2*P.L[0]*P.L[0]*P.n_max; N += 400*P.n_max){
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
            for(unsigned n=0; t < burnin + until; ++n){
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
      else if (output == "number"){
        ofstream outfile;
        string name = "./lars_sim/number/hex_";
        string name_number = name+"number"+"_"+tumb+"_"+dens+"_"+size+"_"+occ_p+txt;
        outfile.open(name_number);


        t = HL.run_until(burnin);
        hist_t hist = HL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i+=2){
          second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
          first_moment += hist[i] * ((i+2)/2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        outfile << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
        for(double n=1; t < burnin + until; n*=1.01) {
          t = HL.run_until(burnin + n * every);
          hist_t hist = HL.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i+=2){
            second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
            first_moment += hist[i] * ((i+2)/2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          outfile << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
        }
      } else if (output == "snapshots"){
        ofstream outfile;
        outfile.open("./lars_sim/gif/hexdir.txt");

        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          // only doing a positional output here
          HL.realise_directions();
          outfile << HexDirectionWriter(HL, outfile) << endl;
        }

      }
      else if (output == "stopping time"){
        ofstream outfile;
        
        outfile.open("./lars_sim/Data/stopping/hex_"+occ_p+".txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          vec_d st = HL.stopping(t);
          for (const auto& m : st) outfile << m << " ";
          outfile << endl;
        }
      }
      else if (output == "function"){
        ofstream outfile;
        outfile.open("./lars_sim/testing/hexagonal.txt");
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
      else if (output == "distribution"){
        ofstream outfile, outfile2;
        outfile.open("./lars_sim/Data/dist/hex_"+occ_p+".txt");
        outfile2.open("./lars_sim/Data/dist/hex_dens_"+occ_p+".txt");
        hist_t dist(4);
        for(double n = 0; t < burnin + until; n++) {
          t = HL.run_until(burnin + n * every);
          vec_d dens = HL.density();
          for (const auto& m : dens) outfile2 << m << " ";
          outfile2 << endl;
          hist_t dr = HL.particle_neighbour_dist();
          dist[std::slice(0,dr.size(),1)] += dr;
        }
        for (const auto& m : dist) outfile << m << " ";
        outfile << endl;
      }

    }
  } else{
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
      std::cout << "# occupation number = " << P.n_max << std::endl;

      double t = 0;


      if(output == "clusters") {
        if(localaverage == 0) {
          hist_t sumhist;
          // We sum the histograms over all measurements
          for(unsigned n=0; t < burnin + until; ++n) {
            t = L.run_until(burnin + n * every);
            hist_t hist = L.cluster_distributions();

            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }

          }
          ofstream outfile;
          outfile.open(txtoutput);
          for(const auto& k: sumhist) outfile << k << " ";
          outfile << endl;
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
            ofstream outfile, outfile_nr;
            outfile.open(txtoutput);
            for(const auto& k: sumhist) outfile << k << " ";
            outfile << endl;
          }
        } 
      }
      else if (output=="area"){
        ofstream surf, part, clust, border;
        surf.open("./lars_sim/Data/surf/square_sv_low"+occ_p+".txt");
        part.open("./lars_sim/Data/surf/square_part_low"+occ_p+".txt");
        clust.open("./lars_sim/Data/surf/square_clust_low"+occ_p+".txt");
        border.open("./lars_sim/Data/surf/square_border_low"+occ_p+".txt");

        // unsigned so only one cluster each gets printed
        unsigned check_1 = 0;
        unsigned check_2 = 0;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          vec a = L.surface_volume();
          for (const auto& k : a) surf << k << " ";
          surf << endl;
          part << ParticleWriter(L, part) << endl;
          for (const auto& n : L.cluster_surface()) border << n << " ";
          border << endl;

          if (L.clust_size(8, 10) == 1 && check_1 == 0){
            vec n =L.single_cluster(8, 10);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_1++;
          }

          if (L.clust_size(3700, 4000) == 1 && check_2 == 0){
            vec n = L.single_cluster(3700, 4000);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_2++;
          }

        }
      }
      else if (output == "weighted"){
         ofstream output, part;
         part.open("./lars_sim/Data/weighted/square_part.txt");
         output.open("./lars_sim/Data/weighted/square.txt");
         
         for(double n=1; t < burnin + until; n++) {
            t = L.run_until(burnin + n * every);
            part << ParticleWriter(L, part) << endl;
            hist_t hist = L.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            double weighted = second_moment / first_moment;
            weighted = weighted / double(P.N);
            output << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
         }
      }
      else if (output == "stable"){
        ofstream part, numb;
        part.open("./lars_sim/Data/stable/square.txt");
        numb.open("./lars_sim/Data/stable/square_number.txt");
        
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          // only doing a positional output here
          numb << t << " " << L.number_cluster() << endl;
          part << ParticleWriter(L, part) << endl;
        }
      }
      else if (output == "stopping time"){
        ofstream outfile;
        
        outfile.open("./lars_sim/Data/stopping/square_"+occ_p+".txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          vec_d st = L.stopping(t);
          for (const auto& m : st) outfile << m << " ";
          outfile << endl;
        }
      }
      else if (output == "number"){
        ofstream outfile;
        string name = "./lars_sim/number/square_";
        string name_number = name+"number"+"_"+tumb+"_"+dens+"_"+size+txt;
        outfile.open(name_number);


        t = L.run_until(burnin);
        hist_t hist = L.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i+=2){
          second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
          first_moment += hist[i] * ((i+2)/2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        outfile << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
        for(double n=1; t < burnin + until; n*=1.01) {
          t = L.run_until(burnin + n * every);
          hist_t hist = L.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i+=2){
            second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
            first_moment += hist[i] * ((i+2)/2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          outfile << t << " " << weighted << " " << L.avg_cluster_size_nr() << endl;
        }
      }
      else if (output == "heatmap"){
        ofstream outfile, outfile_avg;
        outfile.open("./lars_sim/heatmap/square_alpha_N.txt");
        outfile_avg.open("./lars_sim/heatmap/square_alpha_N_avg.txt");
        for (double alp = 1e-6; alp <= 1.0; alp*=1.2){
          for (unsigned N = 100; N <= P.L[0]*P.L[0]; N += 200){
            Parameters P_h;
            P_h.N = N;
            P_h.alpha[0] = P.alpha[1] = alp;
            P_h.L = P.L;
            P_h.n_max = P.n_max;
            std::size_t maxsize = 1;
            Lattice LB(P_h, rng);
            t = 0;

            // also taking mean over all used timesteps
            double mean = 0.0;
            double count = 0;
            for(unsigned n=0; t < burnin + until; ++n){
              t = LB.run_until(burnin + n * every);
              maxsize = std::max(maxsize, LB.max_cluster_size());
              mean += LB.avg_cluster_size();
              count += 1;
            }

            mean = mean / count;

            outfile << maxsize << " ";
            outfile_avg << mean << " ";
          }
          outfile << endl;
          outfile_avg << endl;
        } 
      }else if (output == "motility"){
        ofstream outfile;
        outfile.open("./lars_sim/Data/motility/square_perc_low.txt");
        for (double al = 0.0; al < 0.2 ; al+=0.005){
          // defining lattice for new alpha
          P.alpha[0] = P.alpha[1]  = al;
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
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
      
      }
      else if (output == "distribution"){
        ofstream outfile, outfile2;
        outfile.open("./lars_sim/Data/dist/square_"+occ_p+".txt");
        outfile2.open("./lars_sim/Data/dist/square_dens_"+occ_p+".txt");
        hist_t dist(5);
        for(double n = 0; t < burnin + until; n++) {
          t = L.run_until(burnin + n * every);
          vec_d dens = L.density();
          for (const auto& m : dens) outfile2 << m << " ";
          outfile2 << endl;
          hist_t dr = L.particle_neighbour_dist();
          dist[std::slice(0,dr.size(),1)] += dr;
        }
        for (const auto& m : dist) outfile << m << " ";
        outfile << endl;
      }
      else {
        ofstream outfile;
        outfile.open("./lars_sim/gif/square.txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = L.run_until(burnin + n * every);
          if (output == "particles") std::cout << ParticleWriter(L, outfile) << std::endl;
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
      std::cout << "# occupation number = " << P.n_max << std::endl;

      double t = 0;

      if(output == "clusters") {
        if(localaverage == 0) {
          // We sum the histograms over all measurements
          hist_t sumhist;
          for(unsigned n=0; t < burnin + until; ++n) {
            t = TL.run_until(burnin + n * every);
            hist_t hist = TL.cluster_distributions();
            // Only area
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }

          }
          // output for each of the distributions 
          ofstream outfile;
          outfile.open(txtoutput);
          for(const auto& k: sumhist) outfile << k << " ";
          outfile << endl;
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
            ofstream outfile, outfile_nr;
            outfile.open(txtoutput);
            for(const auto& k: sumhist) outfile << k << " ";
            outfile << endl;
          }
        }
      }
      else if (output=="area"){
        ofstream surf, part, clust, border;
        surf.open("./lars_sim/Data/surf/tri_sv_low"+occ_p+".txt");
        part.open("./lars_sim/Data/surf/tri_part_low"+occ_p+".txt");
        clust.open("./lars_sim/Data/surf/tri_clust_low"+occ_p+".txt");
        border.open("./lars_sim/Data/surf/tri_border_low"+occ_p+".txt");
        // unsigned so only one cluster each gets printed
        unsigned check_1 = 0;
        unsigned check_2 = 0;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          vec a = TL.surface_volume();
          for (const auto& k : a) surf << k << " ";
          surf << endl;
          part << TriangleParticleWriter(TL, part) << endl;
          for (const auto& n : TL.cluster_surface())border << n << " ";
          border << endl;

          if (TL.clust_size(8, 10) == 1 && check_1 == 0){
            vec n = TL.single_cluster(8, 10);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_1++;
          }

          if (TL.clust_size(3400, 3500) == 1 && check_2 == 0){
            vec n = TL.single_cluster(3400, 3500);
            for (const auto& m : n) clust << m << " ";
            clust << endl;
            check_2++;
          }

        }
      }
      else if (output == "weighted"){
         ofstream output, part;
         part.open("./lars_sim/Data/weighted/tri_part.txt");
         output.open("./lars_sim/Data/weighted/tri.txt");
         
         for(double n=1; t < burnin + until; n++) {
            t = TL.run_until(burnin + n * every);
            part << TriangleParticleWriter(TL, part) << endl;
            hist_t hist = TL.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            double weighted = second_moment / first_moment;
            weighted = weighted / double(P.N);
            output << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
         }
      }
      else if (output == "heatmap"){
        ofstream outfile, outfile_avg;
        outfile.open("./lars_sim/heatmap/tri_alpha_N.txt");
        outfile_avg.open("./lars_sim/heatmap/tri_alpha_N_avg.txt");
        for (double alp = 1e-6; alp <= 1.0; alp*=1.2){
          for (unsigned N = 100; N <= P.L[0]*P.L[0]; N += 200){
            Parameters P_h;
            P_h.N = N;
            P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
            P_h.L = P.L;
            P_h.n_max = P.n_max;
            std::size_t maxsize = 1;
            Triangle_lattice LB(P_h, rng);
            t = 0;

            // also taking mean over all used timesteps
            double mean = 0.0;
            double count = 0;
            for(unsigned n=0; t < burnin + until; ++n){
              t = LB.run_until(burnin + n * every);
              maxsize = std::max(maxsize, LB.max_cluster_size());
              mean += LB.avg_cluster_size();
              count += 1;
            }

            mean = mean / count;

            outfile << maxsize << " ";
            outfile_avg << mean << " ";
          }
          outfile << endl;
          outfile_avg << endl;
        } 
      }
      else if (output == "number"){
        ofstream outfile;
        string name = "./lars_sim/number/tri_";
        string name_number = name+"number"+"_"+tumb+"_"+dens+"_"+size+txt;
        outfile.open(name_number);


        t = TL.run_until(burnin);
        hist_t hist = TL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i+=2){
          second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
          first_moment += hist[i] * ((i+2)/2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        outfile << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
        for(double n=1; t < burnin + until; n*=1.01) {
          t = TL.run_until(burnin + n * every);
          hist_t hist = TL.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i+=2){
            second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
            first_moment += hist[i] * ((i+2)/2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          outfile << t << " " << weighted << " " << TL.avg_cluster_size_nr() << endl;
        }
      }
      else if (output == "stable"){
        ofstream part, numb;
        part.open("./lars_sim/Data/stable/triangular.txt");
        numb.open("./lars_sim/Data/stable/triangular_number.txt");
        
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          // only doing a positional output here
          numb << t << " " << TL.number_cluster() << endl;
          part << TriangleParticleWriter(TL, part) << endl;
        }
      }else if (output == "motility"){
        ofstream outfile;
        outfile.open("./lars_sim/Data/motility/triangular_perc_low.txt");
        for (double al = 0.0; al < 0.2 ; al+=0.005){
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
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
      
      }
      else if (output == "distribution"){
        ofstream outfile, outfile2;
        outfile.open("./lars_sim/Data/dist/tri_"+occ_p+".txt");
        outfile2.open("./lars_sim/Data/dist/tri_dens_"+occ_p+".txt");
        hist_t dist(7);
        for(double n = 0; t < burnin + until; n++) {
          
          t = TL.run_until(burnin + n * every);
          vec_d dens = TL.density();
          for (const auto& y : dens) outfile2 << y << " ";
          outfile2 << endl;
          
          hist_t dr = TL.particle_neighbour_dist();
          dist[std::slice(0,dr.size(),1)] += dr;
          
        }
        
        for (const auto& m : dist) outfile << m << " ";
        outfile << endl;
      }
      
      else if (output == "stopping time"){
        ofstream outfile;
        
        outfile.open("./lars_sim/Data/stopping/tri_"+occ_p+".txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = TL.run_until(burnin + n * every);
          vec_d st = TL.stopping(t);
          for (const auto& m : st) outfile << m << " ";
          outfile << endl;
        }
      }
      else {
        ofstream outfile;
        outfile.open("./lars_sim/gif/triangle.txt");
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
      std::cout << "# occupation number = " << P.n_max << std::endl;

      double t = 0;

      if(output == "clusters") {
        if(localaverage == 0) {
          hist_t sumhist;
          // We sum the histograms over all measurements
          //ofstream outfile_part;
          //outfile_part.open("./lars_sim/gif/hexagonal.txt");

          for(unsigned n=0; t < burnin + until; ++n) {
            t = HL.run_until(burnin + n * every);
            hist_t hist = HL.cluster_distributions();
            // particle output for comparison
            //outfile_part << HexagonalParticleWriter(HL, outfile_part) << endl;
            // Only area
            if(hist.size() > sumhist.size()) {
              hist[std::slice(0,sumhist.size(),1)] += sumhist;
              sumhist = std::move(hist);
            } else {
              sumhist[std::slice(0,hist.size(),1)] += hist;
            }

          }
          ofstream outfile;
          outfile.open(txtoutput);
          for(const auto& k: sumhist) outfile << k << " ";
          outfile << endl;
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
            ofstream outfile;
            outfile.open("./lars_sim/hexdist.txt");
            for(const auto& k: sumhist) outfile << k << " ";
            outfile << endl;
          }
        }
      } 
      else if (output=="area"){
        ofstream surf, part, clust, border;
        surf.open("./lars_sim/Data/surf/hex_sv_low"+occ_p+".txt");
        part.open("./lars_sim/Data/surf/hex_part_low"+occ_p+".txt");
        clust.open("./lars_sim/Data/surf/hex_clust_low"+occ_p+".txt");
        border.open("./lars_sim/Data/surf/hex_border_low"+occ_p+".txt");
        // unsigned so only one cluster each gets printed
        unsigned check_1 = 0;
        unsigned check_2 = 0;
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          vec a = HL.surface_volume();
          for (const auto& k : a) surf << k << " ";
          surf << endl;
          part << HexagonalParticleWriter(HL, part) << endl;
          for (const auto& n : HL.cluster_surface()) border << n << " ";
          border << endl;

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

        }
      }
      else if (output == "weighted"){
         ofstream output, part;
         part.open("./lars_sim/Data/weighted/hex_part.txt");
         output.open("./lars_sim/Data/weighted/hex.txt");
         
         for(double n=1; t < burnin + until; n++) {
            t = HL.run_until(burnin + n * every);
            part << HexagonalParticleWriter(HL, part) << endl;
            hist_t hist = HL.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            double weighted = second_moment / first_moment;
            weighted = weighted / double(P.N);
            output << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
         }
      }
      else if (output == "stable"){
        ofstream part, numb;
        part.open("./lars_sim/Data/stable/hexagonal.txt");
        numb.open("./lars_sim/Data/stable/hexagoonal_number.txt");
        
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          // only doing a positional output here
          numb << t << " " << HL.number_cluster() << endl;
          part << HexagonalParticleWriter(HL, part) << endl;
        }
      }
      else if (output == "heatmap"){
        ofstream outfile, outfile_avg;
        outfile.open("./lars_sim/heatmap/hex_alpha_N.txt");
        outfile_avg.open("./lars_sim/heatmap/hex_alpha_N_avg.txt");
        for (double alp = 1e-6; alp <= 1.0; alp*=1.2){
          for (unsigned N = 200; N <= 2*P.L[0]*P.L[0]; N += 400){
            Parameters P_h;
            P_h.N = N;
            P_h.alpha[0] = P.alpha[1] = P.alpha[2] = alp;
            P_h.L = P.L;
            P_h.n_max = P.n_max;
            std::size_t maxsize = 1;
            Hexagonal_lattice LB(P_h, rng);
            t = 0;

            // also taking mean over all used timesteps
            double mean = 0.0;
            double count = 0;
            for(unsigned n=0; t < burnin + until; ++n){
              t = LB.run_until(burnin + n * every);
              maxsize = std::max(maxsize, LB.max_cluster_size());
              mean += LB.avg_cluster_size();
              count += 1;
            }

            mean = mean / count;

            outfile << maxsize << " ";
            outfile_avg << mean << " ";
          }
          outfile << endl;
          outfile_avg << endl;
        }
      
      }
      
      else if (output == "number"){
        ofstream outfile;
        string name = "./lars_sim/number/hex_";
        string name_number = name+"number"+"_"+tumb+"_"+dens+"_"+size+txt;
        outfile.open(name_number);


        t = HL.run_until(burnin);
        hist_t hist = HL.cluster_distributions();
        double second_moment = 0;
        double first_moment = 0;
        for (unsigned i = 0; i < hist.size(); i+=2){
          second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
          first_moment += hist[i] * ((i+2)/2);
        }
        double weighted = second_moment / first_moment;
        weighted = weighted / double(P.N);
        outfile << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
        for(double n=1; t < burnin + until; n*=1.01) {
          t = HL.run_until(burnin + n * every);
          hist_t hist = HL.cluster_distributions();
          double second_moment = 0;
          double first_moment = 0;
          for (unsigned i = 0; i < hist.size(); i+=2){
            second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
            first_moment += hist[i] * ((i+2)/2);
          }
          double weighted = second_moment / first_moment;
          weighted = weighted / double(P.N);
          outfile << t << " " << weighted << " " << HL.avg_cluster_size_nr() << endl;
        }
      }else if (output == "particles") {
        ofstream outfile;
        outfile.open("./lars_sim/gif/hexagonal.txt");

        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          // only doing a positional output here
          
          outfile << HexagonalParticleWriter(HL, outfile) << endl;

        }

      }else if (output == "motility"){
        ofstream outfile;
        outfile.open("./lars_sim/Data/motility/hexagonal_perc_low.txt");
        for (double al = 0.0; al < 0.2 ; al+=0.005){
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
          for(unsigned n=0; t < burnin + until; ++n) {
            t = LB.run_until(burnin + n * every);
            values_mas.push_back(double(LB.max_cluster_size_nr())/double(P.N));
            values_mot.push_back(LB.motility_fraction());
            rel_mass += double(LB.max_cluster_size_nr())/double(P.N);
            mean += LB.motility_fraction();
            count++;
            //outfile << t << " " << HL.motility_fraction() << " " << rel_mass << endl;

            hist_t hist = LB.cluster_distributions();
            double second_moment = 0;
            double first_moment = 0;
            for (unsigned i = 0; i < hist.size(); i+=2){
              second_moment += hist[i] * ((i+2)/2) * ((i+2)/2);
              first_moment += hist[i] * ((i+2)/2);
            }
            //std::cout << second_moment << " " << first_moment << endl;
            values_wei.push_back(second_moment / first_moment * 1.0/double(P.N));
            weighted += second_moment / first_moment * 1.0/double(P.N);
          }
          mean = mean / count;
          rel_mass = rel_mass / count;
          weighted = weighted / count;


          double cov_mot = 0;
          for (auto& val : values_mot){
            cov_mot += pow(val - mean, 2);
          }
          cov_mot = cov_mot/(values_mot.size() - 1);

          double cov_mas = 0;
          for (auto& val : values_mas){
            cov_mas += pow(val - rel_mass, 2);
          }
          cov_mas = cov_mas/(values_mas.size() - 1);

          for (auto& val : values_wei){
            cov_w += pow(val - weighted, 2);
          }
          cov_w = cov_w/(values_wei.size() - 1);

          outfile << al << " " << mean << " " << cov_mot << " " << rel_mass << " " << cov_mas << " " << weighted << " " << cov_w << endl;
        }
      
      } else if (output == "snapshots"){
        ofstream outfile;
        outfile.open("hexdir.txt");

        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          // only doing a positional output here
          HL.realise_directions();
          outfile << HexDirectionWriter(HL, outfile) << endl;
        }

      }
      else if (output == "function"){
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
        }
      }
      else if (output == "distribution"){
        ofstream outfile, outfile2;
        outfile.open("./lars_sim/Data/dist/hex_"+occ_p+".txt");
        outfile2.open("./lars_sim/Data/dist/hex_dens_"+occ_p+".txt");
        hist_t dist(4);
        for(double n = 0; t < burnin + until; n++) {
          t = HL.run_until(burnin + n * every);
          vec_d dens = HL.density();
          for (const auto& m : dens) outfile2 << m << " ";
          outfile2 << endl;
          hist_t dr = HL.particle_neighbour_dist();
          dist[std::slice(0,dr.size(),1)] += dr;
        }
        for (const auto& m : dist) outfile << m << " ";
        outfile << endl;
      }
      else if (output == "stopping time"){
        ofstream outfile;
        
        outfile.open("./lars_sim/Data/stopping/hex_"+occ_p+".txt");
        for(unsigned n=0; t < burnin + until; ++n) {
          t = HL.run_until(burnin + n * every);
          vec_d st = HL.stopping(t);
          for (const auto& m : st) outfile << m << " ";
          outfile << endl;
        }
      }

    }
  }
  // maybe do other lattices as well ?
  return 0;
}
