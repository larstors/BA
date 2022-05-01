"""
First draft of MIPS on different kinds of lattice using Gillespie

Credit goes to Richard Blythe
"""


from sys import maxsize
import matplotlib.pyplot as plt
import numpy as np
import random
import scipy as sci
import time




"""
class Parameters:
    def __init__(self, L=[100], N=50, alpha=[1.0]):
        self.L = L   # dimensions of simulation
        self.N = N      # Number of particles
        self.alpha = alpha # tumble rate
    
    

class Site(Parameters):
    def __init__(self, L, N, alpha):
        super().__init__(L, N, alpha)
        self.id = 0                     # particle id
        self.occ = 0                    # true if particle, else false
        self.act = 0                    # true if event scheduled, else false
        self.neighbours = 0             # number of neighbours. Max number depends on lattice
        self.direction = 0              # Direction of last hop
        self.hoptime = 0                # time of last hop

class Scheduler(Site):
    def __init__(self, _time, L, N, alpha):
        super().__init__(L, N, alpha)
        self._time = _time
        self.queue = []
    
    def schedule(self, delay, event, position):
        self.queue.append([delay+self._time, event, position])                  # Has structure (time, func, lattice position). It gets sorted after time.
        self.queue.sort(key=lambda x: x[0])

    def advance(self, tmax):
        if tmax < self._time:
            return False
        
        if not self.queue:
            self._time = tmax
            return False
        
        te = self.queue[0]

        if te[0] > self._time:
            self._time = tmax
            return False
        
        self.queue.pop(0)
        self._time = te[0]
        te[1](te[2])                     # ! See if this requires more arguments and if it even works
        return True

    def pending(self):
        return len(self.queue)

"""

class site:
    def __init__(self, check, i):
        self.id = i                     # particle id
        self.occ = check                    # true if particle, else false
        self.act = False                    # true if event scheduled, else false
        self.neighbours = 0             # number of neighbours. Max number depends on lattice
        self.direction = 0              # Direction of last hop
        self.hoptime = 0                # time of last hop

class Square_Lattice(site):
    def __init__(self, _time, L, N, alpha, check=False) -> None:
        super().__init__(check, 1)
        # super().__init__(_time, L, N, alpha)
        self.L = L   # dimensions of simulation
        self.N = N      # Number of particles
        self.alpha = alpha # tumble rate
        self._time = _time
        self.queue = []
        self.lattice_size = np.prod(self.L) 
        self.sites = []
        self.piu()
        self.Fill_lattice()
        
    
    def piu(self):
        """filling list of lattice sites
        """
        for i in range(self.lattice_size):
            self.sites.append(site(False, i))

    def schedule(self, delay, event, position):
        """schedule an event and sort eventlist by time

        Args:
            delay (_type_): how far in future is event
            event (_type_): what is the event
            position (_type_): what site is it
        """
        self.queue.append([delay+self._time, event, position])                  # Has structure (time, func, lattice position). It gets sorted after time.
        self.queue.sort(key=lambda x: x[0])

    def advance(self, tmax):
        if tmax < self._time:
            return False
        
        if not self.queue:
            self._time = tmax
            return False
        
        te = self.queue[0]

        if te[0] > self._time:
            self._time = tmax
            return False
        
        self.queue.pop(0)
        self._time = te[0]
        te[1](te[2])                     # ! See if this requires more arguments and if it even works
        return True

    def pending(self):
        return len(self.queue)
    
    def place(self, n, id, direction, time):
        self.sites[n].id = id
        self.sites[n].direction = direction
        self.sites[n].hoptime = time
        if not self.sites[n].occ:
            self.sites[n].occ = True
            for m in self.neighbours_index(n):
                self.sites[m].neighbours += 1
        return



    def neighbours_index(self, n):
        # array with indices
        index = np.zeros(2*len(self.L))
        #parameter used in following calculation
        below = 1
        for d in range(len(self.L)):
            # size of current dimesion
            L = self.L[d]
            # volume of the d-subspace
            above = below * L

            # some magic I have yet to fully understand
            x = int(n / below) % L
            y = n % below + int(n/above) * above

            # the foward and backward indices to neighbours
            index[2*d] = y + ((x + 1) % L) * below
            index[2*d + 1] = y + ((x + L - 1) % L) * below

            #update below
            below = above

        return index.astype(int)

    def neighbours_index_forward(self, n):
        # array with indices
        index = np.zeros(len(self.L))
        #parameter used in following calculation
        below = 1
        for d in range(len(self.L)):
            # size of current dimesion
            L = self.L[d]
            # volume of the d-subspace
            above = below * L

            # some magic I have yet to fully understand
            x = int(n / below) % L
            y = n % below + int(n/above) * above

            # the foward and backward indices to neighbours
            index[d] = y + ((x + 1) % L) * below

            #update below
            below = above

        return index.astype(int)


    def direct(self):
        """ Random Direction generator

        Returns:
            _type_: _description_
        """
        choice = np.array([i for i in range(len(self.alpha))])
        return np.random.choice(choice, 1, self.alpha)[0]
    
    def Fill_lattice(self):
        """Generates filled lattice
        """
        unplaced = self.N
        for n in range(self.lattice_size):
            if np.random.randint(1, self.lattice_size - n+1) <= unplaced:
                self.place(n, self.N - unplaced, self.direct(), 0.0)
                unplaced -= 1
            
            else:
                self.sites[n].id = unplaced + n


        return 
    






    def sched(self, n):
        self.schedule(np.random.exponential(scale=1/np.sum(self.alpha)), self.move, n)
    
    def move(self, n):
        assert(self.sites[n].occ)            # check if occupied
        assert(self.sites[n].act)            # check if active
        if self.sites[n].neighbours == 2 * len(self.L):
            self.sites[n].act = False
        else:
            if np.random.exponential(scale=1/np.sum(self.alpha)) < self._time - self.sites[n].hoptime:
                self.sites[n].direction = self.direct()
            self.sites[n].hoptime = self._time
            dnbs = self.neighbours_index(n)
            if not self.sites[dnbs[self.sites[n].direction]].occ:
                # ! I am not doing the assert  here, cause it seems to be redundant
                # interchanging n with target site
                vid = self.sites[dnbs[self.sites[n].direction]].id
                self.sites[n].occ = False
                self.sites[n].act = False
                self.place(dnbs[self.sites[n].direction], self.sites[n].id, self.sites[n].direction, self.sites[n].hoptime)
                self.sites[n].id = vid
                for m in dnbs:
                    self.sites[m] -= 1
                    if self.sites[m].occ and not self.sites[m].act:
                        self.sched(m)
            else:
                self.sched(n)
        self.sites[n].act = True

    def run_until(self, time):
        while self.advance(time):
            # ! add check for consistency
            a = 1                                   # just a line to stop the editor from showing something is wrong
        return self._time

    
    def cluster_distribution(self):
        memberof = np.array([i for i in range(self.lattice_size)])
        clusters = []
        for n in range(self.lattice_size):
            clusters.append([n])

        maxsize = 1

        for n in range(self.lattice_size):
            
            for m in self.neighbours_index_forward(n):
                large = memberof[n]
                small = memberof[m]
                # merge clusters
                if (self.sites[n].occ == self.sites[m].occ and small != large):
                    # ensure that large actually is larger, if not make it so
                    #if len(clusters[large]) < len(clusters[small]):
                    #    large, small = small, large
                    
                    
                    for k in clusters[small]:
                        memberof[k] = large

                    # add small to large
                    for k in clusters[small]:
                        clusters[large].append(k)

                    #delete small
                    clusters[small] = []


                    maxsize = max(maxsize, len(clusters[large]))
                    

        
        
        # array with hist sizes
        dists = np.zeros(2 * maxsize)
        # fill it
        for k in range(len(clusters)):
            if len(clusters[k]) != 0:
                dists[2 * len(clusters[k]) - 1 - self.sites[clusters[k][0]].occ] += 1
        
        
        return dists
        

class Triangular_Lattice(site):
    def __init__(self, _time, L, N, alpha, check=False) -> None:
        super().__init__(check, 1)
        # super().__init__(_time, L, N, alpha)
        self.L = L   # dimensions of simulation
        self.N = N      # Number of particles
        self.alpha = alpha # tumble rate ! in 2d it should be 3d
        self._time = _time
        self.queue = []
        self.lattice_size = np.prod(self.L) 
        self.sites = []
        self.piu()
        self.Fill_lattice()
        
    
    def piu(self):
        """filling list of lattice sites
        """
        for i in range(self.lattice_size):
            self.sites.append(site(False, i))

    def schedule(self, delay, event, position):
        """schedule an event and sort eventlist by time

        Args:
            delay (_type_): how far in future is event
            event (_type_): what is the event
            position (_type_): what site is it
        """
        self.queue.append([delay+self._time, event, position])                  # Has structure (time, func, lattice position). It gets sorted after time.
        self.queue.sort(key=lambda x: x[0])

    def advance(self, tmax):
        if tmax < self._time:
            return False
        
        if not self.queue:
            self._time = tmax
            return False
        
        te = self.queue[0]

        if te[0] > self._time:
            self._time = tmax
            return False
        
        self.queue.pop(0)
        self._time = te[0]
        te[1](te[2])                     # ! See if this requires more arguments and if it even works
        return True

    def pending(self):
        return len(self.queue)
    
    def place(self, n, id, direction, time):
        self.sites[n].id = id
        self.sites[n].direction = direction
        self.sites[n].hoptime = time
        if not self.sites[n].occ:
            self.sites[n].occ = True
            for m in self.neighbours_index(n):
                self.sites[m].neighbours += 1
        return



    def neighbours_index(self, n):
        # ! so far this is only valid for 2d
        # choosing to go "up right" and "down left" for the diagonal
        # array with indices
        index = np.zeros(2*len(self.L) + 2) 
        #parameter used in following calculation
        
        x = n % self.L[0]
        y = int(n/self.L[0])

        index[0] = (n + 1) % self.L[0]
        index[1] = (n - 1) % self.L[0]

        index[2] = x + ((y + 1) % self.L[1]) * self.L[0]
        index[3] = x + ((y - 1) % self.L[1]) * self.L[0]

        index[4] = (n + 1) % self.L[0] + ((y + 1) % self.L[1]) * self.L[0]
        index[5] = (n - 1) % self.L[0] + ((y - 1) % self.L[1]) * self.L[0]

        return index.astype(int)

    def neighbours_index_forward(self, n):
        # array with indices
        index = np.zeros(len(self.L) + 1) 
        #parameter used in following calculation
        
        x = n % self.L[0]
        y = int(n/self.L[0])

        index[0] = (n + 1) % self.L[0]

        index[1] = x + ((y + 1) % self.L[1]) * self.L[0]

        index[2] = (n + 1) % self.L[0] + ((y + 1) % self.L[1]) * self.L[0]

        return index.astype(int)


    def direct(self):
        """ Random Direction generator

        Returns:
            _type_: _description_
        """
        choice = np.array([i for i in range(len(self.alpha))])
        return np.random.choice(choice, 1, self.alpha)[0]
    
    def Fill_lattice(self):
        """Generates filled lattice
        """
        unplaced = self.N
        for n in range(self.lattice_size):
            if np.random.randint(1, self.lattice_size - n+1) <= unplaced:
                self.place(n, self.N - unplaced, self.direct(), 0.0)
                unplaced -= 1
            
            else:
                self.sites[n].id = unplaced + n


        return 
    






    def sched(self, n):
        self.schedule(np.random.exponential(scale=1/np.sum(self.alpha)), self.move, n)
    
    def move(self, n):
        assert(self.sites[n].occ)            # check if occupied
        assert(self.sites[n].act)            # check if active
        if self.sites[n].neighbours == 2 * len(self.L):
            self.sites[n].act = False
        else:
            if np.random.exponential(scale=1/np.sum(self.alpha)) < self._time - self.sites[n].hoptime:
                self.sites[n].direction = self.direct()
            self.sites[n].hoptime = self._time
            dnbs = self.neighbours_index(n)
            if not self.sites[dnbs[self.sites[n].direction]].occ:
                # ! I am not doing the assert  here, cause it seems to be redundant
                # interchanging n with target site
                vid = self.sites[dnbs[self.sites[n].direction]].id
                self.sites[n].occ = False
                self.sites[n].act = False
                self.place(dnbs[self.sites[n].direction], self.sites[n].id, self.sites[n].direction, self.sites[n].hoptime)
                self.sites[n].id = vid
                for m in dnbs:
                    self.sites[m] -= 1
                    if self.sites[m].occ and not self.sites[m].act:
                        self.sched(m)
            else:
                self.sched(n)
        self.sites[n].act = True

    def run_until(self, time):
        while self.advance(time):
            # ! add check for consistency
            a = 1                                   # just a line to stop the editor from showing something is wrong
        return self._time

    
    def cluster_distribution(self):
        memberof = np.array([i for i in range(self.lattice_size)])
        clusters = []
        for n in range(self.lattice_size):
            clusters.append([n])

        maxsize = 1

        for n in range(self.lattice_size):
            
            for m in self.neighbours_index_forward(n):
                large = memberof[n]
                small = memberof[m]
                # merge clusters
                if (self.sites[n].occ == self.sites[m].occ and small != large):
                    # ensure that large actually is larger, if not make it so
                    #if len(clusters[large]) < len(clusters[small]):
                    #    large, small = small, large
                    
                    
                    for k in clusters[small]:
                        memberof[k] = large

                    # add small to large
                    for k in clusters[small]:
                        clusters[large].append(k)

                    #delete small
                    clusters[small] = []


                    maxsize = max(maxsize, len(clusters[large]))
                    

        
        
        # array with hist sizes
        dists = np.zeros(2 * maxsize)
        # fill it
        for k in range(len(clusters)):
            if len(clusters[k]) != 0:
                dists[2 * len(clusters[k]) - 1 - self.sites[clusters[k][0]].occ] += 1
        
        
        return dists



class ParticleWriter(Square_Lattice):
    def __init__(self, _time, L, N, alpha) -> None:
        super().__init__(_time, L, N, alpha)
    
    def printer(self):
        for n in range(self.lattice_size):
            if self.sites[n].occ:
                print("%d %d " % (self.sites[n].id, n))







if __name__=="__main__":
    # -----------------------------------------
    # Setting a default so I don't forget it
    # -----------------------------------------

    # P = Parameters([100], 50, [1.0])

    lattice = Square_Lattice(_time = 0, L=[100, 100], N=6000, alpha=[1.5, 0.5])
    burning = 1000
    until = 5000
    every = 2.5
    t = 0
    n = 0
    localaverage = 0
    localinterval = 10

    output = "clusters"

    if output == "clusters":
        if (localaverage == 0):
            sumhist = np.array([])
            while n < burning + until:
                t = lattice.run_until(burning + n * every)
                hist = lattice.cluster_distribution()
                
                # if hist larger, we add sumhist to hist and say that sumhist is hist
                if (len(hist) > len(sumhist)):
                    hist[0:len(sumhist)] += sumhist
                    sumhist = hist
                else:
                    sumhist[0:len(hist)] += hist

                n+=1
            # ! here we need to put in the output of the data

            # where the sumhist != 0 for particles
            ar = np.argwhere(sumhist[::2] !=  0)
            total = np.sum(sumhist[::2])
            freq = sumhist[ar * 2] / total 
                        
            fig = plt.figure()
            plt.loglog(ar, freq, "go")
            plt.xlabel(r"Cluster size $k$")
            plt.ylabel(r"Frequency $P(k)$")
            plt.savefig("run_1.pdf", pdi=200)
            print(ar, freq)
            
        
        else:
            
            while n < burning + until:
                t = lattice.run_until(burning + n * every)
                sumhhist = lattice.cluster_distribution()
                for m in range(localaverage):
                    t = lattice.run_until(burning + n * every + m * localinterval)
                    hist = lattice.cluster_distribution()
                
                    # if hist larger, we add sumhist to hist and say that sumhist is hist
                    if (len(hist) > len(sumhist)):
                        hist[0:len(sumhist)] += sumhist
                        sumhist = hist
                    else:
                        sumhist[0:len(hist)] += hist

                n+=1
            # ! here we need to put in the output of the data




    else:
        while t < burning + until:
            t = lattice.run_until(burning + n * every)
            print(lattice.cluster_distribution())
            #prin = ParticleWriter(lattice._time, lattice.L, lattice.N, lattice.alpha)
            #prin.printer()

            n+=1
    