#!/usr/bin/env python2.7

import matplotlib.pyplot as plt
from ascii_graph import Pyasciigraph
import collections
import numpy as np
import itertools
import random
import sys
import os

from species import species

# pylint: disable=C0103
# pylint: disable=R0903

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1000

class Community(object):

    def __init__(self, K=5000, colrate=0.01, mig_clust_size=1, quiet=False):
        self.quiet = quiet

        ## List for storing species objects that have had sequence
        ## simulated and sumstats calculated
        self.species_objects = []

        ## Settings specific to the uniform metacommunity
        ## This is individuals per species
        self.uniform_inds = 10000000
        self.uniform_species = 1000

        ## Variables associated with the metacommunity (these are poorly named)
        ## total_inds has to be set by set_metacommunity
        ## This is total individuals in the metacommunity
        self.total_inds = 0

        ## Metacommunity will be a 2-d array
        ## 1d: uuid
        ## 2d: metacommunity abundance
        ## 3d: immigration probability
        self.metacommunity = np.array([])

        self.maxabundance = 0
        self.colonization_rate = colrate

        ## Set the size of the cluster of colonizers
        ## Setting this to 1 equates to the basic immigration model
        self.mig_clust_size=mig_clust_size

        ## Variables for tracking the local community
        ## Dimensions are: species uuid identifier, local abundance, colonization time
        self.local_community = np.array([])
        self.local_abundances = []
        self.K = K

        self.extinctions = 0
        self.colonizations = 0
        self.current_time = 0


    def set_metacommunity(self, infile):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there
        """
        dt = np.dtype([('uuid', int), ('abund', np.int64), ('col_prob', float)])
        if infile == "uniform":
            #self.metacommunity = np.array([self.uniform_inds] * self.uniform_species, dtype=np.dtype("int"))
            #abunds = np.recarray([self.uniform_inds] * self.uniform_species, dtype=dt)
            abunds = np.array([self.uniform_inds] * self.uniform_species, dtype=np.int64)
        else:
            if os.path.isfile(infile):
                with open(infile, 'r') as inf:
                    abunds = np.array([int(line.split()[0]) for line in inf], dtype=np.int64)
            else:
                raise Exception("Bad metacommunity input - ".format(infile))

        #self.metacommunity = np.recarray(len(abunds), dtype=dt)
        self.metacommunity = np.zeros(len(abunds), dtype=dt)
        self.metacommunity["uuid"] = np.arange(0, len(abunds))
        self.metacommunity["abund"] = abunds
        self.total_inds = sum(self.metacommunity["abund"])
        print("total_inds", self.total_inds)
        self.maxabundance = np.amax(self.metacommunity["abund"])
        self.metacommunity["col_prob"] = self.metacommunity["abund"]/float(self.total_inds)
        print("metacommunity", self.metacommunity)
        """
        self.metacommunity = np.vstack([np.arange(0, len(self.metacommunity)), self.metacommunity])
        self.total_inds = sum(self.metacommunity[1])
        self.metacommunity = np.vstack([self.metacommunity, np.array(self.metacommunity[1]/float(self.total_inds))])
        self.maxabundance = np.amax(self.metacommunity[1])
        self.local_community = np.zeros((4, len(self.metacommunity[0])))
        self.local_community[0] = self.metacommunity[0]
        print(self.local_community)
        """

    def __str__(self):
        return "<Community {}>".format(self.name)


    def prepopulate(self, mode="landbridge"):
        dt = np.dtype([('uuid', int), ('abund', int), ('col_time', int), ('is_founder', np.bool)])
        #self.local_community = np.recarray(len(self.metacommunity), dtype=dt)
        self.local_community = np.zeros(len(self.metacommunity), dtype=dt)
        self.local_community["uuid"] = np.arange(0, len(self.metacommunity))
        self.local_community["col_time"] = 0
        self.local_community["abund"] = 0
        self.local_community["is_founder"] = False

        if mode == "landbridge":
            ## prepopulate the island w/ total_inds individuals sampled from the metacommunity
            print(self.metacommunity.dtype)
            init_community = np.random.multinomial(self.K, self.metacommunity['col_prob'], size=1)        
            print("init_community", init_community)
            print("Initializing local community:")
            for i, x in enumerate(init_community[0]):
                if x:
                    self.local_community['abund'][i] += x
                    self.local_community['col_time'][i] = 1
                    self.local_community['is_founder'][i] = True
                    self.local_abundances.extend([i] * x)
            print("N species = {}".format(np.where(self.local_community['abund'] != 0)[0].size))
            print("N individuals = {}".format(np.sum(self.local_community['abund'])))


        else:
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            new_species = np.where(self.metacommunity['abund'] == self.maxabundance)[0][0]
            self.local_community['abund'][new_species] = self.K
            self.local_community['col_time'][new_species] = 1
            self.local_community['is_founder'][new_species] = True
            self.local_abundances = [new_species] * self.K
        print("local_community", self.local_community)

    def remove_from_local(self, n=1):
        for i in xrange(n):

            idx = self.get_random_individual()
            #if self.local_community['abund'][idx] == 0:
            #    sys.exit("got a bad idx")
            self.local_community['abund'][idx] -= 1
            ## Record local extinction events
            if self.local_community['abund'][idx] == 0:
                self.extinctions += 1
            ## Don't need this if array works
            try:
                self.local_abundances.remove(idx)
            except Exception as inst:
                print(inst)
                raise


    def get_random_individual(self):
        ## Select the individual to die
        ## SLOWWWWWW
        # victim = np.random.multinomial(1, self.local_community["abund"]/float(self.K), size=1)
        # idx = self.local_community['uuid'][victim[0].astype("bool")]

        return random.choice(self.local_abundances)


    def step(self, nsteps=1):
        for step in xrange(nsteps):

            ## Check probability of an immigration event
            if np.random.random_sample() < self.colonization_rate:
                ## Loop until you draw species unique in the local community
                ## The flag to tell 'while when we're done, set when you successfully 
                ## draw a non-local-doop from the metacommunity
                unique = 0
    
                ## If you set your carrying capacity too high relative to the size of your
                ## metacommunity then you'll get stuck drawing duplicates over and over
                idiot_count = 0
                while not unique:
                    ## Sample from the metacommunity
                    migrant_draw = np.random.multinomial(1, self.metacommunity["col_prob"], size=1)[0]
                    new_species = self.metacommunity["uuid"][migrant_draw.astype("bool")]
                    ##TODO: Should set a flag to guard whether or not to allow multiple colonizations
                    if self.local_community["abund"][new_species] > 0:
                        #print(len(np.nonzero(self.local_community["abund"])[0])),
                        #print(np.nonzero(self.local_community["abund"])[0].sum())
                        #print(new_species, self.local_community["abund"][new_species])
                        #print("multiple colonization events are forbidden, for now")
                        unique = 0
    
                        if idiot_count > MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY:
                            msg = """
    Metacommunity is exhausted w/ respect to local
    community. Either expand the size of the metacommunity,
    decrease the carrying capacity, or switch on multiple
    migration (unimplemented)."""
                            sys.exit(msg)
                        idiot_count +=1
                    else:
                        self.remove_from_local(n=self.mig_clust_size)
                        self.local_community["abund"][new_species] += self.mig_clust_size
                        self.local_community["col_time"][new_species] = self.current_time
                        self.local_community["is_founder"][new_species] = False
                        self.colonizations += 1

#                       Don't need this if array works
                        self.local_abundances.append(new_species)
                        unique = 1
            else:
                ## Not a colonization event, remove 1 from local community
                self.remove_from_local(n=1)

                ## Sample from the local community, including empty demes
                ## Sample all available from local community (community grows slow in volcanic model)
                ## Also, lots of early turnover
                ## This all was only true before i implemented the full rosindell/harmon model,
                ## There are no more empty demes in the current config
                
                idx = self.get_random_individual()
                self.local_community['abund'][idx] += 1

                ## Don't need this if array works
                self.local_abundances.append(idx)

                ## Sample only from available extant species (early pops grow quickly in the volcanic model)
                ## If you do this, the original colonizer just overwhelms everything else
                ## This is more similar to the Rosindell and Harmon model, in which they simply
                ## prepopulate the island entirely with one species. This is effectively the same
                #self.local_community.append(random.choice([x for x in self.local_community if not x[0] == 0]))
    
            ## update current time
            self.current_time += 1

    def get_abundances(self, octaves=False):
        ## Get abundance classes and counts all in one swoop
        classes_and_counts = np.unique(self.local_community["abund"], return_counts=True)

        ## Conver to a dict cuz that's what downstream expects
        abundance_distribution = collections.OrderedDict(zip(classes_and_counts[0], classes_and_counts[1]))
        
        ## This is another way to do it
        ## abundances = self.local_community[np.nonzero(self.local_community["abund"])]
        ## classes = np.bincount(abundances["abund"])
        ## counts = np.nonzero(classes)[0]
        ## abundance_distribution = zip(classes, classes[counts])

        if octaves:
            dist_in_octaves = collections.OrderedDict()
            min = 1
            max = 2
            while max/2 < len(abundance_distribution):
                ## Count up all species w/in each octave
                count = 0
                ## Here `i` is the abundance class and
                ## `j` is the count for that class
                for i, j in abundance_distribution.items():
                    if (i < max) and (i >= min):
                        count += j
                dist_in_octaves[min] = count
                min = min * 2
                max = max * 2
            abundance_distribution = dist_in_octaves
        return abundance_distribution


    def simulate_seqs(self):
        self.species_objects = []
        ## Setting colonization_time as a scaling factor rather than as a raw tdiv
        ## The old way of doing this is `self.current_time - tdiv`
        #self.species_objects = [species(UUID=UUID, colonization_time=1/float(tdiv), abundance=self.local_community.count(UUID),\
        for UUID, tdiv in self.divergence_times.items():
            #print(self.local_community)
            if UUID in self.local_community:
                meta_abundance = -1
                for x, y in self.species:
                    if UUID[0] == x:
                        meta_abundance = y
                #meta_abundance = [x[1] for x in self.abundances if x[0] == UUID[0]]
                #meta_abundance = self.abundances[self.species.index(UUID[0])]
                abundance = self.local_community.count(UUID)
                #print(self.local_community)
                self.species_objects.append(species(UUID=UUID, colonization_time=self.current_time - tdiv,\
                                        abundance=abundance,\
                                        meta_abundance=meta_abundance))

        for s in self.species_objects:
            s.simulate_seqs()
            s.get_sumstats()

    def shannon_idex(self):
        ## Get local abundances
        local_abunds = self.local_community[np.nonzero(self.local_community["abund"])]
        ## Make raw counts into proportions
        local_abunds = local_abunds/float(self.K)
        ## Calculate shannon index and return
        local_abunds = local_abunds * np.log(local_abunds)
        return local_abunds.sum() * -1

    def get_species(self):
        return(self.species_objects)


if __name__ == "__main__":
    data = Community()
    #data.set_metacommunity("uniform")
    data.set_metacommunity("metacommunity_LS4.txt")
    #data.prepopulate(mode="landbridge")
    data.prepopulate(mode="volcanic")
    for i in range(100000):
        if not i % 10000:
            print("Done {}".format(i))
            #print(i, np.sum(data.local_community["abund"]), np.sum(data.local_community["abund"].astype("bool")))
            #print(data.local_community)
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    print("Species abundance distribution:\n{}".format(abundance_distribution))
    #print("Colonization times per species:\n{}".format(data.local_community[np.nonzero(data.local_community["col_time"])]))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()
    print("Species:\n{}".format(data.get_species()))
    print("Extinction rate - {}".format(data.extinctions/float(data.current_time)))
    print("Colonization rate - {}".format(data.colonizations/float(data.current_time)))
