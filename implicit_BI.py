#!/usr/bin/env python2.7
""" Sample object """

import matplotlib.pyplot as plt
from ascii_graph import Pyasciigraph
from scipy.stats import logser
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
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1500

class implicit_BI(object):
    """ ipyrad Sample object. Links to files associated
    with an individual sample, used to combine samples 
    into Assembly objects."""

    def __init__(self, K=5000, colrate=0.01, exponential=False, quiet=False):
        self.quiet = quiet

        ## List for storing species objects that have had sequence
        ## simulated and sumstats calculated
        self.species_objects = []
        self.exponential = exponential

        ## Settings specific to the uniform metacommunity
        ## This is individuals per species
        self.uniform_inds = 10000000
        self.uniform_species = 1000

        ## Variables associated with the metacommunity (these are poorly named)
        ## total_inds has to be set by set_metacommunity
        ## This is total individuals in the metacommunity
        self.total_inds = 0
        self.immigration_probabilities = []
        self.abundances = []
        self.species = []

        self.maxabundance = 0
        self.colonization_rate = colrate
        self.allow_multiple_colonizations = False

        ## Variables for tracking the local community
        self.local_community = []
        self.local_inds = K
        self.divergence_times = {}
        ## Track number of secondary colonization events per species
        self.post_colonization_migrants = {}

        self.extinctions = 0
        self.colonizations = 0
        self.current_time = 0
        
        ## Vector for tracking lifetimes of excinct species
        ## We can plot this and average it to figure out how long
        ## species hang around in the local community
        self.extinction_times = []

        ## The invasive species identity
        self.invasive = -1
        ## Track how many invasives differentially survived
        self.survived_invasives = 0
        self.invasion_time = 0

        ###########################################################
        ## Variables for non-neutral assembly processes
        ## Note all non-neutral assembly value defaults result
        ## in neutral processes.
        ###########################################################

        ## An ID by trait value dictionary
        self.species_trait_values = {}

        ## Toggle competitive exclusion
        self.competitive_exclusion = False

        ## Toggle environmental filtering
        self.environmental_filtering = False
        self.environmental_optimum = 0


    def set_metacommunity(self, infile, random=False):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there

        random=True will set random trait values in the range [0-1]
        """
        if infile == "logser":
            ## Parameter of the logseries distribution
            p = .98
            self.abundances = logser.rvs(p, size=self.local_inds)
            if random:
                self.species_trait_values = {x:y for x,y in enumerate(np.random.rand(self.uniform_species))}
            else:
                self.species_trait_values = {x:y for x,y in enumerate(np.zeros(len(self.abundances)))}

        elif infile == "uniform":
            #for i in range(self.uniform_inds):
            self.abundances = [self.uniform_inds] * self.uniform_species
            if random:
                self.species_trait_values = {x:y for x,y in enumerate(np.random.rand(self.uniform_species))}
            else:
                self.species_trait_values = {x:y for x,y in enumerate(np.zeros(len(self.abundances)))}
        else:
            if os.path.isfile(infile):
                with open(infile, 'r') as inf:
                    self.abundances = [int(line.split()[0]) for line in inf]
                    ## Try fetching traits as well. If there are no traits then assume we're not doing trait based
                    try:
                        inf.seek()
                        self.species_trait_values = {x:y for x,y in enumerate([float(line.split()[1]) for line in inf])}
                    except:
                        self.species_trait_values = {x:y for x,y in enumerate(np.zeros(len(self.abundances)))}
            else:
                raise Exception("Bad metacommunity input - ".format(infile))

        ## Actually using uuid is v slow
        #self.species = [uuid.uuid4() for _ in enumerate(self.abundances)]
        self.species = [x for x in enumerate(self.abundances)]
        self.total_inds = sum(self.abundances)
        self.immigration_probabilities = [float(self.abundances[i])/self.total_inds for i in range(len(self.abundances))]

        self.maxabundance = np.amax(self.immigration_probabilities)
        

    def __str__(self):
        return "<implicit_BI {}>".format(self.name)


    def prepopulate(self, mode="landbridge"):
        if mode == "landbridge":
            ## prepopulate the island w/ total_inds individuals sampled from the metacommunity
            init_community = np.random.multinomial(self.local_inds, self.immigration_probabilities, size=1)        
            print("Initializing local community:")
            for i, x in enumerate(init_community[0]):
                if x:
                    self.local_community.extend([self.species[i][0]] * x)
                    ## Set founder flag
            self.local_community = [tuple([x, True]) for x in self.local_community]
            print("N species = {}".format(len(set([x[0] for x in self.local_community]))))
            #self.local_community = list(itertools.chain.from_iterable(self.local_community))
            print("N individuals = {}".format(len(self.local_community)))

            ## All species diverge simultaneously upon creation of the island.
            for taxon in self.local_community:
                self.divergence_times[taxon] = 1
        else:
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            ## This is the old way that acts weird
            #self.local_community = [0] * self.local_inds
            #self.local_community = [((x,0), True) for x in self.local_community]
            new_species = (self.species[self.immigration_probabilities.index(self.maxabundance)][0], True)
            self.local_community.append(new_species)
            for i in range(1,self.local_inds):
                self.local_community.append((None,True))
            self.divergence_times[new_species] = 1

    def death_step(self, invasion_time, invasiveness):
        ## Select the individual to die
        victim = random.choice(self.local_community)
        ## If no invasive hasn't invaded then just do the normal sampling
        if self.invasive == -1:
            self.local_community.remove(victim)
        else:
            ## If invasiveness is less than the random value remove the invasive individual
            ## else choose a new individual
            if victim == self.invasive and np.random.rand() < invasiveness:
                self.survived_invasives += 1
                victim = random.choice(self.local_community)
            self.local_community.remove(victim)
        ## Record local extinction events
        if not victim in self.local_community:
            ## This was supposed to not record "extinctions" of empty deme space
            ## but it fucks up the calculation of extinction rate
            ##if not victim[0] == None:
            if True:
                self.extinctions += 1
                try:
                    ## reset 
                    self.post_colonization_migrants[victim[0]] = 0
                    self.extinction_times.append(self.current_time - self.divergence_times[victim])
                except:
                    ## The empty deme will make this freak
                    pass
            ## If the invasive prematurely goes extinct just pick a new one
            if victim == self.invasive:
                print("invasive went extinct")
                self.invasive = -1

    def migrate_no_dupes_step(self):
        ## Loop until you draw species unique in the local community
        ## The flag to tell 'while when we're done, set when you successfully 
        ## draw a non-local-doop from the metacommunity
        unique = 0

        ## If you set your carrying capacity too high relative to the size of your
        ## metacommunity then you'll get stuck drawing duplicates over and over
        idiot_count = 0
        while not unique:
            ## Sample from the metacommunity
            migrant_draw = np.random.multinomial(1, self.immigration_probabilities, size=1)
            #print("Immigration event - {}".format(np.where(migrant_draw == 1)))
            #print("Immigrant - {}".format(self.species[np.where(migrant_draw == 1)[1][0]]))
            new_species = self.species[np.where(migrant_draw == 1)[1][0]]
            ##TODO: Should set a flag to guard whether or not to allow multiple colonizations
            if new_species[0] not in [x[0] for x in self.local_community]:
                ## Got a species not in the local community
                unique = 1
            else:
                #print("multiple colonization forbidden: sp id {}".format(new_species[0]))
                idiot_count +=1
            if idiot_count > MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY:
               msg = """\nMetacommunity is exhausted w/ respect to local
               community. Either expand the size of the metacommunity,
               decrease the carrying capacity, or switch on multiple
               migration (unimplemented)."""
               sys.exit(msg)

        return new_species

    def migrate_step(self):
        """ Allow multiple colonizations. In this case we return the sampled species
        as well as a bool reporting whether or not this is the first colonization of
        this species into the local community so that the coltime can be recorded.
        multiple colonizations of a species do not update coltime, but we record them
        for migration rate calculation."""

        init_col = True
        migrant_draw = np.random.multinomial(1, self.immigration_probabilities, size=1)
        new_species = self.species[np.where(migrant_draw == 1)[1][0]]
        if new_species[0] in [x[0] for x in self.local_community]:
            #print("Multiple colonization: sp id {}".format(new_species[0]))
            ## This is a post-colonization migrant so record the event and tell downstream
            ## not to update the colonization time.
            self.post_colonization_migrants[new_species[0]] += 1
            init_col = False
        else:
            ## This is a new migrant so init the post-colonization count
            self.post_colonization_migrants[new_species[0]] = 0

        return new_species, init_col


    def step(self, nsteps=1, time=0, invasion_time=100000, invasiveness=0.1):
        for step in range(nsteps):
            ## If there are any members of the local community
            if self.local_community:
                ## Do the magic to remove one individual from the local community
                ## After this function returns K = K - 1
                self.death_step(invasion_time, invasiveness)

            ## Check probability of an immigration event
            if np.random.random_sample() < self.colonization_rate:

                ## Grab the new colonizing species
                init_colonization = True
                if self.allow_multiple_colonizations:
                    new_species, init_colonization = self.migrate_step()
                else:
                    new_species = self.migrate_no_dupes_step()

                ## Only record coltime if this is the first time this species enters the local community 
                if init_colonization:
                    self.divergence_times[(new_species[0], False)] = self.current_time

                ## Only set the invasive species once at the time of next migration post invasion time
                ## If invasion time is < 0 this means "Don't do invasive"
                if not invasion_time < 0:
                    if self.invasive == -1 and time >= invasion_time:
                        self.invasive = (new_species[0], False)
                        print("setting invasive species {} at time {}".format(self.invasive, self.current_time))
                        self.invasion_time = self.current_time

                ## Add the colonizer to the local community, record the colonization time
                self.local_community.append((new_species[0], False))
                self.colonizations += 1
            else:
                ## Sample from the local community, including empty demes
                ## Sample all available from local community (community grows slow in volcanic model)
                ## Also, lots of early turnover
                ## This all was only true before i implemented the full rosindell/harmon model,
                ## There are no more empty demes in the current config
                self.local_community.append(random.choice(self.local_community))
    
                ## Sample only from available extant species (early pops grow quickly in the volcanic model)
                ## If you do this, the original colonizer just overwhelms everything else
                ## This is more similar to the Rosindell and Harmon model, in which they simply
                ## prepopulate the island entirely with one species. This is effectively the same
                #self.local_community.append(random.choice([x for x in self.local_community if not x[0] == 0]))
    
            ## update current time
            self.current_time += 1


    def get_abundances(self, octaves=False):
        ## Make a counter for the local_community, counts the number of
        ## individuals w/in each species
        abundances = collections.Counter([x[0] for x in self.local_community])

        ## If we were doing mode=volcanic then there may be some remaining
        ## space in our carrying capacity that is unoccupied (indicated by
        ## None in the abundances.keys()
        try:
            abundances.pop(None)
        except KeyError:
            pass

        ## Make a set of abundances to get all the unique values
        abundance_classes = set(abundances.values())

        ## Now for each abundance class you have to go through and
        ## count the number of species at that abundance.
        ## This is currently stupid because in python there's no
        ## straightforward way to get keys from values in a dict.
        abundance_distribution = collections.OrderedDict()
        for i in abundance_classes:
            count = 0
            for _, v in abundances.items():
                if v == i:
                    count += 1
            abundance_distribution[i] = count
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

    ## How strong is the bottleneck? Strength should be interpreted as percent of local
    ## community to retain
    def bottleneck(self, strength=1):
        reduction = int(round(self.local_inds * strength))
        self.local_community = self.local_community[:reduction]

        ## First remove the extinct species from the species list
        pre = len(self.species_objects)
        self.species_objects = [s for s in self.species_objects if s.uuid in self.local_community]
        ## Update the extinction counter
        self.extinctions += (pre - len(self.species_objects))

        sp = self.species_objects
        ## Update abundances per species that survived the bottleneck
        for i, s in enumerate(sp):
            if s.uuid in self.local_community:
                abund = self.local_community.count(s.uuid)
                s.update_abundance(abund)
                self.species_objects[i] = s

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
                                        exponential=self.exponential, abundance=abundance,\
                                        meta_abundance=meta_abundance))

        for s in self.species_objects:
            s.simulate_seqs()
            s.get_sumstats()
            ## For debugging invasives
            #if s.abundance > 1000:
            #    print("\n{}".format(s))

    def set_species(self, species_objects):
        self.species_objects = species_objects

    def get_species(self):
        return(self.species_objects)


if __name__ == "__main__":
    data = implicit_BI()
    #data.set_metacommunity("uniform")
    data.set_metacommunity("metacommunity_LS4.txt")
    #data.prepopulate(mode="landbridge")
    data.prepopulate(mode="volcanic")
    for i in range(100000):
        if not i % 10000:
            print("Done {}".format(i))
            #print(i, len(data.local_community), len(set(data.local_community)))
            #print(data.local_community)
        data.step()
    abundance_distribution = data.get_abundances(octaves=False)
    print("Species abundance distribution:\n{}".format(abundance_distribution))
    #print("Colonization times per species:\n{}".format(data.divergence_times))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()
    print("Species:\n{}".format(data.get_species()))
    print("Extinction rate - {}".format(data.extinctions/float(data.current_time)))
    print("Colonization rate - {}".format(data.colonizations/float(data.current_time)))
