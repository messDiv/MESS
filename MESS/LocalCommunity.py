#!/usr/bin/env python2.7

try:
    import matplotlib.pyplot as plt
    from ascii_graph import Pyasciigraph
except:
    print("matplotlib and/or ascii_graph not installed, so plotting is disabled.")
from scipy.stats import logser
from collections import OrderedDict
import collections
import numpy as np
import itertools
import random
import sys
import os
import MESS

from MESS.stats import shannon
try:
    from species import species
except:
    print("Species module failed to load, things probably won't work right.")

import logging
LOGGER = logging.getLogger(__name__)

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1500

class LocalCommunity(object):

    def __init__(self, name=None, K=5000, colrate=0.01, allow_multiple_colonizations=False, \
                mig_clust_size=1, exponential=False, quiet=False):
        self.quiet = quiet

        if name is None:
            raise MESSError("LocalCommunity must be named")
        else:
            self.name = name

        ## If you add a parameter to this dictionary you need
        ## to also add a short description to the LOCAL_PARAMS dict
        ## at the end of this file
        ##
        ## Also be sure to add it to _paramschecker so the type gets set correctly
        self.paramsdict = OrderedDict([
                        ("name", self.name),
                        ("mode", "volcanic"),
                        ("K", K),
                        ("colrate", colrate),
                        ("age", 1e6),
                        ("mig_clust_size", mig_clust_size),
                        ("allow_multiple_colonizations", allow_multiple_colonizations),
        ])

        ## List for storing species objects that have had sequence
        ## simulated and sumstats calculated
        self.species_objects = []
        self.exponential = exponential

        ## Settings specific to the uniform metacommunity
        ## This is individuals per species
        self.uniform_inds = 1000000
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
        self.allow_multiple_colonizations = allow_multiple_colonizations
        self.mig_clust_size = mig_clust_size

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
        self.invasiveness = 0
        ## Track how many invasives differentially survived
        self.survived_invasives = 0
        self.invasion_time = -1

        ###########################################################
        ## Variables for non-neutral assembly processes
        ## Note all non-neutral assembly value defaults result
        ## in neutral processes.
        ###########################################################

        ## Fundamental parameters of the metacommunity tree
        self.metcommunity_tree_height = 0
        self.trait_evolution_rate_parameter = 0

        ## An ID by trait value dictionary
        self.species_trait_values = {}

        ## Toggle competitive exclusion
        self.competitive_exclusion = False
        self.competitive_strength = 0

        ## Toggle environmental filtering
        self.environmental_filtering = False
        self.environmental_optimum = 0
        self.environmental_strength = 0

        ## Toggle params for holeing species death probabilities
        self.species_death_probability = {}
        #species death probability will be an unchanging dictionary for enviormental environmental_filtering
        #but will have to be updated for competitive exlcusion models
        self.individual_death_probabilites = []

        ## Set the default metacommunity and prepopulate the island
        self.set_metacommunity("logser")
        self.prepopulate(mode=self.paramsdict["mode"], quiet=quiet)


    ## Return fraction of equilibrium obtained by the local community
    def _lambda(self):
        founder_flags = [x[1] for x in self.local_community]
        percent_equil = float(founder_flags.count(False))/len(founder_flags)
        return percent_equil

    def _paramschecker(self, param, newvalue, quiet=False):
        """ Raises exceptions when params are set to values they should not be"""
        ## TODO: This should actually check the values and make sure they make sense
        ## TODO: Also check here if you're setting the mode parameter you have to rerun prepopulate
        try:
            LOGGER.debug("set param {} - {} = {}".format(self, param, newvalue))

            ## Cast params to correct types
            if param in ["K", "mig_clust_size", "age"]:
                self.paramsdict[param] = int(float(newvalue))
            elif param in ["colrate"]:
                self.paramsdict[param] = float(newvalue)
            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            ## Do something intelligent here?
            raise

    def write_params(self, outfile=None, append=True):
        """
        Write out the parameters for this island to a file.
        Normally this isn't called directly, but by the main
        simulation engine.

        append 
        """
        if outfile is None:
            raise MESSError("LocalCommunity.write_params outfile must be specified.")

        ## If not appending then we are overwriting
        if append:
            filemode = 'a'
        else:
            filemode = 'w'

        with open(outfile, filemode) as paramsfile:
            ## Only write the full header if not appending
            if not append:
                header = "------- MESS params file (v.{})".format(MESS.__version__)
                header += ("-"*(80-len(header)))
                paramsfile.write(header)
            
            header = "------- LocalCommunity params: {}".format(self.name)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            for key, val in self.paramsdict.iteritems():
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)

                padding = (" "*(20-len(paramvalue)))
                paramkey = self.paramsdict.keys().index(key)
                paramindex = " ## [{}] ".format(paramkey)
                LOGGER.debug("{} {} {}".format(key, val, paramindex))
                #name = "[{}]: ".format(paramname(paramkey))
                name = "[{}]: ".format(key)
                #description = paraminfo(paramkey, short=True)
                description = LOCAL_PARAMS[key]
                paramsfile.write("\n" + paramvalue + padding + \
                                        paramindex + name + description)

            paramsfile.write("\n")        


    def set_metacommunity(self, infile, random=False):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there

        random=True will set random trait values in the range [0-1]
        """
        if infile in ["logser", "uniform"]:
            if infile == "logser":
                ## Parameter of the logseries distribution
                p = .98
                self.abundances = logser.rvs(p, size=self.local_inds)
            else:
                self.abundances = [self.uniform_inds] * self.uniform_species


            self.species = ["t"+str(x) for x in range(0, len(self.abundances))]

            if random:
                self.species_trait_values = {x:y for x,y in zip(self.species, np.random.rand(len(self.species)))}
            else:
                self.species_trait_values = {x:y for x,y in zip(self.species, np.random.rand(len(self.species)))}
        else:
            #infile="SpInfo.txt"
            if os.path.isfile(infile):
                with open(infile, 'r') as inf:
                    lines = inf.readlines()
                    self.metcommunity_tree_height = float(lines[0].split()[0])
                    self.trait_evolution_rate_parameter = float(lines[1].split()[0])

                    for i in range(2,len(lines)):
                        info = lines[i].split()
                        self.species.append(info[0])
                        self.species_trait_values[info[0]] = (float(info[1]))
                        self.abundances.append(int(info[2]))

                    ### Try fetching traits as well. If there are no traits then assume we're not doing trait based
                    #try:
                    #    inf.seek()
                    #    self.species_trait_values = {x:y for x,y in enumerate([float(line.split()[1]) for line in inf])}
                    #except:
                    #    self.species_trait_values = {x:y for x,y in enumerate(np.zeros(len(self.abundances)))}
            else:
                raise Exception("Bad metacommunity input - ".format(infile))

        self.total_inds = sum(self.abundances)
        self.immigration_probabilities = [float(self.abundances[i])/self.total_inds for i in range(len(self.abundances))]

        self.maxabundance = np.amax(self.immigration_probabilities)

        ## Init post colonization migrants counters
        self.post_colonization_migrants = {x:0 for x in self.species}

        if self.environmental_filtering | self.competitive_exclusion:
            self.calculate_death_probabilities()


    def __str__(self):
        return "<LocalCommunity {}: Shannon's Entropy {}>".format(self.name, shannon(self.get_abundances(octaves=False)))


    def calculate_death_probabilities(self):
        if self.environmental_filtering:
            #environmental optimum is a random value between (-rate of BM * tree depth, rate of BM * tree depth)
            scale = float(4 * self.trait_evolution_rate_parameter * self.metcommunity_tree_height)
            self.environmental_optimum = np.random.uniform(-scale, scale)

            #determine strength of environment
            self.environmental_strength = float(np.random.uniform(0.1,100))
            self.weight = self.trait_evolution_rate_parameter * self.metcommunity_tree_height * (self.environmental_strength ** 2)

            for i in range(len(self.species)):
                deathProb = 1 - (np.exp(-((self.species_trait_values[self.species[i]] - self.environmental_optimum) ** 2)/self.weight))
                self.species_death_probability[self.species[i]] = deathProb

        if self.competitive_exclusion:
            #determine strength of competition
            self.competitive_strength = float(np.random.uniform(0.1,100))
            self.weight = self.trait_evolution_rate_parameter * self.metcommunity_tree_height * (self.competitive_strength ** 2)


    def prepopulate(self, mode="landbridge", quiet=False):
        ## Clean up local_community if it already exists
        self.local_community = []
        if mode == "landbridge":
            ## prepopulate the island w/ total_inds individuals sampled from the metacommunity
            init_community = np.random.multinomial(self.local_inds, self.immigration_probabilities, size=1)
            for i, x in enumerate(init_community[0]):
                if x:
                    self.local_community.extend([self.species[i]] * x)
                    ## Set founder flag
            self.local_community = [tuple([x, True]) for x in self.local_community]
            #self.local_community = list(itertools.chain.from_iterable(self.local_community))

            ## All species diverge simultaneously upon creation of the island.
            for taxon in self.local_community:
                self.divergence_times[taxon] = 1
        elif mode == "volcanic":
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            new_species = (self.species[self.immigration_probabilities.index(self.maxabundance)], True)
            self.local_community.append(new_species)
            ## prepopulate volcanic either with all the most abundant species in the metacommunity
            ## or with one sample of this species and a bunch of "emtpy deme space". The empty
            ## demes screw up competition/environmental filtering models
            if True:
               self.local_community.extend([new_species] * (self.local_inds - 1))
            else:
                for i in range(1,self.local_inds):
                    self.local_community.append((None,True))
            self.divergence_times[new_species] = 1
        else:
            raise MESSError(BAD_MODE_PARAMETER.format(mode))

        if not quiet:
            print("  Initializing local community:")
            print("    N species = {}".format(len(set([x[0] for x in self.local_community]))))
            print("    N individuals = {}".format(len(self.local_community)))
                            

    def death_step(self):
        ## Select the individual to die

        ##currently this will fail under volcanic model because the entire local community will go extinct
        if self.environmental_filtering:
            species_inLocal = [x[0] for x in self.local_community if x[0] != None]
            #wont need if statement after volcanic model changed
            #print(species_inLocal)
            death_probabilites = []
            for i in range(len(species_inLocal)):
                death_probabilites.append(self.species_death_probability[species_inLocal[i]])
            #print(death_probabilites)

            self.individual_death_probabilites = death_probabilites / sum(death_probabilites)
            #print(self.individual_death_probabilites)
            #print(sum(self.individual_death_probabilites))
            #index = np.arange(0,len(self.individual_death_probabilites),1)
            sample = np.random.multinomial(1, self.individual_death_probabilites)
            #print(sample)
            #self.local_community[victim_index==1]
            victim_index = int(np.arange(0,len(self.individual_death_probabilites),1)[sample==1])
            #print(victim_index)
            #print(self.local_community[victim_index])
            victim = self.local_community[victim_index]

        elif self.competitive_exclusion:
            species_inLocal = [x[0] for x in self.local_community if x[0] != None]

            local_traits = []
            for i in range(len(species_inLocal)):
                local_traits.append(self.species_trait_values[species_inLocal[i]])
            mean_local_trait = np.mean(local_traits)

            for i in range(len(self.species)):
                deathProb = (np.exp(-((self.species_trait_values[self.species[i]] - mean_local_trait) ** 2)/self.weight))
                self.species_death_probability[self.species[i]] = deathProb
            #print(species_inLocal)
            #print(self.species_death_probability)

            death_probabilites = []
            for i in range(len(species_inLocal)):
                death_probabilites.append(self.species_death_probability[species_inLocal[i]])

            self.individual_death_probabilites = death_probabilites / sum(death_probabilites)
            #print(np.sum(self.individual_death_probabilites))
            sample = np.random.multinomial(1, self.individual_death_probabilites)
            #print(sample)
            victim_index = int(np.arange(0,len(self.individual_death_probabilites),1)[sample==1])
            #print(victim_index)
            #print(self.local_community[victim_index])
            victim = self.local_community[victim_index]

        else:
            ## If not trait based just select one individual randomly (neutral0
            victim = random.choice(self.local_community)
        ## If no invasive has invaded then just do the normal sampling
        if self.invasive == -1:
            ## If no invasive species yet just go on
            pass
        else:
            ## If invasiveness is less than the random value remove the invasive individual
            ## else choose a new individual
            if victim == self.invasive and np.random.rand() < self.invasiveness:
                self.survived_invasives += 1
                victim = random.choice(self.local_community)

        try:
            #print(victim)
            self.local_community.remove(victim)
        except Exception as inst:
            print(victim)
            print(self.local_community)
            raise

        ## Record local extinction events
        if not victim in self.local_community:
            ## This was supposed to not record "extinctions" of empty deme space
            ## but it fucks up the calculation of extinction rate
            ##if not victim[0] == None:
            if True:
                self.extinctions += 1
                try:
                    ## reset post colonization migrant count for this species
                    ## Record the lifetime of this species and remove their record from divergence_times
                    self.post_colonization_migrants[victim[0]] = 0
                    self.extinction_times.append(self.current_time - self.divergence_times[victim])
                    del self.divergence_times[victim]
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
            #print(new_species)
            if new_species not in [x[0] for x in self.local_community]:
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
        migrant_draw = np.random.multinomial(1, self.immigration_probabilities, size=1).argmax()
        new_species = self.species[migrant_draw]
        if new_species in [x[0] for x in self.local_community]:
            #print("Multiple colonization: sp id {}".format(new_species[0]))
            ## This is a post-colonization migrant so record the event and tell downstream
            ## not to update the colonization time.
            self.post_colonization_migrants[new_species] += 1
            init_col = False
        else:
            ## This is a new migrant so init the post-colonization count
            self.post_colonization_migrants[new_species] = 0
            #print("New immigrant {}\t{}".format(new_species, self.post_colonization_migrants))

        return new_species, init_col


    def step(self, nsteps=1):
        for step in range(nsteps):
            ## Check probability of an immigration event
            if np.random.random_sample() < self.colonization_rate:
                ## If clustered migration remove the necessary number of additional individuals
                for _ in range(self.paramsdict["mig_clust_size"]):
                    self.death_step()

                ## Grab the new colonizing species
                ## the init_colonization flag is used to test whether to update the divergence time
                init_colonization = True
                if self.paramsdict["allow_multiple_colonizations"]:
                    new_species, init_colonization = self.migrate_step()
                else:
                    new_species = self.migrate_no_dupes_step()

                ## Only record coltime if this is the first time this species enters the local community
                if init_colonization:
                    self.divergence_times[(new_species, False)] = self.current_time

                ## Only set the invasive species once at the time of next migration post invasion time
                ## If invasion time is < 0 this means "Don't do invasive"
                if not self.invasion_time < 0:
                    if self.invasive == -1 and self.current_time >= self.invasion_time:
                        self.invasive = (new_species, False)
                        LOGGER.info("setting invasive species {} at time {}".format(self.invasive, self.current_time))
                        self.invasion_time = self.current_time

                ## Add the colonizer to the local community, record the colonization time
                self.local_community.extend([(new_species, False)] * self.mig_clust_size)
                self.colonizations += 1
            else:
                self.death_step()
                ## Sample from the local community, including empty demes
                ## Sample all available from local community (community grows slow in volcanic model)
                ## Also, lots of early turnover
                ## This all was only true before i implemented the full rosindell/harmon model,
                ## There are no more empty demes in the current config
                chx = np.random.randint(0, len(self.local_community))
                self.local_community.append(self.local_community[chx])
                ## Much faster than random.choice
                #self.local_community.append(random.choice(self.local_community))

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
        for UUID, tcol in self.divergence_times.items():
            if UUID in self.local_community:
                meta_abundance = -1
                for x in self.species:
                    if UUID[0] == x:
                        meta_abundance = self.abundances[int(x[1:])]
                #meta_abundance = [x[1] for x in self.abundances if x[0] == UUID[0]]
                #meta_abundance = self.abundances[self.species.index(UUID[0])]
                abundance = self.local_community.count(UUID)
                #print(self.local_community)
                try:
                    tdiv = self.current_time - tcol
                    self.species_objects.append(species(UUID=UUID, colonization_time=tdiv,\
                                        exponential=self.exponential, abundance=abundance,\
                                        meta_abundance=meta_abundance,
                                        migration_rate=self.post_colonization_migrants[UUID[0]]/float(tdiv)))
                except:
                    print(UUID)
                    print(self.post_colonization_migrants)
                    raise
        for s in self.species_objects:
            s.simulate_seqs()
            s.get_sumstats()
            ## For debugging invasives
            #if s.abundance > 1000:
            #    print("\n{}".format(s))

    def set_species(self, species_objects):
        self.species_objects = species_objects
        print(self.species_objects)

    def get_species(self):
        return(self.species_objects)


#############################
## Model Parameter Info Dicts
#############################
LOCAL_PARAMS = {
    "name" : "Local community name",\
    "mode" : "Local community formation mode (volcanic/landbridge)",\
    "K" : "Local carrying capacity",\
    "colrate" : "Colonization rate into local community",\
    "mig_clust_size" : "# of individuals per colonization event",\
    "age" : "Local community age",\
    "allow_multiple_colonizations" : "Toggle allowing post colonization migration",
}


#############################
## Error messages
#############################
BAD_MODE_PARAMETER = """
    Unrecognized local community mode. Options are 'landbridge' or'volcanic'.
    You put {}.
"""


if __name__ == "__main__":
    data = LocalCommunity("tmp", K=500, allow_multiple_colonizations=True, quiet=True)
    #data.set_metacommunity("uniform")
    #data.environmental_filtering = True
    #data.competitive_exclusion = True
    data.set_metacommunity("uniform")
    data.prepopulate(mode="landbridge", quiet=True)

    for i in range(10000):
        if not i % 100:
            print("{} ".format(i)),
            #print(i, len(data.local_community), len(set(data.local_community)))
            #print(data.local_community)
        data.step()
    print(data.local_community)
    abundance_distribution = data.get_abundances(octaves=False)
    print("\n")
    print(data)

    data.simulate_seqs()
    print("Species abundance distribution:\n{}".format(abundance_distribution))
    #print("Colonization times per species:\n{}".format(data.divergence_times))
    #plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    #plt.show()
    print("Species:\n{}".format(data.get_species()))
    print("Extinction rate - {}".format(data.extinctions/float(data.current_time)))
    print("Colonization rate - {}".format(data.colonizations/float(data.current_time)))
