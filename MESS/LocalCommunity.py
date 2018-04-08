#!/usr/bin/env python2.7

try:
    import matplotlib.pyplot as plt
    from ascii_graph import Pyasciigraph
except:
    print("matplotlib and/or ascii_graph not installed, so plotting is disabled.")
from scipy.stats import logser
from collections import OrderedDict
import collections
import pandas as pd
import numpy as np
import itertools
import random
import sys
import os
import MESS

from util import MESSError
from stats import shannon
from species import species

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

        ## Dictionary of 'secret' parameters that most people won't want to mess with
        self._hackersonly = dict([
                        ("allow_empty", False),
        ])

        ## list for storing the state of our local community. The list is much faster
        ## than a numpy array, plus lists are mutable, which is convenient.
        ## I have tried making the local community a numpy array twice, and both times
        ## it turned out to suck both times.
        self.local_community = []
        self.founder_flags = []

        ## np array for storing info about each species in the local community. This is for
        ## info that would be annoying or impossible to keep track of "per individual". 
        ## The array is of length == len(Counter(local_community))
        ##  * "ids" - species identifiers present in local community
        ##  * "colonization_times" - Colonization times per species
        ##  * "post_colonization_migrants" - Count of post colonization migrants per species
        self.local_info = []

        ## The regional pool that this local community belongs to
        ## this is updated by Region._link_local(), so don't set it by hand
        self.region = ""

        self.files = dict({
                "full_output": [],
                })

        ## summary stats dict
        self.stats = pd.Series(
            index=["_lambda",
                   "generation",
                   "K",
                   "colrate",
                   "colrate_calculated",
                   "extrate_calculated",
                   "R",
                   "shannon",
                   "mean_pi",
                   "stdv_pi",
                   "median_pi",
                   "iqr_pi",
                   "mean_dxy",
                   "stdv_dxy",
                   "median_dxy",
                   "iqr_dxy",
                   "sgd",
                   "trees"]).astype(np.object)

        ## List for storing species objects that have had sequence
        ## simulated and sumstats calculated
        self.species_objects = []
        self.exponential = exponential

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


    ## Return fraction of equilibrium obtained by the local community
    def _lambda(self):
        percent_equil = 0
        try:
            percent_equil = float(self.founder_flags.count(False))/len(self.founder_flags)
        except ZeroDivisionError:
            pass
        return percent_equil


    def _paramschecker(self, param, newvalue, quiet=False):
        """ Raises exceptions when params are set to values they should not be"""
        ## TODO: This should actually check the values and make sure they make sense
        ## TODO: Also check here if you're setting the mode parameter you have to rerun prepopulate
        try:
            LOGGER.debug("set param {} - {} = {}".format(self, param, newvalue))

            ## Cast params to correct types
            if param == "K":
                self.paramsdict[param] = int(float(newvalue))
                self.local_community = []
                self.prepopulate(mode=self.paramsdict["mode"], quiet=True)

            if param in ["mig_clust_size", "age"]:
                self.paramsdict[param] = int(float(newvalue))

            elif param in ["colrate"]:
                self.paramsdict[param] = float(newvalue)

            elif param == "mode":
                ## Must reup the local community if you change the mode
                self.paramsdict[param] = newvalue
                self.local_community = []
                self.prepopulate(mode=self.paramsdict["mode"], quiet=True)

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


    def __str__(self):
        return "<LocalCommunity {}: Shannon's Entropy {}>".format(self.name, shannon(self.get_abundances(octaves=False)))


    ## Not done
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


    def prepopulate(self, quiet=False):
        LOGGER.debug("prepopulating local_community")

        if not self.region:
            msg = "Skip populating the local community as it is unlinked to a region."
            LOGGER.error(msg)
            if MESS.__interactive__: print("    {}".format(msg))
            return

        ## Clean up local_community if it already exists
        self.local_community = []

        if self.paramsdict["mode"] == "landbridge":
            ## prepopulate the island w/ a random sample from the metacommunity
            ## TODO: The underscore here is ignoring trait values
            self.local_community, _ = self.region.get_nmigrants(self.paramsdict["K"])

        elif self.paramsdict["mode"]  == "volcanic":
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            ## TODO: The _ is a standin for trait values, have to do something with them
            new_species, _ = self.region.get_most_abundant()

            ## prepopulate volcanic either with all the most abundant species in the metacommunity
            ## or with one sample of this species and a bunch of "emtpy deme space". The empty
            ## demes screw up competition/environmental filtering models
            if self._hackersonly["allow_empty"]:
                self.local_community = [None] * self.paramsdict["K"]
                self.local_community[0] = new_species
            else:
                self.local_community = [new_species] * self.paramsdict["K"]
        else:
            raise MESSError(BAD_MODE_PARAMETER.format(mode))

        ## At time 0 all individuals are founders and all colonization times
        ## and post colonization migrant counts are zero
        ids = list(set(self.local_community))
        self.local_info = pd.DataFrame([], columns = ids,\
                                            index = ["colonization_times",\
                                                     "post_colonization_migrants"])
        self.local_info = self.local_info.fillna(0)
        self.founder_flags = [True] * len(self.local_community)
        print(self.local_info)
        
        if self.environmental_filtering | self.competitive_exclusion:
            self.calculate_death_probabilities()

        if not quiet:
            print("  Initializing local community:")
            print("    N species = {}".format(len(set(self.local_community))))
            print("    N individuals = {}".format(len(self.local_community)))
                            

    ## Not updated
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

        ## Clean up local community list and founder flag list
        idx = self.local_community.index(victim)
        self.local_community.pop(idx)
        self.founder_flags.pop(idx)

        ## Record local extinction events
        if not victim in self.local_community:
            self.extinctions += 1
            try:
                ## Record the lifetime of this species and remove their record from divergence_times
                ## Remove the species from the local_info array
                self.extinction_times.append(self.current_time - self.local_info[victim]["divergence_times"])
                self.local_info.drop(columns=[victim])
            except:
                ## The empty deme will make this freak
                pass
            ## If the invasive prematurely goes extinct just pick a new one
            if victim[0] == self.invasive:
                LOGGER.info("invasive went extinct")
                self.invasive = -1


    ## Not updated
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
            ## TODO: _ here is the ignored trait value, for now
            new_species, _ = self.region.get_migrant()
            if new_species not in self.local_community:
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

        ## This is a new migrant so init local_info for it
        self.local_info[new_species] = 0
        self.local_info[new_species]["post_colonization_migrants"] = 0
        self.local_info[new_species]["colonization_times"] = self.current_time
        return new_species


    ## Not updated
    def migrate_step(self):
        """ Allow multiple colonizations. In this case we return the sampled species
        as well as a bool reporting whether or not this is the first colonization of
        this species into the local community so that the coltime can be recorded.
        multiple colonizations of a species do not update coltime, but we record them
        for migration rate calculation."""

        new_species = self.region.get_migrant()
        if new_species in self.local_community:
            ## This is a post-colonization migrant so record the event and tell downstream
            ## not to update the colonization time.
            self.local_info[new_species]["post_colonization_migrants"] += 1
        else:
            ## This is a new migrant so init local_info for it
            self.local_info[new_species] = 0
            self.local_info[new_species]["colonization_times"] = self.current_time

        return new_species


    ## Not updated
    def step(self, nsteps=1):
        for step in range(nsteps):
            ## Check probability of an immigration event
            if np.random.random_sample() < self.paramsdict["colrate"]:
                ## If clustered migration remove the necessary number of additional individuals
                for _ in range(self.paramsdict["mig_clust_size"]):
                    self.death_step()

                ## Grab the new colonizing species
                ## the init_colonization flag is used to test whether to update the divergence time
                if self.paramsdict["allow_multiple_colonizations"]:
                    new_species = self.migrate_step()
                else:
                    new_species = self.migrate_no_dupes_step()

                ## Only set the invasive species once at the time of next migration post invasion time
                ## If invasion time is < 0 this means "Don't do invasive"
                if not self.invasion_time < 0:
                    if self.invasive == -1 and self.current_time >= self.invasion_time:
                        self.invasive = new_species
                        LOGGER.info("setting invasive species {} at time {}".format(self.invasive, self.current_time))
                        self.invasion_time = self.current_time

                ## Add the colonizer to the local community, record the colonization time
                self.local_community.append(new_species)
                self.founder_flags.append(False)
                self.colonizations += 1
            else:
                self.death_step()
                ## Sample all available from local community (community grows slow in volcanic model)
                ## This is the fastest way to sample from a list. >4x faster than np.random.choice
                chx = random.choice(self.local_community)
                self.local_community.append(chx)
                idx = self.local_community.index(chx)
                self.founder_flags.append(self.founder_flags[idx])

                ## Sample only from available extant species (early pops grow quickly in the volcanic model)
                ## If you do this, the original colonizer just overwhelms everything else
                ## This is more similar to the Rosindell and Harmon model, in which they simply
                ## prepopulate the island entirely with one species. This is effectively the same
                #self.local_community.append(random.choice([x for x in self.local_community if not x == None]))

            ## update current time
            self.current_time += 1


    def get_abundances(self, octaves=False):
        ## Make a counter for the local_community, counts the number of
        ## individuals w/in each species
        abundances = collections.Counter(self.local_community)

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
    ## Not done
    def bottleneck(self, strength=1):
        reduction = int(round(self.paramsdict["K"] * strength))
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
                    print(self.local_info["post_colonization_migrants"])
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


    def get_stats(self):
        LOGGER.debug("Entering get_stats()")
        self.simulate_seqs()
        self.stats._lambda = self._lambda()
        self.stats.colrate_calculated = self.colonizations/float(self.current_time)
        self.stats.extrate_calculated = self.extinctions/float(self.current_time)
        self.stats.shannon = shannon(self.get_abundances(octaves=False))

        sp = self.get_species()
        if sp:
            pis = np.array([x.pi_local for x in sp])
            dxys = np.array([x.dxy for x in sp])
            self.stats.mean_pi = np.mean(pis)
            self.stats.mean_dxy = np.mean(dxys)

        return self.stats


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
    data = MESS.Region("tmp")
    loc = LocalCommunity("wat", K=100)
    data._link_local(loc)
    print(loc)
    print(loc.local_info)
    ## Allow for either having or not having empty demes
    assert(len(collections.Counter(loc.local_community)) <= 2)
    loc.paramsdict["mode"] = "landbridge"
    loc.prepopulate()
    assert(len(collections.Counter(loc.local_community)) > 3)

    #print(loc.local_info[:10])
    print(loc.get_abundances())
    loc.step(10000)
    print(loc.local_info)
    print(loc)
