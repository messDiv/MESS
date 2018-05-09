#!/usr/bin/env python2.7

from scipy.stats import logser
from collections import OrderedDict
from scipy.stats import iqr
import collections
import pandas as pd
import numpy as np
import itertools
import random
import sys
import os
import MESS

from util import MESSError, _tuplecheck, sample_param_range
from stats import shannon, SAD, SGD
from species import species

import logging
LOGGER = logging.getLogger(__name__)

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1500


class LocalCommunity(object):

    def __init__(self, name=None, K=5000, colrate=0.01, \
                mig_clust_size=1, quiet=False):
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
                        ("age", 100000),
                        ("mig_clust_size", mig_clust_size),
                        ("filtering_optimum", 100)
        ])

        ## A dictionary for holding prior ranges for values we're interested in
        self._priors = dict([
                        ("K", []),
                        ("colrate", []),
                        ("mig_clust_size", [])
        ])

        ## Dictionary of 'secret' parameters that most people won't want to mess with
        ##  * allow_empty is a switch on whether to fully populate local community
        ##      with one species for 'volcanic' mode or to introduce just one
        ##      individual and populate the rest of K with 'empty' demes.
        ##  * outdir is inherited from the Region.simulate() command, so users
        ##      should normally not mess with this.
        self._hackersonly = dict([
                        ("allow_empty", False),
                        ("outdir", []),
        ])

        ## list for storing the state of our local community. The list is much faster
        ## than a numpy array, plus lists are mutable, which is convenient.
        ## I have tried making the local community a numpy array twice, and both times
        ## it turned out to suck both times.
        self.local_community = []
        self.founder_flags = []

        ## pandas Data Frame for storing info about each species in the local community. This is for
        ## info that would be annoying or impossible to keep track of "per individual".
        ## The column headers are species ids, and the dataframe is always of length ==
        ## len(Counter(local_community)), because we prune extinct.
        ## Fields:
        ##  * "colonization_times" - Colonization times per species
        ##  * "post_colonization_migrants" - Count of post colonization migrants per species
        ##  * "abundance_through_time" - A record of abundance through time for this species
        self.local_info = pd.DataFrame([])

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
                   "trees",
                   "mean_ltr",
                   "var_ltr",
                   "mean_rtr",
                   "var_rtr",
                   "mean_dif",
                   "var_dif",
                   "kurtosis",
                   "skewness"]).astype(np.object)

        ## List for storing species that have had sequence
        ## simulated and sumstats calculated
        self.species = []

        ## Various counters for counting stuff that happens
        self.extinctions = 0
        self.colonizations = 0
        self.current_time = 0

        ## Vector for tracking lifetimes of excinct species
        ## We can plot this and average it to figure out how long
        ## species hang around in the local community
        self.extinction_times = []

        ## Dicts for tracking esoteric shit if _log() is called
        self.lambda_through_time = OrderedDict({})
        self.species_through_time = OrderedDict({})

        ## Track number of rejections per death step
        self.rejections = []

        ## Invasiveness is unimplemented and may not make it into MESS proper.
        ## The invasive species identity
        self.invasive = -1
        self.invasiveness = 0
        ## Track how many invasives differentially survived
        self.survived_invasives = 0
        self.invasion_time = -1


    def _copy(self):
        """ Create a new clean copy of this LocalCommunity."""
        LOGGER.debug("Copying LocalCommunity - {}".format(self.name))
        new = LocalCommunity(self.name)

        new.region = self.region
        new._priors = self._priors

        new.paramsdict = self.paramsdict
        ## Get sample from prior range on params that may have priors
        for param in ["K", "colrate", "mig_clust_size"]:
            ## if _priors is empty then this param is fixed
            if np.any(self._priors[param]):
                self.paramsdict[param] = sample_param_range(new._priors[param])[0]

        new._hackersonly = self._hackersonly
        LOGGER.debug("Copied LocalCommunity - {}".format(self.paramsdict))
        return new


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
            if param in ["K", "mig_clust_size", "age"]:
                tup = _tuplecheck(newvalue, dtype=int)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup)[0]
                else:
                    self.paramsdict[param] = tup

            elif param in ["colrate"]:
                tup = _tuplecheck(newvalue, dtype=float)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup)[0]
                else:
                    self.paramsdict[param] = tup

            elif param == "mode":
                ## Must reup the local community if you change the mode
                self.paramsdict[param] = newvalue

            elif param == "filtering_optimum":
                self.paramsdict[param] = float(newvalue)

            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            LOGGER.debug("Bad parameter - {} {}".format(param, newvalue))
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


    def _log(self, full=False):
        """ A function for occasionally logging a ton of information through time.
        Anything that needs to get recorded through time should happen in here.
        'full' will dump a ton of stuff to the outdir, and is normally only really
        used for fancy plotting."""
        if full:
            self.lambda_through_time[self.current_time] = self._lambda()
            self.simulate_seqs()
            self.species_through_time[self.current_time] = self.species

        ## Every once in a while test to be sure our community is the same size
        ## as we think it should be.
        if not len(self.local_community) == self.paramsdict["K"]:
            msg = "k = {} r = {}".format(len(self.local_community),\
                                                len(set(self.local_community)))
            raise MESSError("  Community size violation - {}".format(msg))

        try:
            ## Always log size changes through time
            abunds = collections.Counter(self.local_community)
            LOGGER.debug("_log \n{}".format(self.local_info))
            LOGGER.debug("abunds \n{}".format(abunds))
            for species in self.local_info:
                self.local_info[species]["abundances_through_time"][self.current_time] = abunds[species]
        except Exception as inst:
            raise MESSError("Error in _log() - {}".format(inst))

    def __str__(self):
        return "<LocalCommunity {}: Shannon's Entropy {}>".format(self.name, shannon(self.get_abundances()))


    def prepopulate(self, quiet=False):
        LOGGER.debug("prepopulating local_community - {}".format(self))
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

            try:
                new_species, _ = self.region.get_most_abundant()
            except Exception as inst:
                raise MESSException("Error in prepopulate - {}".format(inst))

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
                                                     "post_colonization_migrants",\
                                                     "abundances_through_time"])
        self.local_info = self.local_info.fillna(0)
        for sp in self.local_info:
            self.local_info[sp] = [0, 0, OrderedDict()]

        self.founder_flags = [True] * len(self.local_community)

        if not quiet:
            print("    Initializing local community:")
            print("      N species = {}".format(len(set(self.local_community))))
            print("      N individuals = {}".format(len(self.local_community)))
        LOGGER.debug("Done prepopulating - {}".format(self))


    ## Not updated
    def death_step(self):
        ## Select the individual to die

        ##currently this will fail under volcanic model because the entire local community will go extinct
        if self.region.paramsdict["community_assembly_model"] == "filtering":
            death_Probability = 0
            reject = 0

            ## While the death probability is less than random uniform numeber between 0 and 1,
            ## keep selecting a new victime to potentially die
            while death_Probability < np.random.uniform(0,1):
                reject = reject + 1
                victim = random.choice(self.local_community)
                victim_trait = self.region.get_trait(victim)
                death_Probability = 1 - (np.exp(-((victim_trait - self.paramsdict["filtering_optimum"]) ** 2)/self.region.get_weight()))

            self.rejections.append(reject)


        if self.region.paramsdict["community_assembly_model"] == "competition":
            death_Probability = 0
            reject = 0
            mean_local_trait = self.region.get_trait_stats(self.local_community)[0]

            while death_Probability < np.random.uniform(0,1):
                reject = reject + 1
                victim = random.choice(self.local_community)
                victim_trait = self.region.get_trait(victim)
                death_Probability = (np.exp(-((victim_trait - mean_local_trait) ** 2)/self.region.get_weight()))

            self.rejections.append(reject)

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
                self.extinction_times.append(self.current_time - self.local_info[victim]["colonization_times"])
                vic_info = self.local_info.pop(victim)
                LOGGER.debug("Extinction victim info \n{}\n{}".format(victim, vic_info))
            except Exception as inst:
                LOGGER.debug(self.local_info)
                raise MESSError("Exception during recording extinction - {}".format(inst))
                pass
            ## If the invasive prematurely goes extinct just pick a new one
            if victim[0] == self.invasive:
                LOGGER.info("invasive went extinct")
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
        self.local_info[new_species] = [self.current_time,\
                                        0,\
                                        OrderedDict([(self.current_time,self.paramsdict["mig_clust_size"])])]

        return new_species


    def migrate_step(self):
        """ Allow multiple colonizations. In this case we return the sampled species
        as well as a bool reporting whether or not this is the first colonization of
        this species into the local community so that the coltime can be recorded.
        multiple colonizations of a species do not update coltime, but we record them
        for migration rate calculation."""

        new_species, _ = self.region.get_migrant()
        if new_species in self.local_community:
            ## This is a post-colonization migrant so record the event and tell downstream
            ## not to update the colonization time.
            self.local_info[new_species]["post_colonization_migrants"] += 1
        else:
            ## This is a new migrant so init local_info for it
            self.local_info[new_species] = [self.current_time,\
                                            0,\
                                            OrderedDict([(self.current_time,self.paramsdict["mig_clust_size"])])]

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
                ## Removed the if statement because multiple colonizations are always allowed
                #if self.region.paramsdict["allow_multiple_colonizations"]:
                new_species = self.migrate_step()
                #else:
                #    new_species = self.migrate_no_dupes_step()

                ## Only set the invasive species once at the time of next migration post invasion time
                ## If invasion time is < 0 this means "Don't do invasive"
                if not self.invasion_time < 0:
                    if self.invasive == -1 and self.current_time >= self.invasion_time:
                        self.invasive = new_species
                        LOGGER.info("setting invasive species {} at time {}".format(self.invasive, self.current_time))
                        self.invasion_time = self.current_time

                ## Add the colonizer to the local community, record the colonization time
                self.local_community.extend([new_species] * self.paramsdict["mig_clust_size"])
                self.founder_flags.extend([False] * self.paramsdict["mig_clust_size"])
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
        return SAD(self.local_community)


    ## How strong is the bottleneck? Strength should be interpreted as percent of local
    ## community to retain
    ## Not done
    def bottleneck(self, strength=1):
        reduction = int(round(self.paramsdict["K"] * strength))
        self.local_community = self.local_community[:reduction]

        ## First remove the extinct species from the species list
        pre = len(self.species)
        self.species = [s for s in self.species if s.uuid in self.local_community]
        ## Update the extinction counter
        self.extinctions += (pre - len(self.species))

        sp = self.species
        ## Update abundances per species that survived the bottleneck
        for i, s in enumerate(sp):
            if s.uuid in self.local_community:
                abund = self.local_community.count(s.uuid)
                s.update_abundance(abund)
                self.species[i] = s


    def simulate_seqs(self):
        self.species = []
        for name, coltime in self.local_info.loc["colonization_times"].iteritems():
            try:
                meta_abund = self.region.get_abundance(name)
                local_abund = self.local_community.count(name)
                tdiv = self.current_time - coltime
                sp = species(UUID=name,
                             colonization_time = tdiv,\
                             growth = self.region.paramsdict["population_growth"],\
                             abundance = local_abund,\
                             meta_abundance = meta_abund,
                             migration_rate = self.local_info[name]["post_colonization_migrants"]/float(tdiv),\
                             abundance_through_time = self.local_info[name]["abundances_through_time"].values())
                sp.simulate_seqs()
                sp.get_sumstats()
                self.species.append(sp)
                ## For debugging invasives
                #if s.abundance > 1000:
                #    print("\n{}".format(s))
            except Exception as inst:
                print(self.local_community)
                print(self.local_info)
                print(len(set(self.local_community)))
                print(self.local_info.shape)
                msg = "Error in simulate_seqs() - {}\nabundance - {} / meta_abundance {}\n{}\n{}".format(name,
                                                                                                         local_abund,
                                                                                                         meta_abund,
                                                                                                         self.local_info[name],
                                                                                                         inst)
                raise MESSError(msg)



    def get_stats(self):
        if self.rejections:
            LOGGER.debug("Average number of rejections - {}".format(np.mean(self.rejections)))

        LOGGER.debug("Entering get_stats()")
        self.simulate_seqs()
        LOGGER.debug("First 5 species - \n{}".format(self.species[:5]))
        self.stats._lambda = self._lambda()
        self.stats.generation = self.current_time
        self.stats.K = self.paramsdict["K"]
        self.stats.colrate = self.paramsdict["colrate"]
        self.stats.colrate_calculated = self.colonizations/float(self.current_time)
        self.stats.extrate_calculated = self.extinctions/float(self.current_time)
        self.stats.R = len(set(self.local_community))
        self.stats.shannon = shannon(self.get_abundances(octaves=False))

        pis = np.array([x.stats["pi_local"] for x in self.species])
        dxys = np.array([x.stats["dxy"] for x in self.species])
        self.stats.mean_pi = np.mean(pis)
        self.stats.stdv_pi = np.std(pis)
        self.stats.median_pi = np.median(pis)
        self.stats.iqr_pi = iqr(pis)
        self.stats.mean_dxy = np.mean(dxys)
        self.stats.stdv_dxy = np.std(dxys)
        self.stats.median_dxy= np.median(dxys)
        self.stats.iqr_dxy = iqr(dxys)

        self.stats.sgd = SGD(pis, dxys)

        self.stats.mean_ltr = self.region.get_trait_stats(self.local_community)[0]
        self.stats.var_ltr = self.region.get_trait_stats(self.local_community)[1]
        self.stats.mean_rtr = self.region.get_trait_stats(self.local_community)[2]
        self.stats.var_rtr = self.region.get_trait_stats(self.local_community)[3]
        self.stats.mean_dif = self.region.get_trait_stats(self.local_community)[4]
        self.stats.var_dif = self.region.get_trait_stats(self.local_community)[5]
        self.stats.kurtosis = self.region.get_trait_stats(self.local_community)[6][0]
        self.stats.skewness = self.region.get_trait_stats(self.local_community)[7][0]

        ## Log to file
        #statsfile = os.path.join(self._hackersonly["outdir"],
        #                         self.paramsdict["name"] + "-simout.txt")
        #self.stats.to_csv(statsfile, na_rep=0, float_format='%.5f')

        #megalog = os.path.join(self._hackersonly["outdir"],
        #                         self.paramsdict["name"] + "-megalog.txt")
        ## concatenate all species results and transpose the data frame so rows are species
        fullstats = pd.concat([sp.stats for sp in self.species], axis=1).T
        #fullstats.to_csv(megalog, index_label=False)

        ## If you ever want to pretty print results
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
    "filtering_optimum" : "optimum trait value, only used during environmental filtering model",\
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
    loc = LocalCommunity("wat", K=5000)
    data._link_local(loc)
    print(loc)
    print(loc.local_info)
    ## Allow for either having or not having empty demes
    assert(len(collections.Counter(loc.local_community)) <= 2)
    loc.paramsdict["mode"] = "landbridge"
    loc.prepopulate()
    assert(len(collections.Counter(loc.local_community)) > 3)
    print(loc.get_abundances())

    loc.paramsdict["mode"] = "volcanic"
    loc.prepopulate()
    loc.step(100000)
    print(len(set(loc.local_community)))
    print(len(loc.local_community))
    print(loc.local_info.shape)
    print("getting stats")
    print(loc)
    print(loc.paramsdict)
