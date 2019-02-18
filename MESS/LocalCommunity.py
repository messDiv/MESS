#!/usr/bin/env python2.7

from __future__ import print_function

from scipy.stats import logser
from collections import OrderedDict
from scipy.stats import iqr,hmean
import collections
import pandas as pd
import numpy as np
import msprime
import itertools
import random
import sys
import os
import MESS

from .util import MESSError, _tuplecheck, sample_param_range, memoize
from .stats import shannon, SAD
from .species import species
from .SGD import SGD

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
                        ("mig_clust_size", mig_clust_size)
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
                        ("allow_empty", True),
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
        ##
        ## NB: If you add another field you need to change the initialization inside
        ##     prepopulate, and also the behavior in _add_local_info()
        ##
        ## Fields:
        ##  * "colonization_times" - Colonization times per species
        ##  * "post_colonization_migrants" - Count of post colonization migrants per species
        ##  * "abundances_through_time" - A record of abundance through time for this species
        ##  * "ancestor" - The species name of the immediate ancestor if local speciation
        ##                 happened, otherwise "".
        self.local_info = pd.DataFrame([])

        ## The regional pool that this local community belongs to
        ## this is updated by Region._link_local(), so don't set it by hand
        self.region = ""

        self.files = dict({
                "full_output": [],
                })

        ## summary stats dict
        self.stats = pd.Series(
            index = ["_lambda",
                   "generation",
                   "K",
                   "colrate",
                   "speciation_probability",
                   "sigma",
                   "trait_rate",
                   "ecological_strength",
                   "filtering_optimum",
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
                   "trees",
                   "mn_local_traits",
                   "var_local_traits",
                   "mn_regional_traits",
                   "var_regional_traits",
                   "reg_loc_mn_trait_dif",
                   "reg_loc_var_trait_dif",
                   "kurtosis_local_traits",
                   "skewness_local_traits"]).astype(np.object)

        ## Will be initialized as a pd series when the region is linked
        self.SGD = ""

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


    def _add_local_info(self, sname, abundances_through_time=0,\
                        ancestor='', ancestral_abundance=[], speciation_completion=0):
        ## Construct a new local_info record for new_species. The fields are:
        ## colonization time - in the forward time model. This gets converted to
        ##      divergence time for the backward time model.
        ## post_colonization_migrants - The count of migrants that have come in
        ##      that are the same species as this one, since colonization
        ## abundances_through_time - Dictionary containing history of population
        ##      size change through time. Default is 0, which indicates a new colonist.
        ## ancestor - If the species was introduced by speciation rather than
        ##      colonization, then it'll have an ancestral species.
        #######################
        ## Protracted speciation parameters that shouldn't be touched otherwise
        #######################
        ## ancestral_abundance - A list of fluctuating ancestral abundances at
        ##      the time the species split from its sister. Default is empty list which indicates
        ##      a new colonist from the metacommunity.
        ## speciation_completion - The point in forward time when this species will become 'good'.
        ##      before this point it is an incipient species and if the simulation is stopped then
        ##      all the individuals will be thrown back into the parent species pool.
        ##      Default is 0 which indicates this is a good species immediately, either a new
        ##      colonizing lineage or a point mutation species.
        if abundances_through_time == 0: abundances_through_time = OrderedDict([(self.current_time,self.paramsdict["mig_clust_size"])])
        self.local_info[sname] = [self.current_time,\
                                        0,\
                                        abundances_through_time,\
                                        ancestor,\
                                        ancestral_abundance,
                                        speciation_completion]


    def _set_region(self, region):
        self.region = region
        self.SGD = SGD([], ndims=region.paramsdict["sgd_dimensions"], nbins=region.paramsdict["sgd_bins"])


    def _get_stats_header(self):
        return pd.concat([self.stats, self.SGD.to_series()]).keys()


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
                self.paramsdict[param] = newvalue

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

            for key, val in self.paramsdict.items():
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)

                padding = (" "*(20-len(paramvalue)))
                paramkey = list(self.paramsdict.keys()).index(key)
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
            LOGGER.debug("_log - {} lambda {} gen {}\n{}".format(self.name, self._lambda(), self.current_time, self.local_info))
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
                                                     "abundances_through_time",\
                                                     "ancestor",\
                                                     "ancestral_abundance",\
                                                     "speciation_completion"])
        self.local_info = self.local_info.fillna(0)
        for sp in self.local_info:
            self.local_info[sp] = [0, 0, OrderedDict(), "", [], 0]

        self.founder_flags = [True] * len(self.local_community)
        if not quiet:
            print("    Initializing local community:")
            print("      N species = {}".format(len(set(self.local_community))))
            print("      N individuals = {}".format(len(self.local_community)))
        LOGGER.debug("Done prepopulating - {}".format(self))


    def death_step(self):

        ## Select the individual to die
        victim = random.choice(self.local_community)
        done = False

        ## Emtpy niche always dies and neutral model always accepts selection regardless
        if victim == None or self.region.paramsdict["community_assembly_model"] == "neutral":
            done = True

        reject = 0
        survival_scalar = 0.25
        while not done:
            ## If you made it this far then you're doing a trait model.

            ## If reject is > 0, then this is a trait model that has chosen an individual
            ## at time 0, so we will keep looping here and sampling only real individuals
            ## Saves a bunch of time looping and grabbing None over and over, esp in early days
            if reject > 0:
                victim = random.choice([x for x in self.local_community if x != None])

            ## This is the "get on with it" switch. Only for trait models, especially
            ## for when the competition model goes off the rails. Here ratchet up the
            ## allowance for death so it doesn't just sit there spinning forever
            if reject % 5 == 0 and reject > 0:
                survival_scalar += .2
                LOGGER.debug("Survival scalar bump - {}".format(survival_scalar))

            death_thresh = np.random.uniform(0,1)
            victim_trait = self.region.get_trait(victim)

            if self.region.paramsdict["community_assembly_model"] == "filtering":

                ## Call to _get_filter is memoized so results are cached
                death_probability = _get_filtering_death_prob(self.region, victim_trait)
                death_probability = (1 - death_probability) * survival_scalar + death_probability
                target_trait_val = self.region.metacommunity.paramsdict["filtering_optimum"]

            elif self.region.paramsdict["community_assembly_model"] == "competition":
                mean_local_trait = self.region.get_trait_stats(self.local_community, mean_only=True)
                death_probability = _get_competition_death_prob(self.region, victim_trait, mean_local_trait)
                #death_probability = (np.exp(-((victim_trait - mean_local_trait) ** 2)/self.region.metacommunity.paramsdict["ecological_strength"]))
                death_probability = (1 - death_probability) * survival_scalar + death_probability
                target_trait_val = mean_local_trait

            LOGGER.debug("rj {} trait {} dprob {} dthr {} target {}".format(reject, victim_trait, death_probability, death_thresh, target_trait_val))
            if death_probability > death_thresh:
                done = True
            else:
                reject = reject + 1

        self.rejections.append(reject)

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
            ## Clean up local community list and founder flag list
            idx = self.local_community.index(victim)
            self.local_community.pop(idx)
            self.founder_flags.pop(idx)

            ## If the species of the victim went locally extinct then clean up local_info
            self._test_local_extinction(victim)
        except Exception as inst:
            import pdb; pdb.set_trace()
            raise inst


    def _test_local_extinction(self, victim):
        ## Record local extinction events
        ancestor = victim
        ## Don't clean up after emtpy niche space
        if not victim in self.local_community and victim != None:
            self.extinctions += 1
            try:
                ## Record the lifetime of this species and remove their record from divergence_times
                ## Remove the species from the local_info array
                self.extinction_times.append(self.current_time - self.local_info[victim]["colonization_times"])
                vic_info = self.local_info.pop(victim)
                offspring = self.local_info.columns[self.local_info.loc["ancestor"] == victim]
                #LOGGER.debug("\nExtinction victim info \n{}\n{}\nOffspring {}".format(victim, vic_info, offspring))

                ## Speciation process nonsense
                ## Update ancestry and population size change history for any species with this one as
                ## direct ancestor
                ancestor = vic_info["ancestor"]
                anc_size_changes = vic_info["abundances_through_time"]
                self.local_info.loc["ancestor"][self.local_info.loc["ancestor"] == victim] = ancestor
                ## I don't think you can vectorize the update

                for o in offspring:
                    LOGGER.debug("offspring {} {}".format(o, self.local_info[o]["abundances_through_time"]))
                    self.local_info[o]["abundances_through_time"].update(anc_size_changes)
                    LOGGER.debug("offspring {} {}".format(o, self.local_info[o]["abundances_through_time"]))

            except Exception as inst:
                LOGGER.debug(self.local_info)
                raise MESSError("Exception during recording extinction - {}".format(inst))
            ## If the invasive prematurely goes extinct just pick a new one
            if victim[0] == self.invasive:
                LOGGER.info("invasive went extinct")
                self.invasive = -1
        return ancestor


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
        self._add_local_info(sname = new_species)

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
            self._add_local_info(sname = new_species)

        return new_species


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
                try:
                    self.death_step()
                    ## Sample all available from local community (community grows slow in volcanic model)
                    ## This is the fastest way to sample from a list. >4x faster than np.random.choice
                    chx = random.choice(self.local_community)
                    self.local_community.append(chx)
                    idx = self.local_community.index(chx)
                    self.founder_flags.append(self.founder_flags[idx])
                except Exception as inst:
                    import pdb; pdb.set_trace()
                    LOGGER.error("Exception in step() - {}".format(inst))
                    raise inst

            ##############################################
            ## Speciation process
            ##############################################
            if self.region.paramsdict["speciation_model"] != "none" and\
               np.random.random_sample() < self.region.paramsdict["speciation_probability"]:

               self.speciate()

            ## update current time
            self.current_time += 1


    def get_abundances(self, octaves=False):
        return SAD(self.local_community)


    def speciate(self):
        """ Occassionally initiate the speciation process. In all modes, one
            one individual is randomly selected to undergo speciation.
            Speciation does not change the founder_flag state.

            Currently there are 3 modes implemented:
            - point_mutation: The randomly selected individual becomes a new
                species of its own, of abundance 1.
            - random_fission: The species of the randomly selected individual
                undergoes random fission. In this mode the abundance of the
                new species is determined by randomly splitting off a chunk
                of the individuals from the parent species. All fission
                sizes are equally likely.
            - protracted:
        """
        LOGGER.debug("Initiate speciation process")

        ## Sample the individual to undergo speciation
        chx = random.choice(self.local_community)
        idx = self.local_community.index(chx)

        ## Construct the new species name.
        ## We want the new name to be globally unique but we don't
        ## want to waste a bunch of time keeping track of it in the
        ## region, so we can do something like this:
        sname = chx + ":{}-{}".format(self.name, self.current_time)

        ## Fetch the abundance history of the parent species
        parent_abunds = self.local_info[chx]["abundances_through_time"]

        ## Inform the regional pool that we have a new species
        trt = self.region.get_trait(chx)
        self.region._record_local_speciation(sname, trt)

        if self.region.paramsdict["speciation_model"] == "point_mutation":

            ## Replace the individual in the local_community with the new species
            self.local_community[idx] = sname

            ## If the new individual removes the last member of the ancestor
            ## species, then you need to do some housekeeping.
            ## TODO: is this really an "extinction" event? we need to clean up
            ## the local_info regardless, but we may not want to increment the
            ## extinctions counter.
            ancestor = self._test_local_extinction(chx)

            ## In the point mutation speciation model the offspring population
            ## immediately inherits the parent population history of size change.
            ## _test_local_extinction() handles all the ancestor inheritence logic.
            self._add_local_info(sname = sname, abundances_through_time=parent_abunds , ancestor = ancestor)

        elif self.region.paramsdict["speciation_model"] == "random_fission":

            ## Remove all individuals of the target species from the local community.
            ## We'll have to add some back as the original species, but doing it
            ## this way simplifies the housekeeping.
            self.local_community = [sp for sp in self.local_community if sp != chx]
            ## Get abundance of the target species so we can perform
            ## the random fission.
            sp_abund = self.paramsdict["K"] - len(self.local_community)

            ## Get the number of individuals to carve off for the new species.
            ## If sp_abund == 1, or if new_abund == sp_abund then this is
            ## essentially anagenetic speciation, as the initial species will
            ## be removed from the local community and replaced by the new sp.
            new_abund = np.random.randint(1, sp_abund)

            ## Push the appropriate number of individuals of the original species
            ## back into the local community. Don't bother testing for sp_abund > 0
            ## here since the extend will happily ignore the empty list.
            sp_abund = sp_abund - new_abund
            self.local_community.extend([chx] * sp_abund)

            ## Test local extinction after parent species is pushed back to
            ## local community. If no parent species remain at this time then
            ## we need to do some housekeeping to track the ancestor of the
            ## new species.
            ancestor = self._test_local_extinction(chx)

            ## Push the new appropriate number of new individuals of the new
            ## species into the local community
            self.local_community.extend([sname] * new_abund)

            ## Same as for point mutation, the new species inherits the abndance
            ## history of the parent, and also _test_local_extinction() handles
            ## all the ancestor inheritence logic.
            self._add_local_info(sname = sname, abundances_through_time=parent_abunds , ancestor = ancestor)

        elif self.region.paramsdict["speciation_model"] == "protracted":
            pass
        else:
            raise MESSError("Unrecognized speciation model - {}".format(self.region.paramsdict["speciation_model"]))


    ## TODO: Unused and not updated to current MESS structure
    ## How strong is the bottleneck? Strength should be interpreted as percent of local
    ## community to retain
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


    ################################################
    ## Functions for driving the backward time model
    ################################################
    ## I always forget how to do this in pandas: self.local_info.iloc[:, idxs]
    def _get_singleton_species(self):
        ## A function to return any species that are not involved in a speciation
        ## event. These are "easier" to handle individually, so we'll do them separately.
        ## First get the candidate list of species.
        idxs = np.where(self.local_info.loc["ancestor"] == "")[0]
        names = self.local_info.columns[idxs]
        ## Get just the species that don't occur as ancestors in the local community
        singletons = set(names).difference(set(self.local_info.loc["ancestor"]))
        ## Return the dataframe with just those individuals without
        return self.local_info[list(singletons)]


    def _get_clades(self):
        ## We need to get all the groups of species that descended from a common
        ## local anscestor so we can run the backwards time model for all of them
        ## combined in a coherent fashion.

        ## Sort them by colonization time
        self.local_info.sort_values("colonization_times", axis=1, inplace=True)

        ## Clades will all be descended from one unique ancestor in the metacommunity
        ## so we will just group them by this identifier. Species names have this
        ## format: t0:<local_name>-<split_time>, so the zeroth element of the split
        ## on ':' will always be the metacommunity species id.
        clades = {x:[] for x in set([y.split(":")[0] for y in self.local_info])}

        for idx in self.local_info:
            clades[idx.split(":")[0]].append(idx)

        LOGGER.debug("Clades - {}".format(clades))

        return clades


    def sim_seqs(self):
        self.species = []
        try:
            ## Hax. Remove the empty deme from local info. This _might_ break the fancy plots.
            self.local_community = [x for x in self.local_community if x != None]
            self.local_info = self.local_info[self.local_info.columns.dropna()]
        except:
            ## df.drop will raise if it doesn't find a matching label to drop, in which case we're done.
            pass

        self.local_info.loc["meta_abund"] = [self.region.get_abundance(x) for x in self.local_info]
        self.local_info.loc["local_abund"] = [self.local_community.count(x) for x in self.local_info]
        self.local_info.loc["colonization_times"] = self.current_time - self.local_info.loc["colonization_times"]
        for cname, species_list in self._get_clades().items():
            ## Just get a slice of the local_info df that represents the species of interest
            dat = self.local_info[species_list]
            ## Dictionary mapping species names to 0-based index (where 0 is metacommunity)
            sp_idxs = OrderedDict({x:ix for ix, x in zip(range(1, len(dat.columns)+1), dat.columns[::-1])})
            sp_idxs[''] = 0

            pop_cfgs = []
            pop_meta = msprime.PopulationConfiguration(sample_size = 10, initial_size = 10000)
            pop_cfgs.append(pop_meta)
            for sp, idx in sp_idxs.items():
                if not sp:
                    continue
                sizechange_times = sorted(dat[sp]["abundances_through_time"], reverse=True)

                size = hmean([dat[sp]["abundances_through_time"][x]*1000 for x in sizechange_times])
                pop_local = msprime.PopulationConfiguration(sample_size = 10, initial_size = size, growth_rate = 0)
                pop_cfgs.append(pop_local)

            ## sp are added in chronological order of coltime, so the'll be in the right order here for
            ## adding in reverse coltime
            split_events = []
            for col in dat.columns:
                #print(col, dat[col].loc['colonization_times'])
                #print(sp_idxs[col])
                split_event = msprime.MassMigration(time = dat[col].loc['colonization_times'],\
                                                    source = sp_idxs[col],\
                                                    destination = sp_idxs[dat[col]["ancestor"]],\
                                                    proportion = 1)
                split_events.append(split_event)

            try:
                #print(sp_idxs, split_events)
                debug = msprime.DemographyDebugger(population_configurations = pop_cfgs,\
                                                 demographic_events = split_events)

                ## Enable this at your own peril, it will dump a ton of shit to stdout
                #debug.print_history()
            except:
                raise

            tree_sequence = msprime.simulate(length = 600,\
                                            mutation_rate = 1e-7,\
                                            population_configurations = pop_cfgs,\
                                            demographic_events = split_events)

            ## Now we have the tree sequence for the clade and can go back through the species list
            ## and pull out stats per species
            LOGGER.debug("Done with hsitory now simulating seqs")
            metasamps = tree_sequence.samples(population=0)
            for sp, idx in sp_idxs.items():
                ## Skip the metacommunity, which doesn't have a parent :'(
                if sp == '':
                    continue
                try:
                    migration_rate = dat[sp]["post_colonization_migrants"]/float(dat[col]['colonization_times'])
                except ZeroDivisionError as inst:
                    ## This should only happen when coltime is 0, which should be never
                    LOGGER.error("Got bad coltime - {}".format(dat[col]))
                    migration_rate = 0

                sp_obj = species(name = sp,
                         species_params = self.region.get_species_params(),
                         divergence_time = dat[sp].loc['colonization_times'],\
                         growth = self.region.paramsdict["population_growth"],\
                         abundance = dat[sp].loc["local_abund"],\
                         meta_abundance = dat[sp]["meta_abund"],
                         migration_rate = migration_rate,\
                         abundance_through_time = {self.current_time - x:y for x, y in list(dat[sp]["abundances_through_time"].items())})
                sp_obj.tree_sequence = tree_sequence
                samps = tree_sequence.samples(population=idx)
                sp_obj.get_sumstats(samps, metasamps)
                self.species.append(sp_obj)


    def simulate_seqs(self):
        self.sim_seqs()
        ##FIXME Oh boy. This was hacked for the speciation code. Gotta clean this up.
        return
        self.species = []
#        for name, coltime in self._get_singleton_species().loc["colonization_times"].items():
        for name, coltime in self.local_info.loc["colonization_times"].iteritems():
            try:
                ## Get the meta abundance for the original colonizing lineage
                meta_abund = self.region.get_abundance(name.split(":")[0])
                local_abund = self.local_community.count(name)
                tdiv = self.current_time - coltime
                tdiv = tdiv / float(self.paramsdict["K"])
                ## Rescale abundances through time so they are "backwards" values
                abundances_through_time = {self.current_time - x:y for x, y in list(self.local_info[name]["abundances_through_time"].items())}
                try:
                    sp = species(name = name,
                                 species_params = self.region.get_species_params(),
                                 divergence_time = tdiv,\
                                 growth = self.region.paramsdict["population_growth"],\
                                 abundance = local_abund,\
                                 meta_abundance = meta_abund,
                                 migration_rate = self.local_info[name]["post_colonization_migrants"]/float(tdiv),\
                                 abundance_through_time = abundances_through_time)
                except Exception as inst:
                    print("Error in creating species {}\n{}".format(name, inst))
                sp.simulate_seqs()
                self.species.append(sp)
                ## For debugging invasives
                #if s.abundance > 1000:
                #    print("\n{}".format(s))
            except Exception as inst:
                print("tdiv = {}".format(tdiv))
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
        ## Model parameters
        self.stats._lambda = self._lambda()
        self.stats.generation = self.current_time
        self.stats.K = self.paramsdict["K"]
        self.stats.colrate = self.paramsdict["colrate"]
        self.stats.speciation_probability = self.region.paramsdict["speciation_probability"]
        self.stats.sigma = self.region.paramsdict["sigma"]
        self.stats.trait_rate = self.region.metacommunity.paramsdict["trait_rate"]
        self.stats.ecological_strength = self.region.metacommunity.paramsdict["ecological_strength"]
        ## Pseudo-parameters
        self.stats.filtering_optimum = self.region.metacommunity.paramsdict["filtering_optimum"]
        try:
            self.stats.colrate_calculated = self.colonizations/float(self.current_time)
            self.stats.extrate_calculated = self.extinctions/float(self.current_time)
        except ZeroDivisionError as inst:
            raise MESSError("Current time should never be zero, check parameter settings.")
        ## Model sumstats
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

        self.SGD = SGD(pis,\
                       dxys,\
                       nbins = self.region.paramsdict["sgd_bins"],\
                       ndims = self.region.paramsdict["sgd_dimensions"])
        LOGGER.debug("SGD - {}".format(self.SGD))

        try:
            trait_stats = self.region.get_trait_stats(self.local_community)
            self.stats.mn_local_traits = trait_stats[0]
            self.stats.var_local_traits = trait_stats[1]
            self.stats.mn_regional_traits = trait_stats[2]
            self.stats.var_regional_traits = trait_stats[3]
            self.stats.reg_loc_mn_trait_dif = trait_stats[4]
            self.stats.reg_loc_var_trait_dif = trait_stats[5]
            self.stats.kurtosis_local_traits = trait_stats[6]
            self.stats.skewness_local_traits = trait_stats[7]

            ## This line will get you phy stats (mean and var of branch lenghts) from the regional/metacommuntiy tree
            ## Not sure what to do about local tree yet..
            #phy_stats = self.region.get_phy_stats(self.region.metacommunity.metacommunity_tree)

            ## Log to file
            #statsfile = os.path.join(self._hackersonly["outdir"],
            #                         self.paramsdict["name"] + "-simout.txt")
            #self.stats.to_csv(statsfile, na_rep=0, float_format='%.5f')

            ## If you don't actually want all the intermediate files then we won't make them
            if self.region._log_files:
                megalog = os.path.join(self._hackersonly["outdir"],
                                         self.paramsdict["name"] + "-{}-megalog.txt".format(self._lambda()))
                ## concatenate all species results and transpose the data frame so rows are species
                fullstats = pd.concat([sp.stats for sp in self.species], axis=1).T
                fullstats.to_csv(megalog, index_label=False)
        except Exception as inst:
            import pdb; pdb.set_trace()
            LOGGER.error("Error in get_stats() - {}".format(inst))
            raise

        return pd.concat([self.stats, self.SGD.to_series()])


@memoize
def _get_filtering_death_prob(region, victim_trait):
    try:
        val = 1 - (np.exp(-((victim_trait - region.metacommunity.paramsdict["filtering_optimum"]) ** 2)/region.metacommunity.paramsdict["ecological_strength"]))
    except Exception as inst:
        raise MESSError("Error getting death prob using trait - {}".format(victim_trait))
    return val


@memoize
def _get_competition_death_prob(region, victim_trait, mean_local_trait):
    try:
        val = death_probability = (np.exp(-((victim_trait - mean_local_trait) ** 2)/region.metacommunity.paramsdict["ecological_strength"]))
    except Exception as inst:
        raise MESSError("Error in geting death prob using trait & mean - {} {}".format(victim_trait, mean_local_trait))
    return val


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
