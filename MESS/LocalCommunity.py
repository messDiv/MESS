#!/usr/bin/env python2.7

from __future__ import print_function

from scipy.stats import logser
from collections import OrderedDict
import collections
import pandas as pd
import numpy as np
import msprime
import itertools
import random
import sys
import os
import MESS

from .util import MESSError, tuplecheck, sample_param_range, memoize
from .stats import *
from .species import species
from .SGD import SGD

import logging
LOGGER = logging.getLogger(__name__)

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1500


class LocalCommunity(object):

    def __init__(self, name="Loc1", K=1000, colrate=0.01, quiet=False):
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
                        ("K", K),
                        ("colrate", colrate),
                        ("speciation_rate", 0),
                        ("background_death", 0.25)
        ])

        ## A dictionary for holding prior ranges for values we're interested in
        self._priors = dict([
                        ("K", []),
                        ("colrate", []),
        ])

        ## Dictionary of 'secret' parameters that most people won't want to mess with
        ##  * allow_empty is a switch on whether to fully populate local community
        ##      with one species for 'volcanic' mode or to introduce just one
        ##      individual and populate the rest of K with 'empty' demes.
        ##  * outdir is inherited from the Region.simulate() command, so users
        ##      should normally not mess with this.
        ##  * mig_clust_size: number of incoming migrants during colonization events.
        ##      The number of individuals sampled to die is adjusted to equal this.
        ##  * age: This was intended to be the time at which the "island" arises
        ##      during the simulation, to allow for islands popping up at
        ##      at different times, but I never implemented this yet.
        ##  * trait_rate_local: The trait evolution rate parameter for local community,
        ##      You can fix this value to something different if you really want,
        ##      but by default we calculate from the trait rate meta divided by global
        ##      birth rate + death rate.
        ##  * mode: Whether to prepopulate the local community as a 'volcanic' or
        ##      'landbridge' style origin.
        self._hackersonly = dict([
                        ("allow_empty", True),
                        ("outdir", []),
                        ("mig_clust_size", 1),
                        ("age", 100000),
                        ("trait_rate_local", 0),
                        ("mode", "volcanic"),
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
        for param in ["K", "colrate"]:
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
        if abundances_through_time == 0: abundances_through_time = OrderedDict([(self.current_time,self._hackersonly["mig_clust_size"])])
        self.local_info[sname] = [self.current_time,\
                                        0,\
                                        abundances_through_time,\
                                        ancestor,\
                                        ancestral_abundance,
                                        speciation_completion]


    def _set_region(self, region):
        self.region = region
        self.SGD = SGD([], ndims=region._hackersonly["sgd_dimensions"], nbins=region._hackersonly["sgd_bins"])
        self._hackersonly["trait_rate_local"] = _get_trait_rate_local(self.region)


    ## Getting params header and parameter values drops the local
    ## community name (param 0), and adds a bunch of pseudo-parameters
    def _get_params_header(self):
        params_header = list(self.paramsdict.keys())[1:]
        params_header = params_header + ["generation", "_lambda", "colrate_calculated", "extrate_calculated",\
                                            "trait_rate_local", "filtering_optimum"]
        return params_header


    def _get_params_values(self):
        ## Get params and drop name
        params_vals = list(self.paramsdict.values())[1:]
        ## We are reporting generations scaled to WF time
        params_vals = params_vals + [self.current_time * 2 / self.paramsdict["K"],\
                                    self._lambda(),\
                                    self.colonizations/float(self.current_time),\
                                    self.extinctions/float(self.current_time),\
                                    self._hackersonly["trait_rate_local"],\
                                    self.region.metacommunity._hackersonly["filtering_optimum"]]
        params_vals = pd.DataFrame(params_vals, index=self._get_params_header())
        return params_vals


    def _paramschecker(self, param, newvalue, quiet=False):
        """ Raises exceptions when params are set to values they should not be"""
        ## TODO: This should actually check the values and make sure they make sense
        ## TODO: Also check here if you're setting the mode parameter you have to rerun prepopulate
        try:

            ## Cast params to correct types
            if param in ["K"]:
                tup = tuplecheck(newvalue, dtype=int)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup)[0]
                else:
                    self.paramsdict[param] = tup

            elif param in ["colrate"]:
                tup = tuplecheck(newvalue, dtype=float)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup)[0]
                else:
                    self.paramsdict[param] = tup

            elif param == "mode":
                ## Must reup the local community if you change the mode
                self.paramsdict[param] = newvalue

            elif param == "speciation_rate":
                self.paramsdict[param] = float(newvalue)

            elif param == "background_death":
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
        return "<LocalCommunity {}>".format(self.name)


    def prepopulate(self, verbose=False):
        LOGGER.debug("prepopulating local_community - {}".format(self))
        if not self.region:
            msg = "Skip populating the local community as it is unlinked to a region."
            LOGGER.error(msg)
            if MESS.__interactive__: print("    {}".format(msg))
            return

        ## Clean up local_community if it already exists
        self.local_community = []

        if self._hackersonly["mode"] == "landbridge":
            ## prepopulate the island w/ a random sample from the metacommunity
            ## TODO: The underscore here is ignoring trait values
            self.local_community, _ = self.region.get_nmigrants(self.paramsdict["K"])

        elif self._hackersonly["mode"]  == "volcanic":
            ## If not landbridge then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            ## TODO: The _ is a standin for trait values, have to do something with them

            try:
                new_species, _ = self.region.get_most_abundant()
            except Exception as inst:
                raise MESSError("Error in prepopulate - {}".format(inst))

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
        if verbose:
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
        survival_scalar = self.paramsdict["background_death"]
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
                survival_scalar += .05
                LOGGER.debug("Survival scalar bump - {}".format(survival_scalar))

            death_thresh = np.random.uniform(0,1)
            victim_trait = self.region.get_trait(victim)

            ## TODO: There's a bunch of nonsense here that could be optimized.
            ## The trait based stuff is so much slower than the neutral, and
            ## the competition model is especially bad, so optimization could
            ## be useful.
            if self.region.paramsdict["community_assembly_model"] == "filtering":

                ## Call to _get_filter is memoized so results are cached
                death_probability = _get_filtering_death_prob(self.region, victim_trait)
                death_probability = (1 - death_probability) * survival_scalar + death_probability
                target_trait_val = self.region.metacommunity._hackersonly["filtering_optimum"]

            elif self.region.paramsdict["community_assembly_model"] == "competition":
                mean_local_trait = self.region.get_trait_mean(self.local_community)
                death_probability = _get_competition_death_prob(self.region, victim_trait, mean_local_trait)
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
            #import pdb; pdb.set_trace()
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
        ## Convert time in generations to timesteps (WF -> Moran)
        for step in range(nsteps * self.paramsdict["K"] / 2):
            ## Check probability of an immigration event
            if np.random.random_sample() < self.paramsdict["colrate"]:
                ## If clustered migration remove the necessary number of additional individuals
                for _ in range(self._hackersonly["mig_clust_size"]):
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
                self.local_community.extend([new_species] * self._hackersonly["mig_clust_size"])
                self.founder_flags.extend([False] * self._hackersonly["mig_clust_size"])
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
                    #import pdb; pdb.set_trace()
                    LOGGER.error("Exception in step() - {}".format(inst))
                    raise inst

            ##############################################
            ## Speciation process
            ##############################################
            if self.region.paramsdict["speciation_model"] != "none" and\
               np.random.random_sample() < self.paramsdict["speciation_rate"]:

               self.speciate()

            ## update current time
            self.current_time += 1


    def get_abundances(self, octaves=False, raw_abunds=False):
        return SAD(self.local_community, octaves=octaves, raw_abunds=raw_abunds)


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

        ## Sample the individual to undergo speciation. Remove all
        ## empty deme space prior to this, if it exists, since we don't
        ## want empty deme space speciating.
        chx = random.choice([sp for sp in self.local_community if sp != None])
        idx = self.local_community.index(chx)


        ## Construct the new species name.
        ## We want the new name to be globally unique but we don't
        ## want to waste a bunch of time keeping track of it in the
        ## region, so we can do something like this:
        sname = chx + ":{}-{}".format(self.name, self.current_time)

        ## Fetch the abundance history of the parent species
        parent_abunds = self.local_info[chx]["abundances_through_time"]

        ## Inform the regional pool that we have a new species
        ## Identify parent's trait value
        parent_trait = self.region.get_trait(chx)

        ## Trait evolution. Offsprint trait is normally distributed
        ## with mean of parent value, and stdv equal to stdv of BM
        ## process in metacommunity times average lineage lifetime
        trt = np.random.normal(parent_trait, self._hackersonly["trait_rate_local"], 1)

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
            ## The `sp_abund+1` here is because randint samples up to sp_abund-1,
            ## so we need to allow for the case of new_abund == sp_abund.
            new_abund = np.random.randint(1, sp_abund+1)

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
        ## TODO: Unused. All lineages are handled with the _get_clades() function.
        ##       Should just remove this.
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
            split_events = []
            meta_abund = self.region.get_abundance(cname)
            pop_meta = msprime.PopulationConfiguration(sample_size = 10, initial_size = meta_abund)
            pop_cfgs.append(pop_meta)
            species_dict = {}
            for sp, idx in sp_idxs.items():
                if not sp:
                    continue
                try:
                    migration_rate = dat[sp]["post_colonization_migrants"]/float(dat[sp]['colonization_times'])
                except ZeroDivisionError as inst:
                    ## This should only happen when coltime is 0, which should be never
                    LOGGER.error("Got bad coltime - {}".format(dat[col]))
                    migration_rate = 0

                sp_obj = species(name = sp,
                         species_params = self.region.get_species_params(),
                         trait_value = self.region.get_trait(sp),
                         divergence_time = dat[sp].loc['colonization_times'],\
                         growth = self.region._hackersonly["population_growth"],\
                         abundance = dat[sp].loc["local_abund"],\
                         meta_abundance = dat[sp]["meta_abund"],
                         migration_rate = migration_rate,\
                         abundance_through_time = {self.current_time - x:y for x, y in list(dat[sp]["abundances_through_time"].items())})

                pop_local = sp_obj._get_local_configuration()
                pop_cfgs.append(pop_local)
                split_event = sp_obj._get_local_meta_split(source_idx = idx, dest_idx = sp_idxs[dat[sp]["ancestor"]])
                split_events.extend(split_event)
                species_dict[idx] = sp_obj

            ## Make sure demographic events are sorted increasing backward in time
            split_events.sort(key = lambda x: x.time)
            try:
                #print("\n", sp_idxs, split_events)
                debug = msprime.DemographyDebugger(population_configurations = pop_cfgs,\
                                                 demographic_events = split_events)

                ## Enable this at your own peril, it will dump a ton of shit to stdout
                #debug.print_history()
            except:
                raise

            tree_sequence = msprime.simulate(length = self.region.paramsdict["sequence_length"],\
                                            mutation_rate = self.region.paramsdict["mutation_rate"],\
                                            population_configurations = pop_cfgs,\
                                            demographic_events = split_events)

            ## This block of code is for debugging the msprime demography
            ## It might be cute to add a command line flag to optionally
            ## save these somewhere smart, but.... 
            if True:
                tree = tree_sequence.first()
                colour_map = {0:"red", 1:"blue"}
                for idx in range(len(sp_idxs) - 2):
                    r, g, b = np.random.randint(0, 255, 3)
                    colour_map[idx+2] = "rgb({}, {}, {})".format(r, g, b)
                node_colours = {u: colour_map[tree.population(u)] for u in tree.nodes()}
                node_labels = {u: (str(u)) for u in tree.nodes()}
                tree.draw(path="/tmp/sp-{}.svg".format(cname), format='svg',\
                            height=500, width=1000, node_labels=node_labels,\
                            node_colours=node_colours)

            ## Now we have the tree sequence for the clade and can go back through the species list
            ## and pull out stats per species
            metasamps = tree_sequence.samples(population=0)
            for sp, idx in sp_idxs.items():
                ## Skip the metacommunity, which doesn't have a parent :'(
                if sp == '':
                    continue
                sp_obj = species_dict[idx]
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
                                 growth = self.region._hackersonly["population_growth"],\
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

        abunds = np.array([x.stats["abundance"] for x in self.species])
        pis = np.array([x.stats["pi_local"] for x in self.species])
        dxys = np.array([x.stats["dxy"] for x in self.species])
        traits = np.array([x.stats["trait"] for x in self.species])

        dat = pd.DataFrame([], columns=["pis", "dxys", "abunds", "traits"])
        dat["abunds"] = abunds
        dat["pis"] = pis
        dat["dxys"] = dxys
        dat["traits"] = traits

        ss = calculate_sumstats(dat, sgd_bins=self.region._hackersonly["sgd_bins"],\
                                    sgd_dims=self.region._hackersonly["sgd_dimensions"],\
                                    metacommunity_traits=self.region.metacommunity._get_trait_values()) 

        ## If you don't actually want all the intermediate files then we won't make them
        if self.region._log_files:
            megalog = os.path.join(self._hackersonly["outdir"],
                                self.paramsdict["name"] + "-{}-megalog.txt".format(self._lambda()))

            ## concatenate all species results and transpose the data frame so rows are species
            fullstats = pd.concat([sp.stats for sp in self.species] , axis=1).T
            fullstats.to_csv(megalog, index_label=False)

        ## paste on the local parameters and pseudo-parameters
        params = self._get_params_values()

        return params.append(ss.T)


@memoize
def _get_filtering_death_prob(region, victim_trait):
    try:
        val = 1 - (np.exp(-((victim_trait - region.metacommunity._hackersonly["filtering_optimum"]) ** 2)/region.metacommunity.paramsdict["ecological_strength"]))
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

@memoize
def _get_trait_rate_local(region):
    try:
        ext = region.metacommunity.paramsdict["birth_rate"] * region.metacommunity.paramsdict["death_proportion"]
        val = region.metacommunity.paramsdict["trait_rate_meta"]/ (region.metacommunity.paramsdict["birth_rate"] + ext)
    except Exception as inst:
        raise MESSError("Error in geting trait rate for local community")
    return val

#############################
## Model Parameter Info Dicts
#############################
LOCAL_PARAMS = {
    "name" : "Local community name",\
    "K" : "Local carrying capacity",\
    "colrate" : "Colonization rate into local community",\
    "speciation_rate" : "# of new species per forward-time generation",\
    "background_death" : "Baseline death probability in trait-based models",\
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
    loc._hackersonly["mode"] = "landbridge"
    loc.prepopulate()
    assert(len(collections.Counter(loc.local_community)) > 3)
    print(loc.get_abundances())

    loc._hackersonly["mode"] = "volcanic"
    loc.prepopulate()
    loc.step(100000)
    print(len(set(loc.local_community)))
    print(len(loc.local_community))
    print(loc.local_info.shape)
    print("getting stats")
    print(loc)
    print(loc.paramsdict)
