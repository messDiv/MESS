from __future__ import print_function

from scipy.stats import logser
from collections import OrderedDict
from scipy.spatial import distance
from time import time, sleep
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


class LocalCommunity(object):
    """
    Construct a local community.

    :param str name: The name of the LocalCommunity.
    :param int J: The number of individuals in the LocalCommunity.
    :param float m: Migration rate into the LocalCommunity. This is the
        probability per time step that a death is replaced by a migrant
        from the metacommunity.
    :param bool quiet: Print out some info about the local community.
    """
    def __init__(self, name="Loc1", J=1000, m=0.01, quiet=False):
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
                        ("J", J),
                        ("m", m),
                        ("speciation_prob", 0),
        ])

        ## A dictionary for holding prior ranges for values we're interested in
        self._priors = dict([
                        ("J", []),
                        ("m", []),
                        ("speciation_prob", []),
        ])

        ## Dictionary of 'secret' parameters that most people won't want to mess with
        ##  * allow_empty is a switch on whether to fully populate local community
        ##      with one species for 'volcanic' mode or to introduce just one
        ##      individual and populate the rest of J with 'empty' demes.
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
        ##      'continental' style origin.
        self._hackersonly = dict([
                        ("allow_empty", False),
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
        ##  * "split_time" - Time of last divergence event involving this species,
        ##                   for keeping track of internal branch lengths.
        self.local_info = pd.DataFrame([])

        ## The regional pool that this local community belongs to
        ## This is updated by Region._link_local(), so don't set it by hand
        self.region = ""
        ## death_step will be a reference to the function to select the
        ## individual for removal. It is set dependent on which type of
        ## assembly model the region specifies, and is also set by _link_local()
        self.death_step = ""

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

        ## Invasiveness is unimplemented and may not make it into MESS proper.
        ## The invasive species identity
        self.invasive = -1
        self.invasiveness = 0
        ## Track how many invasives differentially survived
        self.survived_invasives = 0
        self.invasion_time = -1
        self.NAME_SEPARATOR = "."

    def _copy(self):
        """
        Create a new clean copy of this LocalCommunity, resampling new values
        for any parameter that were specified with a prior.
        """
        LOGGER.debug("Copying LocalCommunity - {}".format(self.name))
        new = LocalCommunity(self.name)

        new._set_region(self.region)
        new._priors = self._priors

        new.paramsdict = self.paramsdict
        ## Get sample from prior range on params that may have priors
        for param in ["J", "m", "speciation_prob"]:
            ## if _priors is empty then this param is fixed
            if np.any(self._priors[param]):
                loguniform = False
                if param in ["m", "speciation_prob"]:
                    loguniform = True
                self.paramsdict[param] = sample_param_range(new._priors[param], loguniform=loguniform)[0]

        new._hackersonly = self._hackersonly
        LOGGER.debug("Copied LocalCommunity - {}".format(self.paramsdict))
        return new


    def _lambda(self):
        """
        Return fraction of equilibrium obtained by the local community
        """
        percent_equil = 0
        try:
            percent_equil = float(self.founder_flags.count(False))/len(self.founder_flags)
        except ZeroDivisionError:
            pass
        return percent_equil


    def _add_local_info(self, sname, abundances_through_time=0,\
                        ancestor='', ancestral_abundance=[], speciation_completion=0):
        """
        Construct a new local_info record for new_species. The fields are:
        colonization time - in the forward time model. This gets converted to
            divergence time for the backward time model. This is a misnomer,
            held over from the days when speciation wasn't a thing. This should
            more appropriately be called "species age".
        post_colonization_migrants - The count of migrants that have come in
            that are the same species as this one, since colonization
        abundances_through_time - Dictionary containing history of population
            size change through time. Default is 0, which indicates a new colonist.
        ancestor - If the species was introduced by speciation rather than
            colonization, then it'll have an ancestral species.
        ancestral_abundance - A list of fluctuating ancestral abundances at
            the time the species split from its sister. Default is empty list which indicates
            a new colonist from the metacommunity.
        speciation_completion - The point in forward time when this species will become 'good'.
            before this point it is an incipient species and if the simulation is stopped then
             all the individuals will be thrown back into the parent species pool.
            Default is 0 which indicates this is a good species immediately, either a new
            colonizing lineage or a point mutation species.
        split_time = Time of the last speciation event involving a given species.
            During the next speciation event for this species this value will
            will be converted to the length of the internal branch leading to
            this speciation event. When a new species is created colonization_time
            and split_time will be the same, but this value changes when a new
            species is formed _from_ this one.
        """
        if abundances_through_time == 0: abundances_through_time = OrderedDict([(self.current_time,self._hackersonly["mig_clust_size"])])
        self.local_info[sname] = [self.current_time,\
                                        0,\
                                        abundances_through_time,\
                                        ancestor,\
                                        ancestral_abundance,\
                                        speciation_completion,\
                                        self.current_time]


    def _set_region(self, region):
        """
        Connect this LocalCommunity to its region.
        """
        self.region = region
        self.SGD = SGD([], ndims=region._hackersonly["sgd_dimensions"], nbins=region._hackersonly["sgd_bins"])
        self._hackersonly["trait_rate_local"] = self._get_trait_rate_local()
        self._set_death_step()


    def _set_death_step(self):
        """
        Set the death step method based on which community assembly model we're
        running.
        """
        assembly_model = self.region.paramsdict["community_assembly_model"]
        if assembly_model == "neutral":
            self.death_step = self._neutral_death_step
        elif assembly_model == "competition":
            self.death_step = self._competition_death_step
        elif assembly_model == "filtering":
            self.death_step = self._filtering_death_step
            self._filtering_update_death_probs()
        elif assembly_model == "pairwise_competition":
            self.death_step = self._pairwise_competition_death_step
        else:
            raise Exception("unrecognized community assembly model in _set_death_step: {}".format(assembly_model))


    ## Update global distance matrix (and its exp version)
    ## for the pairwise competition model
    def _distance_matrix_init(self):
        loc_inds = [x for x in self.local_community if x != None]
        local_traits = list(map(self.region.get_trait, loc_inds))
        print(len(local_traits))
        local_traits = [[x] for x in local_traits] ## For use in cdist
        es = 1./self.region.metacommunity.paramsdict["ecological_strength"]
        
        nb_ind = len(loc_inds)
        dist_matrix = distance.cdist(local_traits,local_traits,'sqeuclidean')
        self._exp_distance_matrix = np.exp(-(dist_matrix)/es)
        
        
        
    def _distance_matrix_remove(self,idx):
        self._exp_distance_matrix[idx:-1] = self._exp_distance_matrix[idx+1:]
        self._exp_distance_matrix[:,idx:-1] = self._exp_distance_matrix[:,idx+1:]


    def _distance_matrix_add(self):
        loc_inds = [x for x in self.local_community if x != None]
        nb_ind = len(loc_inds)
        local_traits = list(map(self.region.get_trait, loc_inds))
        local_traits = [[x] for x in local_traits] ## For use in cdist
        es = 1./self.region.metacommunity.paramsdict["ecological_strength"]
        new_dist = np.reshape(distance.cdist(local_traits,[local_traits[-1]]),(nb_ind))
        
        self._exp_distance_matrix[-1] = np.exp(-(new_dist/es))
        self._exp_distance_matrix[:,-1] = self._exp_distance_matrix[-1].T
        ## The new individual has been appended at the end of the local community


    ## Update global death probabilities for the filtering model
    def _filtering_update_death_probs(self, perturb=False):
        """
        Update environmental filtering death probabilities. In a model of zero
        speciation these only ever need to get calculated once. With speciation
        they need to get recalculated upon each speciation event. Also,
        we add a subtle perturbation per generation to simulate shifting
        environmental conditions.

        :param bool perturb: Whether to perturb the filtering optimum by +/- a
            small random value.
        """
        fo = self.region.metacommunity._hackersonly["filtering_optimum"]
        es = self.region.metacommunity.paramsdict["ecological_strength"]
        if perturb:
            noise = np.random.normal(1, 0.05)
            fo = fo * noise
#            print("Upstate fo", fo, noise)
            self.region.metacommunity._hackersonly["filtering_optimum"] = fo

        def dprob(trt):
            dp = 1 - (np.exp(-((trt - fo) ** 2)/es))
## Set a minimum death prob?
#            if dp < 0.75:
#                dp = 0.75
            return dp
        self._filt_death_probs = {sp:dprob(trt) for sp, trt in self.region.metacommunity._get_species_traits().items()}

    ## Getting params header and parameter values drops the local
    ## community name (param 0), and adds a bunch of pseudo-parameters
    def _get_params_header(self):
        params_header = list(self.paramsdict.keys())[1:]
        params_header = params_header + ["generation", "_lambda", "migrate_calculated", "extrate_calculated",\
                                            "trait_rate_local", "filtering_optimum"]
        return params_header


    def _get_params_values(self):
        ## Get params and drop name
        params_vals = list(self.paramsdict.values())[1:]
        ## We are reporting generations scaled to WF time
        params_vals = params_vals + [self.current_time * 2 / self.paramsdict["J"],\
                                    self._lambda(),\
                                    self.colonizations/float(self.current_time),\
                                    self.extinctions/float(self.current_time),\
                                    self._hackersonly["trait_rate_local"],\
                                    self.region.metacommunity._hackersonly["filtering_optimum"]]
        params_vals = pd.DataFrame(params_vals, index=self._get_params_header())
        return params_vals


    def _paramschecker(self, param, newvalue, quiet=False):
        """
        Raises exceptions when params are set to values they should not be.
        """
        ## TODO: This should actually check the values and make sure they make sense
        ## TODO: Also check here if you're setting the mode parameter you have to rerun _prepopulate
        try:

            ## Cast params to correct types
            if param in ["J"]:
                tup = tuplecheck(newvalue, dtype=int)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup)[0]
                else:
                    self.paramsdict[param] = tup

            elif param in ["m"]:
                tup = tuplecheck(newvalue, dtype=float)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup, loguniform=True)[0]
                else:
                    self.paramsdict[param] = tup

            elif param == "mode":
                ## Must reup the local community if you change the mode
                self.paramsdict[param] = newvalue

            elif param == "speciation_prob":
                tup = tuplecheck(newvalue, dtype=float)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup, loguniform=True)[0]
                else:
                    self.paramsdict[param] = tup

            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            LOGGER.debug("Bad parameter - {} {}".format(param, newvalue))
            ## Do something intelligent here?
            raise


    def _get_trait_rate_local(self):
        try:
            ext = self.region.metacommunity.paramsdict["speciation_rate"] * self.region.metacommunity.paramsdict["death_proportion"]
            val = self.region.metacommunity.paramsdict["trait_rate_meta"]/ (self.region.metacommunity.paramsdict["speciation_rate"] + ext)
        except Exception as inst:
            raise MESSError("Error in geting trait rate for local community: {}".format(inst))
        return val


    def _write_params(self, outfile=None, full=False):
        """
        Write out the parameters of this LocalCommunity to a file properly
        formatted as input for `MESS -p <params.txt>`.

        :param string outfile: The name of the params file to write to. If not
            specified this will default to `params-<Region.name>.txt`.
        :param bool full: Whether to write out only the parameters of this
            this particular LocalCommunity realization, or to write out the
            parameters including any prior ranges.
        """
        if outfile is None:
            raise MESSError("LocalCommunity._write_params outfile must be specified.")

        with open(outfile, 'a') as paramsfile:
            header = "------- LocalCommunity params: {}".format(self.name)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            for key, val in self.paramsdict.items():
                paramvalue = str(val)

                ## If it's one of the params with a prior, and if the prior is not
                ## empty and if writing out full, then write the prior, and not
                ## the sampled value
                if full:
                    if key in list(self._priors.keys()):
                        if self._priors[key]:
                            paramvalue = "-".join([str(i) for i in self._priors[key]])

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
        """
        A function for occasionally logging a ton of information through time.
        Anything that needs to get recorded through time should happen in here.

        :param bool full: Dump a ton of stuff to the outdir. Normally only really
            used for fancy plotting.
        """
        if full:
            self.lambda_through_time[self.current_time] = self._lambda()
            self._simulate_seqs()
            self.species_through_time[self.current_time] = self.species

        ## Every once in a while test to be sure our community is the same size
        ## as we think it should be.
        if not len(self.local_community) == self.paramsdict["J"]:
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


    def _prepopulate(self, verbose=False):
        LOGGER.debug("prepopulating local_community - {}".format(self))
        if not self.region:
            msg = "Skip populating the local community as it is unlinked to a region."
            LOGGER.error(msg)
            if MESS.__interactive__: print("    {}".format(msg))
            return

        ## Clean up local_community if it already exists
        self.local_community = []

        if self._hackersonly["mode"] == "continental":
            ## prepopulate the island w/ a random sample from the metacommunity
            ## TODO: The underscore here is ignoring trait values
            self.local_community, _ = self.region._get_migrants(self.paramsdict["J"])

        elif self._hackersonly["mode"]  == "volcanic":
            ## If not continental then doing volcanic, so sample just the most abundant
            ## from the metacommunity
            ## TODO: The _ is a standin for trait values, have to do something with them

            try:
                new_species, _ = self.region._get_most_abundant()
            except Exception as inst:
                raise MESSError("Error in _prepopulate - {}".format(inst))

            ## prepopulate volcanic either with all the most abundant species in the metacommunity
            ## or with one sample of this species and a bunch of "emtpy deme space". The empty
            ## demes screw up competition/environmental filtering models
            if self._hackersonly["allow_empty"]:
                self.local_community = [None] * self.paramsdict["J"]
                self.local_community[0] = new_species
            else:
                self.local_community = [new_species] * self.paramsdict["J"]
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
                                                     "speciation_completion",\
                                                     "split_time"])
        self.local_info = self.local_info.fillna(0)
        for sp in self.local_info:
            self.local_info[sp] = [0, 0, OrderedDict(), "", [], 0, 0]

        self.founder_flags = [True] * len(self.local_community)
        if verbose:
            print("    Initializing local community:")
            print("      N species = {}".format(len(set(self.local_community))))
            print("      N individuals = {}".format(len(self.local_community)))
        LOGGER.debug("Done prepopulating - {}".format(self))

        ## Initialize distance matrix if we are in pairwise competition
        if self.region.paramsdict["community_assembly_model"] == "pairwise_competition":
            self._distance_matrix_init()



    def _neutral_death_step(self):
        victim = random.choice(self.local_community)
        self._finalize_death(victim)


    def _competition_death_step(self):
        victim = random.choice(self.local_community)
        if victim == None:
            pass
        else:
            mean_local_trait = self.region.get_trait_mean(self.local_community)
            ## Get local traits for all individuals in the community (remove None first)
            loc_inds = [x for x in self.local_community if x != None]
            local_traits = list(map(self.region.get_trait, loc_inds))
            ## Scale ecological strength for competition to be in the same units
            ## as for filtering
            es = 1./self.region.metacommunity.paramsdict["ecological_strength"]
            ## Apply the competition equation to get fitness per individual
            death_probs = [np.exp(-((x - mean_local_trait) ** 2)/es) for x in local_traits]
            ## Scale all fitness values to proportions
            death_probs = np.array(death_probs)/np.sum(death_probs)
            ## Get the victim conditioning on unequal death probability
            vic_idx = list(np.random.multinomial(1, death_probs)).index(1)
            victim = loc_inds[vic_idx]

        self._finalize_death(victim)


    def _pairwise_competition_death_step(self):
        victim = random.choice(self.local_community)
        if victim == None:
            pass
        else:
            ## Get local traits for all individuals in the community (remove None first)        
            loc_inds = [x for x in self.local_community if x != None]
            nb_ind = len(loc_inds)
          
            death_probs = np.sum(self._exp_distance_matrix,axis=0)-np.diag(self._exp_distance_matrix)

            ## Scale all fitness values to proportions
            death_probs = np.array(death_probs)/np.sum(death_probs)

            ## Get the victim conditioning on unequal death probability
            vic_idx = list(np.random.multinomial(1, death_probs)).index(1)
            victim = loc_inds[vic_idx]
        self._finalize_death(victim)
        ## Update the distance matrix with the new individual
        self._distance_matrix_remove(idx=vic_idx)


    def _filtering_death_step(self):
        victim = random.choice(self.local_community)
        if victim == None:
            pass
        else:
            ## Get local traits for all individuals in the community (remove None first)
            loc_inds = [x for x in self.local_community if x != None]
            ## Fetch pre-calculated death probs for each individual
            death_probs = [self._filt_death_probs[x] for x in loc_inds]
            ## Scale all fitness values to proportions
            death_probs = np.array(death_probs)/np.sum(death_probs)
            ## Get the victim conditioning on unequal death probability
            vic_idx = list(np.random.multinomial(1, death_probs)).index(1)
            victim = loc_inds[vic_idx]

        self._finalize_death(victim)


    def _finalize_death(self, victim):
        ## More old invasiveness code. Should probably just get rid of it eventually.
        ## If no invasive has invaded then just do the normal sampling
        ##if self.invasive == -1:
            ## If no invasive species yet just go on
        ##    pass
        ##else:
            ## If invasiveness is less than the random value remove the invasive individual
            ## else choose a new individual
        ##    if victim == self.invasive and np.random.rand() < self.invasiveness:
        ##        self.survived_invasives += 1
        ##        victim = random.choice(self.local_community)
        try:
            ## Clean up local community list and founder flag list
            idx = self.local_community.index(victim)
            self.local_community.pop(idx)
            self.founder_flags.pop(idx)

            ## If the species of the victim went locally extinct then clean up local_info
            self._test_local_extinction(victim)
        except Exception as inst:
            raise MESSError("Error in _finalize_death(): {}".format(inst))


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
                #LOGGER.debug("\nExtinction victim info \n{}\n{}\nOffspring {}".format(victim, vic_info, offspring))

                ## Speciation process nonsense
                ## Update ancestry and population size change history for any species with this one as
                ## direct ancestor
                ancestor = vic_info["ancestor"]
                anc_size_changes = vic_info["abundances_through_time"]
                ## Get all species in the local community with victim as ancestor
                offspring = self.local_info.loc["ancestor"] == victim
                ## species identities for only the relevant species
                oidx = offspring.index[offspring]
                self.local_info.loc["ancestor", oidx] = ancestor
                # This 'works' but gives the copy of slice pandas warning
                #self.local_info.loc["ancestor"][self.local_info.loc["ancestor"] == victim] = ancestor
                ## I don't think you can vectorize the update

                for o in oidx:
                    LOGGER.debug("offspring {} {}".format(o, self.local_info[o]["abundances_through_time"]))
                    self.local_info[o]["abundances_through_time"].update(anc_size_changes)
                    LOGGER.debug("offspring {} {}".format(o, self.local_info[o]["abundances_through_time"]))

            except Exception as inst:
                LOGGER.debug(self.local_info)
                raise MESSError("Exception during recording extinction - {}".format(inst))
            ## If the invasive prematurely goes extinct just pick a new one
            if victim == self.invasive:
                LOGGER.info("invasive went extinct")
                self.invasive = -1
        return ancestor


    def _migrate_step(self):
        """
        Allow multiple colonizations. In this case we return the sampled species
        as well as a bool reporting whether or not this is the first colonization of
        this species into the local community so that the coltime can be recorded.
        multiple colonizations of a species do not update coltime, but we record them
        for migration rate calculation.
        """

        new_species, _ = self.region._get_migrant()
        if new_species in self.local_community:
            ## This is a post-colonization migrant so record the event and tell downstream
            ## not to update the colonization time.
            self.local_info[new_species]["post_colonization_migrants"] += 1
        else:
            ## This is a new migrant so init local_info for it
            self._add_local_info(sname = new_species)

        return new_species


    def step(self, nsteps=1):
        t0 = time()
        """
        Run one or more generations of birth/death/colonization timesteps. A
        generation is J/2 timesteps (convert from Moran to WF generations).

        :param int nsteps: The number of generations to simulate.
        """
        ## Convert time in generations to timesteps (WF -> Moran)
        for step in range(int(nsteps * self.paramsdict["J"]/2.)):
            chx = ''
            ## Check probability of an immigration event
            if np.random.random_sample() < self.paramsdict["m"]:
                ## If clustered migration remove the necessary number of additional individuals
                for _ in range(self._hackersonly["mig_clust_size"]):
                    self.death_step()

                ## Grab the new colonizing species
                new_species = self._migrate_step()
                chx = new_species

                ## The invasion code "works" in that it worked last time I tried it, but it's
                ## not doing anything right now except slowing down the process. I don't want to dl
                ## it in case we want to resurrect, so it's just commented out for now. 5/2019 iao.
                ## Only set the invasive species once at the time of next migration post invasion time
                ## If invasion time is < 0 this means "Don't do invasive"
                #if not self.invasion_time < 0:
                #    if self.invasive == -1 and self.current_time >= self.invasion_time:
                #        self.invasive = new_species
                #        LOGGER.info("setting invasive species {} at time {}".format(self.invasive, self.current_time))
                #        self.invasion_time = self.current_time

                ## Add the colonizer to the local community, record the colonization time
                self.local_community.extend([new_species] * self._hackersonly["mig_clust_size"])
                self.founder_flags.extend([False] * self._hackersonly["mig_clust_size"])
                self.colonizations += 1
            else:
                try:
                    self.death_step()
                    ## Sample all available from local community (community grows slow in volcanic model)
                    ## This is the fastest way to sample from a list. >4x faster than np.random.choice
                    ## Don't allow empty space to reproduce
                    chx = random.choice([x for x in self.local_community if x != None])
                    self.local_community.append(chx)
                    idx = self.local_community.index(chx)
                    self.founder_flags.append(self.founder_flags[idx])
                except Exception as inst:
                    LOGGER.error("Exception in step() - {}".format(inst))
                    raise inst

            ## WARNING : Pairwise competition assumes that "mig_clust_size" is one !

            if self.region.paramsdict["community_assembly_model"] == "pairwise_competition":
                if self._hackersonly["mig_clust_size"] != 1:
                    raise MESSError("Pairwise competition only handles migration cluster of size 1")
                self._distance_matrix_add()

            ##############################################
            ## Speciation process
            ##############################################
            if self.region.paramsdict["speciation_model"] != "none" and\
               np.random.random_sample() < self.paramsdict["speciation_prob"] and\
               chx != None:

               self._speciate(chx)


            ## update current time
            self.current_time += 1

        ## Perturb environmental optimum for filtering
        self._filtering_update_death_probs(perturb=True)
        t1=time()
        #print("timestep: ",t1-t0)


    def get_abundances(self, octaves=False, raw_abunds=False):
        """
        Get the SAD of the local community.

        :param bool octaves: Return the SAD binned into size-class octaves.
        :param bool raw_abunds: Return the actual list of abundances per
            species, without binning into SAD.

        :return: If `raw_abunds` then returns a list of abundances per species,
            otherwise returns an OrderedDict with keys as abundance classes
            and values as counts of species per class.
        """
        return SAD(self.local_community, octaves=octaves, raw_abunds=raw_abunds)


    def _speciate(self, chx):
        """
        Occassionally initiate the speciation process. In all modes, one
        one individual is randomly selected to undergo speciation.
        Speciation does not change the founder_flag state.

        Currently there are 3 modes implemented::

            point_mutation - The randomly selected individual becomes a new
            species of its own, of abundance 1.

            random_fission -  The species of the randomly selected individual
            undergoes random fission. In this mode the abundance of the
            new species is determined by randomly splitting off a chunk
            of the individuals from the parent species. All fission
            sizes are equally likely.

            protracted - watdo
        """
        LOGGER.debug("Initiate speciation process - {}".format(chx))

        idx = self.local_community.index(chx)

        ## Construct the new species name.
        ## We want the new name to be globally unique but we don't
        ## want to waste a bunch of time keeping track of it in the
        ## region, so we can do something like this:
        sname = chx + "{}{}-{}".format(self.NAME_SEPARATOR, self.name, self.current_time)

        ## Fetch the abundance history of the parent species
        parent_abunds = self.local_info[chx]["abundances_through_time"]

        ## Inform the regional pool that we have a new species
        ## Identify parent's trait value
        parent_trait = self.region.get_trait(chx)

        ## Trait evolution. Offspring trait is normally distributed
        ## with mean of parent value, and stdv equal to stdv of BM
        ## process in metacommunity times average lineage lifetime
        trt = np.random.normal(parent_trait, self._hackersonly["trait_rate_local"], 1)[0]

        ## Get internal branch length leading to this speciation event
        ## and update split_time for the parental
        split_time = self.local_info[chx]["split_time"]
        self.local_info[chx]["split_time"] = self.current_time

        ## split_time is an absolute time, so convert it to relative
        ## time.
        split_time = self.current_time - split_time
        ## scale to branch length in Mya
        ## TODO: This assumes Mya is the right timescale and also
        ##       implicitly assumes generation time is 1 year.
        split_time = split_time/1e6

        self.region._record_local_speciation(sname=sname,\
                                            trait_value=trt,\
                                            ancestor=chx,\
                                            branch_length=split_time)

        ## If filtering then update the death probabilities to record
        ## death probability of the new species
        if self.region.paramsdict["community_assembly_model"] == "filtering":
            self._filtering_update_death_probs()

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
            ## Speciation flips the founder_flag for the new species
            self.founder_flags[idx] = False

        elif self.region.paramsdict["speciation_model"] == "random_fission":
            ## TODO: This doesn't handle the founder_flag housekeeping AT ALL!

            ## Remove all individuals of the target species from the local community.
            ## We'll have to add some back as the original species, but doing it
            ## this way simplifies the housekeeping.
            self.local_community = [sp for sp in self.local_community if sp != chx]
            ## Get abundance of the target species so we can perform
            ## the random fission.
            sp_abund = self.paramsdict["J"] - len(self.local_community)

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


    ##TODO: Unused and not updated to current MESS structure
    def _bottleneck(self, strength=1):
        """
        How strong is the bottleneck? Strength should be interpreted as percent
        of local community to retain
        """
        reduction = int(round(self.paramsdict["J"] * strength))
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
    def _get_clades(self):
        ## We need to get all the groups of species that descended from a common
        ## local anscestor so we can run the backwards time model for all of them
        ## combined in a coherent fashion.

        ## Sort them by colonization time
        self.local_info.sort_values("colonization_times", axis=1, inplace=True)

        ## Clades will all be descended from one unique ancestor in the metacommunity
        ## so we will just group them by this identifier. Species names have this
        ## format: t0.<local_name>-<split_time>, so the zeroth element of the split
        ## on the NAME_SEPARATOR char ('.' at this moment)  will always be the
        ## metacommunity species id.
        clades = {x:[] for x in set([y.split("{}".format(self.NAME_SEPARATOR))[0] for y in self.local_info])}

        for idx in self.local_info:
            clades[idx.split("{}".format(self.NAME_SEPARATOR))[0]].append(idx)

        LOGGER.debug("Clades - {}".format(clades))

        return clades


    def _simulate_seqs(self):
        """
        Simulate genetic variation for each species in the local community.
        """
        self.species = []
        local_info_bak = self.local_info.copy(deep=True)
        try:
            ## Hax. Remove the empty deme from local info. This _might_ break the fancy plots.
            self.local_community = [x for x in self.local_community if x != None]
            self.local_info = self.local_info[self.local_info.columns.dropna()]
        except:
            ## df.drop will raise if it doesn't find a matching label to drop, in which case we're done.
            pass

        self.local_info.loc["meta_abund"] = [self.region._get_abundance(x) for x in self.local_info]
        self.local_info.loc["local_abund"] = [self.local_community.count(x) for x in self.local_info]
        self.local_info.loc["colonization_times"] = self.current_time - self.local_info.loc["colonization_times"]
        for cname, species_list in self._get_clades().items():
            ## Just get a slice of the local_info df that represents the species of interest
            dat = self.local_info[species_list]
            ## Dictionary mapping species names to 0-based index (where 0 is metacommunity)
            sp_idxs = OrderedDict({x:ix for ix, x in zip(list(range(1, len(dat.columns)+1)), dat.columns[::-1])})
            sp_idxs[''] = 0

            pop_cfgs = []
            split_events = []
            meta_abund = self.region._get_abundance(cname)
            try:
                pop_meta = msprime.PopulationConfiguration(sample_size = 10, initial_size = meta_abund)
            except ValueError:
                raise
            pop_cfgs.append(pop_meta)

            # make the tree for this clade
            all_sp = self.region.metacommunity.community["ids"]
            drop = set(all_sp).difference(set(species_list))
            tre = self.region.metacommunity.metacommunity_tree.drop_tips(drop)
            tre = tre.write(tree_format=TREE_FORMAT)

            species_dict = {}
            for sp, idx in sp_idxs.items():
                if not sp:
                    continue
                try:
                    migration_rate = dat[sp]["post_colonization_migrants"]/float(dat[sp]['colonization_times'])
                except ZeroDivisionError as inst:
                    ## This should only happen when coltime is 0, which should be never
                    LOGGER.error("Got bad divergence time - {}".format(dat[sp]['colonization_times']))
                    migration_rate = 0

                sp_obj = species(name = sp,
                         species_params = self.region.get_species_params(),
                         tree = tre,
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
            if False:
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

        ## Return local_info to the original state in case you want to keep running
        ## forward simulations
        self.local_info = local_info_bak


    def get_community_data(self, trees=False):
        """
        Gather the community data and format it in such a way as to prepare it
        for calling MESS.stats.calculate_sumstats(). This is a way of getting
        simulated data that is in the exact format empirical data is required
        to be in. Useful for debugging and experimentation.

        :return: A pandas.DataFrame with 5 columns: "pi", "dxy", "abundance",
            "trait", and "tree", and one row per species.
        """
        abunds = np.array([x.stats["abundance"] for x in self.species])
        pis = np.array([x.stats["pi_local"] for x in self.species])
        dxys = np.array([x.stats["dxy"] for x in self.species])
        traits = np.array([x.stats["trait"] for x in self.species])
        trees = np.array([x.stats["tree"] for x in self.species])

        dat = pd.DataFrame([], columns=["pi", "dxy", "abundance", "trait"])
        dat["abundance"] = abunds
        dat["pi"] = pis
        dat["dxy"] = dxys
        dat["trait"] = traits
        dat["tree"] = trees

        return dat


    def get_stats(self):
        """
        Simulate genetic variation per species in the local community, then
        aggregate abundance, pi, dxy, and trait data for all species
        and calculate summary statistics.

        :return: A pandas.DataFrame including all MESS model parameters and
            all summary statistics.
        """
        LOGGER.debug("Entering get_stats()")
        self._simulate_seqs()
        LOGGER.debug("First 5 species - \n{}".format(self.species[:5]))

        dat = self.get_community_data()

        ## We have to do a bunch of preprocessing for the phylogenetic data

        ## Find all species that speciated in the local community and then
        ## went extinct. We need to prune these out of the metacommunity tree.
        all_local_sp = [x for x in self.region.metacommunity.community["ids"] if self.NAME_SEPARATOR in x]
        extinct_local_sp = set(all_local_sp).difference(set(self.local_community))
        self.extant_meta_toytree = self.region.metacommunity.metacommunity_tree.drop_tips(extinct_local_sp)

        ## We need to update the branch lengths at the tips for extant taxa
        for node in self.extant_meta_toytree.treenode.traverse():
            ## Don't update branch lengths for species that didn't speciate
            ## in the local community
            if node.is_leaf() and node.dist == 0:
                try:
                    node.dist = (self.current_time - self.local_info[node.name]["split_time"])/1e6

                except KeyError:
                    ## If not in local_info, then this is a tip from the
                    ## metacommunity, so the branch length is already fine.
                    pass

        ## The local toytree is only for debugging purposes, as it will be
        ## created independently by the calculate_sumstats method.
        all_non_local = set(self.extant_meta_toytree.get_tip_labels()).difference(set(self.local_community))
        self.extant_local_toytree = self.extant_meta_toytree.drop_tips(all_non_local)

        meta_tree_newick = self.extant_meta_toytree.write(tree_format=TREE_FORMAT)

        ss = calculate_sumstats(dat, sgd_bins=self.region._hackersonly["sgd_bins"],\
                                    sgd_dims=self.region._hackersonly["sgd_dimensions"],\
                                    metacommunity_traits=self.region.metacommunity._get_trait_values(),
                                    metacommunity_tree=meta_tree_newick,
                                    normalize_hills=self.region._hackersonly["normalize_hills"],
                                    verbose=False)

        ## If you don't actually want all the intermediate files then we won't make them
        if self.region._log_files:
            self._log_files()

        ## paste on the local parameters and pseudo-parameters
        params = self._get_params_values()
        return pd.concat([params, ss.T], axis=0)


    def _log_files(self):
        """
        Write a gigantic file of all the simulation to megalog.txt. This is
        mostly used for fancy plotting. Also dump the simulated data to the
        outdir in a form appropriate for MESS.util.import_empirical().
        """
        megalog = os.path.join(self._hackersonly["outdir"],
                            self.paramsdict["name"] + "-{}-megalog.txt".format(self._lambda()))

        ## concatenate all species results and transpose the data frame so rows are species
        fullstats = pd.concat([sp.stats for sp in self.species] , axis=1).T
        fullstats.to_csv(megalog, index_label=False)

        outdir = os.path.join(self._hackersonly["outdir"],
                            self.paramsdict["name"] + "-{}-data".format(self._lambda()))
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        meta_tree_newick = self.extant_meta_toytree.write(tree_format=TREE_FORMAT)
        with open(os.path.join(outdir, "{}-metacommunity.tre".format(self.paramsdict["name"])), 'w') as treefile:
            treefile.write(meta_tree_newick)

        dat = self.get_community_data()
        dat = dat.drop("pi", axis=1)
        dat = dat.drop("dxy", axis=1)
        dat = dat.drop("tree", axis=1)
        dat = pd.concat([pd.DataFrame([x.name for x in self.species], columns=["name"]), dat], axis=1)
        datfile = os.path.join(outdir, "{}.csv".format(self.paramsdict["name"]))
        dat.to_csv(datfile, sep=",", index=False)

        fastadir = os.path.join(outdir, self.paramsdict["name"] + "-fastas")
        if not os.path.exists(fastadir):
            os.mkdir(fastadir)

        # Generate synthetic fastq data per species
        for sp in self.species:
            haps = list(sp.tree_sequence.haplotypes())[:sp.paramsdict["sample_size_local"]]

            length = self.region.paramsdict["sequence_length"]
            nsnps = len(haps[0])

            seq = np.random.choice(["C", "A", "T", "G"], length-nsnps)
            ref_dat = "".join(np.random.choice(["C", "A", "T", "G"], nsnps))
            alt_dat = "".join(np.random.choice(["C", "A", "T", "G"], nsnps))

            idxs = np.random.choice(list(range(0, length-nsnps)), nsnps, replace=False)

            fastas = []
            for sample_idx, hap in enumerate(haps):
                fastas.append(">Sample-{}".format(sample_idx))

                snps = np.array([int(x) for x in list(hap)])

                ref_idxs = idxs[snps==0]
                alt_idxs = idxs[snps==1]
                tmpseq = np.insert(seq, ref_idxs, ref_dat)
                tmpseq = np.insert(tmpseq, alt_idxs, alt_dat)
                fastas.append("".join(tmpseq))

                with open(os.path.join(fastadir, "Species-{}.fasta".format(sp.name)), 'w') as fasta_file:
                    fasta_file.write("\n".join(fastas))


#############################
## Model Parameter Info Dicts
#############################
LOCAL_PARAMS = {
    "name" : "Local community name",\
    "J" : "Number of individuals in the local community",\
    "m" : "Migration rate into local community",\
    "speciation_prob" : "Probability of speciation per timestep in local community",\
}


#############################
## Error messages
#############################
BAD_MODE_PARAMETER = """
    Unrecognized local community mode. Options are 'continental' or'volcanic'.
    You put {}.
"""


if __name__ == "__main__":
    data = MESS.Region("tmp")
    loc = LocalCommunity("wat", J=5000)
    data._link_local(loc)
    print(loc)
    print(loc.local_info)
    ## Allow for either having or not having empty demes
    assert(len(collections.Counter(loc.local_community)) <= 2)
    loc._hackersonly["mode"] = "continental"
    loc._prepopulate()
    assert(len(collections.Counter(loc.local_community)) > 3)
    print(loc.get_abundances())

    loc._hackersonly["mode"] = "volcanic"
    loc._prepopulate()
    loc.step(100000)
    print(len(set(loc.local_community)))
    print(len(loc.local_community))
    print(loc.local_info.shape)
    print("getting stats")
    print(loc)
    print(loc.paramsdict)
