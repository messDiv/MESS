
import dendropy
import collections
import pandas as pd
import numpy as np
import itertools
import random
import sys
import os
import warnings
import MESS

from collections import OrderedDict
from MESS.util import tuplecheck, sample_param_range, MESSError, set_params
from rpy2 import robjects
from scipy.stats import lognorm

import logging
LOGGER = logging.getLogger(__name__)

## id is a variable length string so we set the dtype as "object"
## to allow for reference pointing to string objects
METACOMMUNITY_DTYPE = np.dtype([('ids', object),
                                ('immigration_probabilities', 'f8'),
                                ('abundances', 'i8'),
                                ('trait_values', 'f8')])

class Metacommunity(object):
    """
    The metacommunity from which individuals are sampled for colonization to the
    local community.

    :param str meta_type: Specify the distribution of abundance among species
        in the metacommunity. Can be one of: `logser`, `lognorm`, `uniform`,
        or the file name from which to read metacommunity abundances.
    :param bool quiet: Whether to print info about metacommunity construction.
    """

    def __init__(self, meta_type="logser", quiet=False):
        self.quiet = quiet

        ## If you add a parameter to this dictionary you need
        ## to also add a short description to the LOCAL_PARAMS dict
        ## at the end of this file
        ##
        ## Also be sure to add it to _paramschecker so the type gets set correctly
        self.paramsdict = OrderedDict([
                        ("S_m", 100),
                        ("J_m", 750000),
                        ("speciation_rate", 2),
                        ("death_proportion", 0.7),
                        ("trait_rate_meta", 2),
                        ("ecological_strength", 1),
                        ("intrasp_competition", 1),
                        ("intersp_competition", 1),
                        ("mutualism_proportion", 0.5),
        ])

        ## elite hackers only internal dictionary, normally you shouldn't mess with this
        ##  * metacommunity_type: Options: uniform/logser/<filename>
        ##  * lognorm_shape: Shape parameter of the lognormal distribution, if you
        ##      choose that option, otherwise it does nothing.
        ##  * filtering_optimum: optimum trait value, only used during environmental
        ##      filtering model. It isn't even a real parameter at this point because
        ##      you can't set it, it's constructed during the simulation.
        self._hackersonly= OrderedDict([
                        ("metacommunity_type", meta_type),
                        ("lognorm_shape", 1.98),
                        ("filtering_optimum", 1),
                        ("inter_term_rate", 0.01), #Arbitrary : to be precised
        ])

        ## A dictionary for holding prior ranges for values we're interested in
        self._priors = dict([
                        ("S_m", []),
                        ("J_m", []),
                        ("speciation_rate", []),
                        ("death_proportion", []),
                        ("trait_rate_meta", []),
                        ("ecological_strength", []),
        ])

        ## The Newick formatted tree for the metacommunity
        self.metacommunity_tree = ""

        ## A structured numpy array for holding tip labels, abundances, colonization
        ## probabilities and trait values
        self.community = []

        ## A dictionary for fast trait value lookups
        self.trait_dict = {}

        ## Adictionary for extended interaction terms (when speciation)
        self.added_species = {}

        self._set_metacommunity()
        LOGGER.debug("Metacommunity paramsdict - {}".format(self.paramsdict))



    def __str__(self):
        return "<Metacommunity: {} Richness {}>".format(self._hackersonly["metacommunity_type"],\
                                                        self.paramsdict["S_m"])

    def _resample_priors(self):
        """
        Draw new parameter valuse if priors are specified for any parameters.
        Parameters are sampled uniform, except ecological strength, which is
        sampled loguniform.
        """
        for k,v in list(self._priors.items()):
            if np.array(v).any():
                loguniform = False
                if k in ["ecological_strength"]:
                    loguniform = True
                self.paramsdict[k] = sample_param_range(v, loguniform=loguniform)[0]


    def _simulate_metacommunity(self, J, S_m, speciation_rate, death_proportion, trait_rate_meta):
        """
        Simulate the tree, evolve the traits on it and paste on abundances.
        Calls out to R code.
        """

        make_meta = """## function to make the mainland meta community with a phylo, traits, and abundances
        ## required packages:
        ##    ape,
        ##    TreeSim,
        ##    nleqslv
        ## arguments:
        #' @param J the number of individuals in the meta community
        #' @param S_m the number of species in the meta community
        #' @param speciation_rate the speciation rate
        #' @param death_proportion the proportional extinction rate
        #' @param trait_rate_meta the rate of brownian motion
        makeMeta <- function(Jm, S, lambda, deathFrac, sigma2) {
          ## the tree
          ## suppress an annoying message about geiger overwritting an ape method
          tre <- suppressMessages(TreeSim::sim.bd.taxa(S, numbsim = 1, lambda = lambda, mu = lambda * deathFrac,
                                      complete = FALSE)[[1]])

          ## the traits
          trt <- ape::rTraitCont(tre, sigma = sqrt(sigma2))
          trt <- data.frame(name = names(trt), value = as.numeric(trt))

          ## parameters for the log-series
          nBar <- Jm / S
          p <- 1 - 1/nBar
          b <- -log(p)

          ## the abundances
          sad_logseries = function(S, N) {
              expit = function(x) {
                  return(1 / (1 + exp(-x)))
              }
              
              f = function(ax, S, N) {
                  ax[2] <- expit(ax[2])
                  
                  return(c(-ax[1] * log(1 - ax[2]) - S, prod(ax) / (1 - ax[2]) - N))
              }
              
              Pn = function(n, ax) {
                  o <- ax[1] / n * ax[2]^n
                  S <- -ax[1] * log(1 - ax[2])
                  
                  return(o / S)
              }
              
              axSol = nleqslv::nleqslv(c(20, 2), f, S = S, N = N, control = list(maxit = 1000))$x
              axSol[2] = expit(axSol[2])
              
              PP = cumsum(Pn(1:(N - S + 1), axSol))
              p2match = seq(1, 1/S, length.out = S) - 1 / (2 * S)
              nn = rep(0, S)
              
              for(i in 1:S) {
                  nn[i] = which.min(abs(p2match[i] - PP))
              }
    
              return(nn)
          }
          #abund <- sads::rls(length(trt), length(trt), 0.01)
          #abund <- meteR::sad(meteR::meteESF(S0 = S, N0 = Jm))$r(S)
          abund <- sad_logseries(S, Jm)

          ## return it all in a list
          tre <- ape::write.tree(tre)
          return(list(phylo = tre, traits = trt, abundance = abund))
        }"""

        rmake_meta = robjects.r(make_meta)
        res = rmake_meta(J,
                        S_m,
                        float(speciation_rate),
                        float(death_proportion),
                        float(trait_rate_meta))
        tree = res[0][0]
        traits = pd.DataFrame(np.array(res[1]).T,
                                columns=["name", "value"])
        abunds = np.array(res[2])

        return tree, abunds, traits


    def _paramschecker(self, param, newvalue, quiet=True):
        """
        Raises exceptions when params are set to values they should not be.
        """
        ## TODO: This should actually check the values and make sure they make sense
        try:
            if (not quiet) and MESS.__interactive__:
                print("  Updating Metacommunity parameters requires running _set_metacommunity()"\
                        + " to apply the changes.")

            ## Cast params to correct types 
            if param in ["S_m", "J_m", "speciation_rate", "death_proportion", "trait_rate_meta",
                            "ecological_strength","mutualism_proportion"]:
                dtype = float 
                if param in ["S_m", "J_m"]:
                    dtype = int
                tup = tuplecheck(newvalue, dtype=dtype)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    loguniform = False
                    if param in ["ecological_strength"]:
                        loguniform = True
                    self.paramsdict[param] = sample_param_range(tup, loguniform=loguniform)[0]
                else:
                    self.paramsdict[param] = tup
                LOGGER.debug("{} {}".format(param, tup))

            elif param == "intrasp_competition":
                if newvalue=="*":
                    alpha = np.random.gamma(0.25,4)
                    beta = np.random.gamma(0.25,4)
                    newvalue = (alpha,beta)
                # Aribitrary priors here ! Should be tuneable from parameter file !
                tup = tuplecheck(newvalue, dtype=float)
                self.paramsdict[param] = tup
                #First value : alpha | second value : beta | if only one : fixed

            elif param == "intersp_competition":
                if newvalue=="*":
                    alpha = np.random.gamma(0.25,4)
                    beta = np.random.gamma(0.25,4)
                    newvalue = (alpha,beta)
                # Aribitrary priors here ! Should be tuneable from parameter file !
                tup = tuplecheck(newvalue, dtype=float)
                self.paramsdict[param] = tup
                #First value : alpha | second value : beta | if only one : fixed


        except Exception as inst:
            ## Do something intelligent here?
            raise MESSError("Error {}\n    Bad parameter {} - Bad value {}".format(inst, param, newvalue))


    def _get_params_header(self):
        return list(self.paramsdict.keys())


    def _get_params_values(self):
        return list(self.paramsdict.values())


    def _get_trait_values(self):
        return list(self.community["trait_values"])


    def _get_species_traits(self):
        return self.trait_dict


    def _get_interaction_term(self,species1,species2):
        
        return self.interaction_matrix[self.species_dict[species1]][self.species_dict[species2]]
    #    try:
    #         return self.interaction_matrix[int(float(species1))-1][int(float(species2))-1]
    # #Species from 1 to 100 but matrix from 0 to 99
    #     except ValueError:
    #         if species1 in self.added_species:
    #             species1 = self.added_species[species1]
    #         if species2 in self.added_species:
    #             species2 = self.added_species[species2]
    #         return self.interaction_matrix[int(float(species1))-1][int(float(species2))-1]

    def _create_interaction(self,new, parent):
        ## Extend matrix
        self.interaction_matrix = np.hstack([
            np.vstack([
            self.interaction_matrix, np.zeros((len(self.interaction_matrix)))
            ]),
            np.zeros((1,len(self.interaction_matrix)+1)).T
            ])
        ## Add interspecific interaction
        if type(self.paramsdict["intersp_competition"])==type((0,0)):
            # Draw according to gamma distribution 
            # self.interaction_matrix[-1] = [np.random.normal(
            #     self._get_interaction_term(parent, sp),
            #     self._hackersonly["inter_term_rate"], 1)[0]
            # for sp in range(1,len(self.interaction_matrix)+1)]

            # self.interaction_matrix[:,-1] = [np.random.normal(
            #     self._get_interaction_term(sp, parent),
            #     self._hackersonly["inter_term_rate"], 1)[0]
            # for sp in range(1,len(self.interaction_matrix)+1)]
            rds = np.random.gamma(shape = self.paramsdict["intersp_competition"][0]
                            , scale = self.paramsdict["intersp_competition"][1]
                            , size = (2,len(self.interaction_matrix)))
            self.interaction_matrix[-1] = rds[0]
            self.interaction_matrix[:,-1] = rds[1]

        else:
            ## Expecting just a single value 
            self.interaction_matrix[-1] = np.array([self.paramsdict["intersp_competition"]] * len(self.interaction_matrix))
            self.interaction_matrix[:,-1] = np.array([self.paramsdict["intersp_competition"]] * len(self.interaction_matrix))

        ## Add intraspecific interaction
        ## Write over precedent terms for the diagonal of the matrix
        if type(self.paramsdict["intrasp_competition"])==type((0,0)):
            # ## Draw according to gamma distribution 
            # rd = np.random.normal(
            #     self._get_interaction_term(parent, parent),
            #     self._hackersonly["inter_term_rate"], 3)
            np.random.gamma(shape = self.paramsdict["intrasp_competition"][0]
                                , scale = self.paramsdict["intrasp_competition"][1]
                                , size = 1)
            # self.interaction_matrix[-1][-1] = rd[0]
            # try :
            #     self.interaction_matrix[-1][int(float(parent))-1] = rd[1]
            #     self.interaction_matrix[int(float(parent))-1][-1] = rd[2]
            # except ValueError:
            #     parent = self.added_species[parent]
            #     self.interaction_matrix[-1][int(float(parent))-1] = rd[1]
            #     self.interaction_matrix[int(float(parent))-1][-1] = rd[2]
            # Interaction with parent species derives from intraspecific interaction
        else:
            ## Expecting just a single value 
            self.interaction_matrix[-1][-1] = self.paramsdict["intrasp_competition"]
      

        # ## Nature of the interaction
        # if type(self.paramsdict["intrasp_competition"])==type((0,0)) and type(self.paramsdict["intersp_competition"])==type((0,0)):
        #     factors = np.random.binomial(n = 1,
        #         p = self.paramsdict["mutualism_proportion"],
        #         size =(2,len(self.interaction_matrix)))
        #     factors[factors == 0] = -1 
        #     self.interaction_matrix[-1] = self.interaction_matrix[-1] * factors[0]
        #     self.interaction_matrix[:,-1] = self.interaction_matrix[:,-1] * factors[1].T
        #     # A positive value is a positive interaction,
        #     # A negative value is a negative interaction
        ## NO redraw from type if evolution from parent
        self.species_dict[new] = len(self.interaction_matrix)-1


    def _write_params(self, outfile=None, full=False):
        """
        Write out the parameters of this Metacommunity to a file properly
        formatted as input for `MESS -p <params.txt>`.

        :param string outfile: The name of the params file to write to. If not 
            specified this will default to `params-<Region.name>.txt`.
        :param bool full: Whether to write out only the parameters of this 
            this particular Metacommunity realization, or to write out the 
            parameters including any prior ranges.
        """
        if outfile is None:
            raise MESSError("Metacommunity._write_params outfile must be specified.")

        with open(outfile, 'a') as paramsfile:
            header = "------- Metacommunity params: "
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            for key, val in list(self.paramsdict.items()):
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


    def _set_metacommunity(self, random=False, resample=False):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there

        :param bool random: Whether to generate random trait values in the 
            range [0-1]
        :param bool resample: Whether to resample from any specified priors.
        """
        meta_type = self._hackersonly["metacommunity_type"]
        LOGGER.debug("Enter set_metacommunity - {}".format(meta_type))

        if resample:
            self._resample_priors()

        ## Accumulators for bringing in all the values. These will
        ## eventually all get shoved into self.community
        abundances = np.array([])
        ids = np.array([])
        trait_values = np.array([])

        ## Two distributions are being left in here as hidden options, the ids and trait values
        ## will get populated below
        ##
        ## These are for testing purposes and should normally be ignored
        if meta_type == "lognorm":
            abundances = lognorm.rvs(self._hackersonly["lognorm_shape"],\
                                        loc=1,
                                        size=self.paramsdict["S_m"])

        elif meta_type == "uniform":
            abundances = np.array([self.paramsdict["J_m"] / self.paramsdict["S_m"]]\
                               * self.paramsdict["S_m"])

        ## Get Abundances by simulating a tree, evolving traits on it, and sprinkling abundances
        ## This is the primary means of setting the metacommunity that is driven by Andy's
        ## R code.
        elif meta_type == "logser":
            tree, abunds, traits = self._simulate_metacommunity(self.paramsdict["J_m"],\
                                                                self.paramsdict["S_m"],\
                                                                self.paramsdict["speciation_rate"],\
                                                                self.paramsdict["death_proportion"],\
                                                                self.paramsdict["trait_rate_meta"])
            #handle = Phylo.read(StringIO(tree), "newick")
            handle = dendropy.Tree.get(data=tree, schema="newick")
            self.metacommunity_tree = handle

            abundances = abunds
            ids = traits["name"].values.astype(str)
            trait_values = traits["value"].values

            self._hackersonly["filtering_optimum"] = np.random.normal(loc=np.mean(trait_values), scale=np.std(trait_values), size=1)[0]

        ## Attempt to read tree/ids/abunds/traits from a file. If it fails, fall back to just
        ## try reading the old list of abundances format.
        ##
        ## TODO: This input file format doesn't include a tree, should we just delete this option?
        ## Both of these input file formats should be considered deprecated for the most part.
        elif os.path.isfile(meta_type):
            try:
                with open(meta_type, 'r') as infile:
                    lines = infile.readlines()
                    self.metcommunity_tree_height = float(lines[0].split()[0])
                    self.paramsdict["trait_rate_meta"] = float(lines[1].split()[0])

                    for i in range(2,len(lines)):
                        info = lines[i].split()
                        ## Doing these somewhat out of order because i'm relying
                        ## on the info[1] failing and raising for the old type
                        ## file that is only a list of abundances. If you do the
                        ## ids[0] first it'll succeed and then fuck up the ids col
                        ## downstream.
                        np.append(trait_values, float(info[1]))
                        np.append(abundances, int(info[2]))
                        np.append(ids, info[0])
            except IndexError as inst:
                ## Could be an old style file just containing abundances, one per row
                try:
                    with open(meta_type, 'r') as infile:
                        lines = infile.readlines()
                        abundances = [int(line.split()[0]) for line in lines]
                except Exception as inst:
                    raise MESSError("  Malformed metacommunity specification file - {}\n    {}".format(meta_type, inst))

            ## If reading from a file then the number of species will not correspond
            ## with the value already in the paramsdict, so we need to update the S_m count
            ## and reup the community ndarray
            LOGGER.debug("Read S_m from file - {}".format(len(abundances)))
            self.paramsdict["S_m"] = len(abundances)
        else:
            raise MESSError("  Unrecognized metacommunity input - {}".format(meta_type))

        ## This next set of conditionals is responsible for filling in trait values and ids
        ## for the abundance only file and the 2 distributions specified at the beginning.
        ## TODO: optionally set random trait values?

        if random or not trait_values.size:
            LOGGER.debug("Using random trait values")
            trait_values = np.random.rand(self.paramsdict["S_m"])

        ## If ids haven't been assigned yet, do that here
        if not ids.size:
            ids = np.array(["t"+str(x) for x in range(0, self.paramsdict["S_m"])])

        ## TODO: You can use msprime to simulate a random coalescent tree and get a newick from it.
        ## Is this worth doing?
        if not self.metacommunity_tree:
            pass

        self.community = np.zeros([self.paramsdict["S_m"]], dtype=METACOMMUNITY_DTYPE)
        ## Populate the metacommunity ndarray
        ## FIXME: Sometimes the logser R call returns a number of species not == to the #
        ##        requested, so this will raise. The better thing to do here would be to
        ##        fix the R code.
        try:
            self.community["abundances"] = np.array(abundances)
            self.community["ids"] = ids
            self.community['trait_values'] = np.array(trait_values)

            self.trait_dict = {x["ids"]:x["trait_values"] for x in self.community}
            ## Initiate dictionary for interaction matrix
            try:
                self.species_dict = { s : int(float(s)) - 1 for s in self.community["ids"]}
            except:
                pass

        except ValueError as inst:
            msg = \
"""
  Attempting to set metacommunity size {} with {} species. This can happen sometimes
  with the `logser` metacommunity simulation. Simplest to ignore and rerun it.
"""
            raise MESSError(msg.format(self.paramsdict["S_m"], len(ids)))

        ### Fill the interaction martix ###
        ## Strength of the interaction
        self.interaction_matrix = np.zeros((self.paramsdict["S_m"],self.paramsdict["S_m"]))

            ## Begin with the interspecific interaction
        try:
            if type(self.paramsdict["intersp_competition"])==type((0,0)):
                ## Draw according to gamma distribution
                self.interaction_matrix = np.random.gamma(shape = self.paramsdict["intersp_competition"][0]
                                    , scale = self.paramsdict["intersp_competition"][1]
                                    , size = (self.paramsdict["S_m"],self.paramsdict["S_m"]))
            else:
                ## Expecting just a single value
                self.interaction_matrix = np.array([[self.paramsdict["intersp_competition"] for _ in range(self.paramsdict["S_m"])] for _ in range(self.paramsdict["S_m"])])

        except ValueError as inst:
            msg = \
"""
    Attempting to set interspecific interaction matrix with {} factors but it can only be 1 (exact value) or 2 terms (shape and scale parameters for gamma distribution
"""         
            raise MESSError(msg.format(self.paramsdict["intrasp_competition"]))


        try:
            ## Write over precedent terms for the diagonal of the matrix
            if type(self.paramsdict["intrasp_competition"])==type((0,0)):
                ## Draw according to gamma distribution
                rds = np.random.gamma(shape = self.paramsdict["intrasp_competition"][0]
                                    , scale = self.paramsdict["intrasp_competition"][1]
                                    , size = self.paramsdict["S_m"])
                for i in range(self.paramsdict["S_m"]):
                    self.interaction_matrix[i][i] = rds[i]
            else:
                ## Expecting just a single value
                for i in range(self.paramsdict["S_m"]):
                    self.interaction_matrix[i][i] = self.paramsdict["intrasp_competition"]
        except ValueError as inst:
            msg = \
"""
    Attempting to set intraspecific interaction matrix with {} factors but it can only be 1 (exact value) or 2 terms (shape and scale parameters for gamma distribution
"""         
            raise MESSError(msg.format(self.paramsdict["intrasp_competition"]))

        ## Nature of the interaction
        if type(self.paramsdict["intrasp_competition"])==type((0,0)) and type(self.paramsdict["intersp_competition"])==type((0,0)):
            try:
                factors = np.random.binomial(n = 1,
                    p = self.paramsdict["mutualism_proportion"],
                    size =(self.paramsdict["S_m"],self.paramsdict["S_m"]))
                factors[factors == 0] = -1
                self.interaction_matrix = self.interaction_matrix * factors
                # A positive value is a positive interaction,
                # A negative value is a negative interaction
            except ValueError as inst:
                raise MESSError("Error while drawing the nature of the interactions")

        elif type(self.paramsdict["intrasp_competition"])==type((0,0)):
            raise MESSError("Inter- and intra- specific competition have to be either both fixed or both gamma distributed (TODO)")
        elif type(self.paramsdict["intersp_competition"])==type((0,0)):
            raise MESSError("Inter- and intra- specific competition have to be either both fixed or both gamma distributed (TODO)")
        else:
            pass
            # No random proportion of mutualism if values fixed
                    

        ## Calculate immigration probabilities
        ## Here the true Jm under the logseries model will not equal the true
        ## sum of abundances, since it is a random variable per simulation.
        Jm = np.sum(self.community["abundances"])
        LOGGER.debug("Size of metacommunity - {}".format(Jm))
        self.community['immigration_probabilities'] = self.community["abundances"]/float(Jm)
        LOGGER.debug("Metacommunity info: shape {}\n[:10] {}".format(self.community.shape, self.community[:10]))



    def _update_species_pool(self, sname, trait_value):
        """
        Add a new species to the species pool. This is on account of speciation
        in the local communities and we need to keep track of the trait values
        globally. New species are appended to the end   and given dummy values
        for their immigration probability and regional  abundance. There are no
        dynamics in the metacommunity so these  new species will never act as
        colonists from the metacommunity. The form of this call is ugly for
        stupid reasons.
 
        :param str sname: The name of the new species to add to the
            metacommunity.
        :param float trait_value: The trait value of the new species.
        """
        try:
            LOGGER.debug("Adding species/trait_value - {}/{}".format(sname, trait_value))
            self.community = np.hstack((self.community,\
                                        np.array([tuple([sname, 0, 0, trait_value])], dtype=METACOMMUNITY_DTYPE)))
            self.trait_dict[sname] = trait_value
        except Exception as inst:
            LOGGER.error("Error in Metacommunity._update_species_pool - {}".format(inst))
            LOGGER.error("sname/trait_value - {}/{}".format(sname, trait_value))
            raise


    def _get_migrant(self):
        """
        Return one individual randomly sampled from the entire metacommunity.

        :return: A tuple containing the species ID (str) and the trait
            value (float).
        """
        migrant_draw = np.random.multinomial(1, self.community["immigration_probabilities"], size=1).argmax()
        new_species = self.community["ids"][migrant_draw]
        trait_value = self.community["trait_values"][migrant_draw]

        return new_species, trait_value


    def get_migrants(self, nmigrants=1):
        """
        Sample individuals from the Metacommunity. Each individual is
        independently sampled with replacement from the Metacommunity.

        :return: A tuple of lists of species IDs (str) and trait values (float).
        """
        migrants = np.array((nmigrants), dtype = object)
        trait_vals = np.array((nmigrants), dtype = float)
        for i in range(nmigrants):
            mig, trait = self._get_migrant()
            migrants[i] = mig
            trait_vals[i] = trait
        return migrants, trait_vals


#############################
## Metacommunity Parameter Info Dicts
#############################
LOCAL_PARAMS = {
    "S_m" : "Number of species in the regional pool",\
    "J_m" : "Total # of individuals in the regional pool",\
    "speciation_rate" : "Speciation rate of metacommunity",\
    "death_proportion" : "Proportion of speciation rate to be extinction rate",\
    "trait_rate_meta" : "Trait evolution rate parameter for metacommunity",\
    "ecological_strength" : "Strength of community assembly process on phenotypic change",\
    "intersp_competition" : "Strength of the interspecific competition (fixed single value or parameters for the gamma distribution)",\
    "intrasp_competition" : "Strength of the intraspecific competition (fixed single value or parameters for the gamma distribution)",\
    "mutualism_proportion" : "Percentage of inter- and intraspecific interaction terms that are beneficial for one the species (only used if inter- and intra- competition terms are not fixed)"
    }


if __name__ == "__main__":
    print("Test logser")
    data = Metacommunity("logser")
    print("{} {}".format(data, data.community[:10]))
    print("Test uniform")
    data = Metacommunity("uniform")
    print("{} {}".format(data, data.community[:10]))
    print("Test lognorm")
    data = Metacommunity("lognorm")
    print("{} {}".format(data, data.community[:10]))
    print("Test full file")
    data = Metacommunity("../SpInfo.txt")
    print("{} {}".format(data, data.community[:10]))

    for x in range(10):
        print(data._get_migrant())

    migs, traits = data.get_migrants(5)
    print(migs, traits)
