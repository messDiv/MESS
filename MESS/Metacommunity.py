#!/usr/bin/env python2.7

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
from MESS.util import MESSError,set_params

import logging
LOGGER = logging.getLogger(__name__)

## Limit on the number of redraws in the event of disallowed
## multiple migration, error out and warn if exceeded
MAX_DUPLICATE_REDRAWS_FROM_METACOMMUNITY = 1500

METACOMMUNITY_DTYPE = np.dtype([('ids', '|S10'),
                                ('immigration_probabilities', 'f8'),
                                ('abundances', 'i4'),
                                ('trait_values', 'f4')])

class Metacommunity(object):

    def __init__(self, meta_type="logser", quiet=False):
        self.quiet = quiet

        ## If you add a parameter to this dictionary you need
        ## to also add a short description to the LOCAL_PARAMS dict
        ## at the end of this file
        ##
        ## Also be sure to add it to _paramschecker so the type gets set correctly
        self.paramsdict = OrderedDict([
                        ("metacommunity_type", meta_type),
                        ("nspecies", 1000),
                        ("uniform_abund", 1000000),
                        ("logser_shape", 0.98),
                        ("J", 1000000),
        ])

        self.community = np.zeros([self.paramsdict["nspecies"]], dtype=METACOMMUNITY_DTYPE)

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

        ## Populate community with default values
        self.set_metacommunity()

        LOGGER.debug("Metacommunity paramsdict - {}".format(self.paramsdict))


    def __str__(self):
        return "<Metacommunity: {} Richness {}>".format(self.paramsdict["metacommunity_type"],\
                                                        self.paramsdict["nspecies"])


    def _paramschecker(self, param, newvalue, quiet=False):
        """ Raises exceptions when params are set to values they should not be"""
        ## TODO: This should actually check the values and make sure they make sense
        try:
            LOGGER.debug("set param {} - {} = {}".format(self, param, newvalue))
            if not quiet and MESS.__interactive__:
                print("  Updating Metacommunity parameters requires running set_metacommunity()"\
                        + " to apply the changes.")
            ## Cast params to correct types
            elif param == "metacommunity_type":
                self.paramsdict[param] = newvalue

            elif param == "nspecies":
                self.paramsdict[param] = int(float(newvalue))

            elif param == "uniform_abund":
                self.paramsdict[param] = int(float(newvalue))

            elif param == "logser_shape":
                self.paramsdict[param] = float(newvalue)

            elif param == "J":
                ## Do nothing. J is calculated from the data and not set, for now.
                pass

        except Exception as inst:
            ## Do something intelligent here?
            raise MESSError("Bad parameter {} - Bad value {}".format(param, newvalue))


    def write_params(self, outfile=None, append=True):
        """
        Write out the parameters for this island to a file.
        Normally this isn't called directly, but by the main
        simulation engine.

        append 
        """
        if outfile is None:
            raise MESSError("Metacommunity.write_params outfile must be specified.")

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
            
            header = "------- Metacommunity params: "
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


    def set_metacommunity(self, random=False):
        """
        For setting the metacommunity you can either generate a random
        uniform community or read on in from a file that's basically just
        a long list of abundances (as ints). Abundances are set from one
        of these locations then the species labels and immigration probs
        are calculated from there

        random=True will set random trait values in the range [0-1]
        """
        meta_type = self.paramsdict["metacommunity_type"]

        ## Accumulators for bringing in all the values. These will
        ## eventually all get shoved into self.community
        abundances = []
        ids = []
        trait_values = []

        ## Get Abundances
        if meta_type == "logser":
            abundances = logser.rvs(self.paramsdict["logser_shape"],\
                                        size=self.paramsdict["nspecies"])
        elif meta_type == "uniform":
            abundances = [self.paramsdict["uniform_abund"]]\
                               * self.paramsdict["nspecies"]
        elif os.path.isfile(meta_type):
            try:
                with open(meta_type, 'r') as infile:
                    lines = infile.readlines()
                    self.metcommunity_tree_height = float(lines[0].split()[0])
                    self.trait_evolution_rate_parameter = float(lines[1].split()[0])

                    for i in range(2,len(lines)):
                        info = lines[i].split()
                        ## Doing these somewhat out of order because i'm relying
                        ## on the info[1] failing and raising for the old type
                        ## file that is only a list of abundances. If you do the
                        ## ids[0] first it'll succeed and then fuck up the ids col
                        ## downstream.
                        trait_values.append(float(info[1]))
                        abundances.append(int(info[2]))
                        ids.append(info[0])
            except IndexError as inst:
                ## Could be an old style file just containing abundances, one per row
                try:
                    with open(meta_type, 'r') as infile:
                        lines = infile.readlines()
                        abundances = [int(line.split()[0]) for line in lines]
                except Exception as inst:
                    raise MESSError("  Malformed metacommunity specification file - {}\n    {}".format(meta_type, inst))

            ## If reading from a file then the number of species will not correspond
            ## with the value already in the paramsdict, so we need to update the nspecies count
            ## and reup the community ndarray
            LOGGER.debug("Read nspecies from file - {}".format(len(abundances)))
            self.paramsdict["nspecies"] = len(abundances)
            self.community = np.zeros([self.paramsdict["nspecies"]], dtype=METACOMMUNITY_DTYPE)
        else:
            raise MESSError("  Unrecognized metacommunity input - {}".format(meta_type))

        ## Get trait values per species
        ## TODO: Megan, is this supposed to be doing the same thing whether or not the
        ## random flag is set? seems lol
        ## Don't do random if we read in from the full_monty file
        if random and not trait_values:
            trait_values = np.random.rand(self.paramsdict["nspecies"])
        else:
            trait_values = np.random.rand(self.paramsdict["nspecies"])

        ## If ids haven't been assigned yet, do that here
        if not ids:
            ids = np.array(["t"+str(x) for x in range(0, self.paramsdict["nspecies"])])

        ## Populate the metacommunity ndarray
        self.community["abundances"] = np.array(abundances)
        self.community["ids"] = ids
        self.community['trait_values'] = np.array(trait_values)
        
        ## Calculate immigration probabilities
        self.paramsdict["J"] = np.sum(self.community["abundances"])
        LOGGER.debug("Size of metacommunity - {}".format(self.paramsdict["J"]))
        self.community['immigration_probabilities'] = self.community["abundances"]/float(self.paramsdict["J"])
        LOGGER.debug("Metacommunity info: shape {}\n[:10] {}".format(self.community.shape, self.community[:10]))

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

    def get_migrant(self):
        """ Return one 
        """
        migrant_draw = np.random.multinomial(1, self.community["immigration_probabilities"], size=1).argmax()
        new_species = self.community["ids"][migrant_draw]
        trait_value = self.community["trait_values"][migrant_draw]

        LOGGER.debug("Migrant idx {}\tid {}\t trait_val {}".format(migrant_draw, new_species, trait_value))
        return new_species, trait_value


    def get_nmigrants(self, nmigrants=1):
        migrants = []
        trait_vals = []
        for i in range(nmigrants):
            mig, trait = self.get_migrant()
            migrants.append(mig)
            trait_vals.append(trait)
        return migrants, trait_vals


#############################
## Metacommunity Parameter Info Dicts
#############################
LOCAL_PARAMS = {
    "metacommunity_type" : "Options: uniform/logser/<filename>",\
    "nspecies" : "Number of species in the regional pool",\
    "uniform_abund" : "If uniform: abundance per species",\
    "logser_shape" : "If logser: Shape parameter of the distribution",\
    "J" : "Total # of individuals in the regional pool (calculated)",
}



if __name__ == "__main__":
    data = Metacommunity("logser")
    print("{} {}".format(data, data.community[:10]))
    data = Metacommunity("uniform")
    print("{} {}".format(data, data.community[:10]))
    data = Metacommunity("../SpInfo.txt")
    print("{} {}".format(data, data.community[:10]))
    data = Metacommunity("../metacommunity_LS4.txt")
    print("{} {}".format(data, data.community[:10]))

    for x in xrange(10):
        print(data.get_migrant())

    data = set_params(data, "metacommunity_type", "logser")

    migs, traits = data.get_nmigrants(5)
    print(migs, traits)
