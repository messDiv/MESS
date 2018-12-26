#!/usr/bin/env python2.7

from scipy.stats import hmean
import numpy as np
import pandas as pd
import collections      # For Counter
import itertools
import msprime
import math
import os

from .stats import *

import logging
LOGGER = logging.getLogger(__name__)

class species(object):

    def __init__(self, name = "", species_params={"mutation_rate":0.00000222,
                                                  "sigma":1000,
                                                  "sequence_length":570},
                growth="constant", abundance = 1,
                meta_abundance = 1, divergence_time = 0, migration_rate=0,
                abundance_through_time={}):

        ## I'm calling anything that is invariant across species
        ## a 'parameter' even though I know migration rate is a
        ## parameter too
        self.paramsdict = dict([
                            ("sigma", species_params["sigma"]),
                            ("sequence_length", species_params["sequence_length"]),
                            ("mutation_rate", species_params["mutation_rate"]),
                            ("sample_size_local", 10),
                            ("sample_size_meta", 10),
                            ("abundance_through_time", collections.OrderedDict(sorted(abundance_through_time.items()))),
        ])

        ## Stats for this species in a nice pd series
        ## pi_tot is pi calculated for local and metacommunity combined
        ## da is what was called pi_net in msbayes
        self.stats = pd.Series(
            index=[ "name",
                    "tdiv",
                    "migration_rate",
                    "growth_rate",
                    "Ne_local",
                    "Ne_meta",
                    "segsites_tot",
                    "segsites_local",
                    "segsites_meta",
                    "pi_tot",
                    "pi_local",
                    "pi_meta",
                    "da",
                    "dxy",
                    "tajd_local"]).astype(np.object)

        self.stats["name"] = name
        self.stats["Ne_local"] = abundance * self.paramsdict["sigma"]
        self.stats["Ne_meta"] = meta_abundance * self.paramsdict["sigma"]
        self.stats["tdiv"] = divergence_time * self.paramsdict["sigma"]

        ## Parameters
        self.name = name
        self.tree_sequence = []

        ## Calculate the growth rate and the initial population size
        if growth == "exponential":
            ## TODO: Add a hackersonly param to tune the number of founders
            ## or maybe just use mig_clust_size?
            initial_size = 1.
            self.stats["growth_rate"] = -np.log(initial_size/self.stats["local_Ne"])/self.stats["tdiv"]
        ## Instantaneous expansion to current size and constant size through time
        elif growth == "constant":
            self.stats["growth_rate"] = 0
        ## Take the harmonic mean of the abundance trajectory through time
        elif growth == "harmonic":
            try:
                harmonic = hmean(np.array(abundance_through_time.values()) * self.paramsdict["sigma"])
                LOGGER.debug("sname Ne hmean- {} {} {}".format(self.name, self.stats["Ne_local"], harmonic))
                self.stats["Ne_local"] = harmonic
                self.stats["growth_rate"] = 0
            except Exception as inst:
                LOGGER.debug("Failed harmonic mean for {}".format(self))
                raise inst
        ## Do piecewise constant population size changes through time
        elif growth == "piecewise":
            ##TODO: piecewise population size change is unimplemented
            pass
        else:
            raise MESSError("Unrecognized population growth parameter - {}".format(growth))

        ## Set migration rate here to account for constant vs harmonic mean of Ne
        self.stats["migration_rate"] = migration_rate / self.stats["Ne_local"]

    @classmethod
    def from_df(self, df):
        """ Class method for making creating a species object somewhat less annoying.
            The dataframe passed in must contain at least the following information.
            Colonization times and times of abundances through time need to be in
            backward time units.

            colonization_times               13623
            post_colonization_migrants           0
            abundances_through_time     {13623: 1}
            ancestor                              
            meta_abund                        1000
            local_abund                          2"""

        ## TODO: this doesn't pull in the growth model so will always be constant
        ## I guess you can set it yourself downstream
        ## TODO: This also isn't even close to right. Colonization time is still forward time,
        ## and consequently both divergence time and migration rate aren't scaled / K
        name = df.columns[0]
        sp_info = df[name]
        sp = species(name = name,
                     divergence_time = sp_info["colonization_times"],\
                     abundance = sp_info["local_abund"],\
                     meta_abundance = sp_info["meta_abund"],\
                     migration_rate = sp_info["post_colonization_migrants"]/float(sp_info["colonization_times"]),\
                     abundance_through_time = sp_info["abundances_through_time"])
        return sp


    def __str__(self):
        return self.stats.to_string()


    def __repr__(self):
        return self.__str__()


    def _get_local_configuration(self):
        pop_local = msprime.PopulationConfiguration(\
                        sample_size = self.paramsdict["sample_size_local"],\
                        initial_size = self.stats["Ne_local"],
                        growth_rate = self.stats["growth_rate"])
        return pop_local


    def _get_meta_configuration(self):
        pop_meta = msprime.PopulationConfiguration(\
                        sample_size=self.paramsdict["sample_size_meta"],\
                        initial_size=self.stats["Ne_meta"])
        return pop_meta


    def _get_size_changes_through_time(self, pop_id = 0):
        size_change_events = []

        for time, size in self.paramsdict["abundance_through_time"].items():

            local_size_change = msprime.PopulationParametersChange(\
                                            time = time,\
                                            initial_size = size,\
                                            population_id = pop_id)
            size_change_events.append(local_size_change)

        return size_change_events


    def _get_local_meta_split(self, founder_idx = 0, meta_idx = 1):
        ## Going backwards in time, at colonization time throw all lineages from
        ## the local community back into the metacommunity
        split_event = msprime.MassMigration(time = self.stats["tdiv"] + 1,\
                                            source = founder_idx,\
                                            destination = meta_idx,\
                                            proportion = 1)

        local_rate_change = msprime.PopulationParametersChange(\
                                            time = self.stats["tdiv"],\
                                            growth_rate = 0,\
                                            population_id = 0)

        ## TODO: Could mess with 'initial_size' here, but you don't want
        ## to sample too much from the metacommunity or the local pi
        ## goes way up.
        local_size_change = msprime.PopulationParametersChange(\
                                            time = self.stats["tdiv"],\
                                            initial_size = .01,\
                                            population_id = founder_idx)

        migrate_change = msprime.MigrationRateChange(
                                            time = self.stats["tdiv"],\
                                            rate = 0)

        return [migrate_change, local_size_change, local_rate_change, split_event]


    def simulate_seqs(self):
        LOGGER.debug("Entering simulate_seqs - {}".format(self.name))
        ## TODO: Here we are assuming only one island and that the migration
        ## is only ever unidirectional from the mainland to the island
        migmat = [[0, self.stats["migration_rate"]],
                    [0, 0]]

        pop_local = self._get_local_configuration()
        pop_meta = self._get_meta_configuration()

        split_events = self._get_local_meta_split()

        ## Useful for debugging demographic events.
        if LOGGER.getEffectiveLevel() == 10:
            debug = msprime.DemographyDebugger(population_configurations = [pop_local, pop_meta],\
                                               migration_matrix = migmat,\
                                               demographic_events = split_events)

            ## Enable this at your own peril, it will dump a ton of shit to stdout
            debug.print_history()
            ## I briefly toyed with the idea of logging this to a file, but you really
            ## don't need it that often, and it'd be a pain to get the outdir in here.
            #debugfile = os.path.join(self._hackersonly["outdir"],
            #                    self.paramsdict["name"] + "-simout.txt")
            #with open(debugfile, 'a') as outfile:
            #    outfile.write(debug.print_history())

        LOGGER.debug("Executing msprime - {}".format(self.name))
        self.tree_sequence = msprime.simulate(length = self.paramsdict["sequence_length"],\
                                              mutation_rate = self.paramsdict["mutation_rate"],\
                                              population_configurations = [pop_local, pop_meta],\
                                              migration_matrix = migmat,\
                                              demographic_events = split_events)
        self.get_sumstats()

        tree = self.tree_sequence.first()
        colour_map = {0:"red", 1:"blue"}
        node_colours = {u: colour_map[tree.population(u)] for u in tree.nodes()}
        node_labels = {
            u: (str(u) if u < 8 else "{} (t={:.2f})".format(u, tree.time(u))) 
            for u in tree.nodes()}
        #tree.draw(path="/tmp/{}.svg".format(self.name.replace(" ", "_")), height=500, width=1000, node_labels=node_labels, node_colours=node_colours)


    def get_sumstats(self, samps=np.array([]), metasamps=np.array([])):

        LOGGER.debug("Entering get_sumstats - {}".format(self.name))

        self.stats["segsites_tot"] = len(next(self.tree_sequence.haplotypes()))
        LOGGER.debug("segsites_tot {}".format(self.stats["segsites_tot"]))

        all_haps = self.tree_sequence.haplotypes()
        ## Get population specific haplotypes
        if samps.size > 0 and metasamps.size > 0:
            LOGGER.debug("samps specified - samps and metasamps - {} {}".format(samps, metasamps))
            ## pairwise diversity per base for local and meta samples combined
            self.stats["pi_tot"] = self.tree_sequence.get_pairwise_diversity(np.concatenate([samps, metasamps]))\
                                    / self.paramsdict["sequence_length"]
            all_haps = list(all_haps)
            island_haps = [all_haps[x] for x in samps]
            meta_haps = [all_haps[x] for x in metasamps]
        elif samps.size == 0 and samps.size == 0:
            LOGGER.debug("samps empty - samps and metasamps - {} {}".format(samps, metasamps))
            self.stats["pi_tot"] = self.tree_sequence.get_pairwise_diversity()\
                                    / self.paramsdict["sequence_length"]
            LOGGER.debug("pi_tot {}".format(self.stats["pi_tot"]))
            try:
                island_haps = [next(all_haps) for _ in range(self.paramsdict["sample_size_local"])]
                meta_haps = [next(all_haps) for _ in range(self.paramsdict["sample_size_meta"])]
            except Exception as inst:
                raise MESSError("Error fetching local/meta haplotypes {}".format(self.name))
        else:
            raise Exception("Problem in get_sumstats: Samps and metasamps must either both be empty or both contain a list of samples")

        ## Calculate S for each population
        ## Turn to an np array and transpose it, then test if all elements
        ## are the same by making a set from each row and checking if len > 1
        ## In the transposed array rows are bases and columns are individuals
        ## There's probably a much fucking smarter way to do this
        ihaps_t = np.transpose(np.array([map(int, list(x)) for x in island_haps]))
        mhaps_t = np.transpose(np.array([map(int, list(x)) for x in meta_haps]))

        ## Counter makes a dict, so just get the counts for 2, which indicates
        ## sites segregating in the pop
        ## S will not always == S_local + S_meta. If a site is fixed in one pop and not
        ## present in the other then S will be less than the total. If a site is segragating
        ## in both pops then S will be greater than the total.
        ## There's probably a smarter way to do this too....
        self.stats["segsites_local"] = collections.Counter(\
                                        [len(set(ihaps_t[x])) for x in range(len(ihaps_t))])[2]
        self.stats["segsites_meta"] = collections.Counter(\
                                        [len(set(mhaps_t[x])) for x in range(len(mhaps_t))])[2]

        ## Pass in the transposed arrays, since we already have them
        self.stats["pi_local"] = get_pi(ihaps_t) / self.paramsdict["sequence_length"]
        self.stats["pi_meta"] = get_pi(mhaps_t) / self.paramsdict["sequence_length"]

        ## get pairwise differences between populations while ignoring differences
        ## within populations (Dxy)
        self.stats["dxy"] = get_dxy(ihaps_t, mhaps_t) / self.paramsdict["sequence_length"]

        self.stats["da"] = self.stats["dxy"] - (self.stats["pi_local"] + self.stats["pi_meta"])/2

        ## Forbid biologically unrealistic values of pi
        ## TODO: This is dumb
        if self.stats["pi_meta"] > 0.2 or self.stats["pi_local"] > 0.2:
            LOGGER.error("Bad pi {}".format(self))
            #self.simulate_seqs()
            #self.get_sumstats()

        self.stats["tajd_local"] = tajD_island(island_haps, self.stats["segsites_local"])


    ## This is hackish and is used by the LocalCommunity.bottleneck() function
    ## that is more or less untested
    def update_abundance(self, abund):
        self.stats["Ne_local"] = abund * self.paramsdict["sigma"]


#######################################
## Various summary statistics
#######################################
def tajD_island(haplotypes, S):
    if len(haplotypes) == 0:
        return 0
    if not any(haplotypes):
        return 0
    if S == 0:
        return 0
    d_num = pairwise_diffs(haplotypes) - watt_theta(len(haplotypes), S)
    ddenom = tajD_denom(len(haplotypes), S)
    if ddenom == 0:
        D = 0
    else:
        D = d_num/ddenom
        #print("nhaps {} nuniq {} S {} D {} num {} denom {}".format(len(haplotypes), len(set(haplotypes)), S, D, d_num, ddenom))
    return D


def pairwise_diffs(haplotypes):
    tot = 0
    for count, i in enumerate(itertools.combinations(haplotypes,2)):
        tot += sum(a != b for (a,b) in zip(i[0], i[1]))
    ## Number of comparisons is count + 1 bcz enumerate starts at 0
    return tot/float(count+1)


def watt_theta(n, S):
    return S/sum([1./x for x in xrange(1,n)])


## Fuckin helps if you do it right. This page has a nice worked example with values for each
## subfunction so you can check your equations:
## https://ocw.mit.edu/courses/health-sciences-and-technology/hst-508-quantitative-genomics-fall-2005/study-materials/tajimad1.pdf
def tajD_denom(n, S):
    b1 = (n+1)/float(3*(n-1))
    a1 = sum([1./x for x in xrange(1, n)])
    c1 = b1 - (1./a1)
    e1 = c1/a1
    a2 = sum([1./(x**2) for x in xrange(1, n)])
    b2 = (2.*(n**2 + n + 3))/(9*n*(n-1))
    c2 = b2 - (n+2)/(a1*n) + (a2/(a1**2))
    e2 = c2/(a1**2+a2)
    ddenom = math.sqrt(e1*S + e2*S*(S-1))
    return ddenom


def get_pi(haplotypes):
    ## If no seg sites in a pop then haplotypes will be 0 length
    if haplotypes.size == 0:
        return 0
    n = len(haplotypes[0])
    n_comparisons = float(n) * (n - 1) / 2

    pi = 0
    for hap in haplotypes:
        k = np.count_nonzero(hap)
        pi += float(k) * (n - k) / n_comparisons
    return(pi)


def get_dxy(ihaps_t, mhaps_t):
    ## If no seg sites in a pop then haplotypes will be 0 length
    if ihaps_t.size == 0 or mhaps_t.size == 0:
        return 0

    ## Number of comparisons is == to n_island * n_metacommunity`
    ## n_metacommunity
    n_island = ihaps_t.shape[1]
    n_meta = mhaps_t.shape[1]
    n_comparisons = n_island * n_meta

    dxy = 0
    len_seq = ihaps_t.shape[0]

    ## ibases and mbases are now a list of all bases at a particular
    ## site within each population
    ## There's probably a more elegant way to do this but I was
    ## gunning for readibility. Probably failed.
    for ibases, mbases in zip(ihaps_t, mhaps_t):
        nonzeros_island = np.count_nonzero(ibases)
        nonzeros_meta = np.count_nonzero(mbases)
        zeros_island = n_island - nonzeros_island
        zeros_meta = n_meta - nonzeros_meta

        dxy += (nonzeros_island * zeros_meta \
                + zeros_island * nonzeros_meta) / float(n_comparisons)
    return dxy


####################################################
## Function to build the msprime model for clades of
## species related in the local community.
## TODO: This lives here now because I can't think
## of a better place to put it...
####################################################
def sim_clade(forward_info):
    ## Build the nasty msprime command for species all related
    ## by local speciation. Run the msprime, and gather summary
    ## statistics

    ## Here local info has to be transformed so that all time
    ## values are backards in time.
    ## TODO:

    ## Get species objects for each species in this clade
    clade_species = [species.from_df(pd.DataFrame(forward_info[x])) for x in forward_info]
    LOGGER.debug("Got clade_species:\n{}".format(clade_species))

    


if __name__ == "__main__":
    import MESS
    sp = species("wat", abundance=100, meta_abundance=1000, divergence_time=1000000)
    sp.simulate_seqs()
    print(sp)

    #data = MESS.Region("wat")
    #data.add_local_community("tmp", 500, 0.5)
    #data.simulate(_lambda=0.5)
