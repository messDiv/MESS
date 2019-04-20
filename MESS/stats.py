from __future__ import print_function

import numpy as np
import pandas as pd
import math
from collections import Counter, OrderedDict
from scipy.stats import entropy, kurtosis, iqr, skew
from sklearn.metrics import pairwise_distances
from MESS.SGD import SGD


def generalized_hill_number(abunds, vals=None, order=1):
    ## This is the Chao et al generalized Hill # formula
    ## Get one hill humber from a list of abundances (a column vector from the OTU table)
    ## Generalized function to calculate one hill number from a distribution
    ## of values of some statistic and abundances.
    ## vals:   A list of statistics calculated per species (e.g. pi)
    ## abunds: A list of corresponding abundances per species

    ## Degenerate edge cases can cause all zero values, particulary for pi
    ## in which case we bail out immediately
    if not np.any(abunds):
        return 0

    ## If vals is empty then populate the vector with a list of ones
    ## and this function collapses to the standard Hill number for abundance
    if vals is None:
        vals = np.ones(len(abunds))

    ## Make sure vals is a np array or else order > 2 will act crazy
    vals = np.array(vals)
    ## sum of values weighted by abundance
    V_bar = np.sum(vals*abunds)

    if order == 1:
        proportions = vals*(abunds/V_bar)
        h = np.exp(-np.sum(proportions * np.log(abunds/V_bar)))
    else:
        h = np.sum(vals*(abunds/V_bar)**order)**(1./(1-order))
    h = h/V_bar
    return h


def trait_hill_number(abunds, traits, order=1):
    ## Create the outer product of the abundances, flatten to a vector and divide
    ## by J**2 (essentially converting all pairwise abundance to proportion).
    pij = np.outer(abunds, abunds).flatten()/(np.sum(abunds)**2.)
    ## Reshape the np.array to make sklearn happy about it, then flatten it to a vector
    ## Then get the pairwise euclidean distances between all species trait values
    dij = pairwise_distances(traits.values.reshape(-1, 1)).flatten()
    return generalized_hill_number(abunds=pij, vals=dij, order=order)


## Get one hill humber from a list of abundances (a column vector from the OTU table)
def hill_number(abunds, order=0):
    ## Make sure abunds is a np array or else order > 2 will act crazy
    abunds = np.array(abunds)
    ## Degenerate edge cases can cause all zero values, particulary for pi
    if not np.any(abunds):
        return 0
    if order == 0:
        return len(np.nonzero(abunds)[0])
    if order == 1:
        h1 = np.exp(entropy(abunds))
        return h1
    tot = float(np.sum(abunds))
    proportions = np.array(abunds[abunds > 0])/tot
    prop_order = proportions**order
    h2 = np.sum(prop_order)**(1./(1-order))
    return h2


## Get all hill numbers from 0 to 'orders' from a column vector of abundances
def hill_numbers(abunds, orders, granularity=None, do_negative=False):
    ret = []
    min_order = 0
    if not granularity: granularity = orders + 1
    if do_negative:
        min_order = -orders
        granularity *= 2
    for order in np.linspace(min_order, orders, granularity):
        ret.append(hill_number(abunds, order))
    return np.array(ret)


def SAD(community, from_abundances=False, octaves=False, raw_abunds=False):
    """Generate the species abundance distribution either raw or in 
    octaves of powers of 2. The input here is a simple list of "individuals"
    specified just by their species identifier.
    """

    ## If input data is abundances per species and not the raw community
    ## data then we unpack it first
    if from_abundances:
        tmp = []
        for i, sp in enumerate(community):
            tmp.extend([i] * sp)
        community = tmp

    ## Make a counter for the local_community, counts the number of
    ## individuals w/in each species
    abundances = Counter(community)

    ## If we were doing mode=volcanic then there may be some remaining
    ## space in our carrying capacity that is unoccupied (indicated by
    ## None in the abundances.keys()
    try:
        abundances.pop(None)
    except KeyError:
        pass

    ## If raw_abunds then we just want to get the list of all abundances
    ## of all species as a list, don't do any SAD binning.
    if raw_abunds:
        return abundances.values()

    ## Now for each abundance class you have to go through and
    ## count the number of species at that abundance.
    abundance_distribution = Counter(list(abundances.values()))
    abundance_distribution = OrderedDict(sorted(abundance_distribution.items()))
    if octaves:
        dist_in_octaves = OrderedDict()
        minval = 1
        maxval = 2
        maxabund = max(abundance_distribution.keys())
        while maxval/2 < maxabund:
            ## Count up all species w/in each octave
            count = 0
            ## Here `i` is the abundance class and
            ## `j` is the count for that class
            for i, j in list(abundance_distribution.items()):
                if (i < maxval) and (i >= minval):
                    count += j
            dist_in_octaves[minval] = count
            minval = minval * 2
            maxval = maxval * 2
        abundance_distribution = dist_in_octaves
    return abundance_distribution


## TODO: Replace this and the dxy function with the pi/dxy 
## functions from easyCGD which are much nicer and smarter.
def pi(haplotypes):
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


def dxy(ihaps_t, mhaps_t):
    ## If no seg sites in a pop then haplotypes will be 0 length
    if ihaps_t.size == 0 or mhaps_t.size == 0:
        return 0

    ## Number of comparisons is == to n_island * n_metacommunity`
    ## n_metacommunity
    n_island = ihaps_t.shape[1]
    n_meta = mhaps_t.shape[1]
    n_comparisons = n_island * n_meta

    dxy = 0

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


def _get_sumstats_header(sgd_bins=10, sgd_dims=2, metacommunity_traits=None):
    ## Create some random data so the sumstats calculation doesn't freak out
    pis = np.random.random(10)/10
    dxys = np.random.random(10)/10
    abunds = np.random.randint(1, 100, 10)
    trts = np.random.random(10)*10
    meta_trts = np.random.random(10)*10

    dat = pd.DataFrame([], columns=["pi", "dxy", "abundance", "trait"])
    dat["pi"] = pis
    dat["dxy"] = dxys
    dat["abundance"] = abunds
    dat["trait"] = trts

    header = list(calculate_sumstats(dat, sgd_bins, sgd_dims, metacommunity_traits=metacommunity_traits).columns)
    return header


def calculate_sumstats(diversity_df, sgd_bins=10, sgd_dims=2, metacommunity_traits=None, verbose=False):

    moments = OrderedDict()
    for name, func in zip(["mean", "std", "skewness", "kurtosis", "median", "iqr"],\
                            [np.mean, np.std, skew, kurtosis, np.median, iqr]):
        moments[name] = func

    stat_dict = OrderedDict({})
    stat_dict["S"] = len(diversity_df)

    try:
        ## Abundance Hill #s
        for order in range(1,5):
            stat_dict["abund_h{}".format(order)] = hill_number(diversity_df["abundance"], order=order)
    except KeyError as inst:
        if verbose: print("  No abundance data present")

    try:
        for order in range(1,5):
            stat_dict["pi_h{}".format(order)] = hill_number(diversity_df["pi"], order=order)

        for name, func in moments.items():
            stat_dict["{}_pi".format(name)] = func(diversity_df["pi"])
    except KeyError as inst:
        if verbose: print("  No pi data present")

    try:
        for name, func in moments.items():
            stat_dict["{}_dxys".format(name)] = func(diversity_df["dxy"])
    except KeyError as inst:
        if verbose: print("  No dxy data present")

    ## Tree stats
    stat_dict["trees"] = 0

    try:
        for order in range(1,5):
            stat_dict["trait_h{}".format(order)] = trait_hill_number(diversity_df["abundance"],\
                                                                    diversity_df["trait"],\
                                                                    order=order)
        for name, func in moments.items():
            stat_dict["{}_local_traits".format(name)] = func(diversity_df["trait"])
    except:
        if verbose: print("  No trait data present")

    if np.any(metacommunity_traits):
        for name, func in moments.items():
            val = func(metacommunity_traits)
            stat_dict["{}_regional_traits".format(name)] = val
        ## Double for-loop here so the traits stay together in the sumstats output
        for name, func in moments.items():
            val = func(metacommunity_traits)
            stat_dict["reg_loc_{}_trait_dif".format(name)] = val - stat_dict["{}_local_traits".format(name)]

    try:
        sgd = SGD(diversity_df["pi"],\
                  diversity_df["dxy"],\
                  nbins = sgd_bins, ndims = sgd_dims)
        stat_dict.update(sgd.to_dict())
    except KeyError as inst:
        ## No Dxy data, try just py
        try:
            sgd = SGD(diversity_df["pi"],\
                      nbins = sgd_bins, ndims = 1)
            stat_dict.update(sgd.to_dict())
        except KeyError as inst:
            ## No py data. Skip SGD.
            pass

    return pd.DataFrame(stat_dict, index=[0])


if __name__ == "__main__":
    ## Test SAD()
    dat = [1,2,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]
    dat.extend([9] * 100)
    sad = SAD(dat)
    sad_oct = SAD(dat, octaves=True)
    #assert(cmp(sad.values(), [2, 2, 1, 1, 1, 1, 1]))
    #assert(cmp(sad_oct.values(), [2, 2, 1, 1, 2, 0, 1])) 
    print("SAD - {}".format(sad))
    print("SAD octaves - {}".format(sad_oct))

    pis = np.random.random(10)/10
    dxys = np.random.random(10)/10
    abunds = np.random.randint(1, 100, 10)
    trts = np.random.random(10)*10

    dat = pd.DataFrame([], columns=["pis", "dxys", "abunds", "traits"])
    dat["pis"] = pis
    dat["dxys"] = dxys
    dat["abunds"] = abunds
    dat["traits"] = trts

    ss = calculate_sumstats(dat, sgd_bins=10, sgd_dims=1)
    print("sumstats", len(ss.values), ss)
    header = _get_sumstats_header(sgd_bins=10, sgd_dims=1)
    print("header", len(header), header)
