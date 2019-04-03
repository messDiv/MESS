from __future__ import print_function

import numpy as np
import pandas as pd
import math
from collections import Counter, OrderedDict
from scipy.stats import entropy, kurtosis, iqr, skew
from MESS.SGD import SGD


## Get one hill humber from a list of abundances (a column vector from the OTU table)
def hill_number(abunds, order):
    ## Make sure abunds is a np array or else order > 2 will act crazy
    abunds = np.array(abunds)
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

    dat = pd.DataFrame([], columns=["pis", "dxys", "abunds", "traits"])
    dat["pis"] = pis
    dat["dxys"] = dxys
    dat["abunds"] = abunds
    dat["traits"] = trts

    header = list(calculate_sumstats(dat, sgd_bins, sgd_dims, metacommunity_traits=metacommunity_traits).columns)
    return header


def calculate_sumstats(diversity_df, sgd_bins=10, sgd_dims=2, metacommunity_traits=None):

    moments = OrderedDict({"mean":np.mean,
                "std":np.std,
                "skewness":skew,
                "kurtosis":kurtosis,
                "median":np.median,
                "iqr":iqr})

    stat_dict = OrderedDict({})
    stat_dict["R"] = len(diversity_df)

    ## Abundance Hill #s
    for order in range(1,5):
        stat_dict["abund_h{}".format(order)] = hill_number(diversity_df["abunds"], order=order)

    ## pi/dxy stats
    for order in range(1,5):
        stat_dict["pi_h{}".format(order)] = hill_number(diversity_df["pis"], order=order)

    for name, func in moments.items():
        stat_dict["{}_pi".format(name)] = func(diversity_df["pis"])

    for name, func in moments.items():
        stat_dict["{}_dxys".format(name)] = func(diversity_df["dxys"])

    ## Tree stats
    stat_dict["trees"] = 0

    ## Trait stats
    for name, func in moments.items():
        stat_dict["{}_local_traits".format(name)] = func(diversity_df["traits"])

    if np.any(metacommunity_traits):
        for name, func in moments.items():
            val = func(metacommunity_traits)
            stat_dict["{}_regional_traits".format(name)] = val
            stat_dict["reg_loc_{}_trait_dif".format(name)] = val - stat_dict["{}_local_traits".format(name)]

    sgd = SGD(diversity_df["pis"],\
              diversity_df["dxys"],\
              nbins = sgd_bins, ndims = sgd_dims)

    stat_dict.update(sgd.to_dict())

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
