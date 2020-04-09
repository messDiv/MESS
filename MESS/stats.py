

import numpy as np
import pandas as pd
import math
from collections import Counter, OrderedDict
from itertools import combinations
from scipy.stats import entropy, kurtosis, hmean, iqr, skew, spearmanr
from sklearn.metrics import pairwise_distances
from MESS.SGD import SGD

import MESS

data_axes = ["abundance", "pi", "dxy", "trait"]


def generalized_hill_number(abunds, vals=None, order=1, scale=True, verbose=False):
    """
    This is the Chao et al (2014) generalized Hill # formula. Generalized
    function to calculate one Hill number from a distribution of values of
    some statistic and abundances.

    :param array-like abunds: An `array-like` of abundances per species.
    :param array-like vals: An `array-like` of values per species (e.g. 
        pi values per species, or trait values per species). If this parameter
        is empty then only abundance Hill numbers are calculated.
    :param float order: The Hill number to calculate. 0 is species richness.
        Positive values calculate Hill numbers placing increasing weight on
        the most abundant species. Negative values can also be specified
        (placing more weight on the rare species), but these are uncommonly
        used in practice.
    :param bool scale: Whether to scale to effective numbers of species, or
        return the raw attribute diversity. Equivalent to equation 5c in
        Chao et al 2014. You will almost never want to disable this.

    :return float: The generalized Hill number of order `order` for the given
        data axis using the formula proposed by Chao et al (2014).
    """
    ## Degenerate edge cases can cause all zero values, particulary for pi
    ## in which case we bail out immediately
    if not np.any(abunds):
        return 0

    ## Be sure abundance is scaled to relative abundance and convert to np
    abunds = np.array(abunds)/np.sum(abunds)

    ## If vals is empty then populate the vector with a list of ones
    ## and this function collapses to the standard Hill number for abundance
    if vals is None:
        vals = np.ones(len(abunds))

    ## Make sure vals is an np array or else order > 2 will act crazy
    vals = np.array(vals)
    if verbose: print(("sums:", "dij", np.sum(vals), "pij", np.sum(abunds)))
    ## sum of values weighted by abundance
    V_bar = np.sum(vals*abunds)
    if verbose: print(("vbar", V_bar))

    ## Use the special formula for order = 1
    if order == 1:
        proportions = vals*(abunds/V_bar)
        h = np.exp(-np.sum(proportions * np.log(abunds/V_bar)))
    else:
        h = np.sum(vals*(abunds/V_bar)**order)**(1./(1-order))
    if scale: h = h/V_bar
    return h


def trait_hill_number(abunds, traits, order=1, verbose=False):
    """
    Calculate trait Hill numbers using the generalized Hill number equation.
    For this function `abunds` and `traits` should be array-likes of equal
    length where the abundance of the `i`-th element corresponds to the
    trait value of the `i`-th species.

    :param array-like abunds: Abundances per species in the local community.
    :param array-like traits: Corresponding trait values for each species.
    :param float order: Order of the Hill number to calculate.
    :param bool verbose: Whether to print informational messages.

    :return float: The Hill number of the trait values for this community
        using the generalized forumla proposed by Chao et al (2014).
    """
    ## If there's only one species in the community then the pairwise_distances
    ## function will return 0, and generalized_hill will return inf,
    ## and you'll get nans in the final output. Not good. If only one species
    ## in local, then just return 1 for trait hills. This agrees with the other
    ## hill functions calc'd on just one species.
    if len(traits) == 1:
        return 1
    ## Create the outer product of the abundances, flatten to a vector and divide
    ## by J**2 (essentially converting all pairwise abundance to proportion).
    pij = np.outer(abunds, abunds).flatten()/(np.sum(abunds)**2.)
    ## Reshape the np.array to make sklearn happy about it, then flatten it to a vector
    ## Then get the pairwise euclidean distances between all species trait values
    dij = pairwise_distances(traits.values.reshape(-1, 1)).flatten()
    return generalized_hill_number(abunds=pij, vals=dij, order=order, verbose=verbose)**(1/2.)


def hill_number(abunds, order=0):
    """
    Get one hill humber from an array-like of abundances.

    :param array-like abunds: An `array-like` of counts of individuals per
        species.
    :param float order: The order of the Hill number to calculate.

    :return: The Hill number of order `order`.
    """
    ## Make sure abunds is an np array or else order > 2 will act crazy
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


def hill_numbers(abunds, orders, granularity=None, do_negative=False):
    """
    Get all hill numbers from 0 to 'orders' from an array-like of abundances.
    If `granularity` is specified then calculate fractional Hill numbers
    such that the returned `numpy.array` contains `granularity` count of Hill
    numbers calculated and equally spaced between 0 and `orders`. The
    `do_negative` parameter will include both positive and negative values
    of Hill numbers up to order equal to `orders`.

    :param array-like abunds: An `array-like` of counts of individuals per
        species.
    :param int orders: The max order of Hill numbers to calculate.
    :param int granularity: The number of equally spaced fractional Hill
        numbers to calculate between 0 and `orders`. If not specified then
        returns only integer values of Hill numbers between 0 and `orders`.
    :param bool do_negative: Whether to calculate negative as well as positive
        Hill numbers between 0 and `orders`.

    :return numpy.array: An array of Hill numbers from 0 to (+/-) `orders`
        in intervals of defined by the `granularity` parameter.
    """
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
    """
    Generate the species abundance distribution either raw or in octaves of
    powers of 2. The input here is a simple list of individuals specified
    as species names.

    .. note:: The `from_abundances` and `raw_abunds` arguments are probably
        not useful and may be deprecated in the future.

    :param list community: The list of individual in the local community.
        This list should be of length `J` and should contain only the species
        id of each individual in the local community.
    :param bool from_abundances: If the input community list is already
        counts of individuals per species then unpack it into a list of
        individuals per species.
    :param bool octaves: Whether to return the SAD in octaves of powers of 2.
    :param bool raw_abunds: Wether to return the SAD binned into abundance
        classes or simply return the list of abundances per species.

    :return OrderedDict: An ordered dictionary with keys equal to abundance
        classes and values equal to counts of species per class.
    """

    ## If input data is abundances per species and not the raw community
    ## data then we unpack it first
    if from_abundances:
        tmp = []
        for i, sp in enumerate(community):
            try:
                int(sp)
            except:
                # If sp won't cast to int it's bad data so skip it.
                continue
            tmp.extend([i] * int(sp))
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
        return list(abundances.values())

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


def _hnum(N):
    """
    The N-th harmonic number. Used by the Watterson function.

    :param int N: The harmonic number to calculate.
    """
    return (N)*(1./hmean(list(range(1, N+1))))


def _n_segsites_pooled(seq):
    ## This is pooled data so we are looking at one sequence per
    ## species. Only count IUPAC ambiguity codes.
    c = Counter(seq)
    try:
        for i in ["A", "C", "G", "T", "-", "?", "N", "\n"]:
            c.pop(i)
    except:
        pass
    return sum(c.values())


def Watterson(seqs, nsamples=0, per_base=True):
    """
    Calculate Watterson's theta and optionally average over sequence length.

    :param str/array-like seqs: The DNA sequence(s) over which to 
        calculate the statistic. This parameter can be a single DNA sequence
        as a string, in which case we assume it is pooled data with IUPAC
        ambiguity codes indicating segregating sites. It can also be a string
        indicating the path to a fasta file, or a list of sequences which
        may or may not include the sample names (they will be removed).
    :param int nsamples: The number of samples in the pooled sequence
        (for pooled data only).
    :param bool per_base: Whether to average over the length of the sequence.

    :return float: The value of Watterson's estimator of theta per base.
    """
    if isinstance(seqs, str):
        try:
            segsites = _n_segsites_pooled(seqs)
        except Exception as inst:
            try:
                ## If seqs is a file just read in the data here and do the
                ## calculation in the block below
                seqs = open(seqs).readlines()
            except (IOError, TypeError) as inst:
                ## Not a file
                pass

    if isinstance(seqs, list):
        ## Remove all the fasta names, if any, convert to a numpy array
        ## and transpose, so now it's a list of lists of bases per sample.
        seqs = np.array([list(x) for x in seqs if not ">" in x]).T
        ## Get the number of samples
        nsamples = seqs.shape[1]
        ## Count the number of non-unique bases (the -1 turns monomorphic sites
        ## into zeros.
        #Py2 segsites = np.count_nonzero(np.array(map(len, map(lambda x: set(x), seqs))) - 1)
        segsites = np.count_nonzero(np.array(list(map(len, [set(x) for x in seqs]))) - 1)

    if nsamples == 0:
        raise MESS.util.MESSError("If using pooled data you must include the `nsamples` parameter.")

    w_theta = segsites/_hnum(nsamples)

    if per_base:
        w_theta = w_theta/float(len(seqs))

    return w_theta


def _get_sumstats_header(sgd_bins=10, sgd_dims=1, metacommunity_traits=None):
    """
    Return the full list of supported summary statistics. This is an internal
    function used to generate, among other things, the SIMOUT header.

    :param int sgd_bins: The number of bins per axis of the SGD.
    :param int sgd_dims: The number of dimensions of the SGD (1 or 2)
    :param array-like metacommunity_traits: Metacommunity trait values.
    """
    ## Create some random data so the sumstats calculation doesn't freak out
    pis = MESS.rng.rng.random(10)/10
    dxys = MESS.rng.rng.random(10)/10
    abunds = MESS.rng.rng.integers(1, 100, 10)
    trts = MESS.rng.rng.random(10)*10
    meta_trts = MESS.rng.rng.random(10)*10

    dat = pd.DataFrame([], columns=["pi", "dxy", "abundance", "trait"])
    dat["pi"] = pis
    dat["dxy"] = dxys
    dat["abundance"] = abunds
    dat["trait"] = trts
    header = list(calculate_sumstats(dat, sgd_bins, sgd_dims, metacommunity_traits=metacommunity_traits).columns)
    return header


def calculate_sumstats(diversity_df, sgd_bins=10, sgd_dims=1, metacommunity_traits=None, verbose=False):
    """
    Calculate all summary statistics on a dataset composed of one or more of
    the target data axes. This function will automatically detect the
    appropriate set of summary statistics to calculate based on the columns
    of the dataset. The passed in `diversity_df` may contain one or more of
    the following data axes:

    * ``abundance``: Abundances per species as counts of individuals.
    * ``pi``: Nucleotide diversity per base per species.
    * ``dxy``: Absolute divergence between each species in the local community and the sister species in the metacommunity.
    * ``trait``: The trait value of each species. Trait values are continuous and the distribution of trait values in the local community should be zero centered.

    .. note:: This method should be used for calculating summary statistics for
        all empirical datasets as this is the method that is used to generate
        summary statistics for the simulations. This guarantees that observed
        and simulated statistics are calculated identically.


    :param pandas.DataFrame diversity_df: A DataFrame containing the empirical
        data. This df has a very specific format which is specified above.
    :param int sgd_bins: The number of bins per axis of the constructed SGD.
        This must match the number of bins specified for simulations.
    :param int sgd_dims: The number of dimensions of the constructed SGD.
        This value can be either 1 (pi only) or 2 (both pi and dxy). This
        parameter must match the number of dimensions specified for simulations.
    :param array-like metacommunity_traits: A list or np.array of the trait values
        from the metacommunity. Used for calculating some of the trait based
        summary statistics. These values take the same form as the values of the
        local trait data (i.e. they should be continuous and zero centered).
    :param bool verbose: Whether to print some informational messages.

    :return: Returns a pandas.Dataframe with one row containing all of the
        applicable summary statistic for the input dataframe.
    """
    ## Check diversity_df formatting
    axes = ["abundance", "pi", "dxy", "trait"]
    if len(set(diversity_df.columns).intersection(set(axes))) == 0:
        raise MESS.util.MESSError("\nColumn names must be one or more of exactly: {}\nYou had: {}".format(axes, list(diversity_df.columns)))

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

        for name, func in list(moments.items()):
            stat_dict["{}_pi".format(name)] = func(diversity_df["pi"])
    except KeyError as inst:
        if verbose: print("  No pi data present")

    try:
        for name, func in list(moments.items()):
            stat_dict["{}_dxys".format(name)] = func(diversity_df["dxy"])
    except KeyError as inst:
        if verbose: print("  No dxy data present")

    ## Tree stats
    stat_dict["trees"] = 0

    try:
        for order in range(1,5):
            stat_dict["trait_h{}".format(order)] = trait_hill_number(diversity_df["abundance"],\
                                                                    diversity_df["trait"],\
                                                                    order=order,\
                                                                    verbose=verbose)
    except:
        if verbose: print("  Trait hill numbers ignored since abundances not present.")
    finally:
        try:
            for name, func in list(moments.items()):
                stat_dict["{}_local_traits".format(name)] = func(diversity_df["trait"])
            if np.any(metacommunity_traits):
                for name, func in list(moments.items()):
                    val = func(metacommunity_traits)
                    stat_dict["{}_regional_traits".format(name)] = val
                ## Double for-loop here so the traits stay together in the sumstats output
                for name, func in list(moments.items()):
                    val = func(metacommunity_traits)
                    stat_dict["reg_loc_{}_trait_dif".format(name)] = val - stat_dict["{}_local_traits".format(name)]
        except:
            if verbose: print("  No trait data present")


    ## Calculate correlations between all the distributions we have
    ## This may not make much sense for the trait data at this point.
    ## Here we take the spearmanr()[0] because it returns a tuple of
    ## corellation and p-value.
    ##
    ## TODO: If the variance of either of the variables is 0, e.g.
    ## if all pi values are 0, then spearmanr is undefined and returns nan.
    ## Here we'll just call this 0, even though it's not technically correct
    ## it's close enough to what we want. Would be better to do something smrter.
    ## Only look at the valid columns, in case the input df has some weird shit
    ## going on.
    valid = set(["abundance", "pi", "dxy", "trait"])
    for pair in combinations(sorted(set(diversity_df.columns).intersection(valid)), r=2):
        tmp_df = diversity_df.copy()
        ## If doing traits then transform the trait values into distance from
        ## local trait mean. Should see positive correlation in filtering and
        ## negative in competition.
        if "trait" in pair:
            #idx = pair.index("trait")
            tmp_df["trait"] = np.abs(tmp_df["trait"] - tmp_df["trait"].mean())
        cor = spearmanr(tmp_df[pair[0]], tmp_df[pair[1]])[0]
        if np.isnan(cor): cor = 0
        stat_dict["{}_{}_cor".format(pair[0], pair[1])] = cor

    try:
        sgd = SGD(diversity_df["pi"],\
                  diversity_df["dxy"],\
                  nbins = sgd_bins, ndims = sgd_dims)
        stat_dict.update(sgd.to_dict())
    except KeyError as inst:
        ## No Dxy data, try just pi
        try:
            sgd = SGD(diversity_df["pi"],\
                      nbins = sgd_bins, ndims = 1)
            stat_dict.update(sgd.to_dict())
        except KeyError as inst:
            ## No py data. Skip SGD.
            pass

    return pd.DataFrame(stat_dict, index=[0])


def feature_sets(empirical_df=None):
    """
    Convenience function for getting all the different combinations of summary
    statistics (features) that are relevant for the axes available in the
    empirical data. For example, if you have abundance and pi in your observed
    data then this function will return 3 feature sets relevant to: using only
    abundance as the data, using only pi as the data, and using both abundance
    and pi.

    :param pandas.DataFrame empirical_df: A DataFrame containing the empirical
        data. This df has a very specific format which is documented here. If
        none is provided then summary statistics for all possible combinations
        of data axes are returned

    :return: A dictionary with keys being string descriptors the focal data axes
        and values being a list of summary statistics. Data axes that encompass
        more than one axis will be keyed by a single string with the name of
        both axes concatenated by a '+' sign (e.g. "abundance+pi+trait" is the
        key that will give summary statistics relevant for all three of these
        data axes.
    """
    empirical_axes = ["abundance", "pi", "dxy", "trait"]
    if empirical_df is not None:
        empirical_axes = empirical_df.columns

    ## Pass in a dummy list of metacommnity traits, or else the call to
    ## _get_sumstats_header won't return metacommunity traits.
    sumstats = _get_sumstats_header(metacommunity_traits=MESS.rng.rng.rand(10))
    feature_sets = OrderedDict({})
    feature_sets["abundance"] = ["S"] + [x for x in sumstats if "abund" in x and "_cor" not in x]
    feature_sets["pi"] = [x for x in sumstats if ("pi" in x or "SGD" in x) and "_cor" not in x]
    feature_sets["dxy"] = [x for x in sumstats if "dxy" in x and "_cor" not in x and "_cor" not in x]
    feature_sets["trait"] = [x for x in sumstats if "trait" in x and "_cor" not in x]

    feature_sets["abundance+pi"] =  feature_sets["abundance"] + feature_sets["pi"] + ["abundance_pi_cor"]
    feature_sets["pi+dxy"] = feature_sets["pi"] + feature_sets["dxy"] + ["dxy_pi_cor"]
    feature_sets["pi+trait"] = feature_sets["pi"] + feature_sets["trait"] + ["pi_trait_cor"]
    feature_sets["abundance+trait"] = feature_sets["abundance"] + feature_sets["trait"] + ["abundance_trait_cor"]

    feature_sets["abundance+pi+trait"] = feature_sets["abundance+pi"] + feature_sets["trait"] + ["abundance_trait_cor"] + ["pi_trait_cor"]

    feature_sets["all"] = sumstats

    missing_axes = set(MESS.stats.data_axes).difference(set(empirical_axes))
    if missing_axes:
        feature_sets.pop("all", '')
        for k in feature_sets:
            for m in missing_axes:
                if m in k:
                    feature_sets.pop(k, '')

    return feature_sets


def get_feature_set(data_axes=''):
    """
    TODO: This doesn't exactly work how I want yet.
    """
    if not data_axes:
        data_axes = MESS.stats.data_axes
    return feature_sets()[data_axes]


if __name__ == "__main__":
    ## Test SAD()
    dat = [1,2,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]
    dat.extend([9] * 100)
    sad = SAD(dat)
    sad_oct = SAD(dat, octaves=True)
    #assert(cmp(sad.values(), [2, 2, 1, 1, 1, 1, 1]))
    #assert(cmp(sad_oct.values(), [2, 2, 1, 1, 2, 0, 1])) 
    print(("SAD - {}".format(sad)))
    print(("SAD octaves - {}".format(sad_oct)))

    pis = MESS.rng.rng.random(10)/10
    dxys = MESS.rng.rng.random(10)/10
    abunds = MESS.rng.rng.integers(1, 100, 10)
    trts = MESS.rng.rng.random(10)*10

    dat = pd.DataFrame([], columns=["pis", "dxys", "abunds", "traits"])
    dat["pis"] = pis
    dat["dxys"] = dxys
    dat["abunds"] = abunds
    dat["traits"] = trts

    ss = calculate_sumstats(dat, sgd_bins=10, sgd_dims=1)
    print(("sumstats", len(ss.values), ss))
    header = _get_sumstats_header(sgd_bins=10, sgd_dims=1)
    print(("header", len(header), header))
