from __future__ import print_function

import collections
import numpy as np
import math

## Here abundances is an ordered dict of tuples which are (abundance, count)
## This is the typical format returned by LocalCommunity.get_abundances(octaves=False)
def shannon(abundances):
    ## Unpack the abundance dist
    abunds = [v for v in list(abundances.values())]
    tot = np.sum(abunds)
    return -1 * np.sum([x/float(tot) * math.log(x/float(tot)) for x in abunds  if x > 0])


def SAD(community, from_abundances=False, octaves=False):
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
    abundances = collections.Counter(community)

    ## If we were doing mode=volcanic then there may be some remaining
    ## space in our carrying capacity that is unoccupied (indicated by
    ## None in the abundances.keys()
    try:
        abundances.pop(None)
    except KeyError:
        pass

    ## Now for each abundance class you have to go through and
    ## count the number of species at that abundance.
    abundance_distribution = collections.Counter(list(abundances.values()))
    abundance_distribution = collections.OrderedDict(sorted(abundance_distribution.items()))
    if octaves:
        dist_in_octaves = collections.OrderedDict()
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
