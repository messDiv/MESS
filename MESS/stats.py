
import collections
import numpy as np
import math

## Here abundances is an ordered dict of tuples which are (abundance, count)
## This is the typical format returned by LocalCommunity.get_abundances(octaves=False)
def shannon(abundances):
    ## Unpack the abundance dist
    abunds = [v for v in abundances.values()]
    tot = np.sum(abunds)
    return -1 * np.sum([x/float(tot) * math.log(x/float(tot)) for x in abunds  if x > 0])


def SAD(list_of_inds, octaves=False):
    ## Make a counter for the local_community, counts the number of
    ## individuals w/in each species
    abundances = collections.Counter(list_of_inds)

    ## If we were doing mode=volcanic then there may be some remaining
    ## space in our carrying capacity that is unoccupied (indicated by
    ## None in the abundances.keys()
    try:
        abundances.pop(None)
    except KeyError:
        pass

    ## Now for each abundance class you have to go through and
    ## count the number of species at that abundance.
    abundance_distribution = collections.Counter(abundances.values())
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
            for i, j in abundance_distribution.items():
                if (i < maxval) and (i >= minval):
                    count += j
            dist_in_octaves[minval] = count
            minval = minval * 2
            maxval = maxval * 2
        abundance_distribution = dist_in_octaves
    return abundance_distribution

if __name__ == "__main__":
    ## Test SAD()
    dat = [1,2,3,3,4,4,4,4,4,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,7,7,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8]
    dat.extend([9] * 100)
    sad = SAD([1,2,3,4,4,4,4,5,5,5,5,5,5])
    sad_oct = SAD([1,2,3,4,4,4,4,5,5,5,5,5,5], octaves=True)
    assert(cmp(sad.values(), [2, 2, 1, 1, 1, 1, 1]))
    assert(cmp(sad_oct.values(), [2, 2, 1, 1, 2, 0, 1])) 
    print("SAD - {}".format(sad))
    print("SAD octaves - {}".format(sad_oct))
