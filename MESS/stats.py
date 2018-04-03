
import numpy as np
import math

## Here abundances is an ordered dict of tuples which are (abundance, count)
## This is the typical format returned by LocalCommunity.get_abundances(octaves=False)
def shannon(abundances):
    ## Unpack the abundance dist
    abunds = [v for v in abundances.values()]
    tot = np.sum(abunds)
    return -1 * np.sum([x/float(tot) * math.log(x/float(tot)) for x in abunds  if x > 0])
