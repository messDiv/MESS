#!/usr/bin/env python2.7

import numpy as np
import pandas as pd
import collections
import itertools
import math
import os

import logging
LOGGER = logging.getLogger(__name__)

class SGD(object):

    def __init__(self, pis, dxys=[], ndims=2, nbins=10):
        """ Construct the Species Genetic Diversity histogram. Input is
        a list of pi values and optionally a list of dxy values if 2D
        is desired. nbins is the dimension of one axis of the SGD, and
        the result will be square.

        ## TODO: Possibly experiment with normed/density params
        """

        self.ndims = ndims
        self.nbins = nbins

        flatten = False
        if ndims == 1:
            flatten = True

        if not len(pis):
            shape = (nbins, nbins) if ndims == 2 else (nbins,)
            hist = np.zeros(shape)
        elif not len(dxys):
            hist, xedges = np.histogram(pis, bins=nbins)
        else:
            hist, xedges, yedges = np.histogram2d(pis, dxys, bins=nbins)
            ## numpy 2d histogram is oriented somewhat unintuitively
            ## transposing the array turns y axis into rows
            hist = np.sum(hist.T, axis=0) if flatten else hist.T
        ## Return as int array, since histogram2d uses float64 bins
        self.sgd = hist.astype(int)


    @classmethod
    def from_df(self, df):
        pass   

    def __str__(self):
        return np.array2string(self.sgd).replace('\n', '').replace("[", '').replace("]", "")


    def __repr__(self):
        return self.__str__()


    def to_series(self):
        if self.ndims == 1:
            self._sgd_header_template = "SGD_{}"
            self._sgd_header = [self._sgd_header_template.format(x) for x in range(self.nbins)]
            dat = self.sgd
        else:
            self._sgd_header_template = "SGD_{}_{}"
            self._sgd_header = [self._sgd_header_template.format(x,y) for x in range(self.nbins)\
                                                                      for y in range(self.nbins)]
            dat = self.sgd.ravel()

        return pd.Series(dat, index=self._sgd_header)


    def to_dict(self):
        return collections.OrderedDict(sorted(self.to_series().to_dict().items()))


if __name__ == "__main__":
    ## Test SAD()
    import numpy as np

    pis = MESS.rng.rng.exponential(0.05, size=100)
    dxys = MESS.rng.rng.exponential(0.05, size=100)
    sgd = SGD(pis)
    print(sgd)
    sgd2 = SGD(pis, dxys, nbins=5)
    print(sgd2)

    print((sgd2.to_series()))

    print((SGD([])))
    print((SGD([], ndims=1)))
    print((sgd2.to_dict()))
