
"""
.. module:: rng
   :synopsis: Random generator -- enables repeatability of the experiments

"""

import numpy as np
import logging
LOGGER = logging.getLogger(__name__)

def init(sd=None):
    if sd == None or sd=="*":
        sd = str(np.random.randint(2147483647))
        # R only support 32bits integers
    LOGGER.debug("Intializing random number generator with seed - {}".format(sd))

    global rng
    rng = np.random.default_rng(int(sd))

    global seed
    seed = sd

@property
def rng():
    return rng
