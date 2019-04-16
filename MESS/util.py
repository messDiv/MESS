from __future__ import print_function
import collections
import functools
import subprocess
import shlex
import glob
import sys
import os

import logging
LOGGER = logging.getLogger(__name__)

from .stats import *
from .SGD import SGD

## Custom exception class
class MESSError(Exception):
    """ General MESS exception """
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


def detect_cpus():
    """
    Detects the number of CPUs on a system. This is better than asking
    ipyparallel since ipp has to wait for Engines to spin up.
    """
    # Linux, Unix and MacOS:
    if hasattr(os, "sysconf"):
        if os.sysconf_names.has_key("SC_NPROCESSORS_ONLN"):
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if os.environ.has_key("NUMBER_OF_PROCESSORS"):
        ncpus = int(os.environ["NUMBER_OF_PROCESSORS"])
        if ncpus > 0:
            return ncpus
    return 1 # Default


## quicksort stolen from the internet
## This sorts species by Ne
def qsort(arr):
     if len(arr) <= 1:
          return arr
     else:
          return qsort([x for x in arr[1:] if x.stats["Ne_local"]<arr[0].stats["Ne_local"]])\
                 + [arr[0]] +\
                 qsort([x for x in arr[1:] if x.stats["Ne_local"]>=arr[0].stats["Ne_local"]])


def progressbar(nsims, finished, msg=""):
    """ prints a progress bar """
    progress = 100*(finished / float(nsims))
    hashes = '#'*int(progress/5.)
    nohash = ' '*int(20-len(hashes))
    print("\r  [{}] {:>3}% {} ".format(hashes+nohash, int(progress), msg), end="")
    sys.stdout.flush()


def _expander(namepath):
    """ expand ./ ~ and ../ designators in location names """
    if "~" in namepath:
        namepath = os.path.expanduser(namepath)
    else:
        namepath = os.path.abspath(namepath)
    return namepath


def sample_param_range(param, nsamps=1):
    """ Sample a parameter from a range. This is used to allow parameter
    values in the params file to be specified as a tuple and have simulations
    sample from the range on the fly.
    """
    LOGGER.debug("Sampled from range {}".format(param))
    if isinstance(param, tuple):
        if isinstance(param[0], float):
            param = np.random.uniform(param[0], param[1], nsamps)
        else:
            param = np.random.randint(param[0], param[1], nsamps)
    elif param == 0:
        param = np.round(np.random.random(nsamps), 3)
    else:
        param = [param] * nsamps
    LOGGER.debug("Sampled Value: {}".format(param))
    return param


def tuplecheck(newvalue, dtype=str):
    """
    Takes a string argument and returns value as a tuple.
    Needed for paramfile conversion from CLI to set_params args
    """
    ## TODO: This actually should work

    ## If it's a list then this is probably api mode so the types
    ## of the values should be fine.
    if isinstance(newvalue, tuple):
        ## Already a tuple so we good to go
        pass
    elif isinstance(newvalue, list):
        try:
            newvalue = tuple(newvalue)
        except TypeError:
            raise MESSError("tuplecheck failed for {}, improper list format".format(newvalue))
    else:
        try:
            ## If it's just one value of the proper dtype this should
            ## suffice to catch it.
            newvalue = dtype(newvalue)
        except Exception as inst:
            ## Failed to cast to dtype, so this is probably a prior range
            ## Using comma as seperator here to avoid ugliness of using -,
            ## because we may want negative values in some priors.
            try:
                newvalue = newvalue.rstrip(")").strip("(")
                minval = dtype(newvalue.split("-")[0].strip())
                maxval = dtype(newvalue.split("-")[1].strip())
                newvalue = tuple([minval, maxval])
            except Exception as inst:
                raise MESSError("{}\ttuplecheck() failed to cast to {} - {}"\
                            .format(inst, dtype, newvalue))

    LOGGER.debug("Returning tuple - {}".format(newvalue))
    return newvalue


def import_empirical(input_dir):
    """ Validate, format, and import empirical data from a directory."""
    if not os.path.exists(input_dir):
        raise MESSError("Importing empirical directory that doesn't exist - {}".format(input_dir))

    abundfile = os.path.join(input_dir, "abunds.txt")
    if os.path.isfile(abundfile):
        dat = open(abundfile).read().strip().split(",")
        dat = map(int, dat)
        print("Got empirical abundance Hill_1 - {}".format(hill_number(SAD(dat, from_abundances=True, raw_abunds=True))))

    fastadir = os.path.join(input_dir, "spider-fasta")
    if os.path.exists(fastadir):
        files = glob.glob(os.path.join(fastadir, "*.fasta"))
        pis = []
        for f in files:
            dat = read_sequence(f)
            pis.append(pi(dat))
            print(pi(dat), f)
        dat = SGD(pis)
        print("Empirical SGD - {}".format(dat))


def read_sequence(fastafile):
    """ Read in sequence data from a file and return just the sequencess
    as an array. """
    with open(fastafile) as infile:
        datlines = infile.readlines()[1::2]
        datlines = [x.strip() for x in datlines]
        return(np.array(datlines))


def set_params(data, param, newvalue, quiet=False):
    """
    Set a parameter to a new value. Raises error if newvalue is wrong type.
    This is used to set parameters on both the Region and LocalCommunity
    paramsdicts.

    Parameters
    ----------
    param : str
        The string name (e.g., "project_dir") for the parameter 
        that will be changed.
    
    newvalue : int, str, or tuple
        The new value for the parameter selected for `param`.
        If the wrong type is entered for newvalue
        (e.g., a str when it should be an int), an error will be raised.
        Further information about each parameter is also available
        in the documentation.
    """
    LOGGER.debug("set param: {} {} = {}".format(data, param, newvalue))
    allowed_params = data.paramsdict.keys()
    ## require parameter recognition
    if not param in allowed_params:
        raise MESSError("Parameter key not recognized: {}"\
                                .format(param))
    try:
        data._paramschecker(param, newvalue, quiet)
    except Exception as inst:
        raise MESSError(BAD_PARAMETER.format(param, inst, newvalue))
    return data


class memoize(object):
   '''Decorator. Caches a function's return value each time it is called.
   If called later with the same arguments, the cached value is returned
   (not reevaluated).

   Class poached from https://wiki.python.org/moin/PythonDecoratorLibrary#Memoize
   '''
   def __init__(self, func):
      self.func = func
      self.cache = {}
   def __call__(self, *args):
      if not isinstance(args, collections.Hashable):
         # uncacheable. a list, for instance.
         # better to not cache than blow up.
         return self.func(*args)
      if args in self.cache:
         return self.cache[args]
      else:
         value = self.func(*args)
         self.cache[args] = value
         return value
   def __repr__(self):
      '''Return the function's docstring.'''
      return self.func.__doc__
   def __get__(self, obj, objtype):
      '''Support instance methods.'''
      return functools.partial(self.__call__, obj)


## Error messages
BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """

if __name__ == "__main__":
    import_empirical("empirical_data")
