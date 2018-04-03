from __future__ import print_function
import subprocess
import shlex
import sys
import os

import logging
LOGGER = logging.getLogger(__name__)

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
def qsort(arr):
     if len(arr) <= 1:
          return arr
     else:
          return qsort([x for x in arr[1:] if x.abundance<arr[0].abundance])\
                    + [arr[0]] + qsort([x for x in arr[1:] if x.abundance>=arr[0].abundance])


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


def set_params(data, param, newvalue):
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
    allowed_params = data.paramsdict.keys()
    ## require parameter recognition
    if not param in allowed_params:
        raise MESSError("Parameter key not recognized: {}"\
                                .format(param))
    try:
        data._paramschecker(param, newvalue)
    except Exception as inst:
        raise MESSError(BAD_PARAMETER.format(param, inst, newvalue))
    return data


#############################
## Model Parameter Info Dicts
#############################
REGION_PARAMS = {
    "simulation_name" : "The name of this simulation scenario",\
    "project_dir" : "Where to save files",\
    "metacommunity_type" : "Specify structure of the metacommunity. Options: uniform/logser/<filename>",\
    "data_model" : "Structure of data output to reference table (see docs)",\
    "recording_period" : "Number of forward-time generations for gathering data."\
}

LOCAL_PARAMS = {
    "name" : "Local community name",\
    "mode" : "Local community formation mode (volcanic/landbridge)",\
    "K" : "Local carrying capacity",\
    "colrate" : "Colonization rate into local community",\
    "mig_clust_size" : "# of individuals per colonization event",\
    "age" : "Local community age"
}

## Error messages
BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """
if __name__ == "__main__":
    print(REGION_PARAMS)
