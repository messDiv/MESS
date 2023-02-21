"""
Massive Eco-Evolutionary Synthesis Simulations
"""

__version__ = "0.1.2"
__author__ = "Isaac Overcast, Megan Ruffley"

## Possible values for __loglevel__: "DEBUG"  "INFO"  "WARN"  "ERROR"
__debugflag__ = "./.debug"
__debugfile__ = "./mess_log.txt"

## define state vars
__interactive__ = 1      ## CLI __main__ changes to 0

# pylint: disable=C0103
import os as _os
import atexit as _atexit

## Force matplotlib to behave on headless environments
## If running inside a jupyter notebook and you call %matplotlib inline
## then this call to matplotlib.use will raise a warning, which is safe
## to ignore.
import matplotlib
import warnings
with warnings.catch_warnings(record=True) as w:
    matplotlib.use("agg")

from .Region import Region
from .LocalCommunity import LocalCommunity
from .Metacommunity import Metacommunity
from . import util
from . import plotting
from . import stats

## sklearn throws a DeprecationWarning: sklearn.externals.joblib
## Will be fixed in 0.23: https://github.com/deephyper/deephyper/issues/12
## Very annoying if you don't catch and ignore this.
with warnings.catch_warnings(record=True) as w:
    from . import inference

####################################################################
## create logger for debugging
## this needs to come after __loglevel__ definition
## sets log config and prints warning if __loglevel__ is in hackers mode
import logging as _logging
import logging.config as _lconfig

## debug is set based on whether the flag exists
if _os.path.exists(__debugflag__):
    __loglevel__ = "DEBUG"
else:
    __loglevel__ = "ERROR"

## check that all dependencies exist and are working
import subprocess as _subprocess
import sys as _sys
import socket as _socket


_LOGGER = _logging.getLogger(__name__)
if __loglevel__ == "DEBUG":
    _LOGGER.debug("Engine init")

## Catch warnings and write them to the log file
#def warn(*args, **kwargs):
#    for arg in args:
#        _LOGGER.warn(arg)
#warnings.warn = warn

def cluster_info(ipyclient, spacer=""):
    """ reports host and engine info for an ipyclient """    
    ## get engine data, skips busy engines.    
    hosts = []
    for eid in ipyclient.ids:
        engine = ipyclient[eid]
        if not engine.outstanding:
            hosts.append(engine.apply(_socket.gethostname))

    ## report it
    hosts = [i.get() for i in hosts]
    result = []
    for hostname in set(hosts):
        result.append("{}host compute node: [{} cores] on {}"\
            .format(spacer, hosts.count(hostname), hostname))
    print(("\n".join(result)))



def _debug_on():
    """
    Turns on debugging by creating hidden tmp file
    This is only run by the __main__ engine.
    """
    ## make tmp file and set loglevel for top-level init
    with open(__debugflag__, 'w') as dfile:
        dfile.write("wat")
    __loglevel__ = "DEBUG"
    _LOGGER.info("debugging turned on and registered to be turned off at exit")
    _set_debug_dict(__loglevel__)


def _set_debug_dict(__loglevel__):
    """ set the debug dict """

    _lconfig.dictConfig({
    'version': 1,
    'disable_existing_loggers': False,

    'formatters': {
        'standard': {
            'format': "%(asctime)s \t"\
                     +"pid=%(process)d \t"\
                     +"[%(filename)s]\t"\
                     +"%(levelname)s \t"\
                     +"%(message)s"
        },
    },
    'handlers': {
        __name__: {
            'level':__loglevel__,
            'class':'logging.FileHandler',
            'filename':__debugfile__,
            'formatter':"standard",
            'mode':'a+'
        }
    },
    'loggers':{
        __name__: {
            'handlers': [__name__],
            'level': __loglevel__,
            'propogate': True
        }
    }
    })

_set_debug_dict(__loglevel__)


def _debug_off():
    """ turns off debugging by removing hidden tmp file """
    if _os.path.exists(__debugflag__):
        _os.remove(__debugflag__)
    __loglevel__ = "ERROR"
    _LOGGER.info("debugging turned off")
    _set_debug_dict(__loglevel__)
  

if __name__ == "__main__":
    pass
