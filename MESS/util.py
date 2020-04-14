
import collections
import glob
import functools
import itertools
import logging
import numpy as np
import os
import pandas as pd
import random
import shlex
import subprocess
import sys
import toytree

from .stats import hill_number
from .SGD import SGD

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
        if "SC_NPROCESSORS_ONLN" in os.sysconf_names:
            # Linux & Unix:
            ncpus = os.sysconf("SC_NPROCESSORS_ONLN")
            if isinstance(ncpus, int) and ncpus > 0:
                return ncpus
        else: # OSX:
            return int(os.popen2("sysctl -n hw.ncpu")[1].read())
    # Windows:
    if "NUMBER_OF_PROCESSORS" in os.environ:
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


def sample_param_range(param, nsamps=1, loguniform=False):
    """ Sample a parameter from a range. This is used to allow parameter
    values in the params file to be specified as a tuple and have simulations
    sample from the range on the fly.
    """
    LOGGER.debug("Sampled from range {}".format(param))
    if isinstance(param, tuple):
        if isinstance(param[0], float):
            if loguniform:
                param = np.random.uniform(np.log10(param[0]), np.log10(param[1]), size=nsamps)
                param = np.power(10, param)
            else:
                param = np.round(np.random.uniform(param[0], param[1], nsamps), 5)
        else:
            if loguniform:
                param = np.random.uniform(np.log10(param[0]), np.log10(param[1]), nsamps)
                param = np.int32(np.power(10, param))
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

    ## If it's a list then this is probably api mode so the types
    ## of the values should be fine.
    if isinstance(newvalue, tuple):
        ## Already a tuple so we good to go, just make sure the
        ## values are the right dtype
        newvalue = (dtype(newvalue[0]), dtype(newvalue[1]))
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
            ## In some cases you may get string representations of float
            ## values that _should_ be int. Here we first cast the string
            ## to a float, and then cast to dtype.
            ## https://stackoverflow.com/questions/1841565/valueerror-invalid-literal-for-int-with-base-10
            try:
                newvalue = newvalue.rstrip(")").strip("(")
                minval = dtype(float(newvalue.split("-")[0].strip()))
                maxval = dtype(float(newvalue.split("-")[1].strip()))
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
        dat = list(map(int, dat))
        print("Got empirical abundance Hill_1 - {}".format(hill_number(SAD(dat, from_abundances=True, raw_abunds=True), order=1)))

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


def set_params(data, param, newvalue, quiet=True):
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
    #allowed_params = list(data.paramsdict.keys())
    ## require parameter recognition
    if not param in list(data.paramsdict.keys()):
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


def synthetic_community(model="random", nspecies=10, nspecies_meta=None):
    """
    Generate artificial, toy data in the format required for the
    MESS.stats.calculate_sumstats function. Useful for testing.

    .. note::
        The competition and filtering models *should* generate trait
        distributions that correlate in expected ways with abundance,
        but they don't.

    :param str model: One of: `random`, `correlated`, `filtering` or
        `competition`. Random generates random uniform values for
        all data axes and scales them to look plausible. Max abundance is
        fixed at 1000.
    :param int nspecies: The number of species to return.
    :param int nspecies_meta: Number of metacommunity species to return

    :return: A pandas.DataFrame with 4 columns (abundance, pi, dxy, trait)
        and `nspecies` rows, populated with toy data.
    """
    abunds = np.random.randint(1, 1000, nspecies)

    if model == "random":
        pis = np.random.random(nspecies)/10
        dxys = np.random.random(nspecies)/10
        trts = np.random.random(nspecies)*10
    elif model == "correlated":
        pis = abunds/1000.
        dxys = pis*1.15
        trts = abunds/10.
    elif model == "filtering":
        def filt(strength=10, victim_trait=-3.32923208, filt_opt=0):
            return 1 - (np.exp(-((victim_trait - filt_opt) ** 2)/strength))
        trts = 1-np.array([filt(strength=10, victim_trait=x/100, filt_opt=abunds.max()/100) for x in abunds])
        pis = np.random.random(nspecies)/10
        dxys = np.random.random(nspecies)/10
    elif model == "competition":
        trts = np.array([x - abunds.mean() for x in abunds])/10
        pis = np.random.random(nspecies)/10
        dxys = np.random.random(nspecies)/10
    dat = pd.DataFrame([], columns=["pi", "dxy", "abundance", "trait"],)
    dat["pi"] = pis
    dat["dxy"] = dxys
    dat["abundance"] = abunds
    dat["trait"] = trts

    ## Get the random trees. The simulated local species will be
    ## named l0, l1, l2 ('l' for local). In this way the local
    ## trees can be pasted on to the metacommunity tree without
    ## name collisions (toytree random trees are generated with
    ## species names like r0, r1, r2, r3. Greatly simplifies the
    ## problem of pruning and grafting the local trees.
    for clade in _random_ksubset(dat.index.to_list(),
                                np.random.randint(low=1, high=nspecies)):
        ids = list(clade)
        nw = _random_newick(len(clade))
        nw = _update_newick_names(nw, ["l"+str(x) for x in ids])
        for idx in ids:
            dat.loc[idx, "tree"] = nw

    ## Give them names that will be like the names for the randomly
    ## simulated trees. Plus I think 'ape' will fall over with integer
    ## species ids.
    dat.index = ["l"+str(x) for x in range(0, nspecies)]

    if nspecies_meta:
        mtree = random_metatree(dat["tree"], ntips=nspecies_meta)
        return dat, mtree
    else:
        return dat


## Helper functions for building the reandom trees
def _random_ksubset(ls, k):
    ls = list(ls)
    indices = list(range(k))
    indices.extend([random.choice(list(range(k))) for _ in range(len(ls) - k)])
    random.shuffle(indices)
    return [{x[1] for x in xs} for (_, xs) in
            itertools.groupby(sorted(zip(indices, ls)), lambda x: x[0])]


def _random_newick(ntips, treeheight=1, treefunc=toytree.rtree.unittree):
    if ntips < 2:
        return "(r0:1);"
    else:
        return treefunc(ntips=ntips).write(tree_format=5)


def _update_newick_names(nw, names):
    tre = toytree.tree(nw)
    leaves = tre.get_tip_labels()
    nodes = [tre.treenode.search_nodes(name=x)[0] for x in leaves]
    for node,name in zip(nodes, names):
        node.name = name
    return tre.write(tree_format=5)


def random_metatree(local_trees, ntips=20, treefunc=toytree.rtree.unittree):
    metatree = _random_newick(ntips=ntips, treeheight=10, treefunc=treefunc)

    tre = toytree.tree(metatree)
    tips = tre.get_tip_labels()

    local_trees = list(set(local_trees))

    drop_tips = np.random.choice(tips, len(local_trees), replace=False)
    for tip, ltree in zip (drop_tips, local_trees):
        drop_node = tre.treenode.get_leaves_by_name(name=tip)[0]
        sis = drop_node.get_sisters()[0]
        ttree = toytree.tree(ltree)
        _root = ttree.treenode.get_tree_root()
        height = _root.get_distance(_root.get_farthest_leaf()[0])
        ttree.treenode.dist = drop_node.dist - height
        new_node = drop_node.add_sister(sister=ttree.treenode)

        _ = sis.remove_sister(sister=drop_node)

    return tre.write(tree_format=5)


## Error messages
BAD_PARAMETER = """\
    Error setting parameter '{}'
    {}
    You entered: {}
    """

if __name__ == "__main__":
    import_empirical("empirical_data")
