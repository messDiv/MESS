from util import *
from collections import OrderedDict
import string
import os

import MESS

import logging
LOGGER = logging.getLogger(__name__)

class Region(object):
    """ A Massive eco-evolutionary synthesis simulation object.

    """

    def __init__(self, name, quiet=False, **kwargs):
        if not name:
            raise MESSError(REQUIRE_NAME)

        ## Do some checking here to make sure the name doesn't have
        ## special characters, spaces, or path delimiters. Allow _ and -.
        ## This will raise an error immediately if there are bad chars in name.
        self._check_name(name)
        self.name = name

        self._version = MESS.__version__

        ## stores default ipcluster launch info
        self._ipcluster = {
            "cluster_id" : "",
            "profile" : "default",
            "engines" : "Local",
            "quiet" : 0,
            "timeout" : 120,
            "cores" : 0, #detect_cpus(),
            "threads" : 2,
            "pids": {},
            }

        ## the default params dict
        self.paramsdict = OrderedDict([
                       ("simulation_name", name),
                       ("project_dir", "./default_MESS"),
                       ("trim_loci", (0, 0, 0, 0)),
                       ("output_formats", ['p', 's', 'v']),
                       ("pop_assign_file", ""),
        ])

        self.islands = {}


    def __str__(self):
        return "<MESS.Region object {}: {}>".format(self.name, self.islands.keys())

    ## Test assembly name is valid and raise if it contains any special characters
    def _check_name(self, name):
        invalid_chars = string.punctuation.replace("_", "")\
                                          .replace("-", "")+ " "
        if any(char in invalid_chars for char in name):
            raise MESSError(BAD_MESS_NAME.format(name))


    def _link_local(self, local_community):
        """ Just link a local community that was created externally.
        This is primarily used by __main__ during the initialization process."""
        self.islands[local_community.paramsdict["name"]] = local_community


    def _paramschecker(self, param, newvalue):
        """ Raises exceptions when params are set to values they should not be"""
        ## TODO: This should actually check the values and make sure they make sense
        try:
            LOGGER.debug("set param {} - {} = {}".format(self, param, newvalue))
            self.paramsdict[param] = newvalue
        except Exception as inst:
            ## Do something intelligent here?
            raise

    def write_params(self, outfile=None, force=False):
        """ Write out the parameters of this model to a file properly
        formatted as input for `MESS -p <params.txt>`. A good and
        simple way to share/archive parameter settings for simulations.
        This is also the function that's used by __main__ to
        generate default params.txt files for `MESS -n`
        """
        if outfile is None:
            outfile = "params-"+self.name+".txt"

        ## Test if params file already exists?
        ## If not forcing, test for file and bail out if it exists
        if not force:
            if os.path.isfile(outfile):
                raise MESSError(PARAMS_EXISTS.format(outfile))

        with open(outfile, 'w') as paramsfile:
            ## Write the header. Format to 80 columns
            header = "------- MESS params file (v.{})".format(MESS.__version__)
            header += ("-"*(80-len(header)))
            paramsfile.write(header)

            ## Whip through the current paramsdict and write out the current
            ## param value, the ordered dict index number. Also,
            ## get the short description from paramsinfo. Make it look pretty,
            ## pad nicely if at all possible.
            for key, val in self.paramsdict.iteritems():
                ## If multiple elements, write them out comma separated
                if isinstance(val, list) or isinstance(val, tuple):
                    paramvalue = ", ".join([str(i) for i in val])
                else:
                    paramvalue = str(val)

                padding = (" "*(30-len(paramvalue)))
                paramkey = self.paramsdict.keys().index(key)
                paramindex = " ## [{}] ".format(paramkey)
                LOGGER.debug("{} {} {}".format(key, val, paramindex))
                #name = "[{}]: ".format(paramname(paramkey))
                name = "[{}]: ".format(key)
                #description = paraminfo(paramkey, short=True)
                description = "wat"
                paramsfile.write("\n" + paramvalue + padding + \
                                        paramindex + name + description)

            paramsfile.write("\n")

        for island in self.islands.values():
            island.write_params(outfile)


    def add_local_community(self, name, K, c, quiet):
        loc = MESS.LocalCommunity(name, K, c, quiet)
        ## TODO: Ensure island names aren't dupes
        self.islands[name] = loc


BAD_MESS_NAME = """\
    No spaces or special characters of any kind are allowed in the simulation 
    name. Special characters include all punctuation except dash '-' and 
    underscore '_'. A good practice is to replace spaces with underscores '_'.
    An example of a good simulation name is: hawaiian_arthropods 
    
    Here's what you put:
    {}
    """

PARAMS_EXISTS = """
    Error: Params file already exists: {}
    Use force argument to overwrite.
    """

REQUIRE_NAME = """\
    Simulation scenario name _must_ be set. This is the first parameter in the
    params.txt file, and will be used as a prefix for output files. It should be a 
    short string with no special characters, i.e., not a path (no \"/\" characters).
    If you need a suggestion, name it after the location you're working on.
    """

