from util import *
from collections import OrderedDict
import numpy as np
import time
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
        ## If you add a parameter to this dictionary you need
        ## to also add a short description to the REGION_PARAMS dict
        ## at the end of this file
        ##
        ## Also be sure to add it to _paramschecker so the type gets set correctly
        self.paramsdict = OrderedDict([
                       ("simulation_name", name),
                       ("project_dir", "./default_MESS"),
                       ("metacommunity_type", "logser"),
                       ("data_model", 4),
                       ("generations", 0),
                       ("recording_period", 10000),
                       ("allow_multiple_colonizations", True),
        ])

        ## Track local communities in this model and colonization rates among them
        self.islands = {}
        self.colonization_matrix = []


    #########################
    ## Housekeeping functions
    #########################
    def __str__(self):
        return "<MESS.Region {}: {}>".format(self.name, self.islands.keys())

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


    def _paramschecker(self, param, newvalue, quiet=False):
        """ Raises exceptions when params are set to values they should not be"""
        ## TODO: This should actually check the values and make sure they make sense
        try:
            LOGGER.debug("set param {} - {} = {}".format(self, param, newvalue))
            ## Cast params to correct types
            if param in ["data_model", "recording_period"]:
                self.paramsdict[param] = int(float(newvalue))

            elif param == "project_dir":
                ## If it already exists then just inform the user that we'll be adding
                ## more simulations to the current project directory
                self.paramsdict[param] = newvalue
                if not os.path.exists(self.paramsdict["project_dir"]):
                    os.mkdir(self.paramsdict["project_dir"])
                else:
                    if not quiet:
                        print("    Project directory exists. Additional simulations will be appended.")

            elif param == "generations":
                if "-" in newvalue:
                    low = int(float(newvalue.split("-")[0].strip()))
                    high = int(float(newvalue.split("-")[1].strip()))
                    if low >= high:
                        raise MESSError("Bad parameter for generations: low must be < high value.")
                    self.paramsdict[param] = tuple([low, high])
                    print(self.paramsdict[param])
                else:
                    self.paramsdict[param] = int(float(newvalue))

            elif param == "metacommunity_type":
                ## TODO: Check that the types are ok
                self.paramsdict[param] = newvalue

            elif param == "allow_multiple_colonizations":
                self.paramsdict[param] = newvalue.lower() in ["true"]
                
            else:
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

                padding = (" "*(20-len(paramvalue)))
                paramkey = self.paramsdict.keys().index(key)
                paramindex = " ## [{}] ".format(paramkey)
                LOGGER.debug("{} {} {}".format(key, val, paramindex))
                #name = "[{}]: ".format(paramname(paramkey))
                name = "[{}]: ".format(key)
                #description = paraminfo(paramkey, short=True)
                description = REGION_PARAMS[key]
                paramsfile.write("\n" + paramvalue + padding + \
                                        paramindex + name + description)

            paramsfile.write("\n")

        for island in self.islands.values():
            island.write_params(outfile)


    ########################
    ## Model functions/API
    ########################
    def add_local_community(self, name, K, c, quiet=False):
        loc = MESS.LocalCommunity(name, K, c, quiet)
        ## TODO: Ensure island names aren't dupes
        self.islands[name] = loc


    def update_colonization_matrix(self, matrix):
        """ Set the matrix that describes colonization rate between local communities."""
        ## TODO: Make sure the matrix is the right shape
        self.colonization_matrix = matrix


    ## Main function for managing cluster parallelized simulations
    def run(self, sims, force=False, ipyclient=None):
        """ Do the heavy lifting here"""
        print("    Generating {} simulations.".format(sims))
        if not ipyclient:
            ## Run serially
            for i in xrange(sims):
                progressbar(sims, i)
                gens = self.paramsdict["generations"]
                if type(gens) is tuple:
                    ## Generations is a range so uniform random sample for this sim
                    gens = np.random.randint(gens[0], gens[1])
                    self.simulate(nsteps=gens)
                else:
                    ## If doing lambda we only really care about 3ish signficant figures
                    _lambda = round(np.random.random(), 3)
                    self.simulate(_lambda=_lambda)
            progressbar(100, 100, "Finished {} simulations\n".format(sims))


    ## This is the function that will be run inside the cluster engine, so
    ## everything it does needs to be
    ##
    ## The _lambda getting an underscore here is an artifact of the choice
    ## of symbol for %equilibrium colliding with the python anonymous function
    ## name. Don't confuse them.
    ##
    ## TODO: Need to think about how to quantify lambda if there is more
    ## than one island
    def simulate(self, _lambda=0, nsteps=0):
        LOGGER.debug("Entering simulate(): lambda={}, nsteps={}".format(_lambda, nsteps))
        if _lambda > 0 and nsteps > 0:
            LOGGER.error("simulate accepts only one of either lambda or nsteps args")
            return
        ## Make output directory formatted like <output_dir>/<name>-<timestamp><random 2 digit #>
        ## This will _mostly_ avoid output directory collisions, right?
        outdir = os.path.join(self.paramsdict["project_dir"],\
                              self.paramsdict["simulation_name"]\
                              + "-" + str(time.time()).replace(".", "")[-7:]\
                              + str(np.random.randint(100)))
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        step = 0
        ## Create an temp function to test whether we've reached the end of this simulation
        if _lambda > 0:
            ## In the absence of a better strategy just test for lambda in the first local community
            done_func = lambda: self.islands.values()[0]._lambda() >= _lambda
        else:
            done_func = lambda: step >= nsteps

        while not done_func():
            for island in self.islands.values():
                island.step()
            step += 1


#############################
## Model Parameter Info Dicts
#############################
REGION_PARAMS = {
    "simulation_name" : "The name of this simulation scenario",\
    "project_dir" : "Where to save files",\
    "metacommunity_type" : "Options: uniform/logser/<filename>",\
    "data_model" : "Structure of data output to reference table (see docs)",\
    "generations" : "Duration of simulations. Specify int range or 0 for lambda.",\
    "recording_period" : "Number of forward-time generations between samples for logging",\
    "allow_multiple_colonizations" : "Toggle allowing post colonization migration",\
}


#############################
## Global error messages
#############################
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

if __name__ == "__main__":
    logging.info("wat")
    data = Region("tmp")
    data.add_local_community("tmploc", 1000, 0.01)
    #print("Testing step function.")
    #data.simulate(_lambda=0, nsteps=1000)
    #print(data.islands.values()[0])
    #data.simulate(_lambda=0, nsteps=1000)
    #print(data.islands.values()[0])
    print("Testing lambda function.")
    data.simulate(_lambda=1)
    #print(data.islands.values()[0])
    #data.run(10)
