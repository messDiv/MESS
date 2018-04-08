from util import *
from collections import OrderedDict
import datetime
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
                       ("data_model", 4),
                       ("generations", 0),
                       ("recording_period", 10000),
                       ("allow_multiple_colonizations", True),
        ])

        ## Track local communities in this model and colonization rates among them
        self.metacommunity = MESS.Metacommunity()
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
        LOGGER.debug("Linking local community - {}".format(local_community))
        local_community.region = self
        local_community.prepopulate()
        self.islands[local_community.paramsdict["name"]] = local_community

    def _link_metacommunity(self, metacommunity):
        """ Just import a metacommunity object that's been created externally."""
        LOGGER.debug("Linking metacommunity - {}".format(metacommunity))
        self.metacommunity = metacommunity


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
                else:
                    self.paramsdict[param] = int(float(newvalue))

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

        ## Write parameters of the metacommunity
        self.metacommunity.write_params(outfile)

        ## Write parameters for each island
        for island in self.islands.values():
            island.write_params(outfile)


    ########################
    ## Model functions/API
    ########################
    def add_local_community(self, name, K, c, quiet=False):
        loc = MESS.LocalCommunity(name, K, c, quiet)
        ## TODO: Ensure island names aren't dupes
        self.islands[name] = loc


    def set_metacommunity(self, meta_type):
        pass    


    def set_colonization_matrix(self, matrix):
        """ Set the matrix that describes colonization rate between local communities."""
        ## TODO: Make sure the matrix is the right shape
        self.colonization_matrix = matrix


    ###############################################
    ## Accessor for sampling from the regional pool
    ###############################################
    def get_nmigrants(self, nmigrants=1):
        """Get a sample of inidividuals from the regional pool.
        Returns a list of species ids"""

        ## TODO: This could potentially be used to draw migrants from
        ## the local island pool as well as the metacommunity
        migrants, trait_vals = self.metacommunity.get_nmigrants(nmigrants)
        return migrants, trait_vals


    def get_migrant(self):
        return self.metacommunity.get_migrant()


    def get_most_abundant(self):
        """Just get the most abundant species from the metacommunity"""

        max_idx = self.metacommunity.community["abundances"].argmax()
        new_species = self.metacommunity.community["ids"][max_idx] 
        trait_value = self.metacommunity.community["trait_values"][max_idx]
        return new_species, trait_value


    def get_abundance(self, species=None):
        """Return abundance of a species in the regional pool."""
        ## This is dumb, metacommunity should be pd
        return self.metacommunity.community["abundances"]\
                [np.where(self.metacommunity.community["ids"] == "t1")][0]

    ## Main function for managing cluster parallelized simulations
    def run(self, sims, force=False, ipyclient=None, quiet=False):
        """ Do the heavy lifting here"""
        print("    Generating {} simulations.".format(sims))

        ## Just get all the time values to simulate up front
        ## Doesn't save time really, just makes housekeeping easier
        gens = self.paramsdict["generations"]
        do_lambda = False
        if type(gens) is tuple:
            ## Generations is a range so uniform random sample for this sim
            gens = np.random.randint(gens[0], gens[1], sims)
        elif gens > 0:
            ## Generations is just a value
            gens = np.array([gens] * sims)
        else:
            ## If doing lambda we only really care about 3ish signficant figures
            gens = np.round(np.random.random(sims), 3)
            do_lambda = True
        LOGGER.debug("Sample of durations of simulations: {}".format(gens[:10]))

        ## Run serially. This will be slow for anything but toy models.
        printstr = " Performing Simulations    | {} |"
        if not ipyclient:
            start = time.time()
            for i in xrange(sims):
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(sims, i, printstr.format(elapsed))
                if not do_lambda:
                    res = self.simulate(nsteps=gens[i])
                else:
                    res = self.simulate(_lambda=gens[i], quiet=quiet)
                LOGGER.debug("Finished simulation {} stats:\n{}".format(i, res))
            progressbar(100, 100, " Finished {} simulations\n".format(sims))

        ## Parallelize
        else:
            parallel_jobs = {}
            ## Magic to make the Region() object picklable
            ipyclient[:].use_dill()
            lbview = ipyclient.load_balanced_view()
            for i in xrange(sims):
                parallel_jobs[i] = lbview.apply(simulate, *(self, gens[i], quiet))

            ## Wait for all jobs to finish
            start = time.time()
            while 1:
                fin = [i.ready() for i in parallel_jobs.values()]
                elapsed = datetime.timedelta(seconds=int(time.time()-start))
                progressbar(len(fin), sum(fin),
                    printstr.format(elapsed))
                time.sleep(0.1)
                if len(fin) == sum(fin):
                    print("")
                    break
            progressbar(100, 100, "\n    Finished {} simulations\n".format(sims))

            faildict = {}
            passdict = {}
            ## Gather results
            for result in parallel_jobs:
                if not parallel_jobs[result].successful():
                    faildict[result] = parallel_jobs[result].metadata.error
                else:
                    passdict[result] = parallel_jobs[result].result()
            print(passdict)
            print(faildict)


    ## This is the function that will be run inside the cluster engine, so
    ## everything it does needs to be
    ##
    ## The _lambda getting an underscore here is an artifact of the choice
    ## of symbol for %equilibrium colliding with the python anonymous function
    ## name. Don't confuse them.
    ##
    ## TODO: Need to think about how to quantify lambda if there is more
    ## than one island
    def simulate(self, _lambda=0, nsteps=0, quiet=True):
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
            ## This is not ideal functionality, but it at least gives you a sense of progress
            if not quiet:
                progressbar(100, 100, "{}%".format(self.islands.values()[0]._lambda()))
            for island in self.islands.values():
                island.step()
            step += 1
        ## TODO: Combine stats across local communities if more than on
        outfile = os.path.join(outdir, "simout.txt")
        statsdf = self.islands.values()[0].get_stats()
        statsdf.to_csv(outfile, na_rep=0, float_format='%.5f')
        return statsdf


def simulate(data, time=time, quiet=True):
    import os
    LOGGER.debug("Entering sim - {} on pid {}\n{}".format(data, os.getpid(), data.paramsdict))
    res = data.simulate(_lambda=time, quiet=quiet)
    LOGGER.debug("Leaving sim - {} on pid {}\n{}".format(data, os.getpid(),\
                                                        [str(x) for x in data.islands.values()]))
    return res

#############################
## Model Parameter Info Dicts
#############################
REGION_PARAMS = {
    "simulation_name" : "The name of this simulation scenario",\
    "project_dir" : "Where to save files",\
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
    data.add_local_community("tmploc", 500, 0.01)
    #print("Testing step function.")
    #data.simulate(_lambda=0, nsteps=1000)
    #print(data.islands.values()[0])
    #data.simulate(_lambda=0, nsteps=1000)
    #print(data.islands.values()[0])
    print("Testing lambda function.")
    data.simulate(_lambda=.4)
    #data.run(10)
    print(data.get_most_abundant())
