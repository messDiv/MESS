from util import *
from collections import OrderedDict
import datetime
import numpy as np
from scipy.stats import kurtosis
from scipy.stats import skew
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
                       ("population_growth", "constant"),
                       ("community_assembly_model", "neutral"),
        ])

        ## Track local communities in this model and colonization rates among them
        ## TODO: I think we want the default metacommunity type to be 'logser', but
        ## the new _sim_metacommunity function takes a little time, so startup is laggy.
        self.metacommunity = MESS.Metacommunity()
        ## Populate the default metacommunity
        self.metacommunity.set_metacommunity()

        self.islands = {}
        self.colonization_matrix = []


    #########################
    ## Housekeeping functions
    #########################
    def __str__(self):
        return "<MESS.Region {}: {}>".format(self.paramsdict["simulation_name"], self.islands.keys())

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
        self.islands[local_community.paramsdict["name"]] = local_community
        local_community.prepopulate(quiet=True)


    def _link_metacommunity(self, metacommunity):
        """ Just import a metacommunity object that's been created externally."""
        LOGGER.debug("Linking metacommunity - {}".format(metacommunity))
        self.metacommunity = metacommunity
        self.metacommunity.set_metacommunity()

        for locname in self.islands.keys():
            self.islands[locname].prepopulate()


    def _get_simulation_outdir(self, prefix=""):
        """ Construct an output directory for a simulation run.
        Make output directory formatted like <output_dir>/<name>-<timestamp><random 2 digit #>
        This will _mostly_ avoid output directory collisions, right?"""

        dirname = prefix + self.paramsdict["simulation_name"]
        outdir = os.path.join(self.paramsdict["project_dir"],\
                              dirname\
                              + "-" + str(time.time()).replace(".", "")[-7:]\
                              + str(np.random.randint(100)))
        if not os.path.exists(outdir):
            os.mkdir(outdir)

        ## Push the outdir into each of the local communities.
        for name, island in self.islands.items():
            island._hackersonly["outdir"] = outdir
            self.islands[name] = island

        return outdir


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

            elif param == "community_assembly_model":
                self.paramsdict[param] = newvalue

            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            ## Do something intelligent here?
            raise


    def _reset_local_communities(self):
        LOGGER.debug("_reset_local_community()")
        for name, island in self.islands.items():
            new = island._copy()
            new.prepopulate(quiet=True)
            self.islands[name] = new


    def write_params(self, outfile=None, outdir=None, force=False):
        """ Write out the parameters of this model to a file properly
        formatted as input for `MESS -p <params.txt>`. A good and
        simple way to share/archive parameter settings for simulations.
        This is also the function that's used by __main__ to
        generate default params.txt files for `MESS -n`
        """
        if outfile is None:
            outfile = "params-"+self.paramsdict["simulation_name"]+".txt"

        if not outdir is None:
            if os.path.exists(outdir):
                outfile = os.path.join(outdir, outfile)
            else:
                raise MESSError(NO_OUTDIR).format(outdir)

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
        loc = MESS.LocalCommunity(name=name, K=K, colrate=c, quiet=quiet)
        ## TODO: Ensure island names aren't dupes
        self._link_local(loc)


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

        ## Preserve any parameter state that has changed:
        self.write_params(force=True)

        ## Open output file. If force then overwrite existing, otherwise just append.
        append = 'a'
        if force:
            append = 'w'
        SIMOUT = open(os.path.join(self.paramsdict["project_dir"], "SIMOUT.txt"), append)
        SIMOUT.write("\t".join(self.islands.values()[0].stats.keys()) + "\n")

        ## Just get all the time values to simulate up front
        ## Doesn't save time really, just makes housekeeping easier
        gens = sample_param_range(self.paramsdict["generations"], nsamps=sims)
        ## Check if we're doing steps or lambda
        do_lambda = False
        if isinstance(gens[0], float):
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

                SIMOUT.write("\t".join(map(str, res.values)) + "\n")
                LOGGER.debug("Finished simulation {} stats:\n{}".format(i, res))
            progressbar(100, 100, " Finished {} simulations\n".format(sims))

        ## Parallelize
        else:
            parallel_jobs = {}
            ## Magic to make the Region() object picklable
            ipyclient[:].use_dill()
            lbview = ipyclient.load_balanced_view()
            for i in xrange(sims):
                parallel_jobs[i] = lbview.apply(simulate, *(self, gens[i]))

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
                    res = passdict[result]
                    SIMOUT.write("\t".join(map(str, res.values)) + "\n")

            ## TODO: Do something more intelligent in the event that any of them fails.
            LOGGER.error("Any failed simulations will report here: {}".format(faildict))


    ## This is the function that will be run inside the cluster engine, so
    ## everything it does needs to be
    ##
    ## The _lambda getting an underscore here is an artifact of the choice
    ## of symbol for %equilibrium colliding with the python anonymous function
    ## name. Don't confuse them.
    ##
    ## TODO: Need to think about how to quantify lambda if there is more
    ## than one island
    ## TODO: Need to implement the force flag here to allow or prevent overwriting
    def simulate(self, _lambda=0, nsteps=0, log_full=False, quiet=True):
        ## Run one realization of this Region. simulate() accepts EITHER _lambda
        ## value to run to OR # of steps to perform, and not both.
        ##  * log_full turns on extremely verbose logging to files in the outdir
        ##  * quiet toggles the more fine grained progress bar and normally
        ##    should be set to True, especially on cluster jobs, otherwise it
        ##    floods the stdout pipe
        LOGGER.debug("Entering simulate(): lambda={}, nsteps={}".format(_lambda, nsteps))
        if _lambda > 0 and nsteps > 0:
            LOGGER.error("simulate accepts only one of either lambda or nsteps args")
            return

        ## Not as big of a deal on ipp simulations, but if you're running on a local computer
        ## the local communities need to get reupped for each simulation.
        self._reset_local_communities()

        ## Get an output directory for dumping data
        outdir = self._get_simulation_outdir()

        step = 0
        ## Create an temp function to test whether we've reached the end of this simulation
        if _lambda > 0:
            ## In the absence of a better strategy just test for lambda in the first local community
            done_func = lambda: self.islands.values()[0]._lambda() >= _lambda
        else:
            done_func = lambda: step >= nsteps

        while not done_func():
            ## This is not ideal functionality, but it at least gives you a sense of progress
            ## Especially don't do this on ipyparallel engines, because it floods the pipe.
            if not quiet:
                progressbar(100, 100, "{0:.4f}".format(self.islands.values()[0]._lambda()))
            for island in self.islands.values():
                island.step()
            step += 1
            if not step % self.paramsdict["recording_period"]:
               for island in self.islands.values():
                    island._log(full=log_full)
        ## TODO: Combine stats across local communities if more than one
        for name, island in self.islands.items():
            statsdf = island.get_stats()

        self.write_params(outdir=outdir)
        return statsdf


    def fancy_plots(self, quiet=True):
        LOGGER.debug("Entering fancy_plots()")

        self.simulate(_lambda=1, log_full=True, quiet=quiet)

        outdir = self._get_simulation_outdir(prefix="fancy-")

        for island in self.islands.values():
            try:
                MESS.plotting.plot_rank_abundance_through_time(outdir,
                                                        island.species_through_time,
                                                        island.lambda_through_time,
                                                        verbose=False)
                MESS.plotting.normalized_pi_dxy_heatmaps(outdir,
                                                        island.species_through_time,
                                                        island.lambda_through_time)
                MESS.plotting.normalized_pi_dxy_heatmaps(outdir,
                                                        island.species_through_time,
                                                        island.lambda_through_time,
                                                        one_d=True)
                MESS.plotting.plot_abundance_vs_colonization_time(outdir,
                                                        island.species_through_time,
                                                        island.lambda_through_time)
            except Exception as inst:
                print("    Exception in fancy_plots() - {}".format(inst))

    def get_trait(self, loc_id):
        return self.metacommunity.community['trait_values'][self.metacommunity.community["ids"] == loc_id]
        self.local_info


    def get_trait_stats(self, local_com):
        local_traits = []
        for id in np.unique(local_com):
            local_traits.append(self.metacommunity.community["trait_values"][self.metacommunity.community["ids"] == id])
        return [np.mean(local_traits),
                np.var(local_traits),
                np.mean(self.metacommunity.community["trait_values"]),
                np.var(self.metacommunity.community["trait_values"]),
                np.mean(self.metacommunity.community["trait_values"]) - np.mean(local_traits),
                np.var(self.metacommunity.community["trait_values"]) - np.var(local_traits),
                kurtosis(local_traits),
                skew(local_traits)]

    def get_weight(self):
        weight =   (self.metacommunity.paramsdict["trait_strength"] *
                    self.metacommunity.paramsdict["trait_rate"])
                    #* self.metacommunity.metcommunity_tree_height)
        return weight

    #def get_local_phy(self):

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
    "population_growth" : "Rate of growth since colonization: exponential/constant",\
    "community_assembly_model" : "Model of Community Assembly: neutral, filtering, competition",\
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

NO_OUTDIR = """
    Error: Attempting to write params to a directory that doesn't exist - {}
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
