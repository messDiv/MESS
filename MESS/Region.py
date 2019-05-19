import datetime
import logging
import numpy as np
import os
import string
import time

from collections import OrderedDict

import MESS
from .stats import *
from .util import *

LOGGER = logging.getLogger(__name__)

class Region(object):
    """
    The MESS Region is the fundamental unit of a batch of simulation scenarios.
    A Region encompasses both a Metacommunity and one or more Local Communities,
    and orchestrates the community assembly process.

    :param str name: The name of this MESS simulation. This is used for
        creating output files.
    :param bool quiet: Whether to print progress of simulations or remain silent.
    :param bool log_files: For each community assembly simulation create a
        a directory in the `outdir`, write the exact parameters for the simulation,
        and dump the megalog to a file. The megalog includes all information
        about the final state of the local community, prior to calculating
        summary statistics per species. Primarily for debugging purposes.
    """

    def __init__(self, name, quiet=False, log_files=False):
        if not name:
            raise MESSError(REQUIRE_NAME)

        self._log_files = log_files

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
                       ("generations", 0),
                       ("community_assembly_model", "neutral"),
                       ("speciation_model", "point_mutation"),
                       ("mutation_rate", 2.2e-8),
                       ("alpha", 2000),
                       ("sequence_length", 570),
        ])

        ## A dictionary for holding prior ranges for values we're interested in
        self._priors = dict([
                        ("alpha", []),
                        ("community_assembly_model", []),
                        ("generations", []),
        ])

        ## elite hackers only internal dictionary, normally you shouldn't mess with this
        ##  * population_growth: Rate of growth since colonization: exponential/constant/harmonic.
        ##      'harmonic' is the only sensible one, so we'll use this as the default always.
        ##  * sgd_dimensions: Number of dimensions for simulated SGD: 1 or 2. 1 indicates
        ##      the presence of only local pi information, and 2 indicates pi for
        ##      metacommunity species, so dxy is calculated.
        ##  * sgd_bins: Number of bins per axis for the SGD histogram
        ##  * recording_period: Number of forward-time ticks between samples for logging
        ##  * abundance_speciation_ratio: Relationship between abundance and speciation 
        ##      probability: proportional, inverse, uniform. This is something that
        ##      Rosie mentioned and could be worth trying in the future. Right now
        ##      it is wholly unimplemented. TODO

        self._hackersonly = dict([
                       ("population_growth", "harmonic"),
                       ("sgd_dimensions", 1),
                       ("sgd_bins", 10),
                       ("recording_period", 10),
                       ("abundance_speciation_ratio", "proportional"),
        ])

        ## Add simple default local and metacommunities
        self.metacommunity = MESS.Metacommunity()
        self.islands = {}
        self._link_local(MESS.LocalCommunity())

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


    def _link_local(self, local_community, add=False):
        """
        Just link a local community that was created externally.
        This is primarily used by __main__ during the initialization process.

        :param local_community: The target local community to link.
        :param bool add: Wether to add this local community to the pool, or to
            link and replace.
        """
        LOGGER.debug("Linking local community - {}".format(local_community))
        local_community._set_region(self)
        ## TODO: This is somewhat hax. If we want to link and replace
        ## here we're just blanking whatever's there and replacing.
        ## In the context of a multi island system this might be ugly.
        if not add:
            self.islands = {}
        self.islands[local_community.paramsdict["name"]] = local_community
        local_community.prepopulate()


    def _link_metacommunity(self, metacommunity):
        """
        Just import a metacommunity object that's been created externally.

        :param metacommunity: The metacommunity to link.
        """
        LOGGER.debug("Linking metacommunity - {}".format(metacommunity))

        self.metacommunity = metacommunity

        for locname in self.islands.keys():
            self.islands[locname].prepopulate()


    def _get_simulation_outdir(self, prefix=""):
        """
        Construct an output directory for a simulation run.
        Make output directory formatted like <output_dir>/<name>-<timestamp><random 2 digit #>
        This will _mostly_ avoid output directory collisions, right?

        :param string prefix: The directory within which to create the
            simulation output directory.
        """

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


    def _paramschecker(self, param, newvalue, quiet=True):
        """
        Check and set parameters for the Region. Raises exceptions when params
        are set to values they should not be.

        :param string param: The parameter to set.
        :param newvalue: The value of the parameter.
        :param bool quiet: Whether to print info.
        """
        ## TODO: This should actually check the values and make sure they make sense
        try:
            ## Cast params to correct types
            if param == "project_dir":
                ## If it already exists then just inform the user that we'll be adding
                ## more simulations to the current project directory
                self.paramsdict[param] = newvalue
                if not os.path.exists(self.paramsdict["project_dir"]):
                    os.mkdir(self.paramsdict["project_dir"])
                else:
                    if not quiet:
                        print("  Project directory exists. Additional simulations will be appended.")

            elif param == "generations":
                tup = tuplecheck(newvalue, dtype=float)

                ## If specifying in generations then cast to int, otherwise it's
                ## lambda (between 0 and 1) so leave it as float.
                if isinstance(tup, tuple):
                    if tup[0] > tup[1]:
                        raise MESSError("Generations parameters malformed: {} must be < {}".format(tup[0], tup[1]))
                    if tup[0] > 1:
                        tup = tuple(map(int, tup))
                    self._priors[param] = tup
                    self.paramsdict[param] = tup
                else:
                    if tup > 1:
                        tup = int(tup)
                    self.paramsdict[param] = tup

            elif param == "community_assembly_model":
                if newvalue == "*":
                    self._priors[param] = ["neutral", "filtering", "competition"]
                    newvalue = np.random.choice(self._priors[param])
                self.paramsdict[param] = newvalue

            elif param == "speciation_model":
                self.paramsdict[param] = newvalue

            elif param == "mutation_rate":
                self.paramsdict[param] = float(newvalue)

            elif param == "alpha":
                tup = tuplecheck(newvalue, dtype=int)
                if isinstance(tup, tuple):
                    self._priors[param] = tup
                    self.paramsdict[param] = sample_param_range(tup)[0]
                else:
                    self.paramsdict[param] = tup

            elif param == "sequence_length":
                self.paramsdict[param] = float(newvalue)
            else:
                self.paramsdict[param] = newvalue
        except Exception as inst:
            ## Do something intelligent here?
            raise


    def _record_local_speciation(self, sname, trait_value):
        """
        Local speciation events need to be recorded, for example the trait 
        value needs to get added to the metacommunity trait value list. Maybe
        other stuff needs to happen at speciation time, so we create a function.
        This function is called by LocalCommunity.speciate().

        :param string sname: The ID of the new species.
        :param float trait_value: The trait value of the new species.
        """
        self.metacommunity._update_species_pool(sname, trait_value)


    def _reset_local_communities(self):
        """
        Flip the local community. Basically make a copy, resample any parameters
        that were specified with priors and prepopulate it.
        """
        LOGGER.debug("_reset_local_community()")
        for name, island in self.islands.items():
            new = island._copy()
            new.prepopulate()
            self.islands[name] = new


    def _reset_metacommunity(self):
        """
        Flip the metacommunity. Just regenerate a new metacommunity, resampling
        from any priors that were specified.
        """
        ## Calling set_metacommunity() again will regenerate a new
        ## metacommunity using the same parameters each time.
        self.metacommunity.set_metacommunity(resample=True)


    ## Getting parameters header and parameters carves off
    ## the simulation name and the project directory
    def _get_params_header(self):
        return list(self.paramsdict.keys())[2:]


    def _get_params_values(self):
        return list(self.paramsdict.values())[2:]


    def set_param(self, param, value, quiet=True):
        """
        A convenience function for setting parameters in the API mode, which
        turns out to be a little annoying if you don't provide this. With
        the set_param method you can set parameters on the Region, the
        Metacommunity, or the LocalCommunity. Simply pass the parameter
        name and the value, and this method identifies the appropriate target
        parameter.

        :param string param: The name of the parameter to set.
        :param value: The value of the parameter to set.
        :param bool quiet: Whether to print info to the console.
        """
        try:
            self = set_params(self, param, value, quiet)
        except:
            try:
                self.metacommunity = set_params(self.metacommunity, param, value, quiet)
            except:
                try:
                    name, loc = self.islands.items()[0]
                    self.islands[name] = set_params(loc, param, value, quiet)
                except:
                    raise MESSError("Bad param/value {}/{}".format(param, value))


    def write_params(self, outfile=None, outdir=None, full=False, force=False):
        """
        Write out the parameters of this model to a file properly formatted as
        input for `MESS -p <params.txt>`. A good and simple way to 
        share/archive parameter settings for simulations. This is also the 
        function that's used by __main__ to generate default params.txt files
        for `MESS -n`.

        :param string outfile: The name of the params file to generate. If not
            specified this will default to `params-<Region.name>.txt`.
        :param string outdir: The directory to write the params file to. If not
            specified this will default to the project_dir.
        :param bool full: Whether to write out only the parameters of the
            specific parameter values of this Region, or to write out the
            parameters including prior ranges for parameter values..
        :param bool force: Whether to overwrite if a file already exists.
        """
        if outfile is None:
            outfile = "params-"+self.paramsdict["simulation_name"]+".txt"

        ## If outdir is blank then default to writing to the project dir
        if outdir is None:
            outdir = self.paramsdict["project_dir"]
        elif not os.path.exists(outdir):
            raise MESSError(NO_OUTDIR).format(outdir)

        outfile = os.path.join(outdir, outfile)
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
                paramvalue = str(val)

                ## If it's one of the params with a prior, and if the prior is not
                ## empty and if writing out full, then write the prior, and not
                ## the sampled value
                if full:
                    if key in list(self._priors.keys()):
                        ## The prior on community assembly model is a little goofy
                        ## since it's a list, and not a search range
                        if key == "community_assembly_model" and self._priors[key]:
                            paramvalue = "*"
                        elif self._priors[key]:
                            paramvalue = "-".join([str(i) for i in self._priors[key]])

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
        self.metacommunity._write_params(outfile, full=full)

        ## Write parameters for each island
        for island in self.islands.values():
            island.write_params(outfile, full=full)


    ########################
    ## Model functions/API
    ########################
    def add_local_community(self, name, J, m, quiet=False):
        loc = MESS.LocalCommunity(name=name, J=J, m=m, quiet=quiet)
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
    def _get_nmigrants(self, nmigrants=1):
        """Get a sample of inidividuals from the regional pool.
        Returns a list of species ids"""

        ## TODO: This could potentially be used to draw migrants from
        ## the local island pool as well as the metacommunity
        migrants, trait_vals = self.metacommunity._get_nmigrants(nmigrants)
        return migrants, trait_vals


    def _get_migrant(self):
        return self.metacommunity._get_migrant()


    def _get_most_abundant(self):
        """Just get the most abundant species from the metacommunity"""
        max_idx = self.metacommunity.community["abundances"].argmax()
        new_species = self.metacommunity.community["ids"][max_idx]
        trait_value = self.metacommunity.community["trait_values"][max_idx]

        ## Nudge the trait value of the initial colonizer or else
        ## for filtering models you can get in big trouble (sims run forever)
        if abs(trait_value - self.metacommunity._hackersonly["filtering_optimum"]) < 1:
            trait_value += 1
            self.metacommunity.community["trait_values"][max_idx] = trait_value

        return new_species, trait_value


    def _get_abundance(self, species=None):
        """Return abundance of a species in the regional pool."""
        ## This is dumb, metacommunity should be pd
        return self.metacommunity.community["abundances"]\
                [np.where(self.metacommunity.community["ids"] == species)][0]


    def get_species_params(self):
        return {"mutation_rate":self.paramsdict["mutation_rate"],
                "alpha": self.paramsdict["alpha"],
                "sequence_length":self.paramsdict["sequence_length"]}


    ## Main function for managing cluster parallelized simulations
    def run(self, sims, force=False, ipyclient=None, quiet=False):
        """ Do the heavy lifting here"""
        if not quiet: print("    Generating {} simulation(s).".format(sims))

        if not os.path.exists(self.paramsdict["project_dir"]):
            os.mkdir(self.paramsdict["project_dir"])

        simfile = os.path.join(self.paramsdict["project_dir"], "SIMOUT.txt")
        ## Open output file. If force then overwrite existing, otherwise just append.
        append = 'a'
        if force:
            append = 'w'
            ## Prevent from shooting yourself in the foot with -f
            try:
                os.rename(simfile, simfile+".bak")
            except FileNotFoundError:
                ## If the simfile doesn't exist catch the error and move on
                pass
        ## Decide whether to print the header, if stuff is already in there then
        ## don't print the header, unless you're doing force because this opens
        ## in overwrite mode.
        params = self.metacommunity._get_params_header() +\
                 self._get_params_header() +\
                 list(self.islands.values())[0]._get_params_header()
        header = "\t".join(params + MESS.stats._get_sumstats_header(sgd_bins=self._hackersonly["sgd_bins"],\
                                                                    sgd_dims=self._hackersonly["sgd_dimensions"],
                                                                    metacommunity_traits=self.metacommunity._get_trait_values())) + "\n"
        LOGGER.debug("SIMOUT header - {}".format(header))
        if len(open(simfile, 'a+').readline()) > 0 and not force:
            header = ""
        SIMOUT = open(simfile, append)
        SIMOUT.write(header)

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
            try:
                for i in range(sims):
                    elapsed = datetime.timedelta(seconds=int(time.time()-start))
                    if not quiet: progressbar(sims, i, printstr.format(elapsed))

                    if not do_lambda:
                        res = self.simulate(nsteps=gens[i])
                    else:
                        res = self.simulate(_lambda=gens[i], quiet=quiet)

                    SIMOUT.write(res + "\n")
                    LOGGER.debug("Finished simulation {} stats:\n{}".format(i, res))
            except KeyboardInterrupt as inst:
                print("\n    Cancelling remaining simulations")
            if not quiet: progressbar(100, 100, " Finished {} simulations\n".format(i))

        ## Parallelize
        else:
            parallel_jobs = {}

            ## store ipyclient engine pids to the Assembly so we can
            ## hard-interrupt them later if assembly is interrupted.
            ## Only stores pids of engines that aren't busy at this moment,
            ## otherwise it would block here while waiting to find their pids.
            self._ipcluster["pids"] = {}
            for eid in ipyclient.ids:
                engine = ipyclient[eid]
                if not engine.outstanding:
                    pid = engine.apply(os.getpid).get()
                    self._ipcluster["pids"][eid] = pid

            ## Magic to make the Region() object picklable
            ipyclient[:].use_dill()
            lbview = ipyclient.load_balanced_view()
            for i in range(sims):
                if do_lambda:
                    parallel_jobs[i] = lbview.apply(simulate, self, gens[i], 0)
                else:
                    parallel_jobs[i] = lbview.apply(simulate, self, 0, gens[i])

            ## Wait for all jobs to finish
            start = time.time()
            while 1:
                try:
                    fin = [i.ready() for i in parallel_jobs.values()]
                    elapsed = datetime.timedelta(seconds=int(time.time()-start))
                    if not quiet: progressbar(len(fin), sum(fin),
                        printstr.format(elapsed))
                    time.sleep(0.1)
                    if len(fin) == sum(fin):
                        print("")
                        break
                except KeyboardInterrupt as inst:
                    print("\n    Cancelling remaining simulations.")
                    break
            if not quiet: progressbar(100, 100, "\n    Finished {} simulations\n".format(len(fin)))

            faildict = {}
            passdict = {}
            ## Gather results
            for result in parallel_jobs:
                try:
                    if not parallel_jobs[result].successful():
                        faildict[result] = parallel_jobs[result].metadata.error
                    else:
                        passdict[result] = parallel_jobs[result].result()
                        res = passdict[result]
                        SIMOUT.write(res + "\n")
                except Exception as inst:
                    LOGGER.error("Caught a failed simulation - {}".format(inst))
                    ## Don't let one bad apple spoin the bunch,
                    ## so keep trying through the rest of the asyncs

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
        if _lambda == 0 and nsteps == 0:
            msg = "Either _lambda or nsteps must be specified."
            raise MESSError(msg)

        ## If priors for alpha are set then resample alpha
        if self._priors["alpha"]:
            self.paramsdict["alpha"] = sample_param_range(self._priors["alpha"])[0]
            LOGGER.debug("alpha - {}".format(self.paramsdict["alpha"]))

        if self._priors["community_assembly_model"]:
            self.paramsdict["community_assembly_model"] = np.random.choice(self._priors["community_assembly_model"])

        ## Flip the metacommunity per simulation so we get new draws of trait values.
        ## This is a little slow for logser, and also performance scales with metacommunity size
        self._reset_metacommunity()

        ## Not as big of a deal on ipp simulations, but if you're running on a local computer
        ## the local communities need to get reupped for each simulation.
        self._reset_local_communities()

        if self._log_files:
            ## Get an output directory for dumping data
            outdir = self._get_simulation_outdir()
            self.write_params(outdir=outdir)

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
            if not step % self._hackersonly["recording_period"]:
               for island in self.islands.values():
                    island._log(full=log_full)
        ## TODO: Combine stats across local communities if more than one
        for name, island in self.islands.items():
            statsdf = island.get_stats()

        ## Paste regional parameters on the front of the local community
        ## parameters and simulations
        regional_params = self.metacommunity._get_params_values() +\
                            self._get_params_values()
        tmpsimout = regional_params + list(statsdf.T.values[0])
        simout = []
        for x in tmpsimout:
            try:
                simout.append(np.round(x, 5))
            except:
                simout.append(x)
        simout = "\t".join(map(str, np.array(simout)))

        return simout


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
        try:
            trt = self.metacommunity.trait_dict[loc_id]
            #trt = self.metacommunity.community['trait_values'][self.metacommunity.community["ids"] == loc_id][0]
        except Exception as inst:
            raise MESSError("Species has no trait value: {}".format(loc_id))
        return trt


    def get_trait_mean(self, local_com):
        try:
            sp = list(set(local_com))
            mask = np.isin(self.metacommunity.community["ids"], sp)
            local_traits = self.metacommunity.community["trait_values"][mask]

        except Exception as inst:
            raise MESSError("Problem getting traits from local community: {}".format(inst))

        return np.mean(local_traits)


    def get_phy_stats(self, tree):
        total = []
        for edge in tree.postorder_edge_iter():
            if edge.length is not None:
                total.append(edge.length)
        return [np.mean(total), np.var(total), len(total), sum(total)]


    #def get_local_phy(self):


def simulate(data, _lambda=0, nsteps=0, quiet=True):
    import os
    LOGGER.debug("Entering sim - {} on pid {}\n{}".format(data, os.getpid(), data.paramsdict))
    res = data.simulate(_lambda=_lambda, nsteps=nsteps, quiet=quiet)
    LOGGER.debug("Leaving sim - {} on pid {}\n{}".format(data, os.getpid(),\
                                                        [str(x) for x in data.islands.values()]))
    return res

#############################
## Model Parameter Info Dicts
#############################
REGION_PARAMS = {
    "simulation_name" : "The name of this simulation scenario",\
    "project_dir" : "Where to save files",\
    "generations" : "Duration of simulations. Specify int range or 0 for lambda.",\
    "community_assembly_model" : "Model of Community Assembly: neutral, filtering, competition",\
    "speciation_model" : "Type of speciation process: none, point_mutation, protracted, random_fission",\
    "mutation_rate" : "Mutation rate scaled per base per generation",\
    "alpha" : "Abundance/Ne scaling factor",\
    "sequence_length" : "Length in bases of the sequence to simulate",\
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

