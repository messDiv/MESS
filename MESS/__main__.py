""" the main CLI for calling MESS """



import pkg_resources
import argparse
import logging
import atexit
import time
import sys
import os
import MESS

from MESS.util import *
from MESS.parallel import *

LOGGER = logging.getLogger(__name__)


def parse_params(args):
    """ Parse the params file and return Region params as a dict
    and a list of dicts for params for each island."""

    ## check that params.txt file is correctly formatted.
    try:
        with open(args.params) as paramsin:
            plines = paramsin.read()
    except IOError as _:
        sys.exit("  No params file found")

    ## Params file is sectioned off by region and then by n local communities
    ## do each section independently
    plines = plines.split("------- ")
    ## Get seed for global simulations
    seed = [i.split("##")[0].strip() for i in plines[0].split("\n")[1:] if not i.strip() == ""][0]

    ## Get Region params and make into a dict, ignore all blank lines
    items = [i.split("##")[0].strip() for i in plines[1].split("\n")[1:] if not i.strip() == ""]
    keys = list(MESS.Region('null', quiet=True).paramsdict.keys())
    region_params = {str(i):j for i, j in zip(keys, items)}

    LOGGER.debug("Region params - {}".format(region_params))

    ## Get metacommunity params and make a dict
    items = [i.split("##")[0].strip() for i in plines[2].split("\n")[1:] if not i.strip() == ""]
    keys = list(MESS.Metacommunity(meta_type="uniform", quiet=True).paramsdict.keys())
    meta_params = {str(i):j for i, j in zip(keys, items)}
    LOGGER.debug("Metacommunity params - {}".format(meta_params))

    island_params = []
    for i, params in enumerate(plines[3:]):
        items = [i.split("##")[0].strip() for i in params.split("\n")[1:] if not i.strip() == ""]
        keys = list(MESS.LocalCommunity('null', quiet=True).paramsdict.keys())
        island_dict = {str(i):j for i, j in zip(keys, items)}
        island_params.append(island_dict)

    LOGGER.debug("All island params - {}".format(island_params))

    return seed, region_params, meta_params, island_params


def getregion(args, seed, region_params, meta_params, island_params):
    """ 
    loads simulation or creates a new one and set its params as
    read in from the params file. Does not launch ipcluster. 
    """

    project_dir = MESS.util._expander(region_params['project_dir'])
    LOGGER.debug("project_dir: {}".format(project_dir))
    sim_name = region_params['simulation_name']

    ## make sure the working directory exists.
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)

    data = MESS.Region(sim_name, quiet=args.quiet, log_files=args.log_files)

    ## Populate the parameters of the Region
    for param in region_params:
        try:
            data = set_params(data, param, region_params[param], quiet=args.quiet)
        except Exception as inst:
            print("Error in __main__.getregion(): {}".format(inst))
            sys.exit(-1)

    meta = MESS.Metacommunity()
    for param in meta_params:
        try:
            meta = set_params(meta, param, meta_params[param], quiet=args.quiet)
        except Exception as inst:
            print("Error in __main__.getregion(): {}".format(inst))
            sys.exit(-1)
    data._link_metacommunity(meta)
        
    ## Populate the islands w/in the region
    for island in island_params:
        loc = MESS.LocalCommunity(island["name"], quiet=True)
        for param in island:
            try:
                loc = set_params(loc, param, island[param], quiet=args.quiet)
            except Exception as inst:
                print("    Error setting param {} - {}".format(param, island[param]))
                sys.exit(-1)
        data._link_local(loc)

    return data


def _check_version():
    """ Test if there's a newer version and nag the user to upgrade."""
    ## TODO: This doesn't actually work yet

    from urllib.request import urlopen

    from distutils.version import LooseVersion

    try:
        htmldat = urlopen("https://anaconda.org/MESS/mess").readlines()
        curversion = next((x for x in htmldat if "subheader" in x), None).split(">")[1].split("<")[0]
        if LooseVersion(MESS.__version__) < LooseVersion(curversion):
            msg = """
  A new version of MESS is available (v.{}). To upgrade run:

    conda install -c mess mess\n""".format(curversion)
            print(MESS_HEADER + "\n" + msg)
        else:
            pass
            #print("You are up to date")
    except Exception as inst:
        ## Always fail silently
        pass


def showstats(region_params, meta_params, island_params):
    """ This should be used to display anything that might be useful. """

    project_dir = region_params["project_dir"]
    if not project_dir:
        project_dir = "./"

    print("\n  Summary stats of Region: {}".format(region_params["simulation_name"]) \
         +"\n  ------------------------------------------------")

    print("   TODO: Put interesting stuff here.\n")


def parse_command_line():
    """ Parse CLI args."""

    ## create the parser
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    MESS -n data                       ## create new file called params-data.txt 
    MESS -p params-data.txt            ## run MESS with settings in params file
    MESS -p params-data.txt -f         ## run MESS, overwrite existing data.
    """)

    ## add arguments 
    parser.add_argument('-n', metavar='new', dest="new", type=str, 
        help="create new file 'params-{new}.txt' in current directory")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str,
        help="path to params file simulations: params-{name}.txt")

    parser.add_argument("-s", metavar="sims", dest="sims",
        type=int, default=0,
        help="Generate specified number of simulations")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=-1,
        help="number of CPU cores to use (Default=0=All)")

    parser.add_argument('-r', action='store_true', dest="results",
        help="show status of this simulation run")

    parser.add_argument('-e', metavar='empirical', dest="empirical", type=str, 
        default=None, 
        help="Validate and import empirical data.")

    parser.add_argument('-f', action='store_true', dest="force",
        help="force overwrite of existing data")

    parser.add_argument('-q', action='store_true', dest="quiet",
        help="do not print to stderror or stdout.")

    parser.add_argument('-Q', action='store_true', dest="very_quiet",
        help="do not print anything ever.")

    parser.add_argument('-d', action='store_true', dest="debug",
        help="print lots more info to mess_log.txt.")

    parser.add_argument('-l', action='store_true', dest="log_files",
        help="Write out lots of information in one directory per simulation.")
    
    parser.add_argument('-v', action='version', 
        version=str(pkg_resources.get_distribution('mess')),
        help=argparse.SUPPRESS)

    parser.add_argument("--ipcluster", metavar="ipcluster", dest="ipcluster",
        type=str, nargs="?", const="default",
        help="connect to ipcluster profile")

    parser.add_argument("--fancy-plots", action='store_true', dest="fancy_plots",
        help="Construct fancy plots and animated gifs.")

    ## This is kind of cool, but i'm not super married to it yet
    ## It's a little pissy and would require rearranging a bunch of stuff
    #subparsers = parser.add_subparsers(title="Optional subcommands", metavar="")
    #list_parser = subparsers.add_parser('simulate', help='Generate profuse simulations')
    #list_parser.add_argument('-a', metavar='abundances', dest="abundances",
    #    help='File containing abundance data')
    #list_parser = subparsers.add_parser('import', help='Validate and import empirical data')
    #list_parser.add_argument('-a', metavar='abundances', dest="abundances",
    #    help='File containing abundance data')
    #list_parser.add_argument('-s', metavar="sequence_dir", dest="sequence_dir",
    #    help='Directory containing fasta files, one per species')
    #list_parser = subparsers.add_parser('ABC', help='Perform ABC model selection / parameter estimation')
    #list_parser = subparsers.add_parser('RF', help='Perform Random Forest model selection / parameter estimation')

    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()

    ## very_quiet implies quiet
    if args.very_quiet:
        args.quiet = True

    if not any(x in ["params", "new"] for x in list(vars(args).keys())):
        print("\nBad arguments: must include at least one of"\
                +" `-p` or `-n`\n")
 
        #parser.print_help()
        sys.exit(1)

    return args


def do_sims(data, args):
    ## if ipyclient is running (and matched profile) then use that one
    if args.ipcluster:
        ipyclient = ipp.Client(cluster_id=args.ipcluster)
        data._ipcluster["cores"] = len(ipyclient.ids)
        if not args.quiet:
            print("    Attached to cluster {} w/ {} engines.".format(args.ipcluster, data._ipcluster["cores"]))

    ## if not then we need to register and launch an ipcluster instance
    elif args.cores >= 0:
        ## set CLI ipcluster terms
        ipyclient = None
        data._ipcluster["cores"] = args.cores if args.cores else detect_cpus()
        data._ipcluster["engines"] = "Local"
        ## register to have a cluster-id with "ip- name"
        data = register_ipcluster(data)
        ipyclient = get_client(**data._ipcluster)
        print(cluster_info(ipyclient))
    ## If args.cores is negative then don't use ipcluster
    else:
        ipyclient = None
        if not args.quiet: print("    Parallelization disabled.")
    try:
        ## Do stuff here
        data.run(
            sims=args.sims,
            force=args.force,
            ipyclient=ipyclient,
            quiet=args.quiet)
    except KeyboardInterrupt as inst:
        print("\n  Keyboard Interrupt by user")
        LOGGER.info("assembly interrupted by user.")
    except MESSError as inst:
        LOGGER.error("MESSError: %s", inst)
        print("\n  Encountered an error (see details in ./mess_log.txt)"+\
              "\n  Error summary is below -------------------------------"+\
              "\n{}".format(inst))
    except Exception as inst:
        LOGGER.error(inst)
        print("\n  Encountered an unexpected error (see ./mess_log.txt)"+\
              "\n  Error message is below -------------------------------"+\
              "\n{}".format(inst))
    finally:
        try:
            ## can't close client if it was never open
            if ipyclient:
                print("Clean up ipcluster {}".format(ipyclient))
                ## send SIGINT (2) to all engines
                try:
                    ipyclient.abort(block=False)
                    time.sleep(1)
                    for engine_id, pid in list(data._ipcluster["pids"].items()):
                        LOGGER.debug("  Cleaning up ipcluster engine/pid {}/{}".format(engine_id, pid))
                        if ipyclient.queue_status()[engine_id]["tasks"]:
                            os.kill(pid, 2)
                            LOGGER.info('interrupted engine {} w/ SIGINT to {}'\
                                    .format(engine_id, pid))
                    time.sleep(1)
                except ipp.NoEnginesRegistered as inst:
                    LOGGER.debug("No engines registered: {}".format(inst))

                ## if CLI, stop jobs and shutdown. Don't use _cli here 
                ## because you can have a CLI object but use the --ipcluster
                ## flag, in which case we don't want to kill ipcluster.
                if 'mess-cli' in data._ipcluster["cluster_id"]:
                    LOGGER.info("  shutting down engines")
                    ipyclient.shutdown(hub=True, block=False)
                    ipyclient.close()
                    LOGGER.info("  finished shutdown")
                else:
                    if not ipyclient.outstanding:
                        ipyclient.purge_everything()
                    else:
                        ## nanny: kill everything, something bad happened
                        ipyclient.shutdown(hub=True, block=False)
                        ipyclient.close()
                        print("\nwarning: ipcluster shutdown and must be restarted")

        ## if exception is close and save, print and ignore
        except Exception as inst2:
            print("warning: error during shutdown:\n{}".format(inst2))
            LOGGER.error("shutdown warning: %s", inst2)


def fancy_plots(data, args):
    """ Fancy plots are nice to generate once in a while, but you don't want to
    be making them for all your actual simulations because they take a long time
    and plus one or two will be very representative of a model."""

    try:
        if args.ipcluster or args.cores >=0:
            print("    Generating fancy plots so running locally on one core.")

        data.fancy_plots(quiet=args.quiet)
    except KeyboardInterrupt as inst:
        print("\n  Keyboard Interrupt by user")
        LOGGER.info("assembly interrupted by user.")
    except MESSError as inst:
        LOGGER.error("MESSError: %s", inst)
        print("\n  Encountered an error (see details in ./mess_log.txt)"+\
              "\n  Error summary is below -------------------------------"+\
              "\n{}".format(inst))
    except Exception as inst:
        LOGGER.error(inst)
        print("\n  Encountered an unexpected error (see ./mess_log.txt)"+\
              "\n  Error message is below -------------------------------"+\
              "\n{}".format(inst))

def main():
    """ main function """
    MESS.__interactive__ = 0  ## Turn off API output

    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    if not args.very_quiet: print(MESS_HEADER)

    ## Turn the debug output written to mess_log.txt up to 11!
    ## Clean up the old one first, it's cleaner to do this here than
    ## at the end (exceptions, etc)
    if os.path.exists(MESS.__debugflag__):
        os.remove(MESS.__debugflag__)

    if args.debug:
        if not args.quiet: print("\n  ** Enabling debug mode **\n")
        MESS._debug_on()
        atexit.register(MESS._debug_off)

    ## create new paramsfile if -n
    if args.new:
        ## Create a tmp assembly, call write_params to make default params.txt
        try:
            tmpassembly = MESS.Region(args.new, quiet=True)
            tmpmeta = MESS.Metacommunity()
            tmpassembly._link_metacommunity(tmpmeta)
            tmplocal = MESS.LocalCommunity("island1", quiet=True)
            tmpassembly._link_local(tmplocal)
            tmpassembly.write_params("params-{}.txt".format(args.new), outdir="./", 
                                     force=args.force)
        except Exception as inst:
            print("Error creating new params file: {}".format(inst))
            sys.exit(2)

        print("\n  New file 'params-{}.txt' created in {}\n".\
               format(args.new, os.path.realpath(os.path.curdir)))
        sys.exit(2)


    ## if params then must provide action argument with it
    if args.params:
        if not any([args.results, args.sims, args.fancy_plots]):
            print(MESS_USAGE)
            sys.exit(2)

    if not args.params:
        if any([args.results, args.sims]):
            print(MESS_USAGE)
            sys.exit(2)

    ## Log the current version. End run around the LOGGER
    ## so it'll always print regardless of log level.
    with open(MESS.__debugfile__, 'a') as logfile:
        logfile.write(MESS_HEADER)
        logfile.write("\n  Begin run: {}".format(time.strftime("%Y-%m-%d %H:%M")))
        logfile.write("\n  Using args {}".format(vars(args)))
        logfile.write("\n  Platform info: {}\n".format(os.uname()))

    ## create new Region or load existing Region
    if args.params:
        seed, region_params, meta_params, island_params = parse_params(args)
        LOGGER.debug("region params - {}\nmetacommunity params - {}\nisland params - {}"\
                    .format(region_params, meta_params, island_params))

        ## launch or load Region with custom profile/pid
        data = getregion(args, seed, region_params, meta_params, island_params)

        ## Validate, format, and import empirical data
        if args.empirical:
            import_empirical(args.empirical)

        ## Print results and exit immediately
        if args.results:
            showstats(region_params, meta_params, island_params)

        elif args.fancy_plots:
            fancy_plots(data, args)

        ## Generate numerous simulations
        elif args.sims:
            ## Only blank the log file if we're actually going to do real
            ## work. In practice the log file should never get this big.
            if os.path.exists(MESS.__debugfile__):
                if os.path.getsize(MESS.__debugfile__) > 50000000:
                    with open(MESS.__debugfile__, 'w') as clear:
                        clear.write("file reset")

            if not args.quiet:
                print("\n    {}".format(data))

            try:
                ## Some binary logic trickery here. If very_quiet is False, and
                ## quiet is True, then we still print this information here. If
                ## very_quiet is false, and quiet is false, the Region.run() 
                ## prints this information. This allows for API mode to have 
                ## the power of the quiet flag as well. 
                if not args.very_quiet and args.quiet: print("    Generating {} simulation(s).".format(args.sims))

                do_sims(data, args)
            except Exception as inst:
                print("  Unexpected error - {}".format(inst))


BAD_PARAMETER_ERROR = """
    Malformed params file: {}
    Bad parameter {} - {}
    {}"""

MESS_HEADER = \
"\n -------------------------------------------------------------"+\
"\n  MESS [v.{}]".format(MESS.__version__)+\
"\n  Massive Eco-Evolutionary Synthesis Simulations"+\
"\n -------------------------------------------------------------"


MESS_USAGE = """
    Must provide action argument along with -p argument for params file. 
    e.g., MESS -p params-test.txt -s 10000           ## run 10000 simulations
    e.g., MESS -p params-test.txt -r                 ## shows results
    """

if __name__ == "__main__": 
    main()
