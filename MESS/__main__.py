#!/usr/bin/env python2
""" the main CLI for calling ipyrad """

from __future__ import print_function, division  # Requires Python 2.7+

import pkg_resources
import argparse
import logging
import atexit
import time
import sys
import os
import MESS

def parse_params(args):
    """ Parse the params file args, create and return Region object."""

    ## check that params.txt file is correctly formatted.
    try:
        with open(args.params) as paramsin:
            plines = paramsin.readlines()
    except IOError as _:
        sys.exit("  No params file found")

    ## make into a dict. Ignore blank lines at the end of file
    ## Really this will ignore all blank lines
    items = [i.split("##")[0].strip() for i in plines[1:] if not i.strip() == ""]

    #keys = [i.split("]")[-2][-1] for i in plines[1:]]
    #keys = range(len(plines)-1)
    keys = mess.MESS('null', quiet=True).paramsdict.keys()
    parsedict = {str(i):j for i, j in zip(keys, items)}
    return parsedict



def getassembly(args, parsedict):
    """ 
    loads assembly or creates a new one and set its params from 
    parsedict. Does not launch ipcluster. 
    """

    ## Creating an assembly with a full path in the name will "work"
    ## but it is potentially dangerous, so here we have assembly_name
    ## and assembly_file, name is used for creating new in cwd, file is
    ## used for loading existing.
    ##
    ## Be nice if the user includes the extension.
    #assembly_name = parsedict['0']
    project_dir = ip.core.assembly._expander(parsedict['project_dir'])
    assembly_name = parsedict['assembly_name']
    assembly_file = os.path.join(project_dir, assembly_name)

    ## Region creation will handle error checking  on
    ## the format of the assembly_name

    ## make sure the working directory exists.
    if not os.path.exists(project_dir):
        os.mkdir(project_dir)

    try:
        ## If 1 and force then go ahead and create a new assembly
        if ('1' in args.steps) and args.force:
            data = ip.Assembly(assembly_name, cli=True)
        else:
            data = ip.load_json(assembly_file, cli=True)
            data._cli = True

    except IPyradWarningExit as _:
        ## if no assembly is found then go ahead and make one
        if '1' not in args.steps:
            raise IPyradWarningExit(\
                "  Error: You must first run step 1 on the assembly: {}"\
                .format(assembly_file))
        else:
            ## create a new assembly object
            data = ip.Assembly(assembly_name, cli=True)

    ## for entering some params...
    for param in parsedict:

        ## trap assignment of assembly_name since it is immutable.
        if param == "assembly_name":
            ## Raise error if user tried to change assembly name
            if parsedict[param] != data.name:
                data.set_params(param, parsedict[param])
        else:
            ## all other params should be handled by set_params
            try:
                data.set_params(param, parsedict[param])
            except IndexError as _:
                print("  Malformed params file: {}".format(args.params))
                print("  Bad parameter {} - {}".format(param, parsedict[param]))
                sys.exit(-1)
    return data


def _check_version():
    """ Test if there's a newer version and nag the user to upgrade."""
    import urllib2
    from distutils.version import LooseVersion

    header = \
    "\n -------------------------------------------------------------"+\
    "\n  MESS [v.{}]".format(mess.__version__)+\
    "\n  Massive Eco-Evolutionary Synthesis Simulations"+\
    "\n -------------------------------------------------------------"

    try:
        htmldat = urllib2.urlopen("https://anaconda.org/ipyrad/ipyrad").readlines()
        curversion = next((x for x in htmldat if "subheader" in x), None).split(">")[1].split("<")[0]
        if LooseVersion(mess.__version__) < LooseVersion(curversion):
            msg = """
  A new version of ipyrad is available (v.{}). To upgrade run:

    conda install -c ipyrad ipyrad\n""".format(curversion)
            print(header + "\n" + msg)
        else:
            pass
            #print("You are up to date")
    except Exception as inst:
        ## Always fail silently
        pass

def parse_command_line():
    """ Parse CLI args."""

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 
    MESS -n data                       ## create new file called params-data.txt 
    MESS -p params-data.txt            ## run MESS with settings in params file
    MESS -p params-data.txt -f         ## run MESS, overwrite existing data.
    """)

    ## add arguments 
    parser.add_argument('-v', '--version', action='version', 
        version=str(pkg_resources.get_distribution('mess')))

    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")

    parser.add_argument('-q', "--quiet", action='store_true',
        help="do not print to stderror or stdout.")

    parser.add_argument('-d', "--debug", action='store_true',
        help="print lots more info to ipyrad_log.txt.")

    parser.add_argument('-n', metavar='new', dest="new", type=str, 
        default=None, 
        help="create new file 'params-{new}.txt' in current directory")

    parser.add_argument('-p', metavar='params', dest="params",
        type=str, default=None,
        help="path to params file simulations: params-{assembly_name}.txt")

    parser.add_argument("-c", metavar="cores", dest="cores",
        type=int, default=0,
        help="number of CPU cores to use (Default=0=All)")

    #parser.add_argument("--ipcluster", metavar="ipcluster", dest="ipcluster",
    #    type=str, nargs="?", const="default",
    #    help="connect to ipcluster profile (default: 'default')")

    ## if no args then return help message
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    ## parse args
    args = parser.parse_args()

    if not any(x in ["params", "new"] for x in vars(args).keys()):
        print("Bad arguments: must include at least one of"\
                +"`-p` or `-n`\n")
        parser.print_help()
        sys.exit(1)

    return args


def main():
    """ main function """
    ## parse params file input (returns to stdout if --help or --version)
    args = parse_command_line()

    ## create new paramsfile if -n
    if args.new:
        ## Create a tmp assembly, call write_params to make default params.txt
        try:
            tmpassembly = MESS.Region(args.new, quiet=True, cli=True)
            tmpassembly.add_local_community("island1", K=1000, c=0.01)
            tmpassembly.write_params("params-{}.txt".format(args.new), 
                                     force=args.force)
        except Exception as inst:
            print(inst)
            sys.exit(2)

        print("\n  New file 'params-{}.txt' created in {}\n".\
               format(args.new, os.path.realpath(os.path.curdir)))
        sys.exit(2)


    ## if params then must provide action argument with it
    if args.params:
        if not any([args.branch, args.results, args.steps]):
            print("""
    Must provide action argument along with -p argument for params file. 
    e.g., MESS -p params-test.txt -r              ## shows results
    e.g., MESS -p params-test.txt -s 12           ## runs steps 1 & 2
    """)
            sys.exit(2)

    if not args.params:
        if any([args.branch, args.results, args.steps]):
            print("""
    Must provide params file for doing steps, or getting results.
    e.g., MESS -p params-test.txt -r              ## shows results
    e.g., MESS -p params-test.txt -s 12           ## runs steps 1 & 2
    """)

    ## always print the header when doing steps
    header = \
    "\n -------------------------------------------------------------"+\
    "\n  MESS [v.{}]".format(MESS.__version__)+\
    "\n  Massive Eco-Evolutionary Synthesis Simulations"+\
    "\n -------------------------------------------------------------"

    ## Log the current version. End run around the LOGGER
    ## so it'll always print regardless of log level.
    with open(MESS.__debugfile__, 'a') as logfile:
        logfile.write(header)
        logfile.write("\n  Begin run: {}".format(time.strftime("%Y-%m-%d %H:%M")))
        logfile.write("\n  Using args {}".format(vars(args)))
        logfile.write("\n  Platform info: {}".format(os.uname()))
    sys.exit()


    ## create new Region or load existing Region, quit if args.results
    if args.params:
        parsedict = parse_params(args)

        if args.branch:
            branch_assembly(args, parsedict)

        elif args.steps:
            ## print header
            print(header)

            ## Only blank the log file if we're actually going to run a new
            ## assembly. This used to be in __init__, but had the side effect
            ## of occasionally blanking the log file in an undesirable fashion
            ## for instance if you run a long assembly and it crashes and
            ## then you run `-r` and it blanks the log, it's crazymaking.
            if os.path.exists(MESS.__debugfile__):
                if os.path.getsize(MESS.__debugfile__) > 50000000:
                    with open(MESS.__debugfile__, 'w') as clear:
                        clear.write("file reset")

            ## run Region steps
            ## launch or load assembly with custom profile/pid
            data = getassembly(args, parsedict)

            ## set CLI ipcluster terms
            data._ipcluster["threads"] = args.threads

            ## if ipyclient is running (and matched profile) then use that one
            if args.ipcluster:
                ipyclient = ipp.Client(profile=args.ipcluster)
                data._ipcluster["cores"] = len(ipyclient)

            ## if not then we need to register and launch an ipcluster instance
            else:
                ## set CLI ipcluster terms
                ipyclient = None
                data._ipcluster["cores"] = args.cores if args.cores else detect_cpus()
                data._ipcluster["engines"] = "Local"
                if args.MPI:
                    data._ipcluster["engines"] = "MPI"
                    if not args.cores:
                        raise IPyradWarningExit("must provide -c argument with --MPI")
                ## register to have a cluster-id with "ip- name"
                data = register_ipcluster(data)

            ## set to print headers
            data._headers = 1

            ## run assembly steps
            steps = list(args.steps)
            data.run(
                steps=steps, 
                force=args.force, 
                preview=args.preview, 
                show_cluster=1, 
                ipyclient=ipyclient)
                     
        if args.results:
            showstats(parsedict)


__PARAMSFILE_TPL="""## MESS Model v.0.0.1 Parameters file ##"""


if __name__ == "__main__": 
    main()
