#!/usr/bin/env python2.7

from __future__ import print_function

import numpy as np
import subprocess
import time
import sys
import os


## None of this is true for this file, it's been modified.

## This script will run each combination of col rate and K 100 times.
## So you'll end up with 400 (more) output directories in ./replicates.
##
## On my box the biggest model runs in a little less than 2 minutes, so
## this kind of hackishly runs 4x at a time and sleeps for 2 minutes, then
## starts the next ones.

## To get time to write stderr to an outfile you have to do something like this:
#$ (time ls) > outfile 2>&1
#$ (time ls) > ls_results 2> time_results

if __name__ == "__main__":
    ## 1e6 sims normally just barely gets to equilibrium for most models
    #NSIMS=1000000
    ## Run all the way to double equilibrium
    NSIMS=0
    SIM_DIRECTORY="bottleneck_sims"
    DATA_MODEL=4
    #RECORDING_INTERVAL = NSIMS/100
    if not os.path.exists(SIM_DIRECTORY):
        os.mkdir(SIM_DIRECTORY)

    bot_strengths = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9]
    for j in xrange(100):
        print("iteration {}".format(j))
        FNULL = open(os.devnull, 'w')
        for i in xrange(10):
            print("Doing {}".format(i))
            pids = {}
            procs = {}
            t = str(time.time())
            for c,k in zip([0.03], [10000]):
                for bot in bot_strengths:
                    sim_string = "B-"+str(bot) + "_" + t
                    cursim_dir = os.path.join(SIM_DIRECTORY, sim_string)
                    cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                            + " -c " + str(c)\
                            + " -k " + str(k)\
                            + " -o " + cursim_dir\
                            + " --model={} ".format(DATA_MODEL)\
                            + " -b " + str(bot)\
                            + " -q "
                            #+ " > /dev/null &"
                    proc = subprocess.Popen(cmd.split(), shell=False, stdout=FNULL)
                    pid = proc.pid
                    procs[pid] = proc
                    pids[pid] = True
            done = False
            while not done:
                ## Test all pids have ended
                try:
                    time.sleep(5)
                    for p in pids.keys():
                        try:
                            os.kill(p, 0)
                            procs[p].communicate()
                        except:
                            pids[p] = False
                            ## This one is dead
                    if any(pids.values()):
                        print(pids.keys())
                        print(pids.values())
                    else:
                        done = True
                except:
                    pass
