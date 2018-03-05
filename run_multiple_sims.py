#!/usr/bin/env python2.7

from __future__ import print_function

import numpy as np
import subprocess
import time
import sys
import os

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
    SIM_DIRECTORY="simout"
    #RECORDING_INTERVAL = NSIMS/100
    if not os.path.exists(SIM_DIRECTORY):
        os.mkdir(SIM_DIRECTORY)

    for j in xrange(10):
        print("iteration {}".format(j))
        ## Sample from uniform distribution
        col_rates = [0.03]
        k_vals = [5000] 
    
        FNULL = open(os.devnull, 'w')
        for i in xrange(10):
            print("Doing {}".format(i))
            pids = {}
            procs = {}
            for c in col_rates:
                for k in k_vals:
                    for l in xrange(10):
                        t = str(time.time())
                        sim_string = "K_"+str(k)+"-C_"+str(c) + "_" + t + "_" + str(l)
                        cursim_dir = os.path.join(SIM_DIRECTORY, sim_string)
                        cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                            + " -c " + str(c)\
                            + " -k " + str(k)\
                            + " -o " + cursim_dir\
                            + " --model=4 "\
                            + " -q "
                            #+ " > /dev/null &"
                        proc = subprocess.Popen(cmd.split(), shell=False, stdout=FNULL)
                        pid = proc.pid
                        procs[pid] = proc
                        pids[pid] = True
                        time.sleep(1)
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
