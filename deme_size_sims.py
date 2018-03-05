#!/usr/bin/env python2.7

from __future__ import print_function

import subprocess
import time
import os


## To get time to write stderr to an outfile you have to do something like this:
#$ (time ls) > outfile 2>&1
#$ (time ls) > ls_results 2> time_results


if __name__ == "__main__":
    #NSIMS=10000
    NSIMS = 0
    SIM_DIRECTORY="deme_size_sims"
    #RECORDING_INTERVAL = NSIMS/100
    if not os.path.exists(SIM_DIRECTORY):
        os.mkdir(SIM_DIRECTORY)

    col_rates = [0.005, 0.01, 0.03]
    k_vals = [1000, 5000, 10000]

    for i in range(0, 1):
        for c in col_rates:
            for k in k_vals:
                sim_string = "K_"+str(k)+"-C_"+str(c)+"-"+str(time.time())
                cursim_dir = os.path.join(SIM_DIRECTORY, sim_string+"_x1")
                cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                        + " -c " + str(c)\
                        + " -k " + str(k)\
                        + " -o " + cursim_dir\
                        + " &"
                print(cmd)
                subprocess.call(cmd, shell=True)

        ## For each condition do the same condition but with 4 migrants 
        ## at a time and 4x the size of K
        for c in col_rates:
            for k in k_vals:
                sim_string = "K_"+str(k)+"-C_"+str(c)+"-"+str(time.time())
                cursim_dir = os.path.join(SIM_DIRECTORY, sim_string+"_x4")
                cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                        + " -c " + str(c)\
                        + " -k " + str(k*4)\
                        + " -o " + cursim_dir\
                        + " -i 4 "\
                        + " &"
                print(cmd)
                subprocess.call(cmd, shell=True)

        ## Also do clusters of 10 immigrants and 10x size of K 
        for c in col_rates:
            for k in k_vals:
                sim_string = "K_"+str(k)+"-C_"+str(c)+"-"+str(time.time())
                cursim_dir = os.path.join(SIM_DIRECTORY, sim_string+"_x10")
                cmd = "./gimmeSAD.py -n "+ str(NSIMS)\
                        + " -c " + str(c)\
                        + " -k " + str(k*10)\
                        + " -o " + cursim_dir\
                        + " -i 10 "\
                        + " &"
                print(cmd)
                subprocess.call(cmd, shell=True)
