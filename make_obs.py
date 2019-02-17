#!/usr/bin/env python2.7

import argparse
import glob
import numpy as np
import pandas as pd
import os
from itertools import combinations
from collections import Counter
from MESS.stats import shannon


def pi(file):
    ## Calculate average pi
    pi = 0
    len_seq = 0
    try:
        f = open(file).readlines()
        ## Get just the sequences
        dat = [list(x.strip()) for x in f if ">" not in x]
        len_seq = len(dat[0])

        ## Transpose, so now we have a list of lists of all bases at each
        ## position.
        dat = np.transpose(np.array(dat))

        ## for each position
        for d in dat:
            ## If the position is _not_ monomorphic
            if len(Counter(d)) > 1:
                ## Enumerate the possible comparisons and for each
                ## comparison calculate the number of pairwise differences,
                ## summing over all sites in the sequence.
                base_count = Counter(d)
                ## ignore indels
                del base_count["-"]
                del base_count["N"]
                for c in combinations(base_count.values(), 2):
                    #print(c)
                    n = c[0] + c[1]
                    n_comparisons = float(n) * (n - 1) / 2
                    pi += float(c[0]) * (n-c[0]) / n_comparisons
    except Exception as inst:
        print("Something happenend - {}".format(inst))
        pi = 0
    ## Average over the length of the whole sequence.
    return pi/len_seq


def make_1D_heat(pis):
    max_pi = max(pis)
    np.linspace(0, max(pis), 10)

    pi_island_bins = np.linspace(0, max(pis), 10)
    heat = np.zeros(10, dtype=int)

    for pi_island in pis:
        count_pi_island = 0
        try:
            while not pi_island <= pi_island_bins[count_pi_island]:
                count_pi_island += 1
                ## increment the heatmap point this corresponds to
            heat[count_pi_island] += 1
        except:
            pass
    return("\t".join(map(str, heat)))


def parse_command_line():
    """ Parse CLI args."""

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 


   """)

    ## add model arguments 
    parser.add_argument('-a', dest="abund_file", metavar='abund_file',
        help="File with observed abundances.")

    parser.add_argument('-f', dest="fasta_files", metavar="fasta_files",
        help="Directory containing observed fasta files, one per species.")

    parser.add_argument('-o', dest="outfile", metavar="outfile",
        default="outfile.obs",
        help="Directory containing observed fasta files, one per species.")

    args = parser.parse_args()

    return args

if __name__ == "__main__":

    args = parse_command_line()

    outfile = open(args.outfile, 'w')

    ## Format header
    if args.abund_file:
        outfile.write("shannon\t")
    if args.fasta_files:
        for row in xrange(10):
                outfile.write("\tbin_{}".format(row))
    outfile.write("\n")

    ## Get shannon if abundances are provided
    if args.abund_file:
        dat = open(args.abund_file).read().strip().split(",")
        dat = map(int, dat)
        dat = Counter(dat)
        outfile.write(str(shannon(dat)) + "\t")
        
    ## Get 1D pi vector
    if args.fasta_files:
        files = glob.glob(os.path.join(args.fasta_files, "*.fasta"))

        colname = [args.fasta_files.strip("/").split("/")[-1]]
        pis = {}
        for f in files:
            ## Just use the OTU/species name for the pis file
            pis[f.split("/")[-1].split(".")[0]] = pi(f)
            pis_df = pd.DataFrame.from_dict(pis, orient="index", columns=colname)

        outpis = args.outfile.split(".")[0]+".pis"
        with open(outpis, 'w') as out:
                pis_df.to_csv(out)

        dat = make_1D_heat(pis.values())
        outfile.write(dat)
        

    outfile.close()
