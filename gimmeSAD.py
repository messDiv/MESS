#!/usr/bin/env python2.7

from __future__ import print_function
from ascii_graph import Pyasciigraph
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
from tabulate import tabulate
from scipy.stats import iqr
import pandas as pd
import numpy as np
import subprocess
import collections
import argparse
import datetime
import shutil
import time
import math
import sys
import os

plt.switch_backend('agg')

import implicit_BI
import implicit_CI

# pylint: disable=C0103
# pylint: disable=R0903

## The fuckin profiling function I always forget:
## python -m cProfile -s cumtime gimmeSAD.py -n 0 -k 1000 -o /tmp/wat -f > res.txt

class gimmeSAD(object):

    def __init__(self):
        pass        

    def __str__(self):
        return "<gimmeSAD {}>".format(self.name)


## quicksort stolen from the internet
def qsort(arr): 
     if len(arr) <= 1:
          return arr
     else:
          return qsort([x for x in arr[1:] if x.abundance<arr[0].abundance])\
                    + [arr[0]] + qsort([x for x in arr[1:] if x.abundance>=arr[0].abundance])

def make_outputfile(model, stats):
    """ Make the output file formatted correctly for each model.
        Model 1 - Local pi vector only
        Model 2 - Local pi vector and observed SAD
        Model 3 - Local pi and dxy
        Model 4 - Local pi, dxy, and observed SAD
        """
    ## Write the header, the same for all.
    stats.write("K\tc\tstep\t%equil\tcolrate\textrate\tshannon")
    ## Write the pi or pi/dxy bins
    if model in [1, 2]:
        for row in xrange(10):
            stats.write("\tbin_{}".format(row))
    elif model in [3, 4]:        
        for row in xrange(10):
            for col in xrange(10):
                stats.write("\tbin_{}_{}".format(row, col))
    else:
        raise Exception("Value for --model must be one of [1,2,3,4]")

    stats.write("\n")  

def get_min_max_stats_through_time(sp_through_time):
    ## Normalization routines
    max_pi_island = 0
    max_dxy = 0

    for sp_list in sp_through_time.values():

        ## Get max pi and max dxy
        pis = np.array([(x.dxy, x.pi_island) for x in sp_list])
        ## pis will be empty if this timeslice includes no extant species
        if pis.any():
            my_max_dxy = max([x[0] for x in pis])
            my_max_pi_island = max([x[1] for x in pis])
            if max_dxy < my_max_dxy:
                max_dxy = my_max_dxy
            if max_pi_island < my_max_pi_island:
                max_pi_island = my_max_pi_island

    return max_pi_island, max_dxy

def normalize_heatmap_to_numpy(sp_list, max_pi, max_dxy):
    ## Get the sumstats for this timeslice
    ## Only include extant species in plots
    #pis = np.array([(x.dxy, x.pi_island) for x in sp_list if x.uuid[0] in extant])
    ## include all species at each timeslice
    pis = np.array([(x.dxy, x.pi_island) for x in sp_list])

    ## Empty heatmap we'll write into
    heat = np.zeros((10,10), dtype=np.int)

    ## Make the bins
    dxy_bins = np.linspace(0, max_dxy, 10, endpoint=True)
    pi_island_bins = np.linspace(0, max_pi, 10, endpoint=True)

    ## Now you have the bins each value belongs in, but you need to 
    ## go through and populate the heat matrix
    for dxy, pi_island in pis:
        count_dxy = 0
        count_pi_island = 0
        try:
            while not dxy <= dxy_bins[count_dxy]:
                count_dxy += 1
            while not pi_island <= pi_island_bins[count_pi_island]:
                count_pi_island += 1
            ## increment the heatmap point this corresponds to
            heat[count_dxy][count_pi_island] += 1
        except Exception as inst:
            ## Got a value bigger than our current max pi/dxy. ignore.
            pass
    return heat
        

## This is hackish because it goes through and writes out some stats at runtime
## then goes back at the very end and writes out the heatmaps per line
## of the file. very confusing.
def write_heats_to_outfile(model, stats, data, sp_through_time, equilibria):
    stats.close()
    stats = open(stats.name, 'r')
    lines = stats.readlines()

    max_pi, max_dxy = get_min_max_stats_through_time(sp_through_time)

    for i, line in enumerate(lines[1:]):
        time = int(line.strip().split()[2])
        heat = normalize_heatmap_to_numpy(sp_through_time[time], max_pi, max_dxy)

        ## Models 1 & 2 use the 1D island pi only vector
        if model in [1, 2]:
        ## This is a 2d heatmap so we have to marginalize over dxy
            heat = heat.sum(axis=0)
        ## Models 3 & 4 use the 2D pi x dxy matrix
        elif model in [3, 4]:
            heat = heat.flatten()
        lines[i+1] = line.strip() + "\t" + "\t".join(map(str,heat))

    stats.close()
    lines[0] = lines[0].strip()
    stats = open(stats.name, 'w')
    stats.write("\n".join(lines))

def write_outfile(model, stats, data, eq):
    ## Calculate some crap of interest
    extrate = data.extinctions/float(data.current_time)
    colrate = data.colonizations/float(data.current_time)
    shan = shannon(data.get_abundances(octaves=False))

    ## Write the common data
    stats.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(data.local_inds, data.colonization_rate,\
                                                      data.current_time, eq, colrate, extrate, shan))
    stats.write("\n")
    return
    heat = heatmap_pi_dxy_ascii(data, labels=False).strip()
    heat = np.array([np.array(x.split(" "), dtype=int) for x in heat.split("\n")])

    ## Models 1 & 2 use the 1D island pi only vector
    if model in [1, 2]:
        ## This is a 2d heatmap so we have to marginalize over dxy
        heat = heat.sum(axis=0)
    ## Models 3 & 4 use the 2D pi x dxy matrix
    elif model in [3, 4]:
        heat = heat.flatten()

    stats.write("\t".join(map(str,heat)))
    stats.write("\n")

## Here abundances is an ordered dict of tuples which are (abundance, count)
def shannon(abundances):
    ## Unpack the abundance dist
    abunds = [v for v in abundances.values()]
    tot = np.sum(abunds)
    return -1 * np.sum([x/float(tot) * math.log(x/float(tot)) for x in abunds  if x > 0])


def abundances_from_sp_list(species, octaves=False):
    ## Make a counter for the local_community, counts the number of
    ## individuals w/in each species
    abundances = collections.Counter([x.abundance for x in species])
    ## Get rid of any 0 abundances cuz it fsck us up.
    try:
        abundances.pop(0)
    except KeyError:
        pass

    abundance_distribution = collections.OrderedDict()
    all_abund = abundances.keys()
    max_abund = max(all_abund)
    for ab_class in np.arange(1, max_abund+1):
        if ab_class in all_abund:
            abundance_distribution[ab_class] = abundances[ab_class]
        else:
            abundance_distribution[ab_class] = 0

    if octaves:
        dist_in_octaves = collections.OrderedDict()
        min_oct = 1
        max_oct = 2
        while max_oct/2 < len(abundance_distribution):
            count = 0
            ## Here `i` is the abundance class and
            ## `j` is the count for that class
            for i, j in abundance_distribution.items():
                if (i < max_oct) and (i >= min_oct):
                    count += j
            dist_in_octaves[min_oct] = count
            min_oct = min_oct * 2
            max_oct = max_oct * 2
        abundance_distribution = dist_in_octaves
    return abundance_distribution


def plot_abundances_ascii(abundance_distribution):
    """
    Plot the abundance distribution as returned by the get_abundances
    method. The abundance dist in this case is a dict of abundance
    classes and counts
    """
    ## Cool ascii bar graph
    ## The bar grapher buddy doesn't like ints for description so you have to transform it
    abundance_distribution = [(str(k), v) for k, v in abundance_distribution.items()]
    graph = Pyasciigraph(graphsymbol="|")
    msg = "\n\n"
    for line in  graph.graph('Simulated Species Abundance Distribution', abundance_distribution):
        msg = msg + line + "\n"
    msg = msg + "###############################################################################\n"
    return msg


def plot_abundances_gui(abundance_distribution):
    plt.bar(abundance_distribution.keys(), abundance_distribution.values())
    plt.show()


def tabulate_sumstats(data):
    sp = data.get_species()
    ## Highlight the invasive
    for i,s in enumerate(sp):
        try:
            if s.uuid == data.invasive:
                s.name = s.name + "*"
                sp[i] = s
        except:
            pass
            ## This will barf if invasive isn't a tuple

    #print("Species colonization times (in generations):\n{}".format([x.colonization_time for x in sp]))
    #print("Species Ne:\n{}".format([x.Ne for x in sp]))
    headers = ["Species Name", "Col time", "Loc Abund", "Meta Abund", "pi", "pi_net", "Dxy",  "S", "S_island", "pi_island", "tajD_island", "S_meta", "pi_meta"]
    acc = [[s.name, s.colonization_time, s.abundance, int(s.meta_abundance), s.pi, s.pi_net, s.dxy, s.S, s.S_island, s.pi_island, s.tajD, s.S_meta, s.pi_meta] for s in sp]

    return tabulate(acc, headers, floatfmt=".4f")


def write_megalog(megalogfile, i, percent_equil, data):
    sp = data.get_species()
    acc = [[s.uuid[0], s.colonization_time, s.abundance, int(s.meta_abundance), s.pi, s.pi_net, s.dxy, s.S, s.S_island, s.pi_island, s.tajD, s.S_meta, s.pi_meta] for s in sp]
    megalogfile.write("\n".join(["{}\t{}\t{}".format(percent_equil, i, "\t".join(map(str, s))) for s in acc]) + "\n")

## This actually is doing pi x dxy, but some of the variable
## names are goofy cuz i developed it for pi x pi_w_island 
def heatmap_pi_dxy(data, write="", title=""):
    """ The write flag specifies whether or not to write the image
    to a file. You must pass in a file name. If this file exists
    it gets overwritten. Normally write should be a full path.
    No need to pass in the file extension."""
    sp = data.get_species()

    heat = np.zeros((10,10), dtype=np.int)

    pis = np.array([(x.dxy, x.pi_island) for x in sp])
    max_pi = max([x[0] for x in pis])
    max_pi_island = max([x[1] for x in pis])

    ## Make the bins
    pi_bins = np.linspace(0, max_pi, 10)
    pi_island_bins = np.linspace(0, max_pi_island, 10)

    ## Now you have the bins each value belongs in, but you need to 
    ## go through and populate the heat matrix
    for pi, pi_island in pis:
        count_pi = 0
        count_pi_island = 0
        while not pi <= pi_bins[count_pi]:
            count_pi += 1
        while not pi_island <= pi_island_bins[count_pi_island]:
            count_pi_island += 1
        ## increment the heatmap point this corresponds to
        heat[count_pi][count_pi_island] += 1
    plt.pcolormesh(heat,cmap=plt.cm.Blues)
    plt.xlabel('Dxy')
    plt.ylabel('Pi_w Island')
    plt.colorbar()
    plt.xticks(np.arange(len(pi_bins)), ["{0:.4f}".format(x) for x in pi_bins], rotation='vertical')
    plt.yticks(np.arange(len(pi_bins)), ["{0:.4f}".format(x) for x in pi_island_bins])
    
    ## If writing to a file, don't bother displaying it, plus it hangs the program
    if write:
        plt.title(title)
        plt.savefig(write+".png")
    else:
        plt.show()
    plt.close()


## Doesn't exactly work right, it'll prune extant species but there's 
## a good chance an extant species has colonized, then gone extinct,
## then recolonized (which is allowable), but it fucks up the plotting
def prune_extant(sp_through_time):
    """ This is fucked because at each recording period the species
    get new names generated for them, even though the UUID is the same """

    ## Get the uuids of all extant species
    extant = [x.uuid[0] for x in sp_through_time.values()[-1]]
    
    for time, sp_list in sp_through_time.items(): 
        sp_through_time[time] = np.array([x for x in sp_list if x.uuid[0] in extant])

    return sp_through_time


def plot_rank_abundance_through_time(outdir, sp_through_time, equilibria,\
                                only_extant=False, stats_models=False, as_curve=False,\
                                verbose=False):
    import seaborn
    seaborn.set
    seaborn.set_style(style="white")

    print("Generating abundance distributions through time")
    ## Make the output directory for heatmaps inside the top level output directory
    abund_out = os.path.join(outdir, "abundance_plots")
    if not os.path.exists(abund_out):
        os.mkdir(abund_out)

    times = equilibria.keys()
    equilibria = equilibria.values()

    ## If you only want to see extant species then prune all the extinct ones
    if only_extant:
        sp_through_time = prune_extant(sp_through_time)

    max_n_species, max_abundance, max_octave, max_class_count, max_n_bins, octave_bin_labels = prep_normalized_plots(sp_through_time)
    if verbose:
        print("info:\n\nmax_n_species - {}\nmax_abundance - {}\nmax_octave - {}\nmax_class_count - {}"\
                + "\nmax_n_bins - {}\noctave_bin_labels - {}\n".format(\
                max_n_species, max_abundance, max_octave, max_class_count,\
                max_n_bins, octave_bin_labels))

    ## Get a list of UUIDs of the extant species at time 0 (present)
    extant = [x.uuid[0] for x in sp_through_time.values()[-1]]
    
    ## Make the abundance plots, one for each timeslice
    ## TODO: Should we include all species or just the ones that live to the end?
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the abundance
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))

    tot_plots = len(sp_through_time.values())
    for i, species in enumerate(sp_through_time.values()):
        progressbar(tot_plots, i+1)

        title = "Time_"+str(file_index+i)
        write = os.path.join(abund_out, title)

        ## Make the Plot
        fig = plt.figure(figsize=(12,5))

        ## Make the SAD subplot
        abund = abundances_from_sp_list(species, octaves=True)
        ax1 = plt.subplot(121)
        plot_sad(abund, max_n_species, max_n_bins, max_class_count, octave_bin_labels, verbose)

        ## Make the rank abundance distribution subplot
        plt.subplot(122)
        plot_rank_abundance(species, max_n_species, max_abundance, stats_models, as_curve)
        
        plt.subplots_adjust(bottom=0.15)

        suptitle = "Time: {}".format(times[i])
        fig.get_axes()[0].annotate(suptitle, (0.5, 0.95), 
                            xycoords='figure fraction', ha='center', 
                            fontsize=24
                            )

        plt.savefig(write+".png")
        plt.close()

    progressbar(100, 100, "\n")
    make_animated_gif(abund_out,\
                        os.path.join(outdir, "abundances_through_time.gif"))


def plot_abundance_vs_colonization_time(outdir, sp_through_time, equilibria,\
                                stats_models, as_curve, only_extant=False,\
                                verbose=False):
    import seaborn
    seaborn.set
    seaborn.set_style(style="white")

    print("Generating abundance vs colonization plots through time")
    ## Make the output directory for heatmaps inside the top level output directory
    abund_out = os.path.join(outdir, "abundance_vs_coloniz")
    if not os.path.exists(abund_out):
        os.mkdir(abund_out)

    times = equilibria.keys()
    equilibria = equilibria.values()

    ## If you only want to see extant species then prune all the extinct ones
    if only_extant:
        sp_through_time = prune_extant(sp_through_time)

    max_n_species, max_abundance, max_octave, max_class_count, max_n_bins, octave_bin_labels = prep_normalized_plots(sp_through_time)
    if verbose:
        print("info:\n\nmax_n_species - {}\nmax_abundance - {}\nmax_octave - {}\nmax_class_count - {}"\
                + "\nmax_n_bins - {}\noctave_bin_labels - {}\n".format(\
                max_n_species, max_abundance, max_octave, max_class_count,\
                max_n_bins, octave_bin_labels))

    ## Get a list of UUIDs of the extant species at time 0 (present)
    extant = [x.uuid[0] for x in sp_through_time.values()[-1]]
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the abundance
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))

    tot_plots = len(sp_through_time.values())
    for i, species in enumerate(sp_through_time.values()):
        progressbar(tot_plots, i+1)

        title = "Time_"+str(file_index+i)
        write = os.path.join(abund_out, title)

        ## Make the Plot
        fig = plt.figure(figsize=(12,5))

        ## Make the abundance vs colonization time plot
        ax1 = plt.subplot(121)
        max_coltime = max([s.colonization_time for s in species])
        plot_abund_vs_colon(species, max_coltime, max_abundance)

        ## Make the rank abundance distribution subplot
        plt.subplot(122)
        plot_rank_abundance(species, max_n_species, max_abundance, stats_models, as_curve)

        plt.subplots_adjust(bottom=0.15)

        suptitle = "Time: {}".format(times[i])
        fig.get_axes()[0].annotate(suptitle, (0.5, 0.95),
                            xycoords='figure fraction', ha='center',
                            fontsize=24
                            )

        plt.savefig(write+".png")
        plt.close()

    progressbar(100, 100, "\n")
    make_animated_gif(abund_out,\
                        os.path.join(outdir, "abundances_vs_coltime.gif"))


def plot_abund_vs_colon(species, max_coltime, max_abundance):
    x = [np.log10(s.colonization_time) for s in species]
    y = [np.log10(s.abundance) for s in species]
    plt.scatter(x, y, color="blue", s=100)
    plt.ylim(0, int(math.ceil(np.log10(max_abundance))))
    plt.xlim(0, 8) #int(math.ceil(np.log10(max_coltime))))
    plt.title("Abundance vs Colonization Time", fontsize=24)
    plt.ylabel("Abundance (log10)", fontsize=20)
    plt.xlabel("Colonization Time (log10)", fontsize=20)


def prep_normalized_plots(sp_through_time):
    ## GET MAX values for abundance and num species so we can normalize the plot axes
    max_n_species = max([len(x) for x in sp_through_time.values()])
    max_abundance = max([max([y.abundance for y in sp]) for sp in sp_through_time.values()])

    ## Get max values for abundance class count and abundance octave
    max_octave = 0
    max_class_count = 0
    max_n_bins = 0
    octave_bin_labels = []
    for sp in sp_through_time.values():
        abund = abundances_from_sp_list(sp, octaves=True)
        octave = max(abund.keys())
        class_count = max(abund.values())
        if octave > max_octave:
            max_octave = octave
        if class_count > max_class_count:
            max_class_count = class_count
        if len(abund.keys()) > max_n_bins:
            max_n_bins = len(abund.keys())
            octave_bin_labels = abund.keys()

    return max_n_species, max_abundance, max_octave, max_class_count, max_n_bins, octave_bin_labels


def make_animated_gif(datadir, outfile):
    """ This function will take all png files in a directory and make them
    into an animated gif. The inputs are the directory with all the images
    and the full path including filename of the file to write out"""

    ## Do the imagemagick conversion, if possible
    ## `convert -delay 100 outdir/* anim.gif`
    ## Define the total time to be 10 seconds, total available
    ## timeslots is * 100 bcz convert flag is in 1/100 seconds
    tot_time = 10 * 100
    if i > tot_time:
        delay = 1
    else:
        delay = int(1000./i)
    ## Default half second intervals
    delay = 50
    cmd = "convert -delay {} ".format(delay)\
            + datadir + "/*.png "\
            + outfile
    try:
        subprocess.Popen(cmd.split())
    except Exception as inst:
        print("Trouble creating abundances through time animated gif - {}".format(inst))
        print("You probably don't have imagemagick installed")


def get_max_heat_bin(sp_through_time, max_pi_island, max_dxy):
    max_heat_bin = 0

    for sp_list in sp_through_time.values():
        ## Get the sumstats for this timeslice
        ## Only include extant species in plots
        #pis = np.array([(x.dxy, x.pi_island) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.dxy, x.pi_island) for x in sp_list])

        ## Empty heatmap we'll write into
        heat = np.zeros((20,20), dtype=np.int)

        ## Make the bins
        dxy_bins = np.linspace(0, max_dxy, 20, endpoint=True)
        pi_island_bins = np.linspace(0, max_pi_island, 20, endpoint=True)

        ## Now you have the bins each value belongs in, but you need to 
        ## go through and populate the heat matrix
        for dxy, pi_island in pis:
            count_dxy = 0
            count_pi_island = 0
            try:
                while not dxy <= dxy_bins[count_dxy]:
                    count_dxy += 1
                while not pi_island <= pi_island_bins[count_pi_island]:
                    count_pi_island += 1
                ## increment the heatmap point this corresponds to
                heat[count_dxy][count_pi_island] += 1
            except Exception as inst:
                ## Got a value bigger than our current max pi/dxy. ignore.
                pass
        if np.amax(heat) > max_heat_bin:
            max_heat_bin = np.amax(heat)

    return max_heat_bin


def normalized_pi_dxy_heatmaps(outdir, sp_through_time, equilibria, one_d=False,\
                                only_extant=False, stats_models=False, as_curve=False,\
                                verbose=False):
    """ Normalize x and y axes for the heatmaps. Only take into account extant species.
    Inputs are the output directory to write to and an ordered dict 
    of the species at every recording duration timepoint """

    if one_d:
        msg = "Generating 1-D pi_w heatmap animation"
        datadir = "1D_pi_heatmaps"
        outfile = "1D-pi_anim.gif"
    else:
        msg = "Generating pi x dxy heatmap animation"
        datadir = "normalized_heatmaps"
        outfile = "pi_dxy_anim.gif"

    print(msg)
    ## Make the output directory for heatmaps inside the top level output directory
    heat_out = os.path.join(outdir, datadir)
    if not os.path.exists(heat_out):
        os.mkdir(heat_out)

    times = equilibria.keys()
    equilibria = equilibria.values()

    ## If you only want to see extant species then prune all the extinct ones
    if only_extant:
        sp_through_time = prune_extant(sp_through_time)

    ## GET MAX pi and dxy
    ## Get a list of UUIDs of the extant species at time 0 (present)
    extant = [x.uuid[0] for x in sp_through_time.values()[-1]]
    ## find the max_pi and max_dxy for all extant species through all timepoints
    max_pi_island = 0
    max_dxy = 0
    
    ## For each recorded timeslice
    my_dxys = []
    my_pi_islands = []

    ## Get variables we care about
    max_n_species, max_abundance, _, _, _, _ = prep_normalized_plots(sp_through_time)
    if verbose:
        print("max_n_species - {}\nmax_abundance - {}".format(max_n_species, max_abundance))
    
    ## Normalization routines
    for sp_list in sp_through_time.values():

        ## Get max pi and max dxy
        pis = np.array([(x.dxy, x.pi_island) for x in sp_list if x.uuid[0] in extant])
        ## pis will be empty if this timeslice includes no extant species
        if pis.any():
            my_max_dxy = max([x[0] for x in pis])
            my_max_pi_island = max([x[1] for x in pis])
            if max_dxy < my_max_dxy:
                max_dxy = my_max_dxy
            if max_pi_island < my_max_pi_island:
                max_pi_island = my_max_pi_island
            my_dxys.append(my_max_dxy)
            my_pi_islands.append(my_max_pi_island)

#    max_dxy = np.average(my_dxys)
#    max_pi_island = np.average(my_pi_islands)
    max_dxy = np.median(my_dxys)
    max_pi_island = np.median(my_pi_islands)
    ## However this function is calculating the heatmaps is fucked up, so you have to
    ## hard code max values for pi and dxy here.
    max_dxy = 0.04
    max_pi_island = 0.02
    if verbose:
        print("Got\tmax_dxy - {}\t max_pi_island - {}".format(max_dxy, max_pi_island))
    ## Make the heatmaps, one for each timeslice
    ## TODO: Should we include all species or just the ones that live to the end?
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the heatmap
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))

    max_bin_val = get_max_heat_bin(sp_through_time, max_pi_island, max_dxy)
    ## Make the throwaway plot to get the colorbar
    cbar_min, cbar_max = (0, max_bin_val)
    step = 1
    Z = [[0,0],[0,0]]
    levels = range(cbar_min, cbar_max+step, step)
    my_colorbar = plt.contourf(Z, levels, cmap=plt.cm.jet)
    plt.clf()


    tot_heatmaps = len(sp_through_time.values())
    for i, sp_list in enumerate(sp_through_time.values()):
        progressbar(tot_heatmaps, i+1)
 
        title = "Time_"+str(file_index+i)
        write = os.path.join(heat_out, title)
        #print("Doing", title)
        
        ## Get the sumstats for this timeslice
        ## Only include extant species in plots
        #pis = np.array([(x.dxy, x.pi_island) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.dxy, x.pi_island) for x in sp_list])

        ## Empty heatmap we'll write into
        heat = np.zeros((20,20), dtype=np.int)

        ## Make the bins
        dxy_bins = np.linspace(0, max_dxy, 20, endpoint=True)
        pi_island_bins = np.linspace(0, max_pi_island, 20, endpoint=True)

        ## Now you have the bins each value belongs in, but you need to 
        ## go through and populate the heat matrix
        for dxy, pi_island in pis:
            count_dxy = 0
            count_pi_island = 0
            try:
                while not dxy <= dxy_bins[count_dxy]:
                    count_dxy += 1
                while not pi_island <= pi_island_bins[count_pi_island]:
                    count_pi_island += 1
                ## increment the heatmap point this corresponds to
                heat[count_dxy][count_pi_island] += 1
            except Exception as inst:
                ## Got a value bigger than our current max pi/dxy. ignore.
                pass

        ## Make the Plot
        fig = plt.figure(figsize=(12,5))
        ## Title plot the % of equilibrium
        #plt.suptitle("%equilibrium = {}".format(equilibria[i]), fontsize=25)
        ## Title plot as time
        #plt.suptitle("Time\n{}".format(times[i]), fontsize=25)
        ax1 = plt.subplot(121)

        ## Trick to treat 0 specially in the colorplot (set it to 'white')
        cmap = plt.cm.jet
        cmap.set_under('white')
        eps = np.spacing(0.1)

        if one_d:
            ## Make the 1D pi_within plot
            one_d_array = np.sum(heat, axis=0)
            one_d_array = np.array([one_d_array])
            if verbose:
                print("1D array - {}".format(one_d_array))

            ## Set the axis so it won't stretch the heatmap
            ax1.set_aspect('equal')
            heat = one_d_array
        else:
            if verbose:
                print("Heatmap array - {}".format(heat))

        plt.pcolormesh(heat,cmap=cmap, vmin=eps)

        plt.colorbar(my_colorbar)

        plt.xlabel(u"Nucleotide diversity (\u03c0)", fontsize=20)
        plt.xticks(np.arange(len(pi_island_bins)), ["{0:.4f}".format(x) for x in pi_island_bins], rotation='vertical')

        if one_d:
            pass
        else:
            pass
        ## This is cool, but it makes the plots really busy
        ## Add numbers to the heatplots
        #    for y in range(heat.shape[0]):
        #        for x in range(heat.shape[1]):
        #            if heat[y, x] < 1:
        #                val = ""
        #            else:
        #                val = str(int(heat[y,x]))
        #            plt.text(x + 0.5, y + 0.5, '%s' % val,
        #                     horizontalalignment='center',
        #                     verticalalignment='center',
        #                     path_effects=[PathEffects.withStroke(linewidth=3,foreground="w")])

#        Trying to adjust the stupid x-axis labels
#        ax1.set_xticks(np.arange(len(pi_island_bins) + 0.5))
#        ax1.set_xticklabels(["{0:.4f}".format(x) for x in pi_island_bins], rotation="vertical", ha="center")

        if one_d:
            plt.title("1D-SGD", fontsize=24)
            plt.yticks([])
        else:
            plt.title("2D-SGD", fontsize=24)
            plt.yticks(np.arange(len(dxy_bins)), ["{0:.4f}".format(x) for x in dxy_bins])
            #plt.ylabel('Pairwise differences between \nisland and metacommunity (Dxy)', fontsize=20)
            plt.ylabel(r"Absolute divergence ($D_{xy}$)", fontsize=20)

        ## Pad margins so labels don't get clipped
        plt.subplots_adjust(bottom=0.25)
        #plt.title(title)

        ## Make the rank abundance plot
        plt.subplot(122)
        plot_rank_abundance(sp_list, max_n_species, max_abundance, stats_models, as_curve)
    
        ## Make the super title
        suptitle = "Time: {}".format(times[i])
        fig.get_axes()[0].annotate(suptitle, (0.5, 0.95), 
                            xycoords='figure fraction', ha='center', 
                            fontsize=24
                            )
        plt.savefig(write+".png")
        plt.close()

    progressbar(100, 100, "\n")

    make_animated_gif(heat_out,\
                        os.path.join(outdir, outfile))


def plot_sad(abund, max_n_species, max_n_bins, max_class_count, octave_bin_labels, verbose):
    ax1 = plt.gca()
    ab_class, count = zip(*abund.items())
    if verbose:
        print(ab_class, count)
    df = pd.DataFrame([x for x in count], index = [str(x) for x in ab_class])
    i = -1
    ## This is hax to make seaborn print out the x-axis labels
    ## If the len of the dataframe is < the max len of any df
    ## we print, even though we're trying to set the octave_bin_labels
    ## it just won't print them. It just looks weird w/o this.
    while len(df) < max_n_bins:
        df.loc[i] = 0
        i -= 1
    bar = df.plot(kind = "bar", legend = False, ax = ax1)
    ## Otherwise 1st bar is truncated
    plt.title("SAD", fontsize=24)
    plt.xlim([0.5, max_n_bins]) 
    ax1.set_xticklabels([str(x) for x in octave_bin_labels])
    plt.setp(bar.get_xticklabels(), rotation=0)
    plt.ylim(0, max_class_count)
    plt.xlabel("Abundance Class", fontsize=20)
    plt.ylabel("Count", fontsize=20)
        

def plot_rank_abundance(sp_list, max_n_species, max_abundance, stats_models=False, as_curve=False):

    species = qsort(sp_list)
    species = species[::-1]
    X = np.arange(0,len(species))
    if as_curve:
        Y = [xx.abundance for xx in species]
        plt.semilogy(X, Y, label="simulated")
        ymax = max_abundance
    else:
        Y = [np.log10(xx.abundance) for xx in species]
        plt.scatter(X, Y, color="blue", s=100, label="simulated")
        ymax = int(math.ceil(np.log10(max_abundance)))
        
    plt.title("Rank Abundance", fontsize=24)
    plt.xlim(0, max_n_species)
    plt.ylim(0, ymax)
    plt.ylabel("Abundance (log10)", fontsize=20)
    plt.xlabel("Rank", fontsize=20)

    ## Whether or not to include a couple common statistical models in the plots
    if stats_models:
        import macroeco as meco
        abund = [xx.abundance for xx in species]
        ## Lognormal
        mu, s = meco.models.lognorm.fit_mle(abund)
        lognorm_rad = meco.models.lognorm.rank(len(abund), mu, s)
        if as_curve:
            Y = lognorm_rad[::-1]
            plt.semilogy(X, Y, label="Lognorm RAD")
        else:
            Y = [int(math.ceil(np.log10(x))) for x in lognorm_rad[::-1]]
            plt.scatter(X, Y, s=100, color="green", label="Lognorm RAD")
        ## Logseries
        p = meco.models.logser.fit_mle(abund)
        logser_rad = meco.models.logser.rank(len(abund), p)
        if as_curve:
            Y = logser_rad[::-1]
            plt.semilogy(X, Y, label="Logseries RAD")
        else:
            Y = [int(math.ceil(np.log10(x))) for x in logser_rad[::-1]]
            plt.scatter(X, Y, s=100, color="red", label="Logseries RAD")
        ## Poisson Lognormal
        mu, s = meco.models.plnorm_ztrunc.fit_mle(abund)
        plnorm_rad = meco.models.plnorm_ztrunc.rank(len(abund), mu, s)
        if as_curve:
            Y = plnorm_rad[::-1]
            plt.semilogy(X, Y, label="Logseries RAD")
        else:
            Y = [int(math.ceil(np.log10(x))) for x in plnorm_rad[::-1]]
            plt.scatter(X, Y, s=100, color="red", label="Poisson Lognorm RAD")

        plt.legend()


def heatmap_pi_dxy_ascii(data, labels=False):
    """ This is kind of a toy. It logs to a file in the output directory
    at recording_period interval, but the axes aren't normalized, so it's
    more semi-informational. """
    sp = data.get_species()

    heat = np.zeros((10,10), dtype=np.int)

    pis = np.array([(x.dxy, x.pi_island) for x in sp])
    ## Set a reasonable default
    max_pi = max_pi_island = 0.1
    if pis.any():
        max_pi = max([x[0] for x in pis])
        max_pi_island = max([x[1] for x in pis])
        #print(max_pi, max_pi_island)

    ## Make the bins
    pi_bins = np.linspace(0, max_pi, 10)
    pi_island_bins = np.linspace(0, max_pi_island, 10)

    ## Now you have the bins each value belongs in, but you need to 
    ## go through and populate the heat matrix
    for pi, pi_island in pis:
        count_pi = 0
        count_pi_island = 0
        while not pi <= pi_bins[count_pi]:
            count_pi += 1
        while not pi_island <= pi_island_bins[count_pi_island]:
            count_pi_island += 1
        ## increment the heatmap point this corresponds to
        heat[count_pi][count_pi_island] += 1
    pivals = [x[0] for x in pis]
    dxyvals = [x[1] for x in pis]

    ret = ""
    if labels:
        ret = "mean/stdv/median/iqr pi {}/{}/{}/{}\tmean/stdv/median/iqr dxy {}/{}/{}/{}".format(\
                                            np.mean(pivals), np.std(pivals), np.median(pivals), iqr(pivals),\
                                                            np.mean(dxyvals), np.std(dxyvals), np.median(dxyvals), iqr(dxyvals))
    ## ascii format the data, and make a weak attempt to convey some information
    if labels:
        ret += "\nDxy\n"

    for i, row in enumerate(heat):
        if labels:
            ret += "{0:.4f}".format(pi_bins[i])
        ret += ' '.join(map(str, row)) + "\n"

    if labels:
        ret += "\t" + str(["{0:.4f}".format(x) for x in pi_island_bins])
        ret += "\n\t\tpi_w island"

    return ret


def write_sizechanges(outdir, yoyo):
    """ Write out sizechanges through time for all species extant
    at the end of the simulation. Takes in an ordered dict of collections """

    out = open(os.path.join(outdir, "pop_yoyos.txt"), 'w')
    keys = yoyo.keys()
    keys.reverse()
    ## Get a list of all extant species, the only ones we care about pop size change in.
    extant = yoyo[keys[0]].keys()
    sizes = {}
    for x in extant:
        sizes[x] = []
    ## For each time slice going back in time, record the current pop size of each 
    ## species of interest
    for k in keys:
       for sp, siz in yoyo[k].items():
            if sp in extant:
                sizes[sp].append(siz)
    for k, v in sizes.items():
        out.write("{} - {}\n".format(k, v))
    

def parse_command_line():
    """ Parse CLI args. Only three options now. """

    ## create the parser
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\n
  * Example command-line usage: 

    gimmeSAD -n 100000                  ## Run the simulation for 100000 steps

    gimmeSAD -o outdir                  ## Directory where tons of junk will get written
    gimmeSAD -r recording_period        ## Num steps to sample and write out after

    gimmeSAD -m metacommunity_LS4.txt   ## Source metacommunity from a file of ints
                                        ## the other option is to pass "uniform"
                                        ## which will generate a random metacommunity
                                        ## with all species at uniform frequency

    gimmeSAD -p volcanic                ## Set the mode for prepopulating (or not)
                                        ## the local community.

    gimmeSAD -i 4                       ## Set # of colonizing inds per colonization event
                                        ## if unset basic immigration is assumed (1 individual)
    gimmeSAD -c 0.001                   ## Set colonization rate

    gimmeSAD -a                         ## Output abundances in octaves
    gimmeSAD -q                         ## Don't be so chatty
    gimmeSAD -v                         ## Be very chatty

   """)

    ## add model arguments 
    parser.add_argument('-a', dest="octaves", action='store_true',
        help="Print species abundances in octaves")

    parser.add_argument('-c', metavar='colrate', dest="colrate", type=float,
        default=0.003,
        help="Set colonization rate")

    parser.add_argument('-C', metavar='colonizers', dest="colonizers", type=int,
        default=0,
        help="Switch mode to clustered colonization and set the # of colonizers per event")

    parser.add_argument('-k', metavar='K', dest="K", type=int,
        default=10000,
        help="Carrying capacity of the island (max # individuals in the local community")

    parser.add_argument('-m', metavar='meta', dest="meta", type=str,
        default="logser",
        help="Source metacommunity from file or generate uniform")

    parser.add_argument('-n', metavar='nsims', dest="nsims", type=int,
        default=50000,
        help="Number of demographic events to simulate")

    parser.add_argument('-o', metavar='outdir', dest="outdir", type=str,
        default="output",
        help="Output directory for log/data files and pngs")

    parser.add_argument("--model", metavar="model", dest="model", type=int,
        default=4,
        help="Data model for writing output files.")

    parser.add_argument('-p', metavar='mode', dest="mode", type=str,
        default="volcanic",
        help="Select mode for prepopulating the island (volcanic/landbridge)")

    parser.add_argument('-r', metavar='recording_period', dest="recording_period", type=int,
        default=0,
        help="Length of timeslice between samples for logging")

    parser.add_argument('-i', metavar='invasion_time', dest="invasion_time", type=int,
        default=-1,
        help="Timestep of invasion of invasive species")

    parser.add_argument('-I', metavar='invasiveness', dest="invasiveness", type=float,
        default=0.1,
        help="Invasiveness of the invasive species")

    parser.add_argument('-e', dest="exponential", action='store_true',
        help="Do exponential growth, otherwise constant")
    
    parser.add_argument('-b', metavar='bottleneck', dest="bottleneck", type=float,
        default=1.,
        help="Strength of the bottleneck")

    ## More esoteric params related to changing the way the plots are drawn
    parser.add_argument("--do_plots", action="store_true",
        help="Generate a bunch of plots and animations. Default is not to do plots.")
    parser.add_argument("--curves", action='store_true',
        help="Plot rank abundance as curves rather than points")
    parser.add_argument("--plot_models", action='store_true',
        help="Add expectations under different statistical models to the rank abundance plot")

    ## Add standard quiet/force/version args
    parser.add_argument('-q', "--quiet", action='store_true',
        help="do not print to stderror or stdout.")
    parser.add_argument('-v', "--verbose", action='store_true',
        help="print lots of debug information")
    parser.add_argument('-f', "--force", action='store_true',
        help="force overwrite of existing data")
    parser.add_argument('--version', action='version',
        version="0.0.0.1")

    ## parse args
    args = parser.parse_args()

    ## Check args
    if not args.meta in ["logser", "uniform"]:
        if not os.path.exists(args.meta):
            parser.print_help()
            sys.exit("Metacommunity file doesn't exist - {}".format(args.meta))

    if args.mode not in ["volcanic", "landbridge"]:
        parser.print_help()
        sys.exit("-p option must be either volcanic or landbridge. You said - {}".\
            format(args.mode))

    return args


def progressbar(nsims, finished, msg=""):
    """ prints a progress bar """
    progress = 100*(finished / float(nsims))
    hashes = '#'*int(progress/5.)
    nohash = ' '*int(20-len(hashes))
    print("\r  [{}] {:>3}% {} ".format(hashes+nohash, int(progress), msg), end="")
    sys.stdout.flush()


if __name__ == "__main__":

    ## Parse command line arguments
    args = parse_command_line()

    if args.verbose:
        print(args)

    if args.colonizers:
        ## Implicit space and clustered immigration
        data = implicit_CI.implicit_CI(K=args.K, colrate=args.colrate, mig_clust_size=args.colonizers, quiet=args.quiet)
    else:
        ## Implicit space, one colonizer per event
        data = implicit_BI.implicit_BI(K=args.K, colrate=args.colrate, exponential=args.exponential, quiet=args.quiet)

    ## Set model parameters
    data.set_metacommunity(args.meta)
    data.prepopulate(mode=args.mode)

    ## Setup output directory and files
    try:
        if os.path.exists(args.outdir) and args.force:
            shutil.rmtree(args.outdir)
            os.mkdir(args.outdir)
        elif not os.path.exists(args.outdir):
            os.mkdir(args.outdir)
        else:
            sys.exit("Output directory exists - {}\nUse the force flag -f to overwrite".format(args.outdir))
        out = open(os.path.join(args.outdir, "pi_x_dxy.log"), "w")
        yoyofile = open(os.path.join(args.outdir, "sizechange_through_time.log"), 'w')
        stats = open(os.path.join(args.outdir, "sumstats.txt"), 'w')
        abundacesfile = open(os.path.join(args.outdir, "abundances.txt"), 'w')
        pidxyfile = open(os.path.join(args.outdir, "pidxy.txt"), 'w')
        extfile = open(os.path.join(args.outdir, "extinction_times.txt"), "w")
        megalogfile = open(os.path.join(args.outdir, "megalog.txt"), "w")
        ## write header
        megalogfile.write("\t".join(["%eq", "step", "Species_uuid", "Col_time", "Loc_Abund", "Meta_Abund", "pi", "pi_net", "Dxy",  "S", "S_island", "pi_island", "tajD_island", "S_meta", "pi_meta"]))
        megalogfile.write("\n")

        ## Make the output file properly formatted for this model and return the file object
        make_outputfile(args.model, stats)

    except Exception as inst:
        sys.exit("problem opening output for writing - {}\n{}".format(args.outdir, inst))

    ## If you don't specify a recording period, then just record
    ## 10% of the time, or every 100000 generations if you are
    ## running to equilibrium
    if not args.recording_period:
        if args.nsims < 1:
            args.recording_period = 100000
        else:
            args.recording_period = args.nsims/10

    ## Start the main loop
    start = time.time()
    elapsed = 0
    reached_equilib = False

    ## Track pop size changes through time
    ## Also track species sumstats through time so we can normalize the
    ## heatmap plots
    yoyo = collections.OrderedDict()
    sp_through_time = collections.OrderedDict()
    equilibria = collections.OrderedDict()

    ## Whether or not to care about equilibrium
    ## or just to run n timesteps
    do_equilibrium = True
    single_equilibrium = True

    ## if args.nsims < 0 just run until double equilibrium
    ## or 10e9 steps (effectively forever)
    ## if args.nsims == 0, do single equilibrium
    if not args.nsims == 0:
        single_equilibrium = False
    if args.nsims < 1:
        args.nsims = 100000000
    else:
        do_equilibrium = False
    #for i in range(1, args.nsims):
    i = 0
    while True:
        if i == args.nsims:
            break
        i += 1
        data.step(time=i, invasion_time=args.invasion_time, invasiveness=args.invasiveness)

        ## Print the progress bar every once in a while
        ## Set a flag for equilibrium. If you've reached it, flip all the
        ## founder flags back to True and keep running til next equilibrium
        ## then stop
        percent_equil = 0
        if not i % 10000:
            ## Record abundance of each species through time
            ## Make a counter for the local_community, counts the number of
            ## individuals w/in each species
            abundances = collections.Counter([x[0] for x in data.local_community])
            yoyofile.write(str(abundances)+"\n")
            yoyo[i] = abundances

            ## Test for equilibrium
            founder_flags = [x[1] for x in data.local_community]
            percent_equil = float(founder_flags.count(False))/len(founder_flags)            
            if do_equilibrium:
                if (not any(founder_flags)) and reached_equilib:
                    print("\nReached second equilibrium")
                    break
                elif not any(founder_flags):
                    if single_equilibrium:
                        print("\nReached single equilibrium.")
                        break
                    else:
                        print("\nReached first equilibrium, reset founder flags")
                        data.local_community = [(x[0], True) for x in data.local_community]
                        reached_equilib = True
                        #break

            ## Update progress bar
            secs = time.time()-start
            elapsed = datetime.timedelta(seconds=secs)
            ## Calculate the remaining time
            rate = float(i)/secs
            remaining_secs = (args.nsims - i) / rate
            progressbar(args.nsims, i, " | %equilib - {} | elapsed - {} | remaining - {}".format(percent_equil, elapsed, datetime.timedelta(seconds=int(remaining_secs))))

        ## if %equilibrium < 75 then zoom in and record data more frequently
        if (percent_equil < 0.75 and do_equilibrium):
            recording_period = int(args.recording_period / 10.)
        else:
            recording_period = args.recording_period
        ## Recording data every once in a while
        if not i % recording_period and not i == 0: 
            ## Every once in a while write out useful info
            data.simulate_seqs()
            sp_through_time[i] = data.get_species()
            ## Save a copy of the local_community
            equilibria[i] = percent_equil
            if args.bottleneck < 1:
                tmp_local = data.local_community
                data.bottleneck(args.bottleneck)
            out.write("atEQ {}\tstep {}\tpercent_equil {}\t shannon {} ".format(reached_equilib, i, percent_equil,\
                                                    shannon(data.get_abundances(octaves=False))) + heatmap_pi_dxy_ascii(data, labels=True)+"\n")
            diversity_stats = dict([(s.uuid[0], (s.pi, s.dxy)) for s in data.get_species()])

            abundacesfile.write("{} {}\n".format(percent_equil, data.get_abundances(octaves=False)))
            pidxyfile.write("{} pi {}\n".format(percent_equil, [s.pi_island for s in data.get_species()]))
            pidxyfile.write("{} dxy {}\n".format(percent_equil, [s.dxy for s in data.get_species()]))
            extfile.write("{} {}\n".format(percent_equil, " ".join([str(x) for x in data.extinction_times])))

            write_megalog(megalogfile, i, percent_equil, data)

            ## blank extinction times
            data.extinction_times = []

            for p,c in diversity_stats.values():
                if p > 0.3 or c > 0.3:
                    print(tabulate_sumstats(data))
                    print("one is fuck: pi {} dxy {}".format(p,c))
                    #sys.exit()
            ## Write to the output file
            if reached_equilib:
                eq = 1
            else:
                eq = percent_equil
            write_outfile(args.model, stats, data, eq)

            if args.bottleneck < 1:
                ## after the bottleneck put the original community back to continue the simulations
                data.local_community = tmp_local

            ## Do extra simulations per timestep
            ## for j in xrange(0,1):
            ##    data.simulate_seqs()
            ##    write_outfile(args.model, stats, data, eq)

    progressbar(100, 100, "  |  {}  steps completed  |  Total runtime   {}".format(i, elapsed))

    if not reached_equilib:
        founder_flags = [x[1] for x in data.local_community]
        percent_equil = float(founder_flags.count(False))/len(founder_flags)
        print("How close to equilibrium? {}".format(percent_equil))
    
    ## When finished simulate the final set of sequences
    data.simulate_seqs()
    sp_through_time[i] = data.get_species()
    if args.bottleneck < 1:
        tmp_local = data.local_community
        data.bottleneck(args.bottleneck)
    pidxyfile.write("{} pi {}\n".format(percent_equil, [s.pi_island for s in data.get_species()]))
    pidxyfile.write("{} dxy {}\n".format(percent_equil, [s.dxy for s in data.get_species()]))
    extfile.write("{} {}\n".format(percent_equil, " ".join([str(x) for x in data.extinction_times])))
    equilibria[i] = percent_equil
    print("Extinction rate - {}".format(data.extinctions/float(data.current_time)))
    print("Colonization rate - {}".format(data.colonizations/float(data.current_time)))
    print("Invasive survival # - {}".format(data.survived_invasives))

    ## Print out some informative business
    ## Get all results and write out final sumstats
    abundance_distribution = data.get_abundances(octaves=args.octaves)
    if not args.quiet:
        print(plot_abundances_ascii(abundance_distribution))
        print(tabulate_sumstats(data))

    ## Write out normalized pi_x_dxy heatmaps to the sumstats file
    write_heats_to_outfile(args.model, stats, data, sp_through_time, equilibria)

    ## Write out to log files
    write_sizechanges(args.outdir, yoyo)
    with open(os.path.join(args.outdir, "gimmeSAD.out"), "w") as stats:
        stats.write("Parameters - {}\n".format(args))
        stats.write("Raw abundance dist - {}\n".format(data.get_abundances(octaves=False)))
        stats.write("Abundance in octaves - {}\n".format(data.get_abundances(octaves=True)))
        stats.write("Shannon's entropy - {}\n".format(shannon(data.get_abundances(octaves=False))))
        stats.write(plot_abundances_ascii(abundance_distribution))
        stats.write("\n")
        stats.write(tabulate_sumstats(data))


    if args.do_plots:
        ## Make the normalized pi_x_dxy heatmaps
        plot_rank_abundance_through_time(args.outdir, sp_through_time, equilibria,\
                    stats_models=args.plot_models, as_curve=args.curves, verbose=args.verbose)
        normalized_pi_dxy_heatmaps(args.outdir, sp_through_time, equilibria,\
                    stats_models=args.plot_models, as_curve=args.curves, verbose=args.verbose)
        normalized_pi_dxy_heatmaps(args.outdir, sp_through_time, equilibria,\
                    stats_models=args.plot_models, as_curve=args.curves, one_d=True, verbose=args.verbose)
        plot_abundance_vs_colonization_time(args.outdir, sp_through_time, equilibria,\
                    stats_models=args.plot_models, as_curve=args.curves)

