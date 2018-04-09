
from __future__ import print_function
from ascii_graph import Pyasciigraph
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
from tabulate import tabulate
import numpy as np
import pandas as pd
import math
import os

plt.switch_backend('agg')

import MESS
from MESS.stats import SAD
from MESS.util import *

###################################################################
## These are old, janky plotting methods imported from gimmeSAD
## so most people probably won't want to mess with this as
## they require weird formats of input that are not straightforward
## to accumulate. Best not to try to use them if you're not me.
###################################################################

def plot_rank_abundance_through_time(outdir, sp_through_time, equilibria,\
                                    stats_models=False, as_curve=False,\
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
    #if only_extant:
    #    sp_through_time = prune_extant(sp_through_time)

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
        progressbar(tot_plots, i+1)

        title = "Time_"+str(file_index+i)
        write = os.path.join(abund_out, title)

        ## Make the Plot
        fig = plt.figure(figsize=(12,5))

        ## Make the SAD subplot
        abund = SAD([x.abundance for x in species], from_abundances=True, octaves=True)
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
                                    stats_models=False, as_curve=False, only_extant=False,\
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
    max_pi_local = 0
    max_dxy = 0
    ## For each recorded timeslice
    my_dxys = []
    my_pi_locals = []

    ## Get variables we care about
    max_n_species, max_abundance, _, _, _, _ = prep_normalized_plots(sp_through_time)
    if verbose:
        print("max_n_species - {}\nmax_abundance - {}".format(max_n_species, max_abundance))

    ## Normalization routines
    for sp_list in sp_through_time.values():

        ## Get max pi and max dxy
        pis = np.array([(x.dxy, x.pi_local) for x in sp_list if x.uuid[0] in extant])
        ## pis will be empty if this timeslice includes no extant species
        if pis.any():
            my_max_dxy = max([x[0] for x in pis])
            my_max_pi_local = max([x[1] for x in pis])
            if max_dxy < my_max_dxy:
                max_dxy = my_max_dxy
            if max_pi_local < my_max_pi_local:
                max_pi_local = my_max_pi_local
            my_dxys.append(my_max_dxy)
            my_pi_locals.append(my_max_pi_local)

#    max_dxy = np.average(my_dxys)
#    max_pi_local = np.average(my_pi_locals)
    max_dxy = np.median(my_dxys)
    max_pi_local = np.median(my_pi_locals)
    ## However this function is calculating the heatmaps is fucked up, so you have to
    ## hard code max values for pi and dxy here.
    max_dxy = 0.04
    max_pi_local = 0.02
    if verbose:
        print("Got\tmax_dxy - {}\t max_pi_local - {}".format(max_dxy, max_pi_local))
    ## Make the heatmaps, one for each timeslice
    ## TODO: Should we include all species or just the ones that live to the end?
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the heatmap
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))

    max_bin_val = get_max_heat_bin(sp_through_time, max_pi_local, max_dxy)
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
        #pis = np.array([(x.dxy, x.pi_local) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.dxy, x.pi_local) for x in sp_list])

        ## Empty heatmap we'll write into
        heat = np.zeros((20,20), dtype=np.int)

        ## Make the bins
        dxy_bins = np.linspace(0, max_dxy, 20, endpoint=True)
        pi_local_bins = np.linspace(0, max_pi_local, 20, endpoint=True)

        ## Now you have the bins each value belongs in, but you need to 
        ## go through and populate the heat matrix
        for dxy, pi_local in pis:
            count_dxy = 0
            count_pi_local = 0
            try:
                while not dxy <= dxy_bins[count_dxy]:
                    count_dxy += 1
                while not pi_local <= pi_local_bins[count_pi_local]:
                    count_pi_local += 1
                ## increment the heatmap point this corresponds to
                heat[count_dxy][count_pi_local] += 1
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
        plt.xticks(np.arange(len(pi_local_bins)), ["{0:.4f}".format(x) for x in pi_local_bins], rotation='vertical')

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
#        ax1.set_xticks(np.arange(len(pi_local_bins) + 0.5))
#        ax1.set_xticklabels(["{0:.4f}".format(x) for x in pi_local_bins], rotation="vertical", ha="center")

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
        abund = SAD([x.abundance for x in sp], from_abundances=True, octaves=True)
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


def plot_abund_vs_colon(species, max_coltime, max_abundance):
    x = [np.log10(s.colonization_time) for s in species]
    y = [np.log10(s.abundance) for s in species]
    plt.scatter(x, y, color="blue", s=100)
    plt.ylim(0, int(math.ceil(np.log10(max_abundance))))
    plt.xlim(0, 8) #int(math.ceil(np.log10(max_coltime))))
    plt.title("Abundance vs Colonization Time", fontsize=24)
    plt.ylabel("Abundance (log10)", fontsize=20)
    plt.xlabel("Colonization Time (log10)", fontsize=20)


def get_max_heat_bin(sp_through_time, max_pi_local, max_dxy):
    max_heat_bin = 0

    for sp_list in sp_through_time.values():
        ## Get the sumstats for this timeslice
        ## Only include extant species in plots
        #pis = np.array([(x.dxy, x.pi_local) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.dxy, x.pi_local) for x in sp_list])

        ## Empty heatmap we'll write into
        heat = np.zeros((20,20), dtype=np.int)

        ## Make the bins
        dxy_bins = np.linspace(0, max_dxy, 20, endpoint=True)
        pi_local_bins = np.linspace(0, max_pi_local, 20, endpoint=True)

        ## Now you have the bins each value belongs in, but you need to 
        ## go through and populate the heat matrix
        for dxy, pi_local in pis:
            count_dxy = 0
            count_pi_local = 0
            try:
                while not dxy <= dxy_bins[count_dxy]:
                    count_dxy += 1
                while not pi_local <= pi_local_bins[count_pi_local]:
                    count_pi_local += 1
                ## increment the heatmap point this corresponds to
                heat[count_dxy][count_pi_local] += 1
            except Exception as inst:
                ## Got a value bigger than our current max pi/dxy. ignore.
                pass
        if np.amax(heat) > max_heat_bin:
            max_heat_bin = np.amax(heat)

    return max_heat_bin


def make_animated_gif(datadir, outfile):
    """ This function will take all png files in a directory and make them
    into an animated gif. The inputs are the directory with all the images
    and the full path including filename of the file to write out"""

    ## Do the imagemagick conversion, if possible
    ## `convert -delay 100 outdir/* anim.gif`
    ## TODO: Define the total time to be 10 seconds, total available
    ## timeslots is * 100 bcz convert flag is in 1/100 seconds
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


if __name__ == "__main__":
    import collections
    import numpy as np
    reg = MESS.Region("tmp")
    loc = MESS.LocalCommunity("tmp", K=1000)
    reg._link_local(loc)
    loc.step(1000)
    loc.simulate_seqs()
    print(loc.species_objects)
    print(SAD([x.abundance for x in loc.species_objects], from_abundances=True))
    print(SAD([x.abundance for x in loc.species_objects], from_abundances=True, octaves=True))

