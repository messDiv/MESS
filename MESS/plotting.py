

import matplotlib
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    matplotlib.use("agg")
import math
import matplotlib.patheffects as PathEffects
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import shutil
import subprocess

import MESS
from MESS.stats import SAD
from MESS.util import *
from sklearn.decomposition import PCA
from sklearn.preprocessing import PowerTransformer, StandardScaler

## A dictionary mapping mModel parameters in params/SIMOUT file format
## to unicode/prettier versions for plotting
##
## Helpful matplotlib info for plotting unicode characters:
## https://matplotlib.org/users/mathtext.html#subscripts-and-superscripts
target_labels = {"mutation_rate":"\u03BC",\
                "alpha":"\u03B1",\
                "S_m":r"$S_M$",\
                "J_m":r"$J_M$",\
                "speciation_rate":"\u03BB",\
                "death_proportion":"\u03B5",\
                "trait_rate_meta":r"$\sigma^2_M",\
                "ecological_strength":r"$s_E$",\
                "J":r"$J$",\
                "m":r"$m$",\
                "speciation_prob":"\u03BD",\
                "_lambda":"\u039B",\
                "generation":"generations"}

## Defaults
#model_colors = {"neutral":"blue",\
#                "filtering":"orange",\
#                "competition":"red"}


## Blue/Orange/Red but custom. A little dark.
#model_colors = {"neutral":"#009FE8",\
#                "filtering":"#F5B700",\
#                "competition":"#FC5130"}

## B/O/R but a little Electric
#model_colors = {"neutral":"#7744AE",\
#                "filtering":"#FFDE19",\
#                "competition":"#FF1E57"}

## B/O/R but Too Neon
#model_colors = {"neutral":"#00A4CC",\
#                "filtering":"#EDFF00",\
#                "competition":"#FF3EA5"}

# ## Miami Dolphins
# model_colors = {"neutral":"#FF8800",\
#                 "filtering":"#00B7C4",\
#                 "competition":"#005679"}

## Primary colors
model_colors = {"neutral":"#375E97",\
                "filtering":"#3F681C",\
                "competition":"#FB6542",\
                "pairwise_competition":"#FFBB00"}


def _filter_sims(simfile,\
                    feature_set='',\
                    nsims=0,\
                    normalize_hills=False,\
                    select='',\
                    tol='',\
                    verbose=False):
    """
    Load simulation data and perform filtering and downsampling that's common
    across many plotting routines. Normally you won't call this directly.

    :param str simfile: 
    :param list feature_set:
    :param int nsims:
    :param bool normalize_hills: Whether to divide all Hill numbers by `S` to
        normalize them to a proportion and promote comparison across
        communities of different sizes.
    :param int/float select: 
    :param int/float tol:
    :param bool verbose: Whether to print progress messages.

    :return: Returns a tuple of pd.DataFrame containing the community assembly
        model class labels for retained simulations and a pd.DataFrame of 
        filtered and pruned simulation summary statistics.
    """
    ## Load the simulations
    sim_df = pd.read_csv(simfile, sep="\t", header=0)

    ## Wether to select only specific timepoints bracketed by `tol` or
    ## just plot everything.
    if select is '':
        pass
    else:
        ## Switch on whether to plot based on lambda or based on generations
        ## Set default tolerance to be pretty small
        if select > 1:
            time = "generation"
            if not tol: tol = 1
        else:
            time = "_lambda"
            if not tol: tol = 0.01

        ## Apply the filter +/- `tol`.
        sim_df = sim_df[(sim_df[time] < (select + tol)) & (sim_df[time] > (select - tol))]

    ## If feature_set is unspecified assume we're using all features.
    if not feature_set:
        feature_set = MESS.stats.feature_sets()["all"]

    ## Normalize all Hill #s by dividing by S, so you can compare across
    ## different sized communities.
    if normalize_hills:
        sim_df[[x for x in sim_df.columns if "_h" in x]] =\
                sim_df[[x for x in sim_df.columns if "_h" in x]]\
                .apply(lambda x: x/sim_df["S"])

    ## Prune the simulations based on selected features and number of
    ## simulations to retain.
    if not nsims:
        nsims = len(sim_df)
    labels = sim_df["community_assembly_model"][:nsims]
    sim_df = sim_df[feature_set][:nsims]

    ## Remove invariant targets (save time)
    sim_df = sim_df.loc[:, (sim_df != sim_df.iloc[0]).any()]
    retained = list(sim_df.columns)
    if verbose: print("Removed invariant targets. Retained: {}".format(list(retained)))

    return labels, sim_df


def plot_simulations_hist(simfile,\
                        ax='',\
                        figsize=(20, 20),\
                        feature_set='',\
                        nsims=1000,\
                        normalize_hills=False,\
                        bins=20,\
                        alpha=0.6,\
                        select='',\
                        tol='',\
                        title='',\
                        outfile='',\
                        verbose=False):
    """
    Simple histogram for each summary statistic. Useful for inspecting model
    performance. Invariant summary statistics will be removed.

    :param str simfile: 
    :param tuple figsize:
    :param list feature_set:
    :param int nsims:
    :param bool normalize_hills: Whether to divide all Hill numbers by `S` to
        normalize them to a proportion and promote comparison across
        communities of different sizes.
    :param int bins: The number of bins per histogram.
    :param float alpha: Set alpha value to determine transparency [0-1], larger
        values increase opacity.
    :param int/float select: 
    :param int/float tol:
    :param str title:
    :param str outfile:
    :param bool verbose:

    :return: Return a list of `matplotlib.pyplot.axis` on which the simulated
        summary statistics have been plotted. This list can be _long_ depending
        on how many statistics you plot.
    """

    ## Filter and downsample the simulations
    labels, sim_df = _filter_sims(simfile,\
                            feature_set=feature_set,\
                            nsims=nsims,\
                            normalize_hills=normalize_hills,\
                            select=select,\
                            tol=tol,\
                            verbose=verbose)

    neut_df = sim_df[labels.values == "neutral"]
    filt_df = sim_df[labels.values == "filtering"]
    comp_df = sim_df[labels.values == "competition"]
    pw_df   = sim_df[labels.values == "pairwise_competition"]
    if verbose: print("Nsims\n  neutral\t{}\n  filtering\t{}\n  competition\t{}\n pairwise_competition\t{}"\
                        .format(len(neut_df), len(filt_df), len(comp_df), len(pw_df)))
    
    ## TODO: Would be cool to have an option to plot kde instead of hist.
    ## i.e. neut_df.plot(kind='kde'). Here it is, but it's untested-ish.
    #bw_method=.5
    #axs = neut_df.plot(kind='kde', figsize=figsize, label="neutral", alpha=alpha,
    #                color=MESS.plotting.model_colors["neutral"],  grid=False, bw_method=bw_method)

    #axs = axs.flatten()[:len(sim_df.columns)]
    #_ = filt_df.plot(kind='kde', ax = axs, label="filtering", alpha=alpha,\
    #                color=MESS.plotting.model_colors["filtering"], grid=False, bw_method=bw_method)
    #_ = comp_df.plot(kind='kde', ax = axs, label="competition", alpha=alpha,\
    #                color=MESS.plotting.model_colors["competition"], grid=False, bw_method=bw_method)

    axs = neut_df.hist(figsize=figsize, label="neutral", alpha=alpha, bins=bins,
                    color=MESS.plotting.model_colors["neutral"],  grid=False)

    ## Flatten the list of axes and trim to make sure there's only exactly the
    ## right number to match the number of summary stats retained.
    axs = axs.flatten()[:len(sim_df.columns)]
    _ = filt_df.hist(ax = axs, label="filtering", alpha=alpha, bins=bins,\
                    color=MESS.plotting.model_colors["filtering"], grid=False)
    _ = comp_df.hist(ax = axs, label="competition", alpha=alpha, bins=bins,\
                    color=MESS.plotting.model_colors["competition"], grid=False)

    _ = pw_df.hist(ax = axs, label="pairwise competition", alpha=alpha, bins=bins,\
                    color=MESS.plotting.model_colors["pairwise_competition"], grid=False)

    plt.tight_layout()
    return axs


def plot_simulations_pca(simfile, ax='',\
                            figsize=(8, 8),\
                            feature_set='',\
                            loadings=False,\
                            nsims=1000,\
                            normalize_hills=False,\
                            select='',\
                            tol='',\
                            title='',\
                            outfile='',\
                            verbose=False):
    """
    Plot summary statistics for simulations projected into PC space.

    :param str simfile: 
    :param matplotlib.pyplot.axis ax:
    :param tuple figsize:
    :param list feature_set:
    :param bool loadings: BROKEN! Whether to plot the loadings in the figure.
    :param int nsims:
    :param bool normalize_hills: Whether to divide all Hill numbers by `S` to
        normalize them to a proportion and promote comparison across
        communities of different sizes.
    :param int/float select: 
    :param int/float tol:
    :param str title:
    :param str outfile:
    :param bool verbose:

    :return: Return the `matplotlib.pyplot.axis` on which the simulations are
        plotted.
    """
    if not ax:
        fig, ax = plt.subplots(figsize=figsize)

    ## Filter and downsample the simulations
    labels, sim_df = _filter_sims(simfile,\
                            feature_set=feature_set,\
                            nsims=nsims,\
                            normalize_hills=normalize_hills,\
                            select=select,\
                            tol=tol,\
                            verbose=verbose)

#    sim_df = StandardScaler().fit_transform(sim_df)
    sim_df = PowerTransformer(method='yeo-johnson').fit_transform(sim_df)

    pca = PCA(n_components=2)
    dat = pca.fit_transform(sim_df)

    ax.scatter(dat[:, 0], dat[:, 1], color=[MESS.plotting.model_colors[x] for x in labels])

    ## Remove a bunch of visual noise
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.tick_params(top='off', bottom='off', left='off', right='off')

    if title:
        ax.set_title(title)

    ## TODO: Doesn't work how I'd like.
    ##print("Explained variance", pca.explained_variance_ratio_)
    ##if loadings:
    ##    for i, comp in enumerate(pca.components_.T):
    ##        plt.arrow(0, 0, pca.components_.T[i,0], pca.components_.T[i,1], color = 'r',alpha = 0.5)
    ##        plt.text(pca.components_.T[i,0]* 1.5, pca.components_.T[i,1] * 1.5, dat[i+2], color = 'black', ha = 'center', va = 'center')

    ## If writing to file then don't plot to screen.
    if outfile:
        try:
            plt.savefig(outfile)
            if verbose: print("Wrote figure to: {}".format(outfile))
        except Exception as inst:
            raise Exception("Failed saving figure: {}".format(inst))
        plt.close()

    return ax


def plots_through_time(simfile,\
                        plot_fun,\
                        tmax=1,\
                        ntimes=100,\
                        outdir='',\
                        outgif='',\
                        verbose=False,\
                        **kwargs):
    """
    :param str simfile:
    :param function plot_fun:
    :param int/float tmax: The maximum time to consider when partitioning the
        simulations. If `tmax` <=1 then we plot `_lambda`, and if > 1 we plot
        `generation`.
    :param int ntimes: The number of time slices to carve the simulations up into.
    :param str outdir:
    :param str outfile:
    :param bool verbose: Whether to print progress information.
    :param dict kwargs: Other parameters to pass to the plotting function.
    """

    try:
        subprocess.check_output("convert")
    except OSError:
        ## 'convert' doesn't exist. Warn and raise.
        raise MESSError(REQUIRE_IMAGEMAGICK_ERROR)
    except subprocess.CalledProcessError:
        ## This is okay. 'convert' returns 1 when called with no args.
        pass

    if not outdir:
        ## Clean up previous runs.
        outdir = "/tmp/MESS_plots"
        if os.path.exists(outdir):
            shutil.rmtree(outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ## If not specified, set the default tolerance to be inversely proportional
    ## to the number of selected time slices. Divide the tolerance by 2 so the
    ## spans of simulations don't overlap.
    if not "tol" in kwargs: tol = tmax/ntimes/2

    for i, time in enumerate(np.linspace(0, tmax, ntimes)):
        if verbose: MESS.util.progressbar(ntimes, i, "Generating PCA plots") 
        outfile = os.path.join(outdir, str(time) + ".png")
        _ = plot_fun(simfile,\
                        select=time,\
                        tol=tol,\
                        title="{0:.2f}".format(time),\
                        outfile=outfile,\
                        **kwargs)

    if verbose: MESS.util.progressbar(ntimes, i, "Creating animation") 

    ## Actually generate the animation now.
    if not outgif:
        outgif = os.path.join(outdir, "{}.gif".format(plot_fun.__name__))
    _make_animated_gif(outdir, outgif)

    if verbose: MESS.util.progressbar(100, 100, "Creating animation") 
 
    
def _make_animated_gif(datadir, outfile, delay=50):
    """
    This function will take all png files in a directory and make them
    into an animated gif. The inputs are the directory with all the images
    and the full path including filename of the file to write out

    :param str datadir: Directory that contains all the component files
        for the animation. These should be .png, and should be alpha-sorted
        in the order of the animation.
    :param str outfile: The name of the file to write the animated gif to.
    :param int delay: Time delay between frame changes in 1/100 second
        increments.
    """

    ## Do the imagemagick conversion, if possible
    ## `convert -delay 100 outdir/* anim.gif`
    ## TODO: Define the total time to be 10 seconds, total available
    ## timeslots is * 100 bcz convert flag is in 1/100 seconds
    ## Default half second intervals
    cmd = "convert -delay {} ".format(delay)\
            + datadir + "/*.png "\
            + outfile
    try:
        res = subprocess.check_output(cmd, shell=True)
    except Exception as inst:
        print("Trouble creating abundances through time animated gif - {}".format(inst))
        print("You probably don't have imagemagick installed")


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

    ## Make the output directory for heatmaps inside the top level output directory
    abund_out = os.path.join(outdir, "abundance_plots")
    if not os.path.exists(abund_out):
        os.mkdir(abund_out)

    times = list(equilibria.keys())
    equilibria = list(equilibria.values())

    ## If you only want to see extant species then prune all the extinct ones
    #if only_extant:
    #    sp_through_time = prune_extant(sp_through_time)

    max_n_species, max_abundance, max_octave, max_class_count, max_n_bins, octave_bin_labels = prep_normalized_plots(np.nan_to_num(sp_through_time))
    if verbose:
        print("info:\n\nmax_n_species - {}\nmax_abundance - {}\nmax_octave - {}\nmax_class_count - {}\nmax_n_bins - {}\noctave_bin_labels - {}\n".format(\
                max_n_species, max_abundance, max_octave, max_class_count,\
                max_n_bins, octave_bin_labels))

    ## Get a list of UUIDs of the extant species at time 0 (present)
    extant = [x.stats["name"] for x in list(sp_through_time.values())[-1]]

    ## Make the abundance plots, one for each timeslice
    ## TODO: Should we include all species or just the ones that live to the end?
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the abundance
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))

    tot_plots = len(sp_through_time)
    for i, species in enumerate(sp_through_time.values()):
        progressbar(tot_plots, i+1)
        progressbar(tot_plots, i+1)

        title = "Time_"+str(file_index+i)
        write = os.path.join(abund_out, title)

        ## Make the Plot
        fig = plt.figure(figsize=(12,5))

        ## Make the SAD subplot, convert Ne back to abundance by dividing by alpha
        if any(np.isnan([x.stats["abundance"] for x in species])):
            continue
        abund = SAD([x.stats["abundance"] for x in species], from_abundances=True, octaves=True)        
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
    _make_animated_gif(abund_out,\
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

    times = list(equilibria.keys())
    equilibria = list(equilibria.values())

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
    extant = [x.stats["name"] for x in list(sp_through_time.values())[-1]]
    nslices = len(sp_through_time)
    ## Create a file index that is large enough so we can create all the abundance
    ## files in an order where we can cat them in numerical order w/o too much fuss
    file_index = 10**len(list(str(nslices)))

    tot_plots = len(sp_through_time)
    for i, species in enumerate(sp_through_time.values()):
        progressbar(tot_plots, i+1)

        title = "Time_"+str(file_index+i)
        write = os.path.join(abund_out, title)

        ## Make the Plot
        fig = plt.figure(figsize=(12,5))

        ## Make the abundance vs colonization time plot
        ax1 = plt.subplot(121)
        max_coltime = max([s.stats["tdiv"]for s in species])
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
    _make_animated_gif(abund_out,\
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

    times = list(equilibria.keys())
    equilibria = list(equilibria.values())

    ## If you only want to see extant species then prune all the extinct ones
    if only_extant:
        sp_through_time = prune_extant(sp_through_time)

    ## GET MAX pi and dxy
    ## Get a list of UUIDs of the extant species at time 0 (present)
    extant = [x.stats["name"] for x in list(sp_through_time.values())[-1]]
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
    for sp_list in list(sp_through_time.values()):

        ## Get max pi and max dxy
        pis = np.array([(x.stats["dxy"], x.stats["pi_local"]) for x in sp_list if x.stats["name"] in extant])
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
    ## TODO: This is hackish
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
    levels = list(range(cbar_min, cbar_max+step, step))
    my_colorbar = plt.contourf(Z, levels, cmap=plt.cm.jet)
    plt.clf()


    tot_heatmaps = len(sp_through_time)
    for i, sp_list in enumerate(sp_through_time.values()):
        progressbar(tot_heatmaps, i+1)

        title = "Time_"+str(file_index+i)
        write = os.path.join(heat_out, title)
        #print("Doing", title)

        ## Get the sumstats for this timeslice
        ## Only include extant species in plots
        #pis = np.array([(x.dxy, x.pi_local) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.stats["dxy"], x.stats["pi_local"]) for x in sp_list])

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

        ## Python 2to3 unicode printing incompatibilty that i don't feel
        ## like figuring out because this is old code.
        #plt.xlabel(u"Nucleotide diversity (\u03c0)", fontsize=20)
        plt.xlabel("Nucleotide diversity (pi)", fontsize=20)
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

    _make_animated_gif(heat_out,\
                        os.path.join(outdir, outfile))


def prep_normalized_plots(sp_through_time):
    ## GET MAX values for abundance and num species so we can normalize the plot axes
    max_n_species = max([len(x) for x in list(sp_through_time.values())])
    max_abundance = max([max([y.stats["Ne_local"]/y.paramsdict["alpha"] for y in sp]) for sp in list(sp_through_time.values())])

    ## Get max values for abundance class count and abundance octave
    max_octave = 0
    max_class_count = 0
    max_n_bins = 0
    octave_bin_labels = []
    for sp in list(sp_through_time.values()):
        abund = SAD([x.stats["abundance"] for x in sp], from_abundances=True, octaves=True)
        octave = max(abund.keys())
        class_count = max(abund.values())
        if octave > max_octave:
            max_octave = octave
        if class_count > max_class_count:
            max_class_count = class_count
        if len(abund) > max_n_bins:
            max_n_bins = len(abund)
            octave_bin_labels = list(abund.keys())

    return max_n_species, max_abundance, max_octave, max_class_count, max_n_bins, octave_bin_labels


def plot_sad(abund, max_n_species, max_n_bins, max_class_count, octave_bin_labels, verbose):
    ax1 = plt.gca()
    ab_class = list(abund.keys())
    count = list(abund.values())
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
        Y = [xx.stats["Ne_local"]/xx.paramsdict["alpha"] for xx in species]
        plt.semilogy(X, Y, label="simulated")
        ymax = max_abundance
    else:
        Y = [np.log10(xx.stats["Ne_local"]/xx.paramsdict["alpha"]) for xx in species]
        plt.scatter(X, Y, color="blue", s=100, label="simulated")
        ymax = int(math.ceil(np.log10(max_abundance)))

    plt.title("Rank Abundance", fontsize=24)
    plt.xlim(0, max_n_species)
    plt.ylim(0, ymax)
    plt.ylabel("Abundance (log10)", fontsize=20)
    plt.xlabel("Rank", fontsize=20)

    ## Whether or not to include a couple common statistical models in the plots
##    if stats_models:
##        import macroeco as meco
##        abund = [xx.stats["Ne_local"]/xx.paramsdict["alpha"] for xx in species]
##        ## Lognormal
##        mu, s = meco.models.lognorm.fit_mle(abund)
##        lognorm_rad = meco.models.lognorm.rank(len(abund), mu, s)
##        if as_curve:
##            Y = lognorm_rad[::-1]
##            plt.semilogy(X, Y, label="Lognorm RAD")
##        else:
##            Y = [int(math.ceil(np.log10(x))) for x in lognorm_rad[::-1]]
##            plt.scatter(X, Y, s=100, color="green", label="Lognorm RAD")
##        ## Logseries
##        p = meco.models.logser.fit_mle(abund)
##        logser_rad = meco.models.logser.rank(len(abund), p)
##        if as_curve:
##            Y = logser_rad[::-1]
##            plt.semilogy(X, Y, label="Logseries RAD")
##        else:
##            Y = [int(math.ceil(np.log10(x))) for x in logser_rad[::-1]]
##            plt.scatter(X, Y, s=100, color="red", label="Logseries RAD")
##        ## Poisson Lognormal
##        mu, s = meco.models.plnorm_ztrunc.fit_mle(abund)
##        plnorm_rad = meco.models.plnorm_ztrunc.rank(len(abund), mu, s)
##        if as_curve:
##            Y = plnorm_rad[::-1]
##            plt.semilogy(X, Y, label="Logseries RAD")
##        else:
##            Y = [int(math.ceil(np.log10(x))) for x in plnorm_rad[::-1]]
##            plt.scatter(X, Y, s=100, color="red", label="Poisson Lognorm RAD")
##
##        plt.legend()


def plot_abund_vs_colon(species, max_coltime, max_abundance):
    x = [np.log10(s.stats["tdiv"]) for s in species]
    y = [np.log10(s.stats["abundance"]) for s in species]
    plt.scatter(x, y, color="blue", s=100)
    plt.ylim(0, 1) #int(math.ceil(np.log10(max_abundance))))
    plt.xlim(0, 8) #int(math.ceil(np.log10(max_coltime))))
    plt.title("Abundance vs Colonization Time", fontsize=24)
    plt.ylabel("Abundance (log10)", fontsize=20)
    plt.xlabel("Colonization Time (log10)", fontsize=20)


def get_max_heat_bin(sp_through_time, max_pi_local, max_dxy):
    max_heat_bin = 0

    for sp_list in list(sp_through_time.values()):
        ## Get the sumstats for this timeslice
        ## Only include extant species in plots
        #pis = np.array([(x.dxy, x.pi_local) for x in sp_list if x.uuid[0] in extant])
        ## include all species at each timeslice
        pis = np.array([(x.stats["dxy"], x.stats["pi_local"]) for x in sp_list])

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


def plot_death_probs(death_probs, outdir, model, local_community_record):
    """
    deaths_probs arg should be a dictionnary of time : death_probabilities
    This function aims at generating an animated gif of the distribution of death_probabilities through time
    """
    ## Create a directory to store the individual images
    print("Generating death probabilities plots thourgh time")
    out_dir = outdir+"/plots_death_probs"
    os.mkdir(out_dir)
    maxtime = len(death_probs.keys())
    names = sorted([str(1000000000+i) for i in range(maxtime)])
    values = [x for x in np.nan_to_num(death_probs.values())]
    maxprob = max(np.reshape(values,len(values[0])*len(values)))
    minprob = min(np.reshape(values,len(values[0])*len(values)))

    tot_plots = len(death_probs)
    for i,time in enumerate(death_probs.keys()):
        progressbar(tot_plots, i+1)
        fig = plt.figure(figsize=(12,5))
        plt.plot(sorted(np.nan_to_num(death_probs[time])))
        plt.ylabel('Probability of death')
        plt.ylim(minprob-minprob/100,maxprob+maxprob/100)
        plt.yscale('log')
        plt.xlabel('Rank')
        plt.title('Distribution of death probabilities at time '+str(time))
        fig.savefig(out_dir+'/'+names[i]+'.png')
        plt.close()


    # Actually generate the animation now.
    namegif = model + "_death_probabilities.gif"
    outgif = os.path.join(out_dir,namegif)
    progressbar(100, 100, "\n")
    _make_animated_gif(out_dir, outgif)

    out_dir2 = outdir+"/plots_death_probs_per_species_"+model
    os.mkdir(out_dir2)

    # Have a different view with death probs for species
    for i,time in enumerate(death_probs.keys()):
        progressbar(tot_plots, i+1)
        fig = plt.figure(figsize=(12,5))
        plt.scatter(local_community_record[time],np.nan_to_num(death_probs[time]))
        plt.ylabel('Probability of death')
        plt.ylim(minprob-minprob/100,maxprob+maxprob/100)
        plt.yscale('log')
        plt.xlabel('species')
        plt.title('Distribution of death probabilities at time '+str(time))
        fig.savefig(out_dir2+'/'+names[i]+'.png')
        plt.close()


    # Actually generate the animation now.

    namegif = model + "_death_probabilities_persp.gif"
    outgif = os.path.join(out_dir2,namegif)
    progressbar(100, 100, "\n")
    _make_animated_gif(out_dir2, outgif)



def plot_traits_repartition(outdir, local_traits_through_time, death_probs):
    """
    local_traits_through_time arg should be a dictionnary of time : local_traits
    This function aims at generating an animated gif of the distribution of traits in the local community through time
    """
    ## Create a directory to store the individual images
    print("Generating local traits plots thourgh time")
    trait_out = os.path.join(outdir, "traits_through_time")
    if not os.path.exists(trait_out):
        os.mkdir(trait_out)

    maxtime = len(local_traits_through_time.keys())
    names = sorted([str(1000000000+i) for i in range(maxtime)])
    values = [x for x in np.nan_to_num(death_probs.values())]
    maxprob = max(np.reshape(values,len(values[0])*len(values)))
    minprob = min(np.reshape(values,len(values[0])*len(values)))
    
    values = [x for x in local_traits_through_time.values()]
    maxtrait = max(np.reshape(values,len(values[0])*len(values)))
    mintrait = min(np.reshape(values,len(values[0])*len(values)))


    tot_plots = len(local_traits_through_time)
    for i,time in enumerate(local_traits_through_time.keys()):
        progressbar(tot_plots, i+1)
        fig = plt.figure(figsize=(12,5))
        ax = fig.add_subplot(111, projection='3d')
        x = local_traits_through_time[time]
        y = np.nan_to_num(death_probs[time])
        hist, xedges, yedges = np.histogram2d(x, y, bins=10, range=[[mintrait, maxtrait], [minprob, maxprob]])

        # Construct arrays for the anchor positions of the 16 bars.
        xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25, indexing="ij")
        xpos = xpos.ravel()
        ypos = ypos.ravel()
        zpos = 0

        # Construct arrays with the dimensions for the 16 bars.
        dx = dy = 0.5 * np.ones_like(zpos)
        dz = hist.ravel()

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')


        plt.ylabel('death probability')
        # plt.ylim(minprob-minprob/100,maxprob+maxprob/100)
        plt.xlabel('trait value')
        # plt.xlim(mintrait-abs(mintrait)/100,maxtrait+maxtrait/100)
        plt.title('Death probabilities and trait repartitions'+str(time))
        fig.savefig(trait_out+'/'+names[i]+'.png')
        plt.close()







REQUIRE_IMAGEMAGICK_ERROR = """
The plots_through_time() function requires the image-magick graphics
processing package which may be installed with conda:

    conda install -c conda-forge imagemagick -y
"""


if __name__ == "__main__":
    import collections
    import numpy as np
    reg = MESS.Region("tmp")
    loc = MESS.LocalCommunity("tmp", J=1000)
    reg._link_local(loc)
    loc.step(1000)
    loc._simulate_seqs()
    print(loc.species_objects)
    print(SAD([x.abundance for x in loc.species_objects], from_abundances=True))
    print(SAD([x.abundance for x in loc.species_objects], from_abundances=True, octaves=True))

