Conda package for MESS. The modules required to build msprime
need to be pulled in and package build time like this:

conda build conda.recipe/MESS -c conda-forge

The modules required to run msprime need to be copied from conda-forge
to the anaconda mess channel (this only needs to be done once, and
every so often if the packages get updated). Then when users run

conda install -c mess mess

it'll just work. The 'macroeco' package is an exception, since it's
not in conda-forge we actually have to maintain the conda package too,
but i'm half interested in just ditching it.
