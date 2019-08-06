# MESS - Massive Eco-evolutionary Synthesis Simulations

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/messDiv/MESS/master)

Adding traits and phylogenies to abundances and genetic diversity

## Introduction
A central challenge in understanding the origins of biodiversity is that, while ecological phenomena occur across observably short time periods, the longer-term drivers and outcomes of these ecological processes can often only be indirectly infered. MESS is an inferential model at the interface between macroecology, macroevolution and population-level processes, which can be applied to data from geological or ecological chronosequences that present communities of different ages. Inferences from these snapshots in time thereby allow for model validation and a link between direct observational methods for local communities and models that make indirect inferences underlying community history using genetic and phylogenetic patterns. The MESS model directly links ecological theories and models of community composition and comparative population genomics, all within a temporal framework. Our approach is to build a unified model bridging theory from phylogenetic and comparative population genomics with ecological theory, so as to understand the history underlying patterns of species diversity. This model can be used to make joint predictions of species abundances and genetic diversities over time. This unified approach bridges ecological and evolutionary theory to elucidate processes responsible for origins and maintenance of species diversity and provide a framework for making predictions about biodiversity dynamics.

Critical support for development of the MESS model was provided by the German Centre for Integrative Biodiversity Research (iDiv) and the Santa Fe Institute.

# Installation
MESS is distributed as a conda package, which makes installing all the dependencies
transparent to the user. The classic conda install goes like this:

* Install [conda](https://docs.conda.io/en/latest/miniconda.html) for python 2.7 (I know.... I know... Soon for python 3)
* Create a new conda environment: `conda create --name MESS` & `conda activate MESS`
* Install MESS: `conda install -c conda-forge -c mess mess`

**NB: Hackers only.** You can also install from a clean exported environment 
file available in github repo conda.build directory:

* Install conda
* `wget https://raw.githubusercontent.com/messDiv/MESS/master/util/install/MESS_environment.yaml`
* `conda env create -f MESS_environment.yaml`

# Documentation
Primary MESS documentation lives in the [MESS readthedocs page](https://pymess.readthedocs.io/en/latest/)

### Building the MESS conda package
This is a "note to self", more than end user docs, so don't do this unless you
really want to build the conda package yourself.
* Install conda
* `conda install conda-build anaconda-client`
* `conda build conda.recipe/MESS -c conda-forge -c mess`

