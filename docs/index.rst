.. MESS documentation master file, created by
   sphinx-quickstart on Tue Apr  3 09:55:21 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Massive Eco-evolutionary Synthesis Simulations (MESS)
=====================================================

MESS is a novel comparative phylogeographic model grounded in community
ecological theory. This integrative approach makes use of four data axes
(distributions of traits, abundances, genetic diversities/divergences, and
phylogenetic patterns) to enable testing alternative community assembly
models (neutral vs non-neutral) and estimating parameters underlying different
assembly processes (e.g. dispersal vs in situ speciation). This method
capitalizes on the widespread use of DNA barcoding and meta-barcoding
approaches.

What kind of data it requires
-----------------------------
MESS requires population-level sampling (5-10 individuals per species) from a local community or multiple local communities. This can be at a variety of scales ranging from a microbial community within a host individual, a locally sampled plot targeting everything from a taxonomic group, to a regional assemblage that emerged via disersal and/or local speciation. Currently only single locus data is supported, so community metabarcoding projects would be quite appropriate. Other data types can be included but are not required (abundances, per taxon trait metrics, and phylogenetic information).

Software
--------
MESS is implemented in python and the code is available on `github <https://github.com/messDiv/MESS>`_.

Try it now!
-----------
Experiment with the MESS model now! Launch the binder instance below and you
can open and run the notebooks in the `jupyter-notebooks` directory.

.. image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/messDiv/MESS/master

References
----------
The MESS manuscript is available on `bioarxiv <https://www.biorxiv.org/content/biorxiv/early/2020/01/31/2020.01.30.927236.full.pdf>`_::

    Overcast I, Ruffley M, Rosindell J, Harmon L, Borges PA, Emerson BC, ... &
        Rominger A. (2020). A unified model of species abundance, genetic
        diversity, and functional diversity reveals the mechanisms structuring
        ecological communities. BioRxiv.

MESS is based on previous work on the gimmeSAD joint neutral model which can be
found here::

    Overcast I, Emerson BC, Hickerson MJ. (2019). An integrated model of
        population genetics and community ecology. Journal of Biogeography,
        46: 816-829. https://doi.org/10.1111/jbi.13541

.. toctree::
   :maxdepth: 1

   installation
   tutorial_cli
   tutorial_api
   inference
   parameters
   api


