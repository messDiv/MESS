.. _tutorial_api:

MESS (Massive Eco-Evolutionary Synthesis Simulations) - API Tutorial
====================================================================

This is the second part of the MESS tutorial in which we introduce the
API mode using jupyter notebooks. This is meant as a broad introduction but
does make the assumption you've already completed the :ref:`CLI tutorial <tutorial_cli>`.
some of the parameters and terminology. We will use as an example in this
tutorial the spider community data set from La Reunion published by
Emerson et al (2017). However, you can follow along with one of the
`other example
datasets <https://github.com/messDiv/MESS/tree/master/jupyter-notebooks/empirical>`__
if you like, the procedure will be identical although your results will
vary.

-  `Running and connecting to a Jupyter notebook`_
-  `Create and parameterize a new MESS Region`_
-  `Run MESS serial simulations in API mode`_
-  `Run MESS parallel simulations in API mode`_
-  `Adding more simulations to a region`_
-  `Using simulations to perform statistical inference`_
-  `References`_

Running and connecting to a Jupyter notebook
--------------------------------------------
For the purposes of this tutorial, all command interactions will take place
inside a jupyter notebook running on your personal computer. For the most
part we will be writing and executing python commands. Jupyter is already
installed as a dependency of MESS, but if you need help getting a server
running see the excellent `Jupyter notebook documentation
<https://jupyter-notebook.readthedocs.io/en/stable/notebook.html#starting-the-notebook-server>`__.


Create and parameterize a new MESS Region
-----------------------------------------

MESS API mode lets you dive under the hood of the CLI mode a bit. You
have all the power of the CLI mode, yet more flexibility. The first step
in API mode is to create a MESS ``Region``. A ``Region`` encompasses a
very large Metacommunity and a much smaller local community which is
connected to it by colonization. In creating a ``Region`` the only thing
you’re required to pass in is a name, so lets go with “LaReunion” (No
spaces!), as this is the region the empirical data is drawn from.

.. code:: python

   reunion = MESS.Region("LaReunion")
   print(reunion.get_params())

::

   ------- MESS params file (v.0.1.0)----------------------------------------------
   LaReunion            ## [0] [simulation_name]: The name of this simulation scenario
   ./default_MESS       ## [1] [project_dir]: Where to save files
   0                    ## [2] [generations]: Duration of simulations. Values/ranges Int for generations, or float [0-1] for lambda.
   neutral              ## [3] [community_assembly_model]: Model of Community Assembly: neutral, filtering, competition
   point_mutation       ## [4] [speciation_model]: Type of speciation process: none, point_mutation, protracted, random_fission
   2.2e-08              ## [5] [mutation_rate]: Mutation rate scaled per base per generation
   2000                 ## [6] [alpha]: Abundance/Ne scaling factor
   570                  ## [7] [sequence_length]: Length in bases of the sequence to simulate
   ------- Metacommunity params: --------------------------------------------------
   100                  ## [0] [S_m]: Number of species in the regional pool
   750000               ## [1] [J_m]: Total # of individuals in the regional pool
   2                    ## [2] [speciation_rate]: Speciation rate of metacommunity
   0.7                  ## [3] [death_proportion]: Proportion of speciation rate to be extinction rate
   2                    ## [4] [trait_rate_meta]: Trait evolution rate parameter for metacommunity
   1                    ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change
   ------- LocalCommunity params: Loc1---------------------------------------------
   Loc1                 ## [0] [name]: Local community name
   1000                 ## [1] [J]: Number of individuals in the local community
   0.01                 ## [2] [m]: Migration rate into local community
   0                    ## [3] [speciation_prob]: Probability of speciation per timestep in local community

These are all the parameters of the model. The defaults are chosen to
reflect a typical oceanic island arthropod community. Don’t worry at
this point about all the parameters, lets focus for now on
``community_assembly_model``, the size of the local community (``J``),
and the rate of migration from the metacommunity to the local community
(``m``). We will set parameter ranges for these, and each simulation
will sample a random value from this range. In a new cell use the
``set_param()`` method to change these values:

.. code:: python

   reunion.set_param("community_assembly_model", "*")
   reunion.set_param("J", "1000-10000")
   reunion.set_param("m", "0.001-0.01")

..

   **NB:** Setting the ``community_assembly_model`` to ``*`` indicates
   that we want to sample uniformly among all three of the model types:
   neutral, competition, and environmental filtering.

Print the params again to prove to yourself that the ranges are now set:

.. code:: python

   print(reunion.get_params())

::

   ------- MESS params file (v.0.1.0)----------------------------------------------
   LaReunion            ## [0] [simulation_name]: The name of this simulation scenario
   ./default_MESS       ## [1] [project_dir]: Where to save files
   0                    ## [2] [generations]: Duration of simulations. Values/ranges Int for generations, or float [0-1] for lambda.
   *                    ## [3] [community_assembly_model]: Model of Community Assembly: neutral, filtering, competition
   point_mutation       ## [4] [speciation_model]: Type of speciation process: none, point_mutation, protracted, random_fission
   2.2e-08              ## [5] [mutation_rate]: Mutation rate scaled per base per generation
   2000                 ## [6] [alpha]: Abundance/Ne scaling factor
   570                  ## [7] [sequence_length]: Length in bases of the sequence to simulate
   ------- Metacommunity params: --------------------------------------------------
   100                  ## [0] [S_m]: Number of species in the regional pool
   750000               ## [1] [J_m]: Total # of individuals in the regional pool
   2                    ## [2] [speciation_rate]: Speciation rate of metacommunity
   0.7                  ## [3] [death_proportion]: Proportion of speciation rate to be extinction rate
   2                    ## [4] [trait_rate_meta]: Trait evolution rate parameter for metacommunity
   1                    ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change
   ------- LocalCommunity params: Loc1---------------------------------------------
   Loc1                 ## [0] [name]: Local community name
   1000-10000           ## [1] [J]: Number of individuals in the local community
   0.001-0.01           ## [2] [m]: Migration rate into local community
   0                    ## [3] [speciation_prob]: Probability of speciation per timestep in local community

Run MESS serial simulations in API mode
---------------------------------------

Now we can run community assembly simulations given our new parameterization
using the ``run()`` method. There is one required argument to this method
(``nsims``) which indicates the number of independent community assembly
realizations to perform.

.. code:: python

   reunion.run(sims=1)

::

      Generating 1 simulation(s).
   [####################] 100%  Finished 0 simulations    | 0:01:02 |


Run MESS parallel simulations in API mode
-----------------------------------------

Like the CLI, the MESS API can make use of all the cores you can throw at it
thanks to integration with the very nice :ref:`IPyparallel library <https://ipyparallel.readthedocs.io/>`.
To take a moment to :ref:`launch an ipcluster instance <mess_parallelization>`.

Now we assume you have an ``ipyclient`` object initialized in your notebook.
The ``Region.run()`` method can also an optional argument (``ipyclient``) for
specifying a connection to an ipyparallel backend, allowing for massive
parallelization. Let's check to make sure how many cores our ipyclient is
attached to:

.. code:: python

    len(ipyclient)

::

    40

Now call run and generate 40 simulations on the 40 cores:

.. code:: python

    reunion.run(sims=40, ipyclient=ipyclient)

::

      Generating 40 simulation(s).
    [####################] 100%  Performing Simulations    | 0:01:31 | 
    [####################] 100% 
      Finished 40 simulations

Now we generated 40 simulations in the parallel in the (roughly) the same time
it took to generate 1 simulation in serial. I say 'roughly' here for two
reasons. First, The simulations are stochastic, and the amount of time any
given simulation will take is Poisson distributed, so sometimes you'l get
'unlucky' with one that takes much longer. Second, by default the
``generations`` parameter is ``0``, which indicates to uniformly sample a
``lambda`` value from the half-open interval [0-1). Small ``lambda`` will
(on average) run faster than large ``lambda``, so again, another source of
variability in runtime.

Adding more simulations to a region
-----------------------------------

As with the CLI mode, if you find you need to add more simulations to a
``Region``, for whatever reason, you can simply call ``run()`` again, and this
will append the new simulations to what has already been run.

Using simulations to perform statistical inference
--------------------------------------------------

You can now proceed to the :ref:`MESS Machine Learning Tutorial <ml_inference>`
to learn how to use the simulations to perform model selection and parameter
estimation on real data.

References
----------
::

    Emerson, B. C., Casquet, J., López, H., Cardoso, P., Borges, P. A.,
        Mollaret, N., … & Thébaud, C. (2017). A combined field survey and
        molecular identification protocol for comparing forest arthropod
        biodiversity across spatial scales. Molecular ecology resources, 17(4),
        694-707.

