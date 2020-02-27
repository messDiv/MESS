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
-  `Create and parameterize a new MESS Region <#Create-MESS-Region>`__
-  `Run MESS simulations in API mode <#Simulate-MESS-API>`__
-  `Free time to experiment with other example
   datasets <#Example-Datasets>`__

Running and connecting to a Jupyter notebook
--------------------------------------------
For the purposes of this tutorial, all command interactions will take place
inside a jupyter notebook running on your personal computer. For the most
part we will be writing and executing python commands. Jupyter is already
installed as a dependency of MESS, but if you need help getting a server
running see the excellent `Jupyter notebook documentation
<https://jupyter-notebook.readthedocs.io/en/stable/notebook.html#starting-the-notebook-server>`__.

Parallelization
---------------
Like the CLI, the MESS API can make use of all the cores you can throw at it
thanks to integration with the very nice IPyparallel library. It is ``highly``
recommended to take a moment to :ref:`launch an ipcluster instance <mess_parallelization>`.

Download and examine example data
---------------------------------
We will be using as an example dataset of community-scale COI sequences (~500bp)
and densely sampled abundances for the spider community on the island of La
Reunion, an overseas department of France, which is the largest of the Mascarene
islands, located in the Indian Ocean approximately 1000 km from Madagascar. The
data we will be using was collected and published by Emerson et al (2017). For
this exercise, we will just grab and use the formatted data from the MESS
github repository. For further instruction on properly formatting and converting
raw data into MESS-ready format see the `MESS raw data handling page <MESS_process_raw_data.html>`__.

In a new cell in your notebook you can download the Reunion spider data
like this:

.. code:: bash

   !wget https://raw.githubusercontent.com/messDiv/MESS/master/empirical_data/Reunion_spiders/spider.dat

..

   **NB:** The ``!`` prior to the command here indicates that the
   jupyter notebook should interpret this as a bash command executed at
   the command line. This is a handy way of executing bash commands from
   within the notebook environment, rather than returning to a terminal
   window on the cluster. It’s just a nice shortcut.

Now make a new cell and import MESS and pandas (which is a python
structured data library providing the DataFrame class), and read in the
data you just downloaded.

.. code:: python

   %matplotlib inline
   import MESS
   import pandas as pd
   spider_df = pd.read_csv("spider.dat", index_col=0)
   spider_df[:5]

..

   **Special Note:** The ``%matplotlib inline`` command tells jupyter to
   draw the figures actually inside your notebook environment. If you
   don’t put this your figures won’t show up!

   **NB:** Importing pandas as ``pd`` is pretty cannonical. It makes
   typing out pandas commands shorter because you can reference it as
   ``pd`` rather than ``pandas``.

..

   **NB:** The final line in the above command asks python to display
   the first 5 rows of the ``spider_df`` dataframe. It should look like
   this:

.. figure:: images/Reunion_spider_df.png

 ## Create and parameterize a new MESS Region

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

 ## Run MESS simulations in API mode

Now we can run a bunch of simulations given our new parameterization
using the ``run()`` method. Run takes one argument (``nsims``) which
indicates the number of independent community assembly realizations to
perform.

.. code:: python

   reunion.run(sims=1)

::

      Generating 1 simulation(s).
   [####################] 100%  Finished 0 simulations    | 0:00:00 |

..

   We will not do this now, but the ``run`` method can also accept an
   ``ipyclient`` argument for specifying a connection to an ipyparallel
   backend, allowing for massive parallelization. For more info see the
   `MESS parallelization documentation <MESS_parallelization.html>`__.

Since it can take quite some time to run a number of simulations
sufficient for model selection and parameter estimation we will use a
suite of pre-baked simulations I generated ahead of time. Fetch them
with ``wget`` from the compphylo site:

::

   !wget https://compphylo.github.io/Oslo2019/MESS_files/MESS_simulations/SIMOUT.txt
   !wc -l SIMOUT.txt

::

   100%[======================================>] 14,234,338  50.0M/s   in 0.3s    
   2019-08-11 01:25:27 (50.0 MB/s) - "SIMOUT.txt.1" saved [14234338/14234338]
   24440 SIMOUT.txt

..

   **NB:** The ``wc`` command counts the number of lines if you pass it
   the ``-l`` flag. You can see this series of ~25k simulations is about
   14MB.

You can now proceed to the :ref:`MESS Machine Learning Tutorial <ml_inference>`.

References
----------
::

    Emerson, B. C., Casquet, J., López, H., Cardoso, P., Borges, P. A.,
        Mollaret, N., … & Thébaud, C. (2017). A combined field survey and
        molecular identification protocol for comparing forest arthropod
        biodiversity across spatial scales. Molecular ecology resources, 17(4),
        694-707.

