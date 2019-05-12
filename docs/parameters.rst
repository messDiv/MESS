.. include:: global.rst

.. _parameters:

Community Assembly Parameters
====================
The parameters contained in a params file affect the behavior of various parts
of the forward-time and backward-time assembly process. The defaults that we 
chose are fairly reasonable values as as starting point, however, you will 
always need to modify at least a few of them (for example, to indicate the 
location of your data), and often times you will want to modify many of the 
parameters.

Below is an explanation of each parameter setting, the eco-evolutionary process
that it affects, and example entries for the parameter into a params.txt file.


.. _simulation_name:

0. simulation_name
-----------------
The simulation name is used as the prefix for all output files. It should be a
unique identifier for this particular set of simulations, meaning the set of 
parameters you are using for the current data set. When I run multiple related
simulations I usually use names indicating the specific parameter combinations 
(e.g., filtering_nospeciation, J5000_neutral). 

Example: New community simulations are created with the -n options to MESS:

.. code-block:: bash

    >>> MESS -n J1000_neutral          ## create a new assembly named J1000_neutral


.. _project_dir:

1. project_dir
--------------
A project directory can be used to group together multiple related simulations.
A new directory will be created at the given path if it does not already exist.
A good name for project_dir will generally be the name of the community/system being 
studied. The project dir path should generally not be changed after simulations/analysis
are initiated, unless the entire directory is moved to a different location/machine.

Example entries into params.txt:

.. code-block:: bash

    /home/watdo/MESS/galapagos         ## [1] create/use project dir called galapagos
    galapagos                          ## [1] create/use project dir called galapagos in cwd (./)


.. _generations:

2. generations
--------------
This parameter specifies the amount of time to run forward-time simulations. 
It can be specified in a number of different ways, but overall time can be 
considered either in terms of Wright-Fisher (WF) generations or in terms of Lambda.
For WF generations you should specify an integer value (or a range of integer values)
which will run the forward-time process for WF * J / 2 time-steps (where a time-step
is one birth/death/colonization/speciation event). For Lambda you may select
either an exact Lambda value (a real value between 0 and 1 exclusive), or you
can set `generations` equal to 0, which will draw random Lambda values between
0 and 1 for each simulatoin.

Example entries into params.txt:

.. code-block:: bash

    0                  ## [2] [generations]: Sample random Lambda values for each simulation 
    100                ## [2] [generations]: Run each simulation for 100 WF generations
    50-100             ## [2] [generations]: Sample uniform between 50-100 WF generations for each simulation


.. _community_assembly_model:

3. community_assembly_model
----------------
With this parameter you may specify a neutral or non-neutral scenario for
the forward time process. There are currently three different options for
this parameter: `neutral`, `filtering`, or `competition`. The `neutral`
case indicates full ecological equivalence of all species, so all
individuals have an equal probability of death at each time-step. In the
`filtering` and `competition` models survival probability is contingent
on proximity of species trait values to the environmental optimum, or distance
from the local trait mean, respectively.

Example entries into params.txt:

.. code-block:: bash

    neutral             ## [3] [community_assembly_model]: Select the neutral process forward-time
    filtering           ## [3] [community_assembly_model]: Select the environmental filtering process
    *                   ## [3] [community_assembly_model]: Randomly choose one of the community assembly models


.. _speciation_model:

4. speciation_model
--------------------

Specify a speciation process in the local community. If `none` then no
speciation happens locally. If `point_mutation` then one individual
will instantaneously speciate at rate `speciation_prob` for each forward-time
step. If `random_fission` then one lineage will randomly split into
two lineages at rate `speciation_prob` with the new lineage receiving
`n = U~(1, local species abundance)` individuals, and the parent lineage 
receiving `1 - n` individuals. `protracted` will specify a model of
protracted speciation, but this is as yet unimplemented.

Example entries into params.txt:

.. code-block:: bash

    none                ## [4] [speciation_model]: No speciation in the local community
    point_mutation      ## [4] [speciation_model]: Point mutation specation process


.. _mutation_rate:

5. mutation_rate
----------------
Specify the mutation rate for backward-time coalescent simulation of
genetic variation. This rate is the per base, per genration probability
of a mutation under an infinite sites model.

Example entries into params.txt:

.. code-block:: bash

    2.2e-08             ## [5] [mutation_rate]: Mutation rate scaled per base per generation

.. _alpha:

6. alpha
---------------------


.. code-block:: bash

    2000                 ## [6] [alpha]: Abundance/Ne scaling factor



.. _sequence_length:

7. sequence_length
------------
There are now many forms of restriction-site associated DNA library preparation
methods and thus many differently named data types. Currently, we categorize
these into :ref:`six data types <data_types>`. Follow the link
to deteremine the appropriate category for your data type.


.. code-block:: python

    rad                       ## [7] rad data type (1 cutter, sonication)
    pairddrad                 ## [7] paired ddrad type (2 different cutters)
    pairgbs                   ## [7] paired gbs type (1 cutter cuts both ends)


.. _S_m:

8. S_m
------
S_m specifies the total number of species to simulate in the metacommunity. Larger
values will result in more singletons in the local community and reduced rates
of multiple-colonization.

Example entries to params.txt file:

.. code-block:: bash

    500                     ## [0] [S_m]: Number of species in the regional pool
    100-1000                ## [0] [S_m]: Number of species in the regional pool


.. _J_m:

9. J_m
------
The total number of individuals in the metacommunity.

Example entries to params.txt:

.. code-block:: bash

    0                      ## [9] allow zero low quality bases in a read
    5                      ## [9] allow up to five low quality bases in a read


.. _speciation_rate:

10. speciation_rate
------------------------

Example entries to params.txt:

.. parsed-literal::

    33                 ## [10] default offset of 33, converts to min score=20


.. _death_proportion:

11. death_proportion
--------------------
This is the minimum depth at which statistical base calls will be made during
step 5 consensus base calling. By default this is set to 6, which for most
reasonable error rates estimates is approximately the minimum depth at which a
heterozygous base call can be distinguished from a sequencing error.

Example entries to params.txt

.. parsed-literal::

    6                 ## [11] set mindepth statistical to 6


.. _trait_rate_meta:

12. trait_rate_meta
---------------------

Example entries to params.txt:

.. parsed-literal::

    6                 ## [12] set to relatively high value similar to mindepth_stat


.. _ecological_strength:
13. ecological_strength
-----------------------
This parameter dictates the strength of interactions in the environmental
filtering and competition models. As the value of this parameter approaches
zero, ecological strength is reduced and the assembly process increasingly
resembles neutrality (ecological equivalence). Larger values increasingly
bias probability of death against individuals with traits farther from 
the environmental optimum (in the filtering model).

In the following examples the environmental optimum is `3.850979`, and the 
ecological strength is varied from 0.001 to 100. Column 0 is species ID,
column 1 is trait value, column 2 is unscaled probability of death, and
column 3 is proportional probability of death. Models with strength of
0.001 and 0.01 are essentially neutral. Strength of 0.1 confers a slight 
advantage to individuals very close to the local optimum (e.g. species 't97').

.. image:: images/ecological_strength_0.001.png
    :width: 25 %
.. image:: images/ecological_strength_0.01.png
    :width: 29 %
.. image:: images/ecological_strength_0.1.png
    :width: 30 %

Ecological strength of 1 (below, left panel) is noticeably non-neutral (e.g. 't97' 
survival probability is 10x greater than average). A value of 10 for this 
parameter generates a _strong_ non-neutral process (below, center panel: 't97' is 100x less 
likely to die than average, and the distribution of death probabilities is
more varied). Ecological strength values >> 10 are _extreme_ and will probably
result in degenerate behavior (e.g. strength of 100 (below, right panel) in which
several of the species will be effectively immortal, with survival probability
thousands of times better than average).

.. image:: images/ecological_strength_1.png
    :width: 30 %
.. image:: images/ecological_strength_10.png
    :width: 30 %
.. image:: images/ecological_strength_100.png
    :width: 30 %

Example entries to params.txt:

.. parsed-literal::

    1             ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change
    0.001-1       ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change


.. _name:

14. name
--------------------
This the level of sequence similarity at which two sequences are identified
as being homologous, and thus cluster together. The value should be entered
as a decimal (e.g., 0.90). We do not recommend using values higher than 0.95,
as homologous sequences may not cluster together at such high threshold due
to the presence of Ns, indels, sequencing errors, or polymorphisms.

Affected steps = 3, 6. Example entries to params.txt:

.. parsed-literal::

    0.90              ## [14] clust threshold set to 90%
    0.85              ## [14] clust threshold set to 85%


.. _J:

15. J
---------------
The maximum number of allowed mismatches between the barcodes in the barcodes
file and those found in the sequenced reads. Default is 0. Barcodes usually differ
by a minimum of 2 bases, so I would not generally recommend using a value >2.

Affected steps = 1. Example entries to params.txt:

.. parsed-literal::

    0              ## [15] allow no mismatches
    1              ## [15] allow 1 mismatched base


.. _m:

16. m
--------------------
It is important to remove Illumina adapters from your data if present. 
Depending on the fidelity of the size selection procedure implemented during
library preparation there is often at least some small proportion of sequences
in which the read length is longer than the actual DNA fragment, such that the
primer/adapter sequence ends up in the read. This occurs more commonly in 
double-digest (GBS, ddRAD) data sets that use a common cutter, and can be 
especially problematic for GBS data sets, in which short fragments are sequenced 
from either end. The `filter_adapters` parameter has three settings (0, 1, or 2). 
**If 0**, then reads are only removed if they contain more Ns than allowed by the
`max_low_qual_bases` parameter. **If 1**, then reads are trimmed to the first base
which has a Qscore < 20 (on either read for paired data), and also removed if there
are too many Ns. **If 2**, then reads are searched for the common Illumina adapter, 
plus the reverse complement of the second cut site (if present), plus the barcode
(if present), and this part of the read is trimmed. This filter is applied using 
code from the software `cutadapt`, which allows for errors within the adapter 
sequence. 

Affected steps = 2. Example entries to params.txt:

.. parsed-literal::

    0                ## [16] No adapter filtering
    1                ## [16] filter based on quality scores
    2                ## [16] strict filter for adapters


.. _speciation_prob:

17. speciation_prob
------------------------
During step 2 if `filter_adapters` is > 0 reads may be trimmed to a shorter length
if they are either low quality or contain Illumina adapter sequences. 
By default MESS will keep trimmed reads down to a minimum length of 35bp. 
If you want to set a higher limit you can do so here.

Affected steps = 2. Example entries to params.txt

.. parsed-literal::

    50                ## [17] minimum trimmed seqlen of 50
    75                ## [17] minimum trimmed seqlen of 75
