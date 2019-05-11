.. include:: global.rst

.. _parameters:

Community Assembly Parameters
====================
The parameters contained in a params file affect the behavior of various parts
of the forward-time and backward-time assembly process. The defaults that we 
chose are fairly reasonable values as as startin point, however, you will 
always need to modify at least a few of them (for example, to indicate the 
location of your data), and often times you will want to modify many of the 
parameters.

Below is an explanation of each parameter setting, the eco-evolutionary process
that it effects, and example entries for the parameter into a params.txt file.


.. _simulation_name:

0. simulation_name
-----------------
The Assembly name is used as the prefix for all output files. It should be a
unique identifier for the assembly, meaning the set of parameters you are using
for the current data set. When I assemble multiple data with different parameter
combinations I usually either name them consecutively (e.g., data1, data2), or
with names indicating their parameter combinations (e.g., data_clust90,
data_clust85). The Assembly name cannot be changed after an Assembly is created
with the ``-n`` flag, but a new Assembly with a different name can be created
by branching the Assembly (see :ref:`branching workflow<branching_workflow>`).

Affected steps: 1-7
Example: new Assemblies are created with the -n or -b options to MESS:

.. code-block:: bash

    >>> MESS -n data1                       ## create a new assembly named data1


.. _project_dir:

1. project_dir
--------------
A project directory can be used to group together multiple related Assemblies.
A new directory will be created at the given path if it does not already exist.
A good name for Project_dir will generally be the name of the organism being studied.
The project dir path should generally not be changed after an analysis is initiated,
unless the entire directory is moved to a different location/machine.

Affected steps: 1-7
Example entries into params.txt:

.. code-block:: bash

    /home/watdo/MESS/galapagos         ## [1] create/use project dir called galapagos
    galapagos                          ## [1] create/use project dir called galapagos


.. _generations:

2. generations
-----------------
This is a path to the location of raw (non-demultiplexed) fastq data files. If
your data are already demultiplexed then this should be left blank. The input
files can be gzip compressed (i.e., have name-endings with .gz). If you enter
a path for raw data files then you should also enter a path to a
:ref:`barcodes file<barcodes file>`.
To select multiple files, or all files in a directory, use a wildcard character (*).

Affected steps = 1
Example entries into params.txt:

.. code-block:: bash

    /home/deren/MESS/tests/data/*.fastq.gz     ## [2] select all gzip data files
    ~/MESS/tests/data/*.fastq.gz               ## [2] select all gzip data files
    ./ipsimdata/rad_example*.fastq.gz            ## [2] select files w/ `rad_example` in name


.. _community_assembly_model:

3. community_assembly_model
----------------
This is a path to the location of a barcodes_file_. This is used in step1
for demuliplexing, and can also be used in step2 to improve the detection of
adapter/primer sequences that should be filtered out. If your data are already
demultiplexed the barcodes path can be left blank.

Affected steps = 1-2.
Example entries into params.txt:

.. code-block:: python

    /home/deren/ipsimdata/rad_example_barcodes.txt    ## [3] select barcode file
    ./ipsimdata/rad_example_barcodes.txt              ## [3] select barcode file


.. _speciation_model:

4. speciation_model
--------------------
This is a path to the location of sorted fastq data. If your data are already
demultiplexed then this is the location from which data will be loaded when
you run step 1. A wildcard character can be used to select multiple
files in directory.

Affected steps = 1
Example entries into params.txt:

.. code-block:: python

    /home/deren/MESS/tests/ipsimdata/*.fastq.gz    ## [4] select all gzip data files
    ~/MESS/tests/ipsimdata/*.fastq                 ## [4] select all fastq data files
    ./ipsimdata/rad_example*.fastq.gz                ## [4] select files w/ `rad_example` in name


.. _mutation_rate:

5. mutation_rate
--------------------
There are four :ref:`Assembly_methods<Assembly_methods>` options in MESS:
denovo, reference, denovo+reference, and denovo-reference.
The latter three all require a reference sequence file (param #6) in fasta
format. See the :ref:`tutorials` for an example.

Affected steps = 3, 6
Example entries into params.txt:

.. code-block:: python

    denovo                            ## [5] denovo assembly
    reference                         ## [5] reference assembly
    denovo+reference                  ## [5] reference addition assembly
    denovo-reference                  ## [5] reference subtraction assembly


.. _alpha:

6. alpha
---------------------
The reference sequence file should be in fasta format. It does
not need to be a complete nuclear genome, but could also be any
other type of data that you wish to map RAD data to; for example
plastome or transcriptome data.

.. code-block:: python

    ~/MESS/tests/ipsimdata/rad_example_genome.fa   ## [6] select fasta file
    ./data/finch_full_genome.fasta                   ## [6] select fasta file


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
