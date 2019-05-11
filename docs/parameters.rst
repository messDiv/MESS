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
-----------------------
The restriction overhang is used during demultiplexing (step1) and also to detect
and filter out adapters/primers (in step2), if the `filter_adapters` parameter
is turned on. Identifying the correct sequence to enter for the restriction_overhang
can be tricky. You do not enter the restriction recognition sequence, but rather
the portion of this sequence that is left attached to the sequenced read after
digestion. For example, the enzyme PstI has the following palindromic sequence,
with `^` indicating the cut position.

.. code-block:: python

    5'...C TGCA^G...'3
    3'...G^ACGT C...'5

Digestion with this enzyme results in DNA fragments with the sequence ``CTGCA``
adjacent to the cut site, which when sequenced results in the reverse complement
``TGCAG`` as the restriction overhang at the beginning of each read.
The easiest way to identify the restriction overhang is simply to look at
the raw (or demultiplexed) data files yourself. The restriction overhang
will be the (mostly) invariant sequence that occurs at
the very beginning of each read if the data are already demultiplexed, or right
after the barcode in the case of non-demultplexed data. Use the command
below to peek at the first few lines of your fastQ files to find the invariant
sequence.

.. code-block:: bash

    ## gunzip decompresses the file,
    ## the flag -c means print to screen,
    ## and `less` tells it to only print the first 100 lines
    zless 100 my_R1_input_file.fastq.gz

This will print something like the following. You can see that each of the
lines of sequence data begins with `TGCAG` followed by variable sequence data.
For data that used two-cutters (ddrad), you will likely not be able to see the second
cutter overhang for single-end reads, but if your data are paired-end, then
the `_R2_` files will begin with the second `restriction_overhang`. The second
restriction_overhang is only used to detect adapters/primers if the
`filter_adapters`_ parameter is set > 1. The second `restriction_overhang`
can optionally be left blank.

.. parsed-literal::

    @HWI-ST609:152:D0RDLACXX:2:2202:18249:93964 1:N:0:
    TGCAGCAGCAAGTGCTATTCGTACAGTCATCGATCAGGGTATGCAACGAGCAGAAGTCATGATAAAGGGTCCCGGTCTAGGAAGAGACGCAGCATTA
    +
    BDFFHHHHHJJJHIJJJJJIJIGJJJJJJJJJJJJJJJJDGGIJJJJJHIIIJJJJHIIHIGIGHHHHFFFFEDDDDDACCDCDDDDDDDDDBBDC:
    @HWI-ST609:152:D0RDLACXX:2:2202:18428:93944 1:N:0:
    TGCAGGATATATAAAGAATATACCAATCCTAAGGATCCATAGATTTAATTGTGGATCCAACAATAGAAACATCGGCTCAACCCTTTTAGTAAAAGAT
    +
    ADEFGHFHGIJJJJJJIJJIIJJJIJJIJGIJJJJJJJJIJJIJJJJIIIGGHIEGHJJJJJJG@@CG@@DDHHEFF>?A@;>CCCDC?CDDCCDDC
    @HWI-ST609:152:D0RDLACXX:2:2202:18489:93970 1:N:0:
    TGCAGCCATTATGTGGCATAGGGGTTACATCCCGTACAAAAGTTAATAGTATACCACTCCTACGAATAGCTCGTAATGCTGCGTCTCTTCCTAGACC
    +
    BDFFHHHHHJJJJIJJIJJJJJJJHIJJJJJJJJHIIJJJJIFHIJJJJFGIIJFHIJJIJJIFBHHFFDFEBACEDCDDDDBBBDCCCDDDCDDC:
    @HWI-ST609:152:D0RDLACXX:2:2202:18714:93960 1:N:0:
    TGCAGCATCTGGAAATTATGGGGTTATTTCACAGAAGCTGGAATCTCTTGGGCAATTTCACAGAATCTGGGAATATCTGGGGTAAATCTGCAAGATC
    +
    BDEFHHHHHJJJIJJJJJJJJJJCGIIJJGHJJJJJJJJJJIJJIJJJIHIJJJJJJJHHIIJJJJJJJHGHHHFEFFFDEDABDDFDDEDDDDDDA
    @HWI-ST609:152:D0RDLACXX:2:2202:20484:93905 1:N:0:
    TGCAGAGGGGGATTTTCTGGAGTTCTGAGCATGGACTCGTCCCGTTTGTGCTGCTCGAACACTGACGTTACTTCGTTGATCCCTATGGACTTGGTCA
    +
    ADEEHHHHGJJHIIJJJJJJIJFHJIIIJJJJIIIJJJJEHHIIGHIGGHGHHHHFFEEEDDDD;CBBDDDACC@DC?<CDDCCCCCCA@CCD>>@:
    @HWI-ST609:152:D0RDLACXX:2:2202:20938:93852 1:N:0:
    TGCAGAAGCTGGAGATTCTGGGGCAGCTTTGCAGCAAGCTGAAAATTCTGGGGGTCGATCTGCAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG

In some cases restriction enzymes can bind to more than one specific sequence,
for example ApoI will bind to AATTY (i.e. AATTC and AATTT). If you used an 
enzyme with reduced specificity you can include ambiguity codes in the 
restriction overhang sequence.

Affected steps = 1,2. Example entries to params.txt file:

.. code-block:: bash

    TGCAG                     ## [8] single cutter (e.g., rad, gbs)
    TGCAG, AATT               ## [8] double digest (e.g., ddrad, pairddrad)
    CWGC                      ## [8] single cutter w/ degenerate base


NB: 3RAD and SeqCap data can use up to 4 restriction enzymes. If you have this
kind of data, simply list all the restriction overhangs for all your cutters.

.. code-block:: bash

    CTAGA, CTAGC, AATTC               ## [8] 3rad data (multiple cutters)


.. _J_m:

9. J_m
---------------------
During step 2 bases are trimmed from the 3' end of reads when the quality score
is consistently below 20 (which can be modified by modifying phred_Qscore_offset_). 
However, your reads may still contain some number of ambiguous (N) sites that 
were not trimmed based on quality scores, and these will affect the efficiency
and accuracy of clustering downstream. This parameter sets the upper limit on the 
number of Ns allowed in reads. The default value for
`max_low_qual_bases` is 5. I would generally recommend against increasing 
this value greatly. 

Affected steps = 2. Example entries to params.txt:

.. code-block:: bash

    0                      ## [9] allow zero low quality bases in a read
    5                      ## [9] allow up to five low quality bases in a read


.. _speciation_rate:

10. speciation_rate
------------------------
Bases are trimmed from the 3' end of reads if their quality scores is below
this 20. The default offset for quality scores is 33. Some 
older data use a qscore offset of 64, but this is increasingly rare. You
can toggle the offset number to change the threshold for trimming. 
For example, reducing the offset from 33 to 23 is equivalent to changing the 
minimum quality score from 20 to 10, which is approximately 
95% probability of a correct base call.

Affected steps = 2. Example entries to params.txt:

.. parsed-literal::

    33                 ## [10] default offset of 33, converts to min score=20
    43                 ## [10] offset increased by 10, converts to min score=30
    64                 ## [10] offset used by older data, converts to min score=20.


.. _death_proportion:

11. death_proportion
-------------------------
This is the minimum depth at which statistical base calls will be made during
step 5 consensus base calling. By default this is set to 6, which for most
reasonable error rates estimates is approximately the minimum depth at which a
heterozygous base call can be distinguished from a sequencing error.

Affected steps = 4, 5. Example entries to params.txt

.. parsed-literal::

    6                 ## [11] set mindepth statistical to 6
    10                ## [11] set to 10


.. _trait_rate_meta:

12. trait_rate_meta
---------------------
This is the minimum depth at which majority rule base calls are made during
step 5 consensus base calling. By default this is set to the same value as
mindepth_statistical, such that only statistical base calls are made. This
value must be <= mindepth_statistical. If lower, then sites with coverage
>= mindepth_majrule and < mindepth_statistical will make majority rule calls.
If your data set is very low coverage such that many clusters are excluded due
to low sequencing depth then lowering mindepth_majrule can be an effective way
to increase the amount of usable information in your data set. However, you
should be aware the majority rule consensus base calls will underestimate
heterozygosity.

Affected steps = 4, 5. Example entries to params.txt:

.. parsed-literal::

    6                 ## [12] set to relatively high value similar to mindepth_stat
    2                 ## [12] set below the statistical limit for base calls.


.. _ecological_strength:
13. ecological_strength
-----------------------
This parameter dictates the strength of interactions in the environmental
filtering and competition models.

![ecological_strength_0.001](images/ecological_strength_0.001.png)

Sequencing coverage is often highly uneven among due to differences in the
rate at which fragments are amplified during library preparation, the extent
to which varies across different library prep methods. Moreover, repetitive
regions of the genome may appear highly similar and thus cluster as high depth
clusters. Setting a maxdepth helps to remove the latter problem, but at the
expense of potentially removing good clusters that simply were sequenced
to high depth. The default maxdepth is set quite high (10,000), but you may
change it as you see fit.

Example entries to params.txt:

.. parsed-literal::

    1       ## [5] [ecological_strength]: Strength of community assembly process on phenotypic change


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
