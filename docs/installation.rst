.. _sec-installation:

============
Installation
============

``MESS`` requires Python >= 3.5. Binaries can be downloaded with `conda <https://conda.io/docs/>`_ or built with `pip <https://pip.readthedocs.io/en/stable/>`_.

---------------------------------------------
Method 1: Install pre-built binary from conda
---------------------------------------------

1. Download `anaconda <https://www.anaconda.com/download/>`_ or `miniconda <https://conda.io/miniconda.html>`_.
2. (Optional) create a separate `conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_ to install into:

.. code:: bash

    conda create -n mess-env
    source activate mess-env

3. Install:

.. code:: bash

    conda install mess -c messDiv -c conda-forge

-------------------------------------
Method 2: Build from source using pip
-------------------------------------

Prerequisites:

* pip
* C compiler with OpenMP support

Clone the git repository, then ``pip install`` the project root:

.. code:: bash

    git clone https://github.com/messDiv/MESS
    cd MESS/
    pip install .

Depending on your system, ``pip`` may have trouble installing some
dependencies (such as ``numpy``, ``msprime``, ``pysam``).
In this case, you should manually install these dependencies and try again.

See  `venv <https://docs.python.org/3/tutorial/venv.html>`_ to install into a virtual environment.
