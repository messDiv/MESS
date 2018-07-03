.. _sec-introduction:

==================
Intro CLI Tutorial
==================

**MESS - Massive Eco-Evolutionary Synthesis Simulations**

Example CLI::

  ## Create a new params file and populate with default values
  MESS -n neutral_model

  <edit params-neutral_model.txt>

  ## Generate 10000 simulations 
  MESS -p params-neutral_model.txt -s 10000

  ## Validate format of empirical data from a directory (proper formatting will be checked)
  MESS -i empirical_dir

  ## Perform abc model selection (competition model also previously simulated)
  ## Both neutral and competition models should have been simulated for a similar
  ## number of replicates, MESS will check for this.
  MESS -i empirical_dir --abc params-neutral_model.txt params-competition_model.txt

  ## Perform random forest model selection (competition model also previously simulated)
  MESS -i empirical_dir --RF params-neutral_model.txt params-competition_model.txt

  ## Estimate parameters of empirical data for a given model
  ## TODO: Figure out how to specify which parameters to estimate?
  MESS -i empirical_dir -p params-neutral_model.txt --estimate_params

  ## Generate fancy plots through time for a given model. This will
  ## only run one realization and create several animated gifs <slow>
  MESS -p params-neutral_model.txt --fancy-plots

The code is on `github <https://github.com/messDiv/MESS>`_
