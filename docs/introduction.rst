.. _sec-introduction:

============
Introduction
============

**MESS - Massive Eco-Evolutionary Synthesis Simulations**

Example API::

  import MESS

  ## A mess region contains all universal parameters of the model, metacommunity
  ## information, and one or more local communities that can be connected
  data = MESS.Region("my_first_sim")
  
  ## Define the metacommunity model
  data.set_metacommunity("logseries")
  
  ## Add local communities to the region
  loc1 = data.add_local_community("Island1", K=5000, distance=4500, age=1e6)
  loc2 = data.add_local_community("Island2", K=1000, distance=5000, age=5e5)
  
  ## Define the potentially asymmetric migration matrix
  ## Migration matrix should be square with dimension equal to # of local communities and,
  ## diagonal elements == 0
  data.migration_matrix([[0, 0.05], [0.05, 0]])
  
  ## Run the simulation for x number of generations
  data.step(generations=100000)
  ## Alternatively simulate until some proportion of equilibrium is reached
  results = data.simulate(lambda=0.7)
  
  print(results)
  
Example CLI::

  ## Create a new params file and populate with default values
  MESS -n neutral_model

  <edit params-neutral_model.txt>

  ## Generate 10000 simulations 
  MESS -p params-neutral_model.txt -n 10000

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
