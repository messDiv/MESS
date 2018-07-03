.. _sec-introduction:

==================
Intro API Tutorial
==================

**MESS - Massive Eco-Evolutionary Synthesis Simulations**

Example API::

  import MESS

  ## A mess region contains all universal parameters of the model, metacommunity
  ## information, and one or more local communities that can be connected
  data = MESS.Region("my_first_sim")
  
  ## Define the metacommunity model
  data.set_metacommunity("logseries")
  
  ## Add local communities to the region
  loc1 = data.add_local_community("Island1", K=5000, c=0.01)
  loc2 = data.add_local_community("Island2", K=1000, c=0.05)
  
  ## Define the potentially asymmetric migration matrix
  ## Migration matrix should be square with dimension equal to # of local communities and,
  ## diagonal elements == 0
  data.migration_matrix([[0, 0.05], [0.05, 0]])
  
  ## Run the simulation for x number of generations
  data.simulate(nsims=100000)
  ## Alternatively simulate until some proportion of equilibrium is reached
  results = data.simulate(lambda=0.7)
  
  print(results)
  
The code is on `github <https://github.com/messDiv/MESS>`_
