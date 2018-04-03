.. _sec-introduction:

============
Introduction
============

**MESS - Massive Eco-Evolutionary Synthesis Simulations**

import MESS

## A mess region contains all universal parameters of the model, metacommunity
## information, and one or more local communities that can be connected
mess = MESS.Region("my_first_sim")

## Define the metacommunity model
mess.set_metacommunity("logseries")

## Add local communities to the region
loc1 = mess.add_local_community("Island1", K=5000, distance=4500, age=1e6)
loc2 = mess.add_local_community("Island2", K=1000, distance=5000, age=5e5)

## Define the potentially asymmetric migration matrix
## Migration matrix should be square with dimension equal to # of local communities and,
## diagonal elements == 0
mess.migration_matrix([[0, 0.05], [0.05, 0]])

## Run the simulation for x number of generations, or until a certain lambda is reached
mess.step(generations=100000)
mess.step(lambda=1)

## Some functions to get information about local community status
mess.get_local_community("Island1")

The code is on `github <https://github.com/messDiv/MESS>`_
