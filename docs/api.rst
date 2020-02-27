=================
API Documentation
=================

This is the API documentation for ``MESS``, and provides detailed information
on the Python programming interface. See the :ref:`tutorial_api` for an
introduction to using this API to run simulations.

Simulation model
****************

Region
++++++

.. autoclass:: MESS.Region
    :members:

Metacommunity
+++++++++++++

.. autoclass:: MESS.Metacommunity
    :members:

Local Community
+++++++++++++++

.. autoclass:: MESS.LocalCommunity
    :members:

Inference Procedure
*******************

.. autoclass:: MESS.inference.Ensemble
    :members:

Model Selection (Classification)
++++++++++++++++++++++++++++++++

.. autoclass:: MESS.inference.Classifier
    :members:

Parameter Estimation (Regression)
+++++++++++++++++++++++++++++++++

.. autoclass:: MESS.inference.Regressor
    :members:

Classification Cross-Validation
+++++++++++++++++++++++++++++++

.. autofunction:: MESS.inference.classification_cv

Parameter Estimation Cross-Validation
+++++++++++++++++++++++++++++++++++++

.. autofunction:: MESS.inference.parameter_estimation_cv

Posterior Predictive Checks
+++++++++++++++++++++++++++

.. autofunction:: MESS.inference.posterior_predictive_check

Stats
*****

.. autofunction:: MESS.stats.calculate_sumstats

.. autofunction:: MESS.stats.feature_sets

.. autofunction:: MESS.stats.Watterson

Plotting
++++++++

.. autofunction:: MESS.plotting.plot_simulations_hist

.. autofunction:: MESS.plotting.plot_simulations_pca

.. autofunction:: MESS.plotting.plots_through_time

