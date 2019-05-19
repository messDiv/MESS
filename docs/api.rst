=================
API Documentation
=================

This is the API documentation for ``MESS``, and provides detailed information
on the Python programming interface. See the :ref:`sec_api_tutorial` for an
introduction to using this API to run simulations.

Simulation model
****************

Region
++++++

.. autoclass:: MESS.Region
    :members:

Metacommunity
+++++++++++++

.. autoclass:: MESS.LocalCommunity
    :members:

Local Community
+++++++++++++++

.. autoclass:: MESS.Metacommunity
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

Posterior Predictive Checks
+++++++++++++++++++++++++++

.. autofunction:: MESS.inference.posterior_predictive_check
