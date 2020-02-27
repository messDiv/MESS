.. _mess_parallelization:

MESS Parallelization Docs
-------------------------

Either run this by hand on a compute node or add this line to your
jupyter notebook job submission script before the call to launch the
notebook server:

::

   ipcluster start -n 40 --cluster-id=MESS --daemonize

Here ``-n 40`` specifies the number of cores to use, so you’ll need to
change this to allocate the number you request in the job script.

Now inside the notebook you want to run simulations in, attach to the
running ipyparallel instance like this:

::

   try:
       ipyclient = ipp.Client(cluster_id="MESS")
       print(len(ipyclient))
   except:
       pass

And now you can call use the parallel backend when running simulations
inside a jupyter notebook:

::

   import MESS
   r = MESS.Region("Watdo")
   r.run(nsims=1000, ipyclient=ipyclient)
