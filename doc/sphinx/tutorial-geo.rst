Towards geomechanics
=====================

.. seealso::

	Examples :ref:`gravity-deposition`, :ref:`oedometric-test`, :ref:`periodic-simple-shear`, :ref:`periodic-triaxial-test` deal with topics discussed here.

Parametric studies
-------------------

Input parameters of the simulation (such as size distribution, damping, various contact parameters, …) influence the results, but frequently an analytical relationship is not known. To study such influence, similar simulations differing only in a few parameters can be run and results compared. Yade can be run in *batch mode*, where one simulation script is used in conjunction with *parameter table*, which specifies parameter values for each run of the script. Batch simulation are run non-interactively, i.e. without user intervention; the user must therefore start and stop the simulation explicitly.

Suppose we want to study the influence of :yref:`damping<NewtonIntegrator.damping>` on the evolution of kinetic energy. The script has to be adapted at several places:

#. We have to make sure the script reads relevant parameters from the *parameter table*. This is done using :yref:`yade.utils.readParamsFromTable`; the parameters which are read are created as variables in the ``yade.params.table`` module::

   	utils.readParamsFromTable(damping=.2)      # yade.params.table.damping variable will be created
   	from yade.params import table              # typing table.damping is easier than yade.params.table.damping

   Note that :yref:`yade.utils.readParamsFromTable` takes default values of its parameters, which are used if the script is not run in non-batch mode.

#. Parameters from the table are used at appropriate places::

      NewtonIntegrator(damping=table.damping),

#. The simulation is run non-interactively; we must therefore specify at which point it should stop::

      O.engines+=[PyRunner(iterPeriod=1000,command='checkUnbalancedForce()')]    # call our function defined below periodically

      def checkUnbalancedForce():
         if utils.unbalancedForce<0.05:                      # exit Yade if unbalanced force drops below 0.05
            utils.saveDataTxt(O.tags['d.id']+'.data.bz2')    # save all data into a unique file before exiting   
            import sys
            sys.exit(0)                                      # exit the program

#. Finally, we must start the simulation at the very end of the script::

      O.run()               # run forever, until stopped by checkUnbalancedForce()
      utils.waitIfBatch()   # do not finish the script until the simulation ends; does nothing in non-batch mode

The *parameter table* is a simple text-file, where each line specifies a simulation to run::

	# comments start with # as in python
	damping     # first non-comment line is variable name
	.2
	.4
	.6

Finally, the simulation is run using the special batch command::

   user@machine:~\$ yade-batch parameters.table simulation.py

.. todo:: Parametric studies need to be described better. Perhaps the behavior should be changed so that in batch mode, :yref:`O.run<Omega.run>` and :yref:`yade.utils.waitIfBatch` are run from the main yade script by default. That would however make it only possible to exit the batch via ``sys.exit(0)``, or by a new function like ``utils.exitBatch()`` (which would call ``sys.exit(0)`` for now anyway)

.. rubric:: Exercises

#. Run the gravity deposition script in batch mode, varying :yref:`damping<NewtonIntegrator.damping>` to take values of ``.2``, ``.4``, ``.6``. See the http://localhost:9080 overview page while the batch is running.


Boundary
---------

Particles moving in infinite space usually need some constraints to make the simulation meaningful.

Supports
^^^^^^^^^

So far, supports (unmovable particles) were providing necessary boundary: in the gravity deposition script, :yref:`yade.utils.facetBox` is by internally composed of :yref:`facets <yade.utils.facet>` (triangulation elements), which is ``fixed`` in space; facets are also used for arbitrary triangulated surfaces (see relevant sections of the *User's manual*). Another frequently used boundary is :yref:`yade.utils.wall` (infinite axis-aligned plane).

Periodic
^^^^^^^^^

Periodic boundary is a "boundary" created by using periodic (rather than infinite) space. Such boundary is activated by :yref:`O.periodic=True <Omega.periodic>`, and the space configuration is decribed by :yref:`O.cell <Omega.cell>`. It is well suited for studying bulk material behavior, as boundary effects are avoided, leading to smaller number of particles. On the other hand, it might not be suitable for studying localization, as any cell-level effects (such as shear bands) have to satisfy periodicity as well.

The periodic cell is described by its :yref:`reference size <Cell.refSize>` of box aligned with global axes, and :yref:`current transformation<Cell.trsf>`, which can capture stretch, shear and rotation. Deformation is prescribed via :yref:`velocity gradient<Cell.velGrad>`, which updates the :yref:`transformation<Cell.trsf>` before the next step. :yref:`Homothetic deformation<Cell.homoDeform>` can smear :yref:`velocity gradient<Cell.velGrad>` accross the cell, making the boundary dissolve in the whole cell.
 
Stress and strains can be controlled with :yref:`PeriTriaxController`; it is possible to prescribe mixed strain/stress :yref:`goal<PeriTriaxController.goal>` state using :yref:`PeriTriaxController.stressMask`.

The following creates periodic cloud of spheres and compresses to achieve $\sigma_x$=-10 kPa, $\sigma_y$=-10kPa and $\eps_z$=-0.1. Since stress is specified for $y$ and $z$, :yref:`stressMask<PeriTriaxController.stressMask>` is ``0b011`` (x→1, y→2, z→4, in decimal 1+2=3).

.. ipython::

	Yade [1]: from yade import pack
	
	Yade [1]: sp=pack.SpherePack()

	Yade [1]: sp.makeCloud((1,1,1),(2,2,2),rMean=.2,periodic=True)

	Yade [1]: sp.toSimulation()             # implicitly sets O.periodic=True, and O.cell.refSize to the packing period size

	Yade [1]: O.engines+=[PeriTriaxController(goal=(-1e4,-1e4,-.1),stressMask=0b011,maxUnbalanced=.2,doneHook='functionToRunWhenFinished()')]

When the simulation :yref:`runs<Omega.run>`, :yref:`PeriTriaxController` takes over the control and calls :yref:`doneHook<PeriTriaxController.doneHook>` when :yref:`goal<PeriTriaxController.goal>` is reached. A full simulation with PeriTriaxController might look like the following:

.. literalinclude:: tutorial/peri-triax.py
