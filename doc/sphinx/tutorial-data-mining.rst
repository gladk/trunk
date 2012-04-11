Data mining
=============

Read
-----

Local data
^^^^^^^^^^^

All data of the simulation are accessible from python; when you open the *Inspector*, blue labels of various data can be clicked -- left button for getting to the documentation, middle click to copy the name of the object (use ``Ctrl-V`` or middle-click to paste elsewhere). The interesting objects are among others (see :yref:`Omega` for a full list):

#. :yref:`O.engines<Omega.engines>`
   
   Engines are accessed by their index (position) in the simulation loop::

   	O.engines[0]      # first engine
   	O.engines[-1]     # last engine

   .. note:: The index can change if :yref:`O.engines<Omega.engines>` is modified. *Labeling* introduced below is a better solution for reliable access to a particular engine.

#. :yref:`O.bodies<Omega.bodies>`

   Bodies are identified by their :yref:`id<Body.id>`, which is guaranteed to not change during the whole simulation::

   	O.bodies[0]                                                   # first body
   	[b.shape.radius in O.bodies if isinstance(b.shape,Sphere)]    # list of radii of all spherical bodies
   	sum([b.state.mass for b in O.bodies])                         # sum of masses of all bodies

   .. note:: Uniqueness of :yref:`Body.id` is not guaranteed, since newly created bodies might recycle :yref:`ids<Body.id>` of :yref:`deleted<BodyContainer.erase>` ones.

#. :yref:`O.force<Omega.force>`

   Generalized forces (forces, torques) acting on each particle. They are (usually) reset at the beginning of each step with :yref:`ForceResetter`, subsequently forces from individual interactions are accumulated in :yref:`InteractionLoop`. To access the data, use::

   	O.forces.f(0)     # force on #0
   	O.forces.t(1)     # torque on #1
	
#. :yref:`O.interactions<Omega.interactions>`

   Interactions are identified by :yref:`ids<Body.id>` of the respective interacting particles (they are created and deleted automatically during the simulation)::

   	O.interactions[0,1]   # interactions of #0 with #1
   	O.interactions[1,0]   # the same object
   	O.bodies[0].intrs     # all interactions of body #0

Labels
"""""""

:yref:`Engines<Engine>` and :yref:`functors<Functor>` can be *labeled*, which means that python variable of that name is automatically created.

.. ipython::

	Yade [1]: from yade import *
	
	Yade [1]: O.engines=[
	     ...:    NewtonIntegrator(damping=.2,label='newton')
	     ...: ]
	     ...:

	Yade [1]: newton.damping=.4

	Yade [1]: O.engines[0].damping   # O.engines[0] and newton are the same objects

.. rubric:: Exercises

#. Find meaning of this
   expression::

   	max([b.state.vel.norm() for b in O.bodies])

#. Run the gravity deposition script, pause after a few seconds of simulation. Write expressions that compute

   #. kinetic energy $\sum \frac{1}{2} m_i |v_i| ^2$
   #. average mass (hint: use `numpy.average <http://docs.scipy.org/doc/numpy/reference/generated/numpy.average.html>`__)
   #. maximum $z$-coordinate of all particles
   #. number of interactions of body #1

Global data
^^^^^^^^^^^

Useful measures of what happens in the simulation globally:

unbalanced force
	ratio of maximum contact force and maximum per-body force; measure of staticity, computed with :yref:`yade.utils.unbalancedForce`.
porosity
	ratio of void volume and total volume; computed with :yref:`yade.utils.porosity`.
coordination number
	average number of interactions per particle, :yref:`yade.utils.avgNumInteractions`
stress tensor (periodic boundary conditions)
	averaged force in interactions, computed with :yref:`yade.utils.normalShearStressTensor` and :yref:`yade.utils.stressTensorOfPeriodicCell`
fabric tensor
	distribution of contacts in space (not yet implemented); can be visualized with :yref:`yade.utils.plotDirections`

Energies
""""""""

Evaluating energy data for all components in the simulation (such as gravity work, kinetic energy, plastic dissipation, damping dissipation) can be enabled with ::

	O.trackEnergy=True

Subsequently, energy values are accessible in the :yref:`O.energy<Omega.energy>`; it is a dictionary where its entries can be retrived with ``keys()`` and their values with ``O.energy[key]``.

Save
----

PyRunner
^^^^^^^^^

To save data that we just learned to access, we need to call Python from within the *simulation loop*. :yref:`PyRunner` is created just for that; it inherits periodicy control from :yref:`PeriodicEngine` and takes the code to run as text (must be quoted, i.e. inside ``'...'``) attributed called *command*. For instance, adding this to :yref:`O.engines<Omega.engines>` will print the current step number every second::

	O.engines=O.engines+[ PyRunner(command='print O.iter',realPeriod=1) ]

Writing complicated code inside *command* is awkward; in such case, we define a function that will be called::

	def myFunction():
		'''Print step number, and pause the simulation is unbalanced force is smaller than 0.05.'''
		print O.iter
		if utils.unbalancedForce()<0.05:
			print 'Unbalanced force is smaller than 0.05, pausing.'
			O.pause()
	O.engines=[
		# ...
		PyRunner(command='myFunction()',iterPeriod=100) # call myFunction every 100 steps
	]

.. rubric:: Exercises

#. Run the gravity deposition simulation, but change it such that:

   #. :yref:`yade.utils.unbalancedForce` is printed every 2 seconds.
   #. check every 1000 steps the value of unbalanced force

      * if smaller than 0.2, set :yref:`damping<NewtonIntegrator.damping>` to 0.8 (hint: use labels)
      * if smaller than 0.1, pause the simulation

Keeping history
^^^^^^^^^^^^^^^^^

Yade provides the :yref:`yade.plot` module used for storing and plotting variables (plotting itself will be discussed later). Periodic storing of data is done with :yref:`PyRunner` and the :yref:`yade.plot.addData` function, for instance::

	from yade import plot
	O.engines=[  # ...,
		PyRunner(command='addPlotData()',realPeriod=2)                 # call the addPlotData function every 2 seconds of human time
	]
	def addPlotData():
		# this function adds current values to the history of data, under the names specified
		plot.addData(i=O.iter,t=O.time,Ek=utils.kineticEnergy(),coordNum=utils.avgNumInteractions(),unForce=utils.unbalancedForce())

History is stored in :yref:`yade.plot.data`, and can be accessed using the variable name, e.g. ``plot.data['Ek']``, and saved to text file (for post-processing outside yade) with :yref:`yade.plot.saveTxt`.

Plot
-----

:yref:`yade.plot` provides facilities for plotting history saved with :yref:`yade.plot.addData` as 2d plots. Data to be plotted are specified using dictionary :yref:`yade.plot.plots` ::

	plot.plots={'t':('coordNum','unForce',None,'Ek')}

History of all values is given as the name used for :yref:`yade.plot.addData`; keys of the dictionary are $x$-axis values, and values are sequence of data on the $y$ axis; the ``None`` separates data on the left and right axes (they are scaled independently). The plot itself is created with ::

	plot.plot()         # on the command line, F8 can be used as shorthand

While the plot is open, it will be updated periodically, so that simulation evolution can be seen in real-time.

Energy plots
^^^^^^^^^^^^^

Plotting all energy contributions would be difficult, since names of all energies might not be known in advance. Fortunately, there is a way to handle that in Yade. It consists in two parts:

#. :yref:`yade.plot.addData` is given all the energies that are currently defined::

  		plot.addData(i=O.iter,total=O.energy.total(),**O.energy)

   The :yref:`O.energy.total<EnergyTracker.total>` functions, which sums all energies together. The ``**O.energy`` is special python syntax for converting dictionary (remember that :yref:`O.energy<EnergyTracker>` is a dictionary) to named functions arguments, so that the following two commands are identical::

     function(a=3,b=34)              # give arguments as arguments
     function(**{'a':3,'b':34})      # create arguments from dictionary

#. Data to plot are specified using a *function* that gives names of data to plot, rather than providing the data names directly::

   	plot.plots={'i':['total',O.energy.keys()]}     

   where ``total`` is the name we gave to ``O.energy.total()`` above, while ``O.energy.keys()`` will always return list of currently defined energies.

.. rubric:: Exercises

#. Run the gravity deposition script, plotting unbalanced force and kinetic energy.
#. While the script is running, try changing the :yref:`NewtonIntegrator.damping` parameter (do it from both *Inspector* and from the command-line). What influence does it have on the evolution of unbalanced force and kinetic energy?
#. Think about and write down all energy sources (input); write down also all energy sinks (dissipation).
#. Simulate gravity deposition and plot all energies as they evolve during the simulation.

.. seealso::
	
	Most :ref:`examples` use plotting facilities of Yade, some of them also track energy of the simulation.
