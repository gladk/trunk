.. _hands-on:

Hands-on
========


Shell basics
-------------

Directory tree
^^^^^^^^^^^^^^^

Directory tree is hierarchical way to organize files in operating systems. A typical (reduced) tree looks like this::

   /            Root
   ├──boot        System startup
   ├──bin         Low-level programs
   ├──lib         Low-level libraries
   ├──dev         Hardware access
   ├──sbin        Administration programs
   ├──proc        System information
   ├──var         Files modified by system services
   ├──root        Root (administrator) home directory
   ├──etc         Configuration files
   ├──media       External drives
   ├──tmp         Temporary files
   ├──usr         Everything for normal operation (usr = UNIX system resources)
   │    ├──bin       User programs
   │    ├──sbin      Administration programs
   │    ├──include   Header files for c/c++
   │    ├──lib       Libraries
   │    ├──local     Locally installed software
   │    └──doc       Documentation
   └──home        Contains the user's home directories
        ├──user      Home directory for user
        └──user1     Home directory for user1

Note that there is a single root ``/``; all other disks (such as USB sticks) attach to some point in the tree (e.g. in ``/media``).


Shell navigation
^^^^^^^^^^^^^^^^^

Shell is the UNIX command-line, interface for conversation with the machine. Don't be afraid.


Moving around
"""""""""""""""

The shell is always operated by some ``user``, at some concrete ``machine``; these two are constant. We can move in the directory structure, and the current place where we are is *current directory*. By default, it is the *home directory* which contains all files belonging to the respective user::

	user@machine:~\$                         # user operating at machine, in the directory ~ (= user's home directory)
	user@machine:~\$ ls .                    # list contents of the current directory
	user@machine:~\$ ls foo                  # list contents of directory foo, relative to the dcurrent directory ~ (= ls ~/foo = ls /home/user/foo)
	user@machine:~\$ ls /tmp                 # list contents of /tmp
	user@machine:~\$ cd foo                  # change directory to foo
	user@machine:~/foo\$ ls ~                # list home directory (= ls /home/user)
	user@machine:~/foo\$ cd bar              # change to bar (= cd ~/foo/bar)
	user@machine:~/foo/bar\$ cd ../../foo2   # go to the parent directory twice, then to foo2 (cd ~/foo/bar/../../foo2 = cd ~/foo2 = cd /home/user/foo2)
	user@machine:~/foo2\$ cd                 # go to the home directory (= ls ~ = ls /home/user)
	user@machine:~\$

Users typically have only permissions to write (i.e. modify files) only in their home directory (abbreviated ``~``, usually is ``/home/user``) and ``/tmp``, and permissions to read files in most other parts of the system::

	user@machine:~\$ ls /root    # see what files the administrator has
	ls: cannot open directory /root: Permission denied

Keys
"""""

Useful keys on the command-line are:

============= =========================
<tab>         show possible completions of what is being typed (use abundantly!)
^C (=Ctrl+C)  delete current line
^D            exit the shell
↑↓            move up and down in the command history
^C            interrupt currently running program
^\\           kill currently running program
Shift-PgUp    scroll the screen up (show part output)
Shift-PgDown  scroll the screen down (show future output; works only on quantum computers)
============= =========================

.. Multiple files can be selected using *patterns*.
	user@machine:~/foo2\$ ls *.py            # * replaces any characters (except /)
	params.py  remote.py  timing.py
	user@machine:~/foo2\$ ls ../foo2/*.py
	params.py  remote.py  timing.py

Running programs
"""""""""""""""""

When a program is being run (without giving its full path), several directories are searched for program of that name; those directories are given by ``\$PATH``::

	user@machine:~\$ echo \$PATH     # show the value of $PATH
	/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
	user@machine:~\$ which ls       # say what is the real path of ls

The first part of the command-line is the program to be run (``which``), the remaining parts are *arguments* (``ls`` in this case). It is upt to the program which arguments it understands. Many programs can take special arguments called *options* starting with ``-`` (followed by a single letter) or ``--`` (followed by words); one of the common options is ``-h`` or ``--help``, which displays how to use the program (try ``ls --help``).

Full documentation for each program usually exists as *manual page* (or *man page*), which can be shown using e.g. ``man ls`` (``q`` to exit)

Starting yade
^^^^^^^^^^^^^^

If yade is installed on the machine, it can be (roughly speaking) run as any other program; without any arguments, it runs in the "dialog mode", where a command-line is presented:

::

	user@machine:~\$ yade
	Welcome to Yade bzr2616
	TCP python prompt on localhost:9002, auth cookie `adcusk'
	XMLRPC info provider on http://localhost:21002
	[[ ^L clears screen, ^U kills line. F12 controller, F11 3d view, F10 both, F9 generator, F8 plot. ]]
	Yade [1]:                                            #### hit ^D to exit
	Do you really want to exit ([y]/n)?
	Yade: normal exit.

The command-line is in fact ``python``, enriched with some yade-specific features. (Pure python interpreter can be run with ``python`` or ``ipython`` commands).

Instead of typing commands on-by-one on the command line, they can be be written in a file (with the .py extension) and given as argument to Yade::

	user@machine:~\$ yade simulation.py

For a complete help, see ``man yade``


.. rubric:: Exercises

#. Open the terminal, navigate to your home directory
#. Create a new empty file and save it in ``~/first.py``
#. Change directory to ``/tmp``; delete the file ``~/first.py``
#. Run program ``xeyes``
#. Look at the help of Yade.
#. Look at the *manual page* of Yade
#. Run Yade, exit and run it again.

Python basics
--------------

We assume the reader is familar with `Python tutorial <http://docs.python.org/tutorial/index.html>`__ and only briefly review some of the basic capabilities. The following will run in pure-python interpreter (``python`` or ``ipython``), but also inside Yade, which is a super-set of Python.

Numerical operations and modules:

.. ipython::

	Yade [1]: (1+3*4)**2        # usual rules for operator precedence, ** is exponentiation

	Yade [2]: import math       # gain access to "module" of functions

	Yade [3]: math.sqrt(2)      # use a function from that module

	Yade [4]: import math as m  # use the module under a different name

	Yade [5]: m.cos(m.pi)       

	Yade [6]: from math import *  # import everything so that it can be used without module name

	Yade [7]: cos(pi)

Variables:
	
.. ipython::
	
	Yade [1]: a=1; b,c=2,3       # multiple commands separated with ;, multiple assignment

	Yade [2]: a+b+c


Sequences
^^^^^^^^^

Lists
""""""

Lists are variable-length sequences, which can be modified; they are written with braces ``[...]``, and their elements are accessed with numerical indices:

.. ipython::

	Yade [3]: a=[1,2,3]          # list of numbers

	Yade [4]: a[0]               # first element has index 0

	Yade [5]: a[-1]              # negative counts from the end

	Yade [7]: a[3]               # error

	Yade [5]: len(a)             # number of elements

	Yade [6]: a[1:]              # from second element to the end

	Yade [7]: a+=[4,5]           # extend the list

	Yade [8]: a+=[6]; a.append(7) # extend with single value, both have the same effect

	Yade [9]: 9 in a             # test presence of an element

Lists can be created in various ways:

.. ipython::

	Yade [9]: range(10)

	Yade [1]: range(10)[-1] 

List of squares of even number smaller than 20, i.e. $\left\{a^2\;\forall a\in \{0,\cdots,19\} \;\middle|\; 2 \| a\right\}$ (note the similarity):

.. ipython::

	Yade [1]: [a**2 for a in range(20) if a%2==0]



Tuples
"""""""

Tuples are constant sequences:
	
.. ipython::

	Yade [1]: b=(1,2,3)

	Yade [2]: b[0]

	Yade [3]: b[0]=4              # error

Dictionaries
"""""""""""""

Mapping from keys to values:

.. ipython::

	Yade [1]: czde={'jedna':'ein','dva':'zwei','tri':'drei'}

	Yade [1]: de={1:'ein',2:'zwei',3:'drei'}; cz={1:'jedna',2:'dva',3:'tri'}

	Yade [1]: czde['jedna']         ## access values

	Yade [2]: de[1], cz[2]


Functions, conditionals
^^^^^^^^^^^^^^^^^^^^^^^^

.. ipython::

	Yade [1]: 4==5

	Yade [2]: a=3.1

	Yade [3]: if a<pi: b=0           # conditional statement
	     ...: else: b=1
		  ...:

	Yade [4]: c=0 if a<1 else 1      # conditional expression

	Yade [5] b,c:

	Yade [1]: def square(x): return x**2    # define a new function
	     ...:

	Yade [1]: square(2)                     # and call that function


.. rubric:: Exercises

#. Read the following code and say what wil be the
   values of ``a`` and ``b``::

      a=range(5)
      b=[(aa**2 if aa%2==0 else -aa**2) for aa in a]


Yade basics
------------

Yade objects are constructed in the following manner (this process is also called "instantiation", since we create concrete instances of abstract classes: one individual sphere is an instance of the abstract :yref:`Sphere`, like Socrates is an instance of "man"):

.. ipython::

	Yade [1]: from yade import *

	Yade [1]: Sphere           # try also Sphere?
	
	Yade [1]: s=Sphere()       # create a Sphere, without specifying any attributes

	Yade [1]: s.radius         # 'nan' is a special value meaning "not a number" (i.e. not defined)

	Yade [1]: s.radius=2       # set radius of an existing object

	Yade [1]: s.radius            

	Yade [1]: ss=Sphere(radius=3)   # create Sphere, giving radius directly

	Yade [1]: s.radius, ss.radius     # also try typing s.<tab> to see defined attributes

Particles
^^^^^^^^^

Particles are the "data" component of simulation; they are the objects that will undergo some processes, though do not define those processes yet.

Singles
"""""""""

There is a number of pre-defined functions to create particles of certain type; in order to create a sphere, one has to (see the source of :yref:`yade.utils.sphere` for instance):

#. Create :yref:`Body`
#. Set :yref:`Body.shape` to be an instance of :yref:`Sphere` with some given radius
#. Set :yref:`Body.material` (last-defined material is used, otherwise a default material is created)
#. Set position and orientation in :yref:`Body.state`, compute mass and moment of inertia based on :yref:`Material` and :yref:`Shape`

In order to avoid such tasks, shorthand functions are defined in the :yref:`yade.utils` module; to mention a few of them, they are :yref:`yade.utils.sphere`, :yref:`yade.utils.facet`, :yref:`yade.utils.wall`.

.. ipython::

	Yade [1]: from yade import utils
	
	Yade [1]: s=utils.sphere((0,0,0),radius=1)    # create sphere particle centered at (0,0,0) with radius=1

	Yade [1]: s.shape                       # s.shape describes the geometry of the particle

	Yade [1]: s.shape.radius                # we already know the Sphere class

	Yade [1]: s.state.mass, s.state.inertia # inertia is computed from density and geometry

	Yade [1]: s.state.pos                   # position is the one we prescribed

	Yade [1]: s2=utils.sphere((-2,0,0),radius=1,fixed=True)     # explanation below

In the last example, the particle was fixed in space by the ``fixed=True`` parameter to :yref:`yade.utils.sphere`; such a particle will not move, creating a primitive boundary condition.

A particle object is not yet part of the simulation; in order to do so, a special function is called:

.. ipython::

	Yade [1]: O.bodies.append(s)            # adds particle s to the simulation; returns id of the particle(s) added


Packs
"""""

There are functions to generate a specific arrangement of particles in the :yref:`yade.pack` module; for instance, cloud (random loose packing) of spheres can be generated with the :yref:`yade._packSpheres.SpherePack` class:

.. ipython::

	Yade [1]: from yade import pack

	Yade [1]: sp=pack.SpherePack()                   # create an empty cloud; SpherePack contains only geometrical information

	Yade [1]: sp.makeCloud((1,1,1),(2,2,2),rMean=.2) # put spheres with defined radius inside box given by corners (1,1,1) and (2,2,2)

	Yade [1]: for c,r in sp: print c,r               # print center and radius of all particles (SpherePack is a sequence which can be iterated over)
	     ...:

	Yade [1]: sp.toSimulation()                      # create particles and add them to the simulation

.. note::  :yref:`Vector3<external:miniEigen.Vector3>` is class representing 3D vector; a number of geometry operations are supported by :yref:`Vector3<external:miniEigen.Vector3>`, such as scalar product, dot product, cross product, norm, addition and more. There are other similar classes such as :yref:`Vector2<external:miniEigen.Vector2>`, :yref:`Matrix3<external:miniEigen.Matrix3>` (3×3 matrix), :yref:`Quaterion<external:miniEigen.Quaternion>`, :yref:`Vector6<external:miniEigen.Vector6>` etc. See :yref:`miniEigen<external:miniEigen>` module documentation for details.

Boundaries
""""""""""

:yref:`yade.utils.facet` (triangle :yref:`Facet`) and :yref:`yade.utils.wall` (infinite axes-aligned plane :yref:`Wall`) geometries are typically used to define boundaries. For instance, a "floor" for the simulation can be created like this:
 
.. ipython::

	Yade [1]: O.bodies.append(utils.wall(-1,axis=2))

There are other conveinence functions (like :yref:`yade.utils.facetBox` for creating closed or open rectangular box, or family of :yref:`yade.ymport` functions)

Look inside
^^^^^^^^^^^^

The simulation can be inspected in several ways. All data can be accessed from python directly:

.. ipython::

	Yade [1]: len(O.bodies)

	Yade [1]: O.bodies[1].shape.radius   # radius of body #1 (will give error if not sphere, since only spheres have radius defined)

	Yade [1]: O.bodies[2].state.pos      # position of body #2

Besides that, Yade says this at startup (the line preceding the command-line)::

	[[ ^L clears screen, ^U kills line. F12 controller, F11 3d view, F10 both, F9 generator, F8 plot. ]]

:guilabel:`Controller`
	Pressing ``F12`` brings up a window for controlling the simulation. Although typically no human intervention is done in large simulations (which run "headless", without any graphical interaction), it can be handy in small examples. There are basic information on the simulation (will be used later).
:guilabel:`3d view`
	The 3d view can be opened with F11 (or by clicking on button in the *Controller* -- see below). There is a number of keyboard shortcuts to manipulate it (press ``h`` to get basic help), and it can be moved, rotated and zoomed using mouse.  Display-related settings can be set in the "Display" tab of the controller (such as whether particles are drawn).
:guilabel:`Inspector`
	*Inspector* is opened by clicking on the appropriate button in the *Controller*. It shows (and updated) internal data of the current simulation. In particular, one can have a look at engines, particles (*Bodies*) and interactions (*Interactions*). Clicking at each of the attribute names links to the appropriate section in the documentation.


.. rubric:: Exercises

#. What is this code going to
   do?
   
   .. ipython::

      Yade [1]: O.bodies.append([utils.sphere((2*i,0,0),1) for i in range(1,20)])

#. Create a simple simulation with cloud of spheres enclosed in the box ``(0,0,0)`` and ``(1,1,1)`` with mean radius .1. (hint: :yref:`yade._packSpheres.SpherePack.makeCloud`)

#. Enclose the cloud created above in box with corners ``(0,0,0)`` and ``(1,1,1)``; keep the top of the box open. (hint: :yref:`yade.utils.facetBox`; type ``utils.facetBox?`` or ``utils.facetBox??`` to get help on the command line)

#. Open the 3D view, try zooming in/out; position axes so that $z$ is upwards, $y$ goes to the right and $x$ towards you.

Engines
^^^^^^^

Engines define processes undertaken by particles. As we know from the theoretical introduction, the sequence of engines is called *simulation loop*. Let us define a simple interaction loop:

.. ipython::
	
	Yade [1]: O.reset()

	Yade [1]: O.engines=[                   # newlines and indentations are not important until the brace is closed
	     ...:    ForceResetter(),
	     ...:    InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()]),
	     ...:    InteractionLoop(           # dtto for the parenthesis here
	     ...:        [Ig2_Sphere_Sphere_L3Geom(),Ig2_Wall_Sphere_L3Geom()],
	     ...:        [Ip2_FrictMat_FrictMat_FrictPhys()],
	     ...:        [Law2_L3Geom_FrictPhys_ElPerfPl()]
	     ...:    ),
	     ...:    NewtonIntegrator(damping=.2,label='newton')      # define a name under which we can access this engine easily
	     ...: ]
	     ...:

	Yade [1]: O.engines

	Yade [1]: O.engines[-1]==newton    # is it the same object?

	Yade [1]: newton.damping

Instead of typing everything into the command-line, one can describe simulation in a file (*script*) and then run yade with that file as an argument. We will therefore no longer show the command-line unless necessary; instead, only the script part will be shown. Like this::

	O.engines=[                   # newlines and indentations are not important until the brace is closed
		ForceResetter(),
		InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Wall_Aabb()]),
		InteractionLoop(           # dtto for the parenthesis here
			 [Ig2_Sphere_Sphere_L3Geom_Inc(),Ig2_Wall_Sphere_L3Geom_Inc()],
			 [Ip2_FrictMat_FrictMat_FrictPhys()],
			 [Law2_L3Geom_FrictPhys_ElPerfPl()]
		),
		GravityEngine(gravity=(0,0,-9.81)),              # 9.81 is the gravity acceleration, and we say that
		NewtonIntegrator(damping=.2,label='newton')      # define a name under which we can access this engine easily
	]

Besides engines being run, it is likewise important to define how often they will run. Some engines can run only sometimes (we will see this later), while most of them will run always; the time between two successive runs of engines is *timestep* ($\Dt$). There is a mathematical limit on the timestep value, called *critical timestep*, which is computed from properties of particles. Since there is a function for that, we can just set timestep using :yref:`yade.utils.PWaveTimeStep`::

	O.dt=utils.PWaveTimeStep()

Each time when the simulation loop finishes, time ``O.time`` is advanced by the timestep ``O.dt``:

.. ipython::

	Yade [1]: O.dt=0.01

	Yade [1]: O.time

	Yade [1]: O.step()

	Yade [1]: O.time

For experimenting with a single simulations, it is handy to save it to memory; this can be achieved, once everything is defined, with::

	O.saveTmp()


.. rubric:: Exercises

#. Define *engines* as in the above example, run the *Inspector* and click through the engines to see their sequence.
#. Write a simple script which will

   #. define particles as in the previous exercise (cloud of spheres inside a box open from the top)
   #. define a simple simulation loop, as the one given above
   #. set $\Dt$ equal to the critical P-Wave $\Dt$
   #. save the initial simulation state to memory

#. Run the previously-defined simulation multiple times, while changing the value of timestep (use the :guilabel:`⟳` button to reload the initial configuration).

   #. See what happens as you increase $\Dt$ above the P-Wave value.
   #. Try changing the :yref:`gravity<GravityEngine.gravity>` parameter, before running the simulation.
   #. Try changing :yref:`damping<NewtonIntegrator.damping>`

#. Reload the simulation, open the 3d view, open the *Inspector*, select a particle in the 3d view (shift-click). Then run the simulation and watch how forces on that particle change; pause the simulation somewhere in the middle, look at interactions of this particle.

#. At which point can we say that the deposition is done, so that the simulation can be stopped?

.. seealso::

	The :ref:`bouncing-sphere` example shows a basic simulation.
