**********************
Programmer's manual
**********************

Build system
=============

Yade uses `scons <http://www.scons.org>`__ build system for managing the build process. It takes care of configuration, compilation and installation. SCons is written in python and its build scripts are in python, too. SCons complete documentation can be found in its manual page.


.. _scons-parameters:
Pre-build configuration
-----------------------

We use ``\$`` to denote build variable in strings in this section; in SCons script, they can be used either by writing ``\$variable`` in strings passed to SCons functions, or obtained as attribute of the ``Environment`` instance ``env``, i.e. ``env['variable']``; we use the formed in running text here.

In order to allow parallel installation of multiple yade versions, the installation location follows the pattern ``\$PREFIX/lib/yade\$SUFFIX`` for libraries and ``\$PREFIX/bin/yade\$SUFFIX`` for executables (in the following, we will refer only to the first one). ``\$SUFFIX`` takes the form ``-\$version\$variant``, which further allows multiple different builds of the same version (typically, optimized and debug builds). For instance, the default debug build of version 0.5 would be  installed in ``/usr/local/lib/yade-0.5-dbg/``, the executable being ``/usr/local/bin/yade-0.5-dbg``.

The build process takes place outside the source tree, in directory reffered to as ``\$buildDir`` within those scripts. By default, this directory is ``../build-\$SUFFIX``.

Each build depends on a number of configuration parameters, which are stored in mutually independent *profiles*. They are selected according to the ``profile`` argument to scons (by default, the last profile used, stored in ``scons.current-profile``). Each profile remembers its non-default variables in ``scons.profile-\$profile``.

There is a number of configuration parameters; you can list all of them by ``scons -h``. The following table summarizes only a few that are the most used.

``PREFIX`` [default: ``/usr/local``]
	installation prefix (``PREFIX`` preprocessor macro; ``yade.config.prefix`` in python
``version`` [bzr revision (e.g. bzr1899)]
	first part of suffix (``SUFFIX`` preprocessor macro; ``yade.config.suffix`` in python)]
``variant`` [*(empty)*]
	second part of suffix
``buildPrefix`` [``..``]
	 where to create ``build-\$SUFFIX`` directory 
``debug`` [*False* (0)]
	add debugging symbols to output, enable stack traces on crash
``optimize`` [-1, which means that the opposite of ``debug`` value is used]
	optimize binaries (``#define NDEBUG``; assertions eliminated; ``YADE_CAST`` and ``YADE_PTR_CAST`` are static casts rather than dynamic; LOG_TRACE and LOG_DEBUG are eliminated)
``CPPPATH`` [``/usr/include/vtk-5.2:/usr/include/vtk-5.4``]
	additional colon-separated paths for preprocessor (for atypical header locations). Required by some libraries, such as VTK (reflected by the default)
``LIBPATH`` [*(empty)*]
	additional colon-separated paths for linker
``CXX`` [g++]
	compiler executable
``CXXFLAGS`` [*(empty)*]
	additional compiler flags (may are added automatically)
``jobs`` [4]
	number of concurrent compilations to run 
``brief`` [*True* (1)]
	only show brief notices about what is being done rather than full command-lines during compilation
``linkStrategy`` [monolithic]
	whether to link all plugins in one shared library (``monolithic``) or in one file per plugin (``per-class``); the first option is faster for overall builds, while the latter one makes recompilation of only part of Yade faster; granularity of monolithic build can be changed with the ``chunkSize`` parameter, which determines how many files are compiled at once.
``features`` [log4cxx,opengl,gts,openmp]
	optional comma-separated features to build with (details below; each defines macro ``YADE_\$FEATURE``; available as lowercased list ``yade.config.features`` at runtime


Library detection
^^^^^^^^^^^^^^^^^^
When the ``scons`` command is run, it first checks for presence of all required libraries. Some of them are *essential*, other are *optional* and will be required only if features that need them are enabled.

Essentials
"""""""""""""

compiler
	Obviously c++ compiler is necessary. Yade relies on several extensions of ``g++`` from the `gcc <http://gcc.gnu.org`__ suite and cannot (probably) be built with other compilers.
boost
	`boost <http://www.boost.org>`__ is a large collection of peer-reviewed c++ libraries. Yade currently uses thread, date_time, filesystem, iostreams, regex, serialization, program_options, foreach, python; typically the whole boost bundle will be installed. If you need functionality from other modules, you can make presence of that module mandatory. Only be careful about relying on very new features; due to range of systems yade is or might be used on, it is better to be moderately conservative (read: roughly 3 years backwards compatibility).
python
	`python <http://www.python.org>`__ is the scripting language used by yade. Besides [boost::python]_, yade further requires

	* `ipython <http://www.ipython.org>`__ (terminal interaction)
	* `matplotlib <http://matplotlib.sf.net>`__ (plotting)
	* `numpy <http://www.numpy.org>`__ (matlab-like numerical functionality and accessing numpy arrays from ``c``/``c++`` efficiently)

.. _optional-libraries:
Optional libraries (features)
""""""""""""""""""""""""""""""

The `features` parameter controls optional functionality. Each enabled feature defines preprocessor macro `YADE_FEATURE` (name uppercased) to enable selective exclude/include of parts of code. Code of which compilation depends on a particular features should use ``#ifdef YADE_FEATURE`` constructs to exclude dependent parts.

log4cxx (YADE_LOG4CXX)
	Enable flexible logging system (`log4cxx <http://log4cxx.apache.org>`__), which permits to assign logging levels on per-class basis; doesn't change API, only redefines `LOG_INFO` and other macros accordingly; see :ref:`log4cxx` for details.
opengl (YADE_OPENGL)
	Enable 3d rendering as well as the Qt3-based graphical user interface (in addition to python console).
vtk (YADE_VTK)
	Enable functionality using Visualization Toolkit (`vtk <http://www.vtk.org>`__; e.g. :yref:`VTKRecorder` exporting to files readable with ParaView).
openmp (YADE_OPENMP)
	Enable parallelization using OpenMP, non-intrusive shared-memory parallelization framework; it is only supported for ``g++`` > 4.0. Parallel computation leads to significant performance increase and should be enabled unless you have a special reason for not doing so (e.g. single-core machine). See :ref:`upyade-parallel` for details.
gts (YADE_GTS)
	Enable functionality provided by GNU Triangulated Surface library (`gts <http://gts.sf.net>`__) and build PyGTS, its python interface; used for surface import and construction.
cgal (YADE_CGAL)
	Enable functionality provided by Computation Geometry Algorithms Library (`cgal <http://www.cgal.org>`__); triangulation code in :yref:`MicroMacroAnalyser` and :yref:`PersistentTriangulationCollider` ses its routines.
other
	There might be more features added in the future. Always refer to ``scons -h`` output for possible values.

.. warning:: 
	Due to a long-standing `bug <http://issues.apache.org/jira/browse/LOGCXX-322>`_ in log4cxx, using log4cxx will make yade crash at every exit. We work-around this partially by disabling the crash handler for regular exits, but process exit status will still be non-zero. The batch system (:ref:`yade-multi`) detects successful runs by looking at magic line "Yade: normal exit." in the process' standard output.


Before compilation, SCons will check for presence of libraries required by their respective features [#features]_. Failure will occur if a respective library isn't found. To find out what went wrong, you can inspect ``../build-\$SUFFIX/config.log`` file; it contains exact commands and their output for all performed checks.

.. [#features] Library checks are defined inside the ``SConstruct`` file and you can add your own, should you need it.

.. note::
	Features are not auto-detected on purpose; otherwise problem with library detection might build Yade without expected features, causing specifically problems for automatized builds.

Building
-------------

Yade source tree has the following structure (omiting ``debian``, ``doc``, ``examples`` and ``scripts`` which don't participate in the build process); we shall call each top-level component *module*::

	attic/        ## code that is not currently functional and might be removed unless resurrected
	   lattice/      ## lattice and lattice-like models
	   snow/         ## snow model (is really a DEM)
	core/         ## core simulation building blocks
	extra/        ## miscillanea
	gui/          ## user interfaces
	   qt3/          ## graphical user interface based on qt3 and OpenGL
	   py/           ## python console interface (phased out)
	lib/          ## support libraries, not specific to simulations
	pkg/          ## simulation-specific files
	   common/       ## generally useful classes
	   dem/          ## classes for Discrete Element Method
	py/           ## python modules

Each directory on the top of this hierarchy (except ``pkg``, which is treated specially -- see below) contains file ``SConscript``, determining what files to compile, how to link them together and where should they be installed. Within these script, a scons variable ``env`` (build ``Environment``) contains all the configuration parameters, which are used to influence the build process; they can be either obtained with the ``[]`` operator, but scons also replaces ``\$var`` strings automatically in arguments to its functions::

	if 'opengl' in env['features']:
		env.Install('\$PREFIX/lib/yade\$SUFFIX/',[
			# ...
		])


Header installation
^^^^^^^^^^^^^^^^^^^^
To allow flexibility in source layout, SCons will copy (symlink) all headers into flattened structure within the build directory. First 2 components of the original directory are joind by dash, deeper levels are discarded (in case of ``core`` and ``extra``, only 1 level is used). The following table makes gives a few examples:

============================================================= =========================
Original header location											     Included as     
============================================================= =========================
``core/Scene.hpp``														  ``<yade/core/Scene.hpp>``
``lib/base/Logging.hpp``												  ``<yade/lib-base/Logging.hpp>``
``lib/serialization/Serializable.hpp``								  ``<yade/lib-serialization/Serializable.hpp>``
``pkg/dem/DataClass/SpherePack.hpp``                          ``<yade/pkg-dem/SpherePack.hpp>``
``gui/qt3/QtGUI.hpp``                                         ``<yade/gui-qt3/QtGUI.hpp>``
============================================================= =========================

It is advised to use ``#include<yade/module/Class.hpp>`` style of inclusion rather than ``#include"Class.hpp`` even if you are in the same directory.

What files to compile
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``SConscript`` files in ``lib``, ``core``, ``gui``, ``py`` and ``extra`` explicitly determine what files will be built.

Automatic compilation
""""""""""""""""""""""
In the ``pkg/`` directory, situation is different. In order to maximally ease addition of modules to yade, all ``*.cpp`` files are *automatically scanned* by SCons and considered for compilation. Each file may contain multiple lines that declare features that are necessary for this file to be compiled::

	YADE_REQUIRE_FEATURE(vtk);
	YADE_REQUIRE_FEATURE(gts);

This file will be compiled only if *both* ``vtk`` and ``gts`` features are enabled. Depending on current feature set, only selection of plugins will be compiled.

It is possible to disable compilation of a file by requiring any non-existent feature, such as::

	YADE_REQUIRE_FEATURE(temporarily disabled 345uiysdijkn);

The ``YADE_REQUIRE_FEATURE`` macro expands to nothing during actual compilation.

.. note::
	The source scanner was written by hand and is not official part of SCons. It is fairly primitive and in particular, it doesn't interpret c preprocessor macros, except for a simple non-nested feature-checks like ``#ifdef YADE_*``/``#ifndef YADE_*`` #endif.


.. _linking:

Linking
^^^^^^^^
The order in which modules might depend on each other is given as follows:

=========== ================================== ==============================================
module      resulting shared library           dependencies
=========== ================================== ==============================================
lib         ``libyade-support.so``             can depend on external libraries, may **not** depend on any other part of Yade.
core        ``libcore.so``                     ``yade-support``; *may* depend on external libraries.
pkg         ``libplugins.so`` for monolithic   ``core``, ``yade-support``; may **not** depend on external libraries explcitly (only implicitly, by adding the library to global linker flags in ``SConstruct``)
            builds, ``libClass.so`` for
            per-class (per-plugin) builds.
extra 		(undefined)				    			  (arbitrary)
gui         ``libQtGUI.so``,                   ``lib``, ``core``, ``pkg``
            ``libPythonUI.so``
py          (many files)                       ``lib``, ``core``, ``pkg``, external
=========== ================================== ==============================================

Because ``pkg`` plugins might be linked differently depending on the ``linkStrategy`` option, ``SConscript`` files that need to explicitly declare the dependency should use provided ``linkPlugins`` function which returns libraries in which given plugins will be defined::

	env.SharedLibrary('_packSpheres',['_packSpheres.cpp'],
		SHLIBPREFIX='',
		LIBS=env['LIBS']+[linkPlugins(['Shop','SpherePack']),]
	),

.. note:: ``env['LIBS']`` are libraries that all files are linked to and they should always be part of the ``LIBS`` parameter.

Since plugins in ``pkg`` are not declared in any ``SConscript`` file, other plugins they depend on are again found *automatically* by scannig their ``#include`` directives for the pattern ``#include<yade/module/Plugin.hpp>``. Again, this works well in normal circumastances, but is not necessarily robust.

See scons manpage for meaning of parameters passed to build functions, such as ``SHLIBPREFIX``.

Development tools
=================

Integrated Development Environment and other tools
---------------------------------------------------
A frequently used IDE is Kdevelop. We recommend using this software for navigating in the sources, compiling and debugging. Other usefull tools for debugging and profiling are Valgrind and KCachegrind. A series a wiki pages are dedicated to these tools in the `development section <https://yade-dem.org/wiki/Yade#Development>`__ of the wiki.

Hosting and versioning
----------------------
The Yade project is kindly hosted at `launchpad <https://launchpad.net/yade/>`__, which is used for source code, bug tracking, planning, package downloads and more. Our repository `can be http-browsed <http://bazaar.launchpad.net/~vcs-imports/yade/trunk/files>`__.

The versioning software used is `Bazaar <http://www.bazaar-vcs.org>`__, for which a short tutorial can be found in a `Yade's wiki pages <https://yade-dem.org/wiki/Quick_Bazaar_tutorial>`__. Bazaar is a distributed revision control system. It is available packaged for all major linux distributions.

Build robot
-----------
A build robot hosted at `3SR lab. <http://www.3s-r.hmg.inpg.fr/3sr/?lang=en>`__ is tracking souce code changes.
Each time a change in the source code is commited to the main development branch via bazaar, the "buildbot" downloads and compiles the new version, and start a series of tests.

If a compilation error has been introduced, it will be notified to the yade-dev mailing list and to the commiter, thus helping to fix problems quickly.
If the compilation is successfull, the buildbot starts unit regression tests and "check tests" (see below) and report the results. If all tests are passed, a new version of the documentation is generated and uploaded to the website in `html <https://www.yade-dem.org/doc/>`__ and `pdf <https://yade-dem.org/doc/Yade.pdf>`__ format. As a consequence, those two links always point to the documentation (the one you are reading now) of the last successfull build, and the delay between commits and documentation updates are very short (minutes).
The buildbot activity and logs can be `browsed online <https://yade-dem.org/buildbot/>`__.

Regression tests
----------------
Yade contains two types of regression tests, some are unit tests while others are testing more complex simulations. Altough both types can be considered regression tests, the usage is that we name the first simply "regression tests", while the latest are called "check tests".
Both series of tests can be ran at yade startup by passing the options "test" or "check" ::

	yade --test
	yade --check

Unit regression tests
^^^^^^^^^^^^^^^^^^^^^
Unit regression tests are testing the output of individual functors and engines in well defined conditions. They are defined in the folder :ysrc:`py/tests/`.
The purpose of unit testing is to make sure that the behaviour of the most important classes remains correct during code development. Since they test classes one by one, unit tests can't detect problems coming from the interaction between different engines in a typical simulation. That is why check tests have been introduced. 

Check tests
^^^^^^^^^^^
Check tests perform comparisons of simulation results between different versions of yade, as discussed in http://www.mail-archive.com/yade-dev@lists.launchpad.net/msg05784.html and the whole thread. They differ with regression tests in the sense that they simulate more complex situations and combinations of different engines, and usually don't have a mathematical proof (though there is no restriction on the latest). They compare the values obtained in version N with values obtained in a previous version or any other "expected" results. The reference values must be hardcoded in the script itself or in data files provided with the script. Check tests are based on regular yade scripts, so that users can easily commit their own scripts to trunk in order to get some automatized testing after commits from other developers.

Since the check tests history will be mostly based on standard output generated by "yade ---check", a meaningfull checkTest should include some "print" command telling if something went wrong. If the script itself fails for some reason and can't generate an output, the log will contain "scriptName failure". If the script defines differences on obtained and awaited data, it should print some useful information about the problem and increase the value of global variable resultStatus. After this occurs, the automatic test will stop the execution with error message.

An example check test can be found in checkTestTriax.py. It shows results comparison, output, and how to define the path to data files using "checksPath".
Users are encouraged to add their own scripts into the scripts/test/checks/ folder. Discussion of some specific checktests design in users question is welcome. Note that re-compiling is required before that added scripts can be launched by "yade ---check" (or direct changes have to be performed in "lib" subfolders).
A check test should never need more than a few seconds to run. If your typical script needs more, try and reduce the number of element or the number of steps.

Conventions
============

The following rules that should be respected; documentation is treated separately.

* general

  * C++ source files have ``.hpp`` and ``.cpp`` extensions (for headers and implementation, respectively).
  * All header files should have the ``#pragma once`` multiple-inclusion guard.
  * Try to avoid ``using namespace …`` in header files.
  * Use tabs for indentation. While this is merely visual in ``c++``, it has semantic meaning in python; inadverently mixing tabs and spaces can result in syntax errors.

* capitalization style

  * Types should be always capitalized. Use CamelCase for composed names (``GlobalEngine``). Underscores should be used only in special cases, such as functor names.
  * Class data members and methods must not be capitalized, composed names should use use lowercased camelCase (``glutSlices``). The same applies for functions in python modules.
  * Preprocessor macros are uppercase, separated by underscores; those that are used outside the core take (with exceptions) the form ``YADE_*``, such as :ref:`YADE_CLASS_BASE_DOC`.

* programming style

  * Be defensive, if it has no significant performance impact. Use assertions abundantly: they don't affect performance (in the optimized build) and make spotting error conditions much easier.
  * Use logging abundantly. Again, ``LOG_TRACE`` and ``LOG_DEBUG`` are eliminated from optimized code; unless turned on explicitly, the ouput will be suppressed even in the debug build (see below).
  * Use ``YADE_CAST`` and ``YADE_PTR_CAST`` where you want type-check during debug builds, but fast casting in optimized build. 
  * Initialize all class variables in the default constructor. This avoids bugs that may manifest randomly and are difficult to fix. Initializing with NaN's will help you find otherwise unitialized variable. (This is taken care of by :ref:`YADE_CLASS_BASE_DOC` macros for user classes)


Class naming
-------------

Although for historial reasons the naming scheme is not completely consistent, these rules should be obeyed especially when adding a new class.

GlobalEngines and PartialEngines
	GlobalEngines should be named in a way suggesting that it is a performer of certain action (like :yref:`ForceResetter`, :yref:`InsertionSortCollider`, :yref:`Recorder`); if this is not appropriate, append the ``Engine`` to the characteristics (:yref:`GravityEngine`). :yref:`PartialEngines<PartialEngine>` have no special naming convention different from :yref:`GlobalEngines<GlobalEngine>`.

Dispatchers
	Names of all dispatchers end in ``Dispatcher``. The name is composed of type it creates or, in case it doesn't create any objects, its main characteristics. Currently, the following dispatchers [#opengldispatchers]_ are defined:

.. _dispatcher-names:

	========================== ================ ====================== ============== ===================== ===============
	dispatcher                 arity            dispatch types         created type   functor type          functor prefix
	========================== ================ ====================== ============== ===================== ===============
	:yref:`BoundDispatcher`    1                :yref:`Shape`          :yref:`Bound`  :yref:`BoundFunctor`  ``Bo1``
	:yref:`IGeomDispatcher`    2 (symetric)     2 × :yref:`Shape`      :yref:`IGeom`  :yref:`IGeomFunctor`  ``Ig2``
	:yref:`IPhysDispatcher`    2 (symetric)     2 × :yref:`Material`   :yref:`IPhys`  :yref:`IPhysFunctor`  ``Ip2``
	:yref:`LawDispatcher`      2 (asymetric)    :yref:`IGeom`          *(none)*       :yref:`LawFunctor`    ``Law2``
	                                            :yref:`IPhys`
	========================== ================ ====================== ============== ===================== ===============

	Respective abstract functors for each dispatchers are :yref:`BoundFunctor`, :yref:`IGeomFunctor`, :yref:`IPhysFunctor` and :yref:`LawFunctor`.

Functors
	Functor name is composed of 3 parts, separated by underscore.

	#. prefix, composed of abbreviated functor type and arity (see table above)
	#. Types entering the dispatcher logic (1 for unary and 2 for binary functors)
	#. Return type for functors that create instances, simple characteristics for functors that don't create instances.

	To give a few examples:

	* :yref:`Bo1_Sphere_Aabb` is a :yref:`BoundFunctor` which is called for :yref:`Sphere`, creating an instance of :yref:`Aabb`.
	* :yref:`Ig2_Facet_Sphere_Dem3DofGeom` is binary functor called for :yref:`Facet` and :yref:`Sphere`, creating and instace of :yref:`Dem3DofGeom`.
	* :yref:`Law2_Dem3DofGeom_CpmPhys_Cpm` is binary functor (:yref:`LawFunctor`) called for types :yref:`Dem3Dof (Geom)<Dem3DofGeom>` and :yref:`CpmPhys`.

.. [#opengldispatchers] Not considering OpenGL dispatchers, which might be replaced by regular virtual functions in the future.

Documentation
---------------

**Documenting code properly is one of the most important aspects of sustained development.**

Read it again.

Most code in research software like Yade is not only used, but also read, by developers or even by regular users. Therefore, when adding new class, always mention the following in the documentation:

* purpose
* details of the functionality, unless obvious (algorithms, internal logic)
* limitations (by design, by implementation), bugs
* bibliographical reference, if using non-trivial published algorithms (see below)
* references to other related classes
* hyperlinks to bugs, blueprints, wiki or mailing list about this particular feature.

As much as it is meaningful, you should also

* update any other documentation affected
* provide a simple python script demonstrating the new functionality in ``scripts/test``.


Historically, Yade was using Doxygen for in-source documentation. This documentation is still available (by running ``scons doc``), but was rarely written and used by programmers, and had all the disadvantages of auto-generated documentation. Then, as Python became ubiquitous in yade, python was documented using epydoc generator. Finally, hand-written documentation (this one) started to be written using `Sphinx <http://sphinx.pocoo.org>`__, which was developed originally for documenting Python itself. Disadvantages of the original scatter were different syntaxes, impossibility for cross-linking, non-interactivity and frequently not being up-to-date.

.. _sphinxdocumentation:

Sphinx documentation
^^^^^^^^^^^^^^^^^^^^^
Most c++ classes are wrapped in Python, which provides good introspection and interactive documentation (try writing ``Material?`` in the ipython prompt; or ``help(CpmState)``).

Syntax of documentation is `ReST <http://docutils.sourceforge.net/rst.html`__ (reStructuredText, see `reStructuredText Primer <http://sphinx.pocoo.org/rest.html>`__). It is the same for c++ and python code.

* Documentation of c++ classes exposed to python is given as 3rd argument to :ref:`YADE_CLASS_BASE_DOC` introduced below.

* Python classes/functions are documented using regular python docstrings. Besides explaining functionality, meaning and types of all arguments should also be documented. Short pieces of code might be very helpful. See the :yref:`yade.utils` module for an example.


In addition to standard ReST syntax, yade provides several shorthand macros:

``:yref:``
	creates hyperlink to referenced term, for instance::
		
		:yref:`CpmMat`
	
	becomes :yref:`CpmMat`; link name and target can be different::

		:yref:`Material used in the CPM model<CpmMat>`

	yielding :yref:`Material used in the CPM model<CpmMat>`.

``:ysrc:``
	creates hyperlink to file within the source tree (to its latest version in the repository), for instance :ysrc:`core/Cell.hpp`. Just like with `:yref:`, alternate text can be used with ::
	
		:ysrc:`Link text<target/file>`
		
	like :ysrc:`this<core/Cell.hpp>`.

``|ycomp|``
	is used in attribute description for those that should not be provided by the used, but are auto-computed instead; ``|ycomp|`` expands to |ycomp|.

``|yupdate|``
	marks attributes that are periodically update, being subset of the previous. ``|yupdate|`` expands to |yupdate|.
	
``\$...\$``
	delimits inline math expressions; they will be replaced by::
		
		:math:`...`

	and rendered via LaTeX. To write a single dollar sign, escape it with backslash ``\\$``.
	
	Displayed mathematics (standalone equations) can be inserted as explained in `Math support in Sphinx <http://sphinx.pocoo.org/ext/math.html>`_.



Bibliographical references
^^^^^^^^^^^^^^^^^^^^^^^^^^^

As in any scientific documentation, references to publications are very important. To cite an article, add it to BibTeX file in :ysrc:`doc/references.bib`, using the BibTeX format. Please adhere to the following conventions:

#. Keep entries in the form ``Author2008`` (``Author`` is the first author), ``Author2008b`` etc if multiple articles from one author;
#. Try to fill `mandatory fields <http://en.wikipedia.org/wiki/Bibtex#Entry_Types>`_ for given type of citation;
#. Do not use ``\'{i}`` funny escapes for accents, since they will not work with the HTML output; put everything in straight utf-8.

In your docstring, the ``Author2008`` article can be cited by ``[Author2008]_``; for example::

	According to [Allen1989]_, the integration scheme …

will be rendered as

	According to [Allen1989]_, the integration scheme …

Separate class/function documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some c++ might have long or content-rich documentation, which is rather inconvenient to type in the c++ source itself as string literals. Yade provides a way to write documentation separately in :ysrc:`py/_extraDocs.py` file: it is executed after loading c++ plugins and can set ``__doc__`` attribute of any object directly, overwriting docstring from c++. In such (exceptional) cases:

#. Provide at least a brief description of the class in the c++ code nevertheless, for people only reading the code.
#. Add notice saying "This class is documented in detail in the :ysrc:`py/_extraDocs.py` file". 
#. Add documentation to :ysrc:`py/_etraDocs.py` in this way::
		
		module.YourClass.__doc__='''
			This is the docstring for YourClass.

			Class, methods and functions can be documented this way.

			.. note:: It can use any syntax features you like.

		'''

.. note:: Boost::python embeds function signatures in the docstring (before the one provided by the user). Therefore,  before creating separate documentation of your function, have a look at its ``__doc__`` attribute and copy the first line (and the blank lie afterwards) in the separate docstring. The first line is then used to create the function signature (arguments and return value).

Local documentation
^^^^^^^^^^^^^^^^^^^

.. note:: At some future point, this documentation will be integrated into yade's sources. This section should be updated accordingly in that case.

To generate Yade's documentation locally, get a copy of the `ydoc branch <https://launchpad.net/~eudoxos/+junk/ydoc/>`_ via bzr, then follow instructions in the `README <http://bazaar.launchpad.net/~eudoxos/%2Bjunk/ydoc/annotate/head%3A/README>`_ file.



Internal c++ documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

`doxygen <http://www.doxygen.org>`__ was used for automatic generation of c++ code. Since user-visible classes are defined with sphinx now, it is not meaningful to use doxygen to generate overall documentation. However, take care to document well internal parts of code using regular comments, including public and private data members.


Support framework
=======================

Besides the framework provided by the c++ standard library (including STL), boost and other dependencies, yade provides its own specific services.


Pointers
---------

Shared pointers
^^^^^^^^^^^^^^^^^
Yade makes extensive use of shared pointers ``shared_ptr``. [#sharedptr]_ Although it probably has some performance impacts, it greatly simplifies memory management, ownership management of c++ objects in python and so forth. To obtain raw pointer from a ``shared_ptr``, use its ``get()`` method; raw pointers should be used in case the object will be used only for short time (during a function call, for instance) and not stored anywhere.
 
.. [#sharedptr] Either ``boost::shared_ptr`` or ``tr1::shared_ptr`` is used, but it is always imported with the ``using`` statement so that unqualified ``shared_ptr`` can be used.

Python defines thin wrappers for most c++ Yade classes (for all those registered with :ref:`YADE_CLASS_BASE_DOC` and several others), which can be constructed from ``shared_ptr``; in this way, Python reference counting blends with the ``shared_ptr`` reference counting model, preventing crashes due to python objects pointing to c++ objects that were destructed in the meantime.

Typecasting
^^^^^^^^^^^^

Frequently, pointers have to be typecast; there is choice between static and dynamic casting.

* ``dynamic_cast`` (``dynamic_pointer_cast`` for a ``shared_ptr``) assures cast admissibility by checking runtime type of its argument and returns NULL if the cast is invalid; such check obviously costs time. Invalid cast is easily caught by checking whether the pointer is NULL or not; even if such check (e.g. ``assert``) is absent, dereferencing NULL pointer is easily spotted from the stacktrace (debugger output) after crash. Moreover, ``shared_ptr`` checks that the pointer is non-NULL before dereferencing in debug build and aborts with "Assertion \`px!=0' failed." if the check fails.
 
* ``static_cast`` is fast but potentially dangerous (``static_pointer_cast`` for ``shared_ptr``). Static cast will return non-NULL pointer even if types don't allow the cast (such as casting from ``State*`` to ``Material*``); the consequence of such cast is interpreting garbage data as instance of the class cast to, leading very likely to invalid memory access (segmentation fault, "crash" for short).

To have both speed and safety, Yade provides 2 macros:

``YADE_CAST``
	expands to ``static_cast`` in optimized builds and to ``dynamic_cast`` in debug builds. 

``YADE_PTR_CAST``
	expands to ``static_pointer_cast`` in optimized builds and to ``dynamic_pointer_cast`` in debug builds.



Basic numerics
---------------
The floating point type to use in Yade ``Real``, which is by default typedef for ``double``. [#real]_ 

Yade uses the `Eigen <http://eigen.tuxfamily.org>`_ library for computations. It provides classes for 2d and 3d vectors, quaternions and 3x3 matrices templated by number type; their specialization for the ``Real`` type are typedef'ed with the "r" suffix, and occasionally useful integer types with the "i" suffix:

	* ``Vector2r``, ``Vector2i``
	* ``Vector3r``, ``Vector3i``
	* ``Quaternionr``
	* ``Matrix3r``

Yade additionally defines a class named :yref:`Se3r`, which contains spatial position (``Vector3r Se3r::position``) and orientation (``Quaternionr Se3r::orientation``), since they are frequently used one with another, and it is convenient to pass them as single parameter to functions.

.. [#real] Historically, it was thought that Yade could be also run with single precision based on build-time parameter; it turned out however that the impact on numerical stability was such disastrous that this option is not available now. There is, however, ``QUAD_PRECISION`` parameter to scons, which will make ``Real`` a typedef for ``long double`` (extended precision; quad precision in the proper sense on IA64 processors); this option is experimental and is unlikely to be used in near future, though.

Eigen provides full rich linear algebra functionality. Some code firther uses the [cgal]_ library for computational geometry.

In Python, basic numeric types are wrapped and imported from the ``miniEigen`` module; the types drop the ``r`` type qualifier at the end, the syntax is otherwise similar. ``Se3r`` is not wrapped at all, only converted automatically, rarely as it is needed, from/to a ``(Vector3,Quaternion)`` tuple/list.

.. ipython::

	Yade [0]: from math import pi

	# cross product

	Yade [1]: Vector3(1,2,3).cross(Vector3(0,0,1))  

	# construct quaternion from axis and angle

	Yade [2]: Quaternion(Vector3(0,0,1),pi/2)             

.. note:: Quaternions are internally stored as 4 numbers. Their usual human-readable representation is, however, (normalized) axis and angle of rotation around that axis, and it is also how they are input/output in Python. Raw internal values can be accessed using the ``[0]`` … ``[3]`` element access (or ``.W()``, ``.X()``, ``.Y()`` and ``.Z()`` methods), in both c++ and Python.

.. _rtti:

Run-time type identification (RTTI)
-----------------------------------

Since serialization and dispatchers need extended type and inheritance information, which is not sufficiently provided by standard RTTI. Each yade class is therefore derived from ``Factorable`` and it must use macro to override its virtual functions providing this extended RTTI:

``YADE_CLASS_BASE_DOC(Foo,Bar Baz,"Docstring)`` creates the following virtual methods (mediated via the ``REGISTER_CLASS_AND_BASE`` macro, which is not user-visible and should not be used directly):

* ``std::string getClassName()`` returning class name (``Foo``) as string. (There is the ``typeid(instanceOrType).name()`` standard c++ construct, but the name returned is compiler-dependent.)
* ``unsigned getBaseClassNumber()`` returning number of base classes (in this case, 2).
* ``std::string getBaseClassName(unsigned i=0)`` returning name of *i*-th base class (here, ``Bar`` for i=0 and ``Baz`` for i=1).

.. warning:: RTTI relies on virtual functions; in order for virtual functions to work, at least one virtual method must be present in the implementation (``.cpp``) file. Otherwise, virtual method table (vtable) will not be generated for this class by the compiler, preventing virtual methods from functioning properly.

Some RTTI information can be accessed from python:

.. ipython::
	
	Yade [1]: import yade.system

	Yade [2]: yade.system.childClasses('Shape')

	Yade [3]: Sphere().name            ## getClassName()


Serialization
--------------

Serialization serves to save simulation to file and restore it later. This process has several necessary conditions:

* classes know which attributes (data members) they have and what are their names (as strings);
* creating class instances based solely on its name;
* knowing what classes are defined inside a particular shared library (plugin).

This functionality is provided by 3 macros and 4 optional methods; details are provided below.

``Serializable::preLoad``, ``Serializable::preSave``, ``Serializable::postLoad``, ``Serializable::postSave``
	Prepare attributes before serialization (saving) or deserialization (loading) or process them after serialization or deserialization.
	
	See :ref:`attributeregistration`.
``YADE_CLASS_BASE_DOC_*``
	Inside the class declaration (i.e. in the ``.hpp`` file within the ``class Foo { /* … */};`` block). See :ref:`attributeregistration`.

	Enumerate class attributes that should be saved and loaded; associate each attribute with its literal name, which can be used to retrieve it. See :ref:`YADE_CLASS_BASE_DOC`.                    

	Additionally documents the class in python, adds methods for attribute access from python, and documents each attribute.
``REGISTER_SERIALIZABLE``
	In header file, but *after* the class declaration block. See :ref:`classfactory`.
	
	Associate literal name of the class with functions that will create its new instance (``ClassFactory``).
``YADE_PLUGIN``
	In the implementation ``.cpp`` file. See :ref:`plugins`.

	Declare what classes are declared inside a particular plugin at time the plugin is being loaded (yade startup).

.. _attributeregistration:

Attribute registration
^^^^^^^^^^^^^^^^^^^^^^

All (serializable) types in Yade are one of the following:

* Type deriving from :yref:`Serializable`, which provide information on how to serialize themselves via overriding the ``Serializable::registerAttributes`` method; it declares data members that should be serialzed along with their literal names, by which they are identified. This method then invokes ``registerAttributes`` of its base class (until ``Serializable`` itself is reached); in this way, derived classes properly serialize data of their base classes.

  This funcionality is hidden behind the macro :ref:`YADE_CLASS_BASE_DOC` used in class declaration body (header file), which takes base class and list of attributes::

	YADE_CLASS_BASE_DOC_ATTRS(ThisClass,BaseClass,"class documentation",((type1,attribute1,initValue1,,"Documentation for attribute 1"))((type2,attribute2,initValue2,,"Documentation for attribute 2"));

  Note that attributes are encodes in double parentheses, not separated by commas. Empty attribute list can be given simply by ``YADE_CLASS_BASE_DOC_ATTRS(ThisClass,BaseClass,"documentation",)`` (the last comma is mandatory), or by omiting ``ATTRS`` from macro name and last parameter altogether.

* Fundamental type: strings, various number types, booleans, ``Vector3r`` and others. Their "handlers" (serializers and deserializers) are defined in ``lib/serialization``.

* Standard container of any serializable objects.

* Shared pointer to serializable object.

Yade uses the excellent `boost::serialization <http://www.boost.org/doc/libs/release/libs/serialization/>`_ library internally for serialization of data.

.. note:: ``YADE_CLASS_BASE_DOC_ATTRS`` also generates code for attribute access from python; this will be discussed later. Since this macro serves both purposes, the consequence is that attributes that are serialized can always be accessed from python.

Yade also provides callback for before/after (de) serialization, virtual functions :yref:`Serializable::preProcessAttributes` and :yref:`Serializable::postProcessAttributes`, which receive one ``bool deserializing`` argument (``true`` when deserializing, ``false`` when serializing). Their default implementation in :yref:`Serializable` doesn't do anything, but their typical use is:

* converting some non-serializable internal data structure of the class (such as multi-dimensional array, hash table, array of pointers) into a serializable one (pre-processing) and fill this non-serializable structure back after deserialization (post-processing); for instance, InteractionContainer uses these hooks to ask its concrete implementation to store its contents to a unified storage (``vector<shared_ptr<Interaction> >``) before serialization and to restore from it after deserialization.
* precomputing non-serialized attributes from the serialized values; e.g. :yref:`Facet` computes its (local) edge normals and edge lengths from vertices' coordinates.


.. _classfactory:

Class factory
^^^^^^^^^^^^^^
Each serializable class must use ``REGISTER_SERIALIZABLE``, which defines function to create that class by ``ClassFactory``. ``ClassFactory`` is able to instantiate a class given its name (as string), which is necessary for deserialization.

Although mostly used internally by the serialization framework, programmer can ask for a class instantiation using ``shared_ptr<Factorable> f=ClassFactory::instance().createShared("ClassName");``, casting the returned ``shared_ptr<Factorable>`` to desired type afterwards. :yref:`Serializable` itself derives from ``Factorable``, i.e. all serializable types are also factorable (It is possible that different mechanism will be in place if boost::serialization is used, though.)

.. _plugins:

Plugin registration
^^^^^^^^^^^^^^^^^^^
Yade loads dynamic libraries containing all its functionality at startup. ClassFactory must be taught about classes each particular file provides. ``YADE_PLUGIN`` serves this purpose and, contrary to :ref:`YADE_CLASS_BASE_DOC`, must be place in the implementation (.cpp) file. It simple enumerates classes that are provided by this file::

	YADE_PLUGIN((ClassFoo)(ClassBar));

.. note:: You must use parentheses around the class name even if there is only one (preprocessor limitation): ``YADE_PLUGIN((classFoo));``. If there is no class in this file, do not use this macro at all.

Internally, this macro creates function ``registerThisPluginClasses_`` declared specially as ``__attribute__((constructor))`` (see `GCC Function Attributes <http://gcc.gnu.org/onlinedocs/gcc/Function-Attributes.html>`_); this attributes makes the function being executed when the plugin is loaded via ``dlopen`` from ``ClassFactory::load(...)``. It registers all factorable classes from that file in the :ref:`classfactory`.

.. note::
	Classes that do not derive from ``Factorable``, such as ``Shop`` or ``SpherePack``, are not declared with ``YADE_PLUGIN``.

---------------

This is an example of a serializable class header:

.. code-block:: c++

	/*! Homogeneous gravity field; applies gravity×mass force on all bodies. */
	class GravityEngine: public GlobalEngine{
		public:
			virtual void action();
		// registering class and its base for the RTTI system
		YADE_CLASS_BASE_DOC_ATTRS(GravityEngine,GlobalEngine,
			// documentation visible from python and generated reference documentation
			"Homogeneous gravity field; applies gravity×mass force on all bodies.",
			// enumerating attributes here, include documentation
			((Vector3r,gravity,Vector3r::ZERO,"acceleration, zero by default [kgms⁻²]"))
		);
	};
	// registration function for ClassFactory
	REGISTER_SERIALIZABLE(GravityEngine);

and this is the implementation:

.. code-block:: c++

	#include<yade/pkg-common/GravityEngine.hpp>
	#include<yade/core/Scene.hpp>

	// registering the plugin
	YADE_PLUGIN((GravityEngine));

	void GravityEngine::action(){
		/* do the work here */
	}

We can create a mini-simulation (with only one GravityEngine):

.. ipython::

	Yade [1]: O.engines=[GravityEngine(gravity=Vector3(0,0,-9.81))]

	Yade [2]: O.save('abc.xml')


and the XML looks like this:

.. literalinclude:: abc.xml
	:language: xml


.. warning:: Since XML files closely reflect implementation details of Yade, they will not be compatible between different versions. Use them only for short-term saving of scenes. Python is *the* high-level description Yade uses.

.. _pythonattributeaccess:

Python attribute access
^^^^^^^^^^^^^^^^^^^^^^^^

The macro :ref:`YADE_CLASS_BASE_DOC` introduced above is (behind the scenes) also used to create functions for accessing attributes from Python. As already noted, set of serialized attributes and set of attributes accessible from Python are identical. Besides attribute access, these wrapper classes imitate also some functionality of regular python dictionaries:

.. ipython::
	
	Yade [1]: s=Sphere()

	Yade [2]: s.radius              ## read-access

	Yade [3]: s.radius=4.           ## write access

	Yade [4]: s.dict().keys()              ## show all available keys

	Yade [5]: for k in s.dict().keys(): print s.dict()[k]  ## iterate over keys, print their values
	   ...:

	Yade [5]: s.dict()['radius']             ## same as: 'radius' in s.keys()

	Yade [6]: s.dict()                       ## show dictionary of both attributes and values


.. _YADE_CLASS_BASE_DOC:

YADE_CLASS_BASE_DOC_* macro family
-----------------------------------

There is several macros that hide behind them the functionality of :ref:`sphinxdocumentation`, :ref:`rtti`, :ref:`attributeregistration`, :ref:`pythonattributeaccess`, plus automatic attribute initialization and documentation. They are all defined as shorthands for base macro ``YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY`` with some arguments left out. They must be placed in class declaration's body (``.hpp`` file):

.. code-block:: c++

	#define YADE_CLASS_BASE_DOC(klass,base,doc) \
		YADE_CLASS_BASE_DOC_ATTRS(klass,base,doc,)
	#define YADE_CLASS_BASE_DOC_ATTRS(klass,base,doc,attrs) \
		YADE_CLASS_BASE_DOC_ATTRS_CTOR(klass,base,doc,attrs,)
	#define YADE_CLASS_BASE_DOC_ATTRS_CTOR(klass,base,doc,attrs,ctor) \
		YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(klass,base,doc,attrs,ctor,)
	#define YADE_CLASS_BASE_DOC_ATTRS_CTOR_PY(klass,base,doc,attrs,ctor,py) \
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(klass,base,doc,attrs,,ctor,py)
	#define YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(klass,base,doc,attrs,init,ctor,py) \
		YADE_CLASS_BASE_DOC_ATTRS_DEPREC_INIT_CTOR_PY(klass,base,doc,attrs,,init,ctor,py)


Expected parameters are indicated by macro name components separated with underscores. Their meaning is as follows:

``klass``
	(unquoted) name of this class (used for RTTI and python) 
``base``
	(unquoted) name of the base class (used for RTTI and python)
``doc``
	docstring of this class, written in the ReST syntax. This docstring will appear in generated documentation (such as :yref:`CpmMat`). It can be as long as necessary, but sequences interpreted by c++ compiler must be properly escaped (therefore some backslashes must be doubled, like in :math:`\sigma=\epsilon E`::
	
		":math:`\\sigma=\\epsilon E"
		
	Use ``\n`` and ``\t`` for indentation inside the docstring. Hyperlink the documentation abundantly with ``yref`` (all references to other classes should be hyperlinks).

	See :ref:`sphinxdocumentation` for syntax details.
``attrs``
	Attribute must be written in the form of parethesized list:

	.. code-block:: c++

		((type1,attr1,initValue1,attrFlags,"Attribute 1 documentation"))
		((type2,attr2,,,"Attribute 2 documentation"))  // initValue and attrFlags unspecified

	This will expand to
	
	#. data members declaration in c++ (note that all attributes are *public*):

		.. code-block:: c++

			public: type1 attr1;
				type2 attr2;

	#. Initializers of the default (argument-less) constructor, for attributes that have non-empty ``initValue``:

		.. code-block:: c++

			Klass(): attr1(initValue1), attr2() { /* constructor body */ }

		No initial value will be assigned for attribute of which initial value is left empty (as is for attr2 in the above example). Note that you still have to write the commas.

	#. Registration of the attribute in the serialization system (unless disabled by attrFlags -- see below)

	#. Registration of the attribute in python (unless disabled by attrFlags), so that it can be accessed as ``klass().name1``.
		The attribute is read-write by default, see attrFlags to change that.
	
		This attribute will carry the docstring provided, along with knowledge of the initial value. You can add text description to the default value using the comma operator of c++ and casting the char* to (void):

		.. code-block:: c++

			((Real,dmgTau,((void)"deactivated if negative",-1),,"Characteristic time for normal viscosity. [s]"))

		leading to :yref:`CpmMat::dmgTau`.

		The attribute is registered via ``boost::python::add_property`` specifying ``return_by_value`` policy rather than ``return_internal_reference``, which is the default when using ``def_readwrite``. The reason is that we need to honor custom converters for those values; see note in :ref:`customconverters` for details.

	.. admonition:: Attribute flags

		By default, an attribute will be serialized and will be read-write from python. There is a number of flags that can be passed as the 4th argument (empty by default) to change that:

		* ``Attr::noSave`` avoids serialization of the attribute (while still keeping its accessibility from Python)
		* ``Attr::readonly`` makes the attribute read-only from Python
		* ``Attr::triggerPostLoad`` will trigger call to ``postLoad`` function to handle attribute change after its value is set from Python; this is to ensure consistency of other precomputed data which depend on this value (such as ``Cell.trsf`` and such)
		* ``Attr::hidden`` will not expose the attribute to Python at all
		* ``Attr::noResize`` will not permit changing size of the array from Python [not yet used]

		Flags can be combined as usual using bitwise disjunction ``|`` (such as ``Attr::noSave | Attr::readonly``), though in such case the value should be parenthesized to avoid a warning with some compilers (g++ specifically), i.e. ``(Attr::noSave | Attr::readonly)``.

		Currently, the flags logic handled at runtime; that means that even for attributes with ``Attr::noSave``, their serialization template must be defined (although it will never be used). In the future, the implementation might be template-based, avoiding this necessity.


``deprec``
	List of deprecated attribute names. The syntax is ::

		((oldName1,newName1,"Explanation why renamed etc."))
		((oldName2,newName2,"! Explanation why removed and what to do instaed."))

	This will make accessing ``oldName1`` attribute *from Python* return value of ``newName``, but displaying warning message about the attribute name change, displaying provided explanation. This happens whether the access is read or write.

	If the explanation's first character is ``!`` (*bang*), the message will be displayed upon attribute access, but exception will be thrown immediately. Use this in cases where attribute is no longer meaningful or was not straightforwardsly replaced by another, but more complex adaptation of user's script is needed. You still have to give ``newName2``, although its value will never be used -- you can use any variable you like, but something must be given for syntax reasons).

	.. warning::
		Due to compiler limitations, this feature only works if Yade is compiled with gcc >= 4.4. In the contrary case, deprecated attribute functionality is disabled, even if such attributes are declared.

``init``
	Parethesized list of the form:

	.. code-block:: c++

		((attr3,value3)) ((attr4,value4))

	which will be expanded to initializers in the default ctor:

	.. code-block:: c++

		Klass(): /* attributes declared with the attrs argument */ attr4(value4), attr5(value5) { /* constructor body */ }

	The purpose of this argument is to make it possible to initialize constants and references (which are not declared as attributes using this macro themselves, but separately), as that cannot be done in constructor body. This argument is rarely used, though.
		
``ctor``
	will be put directly into the generated constructor's body. Mostly used for calling createIndex(); in the constructor.

	.. note:: 
		The code must not contain commas ouside parentheses (since preprocessor uses commas to separate macro arguments). If you need complex things at construction time, create a separate init() function and call it from the constructor instead.
``py``
	will be appeneded directly after generated python code that registers the class and all its attributes. You can use it to access class methods from python, for instance, to override an existing attribute with the same name etc:

	.. code-block:: c++

		.def_readonly("omega",&CpmPhys::omega,"Damage internal variable")
		.def_readonly("Fn",&CpmPhys::Fn,"Magnitude of normal force.")

	``def_readonly`` will not work for custom types (such as std::vector), as it bypasses conversion registry; see :ref:`customconverters` for details.


Special python constructors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The Python wrapper automatically create constructor that takes keyword (named) arguments corresponding to instance attributes; those attributes are set to values provided in the constructor. In some cases, more flexibility is desired (such as :yref:`InteractionLoop`, which takes 3 lists of functors). For such cases, you can override the function ``Serializable::pyHandleCustomCtorArgs``, which can arbitrarily modify the new (already existing) instance. It should modify in-place arguments given to it, as they will be passed further down to the routine which sets attribute values. In such cases, you should document the constructor::

	.. admonition:: Special constructor

		Constructs from lists of …

which then appears in the documentation similar to :yref:`InteractionLoop`.

Static attributes
^^^^^^^^^^^^^^^^^^^

Some classes (such as OpenGL functors) are instantiated automatically; since we want their attributes to be persistent throughout the session, they are static. To expose class with static attributes, use the ``YADE_CLASS_BASE_DOC_STATICATTRS`` macro. Attribute syntax is the same as for ``YADE_CLASS_BASE_DOC_ATTRS``: 
	
.. code-block:: c++

	class SomeClass: public BaseClass{
		YADE_CLASS_BASE_DOC_STATICATTRS(SomeClass,BaseClass,"Documentation of SomeClass",
			((Type1,attr1,default1,"doc for attr1"))
			((Type2,attr2,default2,"doc for attr2"))
		);
	};

additionally, you *have* to allocate memory for static data members in the ``.cpp`` file (otherwise, error about undefined symbol will appear when the plugin is loaded):

.. code-block: c++

	/* in the .cpp file */
	#include<yade/pkg-something/HeaderFile.hpp>
	/* allocate memory for static attrs;
	   no need to assign initial value, that is done from
	   (hidden) initialization function at class registration
	   time using default values given to macro in the .hpp file.
	*/
	Type1 SomeClass::attr1;
	Type2 SomeClass::attr2;

There is no way to expose class that has both static and non-static attributes using ``YADE_CLASS_BASE_*`` macros. You have to expose non-static attributes normally and wrap static attributes separately in the ``py`` parameter.


.. _valuereference:

Returning attribute by value or by reference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When attribute is passed from c++ to python, it can be passed either as 

* value: new python object representing the original c++ object is constructed, but not bound to it; changing the python object doesn't modify the c++ object, unless explicitly assigned back to it, where inverse conversion takes place and the c++ object is replaced.
* reference: only reference to the underlying c++ object is given back to python; modifying python object will make the c++ object modified automatically.

The way of passing attributes given to ``YADE_CLASS_BASE_DOC_ATTRS`` in the ``attrs`` parameter is determined automatically in the following manner:

* ``Vector3``, ``Vector3i``, ``Vector2``, ``Vector2i``, ``Matrix3`` and ``Quaternion`` objects are passed by *reference*. For instance::
		O.bodies[0].state.pos[0]=1.33
  will assign correct value to ``x`` component of position, without changing the other ones.
* Yade classes (all that use ``shared_ptr`` when declared in python: all classes deriving from :yref:`Serializable` declared with ``YADE_CLASS_BASE_DOC_*``, and some others) are passed as *references* (technically speaking, they are passed by value of the ``shared_ptr``, but by virtue of its sharedness, they appear as references). For instance::
		O.engines[4].damping=.3
  will change :yref:`damping<NewtonIntegrator.damping>` parameter on the original engine object, not on its copy.
* All other types are passed by *value*. This includes, most importantly, sequence types declared in :ref:`customconverters`, such as ``std::vector<shared_ptr<Engine> >``. For this reason, ::
		O.engines[4]=NewtonIntegrator()
  will *not* work as expected; it will replace 5th element of a *copy* of the sequence, and this change will not propagate back to c++.




.. _multiple-dispatch:

Multiple dispatch
------------------
Multiple dispatch is generalization of virtual methods: a :yref:`Dispatcher` decides based on type(s) of its argument(s) which of its :yref:`Functors<Functor>` to call. Numer of arguments (currently 1 or 2) determines *arity* of the dispatcher (and of the functor): unary or binary. For example:

.. code-block:: python

	InsertionSortCollider([Bo1_Sphere_Aabb(),Bo1_Facet_Aabb()])

creates :yref:`InsertionSortCollider`, which internally contains :yref:`Collider.boundDispatcher`, a :yref:`BoundDispatcher` (a :yref:`Dispatcher`), with 2 functors; they receive ``Sphere`` or ``Facet`` instances and create ``Aabb``. This code would look like this in c++:

.. code-block:: c++

	shared_ptr<InsertionSortCollider> collider=(new InsertionSortCollider);
	collider->boundDispatcher->add(new Bo1_Sphere_Aabb());
	collider->boundDispatcher->add(new Bo1_Facet_Aabb());

There are currenly 4 predefined dispatchers (see `dispatcher-names`_) and corresponding functor types. They are inherit from template instantiations of ``Dispatcher1D`` or ``Dispatcher2D`` (for functors, ``Functor1D`` or ``Functor2D``). These templates themselves derive from ``DynlibDispatcher`` (for dispatchers) and ``FunctorWrapper`` (for functors).

Example: IGeomDispatcher
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let's take (the most complicated perhaps) :yref:`IGeomDispatcher`. :yref:`IGeomFunctor`, which is dispatched based on types of 2 :yref:`Shape` instances (a :yref:`Functor`), takes a number of arguments and returns bool. The functor "call" is always provided by its overridden ``Functor::go`` method; it always receives the dispatched instances as first argument(s) (2 × ``const shared_ptr<Shape>&``) and a number of other arguments it needs:

.. code-block:: c++

	class IGeomFunctor: public Functor2D<
	   bool,                                 //return type
	   TYPELIST_7(const shared_ptr<Shape>&,  // 1st class for dispatch 
	      const shared_ptr<Shape>&,          // 2nd class for dispatch
	      const State&,                      // other arguments passed to ::go
	      const State&,                      // …
	      const Vector3r&,                   // … 
	      const bool&,                       // …
	      const shared_ptr<Interaction>&     // …
	   )
	> 

The dispatcher is declared as follows:

.. code-block:: c++

	class IGeomDispatcher: public Dispatcher2D<	
	   Shape,                       // 1st class for dispatch
	   Shape,                       // 2nd class for dispatch
	   IGeomFunctor,  // functor type
	   bool,                        // return type of the functor

	   // follow argument types for functor call
	   // they must be exactly the same as types
	   // given to the IGeomFunctor above.
	   TYPELIST_7(const shared_ptr<Shape>&,  
	      const shared_ptr<Shape>&,
	      const State&,
	      const State&,
	      const Vector3r&,
	      const bool &,
	      const shared_ptr<Interaction>&
	   ),

	   // handle symetry automatically
	   // (if the dispatcher receives Sphere+Facet,
	   // the dispatcher might call functor for Facet+Sphere,
	   // reversing the arguments)
	   false
	>
	{ /* … */ }

Functor derived from IGeomFunctor must then 

* override the ::go method with appropriate arguments (they must match exactly types given to ``TYPELIST_*`` macro);
* declare what types they should be dispatched for, and in what order if they are not the same.

.. code-block:: c++

	class Ig2_Facet_Sphere_Dem3DofGeom: public IGeomFunctor{
	   public:

	   // override the IGeomFunctor::go
	   //   (it is really inherited from FunctorWrapper template,
	   //    therefore not declare explicitly in the
	   //    IGeomFunctor declaration as such)
	   // since dispatcher dispatches only for declared types
	   //   (or types derived from them), we can do 
	   //   static_cast<Facet>(shape1) and static_cast<Sphere>(shape2)
	   //   in the ::go body, without worrying about types being wrong.
	   virtual bool go(
	      // objects for dispatch
	      const shared_ptr<Shape>& shape1, const shared_ptr<Shape>& shape2,
	      // other arguments
	      const State& state1, const State& state2, const Vector3r& shift2,
	      const bool& force, const shared_ptr<Interaction>& c
	   );
	   /* … */

	   // this declares the type we want to be dispatched for, matching
	   //   first 2 arguments to ::go and first 2 classes in TYPELIST_7 above
	   //   shape1 is a Facet and shape2 is a Sphere
	   //   (or vice versa, see lines below)
	   FUNCTOR2D(Facet,Sphere);

	   // declare how to swap the arguments
	   //   so that we can receive those as well
	   DEFINE_FUNCTOR_ORDER_2D(Facet,Sphere);
	   /* … */
	};


Dispatch resolution
^^^^^^^^^^^^^^^^^^^
The dispatcher doesn't always have functors that exactly match the actual types it receives. In the same way as virtual methods, it tries to find the closest match in such way that:

#. the actual instances are derived types of those the functor accepts, or exactly the accepted types;
#. sum of distances from actual to accepted types is sharp-minimized (each step up in the class hierarchy counts as 1)

If no functor is able to accept given types (first condition violated) or multiple functors have the same distance (in condition 2), an exception is thrown.

This resolution mechanism makes it possible, for instance, to have a hierarchy of :yref:`Dem3DofGeom` classes (for different combination of shapes: :yref:`Dem3DofGeom_SphereSphere`, :yref:`Dem3DofGeom_FacetSphere`, :yref:`Dem3DofGeom_WallSphere`), but only provide a :yref:`LawFunctor` accepting ``Dem3DofGeom``, rather than having different laws for each shape combination.

.. note:: Performance implications of dispatch resolution are relatively low. The dispatcher lookup is only done once, and uses fast lookup matrix (1D or 2D); then, the functor found for this type(s) is cached within the ``Interaction`` (or ``Body``) instance. Thus, regular functor call costs the same as dereferencing pointer and calling virtual method. There is `blueprint <https://blueprints.launchpad.net/yade/+spec/devirtualize-functor-calls>`_ to avoid virtual function call as well.

.. note:: At the beginning, the dispatch matrix contains just entries exactly matching given functors. Only when necessary (by passing other types), appropriate entries are filled in as well.

Indexing dispatch types
^^^^^^^^^^^^^^^^^^^^^^^^

Classes entering the dispatch mechanism must provide for fast identification of themselves and of their parent class. [#rttiindex]_ This is called class indexing and all such classes derive from :yref:`Indexable`. There are ``top-level`` Indexables (types that the dispatchers accept) and each derived class registers its index related to this top-level Indexable. Currently, there are:

==================== ===========================
Top-level Indexable  used by
==================== ===========================
:yref:`Shape`        :yref:`BoundFunctor`, :yref:`IGeomDispatcher`
:yref:`Material`     :yref:`IPhysDispatcher`
:yref:`IPhys`        :yref:`LawDispatcher`
:yref:`IGeom`        :yref:`LawDispatcher`
==================== ===========================

The top-level Indexable must use the ``REGISTER_INDEX_COUNTER`` macro, which sets up the machinery for identifying types of derived classes; they must then use the ``REGISTER_CLASS_INDEX`` macro *and* call ``createIndex()`` in their constructor. For instance, taking the :yref:`Shape` class (which is a top-level Indexable):

.. code-block:: c++

	// derive from Indexable
	class Shape: public Serializable, public Indexable {  
	   // never call createIndex() in the top-level Indexable ctor!
	   /* … */

	   // allow index registration for classes deriving from ``Shape``
	   REGISTER_INDEX_COUNTER(Shape);
	};

Now, all derived classes (such as :yref:`Sphere` or :yref:`Facet`) use this:

.. code-block:: c++

	class Sphere: public Shape{
	   /* … */
	   YADE_CLASS_BASE_DOC_ATTRS_CTOR(Sphere,Shape,"docstring",
	      ((Type1,attr1,default1,"docstring1"))
	      /* … */,
	      // this is the CTOR argument
	         // important; assigns index to the class at runtime
	         createIndex(); 
	   );
	   // register index for this class, and give name of the immediate parent class
	   //    (i.e. if there were a class deriving from Sphere, it would use                             
	   //     REGISTER_CLASS_INDEX(SpecialSphere,Sphere),
	   //     not REGISTER_CLASS_INDEX(SpecialSphere,Shape)!)
	   REGISTER_CLASS_INDEX(Sphere,Shape);   
	};

At runtime, each class within the top-level Indexable hierarchy has its own unique numerical index. These indices serve to build the dispatch matrix for each dispatcher.

.. [#rttiindex] The functionality described in :ref:`rtti` serves a different purpose (serialization) and would hurt the performance here. For this reason, classes provide numbers (indices) in addition to strings.

Inspecting dispatch in python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If there is a need to debug/study multiple dispatch, python provides convenient interface for this low-level functionality.

We can inspect indices with the ``dispIndex`` property (note that the top-level indexable ``Shape`` has negative (invalid) class index; we purposively didn't call ``createIndex`` in its constructor):

.. ipython::
	
	Yade [5]: Sphere().dispIndex, Facet().dispIndex, Wall().dispIndex

	Yade [6]: Shape().dispIndex                  # top-level indexable

Dispatch hierarchy for a particular class can be shown with the ``dispHierarchy()`` function, returning list of class names: 0th element is the instance itself, last element is the top-level indexable (again, with invalid index); for instance:

.. ipython::

	Yade [7]: Dem3DofGeom().dispHierarchy()       # parent class of all other Dem3DofGeom_ classes

	Yade [8]: Dem3DofGeom_SphereSphere().dispHierarchy(), Dem3DofGeom_FacetSphere().dispHierarchy(), Dem3DofGeom_WallSphere().dispHierarchy()

	Yade [8]: Dem3DofGeom_WallSphere().dispHierarchy(names=False)   # show numeric indices instead


Dispatchers can also be inspected, using the .dispMatrix() method:

.. ipython::

	Yade [3]: ig=IGeomDispatcher([
	   ...:    Ig2_Sphere_Sphere_Dem3DofGeom(),
	   ...:    Ig2_Facet_Sphere_Dem3DofGeom(),
	   ...:    Ig2_Wall_Sphere_Dem3DofGeom()
	   ...: ])

	Yade [4]: ig.dispMatrix()

	Yade [5]: ig.dispMatrix(False)          # don't convert to class names

We can see that functors make use of symmetry (i.e. that Sphere+Wall are dispatched to the same functor as Wall+Sphere).

Finally, dispatcher can be asked to return functor suitable for given argument(s):

.. ipython::

	Yade [6]: ld=LawDispatcher([Law2_Dem3DofGeom_CpmPhys_Cpm()])

	Yade [7]: ld.dispMatrix()

	# see how the entry for Dem3DofGeom_SphereSphere will be filled after this request

	Yade [8]: ld.dispFunctor(Dem3DofGeom_SphereSphere(),CpmPhys())       

	Yade [9]: ld.dispMatrix()


OpenGL functors
^^^^^^^^^^^^^^^
OpenGL rendering is being done also by 1D functors (dispatched for the type to be rendered). Since it is sufficient to have exactly one class for each rendered type, the functors are found automatically. Their base functor types are ``GlShapeFunctor``, ``GlBoundFunctor``, ``GlIGeomFunctor`` and so on. These classes register the type they render using the ``RENDERS`` macro:

.. code-block:: c++

	class Gl1_Sphere: public GlShapeFunctor {
	   public :
	      virtual void go(const shared_ptr<Shape>&,
	         const shared_ptr<State>&,
	         bool wire,
	         const GLViewInfo&
	      );
	   RENDERS(Sphere);
	   YADE_CLASS_BASE_DOC_STATICATTRS(Gl1_Sphere,GlShapeFunctor,"docstring",
	      ((Type1,staticAttr1,informativeDefault,"docstring"))
	      /* … */
	   );
	};
	REGISTER_SERIALIZABLE(Gl1_Sphere);

You can list available functors of a particular type by querying child classes of the base functor:

.. ipython::

	Yade [1]: import yade.system

	Yade [2]: yade.system.childClasses('GlShapeFunctor')

.. note:: OpenGL functors may disappear in the future, being replaced by virtual functions of each class that can be rendered.


Parallel execution
------------------

Yade was originally not designed with parallel computation in mind, but rather with maximum flexibility (for good or for bad). Parallel execution was added later; in order to not have to rewrite whole Yade from scratch, relatively non-instrusive way of parallelizing was used: `OpenMP <http://www.openmp.org>`__. OpenMP is standartized shared-memory parallel execution environment, where parallel sections are marked by special ``#pragma`` in the code (which means that they can compile with compiler that doesn't support OpenMP) and a few functions to query/manipulate OpenMP runtime if necessary.

There is parallelism at 3 levels:

* Computation, interaction (python, GUI) and rendering threads are separate. This is done via regular threads (boost::threads) and is not related to OpenMP.
* :yref:`ParallelEngine` can run multiple engine groups (which are themselves run serially) in parallel; it rarely finds use in regular simulations, but it could be used for example when coupling with an independent expensive computation:

	.. code-block:: python

		ParallelEngine([
			[Engine1(),Engine2()],   # Engine1 will run before Engine2
			[Engine3()]              # Engine3() will run in parallel with the group [Engine1(),Engine2()]
			                         # arbitrary number of groups can be used
		])

	``Engine2`` will be run after ``Engine1``, but in parallel with ``Engine3``.

	.. warning:: It is your reponsibility to avoid concurrent access to data when using ParallelEngine. Make sure you understand *very well* what the engines run in parallel do.
* Parallelism inside Engines. Some loops over bodies or interactions are parallelized (notably :yref:`InteractionLoop` and :yref:`NewtonIntegrator`, which are treated in detail later (FIXME: link)):

	.. code-block:: c++

		#pragma omp parallel for
		for(long id=0; id<size; id++){
		   const shared_ptr<Body>& b(scene->bodies[id]);
		   /* … */
		}

	.. note :: OpenMP requires loops over contiguous range of integers (OpenMP 3 also accepts containers with random-access iterators).

	If you consider running parallelized loop in your engine, always evalue its benefits. OpenMP has some overhead fo creating threads and distributing workload, which is proportionally more expensive if the loop body execution is fast. The results are highly hardware-dependent (CPU caches, RAM controller).

Maximum number of OpenMP threads is determined by the ``OMP_NUM_THREADS`` environment variable and is constant throughout the program run. Yade main program also sets this variable (before loading OpenMP libraries) if you use the ``-j``/``--threads`` option. It can be queried at runtime with the ``omp_get_max_threads`` function.

At places which are susceptible of being accessed concurrently from multiple threads, Yade provides some mutual exclusion mechanisms, discussed elsewhere (FIXME):

* simultaneously writeable container for :ref:`ForceContainer`,
* mutex for :yref:`Body::state`.

.. _logging:

Logging
--------

Regardless of whether the :ref:`optional-libraries` log4cxx is used or not, yade provides logging macros. [#log4cxxup]_ If log4cxx is enabled, these macros internally operate on the local logger instance (named ``logger``, but that is hidden for the user); if log4cxx is disabled, they send their arguments to standard error output (``cerr``).

.. [#log4cxxup] Because of (seemingly?) no upstream development of log4cxx and a few problems it has, Yade will very likely move to the hypothetical ``boost::logging`` library once it exists. The logging code will not have to be changed, however, as the log4cxx logic is hidden behind these macros.

Log messages are classified by their *severity*, which is one of ``TRACE`` (tracing variables), ``DEBUG`` (generally uninsteresting messages useful for debugging), ``INFO`` (information messages -- only use sparingly), ``WARN`` (warning), ``FATAL`` (serious error, consider throwing an exception with description instead). Logging level determines which messages will be shown -- by default, ``INFO`` and higher will be shown; if you run yade with ``-v`` or ``-vv``, ``DEBUG`` and ``TRACE`` messages will be also enabled (with log4cxx).

Every class using logging should create logger using these 2 macros (they expand to nothing if log4cxx is not used):

``DECLARE_LOGGER;``
	in class declaration body (in the ``.hpp`` file); this declares static variable ``logger``;
``CREATE_LOGGER(ClassName);``
	in the implementation file; it creates and initializes that static variable. The logger will be named ``yade.ClassName``.

The logging macros are the following:

* ``LOG_TRACE``, ``LOG_DEBUG``, ``LOG_INFO``, ``LOG_WARN``, ``LOG_ERROR``, ``LOG_FATAL`` (increasing severity); their argument is fed to the logger stream, hence can contain the ``<<`` operation:

	.. code-block:: c++

		LOG_WARN("Exceeded "<<maxSteps<<" steps in attempts to converge, the result returned will not be precise (relative error "<<relErr<<", tolerance "<<relTol<<")");

	Every log message is prepended filename, line number and function name; the final message that will appear will look like this::

		237763 WARN  yade.ViscosityIterSolver /tmp/yade/trunk/extra/ViscosityIterSolver.cpp:316 newtonRaphsonSolve: Exceeded 30 steps in attempts to converge, the result returned will not be precise (relative error 5.2e-3, tolerance 1e-3)

	The ``237763 WARN  yade.ViscosityIterSolver`` (microseconds from start, severity, logger name) is added by log4cxx and is completely configurable, either programatically, or by using file ``~/.yade-\$SUFFIX/logging.conf``, which is loaded at startup, if present (FIXME: see more etc user's guide)


* special tracing macros ``TRVAR1``, ``TRVAR2``, … ``TRVAR6``, which show both variable name and its value (there are several more macros defined inside ``/lib/base/Logging.hpp``, but they are not generally in use):

	.. code-block:: c++

		TRVAR3(var1,var2,var3);
		// will be expanded to:
		LOG_TRACE("var1="<<var1<<"; var2="<<var2<<"; var3="<<var3);


.. note:: For performance reasons, optimized builds eliminate ``LOG_TRACE`` and ``LOG_DEBUG`` from the code at preprocessor level.

.. note:: Builds without log4cxx (even in debug mode) eliminate ``LOG_TRACE`` and ``LOG_DEBUG``. As there is no way to enable/disable them selectively, the log amount would be huge.

Python provides rudimentary control for the logging system in ``yade.log`` module (FIXME: ref to docs):

.. ipython::

	Yade [2]: from yade import log

	Yade [3]: log.setLevel('InsertionSortCollider',log.DEBUG)  # sets logging level of the yade.InsertionSortCollider logger

	Yade [4]: log.setLevel('',log.WARN)                        # sets logging level of all yade.* loggers (they inherit level from the parent logger, except when overridden)

As of now, there is no python interface for performing logging into log4cxx loggers themselves.

.. _timing:

Timing
-------

Yade provides 2 services for measuring time spent in different pars of the code. One has the granularity of engine and can be enabled at runtime. The other one is finer, but requires adjusting and recompiling the code being measured.

Per-engine timing
^^^^^^^^^^^^^^^^^^
The coarser timing works by merely accumulating numebr of invocations and time (with the precision of the ``clock_gettime`` function) spent in each engine, which can be then post-processed by associated Python module ``yade.timing``. There is a static bool variable controlling whether such measurements take place (disabled by default), which you can change

.. code-block:: c++

	TimingInfo::enabled=True;            // in c++

.. code-block:: python

	O.timingEnabled=True                 ## in python

After running the simulation, ``yade.timing.stats()`` function will show table with the results and percentages:

.. ipython::

	Yade [1]: TriaxialTest(numberOfGrains=100).load()

	Yade [2]: O.engines[0].label='firstEngine'   ## labeled engines will show by labels in the stats table

	Yade [2]: import yade.timing;

	Yade [2]: O.timingEnabled=True
	
	Yade [1]: yade.timing.reset()                   ## not necessary if used for the first time

	Yade [3]: O.run(50); O.wait()

	Yade [5]: yade.timing.stats()

Exec count and time can be accessed and manipulated through ``Engine::timingInfo`` from c++ or ``Engine().execCount`` and ``Engine().execTime`` properties in Python.

In-engine and in-functor timing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Timing within engines (and functors) is based on :yref:`TimingDeltas` class. It is made for timing loops (functors' loop is in their respective dispatcher) and stores cummulatively time differences between *checkpoints*.

.. note:: Fine timing with ``TimingDeltas`` will only work if timing is enabled globally (see previous section). The code would still run, but giving zero times and exec counts.

#. Engine::timingDeltas must point to an instance of :yref:`TimingDeltas` (prefferably instantiate :yref:`TimingDeltas` in the constructor):

	.. code-block:: c++
		
		// header file
		class Law2_Dem3DofGeom_CpmPhys_Cpm: public LawFunctor {
		   /* … */
		   YADE_CLASS_BASE_DOC_ATTRS_CTOR(Law2_Dem3DofGeom_CpmPhys_Cpm,LawFunctor,"docstring",
		      /* attrs */,
		      /* constructor */
		      timingDeltas=shared_ptr<TimingDeltas>(new TimingDeltas);
		   );
		   // ...
		};

#. Inside the loop, start the timing by calling ``timingDeltas->start();``
#. At places of interest, call ``timingDeltas->checkpoint("label")``. The label is used only for post-processing, data are stored based on the checkpoint position, not the label.

	.. warning:: Checkpoints must be always reached in the same order, otherwise the timing data will be garbage. Your code can still branch, but you have to put checkpoints to places which are in common.

	.. code-block:: c++

		void Law2_Dem3DofGeom_CpmPhys_Cpm::go(shared_ptr<IGeom>& _geom,
		                                      shared_ptr<IPhys>& _phys,
		                                      Interaction* I,
		                                      Scene* scene)
		{
		   timingDeltas->start();                        // the point at which the first timing starts
		   // prepare some variables etc here
		   timingDeltas->checkpoint("setup");
		   // find geometrical data (deformations) here
		   timingDeltas->checkpoint("geom");
		   // compute forces here
		   timingDeltas->checkpoint("material");
		   // apply forces, cleanup here
		   timingDeltas->checkpoint("rest");
		}

The output might look like this (note that functors are nested inside dispatchers and ``TimingDeltas`` inside their engine/functor)::

	Name                                    Count                 Time          Rel. time
	-------------------------------------------------------------------------------------
	ForceReseter                        400               9449μs              0.01%      
	BoundDispatcher                     400            1171770μs              1.15%      
	InsertionSortCollider               400            9433093μs              9.24%      
	IGeomDispatcher       400           15177607μs             14.87%      
	IPhysDispatcher        400            9518738μs              9.33%      
	LawDispatcher                       400           64810867μs             63.49%      
	  Law2_Dem3DofGeom_CpmPhys_Cpm                                                     
	    setup                           4926145            7649131μs             15.25%  
	    geom                            4926145           23216292μs             46.28%  
	    material                        4926145            8595686μs             17.14%  
	    rest                            4926145           10700007μs             21.33%  
	    TOTAL                                             50161117μs            100.00%  
	NewtonIntegrator                    400            1866816μs              1.83%      
	"strainer"                          400              21589μs              0.02%      
	"plotDataCollector"                 160              64284μs              0.06%      
	"damageChecker"                       9               3272μs              0.00%      
	TOTAL                                            102077490μs            100.00%      

.. warning:: Do not use :yref:`TimingDeltas` in parallel sections, results might not be meaningful. In particular, avoid timing functors inside :yref:`InteractionLoop` when running with multiple OpenMP threads.

``TimingDeltas`` data are accessible from Python as list of (*label*,*time*,*count*) tuples, one tuple representing each checkpoint:

.. code-block:: python

	deltas=someEngineOrFunctor.timingDeltas.data()
	deltas[0][0] # 0th checkpoint label
	deltas[0][1] # 0th checkpoint time in nanoseconds
	deltas[0][2] # 0th checkpoint execution count
	deltas[1][0] # 1st checkpoint label
	             # …
	deltas.reset() 

Timing overhead
^^^^^^^^^^^^^^^
The overhead of the coarser, per-engine timing, is very small. For simulations with at least several hundreds of elements, they are below the usual time variance (a few percent).

.. TriaxialTest(numberOfGrains=800).load(); O.saveTmp(); import time; O.timingEnabled=True; t0=time.time(); O.run(2000,True); print(time.time()-t0); O.loadTmp(); O.timingEnabled=False; t0=time.time(); O.run(2000,True); print (time.time()-t0)

The finer :yref:`TimingDeltas` timing can have major performance impact and should be only used during debugging and performance-tuning phase. The parts that are file-timed will take disproportionally longer time that the rest of engine; in the output presented above, :yref:`LawDispatcher` takes almost ⅔ of total simulation time in average, but the number would be twice of thrice lower typically (note that each checkpoint was timed almost 5 million times in this particular case).


OpenGL Rendering
-----------------
Yade provides 3d rendering based on `QGLViewer <http://www.libqglviewer.com>`__. It is not meant to be full-featured rendering and post-processing, but rather a way to quickly check that scene is as intended or that simulation behaves sanely. 

.. note:: Although 3d rendering runs in a separate thread, it has performance impact on the computation itself, since interaction container requires mutual exclusion for interaction creation/deletion. The ``InteractionContainer::drawloopmutex`` is either held by the renderer (:yref:`OpenGLRenderingEngine`) or by the insertion/deletion routine. 

.. warning:: There are 2 possible causes of crash, which are not prevented because of serious performance penalty that would result:

	#. access to :yref:`BodyContainer`, in particular deleting bodies from simulation; this is a rare operation, though.
	#. deleting Interaction::phys or Interaction::geom.

Renderable entities (:yref:`Shape`, :yref:`State`, :yref:`Bound`, :yref:`IGeom`, :yref:`IPhys`) have their associated `OpenGL functors`_. An entity is rendered if

#. Rendering such entities is enabled by appropriate attribute in :yref:`OpenGLRenderingEngine`
#. Functor for that particular entity type is found via the :ref:`dispatch mechanism<multiple-dispatch>`.

``Gl1_*`` functors operating on Body's attributes (:yref:`Shape`, :yref:`State`, :yref:`Bound`) are called with the OpenGL context translated and rotated according to :yref:`State::pos` and :yref:`State::ori`. Interaction functors work in global coordinates.




Simulation framework
======================

Besides the support framework mentioned in the previous section, some functionality pertaining to simulation itself is also provided.

There are special containers for storing bodies, interactions and (generalized) forces. Their internal functioning is normally opaque to the programmer, but should be understood as it can influence performance.

Scene
------

``Scene`` is the object containing the whole simulation. Although multiple scenes can be present in the memory, only one of them is active. Saving and loading (serializing and deserializing) the ``Scene`` object should make the simulation run from the point where it left off. 

.. note::
	All :yref:`Engines<Engine>` and functors have interally a ``Scene* scene`` pointer which is updated regularly by engine/functor callers; this ensures that the current scene can be accessed from within user code.

	For outside functions (such as those called from python, or static functions in ``Shop``), you can use ``Omega::instance().getScene()`` to retrieve a ``shared_ptr<Scene>`` of the current scene.

Body container
---------------
Body container is linear storage of bodies. Each body in the simulation has its unique :yref:`id<Body::id>`, under which it must be found in the :yref:`BodyContainer`. Body that is not yet part of the simulation typically has id equal to invalid value ``Body::ID_NONE``, and will have its ``id`` assigned upon insertion into the container. 
The requirements on :yref:`BodyContainer` are

* O(1) access to elements,
* linear-addressability (0…n indexability),
* store ``shared_ptr``, not objects themselves,
* *no* mutual exclusion for insertion/removal (this must be assured by the called, if desired),
* intelligent allocation of ``id`` for new bodies (tracking removed bodies),
* easy iteration over all bodies.

.. note::
	Currently, there is "abstract" class ``BodyContainer``, from which derive concrete implementations; the initial idea was the ability to select at runtime which implementation to use (to find one that performs the best for given simulation). This incurs the penalty of many virtual function calls, and will probably change in the future. All implementations of BodyContainer were removed in the meantime, except ``BodyVector`` (internally a  ``vector<shared_ptr<Body> >`` plus a few methods around), which is the fastest.

Insertion/deletion
^^^^^^^^^^^^^^^^^^

Body insertion is typically used in :yref:`FileGenerator`'s:

.. code-block:: c++

	shared_ptr<Body> body(new Body);
	// … (body setup)
	scene->bodies->insert(body); // assigns the id

Bodies are deleted only rarely:

.. code-block:: c++

	scene->bodies->erase(id);

.. warning::
	Since mutual exclusion is not assured, never insert/erase bodies from parallel sections, unless you explicitly assure there will be no concurrent access.


Iteration
^^^^^^^^^^
The container can be iterated over using ``FOREACH`` macro (shorthand for ``BOOST_FOREACH``):

.. code-block:: c++

	FOREACH(const shared_ptr<Body>& b, *scene->bodies){
	   if(!b) continue;                      // skip deleted bodies
	   /* do something here */              
	}

Note a few important things:

#. Always use ``const shared_ptr<Body>&`` (const reference); that avoids incrementing and decrementing the reference count on each ``shared_ptr``.
#. Take care to skip NULL bodies (``if(!b) continue``): deleted bodies are deallocated from the container, but since body id's must be persistent, their place is simply held by an empty ``shared_ptr<Body>()`` object, which is implicitly convertible to ``false``.

In python, the BodyContainer wrapper also has iteration capabilities; for convenience (which is different from the c++ iterator), NULL bodies as silently skipped:

.. ipython::

	Yade [0]: O.reset()

	Yade [1]: O.bodies.append([Body(),Body(),Body()])

	Yade [2]: O.bodies.erase(1)

	Yade [3]: [b.id for b in O.bodies]

In loops parallelized using OpenMP, the loop must traverse integer interval (rather than using iterators):

.. code-block:: c++

	const long size=(long)bodies.size();       // store this value, since it doesn't change during the loop
	#pragma omp parallel for
	for(long _id=0; _id<size; _id++){
	   const shared_ptr<Body>& b(bodies[_id]);
	   if(!b) continue;
	   /* … */
	}


InteractionContainer
---------------------
Interactions are stored in special container, and each interaction must be uniquely identified by pair of ids (id1,id2).

* O(1) access to elements,
* linear-addressability (0…n indexability),
* store ``shared_ptr``, not objects themselves,
* mutual exclusion for insertion/removal,
* easy iteration over all interactions,
* addressing symmetry, i.e. interaction(id1,id2)≡interaction(id2,id1)

.. note::
	As with BodyContainer, there is "abstract" class InteractionContainer, and then its concrete implementations. Currently, only InteractionVecMap implementation is used and all the other were removed. Therefore, the abstract InteractionContainer class may disappear in the future, to avoid unnecessary virtual calls. 

	Further, there is a `blueprint <https://blueprints.launchpad.net/yade/+spec/intrs-inside-bodies>`_ for storing interactions inside bodies, as that would give extra advantage of quickly getting all interactions of one particular body (currently, this necessitates loop over all interactions); in that case, InteractionContainer would disappear.

Insert/erase
^^^^^^^^^^^^
Creating new interactions and deleting them is delicate topic, since many eleents of simulation must be synchronized; the exact workflow is described in :ref:`interaction-flow`. You will almost certainly never need to insert/delete an interaction manually from the container; if you do, consider designing your code differently.

.. code-block:: c++

	// both insertion and erase are internally protected by a mutex,
	// and can be done from parallel sections safely
	scene->interactions->insert(shared_ptr<Interaction>(new Interactions(id1,id2)));
	scene->interactions->erase(id1,id2);

Iteration
^^^^^^^^^
As with BodyContainer, iteration over interactions should use the ``FOREACH`` macro:

.. code-block:: c++

	FOREACH(const shared_ptr<Interaction>& i, *scene->interactions){
	   if(!i->isReal()) continue;
	   /* … */
	}

Again, note the usage const reference for ``i``. The check ``if(!i->isReal())`` filters away interactions that exist only *potentially*, i.e. there is only :yref:`Bound` overlap of the two bodies, but not (yet) overlap of bodies themselves. The ``i->isReal()`` function is equivalent to ``i->geom && i->phys``. Details are again explained in :ref:`interaction-flow`.

In some cases, such as OpenMP-loops requiring integral index (OpenMP >= 3.0 allows parallelization using random-access iterator as well), you need to iterate over interaction indices instead:

.. code-block:: c++

	inr nIntr=(int)scene->interactions->size(); // hoist container size
	#pragma omp parallel for
	for(int j=0; j<nIntr, j++){
	   const shared_ptr<Interaction>& i(scene->interactions[j]);
	   if(!->isReal()) continue;
	   /* … */
	}

.. _ForceContainer:

ForceContainer
--------------

:yref:`ForceContainer` holds "generalized forces", i.e. forces, torques, (explicit) dispalcements and rotations for each body.

During each computation step, there are typically 3 phases pertaining to forces:

#. Resetting forces to zero (usually done by the :yref:`ForceResetter` engine)
#. Incrementing forces from parallel sections (solving interactions -- from :yref:`LawFunctor`)
#. Reading absolute force values sequentially for each body: forces applied from different interactions are summed together to give overall force applied on that body (:yref:`NewtonIntegrator`, but also various other engine that read forces)

This scenario leads to special design, which allows fast parallel write access:

* each thread has its own storage (zeroed upon request), and only writes to its own storage; this avoids concurrency issues. Each thread identifies itself by the omp_get_thread_num() function provided by the OpenMP runtime.
* before reading absolute values, the container must be synchronized, i.e. values from all threads are summed up and stored separately. This is a relatively slow operation and we provide ForceContainer::syncCount that you might check to find cummulative number of synchronizations and compare it against number of steps. Ideally, ForceContainer is only synchronized once at each step.
* the container is resized whenever an element outside the current range is read/written to (the read returns zero in that case); this avoids the necessity of tracking number of bodies, but also is potential danger (such as ``scene->forces.getForce(1000000000)``, which will probably exhaust your RAM). Unlike c++, Python does check given id against number of bodies.


.. code-block:: c++

	// resetting forces (inside ForceResetter)
	scene->forces.reset()

	// in a parallel section
	scene->forces.addForce(id,force); // add force

	// container is not synced after we wrote to it, sync before reading
	scene->forces.sync();
	const Vector3r& f=scene->forces.getForce(id);

Synchronization is handled automatically if values are read from python:

.. ipython::

	Yade [0]: O.bodies.append(Body())

	Yade [0]: O.forces.addF(0,Vector3(1,2,3))

	Yade [0]: O.forces.f(0)

	Yade [0]: O.forces.f(100)



.. _interaction-flow:

Handling interactions
----------------------

Creating and removing interactions is a rather delicate topic and number of components must cooperate so that the whole behaves as expected.

Terminologically, we distinguish

potential interactions,
	having neither :yref:`geometry<Interaction::geom>` nor :yref:`physics<Interaction::phys>`. :yref:`Interaction.isReal` can be used to query the status (``Interaction::isReal()`` in c++).

real interactions,
	having both :yref:`geometry<Interaction::geom>` and :yref:`physics<Interaction::phys>`. Below, we shall discuss the possibility of interactions that only have geometry but no physics.

During each step in the simulation, the following operations are performed on interactions in a typical simulation:

#. Collider creates potential interactions based on spatial proximity. Not all pairs of bodies are susceptible of entering interaction; the decision is done in Collider::mayCollide:

	* clumps may not enter interactions (only their members can)
	* clump members may not interact if they belong to the same clump
	* bitwise AND on both bodies' :yref:`masks<Body::groupMask>` must be non-zero (i.e. there must be at least one bit set in common)

#. Collider erases interactions that were requested for being erased (see below).

#. :yref:`InteractionLoop` (via :yref:`IGeomDispatcher`) calls appropriate :yref:`IGeomFunctor` based on :yref:`Shape` combination of both bodies, if such functor exists. For real interactions, the functor updates associated :yref:`IGeom`. For potential interactions, the functor returns

	``false``
		if there is no geometrical overlap, and the interaction will stillremain potential-only
	``true``
		if there is geometrical overlap; the functor will have created an :yref:`IGeom` in such case.
	
	.. note ::
		For *real* interactions, the functor *must* return ``true``, even if there is no more spatial overlap between bodies. If you wish to delete an interaction without geometrical overlap, you have to do this in the :yref:`LawFunctor`.
		
		This behavior is deliberate, since different :yref:`laws<LawFunctor>` have different requirements, though ideally using relatively small number of generally useful :yref:`geometry functors<IGeomFunctor>`.

	.. note::
		If there is no functor suitable to handle given combination of :yref:`shapes<Shape>`, the interaction will be left in potential state, without raising any error.

#. For real interactions (already existing or just created in last step), :yref:`InteractionLoop` (via :yref:`IPhysDispatcher`) calls appropriate :yref:`IPhysFunctor` based on :yref:`Material` combination of both bodies. The functor *must* update (or create, if it doesn't exist yet) associated :yref:`IPhys` instance. It is an error if no suitable functor is found, and an exception will be thrown.

#. For real interactions, :yref:`InteractionLoop` (via :yref:`LawDispatcher`) calls appropriate :yref:`LawFunctor` based on combination of :yref:`IGeom` and :yref:`IPhys` of the interaction. Again, it is an error if no functor capable of handling it is found.

#. :yref:`LawDispatcher` can decide that an interaction should be removed (such as if bodies get too far apart for non-cohesive laws; or in case of complete damage for damage models). This is done by calling

	.. code-block:: c++

		InteractionContainer::requestErase(id1,id2)

	Such interaction will not be deleted immediately, but will be reset to potential state. At next step, the collider will call ``InteractionContainer::erasePending``, which will only completely erase interactions the collider indicates; the rest will be kept in potential state.
	

Creating interactions explicitly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Interactions may still be created explicitly with :yref:`yade.utils.createInteraction`, without any spatial requirements. This function searches current engines for dispatchers and uses them. :yref:`IGeomFunctor` is called with the ``force`` parameter, obliging it to return ``true`` even if there is no spatial overlap.


Associating Material and State types
------------------------------------

Some models keep extra :yref:`state<State>` information in the :yref:`Body.state` object, therefore requiring strict association of a :yref:`Material` with a certain :yref:`State` (for instance, :yref:`CpmMat` is associated to :yref:`CpmState` and this combination is supposed by engines such as :yref:`CpmStateUpdater`).

If a :yref:`Material` has such a requirement, it must override 2 virtual methods:

#. :yref:`Material.newAssocState`, which returns a new :yref:`State` object of the corresponding type. The default implementation returns :yref:`State` itself. 
#. :yref:`Material.stateTypeOk`, which checks whether a given :yref:`State` object is of the corresponding type (this check is run at the beginning of the simulation for all particles).

In c++, the code looks like this (for :yref:`CpmMat`):

.. code-block:: c++

	class CpmMat: public FrictMat {
	   public:
	      virtual shared_ptr<State> newAssocState() const { return shared_ptr<State>(new CpmState); }
	      virtual bool stateTypeOk(State* s) const { return (bool)dynamic_cast<CpmState*>(s); }
	   /* ... */
	};

This allows one to construct :yref:`Body` objects from functions such as :yref:`yade.utils.sphere` only by knowing the requires :yref:`Material` type, enforcing the expectation of the model implementor.


Runtime structure
==================

Startup sequence
-----------------
Yade's main program is python script in :ysrc:`core/main/main.py.in`; the build system replaces a few ``\${variables}`` in that file before copying it to its install location. It does the following:

#. Process command-line options, set environment variables based on those options.
#. Import main yade module (``import yade``), residing in :ysrc:`py/__init__.py.in`. This module locates plugins (recursive search for files ``lib*.so`` in the ``lib`` installation directory). :ysrc:`yade.boot<core/main/pyboot.cpp>` module is used to setup logging, temporary directory, … and, most importantly, loads plugins.
#. Manage further actions, such as running scripts given at command line, opening :yref:`yade.qt.Controller` (if desired), launching the ``ipython`` prompt.


Singletons 
-----------
There are several "global variables" that are always accessible from c++ code; properly speaking, they are `Singletons <http://en.wikipedia.org/wiki/Singleton_pattern>`_, classes of which exactly one instance always exists. The interest is to have some general functionality acessible from anywhere in the code, without the necessity of passing pointers to such objects everywhere. The instance is created at startup and can be always retrieved (as non-const reference) using the ``instance()`` static method (e.g. ``Omega::instance().getScene()``).

There are 3 singletons:

``SerializableSingleton``
	Handles serialization/deserialization; it is not used anywhere except for the serialization code proper.
``ClassFactory``
	Registers classes from plugins and able to factor instance of a class given its name as string (the class must derive from ``Factorable``). Not exposed to python.
``Omega``
	Access to simulation(s); deserves separate section due to its importance.

Omega
^^^^^^
The :yref:`Omega` class handles all simulation-related functionality: loading/saving, running, pausing. 

In python, the wrapper class to the singleton is instantiated [#oinst]_ as global variable ``O``. For convenience, :yref:`Omega` is used as proxy for scene's attribute: although multiple ``Scene`` objects may be instantiated in c++, it is always the current scene that :yref:`Omega` represents.

The correspondence of data is literal: :yref:`Omega.materials` corresponds to ``Scene::materials`` of the current scene; likewise for :yref:`materials<Omega.materials>`, :yref:`bodies<Omega.bodies>`, :yref:`interactions<Omega.interactions>`, :yref:`tags<Omega.tags>`, :yref:`cell<Omega.cell>`, :yref:`engines<Omega.engines>`, :yref:`initializers<Omega.initializers>`, :yref:`miscParams<Omega.miscParams>`.

To give an overview of (some) variables:

======================== ===================================
Python                   c++
======================== ===================================
:yref:`Omega.iter`       ``Scene::iter``
:yref:`Omega.dt`         ``Scene::dt``
:yref:`Omega.time`       ``Scene::time``
:yref:`Omega.realtime`   ``Omega::getRealTime()``
:yref:`Omega.stopAtIter` ``Scene::stopAtIter``
======================== ===================================

``Omega`` in c++ contains pointer to the current scene (``Omega::scene``, retrieved by ``Omega::instance().getScene()``). Using :yref:`Omega.switchScene`, it is possible to swap this pointer with ``Omega::sceneAnother``, a completely independent simulation. This can be useful for example (and this motivated this functionality) if while constructing simulation, another simulation has to be run to dynamically generate (i.e. by running simulation) packing of spheres. 

.. rubric:: Footnotes

.. [#oinst] It is understood that instantiating ``Omega()`` in python only instantiates the wrapper class, not the singleton itself.

Engine loop
------------

Running simulation consists in looping over :yref:`Engines<Engine>` and calling them in sequence. This loop is defined in ``Scene::moveToNextTimeStep`` function in :ysrc:`core/Scene.cpp<core/Scene.cpp#L73>`. Before the loop starts, :yref:`O.initializers<Omega.initializers>` are called; they are only run once. The engine loop does the following in each iteration over :yref:`O.engines<Omega.engines>`:

#. set ``Engine::scene`` pointer to point to the current ``Scene``.
#. Call ``Engine::isActivated()``; if it returns ``false``, the engine is skipped.
#. Call ``Engine::action()``
#. If :yref:`O.timingEnabled<Omega.timingEnabled>`, increment :yref:`Engine::execTime` by the difference from the last time reading (either after the previous engine was run, or immediately before the loop started, if this engine comes first). Increment :yref:`Engine::execCount` by 1.

After engines are processed, :yref:`virtual time<Omega.time>` is incremented by :yref:`timestep<Omega.dt>` and :yref:`iteration number<Omega.iter>` is incremented by 1.


Background execution
^^^^^^^^^^^^^^^^^^^^^

The engine loop is (normally) executed in background thread (handled by :ysrc:`SimulationFlow<core/SimulationFlow.hpp>` class), leaving foreground thread free to manage user interaction or running python script. The background thread is managed by :yref:`O.run()<Omega.run>` and :yref:`O.pause()<Omega.pause>` commands. Foreground thread can be blocked until the loop finishes using :yref:`O.wait()<Omega.wait>`.

Single iteration can be run without spawning additional thread using :yref:`O.step()<Omega.step>`.



Python framework
=================


Wrapping c++ classes
---------------------

Each class deriving from :yref:`Serializable` is automatically exposed to python, with access to its (registered) attributes. This is achieved via :ref:`YADE_CLASS_BASE_DOC`. All classes registered in class factory are default-constructed in ``Omega::buildDynlibDatabase``. Then, each serializable class calls ``Serializable::pyRegisterClass`` virtual method, which injects the class wrapper into (initially empty) ``yade.wrapper`` module. ``pyRegisterClass`` is defined by ``YADE_CLASS_BASE_DOC`` and knows about class, base class, docstring, attributes, which subsequently all appear in boost::python class definition.

Wrapped classes define special constructor taking keyword arguments corresponding to class attributes; therefore, it is the same to write:

.. ipython::

	Yade [1]: f1=ForceEngine()

	Yade [2]: f1.ids=[0,4,5]

	Yade [2]: f1.force=Vector3(0,-1,-2)

and 

.. ipython::

	Yade [1]: f2=ForceEngine(ids=[0,4,5],force=Vector3(0,-1,-2))

	Yade [2]: print f1.dict()

	Yade [3]: print f2.dict()

Wrapped classes also inherit from :yref:`Serializable` several special virtual methods: :yref:`dict()<Serializable::dict>` returning all registered class attributes as dictionary (shown above), :yref:`clone()<Serializable::clone>` returning copy of instance (by copying attribute values), :yref:`updateAttrs()<Serializable::updateAttrs>` and :yref:`updateExistingAttrs()<Serializable::updateExistingAttrs>` assigning attributes from given dictionary (the former thrown for unknown attribute, the latter doesn't).

Read-only property ``name`` wraps c++ method ``getClassName()`` returning class name as string. (Since c++ class and the wrapper class always have the same name, getting python type using ``__class__`` and its property ``__name__`` will give the same value).

.. ipython::

	Yade [1]: s=Sphere()

	Yade [2]: s.__class__.__name__


Subclassing c++ types in python
--------------------------------

In some (rare) cases, it can be useful to derive new class from wrapped c++ type in pure python. This is done in the :ref:`yade.pack` module: :yref:`Predicate<yade._packPredicates.Predicate>` is c++ base class; from this class, several c++ classes are derived (such as :yref:`inGtsSurface<yade._packPredicates.Predicate>`), but also python classes (such as the trivial :yref:`inSpace<yade.pack.inSpace>` predicate). ``inSpace`` derives from python class ``Predicate``; it is, however, not direct wrapper of the c++ ``Predicate`` class, since virtual methods would not work.

``boost::python`` provides special ``boost::python::wrapper`` template for such cases, where each overridable virtual method has to be declared explicitly, requesting python override of that method, if present. See `Overridable virtual functions <http://wiki.python.org/moin/boost.python/OverridableVirtualFunctions>`_ for more details.


Reference counting
------------------
Python internally uses `reference counting <http://en.wikipedia.org/wiki/Reference_counting>`_ on all its objects, which is not visible to casual user. It has to be handled explicitly if using pure `Python/C API <http://docs.python.org/c-api/index.html>`_ with ``Py_INCREF`` and similar functions.

``boost::python`` used in Yade fortunately handles reference counting internally. Additionally, it `automatically integrates <http://wiki.python.org/moin/boost.python/PointersAndSmartPointers>`_ reference counting for ``shared_ptr`` and python objects, if class ``A`` is wrapped as ``boost::python::class_<A,shared_ptr<A>>``. Since *all* Yade classes wrapped using :ref:`YADE_CLASS_BASE_DOC` are wrapped in this way, returning ``shared_ptr<…>`` objects from is the preferred way of passing objects from c++ to python.

Returning ``shared_ptr`` is much more efficient, since only one pointer is returned and reference count internally incremented. Modifying the object from python will modify the (same) object in c++ and vice versa. It also makes sure that the c++ object will not be deleted as long as it is used somewhere in python, preventing (important) source of crashes.

.. _customconverters:

Custom converters
-----------------

When an object is passed from c++ to python or vice versa, then either

#. the type is basic type which is transparently passed between c++ and python (int, bool, std::string etc)
#. the type is wrapped by boost::python (such as Yade classes, ``Vector3`` and so on), in which case wrapped object is returned; [#wrap]_

Other classes, including template containers such as ``std::vector`` must have their custom converters written separately. Some of them are provided in :ysrc:`py/wrapper/customConverters.cpp`, notably converters between python (homogeneous, i.e. with all elements of the same type) sequences and c++ ``std::vector`` of corresponding type; look in that source file to add your own converter or for inspiration.

When an object is crossing c++/python boundary, boost::python's global "converters registry" is searched for class that can perform conversion between corresponding c++ and python types. The "converters registry" is common for the whole program instance: there is no need to register converters in each script (by importing ``_customConverters``, for instance), as that is done by yade at startup already.

.. note::
	Custom converters only work for value that are passed by value to python (not "by reference"): some attributes defined using :ref:`YADE_CLASS_BASE_DOC` are passed by value, but if you define your own, make sure that you read and understand `Why is my automatic to-python conversion not being found? <http://www.boost.org/doc/libs/1_42_0/libs/python/doc/v2/faq.html#topythonconversionfailed>`_.

	In short, the default for ``def_readwrite`` and ``def_readonly`` is to return references to underlying c++ objects, which avoids performing conversion on them. For that reason, return value policy must be set to ``return_by_value`` explicitly, using slighly more complicated ``add_property`` syntax, as explained at the page referenced.

.. [#wrap]
	Wrapped classes are automatically registered when the class wrapper is created. If wrapped class derives from another wrapped class (and if this dependency is declared with the ``boost::python::bases`` template, which Yade's classes do automatically), parent class must be registered before derived class, however. (This is handled via loop in ``Omega::buildDynlibDatabase``, which reiterates over classes, skipping failures, until they all successfully register)
	Math classes (Vector3, Matrix3, Quaternion) are wrapped by hand, to be found in :ysrc:`py/mathWrap/miniEigen.cpp`; this module is imported at startup.


Maintaining compatibility
==========================

In Yade development, we identified compatibility to be very strong desire of users. Compatibility concerns python scripts, *not* simulations saved in XML or old c++ code.

Renaming class
---------------
Script :ysrc:`scripts/rename-class.py` should be used to rename class in ``c++`` code. It takes 2 parameters (old name and new name) and must be run from top-level source directory::

	\$ scripts/rename-class.py OldClassName NewClassName
	Replaced 4 occurences, moved 0 files and 0 directories
	Update python scripts (if wanted) by running: perl -pi -e 's/\bOldClassName\b/NewClassName/g' `ls **/*.py |grep -v py/system.py`

This has the following effects:

#. If file or directory has basename ``OldClassName`` (plus extension), it will be renamed using ``bzr``.
#. All occurences of whole word ``OldClassName`` will be replaced by ``NewClassName`` in c++ sources.
#. An extry is added to :ysrc:`py/system.py`, which contains map of deprecated class names. At yade startup, proxy class with ``OldClassName`` will be created, which issues a ``DeprecationWarning`` when being instantiated, informing you of the new name you should use; it creates an instance of ``NewClassName``, hence not disruting your script's functioning::

	Yade [3]: SimpleViscoelasticMat()
	/usr/local/lib/yade-trunk/py/yade/__init__.py:1: DeprecationWarning: Class `SimpleViscoelasticMat' was renamed to (or replaced by) `ViscElMat', update your code! (you can run 'yade --update script.py' to do that automatically)
	->  [3]: <ViscElMat instance at 0x2d06770>

As you have just been informed, you can run ``yade --update`` to all old names with their new names in scripts you provide::

	\$ yade-trunk --update script1.py some/where/script2.py 

This gives you enough freedom to make your class name descriptive and intuitive.

Renaming class attribute
------------------------

Renaming class attribute is handled from c++ code. You have the choice of merely warning at accessing old attribute (giving the new name), or of throwing exception in addition, both with provided explanation. See ``deprec`` parameter to :ref:`YADE_CLASS_BASE_DOC` for details.

Debian packaging instructions
==============================
In order to make parallel installation of several Yade version possible, we adopted similar strategy as e.g. ``gcc`` packagers in Debian did:

#. Real Yade packages are named ``yade-0.30`` (for stable versions) or ``yade-bzr2341`` (for snapshots).
#. They provide ``yade`` or ``yade-snapshot`` virtual packages respectively.
#. Each source package creates several installable packages (using ``bzr2341`` as example version):

   #. ``yade-bzr2341`` with the optimized binaries; the executable binary is ``yade-bzr2341`` (``yade-bzr2341-multi``, …)
   #. ``yade-bzr2341-dbg`` with debug binaries (debugging symbols, non-optimized, and with crash handlers); the executable binary is ``yade-bzr2341-dbg``
   #. ``yade-bzr2341-doc`` with sample scripts and some documentation (see `bug #398176 <https://bugs.launchpad.net/yade/+bug/398176>`_ however)
   #. (future?) ``yade-bzr2341-reference`` with reference documentation (see `bug #401004 <https://bugs.launchpad.net/yade/+bug/401004>`_)
#. Using `Debian alternatives <http://www.debian-administration.org/articles/91>`_, the highest installed package provides additionally commands without the version specification like ``yade``, ``yade-multi``, … as aliases to that version's binaries. (``yade-dbg``, … for the debuggin packages). The exact rule is:

   #. Stable releases have always higher priority than snapshots
   #. Higher versions/revisions have higher pripority than lower versions/revisions.

Prepare source package
----------------------

Debian packaging files are located in :ysrc:`debian/` directory. They contain build recipe :ysrc:`debian/rules`, dependecy and package declarations :ysrc:`debian/control` and maintainer scripts. Some of those files are only provided as templates, where some variables (such as version number) are replaced by special script.

The script :ysrc:`scripts/debian-prep` processes templates in :ysrc:`debian/` and creates files which can be used by debian packaging system. Before running this script:

#. If you are releasing stable version, make sure there is file named ``RELEASE`` containing single line with version number (such as ``0.30``). This will make :ysrc:`scripts/debian-prep` create release packages. In absence of this file, snapshots packaging will be created instead. Release or revision number (as detected by running ``bzr revno`` in the source tree) is stored in ``VERSION`` file, where it is picked up during package build and embedded in the binary.
#. Find out for which debian/ubuntu series your package will be built. This is the name that will appear on the top of (newly created) ``debian/changelog`` file. This name will be usually ``unstable``, ``testing`` or ``stable`` for debian and ``karmic``, ``lucid`` etc for ubuntu. WHen package is uploaded to Launchpad's build service, the package will be built for this specified release.

Then run the script from the top-level directory, giving series name as its first (only) argument::

	\$ scripts/debian-prep lucid

After this, signed debian source package can be created::

	\$ debuild -S -sa -k62A21250 -I -Iattic

(``-k`` gives GPG key identifier, ``-I`` skips ``.bzr`` and similar directories, ``-Iattic`` will skip the useless ``attic`` directory).

Create binary package
---------------------
Local in-tree build
	Once files in ``debian/`` are prepared, packages can be build by issuing::
	\$ fakeroot debian/rules binary
Clean system build
	Using ``pbuilder`` system, package can be built in a chroot containing clean debian/ubuntu system, as if freshly installed. Package dependencies are automatically installed and package build attempted. This is a good way of testing packaging before having the package built remotely at Launchpad. Details are provided at `wiki page <https://www.yade-dem.org/wiki/DebianPackages>`_.
Launchpad build service
	Launchpad provides service to compile package for different ubuntu releases (series), for all supported architectures, and host archive of those packages for download via APT. Having appropriate permissions at Launchpad (verified GPG key), source package can be uploaded to yade's archive by::

		\$ dput ppa:yade-users/ppa ../yade-bzr2341_1_source.changes

	After several hours (depending on load of Launchpad build farm), new binary packages will be published at https://launchpad.net/~yade-users/+archive/ppa.
	
	This process is well documented at https://help.launchpad.net/Packaging/PPA.
