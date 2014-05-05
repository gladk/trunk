###############
Installation
###############

Yade can be installed from packages (precompiled binaries) or source code. 
The choice depends on what you need: if you don't plan to modify Yade itself, 
package installation is easier. In the contrary case, you must download and 
install source code.

Packages
----------

Pre-built packages are provided for all currently supported Debian and Ubuntu 
versions of distributions and available on `yade-dem.org <http://yade-dem.org/packages/>`_  
server. 

These are ``daily`` versions of packages and are updating regularly and include 
all the newly added features.

To install daily-version one needs to add this repository to your 
/etc/apt/sources.list, add a PGP-key AA915EEB as a trusted and install yadedaily ::

	sudo bash -c 'echo "deb http://www.yade-dem.org/packages/ precise/" >> /etc/apt/sources.list'
	wget -O - http://www.yade-dem.org/packages/yadedev_pub.gpg | sudo apt-key add -
	sudo apt-get update
	sudo apt-get install yadedaily

If you have another distribution, not Ubuntu Precise (Version 12.04), be sure to use the
correct name in the first line (for instance, jessie, trusty or wheezy).

After that you can normally start Yade using "yadedaily" or "yadedaily-batch" command.
yadedaily on older distributions can have some disabled features due to older library
versions, shipped with particular distribution. 

Git-repository for packaging stuff is available on `GitHub <https://github.com/yade/yadedaily/>`_. 
Each branch corresponds to one distribution e.g. precise, jessie etc.
The scripts for building all of this stuff is `here <https://github.com/yade/trunk/tree/master/scripts/ppa>`_. 
It uses pbuilder to build packages, so all packages are building in a clean environment.

If you do not need yadedaily-package any more, just remove the
corresponding line in /etc/apt/sources.list and the package itself::

	sudo apt-get remove yadedaily

To remove our key from keyring, execute the following command::

	sudo apt-key remove AA915EEB

Since 2011 all Ubuntu versions (starting from 11.10, Oneiric) and Debian (starting from Wheezy) 
are having already Yade in their main repositories. There are only stable releases are placed.
To install the program, run the following::

	sudo apt-get install yade

To check, what version of Yade is in specific distribution, visit the links
for `Ubuntu <https://launchpad.net/ubuntu/+source/yade>`_ and 
`Debian <http://packages.qa.debian.org/y/yade.html>`_. 
`Debian-Backports <http://backports.debian.org/Instructions>`_ 
repository is updating regularly to bring the newest Yade to a users of stable 
Debians.

Daily and stable Yade versions can coexist without any conflicts.

Source code
------------

Installation from source code is reasonable, when you want to add or 
modify constitutive laws, engines or functions... Installing the latest 
trunk version allows one to use newly added features, which are not yet 
available in packaged versions. 

Download
^^^^^^^^^^

If you want to install from source, you can install either a release 
(numbered version, which is frozen) or the current developement version 
(updated by the developers frequently). You should download the development 
version (called ``trunk``) if you want to modify the source code, as you 
might encounter problems that will be fixed by the developers. Release 
version will not be modified (except for updates due to critical and 
easy-to-fix bugs), but they are in a more stabilized state that trunk 
generally.

#. Releases can be downloaded from the `download page <https://launchpad.net/yade/+download>`_, as compressed archive. Uncompressing the archive gives you a directory with the sources.
#. developement version (trunk) can be obtained from the `code repository <https://github.com/yade/>`_ at github.

We use `GIT <http://git-scm.com/>`_ (the ``git`` command) for code 
management (install the ``git`` package in your distribution and create a GitHub account)::

		git clone git@github.com:yade/trunk.git

will download the whole code repository of ``trunk``. Check out `Yade on github
<https://www.yade-dem.org/wiki/Yade_on_github>`_ wiki page for more.

Alternatively, a read-only checkout is possible via https without a GitHub account (easier if you don't want to modify the main Yade branch)::

		git clone https://github.com/yade/trunk.git
   
For those behind firewall, you can download `any revision  <https://www.yade-dem.org/source/>`_ of the repository as compressed archive.

Release and trunk sources are compiled in the same way.

Prerequisites
^^^^^^^^^^^^^

Yade relies on a number of external software to run; they are checked before the compilation starts.
Some of them are only optional. The last ones are only relevant for using the fluid coupling module (:yref:`FlowEngine`).

* `cmake <http://www.cmake.org/>`_ build system
* `gcc <http://www.gcc.gnu.org>`_ compiler (g++); other compilers will not work; you need g++>=4.2 for openMP support
* `boost <http://www.boost.org/>`_ 1.35 or later
* `qt4 <http://www.qt.nokia.com>`_ library
* `freeglut3 <http://freeglut.sourceforge.net>`_
* `libQGLViewer <http://www.libqglviewer.com>`_
* `python <http://www.python.org>`_, `numpy <http://numpy.scipy.org>`_, `ipython <http://ipython.scipy.org>`_
* `matplotlib <http://matplotlib.sf.net>`_
* `eigen3 <http://eigen.tuxfamily.org>`_ algebra library
* `gdb <http://www.gnu.org/software/gdb>`_ debugger
* `sqlite3 <http://www.sqlite.org>`_ database engine
* `Loki <http://loki-lib.sf.net>`_ library
* `VTK <http://www.vtk.org/>`_ library (optional but recommended)
* `CGAL <http://www.cgal.org/>`_ library (optional)
* `SuiteSparse <http://www.cise.ufl.edu/research/sparse/SuiteSparse/>`_ sparse algebra library (fluid coupling, optional, requires eigen>=3.1)
* `OpenBLAS <http://www.openblas.net/>`_ optimized and parallelized alternative to the standard blas+lapack (fluid coupling, optional)
* `Metis <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview/>`_ matrix preconditioning (fluid coupling, optional)

Most of the list above is very likely already packaged for your distribution. In case you are confronted
with some errors concerning not available packages (e.g. Package libmetis-dev is not available) it may be necessary 
to add yade external ppa from https://launchpad.net/~yade-users/+archive/external::

	sudo add-apt-repository ppa:yade-users/external 
	sudo apt-get update 

The following commands have to be executed in command line of corresponding 
distributions. Just copy&paste to the terminal. To perform commands you 
should have root privileges

.. warning:: If you have Ubuntu 12.10 or older, you need to install libqglviewer-qt4-dev
 package instead of libqglviewer-dev.

 
* **Ubuntu**, **Debian** and their derivatives::

		sudo apt-get install cmake git freeglut3-dev libloki-dev \
		libboost-all-dev fakeroot dpkg-dev build-essential g++ \
		python-dev ipython python-matplotlib libsqlite3-dev python-numpy python-tk gnuplot \
		libgts-dev python-pygraphviz libvtk5-dev python-scientific libeigen3-dev \
		python-xlib python-qt4 pyqt4-dev-tools gtk2-engines-pixbuf python-argparse \
		libqglviewer-dev python-imaging libjs-jquery python-sphinx python-git python-bibtex \
		libxmu-dev libxi-dev libcgal-dev help2man libbz2-dev zlib1g-dev
		

Some of packages (for example, cmake, eigen3) are mandatory, some of them
are optional. Watch for notes and warnings/errors, which are shown
by cmake during configuration step. If the missing package is optional,
some of Yade features will be disabled (see the messages at the end of configuration).
		
Additional packages, which can become mandatory later::

		sudo apt-get install python-gts python-minieigen
		
For effective usage of direct solvers in the PFV-type fluid coupling, the following libraries are recommended, together with eigen>=3.1: blas, lapack, suitesparse, and metis.
All four of them are available in many different versions. Different combinations are possible and not all of them will work. The following was found to be effective on recent deb-based systems. On ubuntu 12.04, better compile openblas with USE_OPENMP=1, else yade will run on a single core::

		sudo apt-get install libopenblas-dev libsuitesparse-metis-dev

Some packages listed here are relatively new and they can be absent
in your distribution (for example, libmetis-dev or python-gts). They can be 
installed from our `external PPA <https://launchpad.net/~yade-users/+archive/external/>`_
or just ignored. In this case some features can be disabled.

If you are using other distribution, than Debian or its derivatives, you should
install the softwares listed above. Their names can differ from the 
names of Debian-packages.


Compilation
^^^^^^^^^^^

You should create a separate build-place-folder, where Yade will be configured 
and where the source code will be compiled. Here is an example for a folderstructure:

    myYade/           ## base directory
            trunk/      ## folder for sourcecode in which you use github
            build/      ## folder in which sources will be compiled; build-directory; use cmake here
            install/    ## installfolder

Then inside this build-directory you should start cmake to configure the compilation process::

	cmake -DINSTALL_PREFIX=/path/to/installfolder /path/to/sources

For the folder structure given above call the following command in folder "build":

    cmake -DINSTALL_PREFIX=../install ../trunk

Additional options can be configured in the same line with the following 
syntax::

	cmake -DOPTION1=VALUE1 -DOPTION2=VALUE2
	
The following options are available:
	
	* INSTALL_PREFIX: path where Yade should be installed (/usr/local by default)
	* LIBRARY_OUTPUT_PATH: path to install libraries (lib by default)
	* DEBUG: compile in debug-mode (OFF by default)
	* CMAKE_VERBOSE_MAKEFILE: output additional information during compiling (OFF by default)
	* SUFFIX: suffix, added after binary-names (version number by default)
	* NOSUFFIX: do not add a suffix after binary-name (OFF by default)
	* YADE_VERSION: explicitely set version number (is defined from git-directory by default)
	* ENABLE_GUI: enable GUI option (ON by default)
	* ENABLE_CGAL: enable CGAL option (ON by default)
	* ENABLE_VTK: enable VTK-export option (ON by default)
	* ENABLE_OPENMP: enable OpenMP-parallelizing option (ON by default)
	* ENABLE_GTS: enable GTS-option (ON by default)
	* ENABLE_GL2PS: enable GL2PS-option (ON by default)
	* ENABLE_LINSOLV: enable LINSOLV-option (ON by default)
	* ENABLE_PFVFLOW: enable PFVFLOW-option, FlowEngine (ON by default)
	* runtimePREFIX: used for packaging, when install directory is not the same is runtime directory (/usr/local by default)
	* CHUNKSIZE: used, if you want several sources to be compiled at once. Increases compilation speed and RAM-consumption during it (1 by default).

For using an extended parameters of cmake, please, follow the corresponding
documentation on cmake-webpage. 

If the compilation is finished without errors, you will see all enabled 
and disabled options. Then start the standard the compilation process::

	make

Installing performs with the following command::

	make install

The "install" command will in fact also recompile if source files have been modified. 
Hence there is no absolute need to type the two commands separately.

The compilation process can take a long time, be patient. An additional
parameter on many cores systems ``-j`` can be added to decrease compilation time
and split the compilation on many cores. For example, on 4-core machines
it would be reasonable to set the parameter ``-j4``. Note, the Yade requires
approximately 2GB/core for compilation, otherwise the swap-file will be used
and a compilation time dramatically increases. After compilation finished successfully
the new built can be started by navigating to /path/to/installfolder/bin and calling yade via (based on version yade-2014-02-20.git-a7048f4)::
    
    cd /path/to/installfolder/bin 
    ./yade-2014-02-20.git-a7048f4

For building the documentation you should at first execute the command "make install"
and then "make doc" to build it. The generated files will be stored in your current
build directory/doc/sphinx/_build.

"make manpage" command generates and moves manpages in a standard place.
"make check" command executes standard test to check the functionality of compiled program.

Yade can be compiled not only by GCC-compiler, but also by `CLANG <http://clang.llvm.org/>`_ 
front-end for the LLVM compiler. For that you set the environment variables CC and CXX 
upon detecting the C and C++ compiler to use::

	export CC=/usr/bin/clang
	export CXX=/usr/bin/clang++
	cmake -DOPTION1=VALUE1 -DOPTION2=VALUE2

Clang does not support OpenMP-parallelizing for the moment, that is why the 
feature will be disabled.
