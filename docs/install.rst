Installation
============

This section provides a quick guide for compiling and running ASYNCH on a Unix based system, including Iowa's local HPC clusters. Windows users can use the Visual Studio project in the ``ide\msvc`` folder.

Requirements
------------

In addition to the source code, several programs and libraries are needed:

-  A C compiler
-  GNU Make
-  An MPI implementation
-  HDF5, libpq libraries

A brief description of each is provided below.

C Compiler
~~~~~~~~~~

ASYNCH is written in ANSI C89 standard. ASYNCH has been successfully compiled and tested with the `GNU C compiler <https://gcc.gnu.org/>`__ gcc (version 4.6 and later), Microsoft Visual Studio 2015 and Intel C++ compiler. Using Intel's compiler is recommanded if available.

GNU Make
~~~~~~~~

`GNU Make <http://www.gnu.org/software/make/>`__ is a utility for directing the generation of binaries from source code. ``make`` is available in most Linux distro.

MPI Implementation
~~~~~~~~~~~~~~~~~~

The Message Passing Interface (MPI) is a standard for transferring data on parallel computers ASYNCH uses MPI for communication between processes to perform calculations in parallel ASYNCH has been successfully used and tested with ``OpenMPI``.

From the `OpenMPI webpage <http://www.open-mpi.org/>`__ :

  The Open MPI Project is an open source Message Passing Interface implementation that is developed and maintained by a consortium of academic, research, and industry partners. OpenMPI is available on Iowa's HPC clusters In theory, ASYNCH should work properly with any other implementation adhering to the MPI standard.

If installing an MPI implementation on a machine for ASYNCH, be sure to install the development packages of the MPI implementation. In Linux repositories, these packages are usually denoted with a ``-dev`` or similar in the package name. If in doubt, try typing "mpirun" and "mpicc" in a terminal. Both of these should be present to run ASYNCH. If mpirun is not present, you have not installed the MPI binaries (meaning you probably haven't tried installing MPI at all). If mpicc is not present, then you are missing the development package.

HDF5 Library
~~~~~~~~~~~~

``HDF5`` is a format and library widely used for storing binary scientific data in an efficient and portable way.

From the `HDF5 webpage <https://support.hdfgroup.org/HDF5/>`__ :

  HDF5 is a data model, library, and file format for storing and managing data. It supports an unlimited variety of datatypes, and is designed for filexible and efficient I/O and for high volume and complex data. HDF5 is portable and is extensible, allowing applications to evolve in their use of HDF5. The HDF5 Technology suite includes tools and applications for managing, manipulating, viewing, and analyzing data in the HDF5 format.

For ASYNCH build, the 1.8.x branch is used.

Optional Libraries
------------------

PostgreSQL Library
~~~~~~~~~~~~~~~~~~

``libpq`` is a library for communicating with a PostgreSQL database, created by the makers of PostgreSQL.

From the `libpq webpage <http://www.postgresql.org/does/9.1/statie/libpq.html>`__ :

  libpq is the C application programmer's interface to PostgreSQL. libpq is a set of library functions that allow client programs to pass queries to the PostgreSQL back-end server and to receive the results of these queries.

libpq can be ommitted while running ``./configure`` with ``--without-postgresql``. If you want to use libpq on a machine for ASYNCH, be sure to install the PostgreSQL development packages. In Linux repositories, these packages are usually denoted by a "-dev" or similar in the package name.

METIS Library
~~~~~~~~~~~~~

``METIS`` is a set of serial programs for partitioning graphs, partitioning finite element meshes, and producing fill reducing orderings for sparse matrices.

From the `METIS webpage <http://glaros.dtc.umn.edu/gkhome/metis/metis/overview>`__ :

  Experiments on a large number of graphs arising in various domains including finite element methods, linear programming, VLSI, and transportation show that METIS produces partitions that are consistently better than those produced by other widely used algorithms. The partitions produced by METIS are consistently 10% to 50% better than those produced by spectral partitioning algorithms.

METIS can be ommitted while running ``./configure`` with ``--without-metis``.

PETSc Library
~~~~~~~~~~~~~

``PETSc`` is a suite of data structures and routines for the scalable (parallel) solution of scientific applications modeled by partial differential equations.

From the `PETSc webpage <http://www.mcs.anl.gov/petsc/>`__ :

  PETSc is intended for use in large-scale application projects, many ongoing computational science projects are built around the PETSc libraries. PETSc is easy to use for beginners. Moreover, its careful design allows advanced users to have detailed control over the solution process.

PETSc can be ommitted while running ``./configure`` with ``--without-petsc``.

Optional Software
-----------------

Some additional software is available that may be useful, depending upon what the user wishes to do. These include

-  Git
-  dos2unix
-  NoMachine NX

Git
~~~

Git is a distributed revision control system. Although not needed for running ASYNCH, the source code repository does require Git for access. Git is in most Linux repositories. GitHub offers `guides about its usage <https://guides.github.com/activities/hello-world/>`__.

dos2unix
~~~~~~~~

This is a useful utility if editing input text files for ASYNCH from a Windows machine. Unix and Windows use a slightly different format for text documents. Although a file may look the same under both a Linux and Windows text editor, subtle diferences can still exist. In general, editing a text file from Linux on a Windows machine will convert the file to the Windows format. To change the format to Unix, use the utility dos2unix. If a file is already under Unix format, this utility will not modify the file Using a text file in Windows format with ASYNCH will result in errors. This process can be slow, depending upon the size of the text files involved As such, ASYNCH does not automatically check the format of input text files. dox2unix can be found in most Linux repositories.

FastX 2
~~~~~~~

FastX is a program for connecting to the HPC systems with a GUI desktop environment. It is similar to No Machine connections but is newer and a little more robust when using Duo 2 factor authentication. This can be useful for users wishing to access an Iowa HPC resource, though this is not the only way. Information about using and obtaining FastX for Iowa HPC resources can be found at `FastX connections <https://wiki.uiowa.edu/display/hpcdocs/FastX+connections#FastXconnections-fastx2>`__.

Source Code, Compiling, and Running ASYNCH
------------------------------------------

The ASYNCH source code is available in a repository hosted by GitHub. Downloading on of the release version the code from the repository requires the use of Git See `Git`_. The source code can also be downloaded directly from GitHub as a zip file.

If the source code is ever updated, you may want to run ``make clean`` before recompiling. This removes all binaries and object files of the old version. Once compiled, ASYNCH can be run with the command:

.. code-block:: sh

  mpirun -np <number of processes> <path>/asynch < gbl filename>

Updating the package
~~~~~~~~~~~~~~~~~~~~

This operation is only necessary if you cloned the git repository. If you are using a release source tarball, you can skip to the next step.

.. code-block:: sh

  autoreconf --install
  make dist

Installing the package
~~~~~~~~~~~~~~~~~~~~~~

These are the generic instruction for an out of source build (prefered method):

.. code-block:: sh

  mkdir build && cd build
  ../configure CFLAGS="-O3 -DNDEBUG"
  make
  make check
  make install

.. note::

   Newer version of gcc requires to add ``-Wno-format-security`` to ``CFLAGS`` argument so that the configure script should be invoked with ``../configure CFLAGS="-DNDEBUG -Wno-format-security"``.
   
Other typical ``configure`` arguments are also available, such as:

- **--prefix**: defines the folder path of the output compiled runnable.

  - Example: ``--prefix=/home/user/asynch/dist``
  
- **CC**: defines the C compiler to be used.

  - Example: ``CC=gcc`` (GNU compiler) or ``CC=icc`` (intel compiler)

Iowa HPC Clusters
-----------------

Currently, the executable used on Neon and Argon ar maintained by yours truly. Users of Iowa's HPC resources should NOT need to download and compile source code on Neon and Argon. Binaries for ASYNCH are located in ``/Dedicated/IFC/.neon/bin`` on Neon and in ``/Dedicated/IFC/.argon/bin`` on Argon. Libraries for linking `libasynch` with your own software are located in the directory ``/Dedicated/IFC/.<cluster>/lib``. As of the compilation of derived work, all required software should be available. The build system included with the source code should work without modification on these clusters.

Setting up the environment on ARGON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These clusters do use third party software through modules. The module for OpenMPI and HDF5 must be loaded once per login session to run ASYNCH. Refers to the :ref:`Getting Started` section for more information. For Argon:

.. code-block:: sh

    # User specific environment and startup programs for Argon

    export PATH=$PATH:$HOME/.local/bin:/Dedicated/IFC/.argon/bin

    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/Dedicated/IFC/.argon/lib

    # Load module OpenMPI and HDF5
    module load zlib/1.2.11_parallel_studio-2017.1
    module load hdf5/1.8.18_parallel_studio-2017.1
    module load openmpi/2.0.1_parallel_studio-2017.1

These load OpenMPI version 1.8.8 for use with the Intel compiler as well as the HDF5 1.8.18 library. Instead of loading these modules manually, the commands can be added to the end of the file ``.bash_profile`` in the user's home directory. Note that Neon and Argon each have a seperate $HOME hence ``.bash_profile`` file. In addition, if using the Python interface functions on Argon, the appropriate Python module must be loaded. This can be done with a call to:

.. code-block:: sh

  module load python27

This can also be added to the ``.bash_profile`` file to automate the loading process.

Installing the package on ARGON
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, ``git clone`` the repository or ``tar xf`` a release package.

Then run the classic GNU build tool chain:

.. code-block:: sh

  mkdir build && cd build
  ../configure --prefix=/Dedicated/IFC/.argon CFLAGS="-O3 -march=core-avx2 -DNDEBUG" PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/Dedicated/IFC/.argon/lib/pkgconfig
  make
  make check
  make install

Updating the package
--------------------

Whenever the ``autoconf`` or ``automake`` files are modified, the build system needs to be update:

.. code-block:: sh

  # Using 'make dist' with a 32 UID
  export TAR_OPTIONS="--owner=0 --group=0 --numeric-owner"

  autoreconf --install
  mkdir build && cd build
  make dist

Standard Makefile Targets
-------------------------

-  ``make all`` Build programs, libraries, documentation, etc. (Same as ``make``.)
-  ``make install`` Install what needs to be installed.
-  ``make install-strip`` Same as ``make install``, then strip debugging symbols.
-  ``make uninstall`` The opposite of ``make install``.
-  ``make clean`` Erase what has been built (the opposite of ``make all``).
-  ``make distclean`` Additionally erase anything ``./configure`` created.
-  ``make check`` Run the test suite, if any.
-  ``make installcheck`` Check the installed programs or libraries, if supported.
-  ``make dist`` Create PACKAGE-VERSION.tar.gz.
