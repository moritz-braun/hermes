SuperLU
--------

.. _SuperLU home page: http://crd.lbl.gov/~xiaoye/SuperLU/
.. _solvers repository: https://github.com/hpfem/solvers
.. _manual: https://github.com/hpfem/solvers/raw/master/manuals/SuperLU.pdf

Hermes currently supports two versions of the SuperLU library - the sequential
one and the multithreaded one. Support for the MPI version will be added in the 
future. Please visit the `SuperLU home page`_ for more information about the
library.

Linux
~~~~~

Sequential
``````````
Download the software package from the `solvers repository`_ and unpack 
it in some temporary directory::
  
  wget https://github.com/hpfem/solvers/raw/master/packages/superlu-4.0.spkg --no-check-certificate
  tar -jxvf superlu-4.0.spkg
  rm superlu-4.0.spkg
  cd superlu-4.0
  
Let's assume that you want to install the library into ``~/solvers/superlu``. 
You may choose any path you like, provided that you have write access to it
(the target directory will be created if it doesn't exist). There are now three 
typical build scenarios that you may follow:

  - Build the default, optimized version of the library::
    
      SPKG_LOCAL=~/solvers/superlu ./spkg-install
    
  - Build the debug version of the library::
    
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --debug
    
    .. __:
    
    This adds debugging symbols to the library, allowing e.g. to step through
    the library source code from within a Hermes debugging session. Note that 
    you may install the optimized and debugging versions to separate directories
    and switch between them when building your application by adjusting the 
    hermes{1|2|3}d CMake.vars file.
    
  - Build the non-optimized version of the library::

      SPKG_LOCAL=~/solvers/superlu ./spkg-install --no-optimizations

    .. __:
          
    The library has problems with running its own test-suite, although the 
    examples included with it as well as Hermes tests seem to run well.
    The multithreaded version failed its tests when any build options had been
    selected (thus the tests have been disabled), while the sequential version 
    passed them at least when compiler optimizations had been disabled 
    (i.e. this build scenario).

For advanced configuration possibilities, please read
src/README, the `manual`_ or visit the `SuperLU home page`_.

Once the library has been built and installed, you may delete the temporary 
directory with the unpacked package to save some disk space (but note that this 
will prevent you from stepping through the library code during debugging), or 
just remove the object files by running

::

  cd src
  make clean 

Now go to the directory with hermes{1|2|3}d. Create the file CMake.vars with the
following lines (or append to the existing one)::

  set(WITH_SUPERLU YES)
  set(SUPERLU_ROOT ~/solvers/superlu) #(or your own installation destination)
  set(SUPERLU_MT   NO)

Finally execute::
  
  rm CMakeCache.txt
  cmake .
  make
    

Multithreaded
`````````````

The multithreaded library can take advantage of modern multi-core
machines by computing the numerical factorization in parallel. The latest version 
specifically tuned for Hermes is available in the 
`solvers repository`_ and you may download and unpack it using the
following commands::

  wget https://github.com/hpfem/solvers/raw/master/packages/superlu_mt-2.0.spkg --no-check-certificate
  tar -jxvf superlu_mt-2.0.spkg
  rm superlu_mt-2.0.spkg
  cd superlu_mt-2.0

(the list of changes made to the original distribution from `SuperLU home page`_
in order to make it compatible with Hermes may be found in src/MODIFICATIONS).

There are two multithreading models supported by SuperLU on Linux

  - `POSIX threads <https://computing.llnl.gov/tutorials/pthreads/>`__ (or Pthreads) - standard model 
    available in most Linux distributions.
    
  - `OpenMP <http://openmp.org/wp/>`__ - should be included in recent GNU compilers (since GCC 4.3.2);
    if you have an older version, you may install it via the libgomp package, e.g.
    in Ubuntu::
    
      sudo apt-get install libgomp1      

Assuming the intended installation directory is ``~/solvers/superlu``, you may
build a particular version of the multithreaded library by issuing one 
of the following available build commands:

  - Build the default, optimized version of the library::
    
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --with-openmp
      
    or
      
    ::
      
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --with-pthreads
    
  - Build the debug version of the library (see the `description above`__)::
    
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --with-openmp --debug
      
    or
      
    ::
      
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --with-pthreads --debug
    
  - Build the non-optimized version of the library (see the `description above`__)::
  
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --with-openmp --no-optimizations  

    or
      
    ::
    
      SPKG_LOCAL=~/solvers/superlu ./spkg-install --with-pthreads --no-optimizations

You may choose any installation destination you like, provided that you have 
write access to it (the target directory will be created if it doesn't exist).    
Also note that the multithreaded library may coexist with the sequential version
in the same directory. For advanced configuration possibilities, please read
src/README, the `manual`_ or visit the `SuperLU home page`_.

Once the library has been built and installed, you may delete the temporary 
directory with the unpacked package to save some disk space (but note that this 
will prevent you from stepping through the library code during debugging), or 
just remove the object files by running

::

  cd src
  make clean 

Now go to the directory with hermes{1|2|3}d. Create the file CMake.vars with the 
following lines (or append to the existing one)::

  set(WITH_SUPERLU YES)
  set(SUPERLU_ROOT ~/solvers/superlu) # or your own installation destination
  set(SUPERLU_MT   YES)
  set(WITH_OPENMP  YES)   # set to NO to use Pthreads rather than OpenMP

Finally execute::

  rm CMakeCache.txt
  cmake .
  make
    
Hermes{1|2|3}d will now be compiled and linked with the multithreaded SuperLU 
library. Before running the parallel calculation, you just need to set the 
environment variable ``OMP_NUM_THREADS`` to the number of threads you wish to 
employ for solution of your system (this is typically the number of cores in your 
multicore machine). For example, on my dual-core laptop I could run

::

  cd hermes2d/tutorial/03-poisson
  OMP_NUM_THREADS=2 ./poisson

Note that you use the variable ``OMP_NUM_THREADS`` with both OpenMP and Pthreads
versions of SuperLU.

Windows MSVC
~~~~~~~~~~~~

http://crd.lbl.gov/~xiaoye/SuperLU/faq.html

MAC OS
~~~~~~

http://www.bleedingmind.com/index.php/2010/07/31/compiling-superlu-on-os-x/
