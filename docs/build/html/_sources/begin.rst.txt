
Overview
========

* This package contain particle code for simulating laser-electron collision
* Free software: BSD license
* Documentation: *to do*

Features
========

* Radiation Reaction (RR)
   * Landau-Liftshitz (LL)
   * Sokolov (SKL)
   * Stochastic (SCS)
* Classical and QED-assisted model for LL & SKL
* Radiation model
   * Nonlinear Compton Scattering
* Simulating a single electron or an electron bunch
* Calculating radiation spectrum and photon number distribution
* Calculating radiation angular distribution
* Pair production
   * Bethe-Heitler
   * nonlinear Breit-Wheeler (*to do*)  
* Laser field beyond paraxial approximation up to fifth order
* Spatial and temporal Gaussian shape laser pulse (Linear & Circular polarisation)
* Load electron beam energy distribution from experimental data

Dependencies
============

* fortran compiler (e.g. gfortran)
* Massage Passing Interface, MPI (e.g. OpenMPI)
* OpenMP 

QuickStart
==========

Install via csh script:

.. code-block:: csh

   ./LEC.cmp


Running simulation

.. code-block:: csh

   echo Data | mpirun -np 32 ./bin/LEC


where ``Data`` is a folder contains the ``input.dat`` file. To set the number of threads ``export OMP_NUM_THREADS=n``, where ``n`` is the number of thread. To run multiple simulations at one run, prepare `DataN` with each input and execute ``./LEC_multirun.jcf``. Please modify this file for ``N``. 
