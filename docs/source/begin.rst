Overview
########

This simulation was developed to simulate the laser-electron bunch collision and written in Fortran. The electrons are described by weighted computational particles and their motion are governed by Lorentz force equation. The laser energy is assumed to be constant, which is true unless the electron charge is in the order of nanoCoulomb. Therefore, Maxwell's equations are not solved and the laser fields are treated as functions in space and time. The particle and laser have three-dimensional components, i.e., :math:`x, y, z, v_{x,y,z}, E,B_{x,y,z}`. Additional :ref:`features <features>` for laser and electron are available. The simulation setup is depicted as follow:

.. figure:: /figures/laser.png

.. admonition:: Note!

   The electron bunch is treated as a two-dimensional sheath. This approximation is valid for near head-on collision. The longitudinal bunch dimension is not supported at the moment. 

.. _features:

Features
########

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
############

The initial program format was FORTRAN77 (CONTINUE, GO TO etc.) with fortran90 features (modules, allocatable array etc.). Standard functions of Massage Passing Interface (MPI) and OpenMP are used. To install the code, the following dependencies are needed in your PC as well as cluster. Please refer to the links below for the installation of these dependencies.

* fortran compiler (e.g. `gfortran <https://gcc.gnu.org/wiki/GFortran>`_)
   
   * Intel Fortran compiler is not supported yet.

* Massage Passing Interface, MPI (e.g. `OpenMPI <http://www.open-mpi.org>`_)

* OpenMP 

   * OpenMP is shipped along with GCC 4.2 and above.

QuickStart
##########

Download the code with git command:

::

   git clone https://github.com/StevE-Ong/LEC.git

or use the Download Zip button.

This code was written on unix based operating system such as macOS and Linux. Installation and execution are performed by using command line. In the command line, change to the code directory. 

Install via csh script:

.. code-block:: csh

   ./LEC.cmp

This script run the following command:

.. code-block:: csh

   #!/bin/csh -x
   set COMPILER=mpif90
   set dir=$PWD
   set SRC0=$dir/src/LEC.f
   set SRC1=$dir/src/laser.f
   set LOADM=$dir/bin/LEC
   set OBJDIR=$dir/obj

   mkdir $dir/obj/
   mkdir $dir/bin/

   #
   #
   $COMPILER -O3 -fopenmp -fopt-info-all -mcmodel=medium -fconvert=big-endian -g
   -fbacktrace -fbounds-check -J$OBJDIR -o $LOADM $SRC0 $SRC1 >& cmp.lst
   #
   echo " ------------    End of  Compile   --------------   "
   #

The compiler option ``-fopenmp`` is required for programs with OpenMP functions. The option ``-fbounds- check`` is used to detect segmentation errors. This option may be excluded except for debug run. The option ``-O3`` specifies the third level of optimisation. The option ``-o`` specify the name of the executable file. Other compliler options are listed `here <https://gcc.gnu.org/onlinedocs/gcc/Option-Index.html#Option-Index_op_letter-M>`_. If compilation is completed successfully, an executable file ``LEC`` is generated and located in ``/bin``. The file ``cmp.lst`` is the compilation log file. Any code error for an unsuccessful compilation will be written here.  

Running simulation

.. code-block:: csh

   echo Data | mpirun -np 32 ./bin/LEC

where ``Data`` is a folder (can be a folder path, e.g. /examples/Data1) contains the input file ``input.dat``. The ``mpirun`` command is used to run the executable file with MPI library. The option ``-n`` or ``-np`` specifying the number of MPI processes (i.e. 32 processes). The number of threads can be specified as

.. code-block:: csh

   export OMP_NUM_THREADS=4

In this case the number of threads is ``4``. 

To run multiple simulations at one execution, the output files with each input are prepared. In :ref:`this examples <examples>` the output files are located in ``/examples/Data1`` and ``/examples/Data2``. The output directory can be changed to your own/preferred directory followed by ``$i`` without space. For example:

.. code-block:: csh

   #!/bin/csh 

   set i = 1

   while($i<3)
       echo "Running simulation "$i
       echo examples/Data$i | mpirun -np 1 ./bin/LEC
       @ i++
   end 

Then run the simulation.

.. code-block:: csh

   ./LEC_multirun.jcf

The following output will be displayed in the command line. In this case, 1 MPI process is used with 64 threads. Radiation emission is calculated. The calculation will terminate with final runtime specified. If there are more than one simulation, the run will continue with similar output.

::

   Running simulation 1
   Specify output directory
   Output directory:examples/Data1/                                                                                     
  
  
 	###L      ########E    #######C      
 	###L      ########E   ##########C  	 
 	###L      ###E       ###C    ###C    
 	###L      ########E  ###C            
 	###L      ########E  ###C            
 	###L      ###E       ###C    ###C    
 	#######L  ########E   ##########C    
 	LASER###  ELECTRON#    COLLISION     
    
   Welcome to Laser Electron Collision code (v-1.3.0)
  
   *****************************************************
   The code is running on           1  processors
 		                 64  threads
   *****************************************************
   Produce radiation

   Iteration =        1000 ; Time step =   2.3852819683908048E-017
   Iteration =        2000 ; Time step =   4.7705639367816095E-017
   Iteration =        3000 ; Time step =   7.1558459051724140E-017
   .
   .
   .
   Iteration =     9998000 ; Time step =   2.3848049119971265E-013
   Iteration =     9999000 ; Time step =   2.3850434401939656E-013
   Iteration =    10000000 ; Time step =   2.3852819683908046E-013
   Calculating angular distribution...
   calculating radiation
   Final runtime =          36 seconds
   Running simulation 2
   Specify output directory
   Output directory:examples/Data2/ 
   .
   .
