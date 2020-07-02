Basic
=====

Input parameters
----------------

The code was written in ``fortran90``. The input parameters are contained in the file ``input.dat`` in the folder ``Data``.

.. code-block:: fortran

   &PARAM1  jobno=001,ksmax=1000000,div=60.d0,div2=40.d0	             ,&END
   &PARAM2  SL=1.0d22,Ev=10.d6,pw=0.8d-6,pp=5.6d-6,sp=1.44d-6                ,&END
   &PARAM3  alpha=0.0d-6,enum=1.0d9,bin=1.d6,shot=1,inc_ang = 0.d0           ,&END
   &PARAM4  xinit=1.1d0,rmass=1.d0,sigmax=0.01d0,sigmay=0.01d0,sigmaz=0.01d0 ,&END
   &PARAM5  iconR=2,QED=1,ipl=0,shape=0,load=0,OutRad=0,OutPairs=0 	     ,&END

**jobno**  The job numbering. Integer.

**ksmax**  Maximum number of time step. Integer.

**div, div2** The number of division/cells per Lamor radius.

**SL** Laser intensity [:math:`\mathrm{W~cm^{-2}}`].

**Ev** Electron energy [eV].

**pw** Laser wavelength [meter].

**pp** Laser pulse length [meter]. Pulse duration is :math:`\tau_\mathrm{L}=\mathrm{pp}/c`. Pulse duration at FWHM is :math:`\tau_\mathrm{L}\times 1.1744`

**sp** Laser waist radius (at :math:`1/e^2`) [meter].

**alpha** Electron waist radius [meter]. Set to 0.0d-6 for a single electron. If set to a value larger than 0, electron energy distribution in 1D and 2D will be outputted.

**enum** Electron number in a bunch. Used for radiation calculation.

**bin** Radiation spectrum bin size [eV].

**shot** Number of shot. Integer.

**inc_ang** Incident angle of the electron with respect to the x-axis [degree].  

**xinit** Initial position of the electron from the collision point. The collision happens at :math:`t=0`. The laser is one pulse length away from :math:`t=0`.

**rmass** The charge to mass ratio of the colliding particle. Setting to ``rmass=-1`` represents positron, ``rmass=-1836`` represents proton.

**sigma_{x,y,z}** Electron momentum energy spreads in three directions.

**iconR** Specifying the particle pusher used. Integer.
   * ``iconR=0`` Lorentz force
   * ``iconR=1`` Sokolov
   * ``iconR=2`` Reduced Landau-Liftshitz
   * ``iconR=3`` Stochastic

**QED** Specifying whether to use QED for Sokolov and Reduced Landau-Liftshitz. ``QED=1`` is mandatory for Stochastic process.

**ipl** Specifying laser polarisation. ``ipl=0`` for linear polarisation and ``ipl=1`` for circular polarisation.

**shape** Laser spatial and temporal profile. Laser profile can be modified in *laser.f*. Integer.
   * ``shape=0`` 0th order Gaussian beam.
   * ``shape=1`` 5th order paraxial approximation.

.. admonition:: Note!

   The temporal profile for 5th order paraxial approximation is not Gaussian. The temporal profile is :math:`g=\mathrm{cosech}((t-x)/\tau_\mathrm{L})`.

**load** Load particle energy distribution from external file. The default filename is ``f_E_smoothed_new_final.txt``. To change the filename, please edit the subroutine ``manual_load``. Integer.

**OutRad** Specifying whether to calculate radiation. When setting ``OutRad=1``, emission spectrum, photon number distribution, radiation angular distribution will be calculated. This part consumes most of the simulation time. Integer.

**OutRad** Specifying whether to calculate pair production. Currently support for Bethe-Heitler pair production. The cross section for Bremsstrahlung and pair production will be calculated. The Z component of nucleus for specific converter is specify in ``module random_commom``. Integer .

Example
-------

Example of a single electron for ``Lorentz vs Sokolov``. The input parameters are in *examples/Data1(2)*

.. code-block:: fortran

   &PARAM1  jobno=001,ksmax=10000000,div=128.d0,div2=40.d0	            ,&END
   &PARAM2  SL=1.0d22,Ev=100.d6,pw=0.8d-6,pp=5.6d-6,sp=1.44d-6              ,&END
   &PARAM3  alpha=0.0d0,enum=1.0d0,bin=1.d6,shot=1,inc_ang = 0.d0           ,&END
   &PARAM4  xinit=5.d0,rmass=1.d0,sigmax=0.01d0,sigmay=0.01d0,sigmaz=0.01d0 ,&END
   &PARAM5  iconR=0,QED=0,ipl=0,shape=0,load=0,OutRad=1,OutPairs=0 	    ,&END

.. code-block:: fortran

   &PARAM1  jobno=002,ksmax=10000000,div=128.d0,div2=40.d0	            ,&END
   &PARAM2  SL=1.0d22,Ev=100.d6,pw=0.8d-6,pp=5.6d-6,sp=1.44d-6              ,&END
   &PARAM3  alpha=0.0d0,enum=1.0d0,bin=1.d6,shot=1,inc_ang = 0.d0           ,&END
   &PARAM4  xinit=5.d0,rmass=1.d0,sigmax=0.01d0,sigmay=0.01d0,sigmaz=0.01d0 ,&END
   &PARAM5  iconR=1,QED=0,ipl=0,shape=0,load=0,OutRad=1,OutPairs=0 	    ,&END

Multiple run can be performed by executing 

.. code-block:: csh

   ./LEC_multirun.jcf

The number ``i`` in this file is set according to the numbering of the folder.

.. code-block:: csh

   #!/bin/csh 

   set i = 1

   while($i<3)
       echo "Running simulation "$i
       echo examples/Data$i | mpirun -np 1 ./bin/LEC
       @ i++
   end 

The electron trajectories

.. figure:: /figures/trajectories.png

The time evolution of electron energy

.. figure:: /figures/energies.png

The radiation spectrum

.. figure:: /figures/spectra.png

The photon number distribution

.. figure:: /figures/photonnumber.png

The radiation angular distribution

.. figure:: /figures/angular_dist.png
