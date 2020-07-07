Input parameters
================

This section describes the details of the input parameters, output files, and examples. The code was written in ``fortran90``. The input parameters are contained in the file ``input.dat`` in the folder ``Data``.

.. code-block:: fortran

   &PARAM1  jobno=001,ksmax=1000000,div=60.d0,div2=40.d0	             ,&END
   &PARAM2  SL=1.0d22,Ev=10.d6,pw=0.8d-6,pp=5.6d-6,sp=1.44d-6                ,&END
   &PARAM3  alpha=0.0d-6,enum=1.0d9,bin=1.d6,shot=1,inc_ang = 0.d0           ,&END
   &PARAM4  xinit=1.1d0,rmass=1.d0,sigmax=0.01d0,sigmay=0.01d0,sigmaz=0.01d0 ,&END
   &PARAM5  iconR=2,QED=1,ipl=0,shape=0,load=0,OutRad=0,OutPairs=0 	     ,&END

**jobno**    The job numbering. Integer.

**ksmax**    Maximum number of time step. Integer.

**div,div2** The number of division or cells per Lamor radius.

**SL**       Laser intensity [:math:`\mathrm{W~cm^{-2}}`].

**Ev**       Electron energy [eV].

**pw** Laser wavelength [meter].

**pp** Laser pulse length [meter]. Pulse duration is :math:`\tau_\mathrm{L}=\mathrm{pp}/c`. Pulse duration at FWHM is :math:`\tau_\mathrm{L}\times 1.1744`

**sp** Laser waist radius (at :math:`1/e^2`) [meter].

**alpha** Electron waist radius [meter]. Set to 0.0d-6 for a single electron. If set to a value larger than 0, electron energy distribution in 1D and 2D will be outputted.

**enum** Electron number in a bunch. Used for radiation calculation. For single electron, ``enum=1.d0``

**bin** Radiation spectrum bin size [eV].

**shot** Number of shot. Integer.

**inc_ang** Incident angle of the electron with respect to the x-axis [degree].  

**xinit** Initial position of the electron from the collision point. The collision happens at :math:`t=0`. The laser is one pulse length away from :math:`t=0`.

**rmass** The charge to mass ratio of the colliding particle. Setting to ``rmass=-1`` represents positron, ``rmass=-1836`` represents proton.

**sigma_{x,y,z}** Electron momentum energy spreads in three directions. This input is ignored for a single electron

**iconR** Specifying the particle pusher used. Integer.

   * ``iconR=0`` Lorentz force
   * ``iconR=1`` Sokolov
   * ``iconR=2`` Reduced Landau-Liftshitz
   * ``iconR=3`` Stochastic

**QED** Specifying whether to use QED for Sokolov and Reduced Landau-Liftshitz. ``QED=1`` is mandatory for Stochastic process.

**ipl** Specifying laser polarisation. 
   * ``ipl=0`` linear polarisation (LP)
   * ``ipl=1`` circular polarisation (CP).

**shape** Laser spatial and temporal profile. Laser profile can be modified in *laser.f*. Integer.

   * ``shape=0`` 0th order Gaussian beam.
   * ``shape=1`` 5th order paraxial approximation.

.. admonition:: Note!

   The temporal profile for 5th order paraxial approximation is not Gaussian. The temporal profile is :math:`g=\mathrm{cosech}((t-x)/\tau_\mathrm{L})`.

**load** Load particle energy distribution from external file. The default filename is ``f_E_smoothed_new_final.txt``. To change the filename, please edit the subroutine ``manual_load``. Integer.

**OutRad** Specifying whether to calculate radiation. When setting ``OutRad=1``, emission spectrum, photon number distribution, radiation angular distribution will be calculated. This part consumes most of the simulation time. Integer.

**OutPairs** Specifying whether to calculate pair production. The code currently support the Bethe-Heitler pair production. The cross section for Bremsstrahlung and pair production will be calculated if ``OutPairs=1``. The Z component of nucleus for the specific converter is specify in ``module random_commom``. Integer.


Output files
============

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

The outputs are written in ASCII format. The file ``output001.dat`` records the detail parameters of the simulation. For example:

.. code-block:: fortran

   Parameters for pulse laser
  
   Laser polarization: linear
  
   0th order Gaussian beam
  
   Laser Intensity               [W/cm2]   1.0000000000000000E+022
   Peak electric field             [V/m]   274000000000000.00     
   Peak magnetic field           [Gauss]   9280000000.0000000     
   Larmor radius for light speed     [m]   1.8318965517241380E-009
   laser wavelength                  [m]   7.9999999999999996E-007
   pulse length                      [m]   5.5999999999999997E-006
   pulse duration                    [s]   1.8666666666666665E-014
   pulse duration (FWHM)             [s]   2.1978133333333330E-014
   waist radius (1/e2)               [m]   1.4400000000000000E-006
  
   parameters for electron beam
   ...

The file ``orbt1q001.dat`` records the trajectories, energy etc. of the particle. For a single electron, there are 7 files recoding the same output. For example:

.. code-block:: fortran

   -0.466547E-13     0.279928E-04    -0.350081E-14    -0.196692E+03    -0.193187E-03     0.100511E+09     0.491296E-07     0.191798E-03     0.205673E-05
   -0.466428E-13     0.279857E-04    -0.141161E-13    -0.196692E+03    -0.390705E-03     0.100511E+09     0.100948E-06     0.392282E-03     0.208935E-05
   -0.466309E-13     0.279785E-04    -0.319586E-13    -0.196692E+03    -0.590080E-03     0.100511E+09     0.152945E-06     0.596550E-03     0.209571E-05
   -0.466190E-13     0.279714E-04    -0.570502E-13    -0.196692E+03    -0.788800E-03     0.100511E+09     0.202518E-06     0.799480E-03     0.207561E-05
   -0.466070E-13     0.279642E-04    -0.893211E-13    -0.196692E+03    -0.984349E-03     0.100511E+09     0.247101E-06     0.995990E-03     0.202914E-05
   ...

The values of each column from the left to right are: time [s], x [m], y [m], :math:`p_x` [normalized], :math:`p_y` [normalized], kinetic energy [eV], work [eV], radiation energy [eV], :math:`\chi_e` [dimensionless]. 

The file ``phtne001.dat`` records the radiation output. For example:

.. code-block:: fortran

   8333.3333333333321        238.86944345907790        6.2633513616217478     
   25000.000000000000        182.50447244093024        7.5507336932200397     
   41666.666666666664        104.14244180601469        9.1380529004804760     
   58333.333333333328        89.422617071344263        9.8579995567498120     
   75000.000000000000        70.619159337234422        10.697476969356655     
   91666.666666666657        63.363841302196569        11.199199843048401 
   ... 

The values in the column from the left to right are the energy [eV], photon number, photon number :math:`\times` energy.

The file ``phtnTe001.dat`` records the radiation angular distribution. For example:

.. code-block:: fortran

   -0.3139E+01    -0.3139E+01     0.0000E+00
   -0.3132E+01    -0.3139E+01     0.0000E+00
   -0.3126E+01    -0.3139E+01     0.0000E+00
   -0.3120E+01    -0.3139E+01     0.0000E+00

The values in the column from the left to right are :math:`\theta_z` [rad], :math:`\theta_y` [rad], radiated energy [a.u].

For an electron bunch, there are more than 7 outputs, depending on the number of MPI processes. Each output record a sample electron information. On the other hand, file such as ``AveEne(jobno).dat``, ``dist_fn(kstep)(jobno).dat``, ``dist_fn2d(kstep)(jobno).dat`` will be output. 

The file ``AveENE`` record the time [s], average kinetic energy [eV], average radiation energy [eV], average + :math:`\sigma` [eV], average - :math:`\sigma` [eV], where :math:`\sigma` is the standard deviation of electron bunch energy. 

The file ``dist_fn`` records energy [eV], electron number [a.u]. 

The file ``dist_fn2d`` records :math:`p_y` [normalized], :math:`p_z` [normalized], electron number [a.u]. 

Python
------

In this examples, the visualisation is performed by using Python in `Jupyter notebook <https://jupyter.org>`_. The python codes can be found in ``/examples/**.ipynb``. The extension ``.ipynb`` stand for Jupiter notebook. In the Jupyter notebook, there is a python function ``import figformat``. This function output/display figures with selected parameters. The figure width, **fig_width** is set to 3.4 inches, represents a single column width of a double column journal. The figure width can be override to any number by writing ``fig.set_size_inches(fig_width*2,fig_width/1.618)`` at each plot. The number ``1.618`` is the Golden ratio. Multiplying or dividing the **fig_width** by the Golden ratio for figure height ensure the nice appearance of a figure. Other parameters such as font size, plot line width, ticks width and etc. can be changed in the file ``figformat.py``.

Gnuplot
-------

On the other hand, a quick visualisation can be performed by using `gnuplot <http://www.gnuplot.info>`_. For example:

:: 

   > plot “***.dat” using ($1):($4) with lines 
   > replot “***.dat” using ($1):($4) with lines

.. _examples:

Examples
========

Single electron
---------------

In this example, we plot several outputs of a single electron. Details of the plotting code can be referred to the Jupyter notebook. It can be viewed in GitHub.

The electron trajectory 

.. figure:: /figures/trajectories.png

The time evolution of electron energy

.. figure:: /figures/energies.png

The radiation spectrum

.. figure:: /figures/spectra.png

The photon number distribution

.. figure:: /figures/photonnumber.png

The radiation angular distribution

.. figure:: /figures/angular_dist.png

Electron bunch
--------------

.. todo:: To do


Models
======

.. todo:: To do

   Details numerical implementation can be obtained in Ref. :cite:`mypop`.

Landau-Liftshitz
----------------

.. math::

   \frac{ dv^{\mu}}{d\tau}=\frac{e}{mc}F^{\mu\nu}v_{\nu}+\tau_{0}\left( \frac{e}{mc} \dot{F}^{\mu\nu} v_{\nu}+\frac{e^{2}}{m^{2}c^{2}}F^{\mu\nu}F_{\alpha\nu}v^{\alpha}
   \frac{e^{2}}{m^{2}c^{2}}(F^{\alpha\nu}v_{\nu})(F_{\alpha\lambda}v^{\lambda})v^{\mu}\right)

Sokolov
-------

.. math::

   \frac{ dp^{\mu}}{d\tau}=\frac{e}{mc}F^{\mu\nu}v_{\nu}-\frac{I_{QED}}{mc^2}p^{\mu}+\tau_{0}\frac{e^{2}}{(mc)^{2}}\frac{I_{QED}}{I_{E}}F^{\mu\nu}F_{\nu\alpha}p^{\alpha}

Stochastic
----------



Quantum
-------



Emission cross-section
----------------------

.. math::

   dW_{em}=\frac{\alpha mc^{2}}{\sqrt{3}\pi\hbar\gamma}\left[\left(1-\xi+\frac{1}{1-\xi} \right)K_{2/3}(\delta)
   -\int_{\delta}^{\infty}K_{1/3}(s)ds  \right] d\xi

.. math::

   \xi=\frac{\hbar\omega}{\gamma mc^{2}},\:\delta=\frac{2\xi}{3(1-\xi)\chi}

and :math:`K_{\nu}(x)` is modified Bessel function. At classical limit :math:`\chi<<1`

.. math::

   dP&=&\mathcal{E}dW_{em}\nonumber\\ &\rightarrow& \frac{e^{2}\omega_{c}}{ \sqrt{3}\pi c}\frac{1}{\gamma^{2}} 
   \frac{\omega}{\omega}_{c}[2K_{2/3}(\delta)-\int_{\delta}^{\infty}K_{1/3}(s)ds]d\omega

reduced to classical synchrotron radiation where :math:`\omega_{c}` is the critical frequency and :math:`\delta\longrightarrow 2\xi/3\chi`.

.. figure:: /figures/qchi.png

The function :math:`q(\chi_e)~\text{for}~\chi_e\ll 1` (blue)

.. math::

    q(\chi_e\ll 1)\approx 1-\frac{55}{16}\sqrt{3}\chi + 48\chi^2 

The function :math:`q(\chi_e)~\text{for}~\chi_e\gg 1` (green)

.. math::

    q(\chi_e\gg 1)\approx\frac{48}{243}\Gamma(\frac{2}{3})\chi^{-4/3} 
    \left[ 1 -\frac{81}{16\Gamma(2/3)}(3\chi)^{-2/3} \right] 

