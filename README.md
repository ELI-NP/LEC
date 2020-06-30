Particle code for laser-electron collision
==========================================
[![Language](https://img.shields.io/badge/language-Fortran90-blue.svg)](https://fortran-lang.org)
[![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause)


## Overview
- This package contain particle code for simulating laser-electron collision
- Free software: BSD license
- Documentation: to do
## Features
- Radiation Reaction (RR)
  - Landau-Liftshitz (LL)
  - Sokolov (SKL)
  - Stochastic (SCS)
- Classical and QED-assisted model for LL & SKL
- Radiation model
  - Nonlinear Compton Scattering
- Simulating a single electron or an electron bunch
- Calculating radiation spectrum and photon number distribution
- Calculating radiation angular distribution
- Pair production
  - Bethe-Heitler
  - nonlinear Breit-Wheeler (to do)  
- Laser field beyond paraxial approximation up to fifth order
- Spatial and temporal Gaussian shape laser pulse (Linear & Circular polarisation)
- Load electron beam energy distribution from experimental data
## QuickStart
#### Requirements:
- fortran compiler (e.g. gfortran)
- Massage Passing Interface, MPI (e.g. OpenMPI)
- OpenMP 
#### Install via Makefile:
```
Makefile
```
#### Install via csh script:
```
./LEC.cmp
```
#### Running simulation:
```
echo Data | mpirun -np 32 ./bin/LEC
```
where `Data` is a folder contains the `input.dat` file. To set the number of threads `export OMP_NUM_THREADS=n`, where `n` is the number of thread. To run multiple simulations at one run, prepare `DataN` with each input and execute `./LEC_multirun.jcf`. Please modify this file for `N`. 
The code are used for the following publications:
- J. F. Ong, T. Moritaka, and H. Takabe, Radiation Reaction in the interaction of Ultraintense laser with matter and gamma ray source, Physics of Plasmas 23, 053117 (2016).
- J. F. Ong, T. Moritaka, and H. Takabe, The suppression of radiation reaction and laser field depletion in laser-electron beam interaction, Physics of Plasmas 25, 033113 (2018).
- J. F. Ong, T. Moritaka, and H. Takabe, Optimizing the energies conversion in laser- electron beam collision, Physics of Plasmas 26, 033102 (2019).
- T. Moritaka et-al., J. Phys. Conf. Ser. 454 012016 (2013).
