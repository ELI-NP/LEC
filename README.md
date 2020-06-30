# LEC
Particle code for laser-electron collision
## Overview
- This package contain particle code for simulating laser-electron collision
- Free software: BSD license
- Documentation: to do
## Features
- Radiation Reaction (RR)
  - Ladau-Liftshitz
  - Sokolov
  - Stochastic
- Classical and QED-assisted model of RR
- Simulating a single electron or an electron bunch
- Calculating radiation spectrum and photon number distribution
- Calculating radiation angular distribution
## QuickStart
#### Requirements:
- fortran compiler (e.g. gfortran)
- Massage Passing Interface, MPI (e.g. OpenMPI)
#### Install via Makefile:
```Makefile```
#### Install via csh:
```./LEC.cmp```
#### Running simulation:
```echo Data | mpirun -np 32 ./bin/LEC```
