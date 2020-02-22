# Dyniso Source Guide 

## Background

Dyniso solves the incompressible Navier-Stokes equations in a box
with the intent of simulating the decay of homogeneous, isotropic 
turbulence.

File | Purpose
-----|--------
dyniso.F90 | Make source file
misc.F90   | Interface to second function
ranf.c     | Interface to random number generator
cleanup    | Script to remove output files
<name>.inp | Various example input files
<name>.mak | Base makefiles
<name>.inc | Platform specific make includes
