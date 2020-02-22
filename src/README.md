# Dyniso Source Guide 

## Background

Dyniso solves the incompressible Navier-Stokes equations in a box
with the intent of simulating the decay of homogeneous, isotropic 
turbulence.

File         | Purpose
-------------|---------------------------------------------------
`dyniso.F90` | Make source file
`misc.F90`   | Interface to second function
`ranf.c`     | Interface to random number generator
`cleanup`    | Script to remove output files
`*.inc`      | Platform specific make includes
`*.inp`      | Various example input files
`*.mak`      | Base makefiles (out dated, kept for reference)

### Notes
1. The `F90` files require the CPP preprocessor 
2. Just typeing `make` will git the makefile options
3. The normal mode of building is `make opt`

