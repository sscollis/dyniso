Dyniso solves the incompressible Navier-Stokes equations in a box
with the intent of simulating the decay of homogeneous, isotropic 
turbulence.

Building Dyniso:

Assuming a Linux system with the Intel compilers in your path
and a threaded FFTW 2.1.5 available at:  $(HOME)/local/fftw-thread
simply do:

cd src
make opt

Running Dyniso:

./dyniso.exe < test.inp

to select the number of OpenMP threads, do:

env OMP_NUM_THREADS=4 ./dyniso.exe < test.inp

The code will output a number of files including:

1. Time statistics in stat?.dat
2. Spectra at selected timesteps in spectra.*
3. Plot3d visualization files in out.{xyz,q}

Notes:
1. The Plot3d files are in fortran unformated (not C-binary)
2. On MacOS-X, you need to install a version of gfortran.
   I use the version from MacPorts which you can select using

   sudo port select --set gcc gcc46

   and you return to the original Apple compilers using

   sudo port select --set gcc none


S. Scott Collis
Thu Feb  9 14:28:14 MST 2012
