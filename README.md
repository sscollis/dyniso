# Dyniso Installation and Usage Guide

![Dyniso logo](https://github.com/sscollis/dyniso/blob/master/docs/dyniso-logo.png)

## Background

Dyniso solves the incompressible Navier-Stokes equations in a box
with the intent of simulating the decay of homogeneous, isotropic 
turbulence.

## Building Dyniso

Assuming a Linux system with the Intel compilers in your path
and a threaded FFTW 2.1.5 available at:  `$(HOME)/local/fftw-thread`
simply do:

    cd src
    make opt

Note, if your compilers are not available using different names you may 
need to explicitly tell make where to find them, such as:

    env CC=gcc-9 FC=gfortran-9 make opt

Dyniso has not been updated yet to work with FFTW3.

## Running Dyniso

To execute dyniso, you redirect a namelist input file, such as:

    ./dyniso.exe < test.inp

to select the number of OpenMP threads, do:

    env OMP_NUM_THREADS=4 ./dyniso.exe < test.inp

The code will output a number of files including:

  1. Time statistics in `stat?.dat`
  2. Spectra at selected timesteps in `spectra.*`
  3. Plot3d visualization files in `out.{xyz,q}`

### Notes:
  1. The Plot3d files are in Fortran unformated (not C-binary)
  2. Performance is highly dependant on the FFT library used.
     If using FFTW v2 you should read:  
     http://www.fftw.org/fftw2_doc/fftw_6.html
  3. Note that FFTW must be build with `--enable-threads and
     --with-openmp`

### Isosurface of Velocity Magnitude 
![Isosurface of velocity magnitude](https://github.com/sscollis/dyniso/blob/master/docs/homo-iso-mag.png)


### Isosurface of Vorticity
![Isosurface of vorticity](https://github.com/sscollis/dyniso/blob/master/docs/homo-iso-vort.png)
