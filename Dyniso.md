# Introduction to Dyniso #

Dyniso solves the incompressible Navier-Stokes equations in a box
with the intent of simulating the decay of homogeneous, isotropic
turbulence.

## Building Dyniso ##

Assuming a Linux or Darwin (MacOS-X) system with the GCC compilers in your path and a [threaded FFTW 2.1.5](BuildFFTW.md) available at:
`$(HOME)/local/fftw-omp`
simply do:
```
cd src
make opt
```

Note:  FFTW v2 is **not** required to run Dyniso as several system-supplied FFT libraries are also supported.

## Running Dyniso ##
```
./dyniso.exe < test.inp
```
to select the number of OpenMP threads, do:
```
env OMP_NUM_THREADS=4 ./dyniso.exe < test.inp
```
The code will output a number of files including:

  1. Time statistics in stat?.dat
  1. Spectra at selected timesteps in spectra.`*`
  1. Plot3d visualization files in out.{xyz,q}

## Notes ##
  1. The Plot3d files are in Fortran unformatted (not C-binary)
  1. On MacOS-X (Darwin), you need to install a version of gfortran.  I use the version from [MacPorts](http://www.macports.org) which you can select using
```
sudo port select --set gcc gcc46
```
> and you return to the original Apple compilers using
```
sudo port select --set gcc none
```

## Examples ##

![http://dyniso.googlecode.com/svn/trunk/docs/homo-iso-mag.png](http://dyniso.googlecode.com/svn/trunk/docs/homo-iso-mag.png)

Contours of velocity magnitude from Dyniso on a 384<sup>3</sup> mesh (with dealiasing this is really a 576<sup>3</sup> mesh), computed using 8 threads on an Apple `MacBook` Pro with Intel Corei7 processors.  Visualized using the Plot3d reader in [Paraview](http://www.paraview.org/).

![http://dyniso.googlecode.com/svn/trunk/docs/homo-iso-vort.png](http://dyniso.googlecode.com/svn/trunk/docs/homo-iso-vort.png)

Contours of voriticyt magnitude from Dyniso on a 384<sup>3</sup> mesh (with dealiasing this is really a 576<sup>3</sup> mesh), computed using 8 threads on an Apple `MacBook` Pro with Intel Corei7 processors.  Visualized using the Plot3d reader in [Paraview](http://www.paraview.org/).