**Dyniso** is intended to be a very-simple, 3d homogeneous, isotropic turbulence simulator using OpenMP (SMP node-level parallelization) with the primary purpose to serve as an educational example.

**Dyniso** currently does **not** support MPI-type parallelism which effectively limits its use to rather small problems.  There are other, open-source isotropic turbulence codes that have more features and that support MPI domain-decomposition parallelism that are likely better for more serious use and you are encouraged to explore those projects.

The primary features of **Dyniso**:
  * Incompressible Navier-Stokes equations
  * Large-Eddy Simulation (LES) using
    * Constant coefficient Smagorinsky model
    * Dynamic Smagorinsky model
  * Fourier pseudo-spectral method in space
  * RK2, RK3, and explicit Euler time integration
  * Variety of initial conditions
    * Synthetic turbulence initial condition
    * Taylor-Green vortex initial condition
    * Restart capability
  * Fortran90 implementation with:
    * OpenMP loop-level parallelism
    * Support for vectorization (SSE)
    * Single or double precision support
  * Builds with:
    * GCC (with Gfortran) on Linux and Darwin
    * ICC (with ifort) on Linux
    * IRIX (hard to find nowadays)
  * Visualization output in Plot3d format

**Dyniso** does require access to a compatible FFT library on your system.  Options include:
  * FFTW v2
  * Intel MKL
  * SGI Complib
  * Cray FFT libraries
And it is easy to implement interfaces to other FFT libraries.

For more information, see [Introduction to Dyniso](https://code.google.com/p/dyniso/wiki/Dyniso#Introduction).