# Building FFTW #

For best performance, dyniso should use FFTW with threads/openmp support built in.  On a Linux platform, I build FFTW using:

env CFLAGS="-O3 -fopenmp" FFLAGS="-O3 -fopenmp" ./configure \
--enable-shared --enable-threads --with-openmp \
--prefix=$HOME/local/fftw-omp