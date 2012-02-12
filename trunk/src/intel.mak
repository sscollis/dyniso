#=============================================================================
#
#  Makefile for homogeneous turbulence code (Intel)
#
#  Author:  Scott Collis
#
#  Revised: 6-7-98
#
#=============================================================================

# Intel + MKL (dealiasing doesn't work since only 2^n is supported)

#OPT    = -O3 -DRANF -DINTEL_MKL
#LIB    = -Vaxlib -L/usr/local/intel/mkl/lib/32 \
#         /usr/local/intel/mkl/lib/32/libmkl_p4.a /usr/lib/libpthread.a

# Intel + FFTW

#OPT   = -O -r8 -DR8 -DRANF -DFFTW -openmp
OPT    = -O -r8 -DR8 -axW -DRANF -DFFTW
#OPT   = -O -r8 -DR8 -axW -DRANF -DFFTW -openmp
LIB    = -Vaxlib -L../lib -lrfftw -lfftw 

FLAGS  = -cpp -w $(OPT)
LFLAGS = -static $(OPT)
FC     = ifc
CC     = icc

.SUFFIXES: .f90 .F90

dyniso: dyniso.o ranf.o misc.o
	$(FC) $(LFLAGS) dyniso.o ranf.o misc.o -o dyniso $(LIB)

clean:
	/bin/rm *.o *.mod *.d

.f90.o:
	$(FC) $(FLAGS) -c $*.f90 

.F90.o:
	$(FC) $(FLAGS) -c $*.F90 

.f.o:
	$(FC) $(FLAGS) -c $*.f

.c.o:
	$(CC) -O3 -c $*.c
