#=============================================================================
#
#  Makefile for homogeneous turbulence code (Intel)
#
#  Author:  Scott Collis
#
#  Revised: 6-7-98
#
#=============================================================================

# Darwin with GCC 

OPT    = -m64 -O3 -fdefault-integer-8 -fdefault-real-8 -DR8 -DRANF -DINTEL_MKL
LIB    = -L/opt/intel/mkl/10.0.3.020/lib -lmkl_p4.a 
#LIB    = -L$(HOME)/local/fftw/lib -lrfftw -lfftw 

FFLAGS  = $(OPT) -I$(HOME)/local/fftw/include
CFLAGS = -O3
LFLAGS = $(OPT) 
FC     = gfortran 
CC     = gcc 

.SUFFIXES: .F90

dyniso: dyniso.o ranf.o misc.o
	$(FC) $(LFLAGS) dyniso.o ranf.o misc.o -o dyniso $(LIB)

clean:
	/bin/rm *.o *.mod

.f90.o:
	$(FC) $(FFLAGS) -c $*.f90 

.F90.o:
	$(FC) $(FFLAGS) -c $*.F90 

.f.o:
	$(FC) $(FFLAGS) -c $*.f

.c.o:
	$(CC) $(CFLAGS) -c $*.c
