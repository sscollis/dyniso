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

OPT    = -m64 -O3 -fdefault-integer-8 -fdefault-real-8 -DR8 -DRANF -DFFTW
LIB    = -L$(HOME)/local/fftw-m64/lib -lrfftw -lfftw 

FLAGS  = $(OPT) -I$(HOME)/local/fftw-m64/include
CFLAGS = -O3 -m64
LFLAGS = $(OPT) 
FC     = gfortran 
CC     = gcc 

.SUFFIXES: .F90

dyniso: dyniso.o ranf.o misc.o
	$(FC) $(LFLAGS) dyniso.o ranf.o misc.o -o dyniso $(LIB)

clean:
	/bin/rm *.o *.mod

.f90.o:
	$(FC) $(FLAGS) -c $*.f90 

.F90.o:
	$(FC) $(FLAGS) -c $*.F90 

.f.o:
	$(FC) $(FLAGS) -c $*.f

.c.o:
	$(CC) $(CFLAGS) -c $*.c
