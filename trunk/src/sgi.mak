#=============================================================================
#
#  Makefile for homogeneous turbulence code (SGI)
#
#  Author:  Scott Collis
#
#  Revised: 6-7-98
#
#=============================================================================

#.... Use random function (same as vectoral)

OPT      = -r8 -DR8 -O2 -mp -DSGI_FFT -DRANF -OPT:IEEE_arithmetic=3
LIB      = -lcomplib.sgimath

#.... Use drand48 function

#OPT     = -n32 -O2 -DRANF -OPT:IEEE_arithmetic=3:roundoff=3
#LIB     = -lcomplib.sgimath drand48.o

#.... Use the rand function

#OPT     = -n32 -O2 -OPT:IEEE_arithmetic=3:roundoff=3
#LIB     = -lcomplib.sgimath

FLAGS    = -cpp $(OPT)
COMP     = f90

.SUFFIXES: .F90

dyniso: dyniso.o ranf.o
	$(COMP) $(FLAGS) dyniso.o ranf.o -o dyniso $(LIB)

vmsiso: vmsiso.o ranf.o
	$(COMP) $(FLAGS) vmsiso.o ranf.o -o vmsiso $(LIB)

tmp: tmp.o ranf.o
	$(COMP) $(FLAGS) tmp.o ranf.o -o tmp $(LIB)

clean:
	/bin/rm *.o *.mod

ensemble: ensemble.f90
	$(COMP) -r8 -O2 ensemble.F90 -o ensemble

.F90.o:
	$(COMP) $(FLAGS) -c $*.F90 

.F.o:
	f77 $(FLAGS) -c $*.F

.c.o:
	cc -O3 -c $*.c
