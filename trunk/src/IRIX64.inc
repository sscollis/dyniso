
# Makefile Flags for 64bit IRIX

DIST = $(HOME)/dist

COMPATH	= /usr/bin
CXX	= CC
CC	= cc
FC	= f90

# Needed to use ITL interface in IRIX

override ADDONS := $(ADDONS)

# dependency check flag

DEP = -M

# parallel 

ifdef PARALLEL
override ADDONS := $(ADDONS) PARALLEL
#MPIDIR  = $(DIST)/mpi
MPICXX	= CC
MPICC	= cc
MPIINC  = 
MPILIB	= -lmpi
endif

#==============================================================================
# Math libraries
#==============================================================================

# ATLAS (fastest)

#ATLASDIR = $(DIST)/atlas/IRIX_SGIIP30
#ATLASDIR = $(DIST)/atlas/IRIX_SGIIP27
#ATLASDIR = $(DIST)/atlas/IRIX_SGIIP35

#MATHINC  = -I$(ATLASDIR)
#MATHLIB  = -L$(ATLASDIR) -llapack -lcblas -latlas -lblas -lm

# SCSL  (is slowest?)

#override ADDONS := $(ADDONS) SCSL
#override ALL := $(ALL)
#MATHLIB = -lscs -lm

# complib.sgimath  (in-between)

#override ALL := $(ALL) clapack
#MATHLIB = -lcomplib.sgimath -lm

# FFTW 

#FFTDIR  = $(DIST)/fftw
#FFTINC  = -I$(FFTDIR)/include
#FFTLIB  = -L$(FFTDIR)/lib -lrfftw -lfftw

# Boost library

#BOOSTINC = -I$(DIST)/boost

# ITL includes

#ITLINC = -I$(DIST)/itl

# optimization settings

# ip27 R10k sgi-Origin 2000
# ip30 R10k sgi-Octane
# ip35 R12k sgi-Origin 3000

OPTXX	= -n32 -O2 -LANG:std -mips4 -OPT:Olimit=0
MOPTXX	= -n32 -Ofast -IPA -LANG:std -mips4 -OPT:Olimit=0:roundoff=3:IEEE_arithmetic=3:alias=typed:ldx=on:liberal_ivdep:div_split:fast_exp

OPT	= $(OPTXX)
MOPT	= $(MOPTXX)

OPTF	= -n32 -r8 -DR8 -O2 -mp -DSGI_FFT -DRANF -OPT:IEEE_arithmetic=3
MOPTF	= $(OPTF)

DEBUG	= -n32 -LANG:std -g
DEBUGF	= -n32 -g

# loader and loader settings

LD	= $(FC)
DEBUGLD	= $(DEBUG)
MOPTLD	= $(MOPTXX)
OPTLD	= $(OPTXX)

INCS    = $(FFTINC) $(MPIINC) $(MATHINC) $(BOOSTINC) $(ITLINC)
LIBS	= $(FFTLIB) $(MPILIB) $(MATHLIB)

# X stuff

XLIBDIRS = -L/usr/X11/lib
XINCDIRS = -I/usr/include/X11
