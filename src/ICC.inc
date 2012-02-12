
# Makefile Flags for Linux using ICC compiler 

LOCAL = $(HOME)/local

COMPATH	=
CXX	= icpc 
CC	= icc 
FC	= ifort 

override ADDONS := $(ADDONS)

# dependency check flag

DEP = -M

# MPI parallel 

ifdef PARALLEL
override ADDONS := $(ADDONS) PARALLEL
MPIDIR  = $(LOCAL)/mpi
MPICXX	= mpicxx 
MPICC	= mpicc 
MPIINC  = 
MPILIB	= -lmpi
endif

#==============================================================================
# Math libraries
#==============================================================================

ifdef USE_OPENMP
OMP = -omp
endif

ifdef USE_FLOAT
FLOAT = -s
endif

# ATLAS

#ATLASDIR = $(LOCAL)/atlas/IRIX_SGIIP30
#ATLASDIR = $(LOCAL)/atlas/IRIX_SGIIP27
#ATLASDIR = $(LOCAL)/atlas/IRIX_SGIIP35

#MATHINC  = -I$(ATLASDIR)
#MATHLIB  = -L$(ATLASDIR) -llapack -lcblas -latlas -lblas -lm

# SCSL

#override ADDONS := $(ADDONS) SCSL
#override ALL := $(ALL)
#MATHLIB = -lscs -lm

# complib.sgimath

#override ALL := $(ALL) clapack
#MATHLIB = -lcomplib.sgimath -lm

# FFTW 

FFTDIR  = $(LOCAL)/fftw$(OMP)$(FLOAT)
FFTINC  = -I$(FFTDIR)/include
ifdef USE_OPENMP
FFTLIB  = -L$(FFTDIR)/lib -lrfftw_threads -lfftw_threads -lrfftw -lfftw -lm
else
FFTLIB  = -L$(FFTDIR)/lib -lrfftw -lfftw -lm
endif

# Boost library

#BOOSTINC = -I$(LOCAL)/boost

# ITL includes

#ITLINC = -I$(LOCAL)/itl

# optimization settings

ifeq ($(USE_OPENMP),1)
  OPENMP = -openmp
endif

ifdef USE_FLOAT
  DOUBLE =
else
  DOUBLE = -r8 -DR8
endif

OPTXX	= -O2
MOPTXX	= -O3

OPT	= $(OPTXX)
MOPT	= $(MOPTXX)

OPTF	= -O2 $(DOUBLE) -axW -DFFTW -DRANF $(OPENMP)
MOPTF	= -fast $(DOUBLE) -axW -DFFTW -DRANF $(OPENMP)

DEBUG	= -g $(DOUBLE) $(OPENMP)
DEBUGF	= -g $(DOUBLE) $(OPENMP)

# loader and loader settings

LD	= $(FC)
DEBUGLD	= $(DEBUGF)
MOPTLD	= $(MOPTF)
OPTLD	= $(OPTF) 

INCS    = $(FFTINC)
LIBS	= $(FFTLIB)

# X stuff

XLIBDIRS = -L/usr/X11/lib
XINCDIRS = -I/usr/include/X11
