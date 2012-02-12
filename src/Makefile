ifndef ARCH
  ARCH   = $(shell uname)
endif
PROG     = dyniso.exe
VPATH    = ../src
LIBDIR   = ../src
HDRDIR   = ../src
INCDIR   = ../src
DEFINES  =
CPPFLAGS =

# Auto sense OS and Machine type
OSVER  = $(shell uname -r)
MACH = $(shell uname -m)

# Set defaults (use OpenMP with double precision)

ifeq ($(USE_OPENMP),0)
  override undefine USE_OPENMP
else
  USE_OPENMP = 1
endif

ifeq ($(USE_FLOAT),0)
  override undefine USE_FLOAT
endif

ALL = dyniso 
CALL = misc

# ARCHitecture specific options

ifeq ($(ARCH),IRIX)                     # Silicon Graphics IRIX
include $(INCDIR)/IRIX.inc
endif

ifeq ($(ARCH),IRIX64)                   # Silicon Graphics IRIX64
include $(INCDIR)/IRIX64.inc
endif

ifeq ($(ARCH),Linux)                    # Linux
include $(INCDIR)/Linux.inc
endif

ifeq ($(ARCH),ICC)                      # Intel CC compiler
include $(INCDIR)/ICC.inc
endif

ifeq ($(ARCH),Darwin)                   # GCC compiler with MPICH
include $(INCDIR)/Darwin.inc
endif

DEF = $(foreach def, $(DEFINES) $(ADDONS), -D$(def))
override FPPFLAGS := $(DEF) $(FPPFLAGS)
override CPPFLAGS :=  $(DEF) $(CPPFLAGS)

OBJS = $(foreach module, $(ALL) $(CALL) $(SPECIAL), $(module).o)

# make variables available for next make
export

help:
	@echo ''
	@echo ' Dyniso makefiles:'
	@echo ''
	@echo ' Run make with one of the following options:'    
	@echo ' a) dbx.....dbx version'
	@echo ' b) opt.....optimized version for '$(ARCH)
	@echo ' c) mopt....maximally optimized version for '$(ARCH)
	@echo ' e) clean...remove all objects'
	@echo ''
	@echo ' ARCH = '$(ARCH)', OSVER = '$(OSVER)', Machine = '$(MACH)
opt:	
	$(MAKE) LTYPE=opt -f dyniso.mak
mopt:
	$(MAKE) LTYPE=mopt -f dyniso.mak
dbx: 
	$(MAKE) LTYPE=g -f dyniso.mak
clean:
	rm -rf *.o *.mod *.exe
