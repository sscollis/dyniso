#=============================================================================
#
#  Makefile for homogeneous turbulence code
#
#  Author:  Scott Collis
#
#  Revised: 2-12-2012
#
#=============================================================================

# Parser generator
LEX = flex -I
YACC = bison -d

# set the compile flags
ifeq ($(LTYPE),g)	# turn on debug flag 
CXXFLAGS := $(DEBUG) $(CCWOFF)
CFLAGS := $(DEBUG) $(CWOFF)
FFLAGS := $(DEBUGF) $(FWOFF)
LDFLAGS	:= $(DEBUGLD)
else			# Maximal optimization flags
ifeq ($(LTYPE),mopt)
CXXFLAGS := $(MOPTXX)
CFLAGS := $(MOPT)
FFLAGS := $(MOPTF)
LDFLAGS	:= $(MOPTLD)
else			# Regular optimization flags
ifeq ($(LTYPE),opt)
CXXFLAGS := $(OPTXX)
CFLAGS := $(OPT)
FFLAGS := $(OPTF)
LDFLAGS	:= $(OPTLD)
endif
endif
endif

override CPPFLAGS := $(CPPFLAGS) $(INCS)
override FPPFLAGS := $(CPPFLAGS)
override LIBS := $(LIBS)
override PROG := $(PROG)

OBJS = $(foreach module, $(ALL) $(CALL) $(SPECIAL), $(module).o)
SRCS = $(foreach module, $(ALL) $(SPECIAL), $(module).F90)

.SUFFIXES: .F90 .F .cpp .c .d

$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) -o $(PROG).tmp $(OBJS) $(LIBS)
	mv $(PROG).tmp ./$(PROG)

.F90.o:
	$(FC) $(FPPFLAGS) $(FFLAGS) -c $(VPATH)/$*.F90 

.F.o:
	$(FC) $(FPPFLAGS) $(FFLAGS) -c $(VPATH)/$*.F

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $(VPATH)/$*.cpp 

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $(VPATH)/$*.c

#  Automatically make dependencies
%.d: %.cpp
	set -e; $(CXX) $(DEP) $(CPPFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

#  Automatically make dependencies
%.d: %.c
	set -e; $(CC) $(DEP) $(CPPFLAGS) $< \
	| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
	[ -s $@ ] || rm -f $@

.PHONY: clean
clean :
	rm -r $(OBJS)
