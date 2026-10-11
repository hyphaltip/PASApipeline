# -*- makefile -*-
# just make headings for each thing and put the headings
# on the line for ALL
PROJECTDIR = .
BINDIR = .

#CC = g++
#CPPC = ${CC}
#CFLAGS = -g -Wall -O3
#CFLAGS = -ggdb 
#  add -xloopinfo to get feedback about -xautopar
#CFLAGS =  -fast -xtarget=ultra3 -xarch=v9
#LDFLAGS = -lm -fast -xtarget=ultra3 -xarch=v9
#CFLAGS =  -fast -xautopar -xipo=2 -xtarget=ultra3 -xarch=v9 -xlibmil -xlibmopt -xprofile=use:./feedback
#LDFLAGS = -lm -fast -xautopar -xipo=2 -xtarget=ultra3 -xarch=v9 -xlibmil -xlibmopt -xprofile=use:./feedback
#-static

# for the 64-bit opterons (use gcc4):
#CFLAGS = -O3 -m64 -mtune=opteron -march=opteron -mfpmath=sse 

# for everything else
CXX ?= g++

# -msse4.2 is what actually matters here: it is what turns __builtin_popcountll
# into a POPCNT instruction instead of a libgcc __popcountdi2 call.  SSE4.2 is
# present on everything since Nehalem/Bulldozer (~2009), so the binary stays
# portable across a heterogeneous cluster and across conda/container hosts.
# -mtune=native only changes scheduling, never the instruction set, so it is
# safe to combine.  Override ARCHFLAGS=-march=native for a single-machine build.
# -msse4.2 exists only on x86-64; on arm64 (Apple Silicon, aarch64) POPCNT is native.
ifneq ($(filter x86_64 amd64 i386 i686,$(shell uname -m)),)
ARCHFLAGS ?= -msse4.2 -mtune=native
else
ARCHFLAGS ?=
endif

CFLAGS = -O3 $(ARCHFLAGS) -flto -funroll-loops -fopenmp -std=c++11
LDFLAGS = -O3 $(ARCHFLAGS) -flto -fopenmp

# macOS: dyld has no fallback search path for @rpath libraries, and libomp (conda's
# llvm-openmp) is linked as @rpath/libomp.dylib. Without an LC_RPATH the binary builds and
# then dies at start: "Library not loaded: @rpath/libomp.dylib". Add the conda library
# directory (or RPATH_DIRS) as an rpath. Linux finds the library through other paths.
ifeq ($(shell uname -s),Darwin)
RPATH_DIRS ?= $(if $(CONDA_PREFIX),$(CONDA_PREFIX)/lib)
LDFLAGS += $(foreach d,$(RPATH_DIRS),-Wl,-rpath,$(d))
endif


DEPEND_FILES = *.cpp *.c *.h
CLEANABLE_FILES = *.o *~ *\# pasa

ALL = pasa 
OBJECTS = alignment_segment.o cdna_alignment.o common_subs.o  cdna_alignment_assembler.o Lobject.o stringFuncts.o pasa.o argProcessor.o


.SUFFIXES: .cpp .c

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o 

.cpp.o:
	$(CXX) $(CFLAGS) -c $*.cpp -o $*.o 

all:    $(ALL)

pasa: ${OBJECTS} pasa.o
	$(CXX) -o $(BINDIR)/$@ $(OBJECTS) $(LDFLAGS)


clean:
	/bin/rm -f $(CLEANABLE_FILES)



