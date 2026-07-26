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
CXX = g++
CFLAGS = -O3 -march=native -flto -funroll-loops -fopenmp -std=c++11
LDFLAGS = -flto -fopenmp


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



