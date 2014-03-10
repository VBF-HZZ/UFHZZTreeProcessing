CXX=g++
CC=gcc
CFLAGS=-O2 -Wall
LDFLAGS=-Llib
PRFFLAGS=-lProof
THRDFLAGS=-lThread
INS=-I$(ROOTSYS)/include/root
INS2=-I$(ROOFITSYS)/include
INSS=-I./include

LD1=-L$(ROOFITSYS)/lib

CFLAGS += `root-config --cflags`
LIBS += `root-config --glibs`

LDa=-lRooFitCore
LDb=-lRooFit

OBJ1=processTrees.o

.PHONY: clean all main test

all: processTrees

processTrees: processTrees.o
	$(CXX) -o processTrees.exe $(OBJ1) $(LIBS)

makeROC: makeROC.o
	$(CXX) -o makeROC.exe makeROC.o $(LIBS)

clean:
	@rm *.o *.exe *~ 


##############RULES##############
.cc.o:
	$(CXX) $(CFLAGS) $(INS) -c $<
.cpp.o:
	$(CXX) $(CFLAGS) $(INS) -c $<


