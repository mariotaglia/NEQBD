#Unix makefile for fortran-file	

# put the name of the target program here
TARGET = neqbd
SRC = BD.f90 switch.f90  analysis.f90  aux.f90  calcforce.f90  cellneighbor.f90  energy.f90  forces.f90  initial.f90  IO.f90  move.f90  mpi.f90  randoms.f90  barostat.f90 splines.f90 averagepot.f90

# some definitions
SHELL = /bin/bash
FFLAGS= -O3 # -fbounds-check -fbacktrace # ${F90FLAGS}


FF = mpif77 #${F90}
VER = ~/bin/NEQBD_PERIOD2

GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)
LFLAGS=-cpp -D_VERSION=\"$(GIT_VERSION)\"

all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o) $(LFLAGS) $(LDFLAGS)
	cp $(TARGET) $(VER)

$(SRC:.f90=.o): $(SRC)
	${FF} -c ${FFLAGS}  $(SRC) $(LFLAGS) $(LDFLAGS)

install: all
	cp $(TARGET) $(VER)

clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 


ifeq (.depend, $(wildcard .depend))
include .depend
endif








































