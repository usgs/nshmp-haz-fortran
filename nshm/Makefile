CC = gcc
CFLAGS = -g -Wall -O0 -fno-strict-aliasing

# Make sure this is the correct compiler
F90 = gfortran -g
# On gldplone this is /home/webm1/bin/irun/bin/gfortran

OCILIB_HOME = "../OCILIB"
OCI_INCLUDE = -I$(OCILIB_HOME)/include
OCI_LIBRARY = -L$(OCILIB_HOME)/lib
OCI_LINKS   = -locilib
OCI_FLAGS   = $(OCI_INCLUDE) $(OCI_LIBRARY) $(OCI_LINKS)

INCLUDE = -Iinclude
LIBRARY = -Llib

TARGETS = nshm_api nshm_grid nshm_grid_f
OBJECTS = build/nshm_api.o build/nshm_grid.o
FOBJS   = build/nshm_grid_f.o

all: init $(TARGETS)

test: all
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		-o test_api src/tests/test_api.c 
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		-o test_grid src/tests/test_grid.c
	$(F90) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) $(FOBJS) \
		-o test_grid_f src/tests/test_grid.f 

nshm_api:
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

nshm_grid: nshm_api
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

nshm_grid_f:
	$(F90) -o build/$@.o -c src/nshm_grid.f
	mv grid.mod include/grid.mod

init: clean
	mkdir build

clean:
	-rm -rf build
	-rm -f test_api
	-rm -f test_grid
	-rm -f include/grid.mod
	-rm -f test_grid_f
