CC = gcc
CFLAGS = -g -Wall -O0 -fno-strict-aliasing

# Make sure this is the correct compiler
F90 = gfortran -g
# On gldplone this is /home/webm1/bin/irun/bin/gfortran

OCILIB_HOME = "/data/ocilib"
OCI_INCLUDE = -I$(OCILIB_HOME)/include
OCI_LIBRARY = -L$(OCILIB_HOME)/lib
OCI_LINKS   = -locilib
OCI_FLAGS   = $(OCI_INCLUDE) $(OCI_LIBRARY) $(OCI_LINKS)

INCLUDE = -Iinclude
LIBRARY = -Llib

TARGETS = nshm_api.o nshm_grid.o nshm_grid_f.o nshm_hazardcurve.o \
          nshm_geo_grid.o
OBJECTS = build/nshm_api.o build/nshm_grid.o build/nshm_hazardcurve.o \
          build/nshm_geo_grid.o
FOBJS   = build/nshm_grid_f.o

all: init $(TARGETS)

lib: all
	ar rvs lib/libnshm.a $(OBJECTS) $(FOBJS)

test: all
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		-o test_api src/tests/test_api.c 
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		-o test_grid src/tests/test_grid.c
	$(F90) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) $(FOBJS) \
		-o test_grid_f src/tests/test_grid.f 
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		-o test_hazardcurve src/tests/test_hazardcurve.c
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		-o test_geo_grid src/tests/test_geo_grid.c

nshm_api.o: src/nshm_api.c nshm_grid.o
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@ -c src/nshm_api.c

nshm_grid.o: src/nshm_grid.c nshm_api.o
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@ -c src/nshm_grid.c

nshm_grid_f.o: src/nshm_grid.f
	$(F90) -o build/$@ -c src/nshm_grid.f
	mv grid.mod include/grid.mod

nshm_hazardcurve.o: src/nshm_hazardcurve.c
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@ -c \
		src/nshm_hazardcurve.c

nshm_geo_grid.o: src/nshm_geo_grid.c
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@ -c \
		src/nshm_geo_grid.c

init:
	-mkdir build

clean:
	-rm -rf build
	-rm -f test_api
	-rm -f test_grid
	-rm -f include/grid.mod
	-rm -f test_grid_f
	-rm -f lib/libnshm.a
	-rm -f test_hazardcurve
	-rm -f test_geo_grid
