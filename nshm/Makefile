CC = gcc
CFLAGS = -g -Wall -O0 -fno-strict-aliasing

F77 = /home/webm1/bin/irun/bin/gfortran

OCI_INCLUDE = -I$(ORACLE_HOME)/rdbms/public
OCI_LIBRARY = -L$(ORACLE_HOME)/lib
OCI_LINKS   = -lsqlplus -lclntsh -lnnz11
OCI_FLAGS   = $(OCI_INCLUDE) $(OCI_LIBRARY) $(OCI_LINKS)

ADHOC_HOME = "../adhoc"
ADHOC_INCLUDE = -I$(ADHOC_HOME)/include
ADHOC_LIBRARY = -L$(ADHOC_HOME)/lib

INCLUDE = -Iinclude $(ADHOC_INCLUDE)
LIBRARY = -Llib $(ADHOC_LIBRARY)

TARGETS = nshm_api nshm_agrid
OBJECTS = build/nshm_api.o build/nshm_agrid.o

all: init $(TARGETS)

lib: all
	ar rvs lib/libadhoc.a $(OBJECTS)

test: all
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
	-o test src/tests/test_agrid_fetch.c -ladhoc

ftest: all
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_INCLUDE) \
		-o build/test_agrid_fortran.o -c src/tests/test_agrid_fortran.c
	$(F77) -g $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) $(OBJECTS) \
		build/test_agrid_fortran.o -o ftest src/tests/test_agrid.f -ladhoc

nshm_api:
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

nshm_agrid: 
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

init: clean
	mkdir build

clean:
	-rm -rf build
	-rm -f test
	-rm -f ftest
