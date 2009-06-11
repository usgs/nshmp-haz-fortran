CC = gcc
CFLAGS = -g -Wall -O3 -fno-strict-aliasing

OCI_INCLUDE = -I$(ORACLE_HOME)/rdbms/public
OCI_LIBRARY = -L$(ORACLE_HOME)/lib
OCI_LINKS   = -lsqlplus -lclntsh -lnnz11
OCI_FLAGS   = $(OCI_INCLUDE) $(OCI_LIBRARY) $(OCI_LINKS)

INCLUDE = -Iinclude
LIBRARY = -Llib

TARGETS = adhoc_utils adhoc_connection adhoc_stmt adhoc_api

OBJECTS = build/adhoc_utils.o build/adhoc_connection.o \
          build/adhoc_stmt.o build/adhoc_api.o

all: init $(TARGETS)

lib: all
	ar rvs lib/libadhoc.a $(OBJECTS)

test: lib
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) -o test_connect \
		src/tests/test_connect.c -ladhoc
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) -o test_statement \
		src/tests/test_statement.c -ladhoc
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) -o test_insert \
		src/tests/test_insert.c -ladhoc
	

adhoc_utils: init
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

adhoc_api: adhoc_utils
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

adhoc_connection: adhoc_api
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

adhoc_stmt: init adhoc_connection 
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

init: clean
	mkdir build

clean:
	-rm -rf build
	-rm -f test_connect
	-rm -f test_statement
	-rm -f test_insert
	-rm -f lib/libadhoc.a
