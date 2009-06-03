CC = gcc
CFLAGS = -g -Wall -O3 -fno-strict-aliasing

OCI_INCLUDE = -I$(ORACLE_HOME)/rdbms/public
OCI_LIBRARY = -L$(ORACLE_HOME)/lib
OCI_LINKS   = -lsqlplus -lclntsh -lnnz11
OCI_FLAGS   = $(OCI_INCLUDE) $(OCI_LIBRARY) $(OCI_LINKS)

INCLUDE = -Iinclude
LIBRARY = -Llib

TARGETS = adhoc_utils adhoc_check_error adhoc_connection adhoc_api

OBJECTS = build/adhoc_utils.o build/adhoc_check_error.o \
          build/adhoc_connection.o build/adhoc_api.o

all: init $(TARGETS)

lib: all
	ar rvs lib/libadhoc.a $(OBJECTS)

test: lib
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) -o test_simple_db_conn \
		src/tests/test_simple_db_conn.c -ladhoc
	$(CC) $(CFLAGS) $(INCLUDE) $(LIBRARY) $(OCI_FLAGS) -o test_simple_connect \
		src/tests/test_simple_connect.c -ladhoc

adhoc_api: adhoc_check_error adhoc_utils
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

adhoc_connection: init adhoc_check_error
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

adhoc_check_error: init adhoc_utils
	$(CC) $(CFLAGS) $(INCLUDE) $(OCI_INCLUDE) -o build/$@.o -c src/$@.c

adhoc_utils: init
	$(CC) $(CFLAGS) $(INCLUDE) -o build/$@.o -c src/$@.c

init: clean
	mkdir build

clean:
	-rm -rf build
	-rm -f test_simple_db_conn
	-rm -f test_simple_connect
	-rm -f lib/libadhoc.a
