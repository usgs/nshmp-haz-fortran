/**
 * File Name: test_statement.c
 * Makefile Target to compile: test
 * Author: Eric Martinez
 * 
 * This is a simple test application to check the proper functionality of the
 * AdHoc_Stmt structure and corresponding functions. Due to the nature of the
 * test (i.e. connecting to a database etc) this test encapsulates the
 * test_connect.c tests as well but since. Basically if this test passes, then
 * so will the test_connect test since we cannot hope to use a statement on an
 * invalid connection. However we keep the old test available in case someting
 * happens with the statement we can individually test to ensure connectivity
 * still functions.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "adhoc_api.h"
#include "adhoc_stmt.h"
#include "adhoc_connection.h"
#include "adhoc_utils.h"

#define OUTPUT_FIELD_INDEX 0

/**
 * 06/11/09 -- EMM:
 * You can toggle the SQL_STR and NUM_FIELDS between the two options below by
 * commenting one of them out. The difference is one will return a specific
 * agrid name and the other will return all the agrid names (in lexicographical
 * order). This is a contrived example/test of the adhoc statement code so don't
 * put too much thought into how this is working.
 */

/* Select a specific one */ 
#define SQL_STR "SELECT AGRID_NAME FROM AGRID_META WHERE AGRID_ID = :id"
#define NUM_FIELDS 2

/* Select all names */ /*
#define SQL_STR "SELECT AGRID_NAME FROM AGRID_META ORDER BY AGRID_ID"
#define NUM_FIELDS 1
*/

int main (int argc, char ** argv) {
	AdHoc_Connection connection;
	AdHoc_Stmt       stmt;
    int              num_rows = 5;
	int i;

	if ( argc < 4 ) {
		fprintf(stderr, "Usage: %s <username> <password> <tnsalias>\n",
			argv[0]);
		return EXIT_FAILURE;
	}

	AdHoc_Field fields[2] = {
		{
			// Max 10 rows in result set. 50 characters per row
			(char **) calloc(10, sizeof(char)*50), (int *) 50, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_STR, "1",
			(void *) 0
		},
		{
			(void *) 1, (int *) -1, (short) NULL, (unsigned short *) NULL, 
			(unsigned short *) NULL, SQLT_INT, ":id", (void *) 0
		}
	};

	adhoc_connect(argv[1], argv[2], argv[3], &connection);

	stmt = adhoc_stmt_create(SQL_STR, fields, NUM_FIELDS, &connection);

	adhoc_stmt_execute(&stmt, &connection);
	adhoc_stmt_fetch(&stmt, &connection, &num_rows);

	printf("Query returned %d rows.\n", num_rows);
	for (i = 0; i < num_rows; ++i) {
    	printf("%s\n",
			(char *) 
				&(stmt.fields[OUTPUT_FIELD_INDEX].value) + 
				( i * 50 * sizeof(char) )
			);
	}

	adhoc_stmt_free(&stmt, &connection);
	adhoc_close_connection(&connection);

	return EXIT_SUCCESS;

}
