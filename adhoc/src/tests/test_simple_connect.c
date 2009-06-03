#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "adhoc_connection.h"
#include "adhoc_api.h"

int main(int argc, char ** argv) {
	AdHoc_Connection conn;

	if (argc != 4) {
		printf("Usage: %s <username> <password> <sid>\n", argv[0]);
		return EXIT_FAILURE;
	}

	adhoc_connect(argv[1], argv[2], argv[3], &conn);

	if (conn.connected) {
		printf("Connected to database.\n");
		adhoc_close_connection(&conn);
		if (conn.connected) {
			printf("Failed to close connection.\n");
		} else {
			printf("Successfully closed connection.\n");
		}
	} else {
		printf("Connection failed.\n");
	}

	return EXIT_SUCCESS;
}
