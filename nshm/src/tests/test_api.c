#include <stdlib.h>
#include <stdio.h>

#include "nshm_api.h"
#include "ocilib.h"

#define LOG_ERR_MSG "You should see this twice. The 2nd time is log_error.\n"

int main(int argc, char ** argv) {
	char * stamp = nshm_timestamp();

	if ( argc < 4 ) {
		printf("Usage: %s <user> <password> <tnsalias>\n", argv[0]);
		return EXIT_FAILURE;
	}

	printf("INFO -- [%s] Starting run.\n", stamp);
	nshm_initialize(argv[1], argv[2], argv[3]);


	printf("Connection established.\n");
	printf("Cleaning up connection.\n");
	nshm_cleanup();

	printf(LOG_ERR_MSG);
	nshm_log_error(LOG_ERR_MSG);

	printf("Generating an error for the error handler callback.\n");
	printf("\n\n==============================================\n\n");
	nshm_initialize("foo", "bar", "baz");
	printf("\n\n==============================================\n\n");

	free(stamp);
	return EXIT_SUCCESS;
}
