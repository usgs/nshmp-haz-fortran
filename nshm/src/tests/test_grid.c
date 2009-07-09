#include <stdlib.h>

#include "nshm_api.h"
#include "nshm_grid.h"

int main(int argc, char ** argv) {
	NSHM_Grid grid;
	int status;

	if ( argc < 4 ) {
		printf("Usage: %s <user> <pass> <tnsalias>\n", argv[0]);
		return EXIT_FAILURE;
	}

	nshm_initialize(argv[1], argv[2], argv[3]);

	printf("#############################################################\n");
	printf("## Fetching an Agrid                                       ##\n");
	printf("#############################################################\n");
	status = nshm_get_agrid(&grid);
	nshm_print_grid(&grid);
	nshm_free_grid(&grid);

	printf("#############################################################\n");
	printf("## Fetching a Bgrid                                        ##\n");
	printf("#############################################################\n");
	status = nshm_get_bgrid(&grid);
	nshm_print_grid(&grid);

	printf("#############################################################\n");
	printf("## Fetching a Mmax                                         ##\n");
	printf("#############################################################\n");
	status = nshm_get_mmax(&grid);
	nshm_print_grid(&grid);

	nshm_cleanup();
	return EXIT_SUCCESS;
}
