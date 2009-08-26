#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "nshm_api.h"
#include "nshm_grid.h"

void make_small_grid(NSHM_Grid *_grid);

int main(int argc, char ** argv) {
	NSHM_Grid grid;
	int status;

	char * query_delete_meta = "\nDELETE FROM %s WHERE %s = %d;\n";
	char * query_delete_data = "DELETE FROM %s WHERE %s = %d;\n\n";

	if ( argc < 4 ) {
		printf("Usage: %s <user> <pass> <tnsalias>\n", argv[0]);
		return EXIT_FAILURE;
	}

	nshm_initialize(argv[1], argv[2], argv[3]);

	printf("#############################################################\n");
	printf("## Fetching an Agrid                                       ##\n");
	printf("#############################################################\n");
	status = nshm_get_agrid(&grid, "adapt_cy");
	nshm_print_grid(&grid);
	nshm_free_grid(&grid);

	printf("#############################################################\n");
	printf("## Fetching a Bgrid                                        ##\n");
	printf("#############################################################\n");
	status = nshm_get_bgrid(&grid, "BVAL_.1");
	nshm_print_grid(&grid);
	nshm_free_grid(&grid);

	printf("#############################################################\n");
	printf("## Fetching a Mmax                                         ##\n");
	printf("#############################################################\n");
	status = nshm_get_mmax(&grid, "gm_j_6p6_7p1");
	nshm_print_grid(&grid);
	nshm_free_grid(&grid);

	printf("#############################################################\n");
	printf("## Writing small sample grids.                             ##\n");
	printf("#############################################################\n");

	make_small_grid(&grid);

	printf("Inserting an agrid...");
	if (nshm_put_agrid(&grid)==NSHM_RETURN_SUCCESS) {
		printf("success!\n"); } else { printf("failed.\n"); }
	printf("You can delete these rows with the following queries:\n");
	printf(query_delete_meta, "AGRID_META", "AGRID_ID", grid.grid_id);
	printf(query_delete_data, "AGRID", "AGRID_ID", grid.grid_id);
	
	printf("Inserting a bgrid...");
	if (nshm_put_bgrid(&grid)==NSHM_RETURN_SUCCESS) {
		printf("success!\n"); } else { printf("failed.\n"); }
	printf("You can delete these rows with the following queries:\n");
	printf(query_delete_meta, "BVALUE_META", "BVAL_ID", grid.grid_id);
	printf(query_delete_data, "BVALUE", "BVAL_ID", grid.grid_id);

	printf("Inserting an Mmax grid...");
	if (nshm_put_mmax(&grid)==NSHM_RETURN_SUCCESS) {
		printf("success!\n"); } else { printf("failed.\n"); }
	printf("You can delete these rows with the following queries:\n");
	printf(query_delete_meta, "MMAX_META", "MMAX_ID", grid.grid_id);
	printf(query_delete_data, "MMAX", "MMAX_ID", grid.grid_id);

	nshm_free_grid(&grid);
	nshm_cleanup();

	return EXIT_SUCCESS;
}

void make_small_grid(NSHM_Grid *_grid) {
	char * name = "Test Sample Grid";

	// This grid setup will have 12 points in it
	_grid->lat_min = 30.0;
	_grid->lat_max = 40.0;
	_grid->lat_inc = 2.0;
	_grid->lng_min = -120.0;
	_grid->lng_max = -118.0;
	_grid->lng_inc = 2.0;

	_grid->grid_values = calloc(12, sizeof(double));
	_grid->grid_values[0] = 1.0;
	_grid->grid_values[1] = 2.0;
	_grid->grid_values[2] = 3.0;
	_grid->grid_values[3] = 4.0;
	_grid->grid_values[4] = 5.0;
	_grid->grid_values[5] = 6.0;
	_grid->grid_values[6] = 7.0;
	_grid->grid_values[7] = 8.0;
	_grid->grid_values[8] = 9.0;
	_grid->grid_values[9] = 10.0;
	_grid->grid_values[10] = 11.0;
	_grid->grid_values[11] = 12.0;

	_grid->grid_name = calloc(strlen(name) + 1, sizeof(char));
	strncpy(_grid->grid_name, name, strlen(name));
}
