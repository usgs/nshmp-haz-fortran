#include <stdlib.h>
#include <stdio.h>

#include "adhoc_api.h"

#include "nshm_api.h"
#include "nshm_agrid_meta.h"
#include "nshm_agrid.h"


void print_agrid(NSHM_Agrid _agrid);

int main(int argc, char ** argv) {

	NSHM_Agrid agrid;

	/* Initialize the DB connection */
	adhoc_init(NSHM_AUTH[INT_DEV]);

	nshm_get_random_agrid(&agrid);

	print_agrid(agrid);

	adhoc_close();
	return EXIT_SUCCESS;
}

void print_agrid(NSHM_Agrid _agrid) {
	int i;

	printf("########################################");
	printf("########################################\n");
	printf("# Agrid Meta Information\n");
	printf("########################################");
	printf("########################################\n");

	printf("\tID:          %d\n", _agrid.metadata->id);
	printf("\tNumRows:     %d\n", _agrid.metadata->num_rows);
	printf("\tDescription: %s\n", _agrid.metadata->description);

	printf("########################################");
	printf("########################################\n");
	printf("# Agrid Data\n");
	printf("########################################");
	printf("########################################\n");



	for (i = 0; i < 10; ++i) {
		printf("\t%3.1f  %4.1f  %6.5f\n",
			_agrid.latitude[i], _agrid.longitude[i], _agrid.value[i]);
	}
}
