#include <stdlib.h>
#include <stdio.h>

#include "adhoc_api.h"

#include "nshm_api.h"
#include "nshm_agrid_meta.h"
#include "nshm_agrid.h"


void print_agrid(NSHM_Agrid * _agrid);

int main(int argc, char ** argv) {

	NSHM_Agrid * agrid;
	agrid = calloc(1, sizeof(NSHM_Agrid));
	agrid->metadata = calloc(1, sizeof(NSHM_AgridMeta));

	/* Initialize the DB connection */
	adhoc_init(NSHM_AUTH[INT_DEV]);

	nshm_get_random_agrid(agrid);

	print_agrid(agrid);

	adhoc_close();
	free(agrid->metadata);
	free(agrid);
	return EXIT_SUCCESS;
}

void print_agrid(NSHM_Agrid * _agrid) {
	int i = 0;

	printf("########################################");
	printf("########################################\n");
	printf("# Agrid Meta Information\n");
	printf("########################################");
	printf("########################################\n");

	printf("\tID:          %d\n", _agrid->metadata->id);
	printf("\tNumRows:     %d\n", _agrid->metadata->num_rows);
	printf("\tLatitude [min, max, inc] = [%4.2f, %4.2f, %4.2f]\n",
		_agrid->metadata->min_lat, _agrid->metadata->max_lat,
		_agrid->metadata->inc_lat);
	printf("\tLongitude [min, max, inc] = [%5.2f, %5.2f, %4.2f]\n",
		_agrid->metadata->min_lng, _agrid->metadata->max_lng,
		_agrid->metadata->inc_lng);
	printf("\tDescription: %s\n", _agrid->metadata->description);

	printf("########################################");
	printf("########################################\n");
	printf("# Agrid Data\n");
	printf("########################################");
	printf("########################################\n");



	for (i = 0; i < 10; ++i) {
		printf("\t%3.1f  %4.1f  %6.5f\n",
			_agrid->latitude[i], _agrid->longitude[i], _agrid->value[i]);
	}
}
