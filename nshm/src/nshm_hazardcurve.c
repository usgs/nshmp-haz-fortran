#include <stdlib.h>
#include <stdio.h>

#include "nshm_hazardcurve.h"

void nshm_print_hazard_curve(NSHM_HazardCurve *_curve) {
	int i;

	printf("Latitude: %f\tLongitude: %f\n", _curve->lat, _curve->lng);

	printf("+-----------------+---------------------------+\n");
	printf("|  Ground Motion  |  Frequency of Exceedance  |\n");
	printf("+-----------------+---------------------------+\n");

	for (i = 0; i < _curve->metainfo.num_points; i++) {
		printf("|  %13g  |  %23g  |\n", _curve->xvals[i], _curve->yvals[i]);
		printf("+-----------------+---------------------------+\n");
	}
}
