#include <stdlib.h>
#include <stdio.h>

#include "nshm_hazardcurve.h"

int main(int argc, char ** argv) {
	
	int i;
	NSHM_HazardCurve * curve = calloc(1, sizeof(NSHM_HazardCurve));
	double g[19] = {
		0.005, 0.007, 0.010, 0.014, 0.019, 0.027,
		0.038, 0.053, 0.074, 0.103, 0.145, 0.203,
		0.284, 0.397, 0.556, 0.778, 1.090, 1.520,
		2.130
	};
	double fex[19] = {
		8.030E-01, 7.068E-01, 5.880E-01, 4.550E-01, 3.218E-01, 2.068E-01,
		1.216E-01, 6.583E-02, 3.380E-02, 1.681E-02, 7.637E-03, 3.092E-03,
		1.093E-03, 3.557E-04, 1.121E-04, 3.314E-05, 7.738E-06, 1.207E-06,
		6.705E-08,
	};

	curve->metainfo.num_points = 0;
	curve->lat = 35.0;
	curve->lng = -118.0;

	for ( i = 0; i < 19; i++) {
		curve->xvals[i] = g[i];
		curve->yvals[i] = fex[i];
		curve->metainfo.num_points++;
	}

	nshm_print_hazard_curve(curve);

	free(curve);
	return EXIT_SUCCESS;
}
