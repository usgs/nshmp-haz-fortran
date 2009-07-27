#ifndef NSHM_HAZARDCURVE_H
#define NSHM_HAZARDCURVE_H

#define NSHM_CURVE_MAX_NUM_POINTS 20

typedef struct _NSHM_HazardCurveMeta {
	/* Place holder for now. */
	/**
	 * Period
	 * Spectral Acceleration
	 * Name
	 * ID
	 * Other info?
	 */
	 int num_points;
} NSHM_HazardCurveMeta;

typedef struct _NSHM_HazardCurve {
	double lat;
	double lng;
	double xvals[NSHM_CURVE_MAX_NUM_POINTS];
	double yvals[NSHM_CURVE_MAX_NUM_POINTS];
	NSHM_HazardCurveMeta metainfo;
} NSHM_HazardCurve;


void nshm_print_hazard_curve(NSHM_HazardCurve *_curve);

#endif
