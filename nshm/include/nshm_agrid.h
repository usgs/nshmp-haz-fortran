#ifndef AGRID_H
#define AGRID_H

#include "nshm_agrid_meta.h"

#define NSHM_MAX_AGRID_LENGTH 1048576  /* At most a MB of space */

typedef struct _NSHM_Agrid {
	NSHM_AgridMeta *metadata;               /* Agrid metadata */
	float *latitude;  /* Latitude coordinates */
	float *longitude; /* Longitude coordinates */
	float *value;     /* Agrid Values */
} NSHM_Agrid;

int nshm_get_random_agrid(NSHM_Agrid *_agrid);

 
#endif
