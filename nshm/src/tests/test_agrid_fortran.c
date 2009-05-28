#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "oci.h"
#include "adhoc_api.h"
#include "nshm_api.h"
#include "nshm_agrid.h"
#include "nshm_agrid_meta.h"

boolean initialized = FALSE;

void fetchagrid_(float ** values, NSHM_AgridMeta * meta)
{
    NSHM_Agrid agrid;

	if (!initialized) {
		adhoc_init(NSHM_AUTH[INT_DEV]);
		initialized = TRUE;
	}

	nshm_get_random_agrid(&agrid);

	meta->id = agrid.metadata->id;
	meta->num_rows = agrid.metadata->num_rows;
	meta->min_lat = agrid.metadata->min_lat;
	meta->max_lat = agrid.metadata->max_lat;
	meta->min_lng = agrid.metadata->min_lng;
	meta->max_lng = agrid.metadata->max_lng;

	memset(meta->description, 0, NSHM_AGRID_META_DESC_LEN);
	strncpy(meta->description, agrid.metadata->description,
		strlen(agrid.metadata->description)
	);

	memmove(values, agrid.value, agrid.metadata->num_rows * sizeof(float));

}
