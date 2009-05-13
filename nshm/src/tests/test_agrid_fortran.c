#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "oci.h"
#include "adhoc_api.h"
#include "nshm_api.h"
#include "nshm_agrid.h"
#include "nshm_agrid_meta.h"

boolean initialized = FALSE;

void fetchagrid_(float * values[], NSHM_AgridMeta * meta, long int desc_len)
{
    NSHM_Agrid agrid;

	if (!initialized) {
		adhoc_init(NSHM_AUTH[INT_DEV]);
		initialized = TRUE;
	}

	nshm_get_random_agrid(&agrid);

	memmove(values, agrid.value, agrid.metadata->num_rows * sizeof(float));

	meta->id = agrid.metadata->id;
	meta->num_rows = agrid.metadata->num_rows;
	memset(meta->description, 0, sizeof(meta->description));
	strncpy(meta->description, agrid.metadata->description,
		strlen(agrid.metadata->description)
	);
}
