#ifndef AGRID_META_H
#define AGRID_META_H

#define NSHM_AGRID_META_DESC_LEN 50

typedef struct _NSHM_AgridMeta {
	int    id;           /* Unique DB id for corresponding  agrid data */
	float  min_lat;      /* Minimum latitude value for corresponding agrid */
	float  max_lat;      /* Maximum latitude value for corresponding agrid */
	float  inc_lat;      /* Degree spacing between latitude values */
	float  min_lng;      /* Minimum longitude value for corresponding agrid */
	float  max_lng;      /* Maximum longitude value for corresponding agrid */
	float  inc_lng;      /* Degree spacing between longitude values */
	int    num_rows;     /* The number of rows fetched from the database */
	char   description[NSHM_AGRID_META_DESC_LEN]; /* Textual description */
} NSHM_AgridMeta;

#endif
