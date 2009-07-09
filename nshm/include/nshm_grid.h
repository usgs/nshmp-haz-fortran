#ifndef NSHM_GRID_H
#define NSHM_GRID_H

// Corresponds to the order in the queries below
#define GRID_ID_IDX 1
#define GRID_NAME_IDX 2
#define GRID_MIN_LAT_IDX 3
#define GRID_MAX_LAT_IDX 4
#define GRID_INC_LAT_IDX 5
#define GRID_MIN_LNG_IDX 6
#define GRID_MAX_LNG_IDX 7
#define GRID_INC_LNG_IDX 8
#define GRID_VAL_IDX 1

// Query to select Agrid Metadata for a single agrid with known AGRID_ID value
#define QUERY_AGRID_META "SELECT AGRID_ID, AGRID_NAME, MIN_LAT, MAX_LAT, \
	INC_LAT, MIN_LON, MAX_LON, INC_LON FROM AGRID_META WHERE AGRID_ID = %d"

// Query to select Agrid Data for a single agrid with known AGRID_ID value
#define QUERY_AGRID_DATA "SELECT AGRID_VAL FROM AGRID WHERE AGRID_ID = %d \
	ORDER BY LAT DESC, LON ASC"

// Query to select Bgrid Metadata for a single agrid with known BVAL_ID value
#define QUERY_BGRID_META "SELECT BVAL_ID, BVAL_NAME, MIN_LAT, MAX_LAT, \
	INC_LAT, MIN_LON, MAX_LON, INC_LON FROM BVALUE_META WHERE BVAL_ID = %d"

// Query to select Bgrid Data for a single agrid with known BGRID_ID value
#define QUERY_BGRID_DATA "SELECT B_VAL FROM BVALUE WHERE BVAL_ID = %d \
	ORDER BY LAT DESC, LON ASC"

// Query to select MMax Metadata for a single agrid with known MMAX_ID value
#define QUERY_MMAX_META "SELECT MMAX_ID, MMAX_NAME, MIN_LAT, MAX_LAT, \
	INC_LAT, MIN_LON, MAX_LON, INC_LON FROM MMAX_META WHERE MMAX_ID = %d"

// Query to select MMax Data for a single agrid with known MMAX_ID value
#define QUERY_MMAX_DATA "SELECT MMAX_VAL FROM MMAX WHERE MMAX_ID = %d \
	ORDER BY LAT DESC, LON ASC"

typedef struct _NSHM_Grid {
	double lat_min;
	double lat_max;
	double lat_inc;
	double lng_min;
	double lng_max;
	double lng_inc;
	double *grid_values;
	char   *grid_name;
	int    grid_id;
} NSHM_Grid;

int nshm_get_agrid(NSHM_Grid *_grid);
int nshm_get_bgrid(NSHM_Grid *_grid);
int nshm_get_mmax(NSHM_Grid *_grid);

void nshm_print_grid(NSHM_Grid *_grid);
int  nshm_num_rows(NSHM_Grid *_grid);
void nshm_free_grid(NSHM_Grid *_grid);

#endif
