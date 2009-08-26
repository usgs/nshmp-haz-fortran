#ifndef NSHM_GRID_H
#define NSHM_GRID_H

// Number of decimal places to keep
#define NSHM_GRID_PRECISION 3.0

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

// Query to select the Agrid ID based on a given grid name
#define QUERY_AGRID_ID "SELECT AGRID_ID FROM AGRID_META WHERE AGRID_NAME = \
	TRIM('%s')"

// Query to select Agrid Metadata for a single agrid with known AGRID_ID value
#define QUERY_AGRID_META "SELECT AGRID_ID, AGRID_NAME, MIN_LAT, MAX_LAT, \
	INC_LAT, MIN_LON, MAX_LON, INC_LON FROM AGRID_META WHERE AGRID_ID = %d"

// Query to select Agrid Data for a single agrid with known AGRID_ID value
#define QUERY_AGRID_DATA "SELECT AGRID_VAL FROM AGRID WHERE AGRID_ID = %d \
	ORDER BY LAT DESC, LON ASC"

// Query to select the Bgrid ID based on a given grid name
#define QUERY_BGRID_ID "SELECT BVAL_ID FROM BVALUE_META WHERE BVAL_NAME = \
	TRIM('%s')"

// Query to select Bgrid Metadata for a single agrid with known BVAL_ID value
#define QUERY_BGRID_META "SELECT BVAL_ID, BVAL_NAME, MIN_LAT, MAX_LAT, \
	INC_LAT, MIN_LON, MAX_LON, INC_LON FROM BVALUE_META WHERE BVAL_ID = %d"

// Query to select Bgrid Data for a single agrid with known BGRID_ID value
#define QUERY_BGRID_DATA "SELECT B_VAL FROM BVALUE WHERE BVAL_ID = %d \
	ORDER BY LAT DESC, LON ASC"

// Query to select the MMax ID based on a given grid name
#define QUERY_MMAX_ID "SELECT MMAX_ID FROM MMAX_META WHERE MMAX_NAME = \
	TRIM('%s')"

// Query to select MMax Metadata for a single agrid with known MMAX_ID value
#define QUERY_MMAX_META "SELECT MMAX_ID, MMAX_NAME, MIN_LAT, MAX_LAT, \
	INC_LAT, MIN_LON, MAX_LON, INC_LON FROM MMAX_META WHERE MMAX_ID = %d"

// Query to select MMax Data for a single agrid with known MMAX_ID value
#define QUERY_MMAX_DATA "SELECT MMAX_VAL FROM MMAX WHERE MMAX_ID = %d \
	ORDER BY LAT DESC, LON ASC"


// Sequence names for the meta table xxx_ID columns.
#define AGRID_SEQ "AGRID_ID_SEQ"
#define BGRID_SEQ "BVAL_ID_SEQ"
#define MMAX_SEQ "MMAX_ID_SEQ"

// Sequence names for the data table xxx_OR columns.
#define AGRID_OR_SEQ "AGRID_OR_SEQ"
#define BGRID_OR_SEQ "BVAL_OR_SEQ"
#define MMAX_OR_SEQ "MMAX_OR_SEQ"

// DML to insert Agrid metadata
#define DML_AGRID_META "INSERT INTO AGRID_META (AGRID_ID, AGRID_NAME, MIN_LAT, \
	MAX_LAT, INC_LAT, MIN_LON, MAX_LON, INC_LON) VALUES (%s.nextval, :name, \
	:lat_min, :lat_max, :lat_inc, :lng_min, :lng_max, :lng_inc)"

// DML to insert Bgrid metadata
#define DML_BGRID_META "INSERT INTO BVALUE_META (BVAL_ID, BVAL_NAME, MIN_LAT, \
	MAX_LAT, INC_LAT, MIN_LON, MAX_LON, INC_LON) VALUES (%s.nextval, :name, \
	:lat_min, :lat_max, :lat_inc, :lng_min, :lng_max, :lng_inc)"

// DML to insert Mmax metadata
#define DML_MMAX_META "INSERT INTO MMAX_META (MMAX_ID, MMAX_NAME, MIN_LAT, \
	MAX_LAT, INC_LAT, MIN_LON, MAX_LON, INC_LON) VALUES (%s.nextval, :name, \
	:lat_min, :lat_max, :lat_inc, :lng_min, :lng_max, :lng_inc)"

// DML to insert Agrid data
#define DML_AGRID_DATA "INSERT INTO AGRID (AGRID_ID, AGRID_VAL, LAT, LON, \
	AGRID_OR) VALUES (%d, :vals, :lats, :lngs, %s.nextval)"

// DML to insert Bgrid data
#define DML_BGRID_DATA "INSERT INTO BVALUE (BVAL_ID, B_VAL, LAT, LON, \
	BVAL_OR) VALUES (%d, :vals, :lats, :lngs, %s.nextval)"

// DML to insert Mmax data
#define DML_MMAX_DATA "INSERT INTO MMAX (MMAX_ID, MMAX_VAL, LAT, LON, \
	MMAX_OR) VALUES (%d, :vals, :lats, :lngs, %s.nextval)"


typedef struct _NSHM_Grid {
	double lat_min;
	double lat_max;
	double lat_inc;
	double lng_min;
	double lng_max;
	double lng_inc;
	int    grid_id;
	char   *grid_name;
	double *grid_values;
} NSHM_Grid;

int nshm_get_agrid(NSHM_Grid *_grid, char *_name);
int nshm_get_bgrid(NSHM_Grid *_grid, char *_name);
int nshm_get_mmax(NSHM_Grid *_grid, char *_name);

int nshm_put_agrid(NSHM_Grid *_grid);
int nshm_put_bgrid(NSHM_Grid *_grid);
int nshm_put_mmax(NSHM_Grid *_grid);

void nshm_print_grid(NSHM_Grid *_grid);
int  nshm_num_rows(NSHM_Grid *_grid);
void nshm_free_grid(NSHM_Grid *_grid);

#endif
