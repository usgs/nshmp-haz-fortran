#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ocilib.h"

#include "nshm_api.h"
#include "nshm_grid.h"

#define LOWER_BOUND_IDX 0
#define UPPER_BOUND_IDX 1

//------------------------------------------------------------------------------
// Private Function Prototypes
//------------------------------------------------------------------------------

int nshm_get_grid_id(char *_stub, char *_name, int *_id);
int nshm_get_grid_meta(NSHM_Grid *_grid, char *_stub, int _id);
int nshm_get_grid_data(NSHM_Grid *_grid, char *_stub, int _id);
int nshm_put_grid_meta(NSHM_Grid *_grid, char *_stub, char *_seq);
int nshm_put_grid_data(NSHM_Grid *_grid, char *_stub, char *_seq);

//------------------------------------------------------------------------------
// FORTRAN API Function Prototypes
//------------------------------------------------------------------------------

/* Meta getter functions */
void nshm_get_agrid_meta_(NSHM_Grid *_grid, char * _name, int _nlen);
void nshm_get_bgrid_meta_(NSHM_Grid *_grid, char * _name, int _nlen);
void nshm_get_mmax_meta_(NSHM_Grid *_grid, char * _name, int _nlen);

/* Data getter functions */
void nshm_get_agrid_data_(NSHM_Grid *_grid, char * _name, int _nlen);
void nshm_get_bgrid_data_(NSHM_Grid *_grid, char * _name, int _nlen);
void nshm_get_mmax_data_(NSHM_Grid *_grid, char * _name, int _nlen);

/* Meta and Data setter functions */
void nshm_put_agrid_(NSHM_Grid *_grid, int nlen);
void nshm_put_bgrid_(NSHM_Grid *_grid, int nlen);
void nshm_put_mmax_(NSHM_Grid *_grid, int nlen);

/* Utility functions */
void nshm_print_grid_(NSHM_Grid *_grid);
void nshm_num_rows_(NSHM_Grid *_grid, int *_num_rows);
void nshm_free_grid_(NSHM_Grid *_grid);

//------------------------------------------------------------------------------
// Public API Function Implementations
//------------------------------------------------------------------------------
int nshm_get_agrid(NSHM_Grid *_grid, char *_name) {
	
	int status;
	int id;

	// Fetch the grid id
	nshm_get_grid_id(QUERY_AGRID_ID, _name, &id);

	// Fetch the metadata for this agrid
	status = nshm_get_grid_meta(_grid, QUERY_AGRID_META, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Agrid metadata.");
		return status;
	}

	// Allocate space for the grid values
	_grid->grid_values = calloc(nshm_num_rows(_grid), sizeof(double));

	// Fetch the data for this agrid
	status = nshm_get_grid_data(_grid, QUERY_AGRID_DATA, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Agrid data.");
		return status;
	}
	
	return NSHM_RETURN_SUCCESS;
}

int nshm_get_bgrid(NSHM_Grid *_grid, char *_name) {
	
	int status;
	int id;

	// Fetch the grid id
	nshm_get_grid_id(QUERY_BGRID_ID, _name, &id);

	// Fetch the metadata for this agrid
	status = nshm_get_grid_meta(_grid, QUERY_BGRID_META, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Bgrid metadata.");
		return status;
	}

	// Allocate space for the grid values
	_grid->grid_values = calloc(nshm_num_rows(_grid), sizeof(double));

	// Fetch the data for this agrid
	status = nshm_get_grid_data(_grid, QUERY_BGRID_DATA, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Bgrid data.");
		return status;
	}
	
	return NSHM_RETURN_SUCCESS;
}

int nshm_get_mmax(NSHM_Grid *_grid, char *_name) {
	
	int status;
	int id;

	// Fetch the grid id
	nshm_get_grid_id(QUERY_MMAX_ID, _name, &id);

	// Fetch the metadata for this agrid
	status = nshm_get_grid_meta(_grid, QUERY_MMAX_META, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Mmax metadata.");
		return status;
	}

	// Allocate space for the grid values
	_grid->grid_values = calloc(nshm_num_rows(_grid), sizeof(double));

	// Fetch the data for this agrid
	status = nshm_get_grid_data(_grid, QUERY_MMAX_DATA, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Mmax data.");
		return status;
	}
	
	return NSHM_RETURN_SUCCESS;
}

int nshm_put_agrid(NSHM_Grid *_grid) {
	int status;

	// Put the meta information in the database and fetch the generated id
	status = nshm_put_grid_meta(_grid, DML_AGRID_META, AGRID_SEQ);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to insert agrid metadata.");
		return status;
	}

	// Put the data into the agrid table in the database
	status = nshm_put_grid_data(_grid, DML_AGRID_DATA, AGRID_OR_SEQ);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to insert agrid data.");
		return status;
	}

	return NSHM_RETURN_SUCCESS;
}

int nshm_put_bgrid(NSHM_Grid *_grid) {
	int status;

	// Put the meta information in the database and fetch the generated id
	status = nshm_put_grid_meta(_grid, DML_BGRID_META, BGRID_SEQ);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to insert bgrid metadata.");
		return status;
	}

	// Put the data into the bgrid table in the database
	status = nshm_put_grid_data(_grid, DML_BGRID_DATA, BGRID_OR_SEQ);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to insert bgrid data.");
		return status;
	}

	return NSHM_RETURN_SUCCESS;
}

int nshm_put_mmax(NSHM_Grid *_grid) {
	int status;

	// Put the meta information in the database and fetch the generated id
	status = nshm_put_grid_meta(_grid, DML_MMAX_META, MMAX_SEQ);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to insert mmax metadata.");
		return status;
	}

	// Put the data into the mmax table in the database
	status = nshm_put_grid_data(_grid, DML_MMAX_DATA, MMAX_OR_SEQ);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to insert mmax data.");
		return status;
	}

	return NSHM_RETURN_SUCCESS;
}

void nshm_print_grid(NSHM_Grid *_grid) {
	int i = 0;
	int num_rows = nshm_num_rows(_grid);
	double lat; double lng;

	// Print the metadata
	printf("Metadata:\n   ID: %d\n   Name: %s\n", 
		_grid->grid_id, _grid->grid_name);
	printf("   Latitude: [ %f, %f : %f ]\n",
		_grid->lat_min, _grid->lat_max, _grid->lat_inc);
	printf("   Longitude: [ %f, %f : %f ]\n",
		_grid->lng_min, _grid->lng_max, _grid->lng_inc);

	// Print the data
	lat = _grid->lat_max; lng = _grid->lng_min;

	printf("Data:\n  Latitude   Longitude    Value\n");

	if ( num_rows > 10 ) {
		// Just print a snippet of data
		i = 0; // Counter of how many iterations we've printed
		for (lat=_grid->lat_max; lat>=_grid->lat_min; lat-=_grid->lat_inc) {
			for (lng=_grid->lng_min; lng<=_grid->lng_max; lng+=_grid->lng_inc) {
				printf("   %4.2f      %5.2f    %7.5f\n", lat, lng,
					_grid->grid_values[i++]);
					if ( i > 4 ) { break; } // Print only 5 rows
			}
			if ( i > 4 ) { break; } // Print only 5 rows
		}

		printf("     .            .       .\n");
		printf("     .            .       .\n");
		printf("     .            .       .\n");
		printf("  %d more rows of data.\n", num_rows - 5);
		
	} else {
		// Just print all of the data. There isn't that much of it.
		for (lat=_grid->lat_max; lat>=_grid->lat_min; lat-=_grid->lat_inc) {
			for (lng=_grid->lng_min; lng<=_grid->lng_max; lng+=_grid->lng_inc) {
				printf("   %4.2f      %5.2f    %7.5f\n", lat, lng,
					_grid->grid_values[i++]);
			}
		}
	}
}

int nshm_num_rows(NSHM_Grid *_grid) {
	double lats = (_grid->lat_max-_grid->lat_min)/_grid->lat_inc+1;
	double lngs = (_grid->lng_max-_grid->lng_min)/_grid->lng_inc+1;
	int num_rows = rint(lats) * rint(lngs);
	return num_rows;
}

void nshm_free_grid(NSHM_Grid *_grid) {
	// Free allocated resources
	free(_grid->grid_name);
	free(_grid->grid_values);

	// Set values to NULL
	_grid->lat_min = 0.0;
	_grid->lat_max = 0.0;
	_grid->lat_inc = 0.0;
	_grid->lng_min = 0.0;
	_grid->lng_max = 0.0;
	_grid->lng_inc = 0.0;
	_grid->grid_id = 0;
}
//------------------------------------------------------------------------------
// FORTRAN API Function Implementations
//------------------------------------------------------------------------------

void nshm_get_agrid_meta_(NSHM_Grid *_grid, char * _name, int _nlen) {
	int id;
	char *name = calloc(_nlen+1, sizeof(char));
	strncpy(name, _name, _nlen);
	nshm_get_grid_id(QUERY_AGRID_ID, name, &id);
	nshm_get_grid_meta(_grid, QUERY_AGRID_META, id);
}   
    
void nshm_get_bgrid_meta_(NSHM_Grid *_grid, char * _name, int _nlen) {
	int id;
	char *name = calloc(_nlen+1, sizeof(char));
	strncpy(name, _name, _nlen);
	nshm_get_grid_id(QUERY_BGRID_ID, name, &id);
	nshm_get_grid_meta(_grid, QUERY_BGRID_META, id);
}   
    
void nshm_get_mmax_meta_(NSHM_Grid *_grid, char * _name, int _nlen) {
	int id;
	char *name = calloc(_nlen+1, sizeof(char));
	strncpy(name, _name, _nlen);
	nshm_get_grid_id(QUERY_MMAX_ID, name, &id);
	nshm_get_grid_meta(_grid, QUERY_MMAX_META, id);
}   
    

void nshm_get_agrid_data_(NSHM_Grid *_grid, char * _name, int _nlen) {
	int id;
	char *name = calloc(_nlen+1, sizeof(char));
	strncpy(name, _name, _nlen);
	nshm_get_grid_id(QUERY_AGRID_ID, name, &id);
	nshm_get_grid_data(_grid, QUERY_AGRID_DATA, id);
}   

void nshm_get_bgrid_data_(NSHM_Grid *_grid, char * _name, int _nlen) {
	int id;
	char *name = calloc(_nlen+1, sizeof(char));
	strncpy(name, _name, _nlen);
	nshm_get_grid_id(QUERY_BGRID_ID, name, &id);
	nshm_get_grid_data(_grid, QUERY_BGRID_DATA, id);
}   

void nshm_get_mmax_data_(NSHM_Grid *_grid, char *_name, int _nlen) {
	int id;
	char *name = calloc(_nlen+1, sizeof(char));
	strncpy(name, _name, _nlen);
	nshm_get_grid_id(QUERY_MMAX_ID, name, &id);
	nshm_get_grid_data(_grid, QUERY_MMAX_DATA, id);
}   


void nshm_put_agrid_(NSHM_Grid *_grid, int nlen) {
	// Copy the name to a null-terminated string
	char *name = calloc(nlen+1, sizeof(char));
	strncpy(name, _grid->grid_name, nlen);

	// Free the original name
	free(_grid->grid_name);

	// Copy the null-terminated name back to the grid_name
	_grid->grid_name = calloc(nlen+1, sizeof(char));
	strncpy(_grid->grid_name, name, nlen);

	// Call the C-method to actually put the grid in the db
	nshm_put_agrid(_grid);
}

void nshm_put_bgrid_(NSHM_Grid *_grid, int nlen) {
	// Copy the name to a null-terminated string
	char *name = calloc(nlen+1, sizeof(char));
	strncpy(name, _grid->grid_name, nlen);

	// Free the original name
	free(_grid->grid_name);

	// Copy the null-terminated name back to the grid_name
	_grid->grid_name = calloc(nlen+1, sizeof(char));
	strncpy(_grid->grid_name, name, nlen);

	// Call the C-method to actually put the grid in the db
	nshm_put_bgrid(_grid);
}

void nshm_put_mmax_(NSHM_Grid *_grid, int nlen) {
	// Copy the name to a null-terminated string
	char *name = calloc(nlen+1, sizeof(char));
	strncpy(name, _grid->grid_name, nlen);

	// Free the original name
	free(_grid->grid_name);

	// Copy the null-terminated name back to the grid_name
	_grid->grid_name = calloc(nlen+1, sizeof(char));
	strncpy(_grid->grid_name, name, nlen);

	// Call the C-method to actually put the grid in the db
	nshm_put_mmax(_grid);
}

void nshm_print_grid_(NSHM_Grid *_grid) { nshm_print_grid(_grid); }

void nshm_num_rows_(NSHM_Grid *_grid, int *_num_rows) {
	*_num_rows = nshm_num_rows(_grid);
}

void nshm_free_grid_(NSHM_Grid *_grid) { nshm_free_grid(_grid); }

//------------------------------------------------------------------------------
// Private Function Implementations
//------------------------------------------------------------------------------

int nshm_get_grid_id(char *_stub, char *_name, int *_id) {
	char *query;
	OCI_Resultset *rs;
	OCI_Statement *st = OCI_StatementCreate(CONNECTION);

	// Build the query
	asprintf(&query, _stub, _name);

	// Execute the query
	OCI_ExecuteStmt(st, query);
	rs = OCI_GetResultset(st);

	// Fetch the results
	if (!OCI_FetchNext(rs)) { return NSHM_RETURN_ERROR; }
	*_id = OCI_GetInt(rs, GRID_ID_IDX);

	// Clean up our resources and return
	free(query);
	OCI_StatementFree(st);
	return NSHM_RETURN_SUCCESS;
}

int nshm_get_grid_meta(NSHM_Grid *_grid, char *_stub, int _id) {
	char *query;
	int name_len;
	double precision;
	OCI_Resultset *rs;
	OCI_Statement *st = OCI_StatementCreate(CONNECTION);

	// How many decimal places of precision to save.
	precision = pow(10, NSHM_GRID_PRECISION);

	// Build the query string
	asprintf(&query, _stub, _id);

	// Execute the query
	OCI_ExecuteStmt(st, query);
	rs = OCI_GetResultset(st);
	
	// Fetch the information
	if (!OCI_FetchNext(rs)) { return NSHM_RETURN_ERROR; }

	// Assign the grid information from our resultset
	_grid->grid_id = OCI_GetInt(rs, GRID_ID_IDX);
	name_len = strlen(OCI_GetString(rs, GRID_NAME_IDX));
	_grid->grid_name = calloc(name_len+1, sizeof(char));
	strncpy(_grid->grid_name, OCI_GetString(rs, GRID_NAME_IDX), name_len);

	_grid->lat_min = OCI_GetDouble(rs, GRID_MIN_LAT_IDX);
	_grid->lat_max = OCI_GetDouble(rs, GRID_MAX_LAT_IDX);
	_grid->lat_inc = OCI_GetDouble(rs, GRID_INC_LAT_IDX);

	_grid->lng_min = OCI_GetDouble(rs, GRID_MIN_LNG_IDX);
	_grid->lng_max = OCI_GetDouble(rs, GRID_MAX_LNG_IDX);
	_grid->lng_inc = OCI_GetDouble(rs, GRID_INC_LNG_IDX);

	free(query);
	OCI_StatementFree(st);
	return NSHM_RETURN_SUCCESS;
}

int nshm_get_grid_data(NSHM_Grid *_grid, char *_stub, int _id) {
	char *query;
	int i, num_rows = nshm_num_rows(_grid);
	OCI_Resultset *rs;
	OCI_Statement *st = OCI_StatementCreate(CONNECTION);

	// Computed the expected number of rows based on the metadata
	// 08/25/09 -- EMM: Do not do this anymore. Memory is allocated in Fortran
	//                  or as part of the public C-interface specification.
	//_grid->grid_values = calloc(num_rows, sizeof(double));

	// Finish out the query with the give _id
	asprintf(&query, _stub, _id);

	// Execute the query
	OCI_ExecuteStmt(st, query);
	rs = OCI_GetResultset(st);

	// Fetch the information into our values array
	for (i=0; i < num_rows && OCI_FetchNext(rs); i++) {
		_grid->grid_values[i] = OCI_GetDouble(rs, GRID_VAL_IDX);
	}

	free(query);
	OCI_StatementFree(st);
	return NSHM_RETURN_SUCCESS;
}

int nshm_put_grid_meta(NSHM_Grid *_grid, char *_stub, char *_seq) {
	char *query;
	OCI_Resultset *rs;
	OCI_Statement *st = OCI_StatementCreate(CONNECTION);

	// Prepare the statement and insert the row into the meta table.
	asprintf(&query, _stub, _seq);
	OCI_Prepare(st, query);
	OCI_BindString(st, ":name", _grid->grid_name, strlen(_grid->grid_name)+1);
	OCI_BindDouble(st, ":lat_min", &_grid->lat_min);
	OCI_BindDouble(st, ":lat_max", &_grid->lat_max);
	OCI_BindDouble(st, ":lat_inc", &_grid->lat_inc);
	OCI_BindDouble(st, ":lng_min", &_grid->lng_min);
	OCI_BindDouble(st, ":lng_max", &_grid->lng_max);
	OCI_BindDouble(st, ":lng_inc", &_grid->lng_inc);
	OCI_Execute(st);
	OCI_StatementFree(st);
	free(query);

	// Now check the sequence number to fetch the curr val to link the meta and
	// data rows together.
	st = OCI_StatementCreate(CONNECTION);
	asprintf(&query, "SELECT %s.currval FROM DUAL", _seq);
	OCI_ExecuteStmt(st, query);
	rs = OCI_GetResultset(st);
	free(query);

	// Fetch the value
	if ( OCI_FetchNext(rs) ) {
		// The grid->grid_id is now set to the inserted row id
		_grid->grid_id = OCI_GetInt(rs, 1);
	} else {
		nshm_log_error("Failed to fetch sequence iterator.");

		OCI_StatementFree(st);
		return NSHM_RETURN_ERROR;
	}

	OCI_Commit(CONNECTION); // Commit changes.
	OCI_StatementFree(st);  // Frees the statement and result set at once.

	return NSHM_RETURN_SUCCESS;
}

int nshm_put_grid_data(NSHM_Grid *_grid, char *_stub, char *_seq) {
	int idx = 0;
	int num_rows = nshm_num_rows(_grid);
	int i, j, latmin, latmax, latinc, lngmin, lngmax, lnginc; // Iterators
	double *latvals, *lngvals;
	double precision;
	char *error, *query;

	OCI_Statement *st;

	precision = pow(10.0, NSHM_GRID_PRECISION);

	// 07/16/09 -- EMM: Using ints for loop iteration  increments and bounds are
	// faster than doubles and also do not have the precision errors.
	latmin = (int) rint(precision * _grid->lat_min);
	latmax = (int) rint(precision * _grid->lat_max);
	latinc = (int) rint(precision * _grid->lat_inc);
	lngmin = (int) rint(precision * _grid->lng_min);
	lngmax = (int) rint(precision * _grid->lng_max);
	lnginc = (int) rint(precision * _grid->lng_inc);

	latvals = calloc(num_rows, sizeof(double));
	lngvals = calloc(num_rows, sizeof(double));

	// Build an array of latitude and longitude values
	for ( i = latmax; i >= latmin; i -= latinc ) {
		for ( j = lngmin; j <= lngmax; j += lnginc ) {
			latvals[idx] = ((double) i) / ((double) precision);
			lngvals[idx] = ((double) j) / ((double) precision);
			idx++;
		}
	}

	if ( idx != num_rows ) {
		asprintf(
			&error,
			"idx and num_rows differed!\n\tidx: %d\n\tnum_rows: %d\n\n",
			idx, num_rows
		);
		nshm_log_error(error);
		free(error); free(latvals); free(lngvals);
		return NSHM_RETURN_ERROR;
	}

	// Prepare the statement and bind our arrays
	asprintf(&query, _stub, _grid->grid_id, _seq);
	st = OCI_StatementCreate(CONNECTION);
	OCI_Prepare(st, query);
	OCI_BindArraySetSize(st, num_rows);
	OCI_BindArrayOfDoubles(st, ":lats", latvals, num_rows);
	OCI_BindArrayOfDoubles(st, ":lngs", lngvals, num_rows);
	OCI_BindArrayOfDoubles(st, ":vals", _grid->grid_values, num_rows);

	// Execute the statement and check if any errors occurred.
	if ( !OCI_Execute(st) ) {
		asprintf(
			&error,
			"There were %d errors executing query.\n",
			OCI_GetBatchErrorCount(st)
		);
		nshm_log_error(error);
		free(error); free(query); free(latvals); free(lngvals);
		return NSHM_RETURN_ERROR;
	}

	OCI_Commit(CONNECTION);

	// Clean up our memory and stuff.
	free(query); free(latvals); free(lngvals);
	OCI_StatementFree(st);

	return NSHM_RETURN_SUCCESS;
}
