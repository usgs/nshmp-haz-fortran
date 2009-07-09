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

int nshm_get_grid_meta(NSHM_Grid *_grid, char *_stub, int _id);
int nshm_get_grid_data(NSHM_Grid *_grid, char *_stub, int _id);

//------------------------------------------------------------------------------
// FORTRAN API Function Prototypes
//------------------------------------------------------------------------------

void nshm_get_agrid_(NSHM_Grid *_grid);
void nshm_get_bgrid_(NSHM_Grid *_grid);
void nshm_get_mmax_(NSHM_Grid *_grid);

void nshm_print_grid_(NSHM_Grid *_grid);
void nshm_num_rows_(NSHM_Grid *_grid, int *_num_rows);
void nshm_free_grid_(NSHM_Grid *_grid);

//------------------------------------------------------------------------------
// Public API Function Implementations
//------------------------------------------------------------------------------
int nshm_get_agrid(NSHM_Grid *_grid) {
	
	int status;
	int id = 1; // Randomize this at some point

	// Fetch the metadata for this agrid
	status = nshm_get_grid_meta(_grid, QUERY_AGRID_META, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Agrid metadata.");
		return status;
	}

	// Fetch the data for this agrid
	status = nshm_get_grid_data(_grid, QUERY_AGRID_DATA, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Agrid data.");
		return status;
	}
	
	return NSHM_RETURN_SUCCESS;
}

int nshm_get_bgrid(NSHM_Grid *_grid) {
	
	int status;
	int id = 1; // There is only 1 bgrid. We always use this one.

	// Fetch the metadata for this agrid
	status = nshm_get_grid_meta(_grid, QUERY_BGRID_META, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Bgrid metadata.");
		return status;
	}

	// Fetch the data for this agrid
	status = nshm_get_grid_data(_grid, QUERY_BGRID_DATA, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Bgrid data.");
		return status;
	}
	
	return NSHM_RETURN_SUCCESS;
}

int nshm_get_mmax(NSHM_Grid *_grid) {
	
	int status;
	int id = 1; // Randomize this at some point

	// Fetch the metadata for this agrid
	status = nshm_get_grid_meta(_grid, QUERY_MMAX_META, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Mmax metadata.");
		return status;
	}

	// Fetch the data for this agrid
	status = nshm_get_grid_data(_grid, QUERY_MMAX_DATA, id);
	if ( status >= NSHM_RETURN_WARNING ) {
		nshm_log_error("Failed to get Mmax data.");
		return status;
	}
	
	return NSHM_RETURN_SUCCESS;
}

void nshm_print_grid(NSHM_Grid *_grid) {
	int i;
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
	if (num_rows > 10) {
		// Just print a snippet of data
		for (i = 0; i < 5; i++) {
			printf("   %4.2f      %5.2f    %7.5f\n", lat, lng,
				_grid->grid_values[i]);
			lng += _grid->lng_inc;
			if ( lng > _grid->lng_max ) {
				lng = _grid->lng_max;
				lat -= _grid->lat_inc;
			}
		}
		printf("     .            .       .\n");
		printf("     .            .       .\n");
		printf("     .            .       .\n");
		printf("  %d more rows of data.\n", num_rows - 5);
		
	} else {
		// Just print all of the data
		for (i = 0; i < num_rows; i++) {
			printf("   %4.2f   %5.2f   %7.5f\n",lat,lng,_grid->grid_values[i]);
			lng += _grid->lng_inc;
			if ( lng > _grid->lng_max ) {
				lng = _grid->lng_max;
				lat -= _grid->lat_inc;
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

void nshm_get_agrid_(NSHM_Grid *_grid) { nshm_get_agrid(_grid); }

void nshm_get_bgrid_(NSHM_Grid *_grid) { nshm_get_bgrid(_grid); }

void nshm_get_mmax_(NSHM_Grid *_grid) { nshm_get_mmax(_grid); }

void nshm_print_grid_(NSHM_Grid *_grid) { nshm_print_grid(_grid); }

void nshm_num_rows_(NSHM_Grid *_grid, int *_num_rows) {
	*_num_rows = nshm_num_rows(_grid);
}

void nshm_free_grid_(NSHM_Grid *_grid) { nshm_free_grid(_grid); }

//------------------------------------------------------------------------------
// Private Function Implementations
//------------------------------------------------------------------------------

int nshm_get_grid_meta(NSHM_Grid *_grid, char *_stub, int _id) {
	char *query;
	int name_len;
	OCI_Resultset *rs;
	OCI_Statement *st = OCI_StatementCreate(CONNECTION);

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
	_grid->grid_name = calloc(1, (sizeof(char) * name_len) + 1);
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
	_grid->grid_values = calloc(num_rows, sizeof(double));

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
