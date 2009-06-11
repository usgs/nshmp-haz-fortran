#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "oci.h"

#include "adhoc_api.h"
#include "adhoc_stmt.h"
#include "adhoc_connection.h"
#include "adhoc_utils.h"

//------------------------------------------------------------------------------
// Private Method Prototypes (i.e. Published in the header file)
//------------------------------------------------------------------------------

int adhoc_associate_field(AdHoc_Stmt * _stmt, AdHoc_Field * _field);

ub2 adhoc_oci_type(int _type);

//------------------------------------------------------------------------------
// Public Methods Implementations (i.e. Published in the header file)
//------------------------------------------------------------------------------

/*
AdHoc_Stmt adhoc_stmt_create(char * _sql, AdHoc_Field _fields[], int _num) {
	AdHoc_Stmt stmt;

	// Copy the SQL string into place
	stmt.sql = calloc(1, strlen(_sql) + 1);
	strncpy(stmt.sql, _sql, strlen(_sql));


	if (_num > 0) {
		// Set the current field list capacity values
		if (_num > ADHOC_NUM_FIELDS_DEFAULT) {
			stmt.fields = calloc(_num, sizeof(AdHoc_Field));
			stmt.max_fields = _num;
		} else {
			stmt.fields = calloc(ADHOC_NUM_FIELDS_DEFAULT, sizeof(AdHoc_Field));
			stmt.max_fields = ADHOC_NUM_FIELDS_DEFAULT;
		}

		// Set the current number of used fields
		memmove(stmt.fields, _fields, sizeof(AdHoc_Field) * _num);
		stmt.num_fields = _num;

	} else {
		// Create the default allocation
		stmt.fields = calloc(ADHOC_NUM_FIELDS_DEFAULT, sizeof(AdHoc_Field));
		stmt.max_fields = ADHOC_NUM_FIELDS_DEFAULT;

		// Set the current number of fields
		stmt.num_fields = 0; 
	}

	return stmt;
}

int adhoc_stmt_prepare(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn) {
	int i;
	int status;

	// Create the OCI Statement
	status = adhoc_check_error(_conn->errhp, OCIStmtPrepare2(
		_conn->svchp, &_stmt->stmthp, _conn->errhp,
		(const OraText *) _stmt->sql, strlen((char *) _stmt->sql),
		(OraText *) NULL, 0, OCI_NTV_SYNTAX, OCI_DEFAULT)
	);

	// Bind/Define input/output fields.
	for (i = 0; i < _stmt->num_fields; ++i) {
		if (atoi(_stmt->fields[i].name) != 0) {
			// Numeric field name indicates an output field. Define by position.
			status = adhoc_stmt_define(_stmt, &_stmt->fields[i], _conn);
		} else {
			// ASCII field name indicates an input field. Bind by name.
			status = adhoc_stmt_bind(_stmt, &_stmt->fields[i], _conn);
		}
		if (status != ADHOC_RETURN_SUCCESS) { break; }
	}

	return status;
}
*/

AdHoc_Stmt adhoc_stmt_create(char * _sql, AdHoc_Field _fields[], int _num,
		AdHoc_Connection * _conn) {

	int status, i;
	char * stamp;
	AdHoc_Stmt stmt;

	// Copy the SQL string into the statement structure
	stmt.sql = calloc(1, strlen(_sql) + 1);
	strncpy(stmt.sql, _sql, strlen(_sql));

	// Prepare the OCI statement
	status = adhoc_check_error(_conn->errhp, OCIStmtPrepare2(
		_conn->svchp, &stmt.stmthp, _conn->errhp,
		(const OraText *) stmt.sql, strlen((char *) stmt.sql),
		(OraText *) NULL, 0, OCI_NTV_SYNTAX, OCI_DEFAULT)
	);

	// Allocate space for the fields
	stmt.max_fields = ( _num > ADHOC_NUM_FIELDS_DEFAULT ) ? 
			_num:ADHOC_NUM_FIELDS_DEFAULT;
	
	stmt.fields = calloc(stmt.max_fields, sizeof(AdHoc_Field));
	stmt.num_fields = 0;

	// Associate and bind each field given to us
	for (i = 0; status == ADHOC_RETURN_SUCCESS && i < _num; i++) {
		if (atoi(_fields[i].name) != 0) {
			// Numeric field name. Associate with a define by position.
			status = adhoc_stmt_define(&stmt, &_fields[i], _conn);
		} else {
			// ASCII field name. Associate with a bind by name.
			status = adhoc_stmt_bind(&stmt, &_fields[i], _conn);
		}
	}

	if (status != ADHOC_RETURN_SUCCESS) {
		get_timestamp(&stamp);
		fprintf(stderr, "Error: %s [%d] -- Statement creation failed.",
			stamp, status);
		adhoc_stmt_free(&stmt, _conn);
	}

	return stmt;
}

int adhoc_stmt_bind(AdHoc_Stmt * _stmt, AdHoc_Field * _field,
		AdHoc_Connection * _conn) {
	int status;
	int idx;

	idx = adhoc_associate_field(_stmt, _field);

	status = adhoc_check_error(_conn->errhp, OCIBindByName(
		_stmt->stmthp, (OCIBind **) &_stmt->fields[idx].handle, _conn->errhp,
		(text *) _stmt->fields[idx].name, (sb4) _stmt->fields[idx].val_lens,
		(void *) &_stmt->fields[idx].value,
		sizeof(_stmt->fields[idx].value), _stmt->fields[idx].type,
		(void *) NULL, (ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL, OCI_DEFAULT)
	);

	return status;
}

int adhoc_stmt_define(AdHoc_Stmt * _stmt, AdHoc_Field * _field,
		AdHoc_Connection * _conn) {
	int status;
	int idx;

	idx = adhoc_associate_field(_stmt, _field);

	status = adhoc_check_error(_conn->errhp, OCIDefineByPos(
		_stmt->stmthp, (OCIDefine **) &_stmt->fields[idx].handle, _conn->errhp,
		atoi(_stmt->fields[idx].name), (void *) &_stmt->fields[idx].value,
		(sb4) _stmt->fields[idx].val_lens, _stmt->fields[idx].type,
		(void *) NULL, (ub2 *) NULL, (ub2 *) NULL, OCI_DEFAULT)
	);

	return status;
}

int adhoc_stmt_execute(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn) {
	int status;

	status = adhoc_check_error(_conn->errhp, OCIStmtExecute(
		_conn->svchp, _stmt->stmthp, _conn->errhp, 0, 0,
		(OCISnapshot *) NULL, (OCISnapshot *) NULL, OCI_DEFAULT)
	);

	return status;	
}

int adhoc_stmt_fetch(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn, int * _rows){
	int status;

	status = OCIStmtFetch2(_stmt->stmthp, _conn->errhp, *_rows,
			OCI_DEFAULT, 0, OCI_DEFAULT);

	if (status == OCI_SUCCESS) {
		status = ADHOC_RETURN_SUCCESS;
	} else if (status == OCI_NO_DATA) {
		// Didn't get as many rows as we asked for
		status = ADHOC_RETURN_WARNING;
		// Find out how many we actually got
		adhoc_check_error(_conn->errhp, OCIAttrGet(
			_stmt->stmthp, OCI_HTYPE_STMT, _rows, (ub4 *) NULL,
			OCI_ATTR_ROWS_FETCHED, _conn->errhp)
		);
	} else {
		// Something went wrong
		status = adhoc_check_error(_conn->errhp, status);
	}

	return status;
}

int adhoc_stmt_free(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn) {
	OCIStmtRelease(_stmt->stmthp,_conn->errhp,(OraText *) NULL,0,OCI_DEFAULT);
	free(_stmt->fields);
	free(_stmt->sql);

	return ADHOC_RETURN_SUCCESS;
}

//------------------------------------------------------------------------------
// Private Methods Implementations (i.e. Not published in the header file)
//------------------------------------------------------------------------------

int adhoc_associate_field(AdHoc_Stmt * _stmt, AdHoc_Field * _field) {
	int idx;

	if (_stmt->num_fields >= _stmt->max_fields) {
		_stmt->max_fields += ADHOC_NUM_FIELDS_INC;
		realloc(_stmt->fields, _stmt->max_fields);
	}

	idx = _stmt->num_fields++;
	memmove(&_stmt->fields[idx], _field, sizeof(AdHoc_Field));

	return idx;
}

ub2 adhoc_oci_type(int  _type) {
	switch (_type) {
		case ADHOC_TYPE_INT:
			return (ub2) SQLT_INT;
			break;
		case ADHOC_TYPE_FLT:
			return (ub2) SQLT_FLT;
			break;
		case ADHOC_TYPE_STR:
			return (ub2) SQLT_STR;
			break;
		default:
			fprintf(stderr, "Unknown OCI type conversion (%d).\n", _type);
			break;
	}
	return (ub2) NULL;
}
