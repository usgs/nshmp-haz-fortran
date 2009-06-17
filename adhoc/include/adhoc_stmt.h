#ifndef ADHOC_STMT_H
#define ADHOC_STMT_H

#include "adhoc_connection.h"

/**
 * These are the curently supported field types that we can bind to or select
 * from in the database. Other field types can be implemented as the need
 * arises.
 */
#define ADHOC_TYPE_INT 0
#define ADHOC_TYPE_FLT 1
#define ADHOC_TYPE_STR 2

/* The initial capacity of in/out fields to bind to a statement */
#define ADHOC_NUM_FIELDS_DEFAULT 8
/* The amount to increate the field array laength for a statement when the
 * current capacity is used */
#define ADHOC_NUM_FIELDS_INC 8

typedef struct _AdHoc_Field {
	void           *value;     /* The input/output buffer for the field */
	int            *val_lens;  /* Length of value buffer(s) */
	short          indicator;  /* Used when fetching strings. Max col width */
	unsigned short *ret_lens;  /* The column level data length from a defn */
	unsigned short *ret_codes; /* The column level return codes from a bind */
	int            type;       /* The AdHoc data type to bind or define */
	char           *name;      /* The name/index of the field/parameter */
	void           *handle;    /* The OCIBind/OCIDefn pointer for this field */
} AdHoc_Field;

typedef struct _AdHoc_Stmt {
	char        *sql;       /* The SQL string for this statement */
	OCIStmt     *stmthp;    /* The OCIStmt object */
	int         num_fields; /* How many does this statement currently hold */
	int         max_fields; /* How many can this statement currently hold */
	AdHoc_Field *fields;    /* The fields for this statement */
} AdHoc_Stmt;

AdHoc_Stmt adhoc_stmt_create(char * _sql, AdHoc_Field _fields[], int _num,
		AdHoc_Connection * _conn);

int adhoc_stmt_bind(AdHoc_Stmt * _stmt, AdHoc_Field * _field,
        AdHoc_Connection * _conn);

int adhoc_stmt_define(AdHoc_Stmt * _stmt, AdHoc_Field * _field,
        AdHoc_Connection * _conn);

int adhoc_stmt_execute(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn);

int adhoc_stmt_fetch(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn,
        int * _rows);

int adhoc_stmt_free(AdHoc_Stmt * _stmt, AdHoc_Connection * _conn);

#endif
