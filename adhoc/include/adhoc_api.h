#ifndef ADHOC_API_H
#define ADHOC_API_H

#define ADHOC_RETURN_ERROR 2
#define ADHOC_RETURN_WARNING 1
#define ADHOC_RETURN_SUCCESS 0

#include "oci.h"

/* These handles are a stop-gap until full abstraction is accomplished. */
extern OCIEnv    *adhoc_env_hp; /* The OCI Environment handle */
extern OCIError  *adhoc_err_hp; /* The OCI Error handle */
extern OCISvcCtx *adhoc_svc_hp; /* The OCI Service Context handle */

/**
 * This structure contains the authentication information to access the
 * database. Initial implementation requires the use of the TNSNAMES type of
 * Oracle connection however this could be expanded to connect in different ways
 * (i.e. EZCONNECT) etc...
 *
 * 05/07/09 -- EMM: Initial implementation.
 */
typedef struct _AdHoc_AuthInfo {
	char * sid;         /* The schema SID */
	char * username;    /* The database user name */
	char * password;    /* The database password for the user */
	char * tns_alias;   /* The TNSNAMES alias for the connection info */
	char * ezconnect;   /* The EZCONNECT string for the connection info */
} AdHoc_AuthInfo;

/**
 * This structure contains all the required handles to connect and interact with
 * the database server through OCI. There will probably be one of these passed
 * around the program runtime or a single global instance unless we choose to
 * implement multithreading etc...
 */
typedef struct _AdHoc_Handles {
	OCIEnv    *env;   /* the OCI environment handle */
	OCIError  *err;   /* The OCI error handle  */
	OCISvcCtx *svc;   /* The OCI service context */
	int       status; /* The adhoc return status of most recent call */
} AdHoc_Handles;

/**
 * This structure represents a simple cursor or interaction with the database.
 * Cursors are used to interact with the database.
 */
typedef struct _AdHoc_Cursor {
	AdHoc_Handles *handles; /* Environment handles for this cursor */
	OCIStmt       *stmt;    /* The OCI Statement for this cursor */
	int           status;   /* The adhoc return status of most recent call */
} AdHoc_Cursor;

/**
 * This structure represents the logical idea of a column-field in the database.
 * each statement can contaim several fields. A field can either be an input
 * (BIND) field, or an output (DEFINE) field. The type of field is indicated by
 * the naming convention for the field's "field" (name) value. BIND fields are
 * named with the ":fieldname" syntax while define fields are named with the
 * "Nfieldname" where "N" is the Nth output column and the column to which you
 * wish to define.
 */
typedef struct _AdHoc_DbField {
	void *value;                /* The value to bind or define */
	short indicator;            /* Used when fetching strings. Max col width */
	unsigned short * ret_lens;  /* The column level data length from a defn */
	unsigned short * ret_codes; /* The column level return code from a bind */
	int   type;                 /* The AdHoc data type to bind or define */
	char *field;                /* The name of the field or parameter */
	void *bind;                 /* The OCIBind pointer for the field */
} AdHoc_DbField;

/**
 * This structure represents the logical idea of an SQL statement. The statement
 * contains a number of fields along with the SQL string and resulting OCIStmt
 * object. This statement can then be executed.
 */
typedef struct _AdHoc_Statement {
	AdHoc_Cursor  cursor;      /* The Cursor information (connection info) */
	char          *sql;        /* The SQL Query to send to the Database */
	int           num_fields;  /* The number of fields in this statement */
	AdHoc_DbField *fields;     /* Fields associated with this statement */
	boolean       initialized; /* Indecates if the statement is ready */
} AdHoc_Statement;



/**
 * Initializes the connection to the database. Creates and associates all
 * required handles. This method must be called before one can interact with the
 * database.
 *
 * AdHoc_AuthInfo _auth_info (IN) -- The DB authentication credentials.
 *
 */
int adhoc_init(AdHoc_AuthInfo _auth_info);

/**
 * Prepares the statement to execute. The given _stmt structure should be
 * partially initialized by the caller before calling this method. That is,
 * before calling this method, the caller should define the SQL and FIELDS of
 * the AdHoc_Statement. In this method call, the AdHoc_Cursor will be 
 * initialized and each of the using the SQL attribute to create an OCIStmt and
 * then binding/defining each of the AdHoc_DbField structures as appropriate.
 *
 * AdHoc_Statement _stmt (IN/OUT) -- The partially initialized AdHoc statement
 *                                   structure.
 */
int adhoc_statement_prepare(AdHoc_Statement * _stmt);

/**
 * Executes a statement agaist the database. Statements executed in this manner
 * are auto-commited upon success.
 *
 * AdHoc_Statement _stmt (IN/OUT) -- The statement to execute. All handles
 * should be bound/defined before calling and after calling the defined buffers
 * will hold the selected data (in the case of a select statement).
 */
int adhoc_query(AdHoc_Statement * _stmt);

/**
 * Closes the connection to the database and frees the memory allocated for the
 * corresponding handles that are no longer used.
 */
int adhoc_close();

#endif
