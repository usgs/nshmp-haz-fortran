#include <stdlib.h>
#include <string.h>

#include "oci.h"
#include "adhoc_connection.h"
#include "adhoc_api.h"
#include "adhoc_utils.h"

int adhoc_connect(char * _username, char * _password, char * _tnsalias,
		AdHoc_Connection * _connection) {
	AdHoc_Connection conn;
	conn.connected = FALSE;
	int status;

	status = adhoc_check_env(conn.envhp, OCIEnvCreate(
		&conn.envhp, OCI_DEFAULT, NULL, NULL, NULL, NULL,
		(size_t) 0, (void **) NULL)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }

	status = adhoc_check_error(conn.errhp, OCIHandleAlloc(
		conn.envhp, (void **) &conn.errhp, OCI_HTYPE_ERROR,
		(size_t) 0, (void **) NULL)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }

	status = adhoc_check_error(conn.errhp, OCIHandleAlloc(
		conn.envhp, (void **) &conn.svchp, OCI_HTYPE_SVCCTX,
		(size_t) 0, (void **) NULL)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }

	status = adhoc_check_error(conn.errhp, OCILogon2(
		conn.envhp, conn.errhp, &conn.svchp,
		(text *) _username, (ub4) strlen(_username),
		(text *) _password, (ub4) strlen(_password),
		(text *) _tnsalias, (ub4) strlen(_tnsalias),
		OCI_DEFAULT)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }
	
	if (status == ADHOC_RETURN_SUCCESS) {
		conn.username = calloc(1, strlen(_username) + 1);
		conn.password = calloc(1, strlen(_username) + 1);
		conn.tnsalias = calloc(1, strlen(_tnsalias) + 1);

		strncpy(conn.username, _username, strlen(_username));
		strncpy(conn.password, _password, strlen(_password));
		strncpy(conn.tnsalias, _tnsalias, strlen(_tnsalias));

		conn.connected = TRUE;
	}

	memmove(_connection, &conn, sizeof(AdHoc_Connection));

	return status;
}

int adhoc_close_connection(AdHoc_Connection * _connection) {
	int status;

	status = adhoc_check_error(_connection->errhp, OCILogoff(
		_connection->svchp, _connection->errhp)
	);

	_connection->connected = FALSE;

	return status;
}

boolean adhoc_is_connected(AdHoc_Connection _connection) {
	return _connection.connected;
}
