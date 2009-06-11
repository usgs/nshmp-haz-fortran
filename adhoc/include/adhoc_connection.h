#ifndef ADHOC_CONNECTION_H
#define ADHOC_CONNECTION_H

#include "oci.h"

typedef struct _AdHoc_Connection {
	int       connection_id;
	char      *username;
	char      *password;
	char      *tnsalias;
	OCIEnv    *envhp;
	OCIError  *errhp;
	OCISvcCtx *svchp;
	boolean   connected;
} AdHoc_Connection;

extern AdHoc_Connection ADHOC_DEFAULT_CONNECTION;

int adhoc_connect(char * _username, char * _password, char * _tnsalias,
		AdHoc_Connection * _connection);

int adhoc_close_connection(AdHoc_Connection * _connection);

boolean adhoc_is_connected(AdHoc_Connection _connection);
#endif
