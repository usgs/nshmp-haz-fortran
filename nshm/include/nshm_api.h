#ifndef NSHM_API_H
#define NSHM_API_H

#include "ocilib.h"

// NSHM Return Codes
#define NSHM_RETURN_SUCCESS 0
#define NSHM_RETURN_WARNING 1
#define NSHM_RETURN_ERROR   2

// NSHM Error Codes
#define NSHM_ERR_CONN_EXISTS_CODE 1981

// NSHM Error Messages
#define NSHM_ERR_CONN_EXISTS_MSG "The connection already exists."


extern OCI_Connection *CONNECTION;

int nshm_initialize(char *_user, char *_pass, char *_tns);

int nshm_cleanup();

char * nshm_timestamp();

void nshm_log_error(char *_message);

#endif
