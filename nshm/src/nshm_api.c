#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "nshm_api.h"
#include "ocilib.h"

OCI_Connection *CONNECTION = NULL;

//------------------------------------------------------------------------------
// Private Function Prototypes
//------------------------------------------------------------------------------

void nshm_error_handler(OCI_Error *_error);

//------------------------------------------------------------------------------
// FORTRAN API Function Prototypes
//------------------------------------------------------------------------------

void nshm_initialize_(char *_user, char *_pass, char *_tns, int ulen, int plen,
		int tlen);

void nshm_cleanup_();

//------------------------------------------------------------------------------
// Public API Function Implementations
//------------------------------------------------------------------------------

int nshm_initialize(char *_user, char *_pass, char *_tns) {
	// Check if we already have an active connection
	if ( CONNECTION != NULL ) {
		char *message;
		char *stamp = nshm_timestamp();
		
		// Generate a log message
		asprintf(&message, "INFO -- [%s] -- %05i\n  %s\n", stamp,
			NSHM_ERR_CONN_EXISTS_CODE, NSHM_ERR_CONN_EXISTS_MSG);

		// Log the generated message
		nshm_log_error(message);

		// Free resources
		free(message);
		free(stamp);

		return NSHM_RETURN_WARNING;
	} // END: if (CONNECTION==NULL)

	// Attempt to initialize the environment
	if (!OCI_Initialize(nshm_error_handler, NULL, OCI_ENV_DEFAULT))
		return NSHM_RETURN_ERROR;


	// Establish the connection
	CONNECTION = OCI_ConnectionCreate(_tns, _user, _pass, OCI_SESSION_DEFAULT);

	return NSHM_RETURN_SUCCESS;
}

int nshm_cleanup() {
	OCI_ConnectionFree(CONNECTION);
	OCI_Cleanup();
	CONNECTION = NULL;
	return NSHM_RETURN_SUCCESS;
}

char * nshm_timestamp() {
	char * stamp;
	char * WEEK_DAYS[] = {"Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"};
   	char * MONTHS[] = {
   		"Jan", "Feb", "Mar", "Apr",
   		"May", "Jun", "Jul", "Aug",
   		"Sep", "Oct", "Nov", "Dec"
   	};
	
   	time_t t;
   	struct tm *lt;
	
   	t = time(NULL);
   	lt = localtime(&t);
	
    asprintf(&stamp, "%s %s %02d, %04d %02d:%02d:%02d",
   		WEEK_DAYS[lt->tm_wday],
   		MONTHS[lt->tm_mon],
   		lt->tm_mday,
   		lt->tm_year+1900,
   		lt->tm_hour,
   		lt->tm_min,
   		lt->tm_sec
   	);	

	return stamp;
}

void nshm_log_error(char *_message) {
	fprintf(stderr, _message);
}


//------------------------------------------------------------------------------
// FORTRAN API Function Implementations
//------------------------------------------------------------------------------

void nshm_initialize_(char *_user, char *_pass, char *_tns, int ulen, int plen,
		int tlen) {
	// FORTRAN does not use null-terminated strings. Copy the FORTRAN strings
	// using their indicator variables into proper null-terminated C-strings.
	char *user = calloc(ulen+1, sizeof(char));
	char *pass = calloc(plen+1, sizeof(char));
	char *tns  = calloc(tlen+1, sizeof(char));

	strncpy(user, _user, ulen);
	strncpy(pass, _pass, plen);
	strncpy(tns, _tns, tlen);
	
	// Now we can call our C routine
    nshm_initialize(user, pass, tns);
}

void nshm_cleanup_() { nshm_cleanup(); }

//------------------------------------------------------------------------------
// Private Function Implementations
//------------------------------------------------------------------------------

void nshm_error_handler(OCI_Error *_error) {
	char *message;
	char *stamp = nshm_timestamp();
	
	
	// Generate a log message
	asprintf(&message, "ERROR -- [%s] -- ORA-%05i\n  %s\n  %s\n",
		stamp,
		OCI_ErrorGetOCICode(_error),
		OCI_ErrorGetString(_error),
		OCI_GetSql(OCI_ErrorGetStatement(_error))
	);
		
	// Log the generated message
	nshm_log_error(message);
	
	// Free resources
	free(message);
	free(stamp);
}

