#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

#include "oci.h"

#include "adhoc_api.h"
#include "adhoc_utils.h"

void get_timestamp(char ** stamp) {
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

	asprintf(stamp, "%s %s %02d, %04d %02d:%02d:%02d", 
		WEEK_DAYS[lt->tm_wday],
		MONTHS[lt->tm_mon],
		lt->tm_mday,
		lt->tm_year+1900,
		lt->tm_hour,
		lt->tm_min,
		lt->tm_sec
	);
}

int adhoc_oci_check(void *_handle_pointer, int _handle_type, int _status) {
	text errbuf[512];
	int errcode = ADHOC_RETURN_SUCCESS;
	char * stamp; get_timestamp(&stamp);

	switch(_status) {
		case OCI_SUCCESS:
			break;
		case OCI_SUCCESS_WITH_INFO:
			OCIErrorGet(
				(void *) _handle_pointer, (ub4) 1, (text *) NULL,
				&errcode, errbuf, (ub4) sizeof(errbuf), _handle_type
			);
			fprintf(stderr, "Warning: %s [%d] -- %s", stamp, errcode, errbuf);
			errcode = ADHOC_RETURN_WARNING;
			break;
		case OCI_NEED_DATA:
			fprintf(stderr, "Error: %s [%d] -- Need Data\n",
				stamp, OCI_NEED_DATA);
			errcode = ADHOC_RETURN_ERROR;
			break;
		case OCI_NO_DATA:
			fprintf(stderr, "Warning: %s [%d] -- No Data\n",
				stamp, OCI_NO_DATA);
			errcode = ADHOC_RETURN_ERROR;
			break;
		case OCI_ERROR:
			OCIErrorGet(
				(void *) _handle_pointer, (ub4) 1, (text *) NULL,
				&errcode, errbuf, (ub4) sizeof(errbuf), _handle_type
			);
			fprintf(stderr, "Error: %s [%d] -- %s", stamp, errcode, errbuf);
			errcode = ADHOC_RETURN_ERROR;
			break;
		case OCI_INVALID_HANDLE:
			fprintf(stderr, "Error: %s [%d] -- Invalid Handle\n",
				stamp, OCI_INVALID_HANDLE);
			errcode = ADHOC_RETURN_ERROR;
			break;
		case OCI_STILL_EXECUTING:
			fprintf(stderr, "Error: %s [%d] -- Still Executing\n",
				stamp, OCI_STILL_EXECUTING);
			errcode = ADHOC_RETURN_ERROR;
			break;
		case OCI_CONTINUE:
			fprintf(stderr, "Error: %s [%d] -- Continue\n",
				stamp, OCI_CONTINUE);
			errcode = ADHOC_RETURN_ERROR;
			break;
		default:
			fprintf(stderr, "Error: %s [%d] -- Unknown OCI Return Status",
				stamp, _status
			);
			errcode = ADHOC_RETURN_ERROR;
			break;
	}

	return errcode;
}

#define addhoc_check_error(_error_hp, _status) adhoc_oci_check((_error_hp), OCI_HTYPE_ERROR, (_status))

#define adhoc_check_env(_env_hp, _status) adhoc_oci_check((_env_hp), OCI_HTYPE_ENV, (_status))
