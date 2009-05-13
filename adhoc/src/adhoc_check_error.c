#include <stdio.h>

#include "oci.h"
#include "adhoc_api.h"
#include "adhoc_utils.h"
#include "adhoc_check_error.h"

int adhoc_oci_check(void *_handle_pointer, ub4 _handle_type, sword _status) {
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
			fprintf(stderr, "Error: %s [%d] -- No Data\n",
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
