#include <stdio.h> /* Might not need this for production but good for dev. */
#include <stdlib.h>
#include <string.h>

#include "oci.h"
#include "adhoc_api.h"
#include "adhoc_check_error.h"


OCIEnv    *adhoc_env_hp;
OCIError  *adhoc_err_hp;
OCISvcCtx *adhoc_svc_hp;

int adhoc_init(AdHoc_AuthInfo _auth_info) {
	int status; 
	/* Create the environment handle */
	status = adhoc_check_env(adhoc_env_hp, OCIEnvCreate(
		&adhoc_env_hp, OCI_DEFAULT, 
		NULL, NULL, NULL, NULL,
		(size_t) 0, (void **) NULL)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }

	/* Create the error handle */
	status = adhoc_check_error(adhoc_err_hp, OCIHandleAlloc(
		adhoc_env_hp, (void **) &adhoc_err_hp, OCI_HTYPE_ERROR,
		(size_t) 0, (void **) NULL)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }

	/* Create the service context */
	status = adhoc_check_error(adhoc_svc_hp, OCIHandleAlloc(
		adhoc_env_hp, (void **) &adhoc_svc_hp, OCI_HTYPE_SVCCTX, 
		(size_t) 0, (void **) NULL)
	);
	if (status >= ADHOC_RETURN_WARNING) { return status; }

	/* Connect to server and create user session */
	status = adhoc_check_error(adhoc_err_hp, OCILogon2(
		adhoc_env_hp, adhoc_err_hp, &adhoc_svc_hp,
		(text *) _auth_info.username, (ub4) strlen(_auth_info.username),
		(text *) _auth_info.password, (ub4) strlen(_auth_info.password),
		(text *) _auth_info.tns_alias, (ub4) strlen(_auth_info.tns_alias),
		OCI_DEFAULT)
	);
	if (status > ADHOC_RETURN_WARNING) { return status; }
	/**
	 * There are actually other ways we could try to connect but for now if this
	 * method fails we will just assume failure.
	 */
	
	return ADHOC_RETURN_SUCCESS;
}

int adhoc_query(AdHoc_Statement * _stmt) {
	return ADHOC_RETURN_WARNING;
}


int adhoc_close() {
	int status;
	status = adhoc_check_error(adhoc_err_hp, OCILogoff(
		adhoc_svc_hp, adhoc_err_hp)
	);
	if (status > ADHOC_RETURN_WARNING) { return status; }

	return ADHOC_RETURN_SUCCESS;
}
