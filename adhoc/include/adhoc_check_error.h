#ifndef ADHOC_CHECK_ERROR_H
#define ADHOC_CHECK_ERROR_H

/**
 * Checks the return status and logs the errors. Returns an AdHoc return status
 * rather than the OCI status. One should defer calling this method if the OCI
 * return status is required but should always call this method following and
 * OCI call in order to keep valid logging.
 *
 * void * _handle_pointer (IN) -- The handle to be inspected in the case of an
 * error. In generall this is an OCIError handle but in the case of a call to
 * OCIEnvCreate (or similar methods) this may be an OCIEnv handle.
 *
 * ub4 _handle_type (IN) -- The OCI handle type identifier. See OCI
 * documentation for details.
 *
 * sword _status (IN) -- The OCI return status from an OCI method call.
 */
int adhoc_oci_check(void *_handle_pointer, ub4 _handle_type, sword _status);

// Macro wrapping our error checker for the OCIError handle case
#define adhoc_check_error(_error_hp, _status) adhoc_oci_check((_error_hp), OCI_HTYPE_ERROR, (_status))

// Macro wrapping our error checker for the OCIEnv handle case
#define adhoc_check_env(_env_hp, _status) adhoc_oci_check((_env_hp), OCI_HTYPE_ENV, (_status))

#endif
