#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "adhoc_api.h"
#include "adhoc_utils.h"

const char * QUERY = "INSERT INTO TEST_INSERT (TEST_VARCHAR, TEST_INT, \
                      TEST_NUMBER) VALUES  ( :str, :int, :flt )";

void do_single_insert();

int main(int argc, char ** argv) {

	// Initialize the connection
	AdHoc_AuthInfo  info = {"teamt", "haz_owner", "ham23*ret", "teamt"};
	adhoc_init(info);

	// Do a sample insert
	do_single_insert();

	return EXIT_SUCCESS;
}

void do_single_insert() {
	char * str_val = "Hello World";
	int    int_val = 5;
	float  flt_val = 10.5;

	OCIStmt *stmt_hp;
	OCIBind *str_val_hp;
	OCIBind *int_val_hp;
	OCIBind *flt_val_hp;

	adhoc_check_error(adhoc_err_hp, OCIStmtPrepare2(
		adhoc_svc_hp, &stmt_hp, adhoc_err_hp, (const OraText *) QUERY,
		strlen((char *) QUERY), (OraText *) NULL, 0, OCI_NTV_SYNTAX,
		OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIBindByPos(
		stmt_hp, &str_val_hp, adhoc_err_hp, 1, (void *) str_val,
		(sb4) strlen(str_val)+1, SQLT_STR, (void *) 0, 
		(ub2 *) 0, (ub2 *) 0, (ub4) 0, (ub4 *) 0,
		OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIBindByPos(
		stmt_hp, &int_val_hp, adhoc_err_hp, 2, (void *) &int_val,
		(sword) sizeof(int_val), SQLT_INT, (void *) NULL, 
		(ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL,
		OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIBindByPos(
		stmt_hp, &flt_val_hp, adhoc_err_hp, 3, (void *) &flt_val,
		(sword) sizeof(flt_val), SQLT_FLT, (void *) NULL,
		(ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL,
		OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIStmtExecute(
		adhoc_svc_hp, stmt_hp, adhoc_err_hp, 1, 0, (OCISnapshot *) NULL,
		(OCISnapshot *) NULL, OCI_COMMIT_ON_SUCCESS)
	);

}
