#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "adhoc_api.h"
#include "adhoc_check_error.h"

#define RESULT_SET_SIZE 10

void query_agrid_name(int _agrid_id);
void query_agrid_names();

int main(int argc, char ** argv) {
	
	AdHoc_AuthInfo int_dev;

	int_dev.username  = "haz_owner";
	int_dev.password  = "ham23*ret";
	int_dev.tns_alias = "teamt";

	/* Initialize connection to database */
	adhoc_init(int_dev);

	/* Do stuff with the connection */
	printf("Checking two rows individually (we know there are only two)...\n");
	query_agrid_name(1);
	query_agrid_name(2);
	printf("\n");

	printf("Checking all rows with one call.\n");
	query_agrid_names();
	printf("\n");

	/* Close connection to database */
	adhoc_close();

	printf("Done!\n");
	return 0;
}

void query_agrid_names() {
	char * query = "SELECT AGRID_NAME FROM AGRID_META ORDER BY AGRID_ID";
	OCIDefine *defn_hp;
	OCIStmt   *stmt_hp;
	char      name[RESULT_SET_SIZE][50];
	sb2       name_idx[RESULT_SET_SIZE];
	ub2       name_len[RESULT_SET_SIZE];
	boolean done = FALSE;
	ub4 rows = 0;
	ub4 i = 0;
	sb4 status;

	/* Create the statement */
	adhoc_check_error(adhoc_err_hp, OCIStmtPrepare2(
		adhoc_svc_hp, &stmt_hp, adhoc_err_hp, (const OraText *) query, 
		strlen((char *) query), (OraText *) NULL, 0, OCI_NTV_SYNTAX, 
		OCI_DEFAULT)
	);

	/* Execute the statement */
	adhoc_check_error(adhoc_err_hp, OCIStmtExecute(
		adhoc_svc_hp, stmt_hp, adhoc_err_hp, 0, 0, (OCISnapshot *) NULL,
		(OCISnapshot *) NULL, OCI_DEFAULT)
	);

	/* Define the output buffer */
	adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
		stmt_hp, &defn_hp, adhoc_err_hp, 1, (void *) name[0], 
		(sb4) sizeof(name[0]), SQLT_STR, (void *) name_idx, (ub2 *) name_len,
		(ub2 *) NULL, OCI_DEFAULT)
	);

	/* Fetch results */
	while (!done) {
		status = OCIStmtFetch(stmt_hp, adhoc_err_hp, RESULT_SET_SIZE,
			OCI_FETCH_NEXT, OCI_DEFAULT);

		if (status == OCI_SUCCESS || status == OCI_NO_DATA) {
			if (status == OCI_SUCCESS) {
				rows = RESULT_SET_SIZE;
			} else if (status == OCI_NO_DATA) {
				adhoc_check_error(adhoc_err_hp, OCIAttrGet(
					stmt_hp, OCI_HTYPE_STMT, &rows, (ub4*) NULL, 
					OCI_ATTR_ROWS_FETCHED, adhoc_err_hp)
				);
				done = TRUE;
			}
		} else {
			adhoc_check_error(adhoc_err_hp, status);
			done = TRUE;
		}

		for ( i=0; i<rows; ++i) {
			printf("%s\n", name[i]);
		}
	} 

	adhoc_check_error(adhoc_err_hp, OCIStmtRelease(
		stmt_hp, adhoc_err_hp, (OraText *) NULL, 0, OCI_DEFAULT)
	);
	
}

void query_agrid_name(int _agrid_id) {
	char * query = "SELECT AGRID_NAME FROM AGRID_META WHERE AGRID_ID = :id";
	OCIBind   *bind_hp;
	OCIDefine *defn_hp;
	OCIStmt   *stmt_hp;
	char      name[50];
	ub4       num_rows = 1;

	adhoc_check_error(adhoc_err_hp, OCIStmtPrepare2(
		adhoc_svc_hp, &stmt_hp, adhoc_err_hp, (const OraText *) query, 
		strlen((char *) query), (OraText *) NULL, 0, OCI_NTV_SYNTAX, 
		OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIBindByName(
		stmt_hp, &bind_hp, adhoc_err_hp, (text *) ":id", -1, 
		(void *) &_agrid_id, sizeof(_agrid_id), SQLT_INT, (void *) NULL, 
		(ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL, OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
		stmt_hp, &defn_hp, adhoc_err_hp, 1, (void *) &name, (sb4) sizeof(name), 
		SQLT_STR, (void *) NULL, (ub2 *) NULL, (ub2 *) NULL, OCI_DEFAULT)
	);

	adhoc_check_error(adhoc_err_hp, OCIStmtExecute(
		adhoc_svc_hp, stmt_hp, adhoc_err_hp, num_rows, 0, (OCISnapshot *) NULL,
		(OCISnapshot *) NULL, OCI_DEFAULT)
	);

	printf("%s\n", name);

	adhoc_check_error(adhoc_err_hp, OCIStmtRelease(
		stmt_hp, adhoc_err_hp, (OraText *) NULL, 0, OCI_DEFAULT)
	);
}
