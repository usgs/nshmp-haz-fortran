#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "adhoc_utils.h"
#include "adhoc_api.h"
#include "adhoc_connection.h"
#include "adhoc_stmt.h"

#define QUERY "INSERT INTO TEST_INSERT (TEST_VARCHAR, TEST_INT, TEST_NUMBER) VALUES  ( :str , :int , :flt )"

void do_single_insert();
void _do_single_insert(char * username, char * password, char * tnsalias);
void __do_single_insert(char * username, char * password, char * tnsalias);

int main(int argc, char ** argv) {
	
	if (argc < 4) {
		fprintf(stderr, "Usage: %s <username> <password> <tnsalias>\n",
			argv[0]
		);
		return EXIT_FAILURE;
	}

	// Initialize the connection
	AdHoc_AuthInfo  info = {argv[3], argv[1], argv[2], argv[3]};
	adhoc_init(info);

	// Do a sample insert
	//do_single_insert();
	__do_single_insert(argv[1], argv[2], argv[3]);

	return EXIT_SUCCESS;
}

void __do_single_insert(char * username, char * password, char * tnsalias) {
	char * str_val = "Hello World";
	int    int_val = 5;
	float  flt_val = 10.5;

	AdHoc_Connection conn;
	AdHoc_Stmt stmt;

	OCIBind *str_val_hp;
	OCIBind *int_val_hp;
	OCIBind *flt_val_hp;

	AdHoc_Field fields[3] = {
		{
			(char *) str_val, (int *) -1, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_STR, ":str",
			(void *) 0
		},
		{
			(int *) &int_val, (int *) -1, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_INT, ":int",
			(void *) 0
		},
		{
			(float *) &flt_val, (int *) -1, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_FLT, ":flt",
			(void *) 0
		}
	};

	adhoc_connect(username, password, tnsalias, &conn);
	stmt = adhoc_stmt_create(QUERY, NULL, 0, &conn);

	adhoc_check_error(conn.errhp, OCIBindByName(
		stmt.stmthp, (OCIBind **) &str_val_hp, conn.errhp,
		(const OraText *) fields[0].name, -1, (void *) fields[0].value,
		(ub4) strlen(fields[0].value)+1, fields[0].type, (void *) 0, (ub2 *) NULL,
		(ub2 *) NULL, (ub4) 0, (ub4 *) 0, OCI_DEFAULT)
	);
	// This inserts the row but the test_varchar value is incorrect (junk)
	//adhoc_stmt_bind(&stmt, &fields[0], &conn);

	adhoc_check_error(conn.errhp, OCIBindByPos(
		stmt.stmthp, &int_val_hp, conn.errhp, 2, (void *) &int_val,
		(sword) sizeof(int_val), SQLT_INT, (void *) NULL, 
		(ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL,
		OCI_DEFAULT)
	);
	// This inserts the row but the test_int value is incorrect (junk)
	//adhoc_stmt_bind(&stmt, &fields[1], &conn);

	adhoc_check_error(conn.errhp, OCIBindByPos(
		stmt.stmthp, &flt_val_hp, conn.errhp, 3, (void *) &flt_val,
		(sword) sizeof(flt_val), SQLT_FLT, (void *) NULL,
		(ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL,
		OCI_DEFAULT)
	);
	// This inserts the row but the test_number value is incorrect (junk)
	//adhoc_stmt_bind(&stmt, &fields[2], &conn);

	adhoc_check_error(conn.errhp, OCIStmtExecute(
		conn.svchp, stmt.stmthp, conn.errhp, 1, 0, (OCISnapshot *) NULL,
		(OCISnapshot *) NULL, OCI_COMMIT_ON_SUCCESS)
	);

}

void _do_single_insert(char * username, char * password, char * tnsalias) {

	char * str_val = "Hello World";
	int    int_val = 5;
	float  flt_val = 10.5;

	AdHoc_Connection conn;
	AdHoc_Stmt stmt;
	AdHoc_Field fields[3] = {
		{
			(char *) str_val, (int *) -1, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_STR, ":str",
			(void *) 0
		},
		{
			(int *) &int_val, (int *) -1, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_INT, ":int",
			(void *) 0
		},
		{
			(float *) &flt_val, (int *) -1, (short) NULL,
			(unsigned short *) NULL, (unsigned short *) NULL, SQLT_FLT, ":flt",
			(void *) 0
		}
	};


	adhoc_connect(username, password, tnsalias, &conn);
	stmt = adhoc_stmt_create((char*) QUERY, fields, 3, &conn);
	adhoc_stmt_execute(&stmt, &conn);
	adhoc_stmt_free(&stmt, &conn);
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
