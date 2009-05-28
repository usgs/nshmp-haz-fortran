#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "oci.h"

#include "adhoc_api.h"
#include "adhoc_check_error.h"

#include "nshm_agrid_meta.h"
#include "nshm_agrid.h"

int nshm_get_agrid_meta(int id, NSHM_Agrid * _agrid);
int nshm_get_agrid_data(int id, NSHM_Agrid * _agrid);

/**
 * Grabs a random agrid from the database.
 *
 * 05/11/09 -- EMM: In this initial "proof of concept", we don't randomize but
 *                  rather just choose to return "agrid #1".
 */
int nshm_get_random_agrid(NSHM_Agrid * _agrid) {
    int status   = ADHOC_RETURN_SUCCESS;
    int agrid_id = 1; /* Will randomize later */

    status = nshm_get_agrid_meta(agrid_id, _agrid);
    if ( status > ADHOC_RETURN_WARNING ) {
        // Maybe print/log something?
        return status;
    }

    status = nshm_get_agrid_data(agrid_id, _agrid);
    if ( status > ADHOC_RETURN_WARNING ) {
        // Maybe print/log something?
        return status;
    }

    return ADHOC_RETURN_SUCCESS;
}

int nshm_get_agrid_meta(int _id, NSHM_Agrid * _agrid) {
    
    char * query = "SELECT AGRID_NAME, MIN_LAT, MAX_LAT, MIN_LON, MAX_LON \
    	FROM AGRID_META WHERE AGRID_ID = :id ";

    OCIBind   *bind_hp;
    OCIDefine *defn_hp;
    OCIDefine *defn_min_lat_hp;
    OCIDefine *defn_max_lat_hp;
    OCIDefine *defn_min_lng_hp;
    OCIDefine *defn_max_lng_hp;
    OCIStmt   *stmt_hp;
    char      name[50];
    float     min_lat, max_lat;
    float     min_lng, max_lng;
    ub4       num_rows = 1;
    NSHM_AgridMeta info;

    adhoc_check_error(adhoc_err_hp, OCIStmtPrepare2(
        adhoc_svc_hp, &stmt_hp, adhoc_err_hp, (const OraText *) query,
        strlen((char *) query), (OraText *) NULL, 0, OCI_NTV_SYNTAX,
        OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIBindByName(
        stmt_hp, &bind_hp, adhoc_err_hp, (text *) ":id", -1,
        (void *) &_id, sizeof(_id), SQLT_INT, (void *) NULL,
        (ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL, OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_hp, adhoc_err_hp, 1, (void *) &name, (sb4) sizeof(name),
        SQLT_STR, (void *) NULL, (ub2 *) NULL, (ub2 *) NULL, OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_min_lat_hp, adhoc_err_hp, 2, (void *) &min_lat,
        (sb4) sizeof(min_lat), SQLT_FLT, (void *) NULL, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_max_lat_hp, adhoc_err_hp, 3, (void *) &max_lat,
        (sb4) sizeof(max_lat), SQLT_FLT, (void *) NULL, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );
    
    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_min_lng_hp, adhoc_err_hp, 4, (void *) &min_lng,
        (sb4) sizeof(min_lng), SQLT_FLT, (void *) NULL, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );
    
    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_max_lng_hp, adhoc_err_hp, 5, (void *) &max_lng,
        (sb4) sizeof(max_lng), SQLT_FLT, (void *) NULL, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );
    
    adhoc_check_error(adhoc_err_hp, OCIStmtExecute(
        adhoc_svc_hp, stmt_hp, adhoc_err_hp, num_rows, 0, (OCISnapshot *) NULL,
        (OCISnapshot *) NULL, OCI_DEFAULT)
    );

    /* Set the agrid meta data information */
    info.id = _id;
    info.min_lat = min_lat;
    info.max_lat = max_lat;
    info.min_lng = min_lng;
    info.max_lng = max_lng;
	memset(info.description, 0, NSHM_AGRID_META_DESC_LEN);
    strncpy(info.description, name, strlen(name));

	_agrid->metadata = calloc(1, sizeof(NSHM_AgridMeta));
    memmove(_agrid->metadata, &info, sizeof(NSHM_AgridMeta));

    adhoc_check_error(adhoc_err_hp, OCIStmtRelease(
        stmt_hp, adhoc_err_hp, (OraText *) NULL, 0, OCI_DEFAULT)
    );
    
    return ADHOC_RETURN_SUCCESS;
}

int nshm_get_agrid_data(int _id, NSHM_Agrid * _agrid) {

    char * query = "SELECT LAT, LON, AGRID_VAL FROM AGRID WHERE AGRID_ID = \
		:id ORDER BY LAT DESC, LON ASC";
    OCIBind   *bind_hp;
    OCIDefine *defn_lat_hp;
    OCIDefine *defn_lng_hp;
    OCIDefine *defn_val_hp;
    OCIStmt   *stmt_hp;
    float lats[NSHM_MAX_AGRID_LENGTH];
    float lngs[NSHM_MAX_AGRID_LENGTH];
    float vals[NSHM_MAX_AGRID_LENGTH];
    float lats_ind[NSHM_MAX_AGRID_LENGTH];
    float lngs_ind[NSHM_MAX_AGRID_LENGTH];
    float vals_ind[NSHM_MAX_AGRID_LENGTH];
    boolean done = FALSE;
    ub4 rows = 0;
    sb4 status;

    adhoc_check_error(adhoc_err_hp, OCIStmtPrepare2(
        adhoc_svc_hp, &stmt_hp, adhoc_err_hp, (const OraText *) query,
        strlen((char *) query), (OraText *) NULL, 0, OCI_NTV_SYNTAX,
        OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIBindByName(
        stmt_hp, &bind_hp, adhoc_err_hp, (text *) ":id", -1,
        (void *) &_id, sizeof(_id), SQLT_INT, (void *) NULL,
        (ub2 *) NULL, (ub2 *) NULL, 0, (ub4 *) NULL, OCI_DEFAULT)
    );
    
    adhoc_check_error(adhoc_err_hp, OCIStmtExecute(
        adhoc_svc_hp, stmt_hp, adhoc_err_hp, 0, 0, (OCISnapshot *) NULL,
        (OCISnapshot *) NULL, OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_lat_hp, adhoc_err_hp, 1, (void *) lats,
        (sb4) sizeof(lats[0]), SQLT_FLT, (void *) lats_ind, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );
    
    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_lng_hp, adhoc_err_hp, 2, (void *) lngs,
        (sb4) sizeof(lngs[0]), SQLT_FLT, (void *) lngs_ind, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );

    adhoc_check_error(adhoc_err_hp, OCIDefineByPos(
        stmt_hp, &defn_val_hp, adhoc_err_hp, 3, (void *) vals,
        (sb4) sizeof(vals[0]), SQLT_FLT, (void *) vals_ind, (ub2 *) NULL,
        (ub2 *) NULL, OCI_DEFAULT)
    );

    while (!done) {
        status = OCIStmtFetch(stmt_hp, adhoc_err_hp, NSHM_MAX_AGRID_LENGTH,
            OCI_FETCH_NEXT, OCI_DEFAULT);

        if ( (status == OCI_SUCCESS) || (status == OCI_NO_DATA) ) {
            if (status == OCI_SUCCESS) {
                rows = NSHM_MAX_AGRID_LENGTH;
            } else if (status == OCI_NO_DATA) {
                /* Might have gotten fewer than allowed. Find out how many*/
                adhoc_check_error(adhoc_err_hp, OCIAttrGet(
                    stmt_hp, OCI_HTYPE_STMT, &rows, (ub4 *) NULL,
                    OCI_ATTR_ROWS_FETCHED, adhoc_err_hp)
                );
                done = TRUE;
            }

            /* Allocate memory for these. Never freed. */
            _agrid->latitude  = malloc(sizeof(float) * rows);
            _agrid->longitude = malloc(sizeof(float) * rows);
            _agrid->value     = malloc(sizeof(float) * rows);

            /* Assign the computed results */
            _agrid->metadata->num_rows = rows;
            memmove(_agrid->latitude, lats, sizeof(float) * rows);
            memmove(_agrid->longitude, lngs, sizeof(float) * rows);
            memmove(_agrid->value, vals, sizeof(float) * rows);
        } else {
            status = adhoc_check_error(adhoc_err_hp, status);
            return status;
        }
    } // END: while(!done)

    adhoc_check_error(adhoc_err_hp, OCIStmtRelease(
        stmt_hp, adhoc_err_hp, (OraText *) NULL, 0, OCI_DEFAULT)
    );

    return ADHOC_RETURN_SUCCESS;
}
