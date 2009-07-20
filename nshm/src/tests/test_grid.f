      PROGRAM MAIN
      USE Grid
      IMPLICIT NONE

C --------------------------------------------------------------------
C --- Define some variables to be used throughout.                 ---
C --------------------------------------------------------------------

      INTEGER :: num_args
      CHARACTER*25 :: username, password, tnsalias

      REAL*8 :: lat_min, lat_max, lat_inc
      REAL*8 :: lng_min, lng_max, lng_inc
      REAL*8 :: grid_values(MAX_GRID_ROWS)
      CHARACTER(MAX_GRID_NAME_LEN) :: grid_name
      INTEGER :: grid_id, num_rows

      TYPE(NSHM_Grid) :: small_grd
      TYPE(NSHM_Grid) :: grd
      num_args = IARGC()

C --- Check for proper usage
      IF (num_args.LT.3) THEN
         WRITE(*,100)
         STOP
      ENDIF

C --- Read in command line arguments
      CALL GETARG(1, username)
      username = TRIM(username)

      CALL GETARG(2, password)
      password = TRIM(password)

      CALL GETARG(3, tnsalias)
      tnsalias = TRIM(tnsalias)
        
C --------------------------------------------------------------------
C --- Begin main body of program.                                  ---
C --------------------------------------------------------------------


C --- You must always initialize before using the API
      CALL nshm_initialize(username, password, tnsalias)

      WRITE(*,101)
      WRITE(*,102)"Fetching an Agrid"
      WRITE(*,101)
C --- Try to get an agrid
      CALL nshm_get_agrid(grd)
      CALL nshm_print_grid(grd)

C --- Try to write a copy of it back like it was a new grid
      CALL nshm_put_agrid(grd)

      WRITE(*,101)
      WRITE(*,102)"Fetching a Bgrid"
      WRITE(*,101)
C --- Try to get a bgrid
      CALL nshm_get_bgrid(grd)
      CALL nshm_print_grid(grd)

C --- Try to write a copy of it back like it was a new grid
      CALL nshm_put_bgrid(grd)

      WRITE(*,101)
      WRITE(*,102)"Fetching a mmax"
      WRITE(*,101)
C --- Try to get a mmax
      CALL nshm_get_mmax(grd)
      CALL nshm_print_grid(grd)
 
C --- Try to write a copy of it back like it was a new grid
      CALL nshm_put_mmax(grd)

C --- Check the number of rows returned from the last grid.
      CALL nshm_num_rows(grd, num_rows)
      WRITE(*,103)num_rows

C --- You should always cleanup when you are done with the API
      CALL nshm_cleanup()

C --------------------------------------------------------------------
C --- End main body of program. Below are format statements only.  ---
C --------------------------------------------------------------------

 100  FORMAT("Usage: test_grid_f <user> <pass> <tns>")
 101  FORMAT("######################################################")
 102  FORMAT("## ",A)
 103  FORMAT("The previous grid fetch returned ",I6," rows of data.")

      END PROGRAM MAIN
