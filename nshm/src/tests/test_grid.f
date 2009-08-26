      PROGRAM MAIN
      USE Grid
      IMPLICIT NONE

C --------------------------------------------------------------------
C --- Define some variables to be used throughout.                 ---
C --------------------------------------------------------------------

      INTEGER :: num_args
      CHARACTER*25 :: username, password, tnsalias
      CHARACTER*50 :: agrid_name, bgrid_name, mmax_name

      INTEGER :: num_rows, i

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

C --- Set agrid, bgrid, and mmax grid names to fetch
      agrid_name = "adapt_cn"
      bgrid_name = "BVAL_.1"
      mmax_name  = "gm_j_6p6_7p1"

C --------------------------------------------------------------------
C --- Begin main body of program.                                  ---
C --------------------------------------------------------------------


C --- You must always initialize before using the API
      CALL nshm_initialize(username, password, tnsalias)

      WRITE(*,101)
      WRITE(*,102)"Fetching an Agrid"
      WRITE(*,101)
C --- Try to get an agrid
      CALL nshm_get_agrid(grd, agrid_name)
      CALL nshm_print_grid(grd)

C --- Try to write a copy of it back like it was a new grid
C      CALL nshm_put_agrid(grd)

      WRITE(*,101)
      WRITE(*,102)"Fetching a Bgrid"
      WRITE(*,101)
C --- Try to get a bgrid
      CALL nshm_get_bgrid(grd, bgrid_name)
      CALL nshm_print_grid(grd)

C --- Try to write a copy of it back like it was a new grid
C      CALL nshm_put_bgrid(grd)

      WRITE(*,101)
      WRITE(*,102)"Fetching a mmax"
      WRITE(*,101)
C --- Try to get a mmax
      CALL nshm_get_mmax(grd, mmax_name)
      CALL nshm_print_grid(grd)
 
C --- Try to write a copy of it back like it was a new grid
C      CALL nshm_put_mmax(grd)

C --- Check the number of rows returned from the last grid.
      CALL nshm_num_rows(grd, num_rows)
      WRITE(*,103)num_rows

C --- Natively in fortran echo the first 20 rows of the last grid
      WRITE(*,104)
      DO i=1, 20
         WRITE(*,105)grd%grid_values(i)
      END DO

C --- You should always cleanup when you are done with the API
      CALL nshm_cleanup()

C --------------------------------------------------------------------
C --- End main body of program. Below are format statements only.  ---
C --------------------------------------------------------------------

 100  FORMAT("Usage: test_grid_f <user> <pass> <tns>")
 101  FORMAT("######################################################")
 102  FORMAT("## ",A)
 103  FORMAT("The previous grid fetch returned ",I6," rows of data.")
 104  FORMAT("They are...")
 105  FORMAT(F10.7)

      END PROGRAM MAIN
