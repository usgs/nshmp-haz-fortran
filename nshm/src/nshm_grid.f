      MODULE GRID
      IMPLICIT NONE

C --- 07/08/09 -- EMM:
C --- FORTRAN does not play nice with C when dealing with dynamically
C --- allocated values. For this reason we choose to statically allocate
C --- the size of 
      INTEGER, PARAMETER :: MAX_GRID_ROWS = 524288

C --- 07/08/09 -- EMM: 
C --- FORTRAN does not play nicely with dynamically allocated strings when
C --- attempting to pass them between C and FORTRAN. For this reason we
C --- instead use a static allocation of this length. This value is based
C --- on the maximum field lenght for the grid_name column in the database.
C --- If the database changes its column width, so too should we change this
C --- value as well.
      INTEGER, PARAMETER :: MAX_GRID_NAME_LEN  = 50


      TYPE NSHM_Grid
         REAL*8 :: lat_min
         REAL*8 :: lat_max
         REAL*8 :: lat_inc
         REAL*8 :: lng_min
         REAL*8 :: lng_max
         REAL*8 :: lng_inc
         REAL*8 :: grid_values(MAX_GRID_ROWS)
         CHARACTER(LEN=MAX_GRID_NAME_LEN) :: grid_name
         INTEGER :: grid_id
      END TYPE NSHM_Grid

      CONTAINS

C ---------------------------------------------------------------------
C --- FORTRAN subroutines go here.
C ---------------------------------------------------------------------

      END MODULE GRID
