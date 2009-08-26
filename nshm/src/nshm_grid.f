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
         INTEGER :: grid_id
         CHARACTER, POINTER :: grid_name
         REAL*8, POINTER, DIMENSION(:) :: grid_values
      END TYPE NSHM_Grid

      CONTAINS

C ---------------------------------------------------------------------
C --- FORTRAN subroutines go here.
C ---------------------------------------------------------------------

      SUBROUTINE nshm_get_agrid(grid, gname)
         TYPE(NSHM_Grid), INTENT(INOUT) :: grid
         CHARACTER(LEN=*), INTENT(IN) :: gname

         INTEGER*4 :: num_points

         CALL nshm_get_agrid_meta(grid, gname)
         CALL nshm_num_rows(grid, num_points)
         ALLOCATE(grid%grid_values(num_points))
         CALL nshm_get_agrid_data(grid, gname)

      END SUBROUTINE nshm_get_agrid

      SUBROUTINE nshm_get_bgrid(grid, gname)
         TYPE(NSHM_Grid), INTENT(INOUT) :: grid
         CHARACTER*50, INTENT(IN) :: gname

         INTEGER*4 :: num_points

         CALL nshm_get_bgrid_meta(grid, gname)
         CALL nshm_num_rows(grid, num_points)
         ALLOCATE(grid%grid_values(num_points))
         CALL nshm_get_bgrid_data(grid, gname)

      END SUBROUTINE nshm_get_bgrid


      SUBROUTINE nshm_get_mmax(grid, gname)
         TYPE(NSHM_Grid), INTENT(INOUT) :: grid
         CHARACTER*50, INTENT(IN) :: gname

         INTEGER*4 :: num_points

         CALL nshm_get_mmax_meta(grid, gname)
         CALL nshm_num_rows(grid, num_points)
         ALLOCATE(grid%grid_values(num_points))
         CALL nshm_get_mmax_data(grid, gname)

      END SUBROUTINE nshm_get_mmax

      END MODULE GRID
