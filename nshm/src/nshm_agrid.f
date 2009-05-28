      MODULE AGRID
      IMPLICIT NONE
C
C
      INTEGER, PARAMETER :: MAX_GRID_ROWS = 524288
      INTEGER, PARAMETER :: MAX_DESC_LEN = 50
C
C
      TYPE NSHM_AgridMeta
        INTEGER :: id
        REAL :: min_lat
        REAL :: max_lat
        REAL :: inc_lat
        REAL :: min_lon
        REAL :: max_lon
        REAL :: inc_lon
        INTEGER :: num_rows
        CHARACTER(LEN=MAX_DESC_LEN) :: description
      END TYPE NSHM_AgridMeta
C
C
      CONTAINS
C
C
C
    
      SUBROUTINE get_agrid ( desc, values, min_lat, max_lat,
     +      inc_lat, min_lon, max_lon, inc_lon, num_rows )
         CHARACTER*50, INTENT(OUT) :: desc
         REAL, INTENT(OUT) :: values(MAX_GRID_ROWS)
         REAL, INTENT(OUT) :: min_lat, max_lat, inc_lat
         REAL, INTENT(OUT) :: min_lon, max_lon, inc_lon
         INTEGER, INTENT(OUT) :: num_rows
C
C
         INTEGER :: i
         TYPE(NSHM_AgridMeta) :: meta
C
C --  Set the description to a null character for the benefit of C
c         meta%description = ""
         CALL fetchagrid(values, meta)

         desc = TRIM(meta%description)
         min_lat = meta%min_lat
         max_lat = meta%max_lat
         inc_lat = meta%inc_lat
         min_lon = meta%min_lon
         max_lon = meta%max_lon
         inc_lon = meta%inc_lon
         num_rows = meta%num_rows
C
C
C
      END SUBROUTINE get_agrid
C
C
      END MODULE AGRID
