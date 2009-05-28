      PROGRAM main
      USE Agrid
      IMPLICIT NONE
C
C
C
C --  Agrid meta parameter values
C
      CHARACTER(MAX_DESC_LEN) :: desc
      REAL :: values(MAX_GRID_ROWS)
      REAL :: min_lat, max_lat, inc_lat
      REAL :: min_lon, max_lon, inc_lon
      INTEGER :: num_rows
C
C --  This is just an iterator for a do-loop later
C
      INTEGER :: i
c      TYPE(NSHM_AgridMeta) :: meta
C
C
 100  FORMAT('Name: ',A50,/,'Rows: ',I6,/,'Latitude [min, max, inc]: [',
     + F5.2,', ',F5.2,', ',F5.2,']',/,'Longitude [min, max, inc]: [',
     + F7.2,', ',F7.2,', ',F5.2,']')
 101  FORMAT('     ',F7.5)
C
C
      CALL get_agrid( desc, values, min_lat, max_lat,
     +      inc_lat, min_lon, max_lon, inc_lon, num_rows )
C
c      meta%description = ""
c      CALL fetchagrid(values, meta)

c      WRITE(*, 100) TRIM(meta%desc), meta%num_rows
      WRITE(*, 100) desc, num_rows, min_lat, max_lat, inc_lat,
     +  min_lon, max_lon, inc_lon
C
      DO i = 1, 10 !num_rows
          WRITE(*,101) values(i)
      END DO
C
      STOP
      END PROGRAM main
