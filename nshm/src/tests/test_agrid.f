      PROGRAM main
      USE Agrid
      IMPLICIT NONE



C --  Agrid meta parameter values

      CHARACTER(MAX_DESC_LEN) :: desc
      REAL :: values(MAX_GRID_ROWS)
      REAL :: min_lat, max_lat, inc_lat
      REAL :: min_lon, max_lon, inc_lon
      INTEGER :: num_rows

C --  This is just an iterator for a do-loop later

      INTEGER :: i


 100  FORMAT('Name: ',A50,/,'Rows: ',I6,/,'Latitude [min, max, inc]: [',
     + F5.2,', ',F5.2,', ',F5.2,']',/,'Longitude [min, max, inc]: [',
     + F7.2,', ',F7.2,', ',F5.2,']')
 101  FORMAT('     ',F7.5)


      CALL get_agrid( desc, values, min_lat, max_lat,
     +      inc_lat, min_lon, max_lon, inc_lon, num_rows )

      WRITE(*, 100) desc, num_rows, min_lat, max_lat, inc_lat,
     +  min_lon, max_lon, inc_lon

      DO i = 1, 10 !num_rows
          WRITE(*,101) values(i)
      END DO

      STOP
      END PROGRAM main
