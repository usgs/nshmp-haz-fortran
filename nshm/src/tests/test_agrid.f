      PROGRAM EXAMPLE
C
C
      TYPE metainfo
        integer id
        real min_lat
        real max_lat
        real lat_inc
        real min_lng
        real max_lng
        real lng_inc
        integer num_rows
        character*50 description ! Don't use or depend on this!!!
      END TYPE metainfo
C
      REAL values(1048576)
      INTEGER i
      TYPE(metainfo) :: meta
C
C
      meta%description = "Do not trust strings between Fortran and C"
      CALL fetchagrid(values, meta)
C
      write(*, 100) meta%id, meta%description, meta%num_rows
 100  format('ID: ',I2,' Name: ',A50,' Rows: ',I6)
C
      do 20 i = 1, meta%num_rows
          write(*, 101) values(i)
 20   continue
 101  format(F6.5)
C
      END
