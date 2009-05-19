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
        integer desc_len
        character, pointer, dimension (:) :: description
C        character*50 description ! Don't use or depend on this!!!
      END TYPE metainfo
C
      REAL values(1048576)
      INTEGER i
      TYPE(metainfo) :: meta
C
C
      allocate(meta%description(0))
      CALL fetchagrid(values, meta)
C
      write(*, 100, advance="no")meta%id
      do i = 1, meta%desc_len
         write(*,103,advance="no")meta%description(i)
      end do
      write(*,102)meta%num_rows
 100  format('ID: ',I2,' Name: ')
 101  format(F6.5)
 102  format(' Rows: ',I6)
 103  format(A1)
C
      do 20 i = 1, 10 !meta%num_rows
          write(*, 101) values(i)
 20   continue
C
      END
