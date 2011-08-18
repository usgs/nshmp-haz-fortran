c--- gethead: reads header of PSHA hazcurve file
c-- to use: gethead.nga filename or names
c -- to compile: f95 gethead.nga.f -o gethead.nga iosubs.o
      type header 
        character*128 region,filin(5)
        real period
        integer nlev
        real xlev(20)
	real xtra(10)
      end type header
      type(header):: headr
      integer readn,ifi,i,nfi
      real v30o
      character name*88
      nfi=iargc()
      if(nfi.eq.0)stop'Usage gethead.nga binary.files'
  	do ifi=1,nfi
      call getarg(ifi,name)
      call openr(name)
      write(6,909)ifi,name
 909	format('File ',i2,': ',a)     
c 30 character name
c      ndata= 308
c 128 character name
      ndata= 896
      call gethead(headr,ndata,readn)
 900  format(a)
	write(6,900)headr%region
      do 1 i=1,5
  1   write(6,900) headr%filin(i)
      write(6,*) headr%period, ' spectral period(s)'
	write(6,*)' Ground motion levels (g):'
	write(6,26)(headr%xlev(j),j=1,headr%nlev)
26	format(6f10.6)
	print *,'Location of grid and other tidbits'
	write(6,*)(headr%xtra(j),j=2,7)
	write(6,*)int(headr%xtra(8)),' number of receivers in the grid'
	write(6,*)headr%xtra(9),' is default Vs30 in m/s'
	write(6,*)headr%xtra(10),
     + ' depth at which Vs30 is 2500 m/s (km)'
	if(ifi.gt.1 .and. v30o .ne. headr%xtra(9))
     + write(6,900) 'warning: vs30 different'
	v30o=headr%xtra(9)
	write(6,900)
	enddo
      end
