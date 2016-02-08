c--- bin2xyzX: converts from C-binary file to ascii file 
c Usage bin2xyzX.v2 input-binary-file  output.file
c v2: Input file may have a header record befor grid info.
c   SHarmsen jan 12 2007
c The gridded mmax, agrid, and bgrid however do not have header records.
c
c SUN compile: f95 bin2xyzX.v2.f iosubs.o -o bin2xyz -e
c PC Compile: gfortran -static -O -ffixed-line-length-none -ffpe-trap= -o bin2xyzX.v2.exe
c      bin2xyzX.v2.f iosubs.o
c
	parameter (NBIG=3200000)
      type  header 
        character*30 name(6)
        real*4 period
        integer*4 nlev
        real*4 xlev(20)
        real*4 extra(10)
      end type header
      type (header) :: headr
      integer*4 readn
      logical ishere
	character*1 ch1
      dimension z(NBIG)
      character name*72,name2*72,name1*72
      if(iargc().lt.2)stop 'Usage bin2xyzX.v2 infile(binary) outfile (asc)'
      call getarg(1,name1)
      call getarg(2,name2)
      inquire(file=name2,exist=ishere)
      if(ishere)then 
      print *,'File ',name2(1:30),' exists'
	print 50,'Overwrite y/n? '
50	format(a,$)
5	format(a)
	read 5,ch1
	if(ch1.eq.'n')then
      print 50,'Please enter another name : '
      read *,name2
	endif
      endif
      open(2,file=name2,status='unknown')
      write(2,25)name1
25	format('#bin2xyzX.v2 converts ',a24,' to ascii')
C       ymin=24.6; ymax=50.0    !latitude range WUS. In 2007 ymax=51.0
	dy =0.1
	print 50,'Enter ymin and ymax in degrees: '
	read *, ymin,ymax
c       extended grid has xmin = -126.5 2007 PSHA
	print 50,'Enter xmin and xmax in degrees: '
	read *,xmin,xmax
	print 50,'Enter dx and dy (degrees) [0.1]: '
	read *,dx,dy
      nx= nint((xmax-xmin)/dx ) +1
      ny= nint((ymax-ymin)/dy) +1
      nrec= nx*ny
      if(nrec.gt.NBIG)then
      print *,'nrec, NBIG = ',nrec, NBIG
      STOP' Please  increase parameter NBIG and retry'
      endif
      print 60,'Below "std" means standard or 308-byte header.'
	write(6,50)'Does file have a header record: 0=no; 1=yes,std; 2=yes,nonstd: '
60	format(a)
	read *,i
      call openr(name1)
      if(i.eq.1)then
      ndata= 308
      elseif(i.eq.2)then
      print 50,'Enter the size (number of bytes ) of the header record: '
      read *,ndata
      else
      i=0
      endif
      if(i.ge.1)call gethead(headr,ndata,readn)
      call getbuf2(z,nrec,readn)
      do 10 i=1,nrec
      iy= (i-1)/nx 
      ix= i- 1- iy*nx
      xlat= ymax- iy*dy
      xlon= xmin+ ix*dx
      write(2,20) xlon,xlat,z(i)
20	format 	(f9.2,1x,f8.2,1x,e11.5)
 10   continue
 	print 88,'Output file is ',name2
88	format(a,a30)
      end
