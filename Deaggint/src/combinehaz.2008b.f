c----program combinehaz.2008b.f;  from  haz2asc.f (Frankel)
c Compile with f95: 
c f95 combinehaz.2008b.f -o ../bin/combinehaz.2008b iosubs.o
c--- 
c---- revised from combinehaz2007.f 7/23/08 for binary write
c---Read header for geographic bounds. If header is empty, use
c   hardwired lat, lon, dlat, dlon
c-- outputs hazard curves, prob. gnd motion  for given lat, lon
c this program has been tested on SOlaris machines. Not tested on Windows.
c Does limited checking for compatibility of input files.
c Steve Harmsen USGS april 2008
      type header_type
        character*30  name(6)
        real  period
        integer  nlev
        real  xlev(20)
        real  extra(10)
      end type header_type
      type(header_type) :: head
      integer*4 readn
      dimension out(14000000),xlev(20)
      character*80 name,nameout,namout,namehaz
      character namepar*80,name2*80,line*78
      dimension pout(1210,530,20),ceus(530,270,20)
      dimension ceus2(1210,530,20),wus(1210,530,20)
c
      write(6,5) "enter name of CEUS prob. file: "
5	format(a,$)
      read(5,900) name
 900  format(a)
      call openr(name)
      ndata= 308
      call gethead(head,ndata,readn)
c------ get map parameters
c         open(unit=1,file=head%name(1),status='old')
c         read(1,*) idum
c         read(1,*) ymin,ymax,dy
c        read(1,*) xmin,xmax,dx
c information should be in the header records.

	if (head%extra(5).eq.0.0.and.head%extra(6).eq.0.)then
c standard ceus grid location.
        ymin= 24.6 
        ymax= 50.0
        dy =0.1
        xmin= -115.0
        xmax= -65.0
        dx= 0.1
	   print *,'Default CEUS grid is defined'
	else
      ymin=head%extra(5)
      print *,'CEUS grid location is read from header record'
      ymax=head%extra(6)
      xmin=head%extra(2)
      xmax=head%extra(3)
        dx=head%extra(4)
        dy=head%extra(7)
      endif
        write(6,*) ymin,ymax,xmin,xmax
        ymaxe=ymax; ymine=ymin
      print *,'These coordinates are for Eastern or CEUS.'
        close(1)
        nx= nint((xmax-xmin)/dx) + 1
        ny= nint((ymax-ymin)/dy) + 1
        write(6,*) nx,ny
        nx2= 2*nx
        ny2= 2*ny
        dx2= dx/2.
        dy2= dy/2.
        nrec= nx*ny
cccc
      nlev= head%nlev
      write(6,*) "nlev ceus=",nlev
      xmaxc=xmax
      ndata=nlev*nrec
      write(6,*) ymin
      ymin0=ymin
      do 300 i=1,nlev
300   xlev(i)= head%xlev(i)
      call getbuf2(out,nlev*nrec,readn)
      write(6,*) ndata,readn
      do 10 iy=1,ny
      do 10 ix=1,nx
      j0= ix-1+ (iy-1)*nx
      j0= j0*nlev
      do 11 j=1,nlev
      ceus(ix,iy,j)=out(j+j0)
  11  continue
  10  continue
      ymin=ymin0
      write(6,*) ymin
c---- interpolate CEUS matrix to .05 degrees
      do 105 ilev=1,nlev
c---- place existing values into new matrix
      do 1 iy=1,ny
      iy2= (2*iy)-1
      do 1 ix=1,nx
      ix2= (2*ix)-1
   1  pout(ix2,iy2,ilev)= ceus(ix,iy,ilev)
c-- fill in even x-values for odd rows
      do 2 iy=1,ny2,2
      do 2 ix=2,nx2-1,2
   2  pout(ix,iy,ilev)= 0.5*(pout(ix-1,iy,ilev)+pout(ix+1,iy,ilev))
c-- fill in even rows
      do 3 iy=2,ny2-1,2
      do 3 ix=1,nx2
  3   pout(ix,iy,ilev)= 0.5*(pout(ix,iy-1,ilev)+pout(ix,iy+1,ilev))
      ioffset= (xmin+125.0)/dx2
      do 20 ix=1,nx2
      do 20 iy=1,ny2
  20  ceus2(ix+ioffset,iy,ilev)= pout(ix,iy,ilev)
 105  continue
      call close(name)
      write(6,*) ymin
      nxbig= (xmax+125.)/dx2  +1
      nybig= (50.-ymin)/dy2  +1
      nrecbig= nxbig*nybig
      write(6,*) nxbig,nybig,nrecbig
cc------read WUS probs
      write(6,5) "Enter name of WUS prob. file: "
      read(5,900) name2
      call openr(name2)
      ndata= 308
      call gethead(head,ndata,readn)
c------ get map parameters
	if (head%extra(5).eq.0.0.and.head%extra(6).eq.0.)then
c standard WUS grid location.
       ymin= 24.6
       ymax= 50.0
       xmin= -125.0
       xmax= -100.
       dy= 0.05
       dx= 0.05
	   print *,'Default WUS grid is defined'
	else
C TRY TO USE INFORMATION STORED IN HEADER RECORD
      xmin=head%extra(2)
      xmax=head%extra(3)
        dx=head%extra(4)
        dy=head%extra(7)
      ymin=head%extra(5)
      print *,'WUS grid location is read from header record'
      ymax=head%extra(6)
      nlev= head%nlev
      write(6,*) "nlev wus=",nlev
        if(dy.ne.dy2.or.dx.ne.dx2)then
        write(6,*)'Spatial sampling not compatible between E and W'
        stop'sampling W has to be double that of sampling E'
        endif
      endif
        write(6,*) ymin,ymax,xmin,xmax
        if(ymin.ne.ymine.or.ymax.ne.ymaxe)then
	write(6,*)'Eastern lats ',ymine,ymaxe,' western lats ',ymin,ymax
        stop' region incompatibility'
        endif
        nx= nint((xmax-xmin)/dx) + 1
        ny= nint((ymax-ymin)/dy) + 1
        nx2= 2*nx
        ny2= 2*ny
        nrec= nx*ny
cccc
      ndata= nlev
      do 43 i=1,nlev
  43  if(xlev(i).ne. head%xlev(i))
     + print *,'warning SA samples differ E&W',xlev(i),head%xlev(i)
      call getbuf2(out,nlev*nrec,readn)
      do 30 iy=1,ny
      do 30 ix=1,nx
      j0= ix-1+(iy-1)*nx
      j0 = j0*nlev
      do 31 j=1,nlev
  31  wus(ix,iy,j)=out(j+j0)
  30  continue
  	i=0
c Send strong to outer, weak to inner orbits.
      do 40 iy=1,nybig
      do 40 ix=1,nxbig
      do 40 j=1,nlev
      i=i+1
      out(i)= wus(ix,iy,j)+ceus2(ix,iy,j)
 40   continue
c----
	head%name(1)='combinehaz.2008b.f USGS'
	head%name(2)='PSHA data for NSHMP 2008'
	head%name(3)='Data to be used on WEB Int. Deagg'
	head%name(4)='Ground motion units: g or gravity'
	head%extra(2)=xmin
	head%extra(3)=xmaxc
	 head%extra(8)=float(nrecbig)
      write(6,*) "enter name of binary output file for haz. curves"
      read(5,900) nameout
      write(6,*) nrecbig
      nhd=308
      ndata=nrecbig*nlev
      call openwx (ip,nameout)
       call puthead(ip,head,nhd,nput)
       call putbufx(ip,out,ndata,nput)
       print *,'binary file written ',nameout
 99   continue
      end
