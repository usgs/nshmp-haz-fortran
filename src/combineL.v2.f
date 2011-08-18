c----program combineL.v2.f
c Harmsen mods of Frankel combineL.f. QA limited checking, feb 21 2007.
c
c--- log -log interpolation
c--- combines eastern and western rate-of-exceedance curves onto 1 grid
c--- finds ground motions for specified annual probability of exceedance
c-- also can output hazard curves for given lat, lon
c try this to compile on Solaris:   f95 combineL.v2.f iosubs.o -o combineL.v2 -e
c
	parameter (nxw=1210,nyw=1530,nxe=1530,nye=1530)
      integer readn
      logical asci,empty_hd
      real prob(14000000),out(4000000),xlev(20)
      type header_type
        character*128  name(6)
        real  period
        integer  nlev
        real  xlev(20)
        real  extra(10)
      end type header_type
      type(header_type) :: head, headr
      character*36 name,nameout,namout
      character namepar*30,name2*30
      dimension pout(nxw,nyw,20),ceus(nxe,nye,20)
      dimension ceus2(nxw,nyw,20),wus(nxw,nyw,20)
c
      write(6,909) "Enter name of CEUS prob. file: "
      read(5,900) name
 900  format(a)
	print 900, name
      call openr(name)
      ndata= 308
      call gethead(head,ndata,readn)
c------ get map parameters
c         open(unit=1,file=head%name(1),status='old')
         headr%name(1)=name
c         read(1,*)igrid	!new
c         read(1,*) ymin,ymax,dy
c get grid loc from header (Harmsen style)
         xmin=head%extra(2)
         xmax=head%extra(3)
         dx=head%extra(4)
         ymin=head%extra(5)
         ymax=head%extra(6)
         dy=head%extra(7)
         empty_hd=xmin.eq.0.
        write(6,*) ymin,ymax,xmin,xmax,dy,dx
        ymine=ymin
        ymaxe=ymax
        if(empty_hd)stop'pgm expecting grid loc data in header. Sorry'
c        close(1)
        nx= nint((xmax-xmin)/dx) +1
        ny= nint((ymax-ymin)/dy) +1
        write(6,*) nx,ny
        nx2= 2*nx
        ny2= 2*ny
        print *,' required CEUS X and y dims ',nx2,ny2
        if(nx2.gt.nxe)stop' increase param nxe '
        if(ny2.gt.nye)stop' increase param nye'
        dx2= dx/2.
        dy2= dy/2.
        nrec= nx*ny
cccc
	score=100.	!start out with a high score. 
c Mark down for minor major incompatibilities between east and west source files.
      nlev= head%nlev
      write(6,*) "nlev=",nlev
      ndata=nlev*nrec
      write(6,*) ymin
      ymin0=ymin
      do 300 i=1,nlev
300   xlev(i)= head%xlev(i)
      call getbuf2(prob,nlev*nrec,readn)
      write(6,*) ndata,readn
      do 10 iy=1,ny
      do 10 ix=1,nx
      j0= ix-1+ (iy-1)*nx
      j0= j0*nlev
      do 11 j=1,nlev
      ceus(ix,iy,j)=prob(j+j0)
  11  continue
  10  continue
      ymin=ymin0
      write(6,*) ymin
c---- interpolate CEUS matrix to .1 degrees
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
      ioffset= (xmin+125)/dx2
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
      write(6,*) "enter name of WUS prob. file"
      read(5,900) name2
	print 900,name2
      call openr(name2)
      ndata= 308
      call gethead(head,ndata,readn)
      headr%period=head%period
c------ get map parameters
c         open(unit=1,file=head%name(1),status='old')
         headr%name(2)=name2
         xmin=head%extra(2)
         xmax=head%extra(3)
         dx=head%extra(4)
         ymin=head%extra(5)
         ymax=head%extra(6)
         dy=head%extra(7)
         empty_hd=xmin.eq.0.
c         read(1,*)igrid	!new
c         read(1,*) ymin,ymax,dy
c        read(1,*) xmin,xmax,dx
        write(6,*) ymin,ymax,dy,xmin,xmax,dx
        if(ymin.ne.ymine)then
        score=score-20.	!remedial action is called for 
        print *,' min lats for east and west do not match '
        endif
        if (ymax.ne.ymaxe)then
        score=score-20.
        print *,'max lats for east and west do not match'
        endif
        close(1)
        nx= nint((xmax-xmin)/dx) +1
        ny= nint((ymax-ymin)/dy) +1
        nx2= 2*nx
        ny2= 2*ny
        nrec= nx*ny
        print *,' required WUS X and y dims ',nx2,ny2
        if(nx2.gt.nxw)stop' increase param nxw '
        if(ny2.gt.nyw)stop' increase param nyw'
cccc
      if(nlev.ne. head%nlev)then
      score=score-30.	!An F for this run would be too generous
      print *,'nlev for east and west not same ',nlev,head%nlev
      endif
      ndata= nlev
      write(6,*) "nlev=",nlev
      do  i=1,nlev
      if(xlev(i).ne. head%xlev(i))then
  	print *,'CEUS and WUS xlev differ ',xlev(i),head%xlev(i)
  	print *,'Are these files compatible?'
  	score = score - 100.*abs(xlev(i)-head%xlev(i))
  	endif
	enddo
  	print *,'Your score so far is ',score,' 100 is an A'
      call getbuf2(prob,nlev*nrec,readn)
      do 30 iy=1,ny
      do 30 ix=1,nx
      j0= ix-1+(iy-1)*nx
      j0 = j0*nlev
      do 31 j=1,nlev
  31  wus(ix,iy,j)=prob(j+j0)
  30  continue
      do 40 j=1,nlev
      do 40 ix=1,nxbig
      do 40 iy=1,nybig
      wus(ix,iy,j)= wus(ix,iy,j)+ceus2(ix,iy,j)
 40   continue
c----
      writE(6,909) "enter 1 to just output hazard curves for sites: "
909 	format(a,$)
      read(5,*) ihaz
	print *,ihaz
      if(ihaz.eq.1) then
 51     write(6,*) "enter lat, lon of site (zero for end)"
        read(5,*) xlat,xlon
        if(xlat.eq.0.) go to 99
        write(6,909) "enter name of output file for haz. curve: "
        read(5,900) nameout
	print *,trim(nameout)
        open(unit=8,file=nameout,status='new')
        ix= (xlon-xmin)/dx2 +1.5
        iy= (ymax-xlat)/dy2 +1.5
        write(8,*) xlat,xlon,iy,ix
        write(8,900) headr%name(1)
        write(8,900) headr%name(2)
        do 50 j=1,nlev
  50    write(8,*) xlev(j),wus(ix,iy,j)
        close(8)
        go to 51
        endif
c----
      write(6,909) "Enter name of output file: "
      read(5,900) namout
	print *,namout
      write(6,699)'Do you want this to be ascii 1 = yes: '
699	format(a,$)
      read *,ii
	print *,ii
      asci=ii.eq.1
        if(asci)then
        open(unit=8,file=namout,status='unknown')
        else
      call openwx(ifp,namout)
      endif
      write(6,909) "Enter annual probability of exceedance: "
      read(5,*) pex
	print *,pex
      pex= alog(pex)
c      write(6,*) "enter rec number"
c      read(5,*) irec
c      write(6,*) prob(irec)
c      do 1 i=1,ndata
c 1    prob(i)= 1.-exp(-p(i)*t)
c
c
      write(6,909) "Enter 0 to output prob gnd motion; 1 for Poisson prob: "
      read(5,*) ip
	print *,ip
      if(ip.eq.1) then
         write(6,*) "enter gnd motion level index,t"
         read(5,*) ilev,t
         if(ilev.lt.1)then
         ilev=1
         print *,'reset ilev to 1'
         score=score-10.
         elseif(ilev.gt.nlev)then
         ilev=nlev
         print *,'reset ilev to nlev'
         score=score-10.
         endif
   	print *,'Your current score is ',score,' 100 is an A'
         ymax= 0.
         do 100 iy=1,nybig
         y = ymax-(iy-1)*dy2
         do 100 ix=1,nxbig
         x= xmin +(ix-1)*dx2
         i= ix+(iy-1)*nxbig
         out(i)= 1.- exp(-wus(ix,iy,ilev)*t)
         out(i)= 100.*out(i)
c             out(i)= wus(ix,iy,ilev)
	if(asci)write(8,200)x,y,out(i)
200	format(f8.2,1x,f7.2,1x,e11.5)
        if(out(i).gt.samax) then
        samax= out(i)
	xloc=x
	yloc=y
	endif
100	continue
  	if(asci)then
  	close(8)
  	print *,'Location of SA or PGAmax is ',xloc,yloc,' max= ',samax
  	stop
  	endif
         call openwx(ifp,namout)
         ndata=308
         headr%extra(9)= ilev
         headr%extra(10)= t
         call puthead(ifp,headr,ndata,ndata)
         call putbufx(ifp,out,nrecbig,readn)
         write(6,*) nrec,readn
         write(6,*) ymax
         go to 99
         endif
      nlev2= nlev-1
c--- interpolate to find ground motions with specified pex
      write(6,699) "Enter scale factor: "
      read(5,*) scale
	print *,scale
      do 60 iy=1,nybig
         y = ymax-(iy-1)*dy2
      do 60 ix=1,nxbig
         x= xmin +(ix-1)*dx2
      i= ix+(iy-1)*nxbig
c      write(6,*) i, out(i-1)
      do 61 j=1,nlev2
      if(wus(ix,iy,j).eq.0.) then
          out(i)= 0.
          go to 61
          endif
      y1= alog(wus(ix,iy,j))
      y2= alog(wus(ix,iy,j+1))
      x1= alog(xlev(j))
      x2= alog(xlev(j+1))
      if((j.eq.1).and.(pex.gt.y1)) then
        out(i)= x1
        go to 59
      elseif((j.eq.nlev2).and.(pex.lt.y2)) then
        out(i)= x2
        go to 59
       elseif((y1.ge.pex).and.(y2.le.pex)) then
        out(i)= x1+ (x2-x1)*(pex-y1)/(y2-y1)
        go to 59
        endif
 61   continue
59	if(out(i).ne.0.)then
	 out(i)= scale*exp(out(i))
	 else
	 out(i)= scale*xlev(1)
	 endif
      if(out(i).gt.samax) then
      samax= out(i)
	xloc=x
	yloc=y
	endif
	if(asci)write(8,200)x,y,out(i)
 60   continue
      ymax= 0.
  	print 339,'Location of SA or PGAmax is ',xloc,yloc,' max= ',samax, namout
339	format(a,f7.2,1x,f6.2,a,f7.4,1x,a20)
 800  format(' ymax=',f12.1)
 	if(asci)then
 	close(8)
 	stop
 	endif
      headr%extra(8)= pex
      headr%extra(9)= scale
      headr%extra(10)= samax
c---- output probabilistic ground motions
      ndata=308
      call puthead(ifp,headr,ndata,readn)
      call putbufx(ifp,out,nrecbig,readn)
      write(6,*) nrecbig,readn
 99   continue
      end
