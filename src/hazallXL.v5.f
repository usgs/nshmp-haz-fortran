c----program hazallXL.v5.f Last rev Dec 19 2012. Stephen Harmsen
c Use either the 308- or the 128-char file name (header rec) in this version,
c Here, the length of the header record is specified on input.
c Solaris:
c f95 hazallXL.v5.f iosubs_128.o -e -o hazallXL.v5
c    -e to extend line length 
c Linix cluster, gldlabm:
c	ifort hazallXL.v5.f iosubs_128.o -O -132 -o ../Bin/hazallXL.v5
c What about -static flag? or -i-static (Intel lib routines)
c
c based on hazallXL.f (frankel). THe big gain in v5 is increased flexibility of
c input files. No longer need to fill up large arrays of zero or 10**-21 curves. Only
c compute for sites within Rmax km of the sources = faults. Previously, all input grids
c had to be defined the same. There were a few code versions (EW) that allowed moderate flexibility.
c --- The "malloc" calls may fail with gfortran compilers. Arrays could be declared
c      with sufficient length to bypass this problem. -fcray-pointer flag should make
c	the gfortran compile a big success.
c This code works when inputting multiple files, all w/different  sampling boxes.
c Can have different NS range as well as EW range, with total, partial or no overlap.
c However, the sampling intervals dx and dy must be the same for all I/O files
c
c usage hazallXL.v5 input.file 
c Note: no longer assumes terminal input, no < input file unlike the earliest version
c Note 2: has been modified for possible ascii output SH.
c  Get data grid info from headers of the input files.
c---- Hazard-curve to point estimate at spec. PE: use log - log interpolation
c--- sums annual prob. of exceed.
c--- from output of hazgridX and hazFX
c--- 
c--- Finds ground motions for specified annual rate of exceedance
c --- Or, finds probability of exceeding an input ground motion level (g) for time t 
c---- Or, can output file with hazard curves
c--- to run: hazallXL.v5  control.script 
c --- control script describes a grid and names several input files, all within/on bdry of grid
c --- File order does not matter. Could be a northern and southern or any reasonably nearby
c --- grids. Continent-scale ok but dont try to grid up the world with this pgm. 
c You will be disappointed if you ask for too much memory.
c--- this version gets most required control info from header records.
c--- input files have hazard curves for each site concatenated
c--- So does the binary output file. 
c --- v5 : first two lines of control file describe the grid of sites to output.
c --- Both I/O assume a rectangular grid  of receivers
c --- The data in the header records are checked to
c determine if different input files in the control script are compatible.
c If not compatible, the reason often is that dx and dy are different.
c If they are 2x coarser you can use hazinterp.nga.f to make dx and dy 1/2 their
c previous values, gaining the compatibility that is needed here.
c Different longitude sampling is  acceptable in this program as long as dx is same
c and overlap occurs or would occur if extended in one direction or another (E or W).
c Arrays are dynamicallly allocated. If you demand too much from your computer, it
c should give you a painful kick and stop running.
c Steve Harmsen, Sept 2009. Oct 2010. 
c Minor rev Mar 25 2008:
c if two supposedly same SA or PGA levels differ by more than 1e-6 in a relative sense
c print a message but otherwise don't. If they differ by more than 0.0001 in an absolute
c sense, the code will shut down. this is same as previous versions.
c Mar 25: Also, identify PGD which has a period index of -2
c May 2009: Add probability of exceeding a specified uniform gm level. Previous
c versions only output a given gm index. This version interpolates (ln-ln).
      type header 
        character*128 name(6)
        real period
        integer nlev
        real xlev(20)
      real extra(10)
      end type header
      type(header):: headr,head
c below, the traditional short header. You can mix 'em up in hazallXL.v5.
      type sheader 
        character*30 name(6)
        real period
        integer nlev
        real xlev(20)
      real extra(10)
      end type sheader
      type(sheader):: shead
      integer readn,nshort/0/,nbad/0/,nmal/0/,nvschang/0/
      integer n2read,nread,i0,nlev,nlevm
      logical ok,new_model,vsmatch/.true./,first,ascii,lastread,short
c      logical, dimension(nfmax) :: eof
c vsmatch is a logical variable that will be false if one or more input files
c have discrepant vs30 compared to the first input file. Vs30 is in new headers.
      real dx,dy,scale,gymax,gymin,gxmax,gxmin,gdx,gdy,tiny,tinytim
c      parameter (nfmax=60,tinytim=1.e-20,tiny=1.e-12)   !max number of input files
c      parameter (nfmax=60,tinytim=1.e-20,tiny=1.e-6)   !max number of input files
      pointer (pp1,prob1),(pp2,prob)
      real, dimension(*):: prob, prob1
      dimension xlev(20),ylev(20)
c
      parameter (nfmax=60,tinytim=1.e-20,tiny=1.e-6)   !max number of input files
      character*80 name(nfmax),namout,nom_de_plume*30      
      integer, dimension(nfmax) :: ifp,nx,ny,recstart,ndata,nrec
      real, dimension (nfmax):: wt
c File names can be long, but over 30 chars no good for header storage.
c
      if(iargc().lt.1)stop 'Enter input file of instructions on command 
     * line, no <'
      call getarg(1,name(1))
      infi=2
      pex=1.	!must initialize >0
c      eof=.false.
      lastread=.false.
      open(infi,file=name(1),status='old',err=1066)
	recstart=1	!starting record
	kstart=1	!starting block of reads
      nom_de_plume=name(1)
       print 900,"This program is hazallXL.v5 (Harmsen, Dec 2012)"
       print 900,'Different W&E long and N&S lat ranges are OK. All 
     * data must be inside master grid'
       print 5,"Enter output latmin,latmax,dlat : "
       read (infi,*)gymin,gymax,gdy
       print *,gymin,gymax,gdy
       print 5,"Enter output longmin,longmax,dlong: "
       read (infi,*)gxmin,gxmax,gdx
       print *, gxmin,gxmax,gdx
       nxout=nint((gxmax-gxmin)/gdx)+1
       nyout=nint((gymax-gymin)/gdy)+1
       nrecout=nxout*nyout
      write(6,5) "enter number of input files: "      
      read(infi,*) nfile
      print *,nfile
      if(nfile.gt.nfmax)stop'increase param nfmax'
      do 90 ifil=1,nfile
3      write(6,*)
      write(6,5) "Enter name of hazard-curve file:  "
5      format(a,$)
50      format(a,$)
      read(infi,900) name(ifil)
      write(6,900)name(ifil)
 900  format(a)
      inquire(file=name(ifil),exist=ok)
      if(ok)then
c      call openr(name(ifil))
      call openr1(name(ifil))
      else
      print *,name(ifil),' not found. Fix your script.'
      stop 'hazallXL.v5 did not finish'
      endif
       write(6,5) "Enter header size(308 or 896) and weight for these 
     * hazard curves: "
      read(infi,*) nhd,wt(ifil)
      print *,nhd,wt(ifil)	!echo input
      if(nhd.eq.308)then
c      print *,'short head read?'
      call getshead(shead,nhd,readn)
c      print *,'yes ',readn
      short=.true.
      nameln=30
       nlev= shead%nlev
      xmin=shead%extra(2) 
      ymin=shead%extra(5) 
      xmax=shead%extra(3)
      ymax=shead%extra(6)
        dxh=shead%extra(4)
        dyh=shead%extra(7)
	vs30o=shead%extra(9)
c	print *,xmin,ymin,xmax,ymax,dxh,dyh,vs30o
	elseif (nhd .eq. 896)then
	short=.false.
	nameln=96
c      print *,'Long head read?'
      call gethead(head,nhd,readn)
c      print *,'yes ',readn
       nlev= head%nlev
      xmin=head%extra(2) 
      ymin=head%extra(5) 
      xmax=head%extra(3)
      ymax=head%extra(6)
        dxh=head%extra(4)
        dyh=head%extra(7)
	vs30o=head%extra(9)
      else
      stop 'Please confine your header sizes to 308 and 896'
      endif
c------ get or check map parameters
      if(ifil.lt.7)headr%name(ifil)= name(ifil)(1:nameln)
        write(*,*) 'ymin = ', ymin
        if(ymax.gt.gymax+tiny) stop 'ymax gt gymax'
        if(ymin.lt.gymin-tiny) stop 'ymin lt gymin'
        if(xmin.lt.gxmin-tiny) stop 'xmin lt gxmin'
        if(xmax.gt.gxmax+tiny) stop 'xmax gt gxmax'
      dx=dxh; dy=dyh
        nx(ifil)= nint((xmax-xmin)/dx) +1
        ny(ifil)= nint((ymax-ymin)/dy) +1
        ndata(ifil)=nx(ifil)*ny(ifil) *nlev    !each file may have a diff. length, ndata
        nrec(ifil)=nx(ifil)*ny(ifil)
        if(dxh.ne.gdx)stop 'dx of file does not match dx of output. 
     * Unacceptable'
        if(dyh.ne.gdy)stop 'dy of file does not match dy of output. 
     * Unacceptable'
       if(ifil.eq.1)then
c some checks will be made as new files are read in to determine data compatibility.
c the first file can extend further south than others. This is not necessarily an incompatibility.
c the output will go as far south as the first file. (based on nrecin)
      nlevm=nlev-1
c set xlev based on 1st input file.
	if(nhd.eq.896)then
      do  i=1,nlev
      xlev(i)= head%xlev(i)
      headr%xlev(i)=xlev(i)
      ylev(i)=alog(xlev(i))
      enddo
        do i=1,10
         headr%extra(i)=head%extra(i)
         enddo      !extra info is not extra.
      else
      head%nlev=shead%nlev
      do  i=1,nlev
      xlev(i)= shead%xlev(i)
      head%xlev(i)= shead%xlev(i)
      headr%xlev(i)=xlev(i)
      ylev(i)=alog(xlev(i))
      enddo
        do i=1,10
         headr%extra(i)=shead%extra(i)
         head%extra(i)=shead%extra(i)
         enddo      !extra info is not extra.
	endif      
	ndataout=nxout*nyout*nlev
c allocate space for output
		pp2=malloc(ndataout*4)
		prob(1:ndataout) = 1.e-21	!initialize. Note explicit mention of all elements
	nrow_out=nxout*nlev
	write(12,*)ndataout,nrow_out,' ndataout nrow_out'
c initialize a few output variables
	endif	!ifil = 1
	nin= ndata(ifil)
	pp1=malloc(nin*4)
	call getbuf2(prob1,nin,nrd)
	recstart(ifil)=1+ (nxout*nint((gymax-ymax)/dyh)+nint((xmin-gxmin)/dxh))*nlev
	print *,' file ndata ',ndata(ifil),' recstart ',recstart(ifil)
	if(ifil.gt.1)then
c compare with first file
      do i=1,nlev
      if(short)head%xlev(i) = shead%xlev(i)
      if(abs(xlev(i)-head%xlev(i)).gt.0.0001)then
      write(6,*) head%xlev(i),xlev(i),' not the same for sample ',i
      if(i.lt.19)then
      nbad=nbad+1
      else
      nmal=nmal+1
      endif
      elseif(abs(1.-head%xlev(i)/xlev(i)).ge.1.e-6)then
c very often the representation of floating point numbers will differ on diff machines
c by about 1e-7 but this is irrelevant for hazard estimation. Assume differences of 1e-6 g
c are worthy of note but not different enough to shut down the process.
      write(6,*)head%xlev(i),xlev(i),' same number?'
      print *,'They are close enough to continue'
      endif
      enddo
        endif	!ifil > 1
c      write(6,*)wt(ifil)
	weight=wt(ifil)
	prob1(1:nin) = prob1(1:nin)*weight	!note explicit mention of all elements.
c the compiler thinks this array has dim(1) even though dynamic allocation is much bigger.
c You have to be explicit about which elements get multiplied. I learned the hard way.
	if(short)then
	 period=shead%period
	 do i=1,10
	 head%extra(i)=shead%extra(i)
	 enddo
	 else
         period=head%period
         endif
         if(period.gt.0.)then
         print *,'Spectral period is ',period,' s'
         elseif (period.eq.0.)then
         print *,"PGA is indicated by zero period"
         elseif(period.gt.-2.)then
         print *,'PGV (units cm/s) is indicated by -1 period'
         elseif(period.gt.-3.)then
         print *,'PGD (units cm) is indicated by -2 period (CB only?)'
         else
c add Cumulative Absolute Velocity label, CAV 10/20/2010 SH.
         print *,'CAV (units g-s) is indicated by -3 period (CB only?)'
         endif
      write(6,*) (head%extra(i),i=1,10)
      if (head%extra(9).ne.0.)vs30n=head%extra(9)
      if(ifil.gt.1.and.vs30o.ne.vs30n)then
      nvschang=nvschang+1
      print *,'*** Above Vs30 is different from Vs30 of previous file ***'
      print *,'********************  Data compatible? *****************'
      vsmatch=.false.
      endif
      vs30o=vs30n
66      format('Input file period mismatch. Now ',f4.2,' was ',f4.2,' s.',/,
     + '*** This is by most standards a serious error. ***')
c      write(6,*) head%nlev
      if(nbad.gt.0)stop 'mismatch of sampled pga or sa.'
      if(nmal.gt.0)print *,'Continuing, should be ok for gm<',
     + head%xlev(19)
	readn=0
	m=1
	nrow_in=nx(ifil)*nlev
	n=nrow_in
	k=recstart(ifil)
	l=k+nrow_out
	write(12,*)'#File number/weight: ',ifil,weight,' start point ',k
	do j=1,ny(ifil)
c	k=k+nrow_out
	do ii=1,n
c	write(12,*)k,l,m,n
c	prob(k:l)= prob(k:l)+ prob1(m:n)	!this was ok but was replaced with explicit loop
	prob(k)=prob(k) + prob1(m)
	k=k+1
	m=m+1
	enddo	!ii loop
c	m=m+nrow_in
c	n=n+nrow_in
c left hand side index must move to same start posn in next row down
	k=l
	l=l+nrow_out
	enddo
10	continue
	call free (pp1)
90	continue
       headr%nlev= head%nlev
      headr%extra(2)=gxmin
      headr%extra(5)=gymin
      headr%extra(3)=gxmax
      headr%extra(6)=gymax
	headr%extra(8)=float(nrecout)
      headr%period=head%period
      write(6,*)
      write(6,50) "Enter 1 to just output summed prob file"
      read(infi,*) iout
      print *,iout
      if(iout.eq.1) then
        write(6,5) "enter name of output file: "
        read(infi,900) namout
        write(6,900)namout
      print 50,'Enter a 1 for ascii output; 0 for binary: '
      read (infi,*)ii
      write(6,*)ii
      if(ii.eq.0)then
       ascii=.false.
        call openwx(ifw,namout,len)
        nhd=896
	write(6,*)'Output binary file has long header structure 896 bytes'
        call puthead(ifw,headr,nhd,readn)
        else      !ascii output
        ascii=.true.
        open(3,file=namout,status='unknown')
        write(3,32)nfile,nom_de_plume,head%period,nlev
32      format('#Pgm hazallXL.v5.f (harmsen) sums ',i2,' hazard curves from ',
     + a,/,'#Lat Long   Rex for spectral period ',f4.2,/,'#GM set(g)',i3)      
        write(3,39)(headr%xlev(i),i=1,nlev)
        endif	!ascii or bin
        elseif(iout.eq.0)then
         write(6,5) "enter name of output file: "
        read(infi,900) namout
        write(6,900)namout
      print 50,'Enter a 1 for ascii output; 0 for binary: '
      read (infi,*)ii
      write(6,*)ii
             if(ii.eq.0)then
      call openwx(ifp,namout)
      ascii=.false.
      else
      open(3,file=namout,status='unknown')
33      format('#Pgm hazallXL.v5.f (harmsen) sums ',i2,' hazard curves. Site Vs30',f6.1,
     + /,'#Long    Lat  SA(g) for spectral period ',f4.2,' Annual Rex ',1pe11.5) 
34      format('#Pgm hazallXL.v5.f (harmsen) sums ',i2,' hazard curves. site Vs30',f6.1,
     + /,'#Long    Lat  PGV(cm/s) or PGD for index ',f4.1,' Rex ',1pe11.5) 
35      format('#Pgm hazallXL.v5.f (harmsen) sums ',i2,' hazard curves. Site Vs30',f6.1,
     + /,'#Long        Lat    CAV(g-s). Annual Rex ',1pe11.5) 
      ascii=.true.     
      endif
      samax= 0.
      nread=0
      i0=1
      write(6,50) "Enter 0 to output prob gnd motion; 1 for Poisson prob: "
      read(infi,*) ip
      print *,ip
	if(ip.lt.1)then
      write(6,50) "Enter mean ann. rate of exceedance: "
      read(infi,*) pex
      print *,pex
      if(ascii.and.head%period.ge.0.)then
      write(3,33)nfile,vs30o, head%period,pex
      elseif(ascii.and.head%period.ge.-2.0)then
      write(3,34)nfile,vs30o, head%period,pex
	else
      write(3,35)nfile,vs30o,pex
      endif	!different header labels
      pex= alog(pex)
      write(6,50) "Enter scale factor: "
      read(infi,*) scale
      print *,scale
      endif	!ip < 1
        endif	!iout=0 or 1
39      format(e11.5)
	if(iout.eq.1.and.ii.eq.0)then
        call putbufx(ifw,prob,ndataout,readn)
        write(6,*) 'writing binary ',ndata(3),readn
        elseif(iout.eq.1)then
 	       k1=1
 	       k2=nlev
	      do i=1,nrecout
	      xi= gxmin + mod(i-1,nxout)*dx
	      yi= gymax - (i-1)/nxout*dy
	      write(3,49)yi,xi,(prob(j), j=k1,k2)
	      k1=k1+nlev
	      k2=k2+nlev
	      enddo
	      elseif(iout.eq.0.and.ii.eq.0)then
	      stop'only ascii output currently permitted with this option'
	 elseif(iout.eq.0.and.ip.eq.0)then
	 do 85 i=1,nrecout			!	***
	      xi=gxmin + mod(i-1,nxout)*dx
	      yi=gymax - (i-1)/nxout*dy
      do 20 j=1,nlevm
      ind= (i-1)*nlev +j
      if(prob(ind).le.tinytim) then
        out= ylev(1)
        go to 80
        endif
      y1= alog(prob(ind))
      y2= alog(prob(ind+1))
      x1= ylev(j)
      x2= ylev(j+1)
      if((j.eq.1).and.(pex.gt.y1)) then
        out= x1
        go to 80
        endif
      if((j.eq.nlevm).and.(pex.lt.y2)) then
        out= x2
        go to 80
        endif
      if((y1.ge.pex).and.(y2.le.pex)) then
        tmp= y2-y1
        if(tmp.eq.0.) then
        out=x1
        else
        out= x1+ (x2-x1)*(pex-y1)/tmp
        endif
        go to 80
        endif
 20   continue
 80   continue
 	out=scale*exp(out)
 	write(3,49)xi,yi,out
       if(out.gt.samax)then
        samax= out
        xistar=xi
        yistar=yi
        endif
85 	continue		!the current section of the data ***
	elseif(iout.eq.0.and.ip.eq.1)then
c probability of exceeding a specific ground motion.
         write(6,50) "enter gnd motion in g and t (yrs): "
         read(infi,*) gm,t
         write(3,349)gm,t
349	format('#Long  Lat  PROB. of SA or PGA>',f6.3,' g in ',f5.1,' yrs.')
         write(6,*)gm, t
        gml=alog(gm)
c log-log interpolation tends to be more accurate.
         j=2
         dowhile(gml.gt.ylev(j).and.j.lt.nlev)
         j=j+1
         enddo
         jm=j-1
      x1= ylev(jm)
      x2= ylev(j)
      frac1=(x2-gml)/(x2-x1)
      frac2=1.0-frac1
      prmax=0.0
         do 100 i=1,nrecout
         ind= (i-1)*nlev + jm
         indp=ind+1
c log(ann rate)
         rate=frac1*alog(prob(ind))+frac2*alog(prob(indp))
         rate=exp(rate)
c back to annual rate space, then onwards to probability.
      prx= 1.- exp(-rate*t)
	      xi=gxmin + mod(i-1,nxout)*dx
	      yi=gymax - (i-1)/nxout*dy
      prmax=max(prmax,prx)
      if(prx.eq.prmax)then
      xistar=xi
      yistar=yi
      endif
100      write(3,49)xi,yi,prx
	
      endif	!iout=0 or 1
49      format(f9.3,1x,f9.3,1p,20(1x,e11.5))
 11   continue
 	if(iout.eq.0.and.ip.eq.0)then
 	print *,'Location of max motion ',xistar,yistar,' max ',samax
      close(3)
      stop'Ascii map file written'
      elseif(iout.eq.0.and.ip.gt.0)then
      	print *,'Location of max probability ',xistar,yistar,' max ',prmax
      close(3)
      stop'Ascii prob. map file written'
      elseif(iout.eq.1.and.ii.eq.0)then
      stop'Binary hazard curve file written'
      elseif(iout.eq.1.and.ii.eq.1)then
      close(3)
      stop'Ascii hazard curve file written'
      endif       !iout=0
c how on earth did we arrive here?
       stop 'hazallXL.v5 normal end'
 1066      stop 'input file on command line not found'
      end
