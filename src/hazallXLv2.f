c----program hazallXL.v2.f
c f95 hazallXL.v2.f iosubs.o -e -o hazallXL.v2
c based on hazallXL.f (frankel). -e to extend line length
c Should be compatible with gfortran on PC Linux cluster
c usage hazallXL.v2 input.file (note, no longer assume terminal input, no < input file)
c modified for possible ascii output SH.
c doesn't necessarily need an ascii control file. Can get info needed from header
c---- Hazard-curve to point estimate at spec. PE: use log - log interpolation
c--- sums annual prob. of exceed.
c--- from output of hazgridX and hazFX
c--- finds Poisson prob. of exceed. for time t or
c--- finds ground motions for specified annual probability of exceedance
c --- Or, finds probability of exceeding an input ground motion level (g)
c--- also can output file with hazard curves
c--- to run: hazallXL.v2  control.script 
c--- this version gets parameter file from header of 1st input file
c--- input files have hazard curves for each site concatenated
c---parfile= parameter file as used in hazgrid, etc. 
c --- v2 : param file has an extra record at beginning with a 0 for
c --- gridded set of receivers, or a 1 for a list of receivers.
c --- If new version of input file, then the data in the header can be checked to
c determine if different input files in the control script are compatible.
c If not compatible, the reason often is that dx and dy are different.
c If they are 2x coarser you can use hazinterp.nga.f to make dx and dy 1/2 their
c previous values, gaining the compatibility that is needed here.
c If max or min latitude is the only difference program tries to accomodate.
c The first file should have min(ymax) and max(ymin) of all regions sampled.
c If so, program will just toss out the extra data from subsequent files.
c Different longitude sampling is not acceptable at this time.
c Steve Harmsen, May 17 2006. Revised June 30 2006. 
c Minor rev Mar 25 2008:
c if two supposedly same SA or PGA levels differ by more than 1e-6 in a relative sense
c print a message but otherwise don't. If they differ by more than 0.0001 in an absolute
c sense, the code will shut down. this is same as previous versions.
c Mar 25: Also, identify PGD which has a period index of -2
c May 2009: Add probability of exceeding a specified uniform gm level. Previous
c versions only output a given gm index. This version interpolates (ln-ln).
      parameter (nbig=12615000)   !how many data can be read in?
      integer nmin/5120000/      !how many data are actually read in?
      type header 
        character*128 name(6)
        real period
        integer nlev
        real xlev(20)
      real extra(10)
      end type header
      type(header):: headr,head
      integer readn,nshort/0/,nbad/0/,nmal/0/,nvschang/0/,nboob/0/
      logical ok,new_model,vsmatch/.true./,first,ascii
c vsmatch is a logical variable that will be false if one or more input files
c have discrepant vs30 compared to the first input file. Vs30 is in new headers.
c short_s becomes true if one or more files didnot compute to southernmost
c latitude, but came close. This is not too bad because Mexico lats for example
c do not need to be computed for WUS hazard maps      
      real latmax,latmin,dx,dy
      real, dimension(nbig):: prob, prob0
      dimension out(300000),xlev(20)
      character name*80,namout*80,nom_de_plume*30      
c File names can be long, but over 30 chars no good for header storage.
c
      if(iargc().lt.1)stop 'Enter input file of instructions on command line, no <'
      call getarg(1,name)
      infi=2
      open(infi,file=name,status='old',err=1066)

      nom_de_plume=name
 9      prob0= 0.
       print *,"This program is hazallXL.v2 (Harmsen)"
       print 5,"Enter a 1 if new_style input with 0 for grid: "
       read (infi,*)i
       print *,i
       new_model=i.eq.1
      write(6,5) "enter number of input files "      
c nfile could be 1 if you want to compute the ground motion for a given
c mean rate of exceedance using a pre-existing hazard curve.
      read(infi,*) nfile
      print *,nfile
      do 11 ifil=1,nfile
c      write(6,*) ifil
3      write(6,*)
      write(6,5) "Enter name of hazard-curve file:  "
5      format(a,$)
      read(infi,900) name
      write(6,900)name
 900  format(a)
      inquire(file=name,exist=ok)
      if(ok)then
      call openr(name)
      else
      print *,name,' not found. Fix your script.'
      stop 'hazallXL.v2 did not finish'
      endif
c  for 30 characters
c     ndata= 308
c  for 128 characters
      ndata= 896
      call gethead(head,ndata,readn)
c------ get or check map parameters
      if(new_model)then
        nxh= nint((head%extra(3)-head%extra(2))/head%extra(4)) +1
        nyh= nint((head%extra(6)-head%extra(5))/head%extra(7)) +1
        nrec1=nxh*nyh      !nrec1 will be a model size for I/O
        dxh=head%extra(4)
        dyh=head%extra(7)
      xmin=head%extra(2)
      ymin=head%extra(5)
      xmax=head%extra(3)
      ymax=head%extra(6)
      dx=dxh
      dy=dyh
      elseif(ifil.eq.1) then
         open(unit=1,file=head%name(1),status='old')
         headr%name(1)=head%name(1)
         if(new_model)then
         read(1,*)igrid
         if(igrid.ne.0)stop'grid of sites was expected'
         endif
        read(1,*) ymin,ymax,dy
        read(1,*) xmin,xmax,dx
c for WUS, ymin,ymax are typically 24.6 and 50 degrees
c for WUS, xmin,xmax are -125 and -100 degrees. dx=dy=0.1 degrees (1996)
        close(1)
      endif
       if(ifil.eq.1)then
c some checks will be made as new files are read in to determine data compatibility
      nlev= head%nlev
c set xlev based on 1st input file.
      do 1 i=1,nlev
  1   xlev(i)= head%xlev(i)
        nx= nint((xmax-xmin)/dx) +1      !improved grid size calc
        ny= nint((ymax-ymin)/dy) +1
        nrecin= nx*ny
       latmax=ymax
      latmin=ymin
        do i=1,10
         headr%extra(i)=head%extra(i)
         enddo      !extra info is not extra.
         period1=head%period
         endif
         period=head%period
         if(period.gt.0.)then
         print *,'Spectral period is ',period,' s'
         elseif (period.eq.0.)then
         print *,"PGA is indicated by zero period"
         elseif(period.gt.-2.)then
         print *,'PGV (units cm/s) is indicated by -1 period'
         else
         print *,'PGD (units cm) is indicated by -2 period (CB only?)'
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
      nrec=nrecin
cccc check for compatibility. 2005 and earlier files might not have some useful header info.
c 2007&later files do if written from a run of the 2007 -2008 code updates.
      if(new_model.and.
     + (nxh.ne.nx.or.nyh.ne.ny.or.dxh.ne.dx.or.dyh.ne.dy))then
      print *,'Header says geographic region not same or not sampled same as 1st'
      print *,'LatitudeN ',latmax,head%extra(6)
      print *,'LatitudeS ',ymin,head%extra(5)
      print *,'LongitudeW ',xmin,head%extra(2)
      print *,'LongitudeE ',xmax,head%extra(3)
      if(head%extra(6).gt.latmax)then
      print *,'Current file has receiver data to north of first file.'
      nyomit = nint((head%extra(6)-latmax)/dyh)
      nomit= nyomit*nxh*nlev
      call getbuf2(prob,nomit,nread)
      print *,'Number of northern rows latitude omitted ',nyomit
      nyh=nyh-nyomit
      endif
      print *,dxh,dx,' dx headr vs dx ascii'
      print *,dyh,dy,' dy headr vs dy ascii'
      print *,nxh,nx,' nx computed from current header versus nx from first'
      print *,nyh,ny,' ny computed from current header versus ny from first'
      if(nxh.eq.nx.and.dyh.eq.dy.and.latmax.eq.head%extra(6).and.nyh.ge.ny)then
      print *,'Current file may go further south. Extra data are omitted'
      nrec=nrec1
      elseif(nxh.eq.nx.and.dyh.eq.dy.and.latmin.gt.head%extra(5)-3.)then
c within 3 degrees: for USA scale, might be good enough
           print *,'Current file cuts off before minimum latitude of others'
           print *,'We will proceed but you should be aware of possible loss'
           print *,'Current ',head%extra(5),' general min ',latmin
           nrec=nx*nint((latmax-head%extra(5))/dy)
      elseif(nxh.eq.nx.and.dyh.eq.dy.and.latmax.lt.head%extra(6)
     + .and.nyh.ge.ny)then
      print *,' Extra data are omitted'
      else
      stop 'Data files may be incompatible. hazallXL.v2 not ready '
      endif
      endif
      if(period1.ne.head%period)then
      write(6,66) head%period,period1
      nboob=nboob+1
      endif
66      format('Input file period mismatch. Now ',f4.2,' was ',f4.2,' s.',/,
     + '*** This is by most standards a serious error. ***')
c      write(6,*) head%nlev
      do i=1,nlev
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
      ndata= nrec*nlev
      if(nbad.gt.0)stop 'mismatch of sampled pga or sa.'
      if(nmal.gt.0)print *,'Continuing, should be ok for gm<',
     + head%xlev(19)
      print *,'ndata ',ndata,' nbig ',nbig
      if(nbig.lt.ndata)stop 'You must redimension nbig>ndata '
      call getbuf2(prob,ndata,readn)
      if(ndata.ne.readn)then
      write(6,*) ndata,readn
         nmin=min(readn,nmin)
         nshort=nshort+1
         if (readn.lt.0.8*ndata)then
         print *,'Warning: Input array size way short of expected'
         else
c sometimes the data in Southernmost regions is less important. Less severe message
         print *,'Input array size somewhat short of expected'
         endif
         endif
      write(6,5) "Enter weight for hazard curves: "
      read(infi,*) factor
      write(6,*)factor
      do 60 i=1,readn
      prob(i)=max(prob(i),1e-20)
 60   prob0(i)= prob0(i)+prob(i)*factor
c      do 61 i=1,ndata
c 61   prob0(i)= prob(i)
      headr%nlev= head%nlev
      do 95 i=1,nlev
 95   headr%xlev(i)= head%xlev(i)
       ifilm=min(6,ifil)
      headr%name(ifilm)= name
c      headr%extra(ifil)= factor      !this space is already rented out.
      headr%period=head%period
 11   continue
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
        call openwx(ifp,namout)
c 30 characters
c        ndata=308
c 128 characters
        ndata= 896
        call puthead(ifp,headr,ndata,readn)
        ndata= nrec*nlev
        call putbufx(ifp,prob0,ndata,readn)
        write(6,*) ndata,readn
         write(6,*) nrec,readn
         if(nboob.gt.0)print *,'Number of period mismatches ',nboob
        stop'binary file was written'
        else      !ascii output
        open(3,file=namout,status='unknown')
        write(3,32)nfile,nom_de_plume,head%period,nlev
32      format('#Pgm hazallXL.v2.f (harmsen) sums ',i2,' hazard curves from ',
     + a,/,'#Lat Long   Pex for spectral period ',f4.2,/,'#GM set(g)',i3)      
        write(3,39)(headr%xlev(i),i=1,nlev)
39      format(e11.5)
      k1=1
      k2=nlev
      do i=1,nrec
      xi=xmin + mod(i-1,nx)*dx
      yi=latmax - (i-1)/nx*dy
c the 2D 1D equivalence. Hats off to George Working Class Boole
      write(3,49)yi,xi,(prob0(j), j=k1,k2)
      k1=k1+nlev
      k2=k2+nlev
      enddo
49      format(f9.3,1x,f9.3,20(1x,e11.5))
      close(3)
      stop'Ascii hazard-curve file written'
      endif      !ii=0
      endif       !iout=1
      write(6,50) "Enter name of output file: "
      read(infi,900) namout
      print *,namout
      print 50,'Enter a 1 for ascii output; 0 for binary: '
50      format(a,$)
      read (infi,*)i
      print *,i
      ascii=i.eq.1
      if(i.eq.0)then
      call openwx(ifp,namout)
      else
      open(3,file=namout,status='unknown')
33      format('#Pgm hazallXL.v2.f (harmsen) sums ',i2,' hazard curves',
     + /,'#Long  Lat  SA(g) for spectral period ',f4.2,' Pex ',e11.5)      
      endif
      write(6,50) "Enter annual probability of exceedance: "
      read(infi,*) pex
      print *,pex
      if(ascii)write(3,33)nfile, head%period,pex
      pex= alog(pex)
      write(6,50) "Enter 0 to output prob gnd motion; 1 for Poisson prob: "
      read(infi,*) ip
      print *,ip
c change from ground-level index to actual ground level. may 2009.
      if(ip.eq.1) then
         write(6,50) "enter gnd motion in g and t (yrs): "
         read(infi,*) gm,t
         write(3,349)gm,t
349	format('#Long  Lat  PROB. of SA or PGA>',f6.3,' g in ',f5.1,' yrs.')
         write(6,*)gm, t
         do j=1,nlev
         xlev(j)=alog(xlev(j))
         enddo
      gml=alog(gm)
c log-log interpolation tends to be more accurate.
         j=2
         dowhile(gml.gt.xlev(j).and.j.lt.nlev)
         j=j+1
         enddo
         jm=j-1
      x1= xlev(jm)
      x2= xlev(j)
      frac1=(x2-gml)/(x2-x1)
      frac2=1.0-frac1
      
         do 100 i=1,nrec
         ind= (i-1)*nlev + jm
         indp=ind+1
c log(ann rate)
         rate=frac1*alog(prob0(ind))+frac2*alog(prob0(indp))
         rate=exp(rate)
c back to annual rate space.
      prx= 1.- exp(-rate*t)
      xi=xmin + mod(i-1,nx)*dx
      yi=latmax - (i-1)/nx*dy
100      write(3,49)xi,yi,prx
c         call openwx(ifp,namout)
c         ndata=308
c         headr%extra(9)= ilev
c         headr%extra(10)= t
c         call puthead(ifp,headr,ndata,ndata)
c         call putbufx(ifp,out,nrec,readn)
c         write(6,*) nrec,readn
         if(nboob.gt.0)print *,'Number of period mismatches ',nboob
         stop 'Ascii file written'
         
         endif
      nlev2= nlev-1
c--- interpolate to find ground motions with specified pex
      do 10 i=1,nrec
c      if(prob(i).eq.0.) then
c        out(i)= alog(xlev(1))
c        go to 10
c        endif
      do 20 j=1,nlev2
      ind= (i-1)*nlev +j
      if(prob0(ind).le.0.) then
        out(i)= 1e-20
        go to 20
        endif
      y1= alog(prob0(ind))
      y2= alog(prob0(ind+1))
      x1= alog(xlev(j))
      x2= alog(xlev(j+1))
      if((j.eq.1).and.(pex.gt.y1)) then
        out(i)= x1
        go to 10
        endif
      if((j.eq.nlev2).and.(pex.lt.y2)) then
        out(i)= x2
        go to 10
        endif
      if((y1.ge.pex).and.(y2.le.pex)) then
        tmp= y2-y1
        if(tmp.eq.0.) write(6,*) i,j,y2,y1
        out(i)= x1+ (x2-x1)*(pex-y1)/tmp
        go to 10
        endif
 20   continue
 10   continue
      write(6,50) "Enter scale factor: "
      read(infi,*) scale
      print *,scale
      ymax= 0.
      do 90 i=1,nrec
      if(out(i).ne.0.) out(i)= exp(out(i))
      if(out(i).eq.0.) out(i)=  xlev(1)
      out(i)= scale*(out(i))
c this is the standard order of data NW to SE. As long as this order is preserved
c more southerly parts of model will correspond to larger indexes. These can be ignored.
      xi=xmin + mod(i-1,nx)*dx
      yi=latmax - (i-1)/nx*dy
      if(ascii)write(3,36)xi,yi, out(i)
36      format(f9.3,1x,f7.3,1x,e11.5)
       if(out(i).gt.ymax)then
        ymax= out(i)
        xistar=xi
        yistar=yi
        endif
90      continue         
      write(6,800) ymax, xistar,yistar
      if(nshort.gt.0)print *,'Number of input files with nread<expected ',nshort
      if(nvschang.gt.0)print *,'Warning: inconsistent Vs30 among input files'
      if(ascii)stop
 800  format(' ymax=',f12.4,' at ',f7.2,1x,f6.2)
      headr%extra(8)= pex
      headr%extra(9)= scale
      headr%extra(10)= ymax
c---- output probabilistic ground motions
      if(i.eq.0)then
      print *,' binary output file name ',namout
c 30 characters
c      ndata=308
c 128 characters
      ndata= 896
      call puthead(ifp,headr,ndata,readn)
      call putbufx(ifp,out,nrec,readn)
      write(6,*) nrec,readn
      endif
 99   continue
       stop 'hazallXL.v2 normal end'
 1066      stop 'input file on command line not found'
      end
