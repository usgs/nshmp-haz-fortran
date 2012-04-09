c  hazpoint.f a Fortran-95 program for psha hazard curve preparation
c----program hazpoint:outputs hazard curve data with given weights at a loc.
c--- Good for cross-checking various programs' output. 
c--- from output of hazgridX and hazFX
c --- Program assumes that hazard curve data files are demultiplexed, 
c   i.e., site1's data, then
c --- site2's data, etc. Data are nlev mean ann. rates per site.
c --- Site grids may differ for various input files. Site should be within
c --- each grid
c -- To compile: f95 hazpoint.f iosubs.o -o hazpoint -e
c--- to run: hazpoint  CAsum1hz (or whatever control file)
c--- this version gets parameter file from header of 1st input file
c--- input files have hazard curves for each site concatenated
c---parfile= parameter file as used in hazgrid, etc. 
c === 2005: header may have important information in the "extra" fields.
c --- Grid location of sites might be contained in extra(2) to extra(7).
c --- Assumed Vs30 (m/s) may be in extra(9) (NGA relations all use Vs30)
c --- Sediment depth may be stored in extra(10) (Campbell-Boz. uses this)
c
      type  header 
        character*128 :: name(6)
        real*4 :: period
        integer*4 :: nlev
        real*4 :: xlev(20)
        real*4 :: extra(10)
      end type header
      type(header) :: head
	character date*8,time*10,zone*5
	integer ival(8)
      integer readn, nhead, idt(3)
	integer*8 nread8
      character*20 site
	integer*8 nlev, iskip
      real xls,yls
      logical grid, pgv
      logical empty_hd
c header can have useful information in extra(). As of
c Oct 20+ 2005 extra(9) contains the vs30 for the run. User should be
c aware of the any variation in vs30 for diff. input files. This
c information, if available, is broadcast on standard out (terminal). However,
c older files have basically empty header in which case empty_head = .true.
c
      logical ok
      dimension prob(100),erate(100)
c erate and prob will differ based on treating weights as aleatory vs epistemic
c respectively. SH Feb 7 2005. erate=estimated rate
c 
      character name*128,namepar*30
c
      if(iargc().lt.1)stop'Put parameter file on command line'
       	write(6,50)'Enter site lat, long (10ths of a degree):'
       read *,yls,xls
      write(6,50)'Enter site name, 2<= string < 20 char: '
      read 5,site
      ms=max(2,index(site,' ')-1)
      site=site(1:ms)//'.'
      ms=ms+1
c       call idate(idt)
        call date_and_time(date,time,zone,ival)
       if(yls.lt.-80.)then
c if Antarctica could be what the dude really wanted. 
       write(6,*)'Switching your so-called lat and long...'
       x=yls
       yls=xls
       xls=x
       endif
c     read(5,900) namepar
      call getarg(1,namepar)
      inquire(file=namepar,exist=ok)
      if(.not.ok)then
      write(6,*)'Command file needed: ',namepar
      write(6,*)'Put in working dir and retry this prog.'
      stop
      endif
	i=index(namepar,' ')-1
       open(3,file=site(1:ms)//namepar(1:i),status='unknown')
       write(6,*)' output will go to file named site.'//namepar
       write(3,100)xls,yls,date,namepar
 100      format('#site coords ',2f9.2,' Program: hazpoint.f (Harmsen) ',
     + /,'# Extracts PSHA spectral values. Run on ',a,
     + ' instructions from ',a20)
      open(unit=1,file=namepar,status='old')
      read(1,*,err=2006) nfile
      do 11 ifil=1,nfile
      read(1,900) name
      read(1,*) factor
 900  format(a)
      inquire(file=name,exist=ok)
      if(.not.ok)then
      write(6,*)'File missing ',name
      stop 'should be in working dir.'
      endif
      call openr(name)
      nhead= 896
      call gethead(head,nhead,readn)
c------ get map parameters. From ascii file or from header
      empty_hd=head%extra(2).eq.0..and.head%extra(6).eq.0.
      if(ifil.eq.1.and.empty_hd) then
      inquire(file=head%name(1),exist=ok)
      if(.not.ok)then
      write(6,*)'file needed: ',head%name(1)
      write(6,*)'Put in working dir and retry this prog.'
      stop
      else
c      write(6,*)'File found for parameters: ',head%name(1)
         open(unit=2,file=head%name(1),status='old')
         read(2,*)nsta
         grid=nsta.eq.0
         if(nsta.gt.20)then
         rewind(2)
         endif
         read(2,*) ymin,ymax,dy
        read(2,*) xmin,xmax,dx
      write(6,*)'Presumed region from which site info will be drawn:'
      write(6,*)ymin,ymax,dy,' lat range &dlat'
      write(6,*)xmin,xmax,dx,' long range &dlong'
        close(2)
        endif
        else
c almost always want data from header        
        dx=head%extra(4)
        dy=head%extra(7)
        xmin=head%extra(2)
        xmax=head%extra(3)
        ymin=head%extra(5)
        ymax=head%extra(6)
        grid=dx.gt.0..and.dy.gt.0.
        endif !ifil .eq. 1
        nx= nint((xmax-xmin)/dx) +1
        ny= nint((ymax-ymin)/dy) +1
        nrec= nx*ny
        nxs=nint((xls-xmin)/dx) + 1
        nys= nint((ymax-yls)/dy)
        nsite=nys*nx+nxs
	write(6,*)nxs,nys,ymax,xmin,dy,dx
        write(6,*)nx,ny,nsite,' nx  ny site_index'
        if(nxs.lt.0.or.nys.lt.0)stop'site outside of study area'
      if(nxs.gt.nx)stop'Your site longitude is east of input-file region'
      if(nsite.gt.nrec)stop 'site outside of study area'
c there are other bad cases but we only catch the above.
c
c I encourage colleagues to use nint in the below type calcs.
c Errors arise when you dont use nint for several xmax   xmin combinations 
c      xmax-xmin ending in .1, .2, .6, .7 are known examples on SUN systems. 
c      ditto for ymax-ymin
cccc
      write(6,909) 'input file: ',head%name(1)
      write(6,909) 'program: ',head%name(3)
      if(head%name(2)(1:1).ne.' ')
     + write(6,909) 'attn model OR agrid file: ',head%name(2)
909      format(a,1x,a)
      if(head%period.gt.0)then
      write(6,*) head%period,' spectral period s'
      pgv=.false.
      elseif(head%period.eq.0.)then
      write(6,*) 'header period indicates PGA'
      else
      write(6,*)' Peak ground velocity cm/s'
      pgv=.true.
      endif
      if(head%extra(9).ne.0.)then
c storing vital info in header records begins Oct 2005 (harmsen only)
      vs30=head%extra(9)
      write(6,*)'Vs30 for above file is ',vs30,' m/s'
      write(6,*)'  Sediment depth is ',head%extra(10),' km.'
      else
      vs30=760.	!assume bc boundary condition as default
      endif
c information on dtor and atten models:
      do i=4,6
      if(head%name(i).ne.'')write(6,*)head%name(i)
      enddo
      if(xmin.ne.head%extra(2).or.xmax.ne.head%extra(3))
     +write(6,*)'First Xmin,xmax ',xmin,xmax,' d. Current ',
     +head%extra(2),head%extra(3)
       nlev= head%nlev
       iskip= 4*(nsite-1)*nlev+nhead
c Skip (bytes) to the site whose location was given. From beginning of
c file, thus header as well as data need to be skipped.
c
c For skipping (lseek) the less frequently used getbuf3 is appropriate here
c
c      write(6,*)ifil,nlev
      do m=1,nlev
      erate(m)=1e-21
      enddo
      call getbuf3(erate,nlev,nread8,iskip)
      if(nread8.lt.nlev)write(6,*)'Read failed ',readn,' expected ',nlev
      call closeio(name)
      inm=index(name,'  ')
      if(pgv)then
      write(3,122)factor,name(1:inm),vs30
122      format(/,'#PGV(cm/s)  rate_exc  weight: ',f5.3,' for ',a,
     + ' vs30 ',f6.1,' m/s')
      else
      write(3,120)factor,name(1:inm),vs30
120      format(/,'#SA(g)  rate_exc  weight: ',f5.3,' for ',a,
     + ' vs30 ',f6.1,' m/s')
      endif
      np=0      
      do 95 i=1,nlev
      rate=erate(i)*factor
        if(rate.gt.0.6e-11 .and..not.pgv)then
        write(3,123) head%xlev(i),rate
        np=np+1
        elseif(rate.gt.0.6e-11)then
        write(3,126)head%xlev(i),rate
        np=np+1
        endif
5      format(a)
50      format(a,$)
 95      continue

123      format(f9.6,1x,e12.5)
126      format(f9.4,1x,e12.5)
      if(np.eq.0)write(6,*)'No hazard curve for these site coordinates. Check file status'
 11   continue
       close(3)
       write(6,*)'Note: data w/ weight*rate < 10**-10 have been omitted'
      stop
2006      write(6,*)'instruction file must begin with an integer, number of files'
      stop 'hazpoint: please check your command line input file'
      end
