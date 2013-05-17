c-----fltrate.2013.f: calculates rates for characteristic EQ's and
c---- incremental a-values, given slip rates and points of faults
c --- Correction feb 23 2010: only increase sliprate on dipping normal or reverse flts,
c --- not on dipping Strike-slip faults. 
c Jan 11 2013. Add branches on dip uncertainty which preserve event rates but vary moment
c rate because of different areas.
c Compiile on Solaris:
c f95 fltrate.2013.f -e -o ../bin/fltrate.2013
c gnu fortran:
c  gfortran fltrate.2013.f -ffixed-line-length-none -finit-local-zero -ffpe-trap= -o fltrate.2013
c The finit-local-zero initializes variables to 0 (logicals to .false.). THis is
c an essential flag when compiling this code on PCs. (Thanks Jennifer Haase)
c However, this flag is not recognized on the Golden cluster.
c --- repair feb 21 2007: magmin for char. with floating cmag needs to be same
c as char mag. Otherwise you will get two gr distributions for those faults.
c Look at faults like steens, sevier/toroweap, and several others.
c --- 
c --- following Advisory Panel meeting of May 3-4, make min mag for NV faults 6.5
c --- check what happens
c --- mod may 4 2007 S Harmsen
c Mod feb 28, include another output file with this summary info:
c ID,name,M_char,Recurr_char,weight,A_GR,B_GR,M_min,MMax,Delta_M
c from Haller request. 2/26. This file is called table.tmp
c----- v1: this version used for May 2006 sensitivity study maps, reported
c to IMW workshop in Reno May 31-June 1 2006. From Frankel fltrate
c --- Input differs from Frankel's in that an additional flag, img, is input.
c --- img is set to 0 to indicate that input Mmax is to be used rather than
c --- the computed value. Geologic constraints historical events,...
c	A global Mmax can override any computed value. We often use 7.5
c	to cap BRP Mmax for long faults such as steens.
c Mod june 7. Different Global Mmax for SS and Dip Slip
c
c --- v2: this version does not round M to nearest tenth but to
c --- nearest .01. Tenths may be too coarse if we wish to determine epistemic
c --- variability due to such things as 10% length uncertainty. 
c shift definition corrected Dec 4 2006.
c ---- v2 affects GR calculations as well as char. Testing May 31 06.
c --- v2 uses event rate below for Bear River and Eglington flts to save
c --- user the effort (memory lapse?) of editing in these values.
c --- Added Oct 12 2006: Last instruction, img, specifies way to calculate
c --- M  from other information. Can use M(L) for normal slip as well, by
c --- making img = 3. This needs further work. In the output record
c --- having to do with iflt, add a 1 (number of mags) and img. Img instruction
c --- is needed by hazFX so that when GR relation is done, same L is used
c --- when floating the ruptures. L= A(M)/downdip width. hazFXnga5: not yet using
c --- this img instruction as of Oct 13 2006.
c ---
c---- All versions of fltrate write char and G-R files for hazFX
c ---- mod may 3 2006: if slip rate < 0 dont use this fault but report
c ---- in the log file. 
c --- mod may 22: for normal slip sources, may compute M(area) instead of M(length) to
c --- study effect of fault dip on rate-of-event estimates and other hazard-related
c --- stuff.
c ---- always report Mmax on input compared to Mmax(L) computed here in the
c ---- log file, but use Mmax of input in the calcs and output. (problem found
c ---- with Quaternary db has Mmax always but is this the one we want?). Emails. 
c --- open item.
c---- Base slip rate onn uplift rate for BRP faults.
c---- magnitude from fault length generally, go to rupture area for normal.
c---- shifts magmin and magmax by dmag/2. for GR
c---- uses seismogenic depth and dip to calculate width
c---- makes separate file for characteristic M<= 6.5
c--- for char M>7.5 makes Mmax=7.5 and uses floating M7.5 for characteristic
c--- modified to accept input cmagmax instead of M7.5
c---- Note: No longer rounds off char M to nearest 0.1 before rate calculation.
c -- writes a message if strike changes by more than 90d on adjacent segments
      implicit none
      character*50 namein,nameout,nameflt,nameout2,namelog,nameout3
      character*50 coord
      character*30 tname(100,3)
      character*3 blank
      logical bogus/.false./,bear/.false./,egg/.false./,fok/.false./
c      logical bogus,bear,egg,fok
      real, dimension (1000):: x,y
      real bear_r/4.35e-4/,egg_r/7.1e-5/
      real fac, dmag, dmaggr, rate, rate1, rate2, rate3, rategr, rate65
      real dip, shift, sdep, ry, rx, rategr3, rategr2, dipo, diff, depth, delta
      real d, d2r, ctime, crate, charmin, bot, baz, b, avglat, arate, allmo, a, a2
      real ymag, ybar, xmu, xmo, xmagl, xmag, xlen, xln, wt, widtho, width, trate, totmo, top, sx, sy
      real sumrate2, sum2, sum1, avgaz, wtin
      integer ic/0/
      integer iuplift, ndip, itype, imax, img, ichar 
c ic = number of faults where mmax is determined from area or length
      real magmax,magmin,magmax0,morate,morateo
c add hslip at each 0.25 d increment in latitude. this is the summed horizotal slip rate.
      real cmag,cmagmax,az,azo,hslip(100),wlist(100)
      real, dimension(5):: dip_br,wt_dip
      integer list(100)
      real mmax(3)
 900  format(a)
 50      format(a,$)
 	list=0
 	wlist=0.0
 	hslip=0.0
 	ilatmn=100
 	ilatmx=0
       if(iargc().gt.0)then
       call getarg(1,namein)
       else
1      write(6,50) "Enter name of input file: "
      read(5,900) namein
      inquire(file=namein,exist=fok)
      if(.not.fok)then
      print *,'File not found. Try again'
      goto 1
      endif
      endif
      open(unit=1,file=namein,status='old')
      write(6,50) "Enter 1 for latlon, 0 for lonlat: "
      read(5,*) latlon
      write(6,50) "enter 1 for uplift rate input: "
      read(5,*) iuplift
      write(6,900) "Name of output file for log is logjam.log "
      namelog='logjam.log'
      write(6,50)'Enter number of dip uncert branches: '
      read *,ndip
	if(ndip.gt.1)then
	      write(6,900)'Enter list of deltadip and weight: '
      do i=1,ndip
      write(6,56)i
56	format('For branch ',i1,': ',$)
	read *,dip_br(i),wt_dip(i)
	enddo
	else
	dip_br(1)=0.0
	wt_dip(1)=1.0
	endif
      open(unit=2,file=namelog,status='unknown')
      open(10,file='normal.flt.tmp',status='unknown')
      print *,'Data on normal fault Mmax going to normal.flt.tmp'
      write(10,1009)
c1009      format('#Name fault',t25,'Xlen width Mmax(area)')
1009      format('#Xlen width Mmax(area) ID Name of fault')
      write(6,*) "enter name of output file for char."
      read(5,900) nameout
      open(unit=3,file=nameout,status='unknown')
      write(6,*) "enter name of output file for G-R"
      read(5,900) nameout2
      open(unit=4,file=nameout2,status='unknown')
      write(6,*) "enter name of output file for M<6.5"
      read(5,900) nameout3
      open(unit=8,file=nameout3,status='unknown')
      write(6,*) "Enter seismogenic thickness (km): "
      read(5,*) sdep
      writE(6,50) "enter max mags for characteristic SS, DS resp: "
      read(5,*) mmax(1),mmax(2)
      mmax(3) = mmax(2)
      print *,'Mmax array ',mmax
      write(6,*)'ID,name,M_char,Recurr_char,weight,A_GR,B_GR,M_min,MMax,Delta_M'
      write(6,*)'will go to a file called TABLE.TMP'
      open(20,file='TABLE.TMP',status='unknown')
      write(20,900)'# DATA FROM '//namein
      write(20,900)'#ID,name      Length(km) M_char  Recurr_char Dip weight  A_GR   B_GR  M_min,MMax,Delta_M'
       read(1,900)nameflt
      dowhile(nameflt(1:1).eq.'#')
      print 900,nameflt
c      write(20,900)nameflt
      read(1,900)nameflt
      enddo
      backspace(1)
      write(6,*) "enter bval,magmin,dmag"
      read(1,*) b,magmin,dmag
      if(magmin.lt.5.05)then
      print *,'magmin is ',magmin,' looks unreasonable for USGS PSHA'
      stop 'fltrate.2013: please check input'
      endif
c      write(3,*) b,magmin,dmag,summag
c--above for G-R calculation of a-value only
c
      blank='   '
      sum1=0.
      sum2=0.
      icharv=1
      igr=2
209      write(6,50)'Enter imag 0 for fixed mag, 1 for calculated: '
      read *,img
      write(6,50)'Enter minimum mag for char events (6.5) : '
      read *,charmin
      if(img.ne.0.and.img.ne.1)goto 209
      rate1=0.
      rate2=0.
      rate3=0.
      rategr= 0.
      rategr2= 0.
      rategr3= 0.
      rate65=0.
      allmo= 0.
      ncnt=0
      d2r= 3.14159/180.
      do 1000 iflt=1,1000
      read(1,900,end=999)nameflt
      egg=index(nameflt,'Eglington').gt.0
c      bear=index(nameflt,'Bear Riv').gt.0
	print 900,nameflt
c      write(6,*) "enter slip rate (mm/yr),type,dip,width,nmg,magmax0,wt"
      read(1,*,end=999,err=1999) rate,itype,xln,dip,top,bot,width,
     + nmg,magmax0,wtin      
c img is new, a flag to determin how mmax is computed
      cmagmax=mmax(itype)
      print *,'cmagmax now ',cmagmax
c rate above can be slip rate (if iuplift .ne. 1) or 
c rate can be uplift rate (if iuplift .eq.1)      just to keep you on edge of seat
      bogus=rate.le.0.0      !some faults have negative or zero rates these will be
      if(bogus)then
      write(6,*)'Please enter event rate for ',nameflt
      read *,trate
      if(trate.gt.0.)bogus=.false.
      endif
c identified in log file but skipped in output    
      if(img.ne.0.and.img.ne.1.and.img.ne.3)
     +  stop'fltrate: this version requires an img flag after mmax0'
      if(bear)rate=bear_r
c some faults may have different "sdep" than the regional value entered,
c example Calaveras C with sdepth 10 km but regional sdep might be 15 km. Below
c effectively kills this fault-specific information. SH may 3 2006.      
      if(sdep.ne.0.) width= sdep/sin(abs(dip)*d2r)
      if((iuplift.eq.1).and.(dip.lt.89.).and.itype.gt.1) then
        rateorig= rate
        rate= rate/sin(abs(dip)*d2r)
c        write(*,*)'rateOrig/rate: ',rateorig,rate
      endif
c--itype= 1 for strike slip; 2 for reverse; 3 for normal
c    write(6,*) "enter number of segment points"
       write(2,*) nameflt
       write(2,*) rate,width,itype
       read(1,*) npts
      if(latlon.eq.1)then
      write(6,*) "enter lat,lon for each point"
      else
      write(6,*)'Enter long,lat for each point'
      endif
      imax=0
      xlen= 0.
      avglat=0.
      do 2 j=1,npts
      read(1,900,end=99) coord 
      if(latlon.eq.0) read(coord,*,err=99) x(j),y(j)
      if(latlon.eq.1) read(coord,*,err=99) y(j),x(j)
      if(j.gt.1) then
        sy= y(j-1)
        sx= x(j-1)
        ry= y(j)
        rx= x(j)
        ybar=0.5*(y(j-1)+y(j))
        call delaz(sy,sx,ry,rx,delta,az,baz)
        xlen= xlen+ delta
        avglat=avglat+ybar*delta
        avgaz=avgaz+az*delta
        if(delta.gt.100.)print *,sx,sy,rx,ry,' adjacent x.y on ',
     + nameflt
c      if(tule)print *,'Length tule spr rim flt now ',xlen
	if(az.lt.0.)az=az+360.
       if(j.gt.2)then
c QA check on fault strike. Big change is flagged but its occurrence does not shut down process.
c compare strike with previous to see if change is reasonable. Dont flag 10 vs 340 for example.
	diff=abs(az-azo)
	if(diff.gt.90.and.diff.lt.270.)print *,' *** Rapid change in strike ',diff,' at j ',j,' *****'
	endif	!j>2
        azo=az
        endif
   2  continue
  99  continue
        if(.not.bogus)then
        xlen=max(1.,xlen)
      write(6,*) xlen,avglat/xlen
c--- for characteristic model, find cmag, crate
c make cmag == charmin if less than charmin. 
c New mod of may 4
       xmagl= 5.08+1.16*alog10(xlen)
      xmag= xmagl
      if(xmag.lt.charmin)then
      write(2,303)xmag,charmin
303      format('xmag from len ',f6.3,' reset to ',f6.3)
      xmag=charmin
      endif
c      write(10,1001)xlen,width,xmag
1001      format(1x,f6.2,1x,f5.2,1x,f5.2,1x,a30)
c Here is the primary v2 change: discretize Mmax to 0.01
      xmag= xmag*100.
      imag= int(xmag+0.5)
      xmag= float(imag)/100.
      if(itype.eq.3)then
      write(10,1001)xlen,width,xmag,nameflt
      avglat=avglat/xlen 
      ilat=(avglat-31.)*4+1
c      print *,avglat,ilat,xlen,nameflt
      ilatmx=max(ilatmx,ilat)
      ilatmn=min(ilatmn,ilat)
      hslip(ilat)=hslip(ilat)+rate/tan(dip*d2r)*wt	!wt can be less than 1: diff epistemic br
      list(ilat)=list(ilat)+1;wlist(ilat)=wlist(ilat)+wt
      if(list(ilat).lt.4)tname(ilat,list(ilat))=nameflt(1:30)
      elseif(itype.eq.1)then
c EW extension from SS movement. should be discussed.
      avglat=avglat/xlen
      ilat=(avglat-31.)*4+1
      avgaz=avgaz/xlen
c      print *,avglat,ilat,xlen,nameflt
      ilatmx=max(ilatmx,ilat)
      ilatmn=min(ilatmn,ilat)
	hslip(ilat)=hslip(ilat)+rate*abs(sin(avgaz*d2r))*wt      
      list(ilat)=list(ilat)+1;wlist(ilat)=wlist(ilat)+wt
      endif
c--- magmax0 ne 0 for floating characteristic EQ
      if(magmax0.gt.0..and.(img.eq.0.or.nmg.eq.0)) then
        xmag=magmax0
        imax=0      !if 1 does a gr for char. this seems wrong. already have gr.
      else
      ic=ic+1
      magmax0=min(xmag,cmagmax)	!redefine according to current xmag
c      magmax0=min(xmag,mmax(itype))	!redefine according to current xmag
      imax=0		!this wasn't in Frankel code. Some big 'uns got thru
        endif
      if((xmag.gt.cmagmax).and.(imax.ne.1)) then
        xmag=cmagmax
        imax=1
        ymag=cmagmax      !new feb21. magmin also needs to be fixed for char.
        else
        ymag=xmag
        endif
      xmu= 3.0e10
      xmo= 1.5*xmag + 9.05
      xmo= 10.**xmo
      rate= rate/1000.
      xlen= xlen*1000.
      width= width*1000.
      morate= xmu*width*xlen*rate
      if(morate.lt.0..and.trate.gt.0.)then
      crate=trate
      morate=crate*xmo
      print *,' hand edited crate and morate '
      else
      crate= morate/xmo
      if(egg)then
      print *,'Eglington may need a special branch hand edit'
      print *,'Computed ann rate is ',crate,' 2002 val ',egg_r
      write(2,*)'Eglington may need a special branch hand edit'
      write(2,*)'Computed ann rate is ',crate,' 2002 val ',egg_r
c      crate=egg_r	!these are edited-in values
      endif	!eglington dawg 
      endif
      if(xmag.gt.6.5) rate1=rate1+crate
      if(xmag.ge.7.0) rate2=rate2+crate
      if(xmag.ge.7.5) rate3=rate3+crate
      if(xmag.eq.6.5) rate65= rate65+ crate
      xlen= xlen/1000.
      cmag=xmag
      write(6,*) xlen,cmag,crate
      magmax= xmag
c---- for floating characteristic EQ; nmag=1
      if(imax.eq.1) then
        xmo= 9.05+1.5*magmax
        xmo= 10.**xmo
        rate= morate/xmo
        a2= alog10(rate)+b*magmax
        endif
ccc
c------for Gutenberg-Richter,  find a-value
      ineq=0
      if(magmax0.ne.0.) magmax=magmax0
c--- shift magmin and magmax to center mag. bins. 
c v2: because magmax may be less than magmin+dmag, we need a better shift value
      if(magmax.ne.magmin) then
      d=magmax-magmin
      if(d.le.0.1)then
      shift=0.5*d
      nmag=1
      elseif(d.le.0.3)then
      nmag=2
      shift=d/4.
      elseif(d.le.0.5)then
      nmag=4
      shift=d/8.
      elseif(d.le.0.7)then
      nmag=6
      shift=d/12.
      else
      nmag=8
      shift=d/16.
      endif
         magmin= magmin+shift
         magmax= magmax-shift
         dmaggr=shift*2.0      ! dmaggr can be smaller than dmag input
         ineq=1
         else
         nmag=1
         endif
c      nmag= (magmax-magmin)/dmaggr +1.4      XXX use above
      totmo= 0.
      do 10 i=1,nmag
      xmag= magmin+(i-1)*dmaggr
      xmo= 1.5*xmag+9.05
      xmo =10.**xmo
      rate= -b*xmag
      rate= 10.**rate
  10  totmo= totmo+xmo*rate
      rate= morate/totmo
      a= alog10(rate)
      write(6,*) nmag
      totmo= 0.
      sumrate2=0.
      do 11 i=1,nmag
      xmag= magmin+(i-1)*dmaggr
      xmo= 1.5*xmag+9.05
      xmo =10.**xmo
      rate= a- b*xmag
      rate=10.**rate
      if(xmag.ge.6.5) rategr= rategr+rate
      if(xmag.ge.7.0) rategr2= rategr2 + rate
      if(xmag.ge.7.5) rategr3= rategr3 + rate
  11  totmo= totmo+ xmo*rate
      sum1= sum1+totmo
      sum2= sum2+ morate
      allmo= allmo+morate
      write(6,*) morate,totmo
c----"a" is incremental a-value  for dmag interval
c-----write stuff to output files
      width= width/1000.
      icharv=1
      igr=2
c--- default min depth
      depth= 0.
c--- re-shift magmin and magmax for output: there is a v2 change here
c from dmag/2 to shift, defined above.
      if(ineq.eq.1) then
        magmin= magmin-shift
        magmax= magmax+shift
        endif
      if((cmag.gt.6.5).and.(imax.eq.0)) then
      do ka=1,ndip
      wt=wt_dip(ka)*wtin
      dipo=dip+dip_br(ka)
      fac=sin(dip*d2r)/sin(dipo*d2r)
      morateo=morate*fac
      widtho=width*fac
      write(3,*) icharv,itype,1,img,blank,nameflt
      write(3,233) cmag,crate,wt,morateo
 233  format(f6.3,1x,f9.7,1x,f6.4,1x,e11.5)
      write(3,*) dipo, widtho, top, xlen
      write(3,*) npts
      write(20,202)nameflt(1:20),xlen,cmag,1./crate,dipo,wt,a,b,magmin,magmax,dmaggr
      do 100 i=1,npts
 100  write(3,905) y(i),x(i)
       enddo	!ka
 905  format(f14.5,1x,f14.5)
      endif
      if((cmag.gt.6.5).and.(imax.eq.1)) then
      do ka=1,ndip
      wt=wt_dip(ka)*wtin
      dipo=dip+dip_br(ka)
      print *,wt,dipo,dip,sin(dipo*d2r)
      fac=sin(dip*d2r)/sin(dipo*d2r)
      morateo=morate*fac
      widtho=width*fac
      write(3,*) igr,itype,1,img,blank,nameflt
      write(3,*) a2,b,ymag,magmax,dmag,wt,morateo      !ymag replaces magmin for char
      arate=10**(a2-b*magmax)
      write(20,202)nameflt(1:20),xlen,magmax,1./arate,dipo,wt,a,b,magmin,magmax,dmaggr
202      format(a20,1x,f6.1,1x,f5.2,1x,f8.1,1x,f6.1,
     +1x,f5.3,1x,f8.5,1x,f5.2,1x,f5.2,1x,f5.2,1x,f6.4)
      write(3,*) dipo, widtho, top, xlen
      write(3,*) npts
      do 101 i=1,npts
 101  write(3,905) y(i),x(i)
      enddo
      endif
      if(cmag.gt.6.5) then
      do ka=1,ndip
      wt=wt_dip(ka)*wtin
      dipo=dip+dip_br(ka)
      fac=sin(dip*d2r)/sin(dipo*d2r)
      morateo=morate*fac
      widtho=width*fac
      write(4,*) igr,itype,1,img,blank,nameflt
c v2 change: try writing dmaggr as the delta-M. Does this work in hazFX? Yes,
c for the nga versions such as hazFXnga5. No for the older versions which 
c precompute tables of Pex.
      write(4,444) a,b,magmin,magmax,dmaggr,wt,morateo
444     format(f8.5,1x,f5.2,1x,f6.3,1x,f6.3,1x,f6.4,1x,f7.5,1x,e11.5)
      write(4,*) dipo, widtho, depth, xlen
      write(4,*) npts
      do 106 i=1,npts
 106  write(4,905) y(i),x(i)
 	enddo	!dip branches new jan 11 2013.
        endif
      if(cmag.le.6.5) then
      do ka=1,ndip
      wt=wt_dip(ka)*wtin
      dipo=dip+dip_br(ka)
      fac=sin(dip*d2r)/sin(dipo*d2r)
      morateo=morate*fac
      widtho=width*fac
      write(8,*) icharv,itype,1,img,blank,nameflt
      write(8,*) cmag,crate,wt,morateo
      write(8,*) dipo, widtho, depth, xlen
      write(20,202)nameflt(1:20),xlen,cmag,1./crate,dipo,wt,0.,b,0.,0.,0.
      write(8,*) npts
      do 105 i=1,npts
 105  write(8,905) y(i),x(i)
       enddo	!do ka loop over dip uncert.
        endif
      ctime= 1./crate
      do ka=1,ndip
      wt=wt_dip(ka)*wtin
      dipo=dip+dip_br(ka)
      fac=sin(dip*d2r)/sin(dipo*d2r)
      morateo=morate*fac
      write(2,*) cmag,crate,ctime,xlen,wt,dipo
      write(2,*) a,b,magmin,magmax,dmag,morateo,xmagl
      enddo	!ka branches.
      else
      write(2,*)' this fault was skipped negative slip rate'
        endif      !if not bogus
 1000 continue     
1999      print *,' LOOK AT NEXT LINE !'
      print *,'err reading file current rate itype ',rate,itype
      print *,'#!@# FOO FOO  PLEASE READ ABOVE NOTE FOO FOO #!@#'
 999  continue
      ncnt= iflt-1
      write(6,*) ncnt
      write(2,*) ncnt
      write(6,*) "char rate > 6.5", rate1
      write(6,*) "char rate >= 7.0", rate2
      write(6,*) "char rate >= 7.5", rate3
      write(6,*) "gr rate >= 6.5=",rategr
      write(6,*) "gr rate >= 7.0=",rategr2
      write(6,*) "gr rate >= 7.5=",rategr3
      write(6,*) "M6.5 rate=",rate65
      sum1= 0.5*rate1 + 0.5*rategr + rate65
      write(6,*) "weighted rate M>=6.5", sum1
      sum1= 0.5*rate2 + 0.5*rategr2
      write(6,*) "weighted rate M>=7.0", sum1
      write(6,*) allmo
      write(2,*)'number of faults with Mmax from area or length ',ic
      close(2)
      open(2,file='horiz.sliprate.dat',status='unknown')
      write(2,400)namein
400	format('#fltrate.2013 determination of horizontal slip rate from file ',a30,
     +/,'#Lat_d hsrate(mm/yr) #flts Weight  names...')
      do i=1,ilatmx
      if(list(i).gt.0)write(2,401)31.+(i-1)*0.25,hslip(i),list(i),wlist(i),
     + tname(i,1),tname(i,2),tname(i,3)
401	format(f6.2,1x,f8.5,6x,i3,1x,f5.2,3(1x,a30))
      enddo
      close(2)
      print *,'look at file horiz.sliprate.dat'
      end
c
c
      subroutine delaz(sorlat,sorlon,stnlat,stnlon,delta,az,baz)
      d2r= 3.14159/180.
      xlat= sorlat*d2r
      xlon= sorlon*d2r
      st0= cos(xlat)
      ct0= sin(xlat) 
      phi0= xlon
      xlat= stnlat*d2r
      xlon= stnlon*d2r
      ct1= sin(xlat)
      st1= cos(xlat)
      sdlon= sin(xlon-phi0)
      cdlon= cos(xlon-phi0)
      cdelt= st0*st1*cdlon+ct0*ct1
      x= st0*ct1-st1*ct0*cdlon
      y= st1*sdlon
      sdelt= sqrt(x*x+y*y)
      delta= atan2(sdelt,cdelt)
      delta= delta/d2r
      az= atan2(y,x)
      az= az/d2r
      x= st1*ct0-st0*ct1*cdlon
      y= -sdlon*st0
      baz= atan2(y,x)
      baz= baz/d2r
      delta= delta*111.2
      return
      end         
