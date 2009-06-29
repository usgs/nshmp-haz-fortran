c sum_haz_2008b.f the WWW  version. Nov 10 2008. minor revs to initialization
c --- this pgm is prelim for 2008 update; it does this:
c --- (1) finds marginal M distribution for 6 levels of PGA (could use other Sa)
c --- (2)  makes plots and a report that should help to perform
c  liquefaction analysis
c     from sum_haz_2008.f
cc Primary purpose: sums w(i)*haz(i,j), i=1,nfil, j=1,...,nlev
c   in (distance, magnitude, epsilon) [ M] also R|M and eps0|M bins. nfil=number of files
c   of PSHA data, can be up to 80? (or modify parameter nfil). haz() is mean
c   annual rate of ground motion exceedances. w(i) is pre-assigned weight.
c   Nsig is number of epsilon bins, where sig represents aleatory and possibly
c   epistemic GM variability. (epistemic if more than one attenuation model used)
c Second purpose: generates a gmt script (not prepared initially)
c PC Linux compile (In 2008, this is gldplone): 
c gfortran ./sum_haz_2008b.f -static -O  -o ./sum_haz_2008b -ffixed-line-length-none
c
c Solaris compile: 
c f77 sum_haz_2008b.f  -e -fast  iosubs.o  -o sum_haz_2008b
c  note: -Bstatic , Solaris System objects to -Bstatic. the
c /usr/lib/libmvec.so is shared. 
c
c The -static argument invokes static library linking, needed for Web computers.
c
c Usage: this program is called from a script written by combine_2008b
c
c arg1 = site name (no blanks!)
c arg2 = site longitude (E+, dec degrees)
c arg3 = site latitiude (dec degrees)
c arg4 = site Vs30 (m/s)
c arg5 -10 = Six PGA levels, units: g
c arg? = computed annual rate of exceedances (exponential notation please)
c arg? = pid (process id, for names)
c arg? = absolute path name of data files  Out
c 
c When OKed, sum_haz_2008b is for automated web-site PSHA deaggregations.
c Input Magnitudes assumed to be greater than 4.6 (Mw). For lower magnitudes,
c some code modifications are required. The 4.6 limit corresponds to a 5.0(MbLg)
c limit in the 2002 National Seismic Hazard Maps. .
c Two output files: name specified and gmt.period (period is 4 char SA period)
c Steve Harmsen updated 11/2008
c
      parameter ( deltam=0.2,nsrcmx=35,tiny=1.e-11)
c We want M=5.2 to go in next higher bin than M=5.19. See amagmin below.
c
      parameter (nfmx=80,nmax=10000, nlev=6,nmmax=50)
c nmax=max number of records per file. Be careful here!
      logical avgctr,epson,getwt,header,ispga,aok,bypass
      logical fltok,ice
      real rbarr,mbarr
      real d,eps,m,mag,m3d,hmax/0.0/,wt
      real rmax/10./,deltar/25./,h_int,h_intm/0.0/
      character*32 name1,name2,a32
      character*1 ch1
      character*10 dateis,avs30
      character*7 faults/'faults/'/
      character*14 e3dmsg
      character*16 comm/'2008 edition. '/
      character*12 pid
      character*10 psa,rate, slong,slat,metric*20
      character*4 labs,label,per,retrnt*5,units/'yrs.'/,anum      !temp hold record fault number
      character*24 msg
      character*6 TSA,SA5*8
      common/lab/msg,label,per,retrnt,units,comm
      common/vs30/vs30
      common/pid/pid,ipid
      common/loca/loc,iloc
      common/lab2/psa,rate,xlong,ylat,metric
       common/lowlim/aratemn,ratett
      common/e3d/r3d,m3d,e3dmsg
       character*4 epr(8)
     1/'3.00','2.00','1.00','0.50','0.30','0.20','0.10','0.00'/  
      character*80 loc,rec,logdir
      character*32 fltnm(40)
c fltnm will hold short descriptive message (Char or GR; dip of fault)
c about significant faults.
       character*1 cseis,flg*2,flg3*3
      logical maybe
      real  amagmin/4.5999/
      real amag(nlev),rsav(nlev),esav(nlev)
      integer, dimension(nlev):: ireturn
      integer iloc,ireg,kf,kfmin,kfmax,kflim/1/,nocc,nvac
	real, dimension(nlev) :: avgr,avgm,ebar0,ratesum,ratepc,pga,arate,return
c fpr this analysis h (rate) is a scalar conditional  rate of exceed. for
c a given source & atten model and pga level.
      real h
      real hsum(nlev),emmin/10./,emmax/1./
c for each gm level find (M,f) distribution and conditional rbar, eps0 at those
c magnitudes. Kramer is not asking for those latter values, but may want to
c show them in plots.
      real, dimension(nmmax,nlev):: htot, rbar,eb0, mbar,ebar
      real, dimension(3,nmmax,nlev):: gmt
c ireturn (units: years) is a preset array. It needs to be revised if this set changes
      ireturn=(/4975, 2475, 975, 475, 224, 72/)
      rbarr=0.0
      ratett=0.0
      ebar0=0.0
      ebarr=0.0
      mbarr=0.0
      label='PSA '
c 2008: all deagg distances hypocentral      or r_cd
c      metric='Rrup(Rcd) and Rjb '
      bypass=.true.
c vector initialization
      ratesum=0.0
c array initializations
      rbar=0.
      ebar=0.0
      eb0=0.0
      mbar=0.
      htot=0.0
c Also initialze 3-comp arrays to zero      
      hnew=0.
      mnew=0.
      rnew=0.
c
      nocc=0
      nvac=0
      imax=1
      jmax=1
      jmin=1
      header=.true.
      call getarg(1,msg)
      icit = min(20,index(msg,' ')-1)
      call getarg(2,slong)
      read(slong,'(f9.4)',err=44)xlong
      if(xlong.lt.-114.95)amagmin=4.999
c      if(xlong.lt.-114.95)metric='Rrup and Rjb '
      comm='2008 UPDATE'
      if(xlong.lt.-129.)comm='From OFR 99-36. '
44      call getarg(3,slat)
      read(slat,'(f7.4)',err=46)ylat
      fac=(ylat-38.)/(xlong+121.75)
46     continue
c new 2008: variable vs30 allowed
      call getarg(4,avs30)
      read(avs30,'(f6.1)')vs30
      if(vs30.ge.1500)then
      iclass=1
      elseif(vs30.ge.780.)then
      iclass=2
      elseif(vs30.ge.740.)then
      iclass=3
      elseif(vs30.ge.380.)then
      iclass=4
      elseif(vs30.ge.340.)then
      iclass=5
      elseif(vs30.ge.190.)then
      iclass=6
      else
      iclass=7	!DE low limit available 180 m/s
      endif
       label='PGA '
       ispga=.true.
       do k=1,nlev
      call getarg(4+k,psa)
       read(psa,'(f11.7)')pga(k)
       enddo
       do k=1,nlev      
      call getarg(10+k,rate)
      read(rate,'(e10.5)')arate(k)
      if(arate(k).gt.0.)return(k)=1./arate(k)
      enddo
c return = mean return time (years) of an exceedance. Should be close to
c the ascii arg retrnt      
      call getarg(17,pid)
      ipid=index(pid,' ')-1
      call getarg(18,loc)
       iloc=index(loc,' ')-1
          call getarg(19,logdir)
       iloq=index(logdir,' ')-1
5      format(a)
50      format(a,$)
      read *,nfil
c      print 50,' Uniform weights summing to 1 y/n? '
c      read 5,ch1
      getwt= .true.
c      print 50,' Use mbar,rbar to represent bin center y/n? '
c      read 5,ch1
c      epson=ch1.eq.'z'
c      print *,'Logical variable epson is ',epson
c opaque way to say that data contain epsiln information in last column
      avgctr=epson.or.ch1.eq.'y'.or.ch1.eq.'Y'
      do i=1,nfil
1      continue
c      print 50,' Enter file name : '
      read 5,name1
c      print 5,name1
      ifile=i
c      print 50,' Enter weight for these data: '
      read *,wt,key
c key is a source-category key, added Jan 2003 to combine similar-source hazard.
c in 2008 there are 35 such categories.
      key=min(key,35)       
c      print *,wt
c	write(2,5)logdir(1:iloq)//name1
      open(1,file=logdir(1:iloq)//name1,status='old',err=20)
c eps1 is the avg epsilon for that file (set of sources of a type)
      ratepc=0.0
      read (1,5,end=20)rec
c there are files that are empty: they correspond to 0-hazard cases.
c these files may not even have a header record.
c      print 5,rec
      read(1,5,end=20)rec
c      endif
      wtsum=wtsum+wt
      do j=1,nmax
      read(1,*,err=20,end=20)d,m,h,eps,k      !eps is epsilon0 averged for this bin
      jbin=1+(m-amagmin)/deltam
      jmax=max(jmax,jbin)
      if(jbin.gt.nmmax)stop 'Sum_haz_2008b: a mag exceeds array bound'
      if(k.gt.6)stop 'Sum_haz_2008b: PGA index exceeds 6. 6 is limit'
      tmp=h*wt
      rbar(jbin,k)=rbar(jbin,k)+d*tmp
      mbar(jbin,k)=mbar(jbin,k)+m*tmp
      eb0(jbin,k)=eb0(jbin,k) + eps*tmp
      htot(jbin,k)=htot(jbin,k)+tmp
      ebar0(k)=ebar0(k) + eps*tmp
      avgr(k)=avgr(k)+d*tmp
      avgm(k)=avgm(k)+m*tmp
      ratesum(k)=ratesum(k)+tmp
101      continue
      enddo
20      close(1)
      enddo      !do i=1,nfil
      ratepc=ratesum/arate*100.0
      read 5,name2
      open(2,file=logdir(1:iloq)//name2,status='unknown')
	k=index(name2,' ')-1
	name2=name2(1:k)//'.dat'
      open(1,file=logdir(1:iloq)//name2,status='unknown')
	avgm=avgm/ratesum
      if(slong(1:1).eq.'-')slong=slong(2:10)
c lose the negative and label west. ratett is the mean rate of exceedance
c computed here. It might differ from the rate on command line.
c ireturn is a preset array. It needs to be revised if this set changes
       write(1,704)msg(1:icit),abs(xlong),ylat,avs30,ireturn,
     1 pga,ratesum
704      format('PSHA Deaggregation: Marginal (M,f). site: ',a,
     1' long: ',f7.3,' W., lat: ',f6.3,' N., Vs30(m/s): ',a,
     1/,'NSHMP 2008  See USGS OFR 2008-1128. dM=0.2 below'
     1/,'Return periods: ',6i5,' yrs.',/,'Exceedance PGA =',
     16f8.5,' g.',/,'Computed Rate_Exc ',6(1x,e9.3),', respectively.')
c aratemn is 0.05% of the hazard (1/2-mil). Don't output if less than that.
      do k=1,nlev
	write(1,12)pga(k),ratesum(k),ireturn(k)
12	format(/,'#PGA g = ',f8.5,1x,' computed rate ',e10.5,' expected return ',
     + i5,' years. ',/,'#Relative_frq. M(Mw) R(km) eps0')
      do j=1,jmax
      m=amagmin+(j-0.5)*deltam
	h=htot(j,k)
	if(h.gt.tiny)then
      rbar(j,k)=rbar(j,k)/h
      eb0(j,k)=eb0(j,k)/h
      rmax=max(dist,rmax)
      mbar(j,k)=mbar(j,k)/h
      h=h/ratesum(k)
c convert to percent contribution
c      gmt(3,j,k)=rbar(j,k)
c      gmt(1,j,k)=htot(j,k)
c      gmt(2,j,k)=mbar(j,k)
      write(1,10)h,mbar(j,k),rbar(j,k),eb0(j,k)
      nocc=nocc+1
      else
      nvac=nvac+1
      write(1,10)0.,m,0.,0.
      endif

10      format(1x,f7.3,2x,f4.2,3x,f6.1,2x,f6.3)
      enddo
      enddo
      htot=1.-exp(-htot)	!rate to probability
	write(2,119)jmax,amagmin
119	format(i2,1x,f5.2,' number of rows of Pex values and Mmin, resp. dM=0.2')
	write(2,122)(pga(k),k=nlev,1,-1)
	write(2,121)(avgm(k),k=nlev,1,-1)
121	format(6(1x,e11.5),' avg Mw')	
122	format(6(1x,e11.5),' PGA(g)')
123	format(6(1x,e11.5),' Pex per yr, M: [',f5.2,1x,f5.2,' )')
	m=amagmin
      do j=1,jmax
       write(2,123)(htot(j,k),k=nlev,1,-1),m,m+deltam
       m=m+deltam
       enddo
       close(1)
       close(2)
c      call gmtscript(gmt,rmax,emmin,emmax,pcmax,rbarr,mbarr,rmode,
c     1      mmode,deltar,imax,jmin,
c     1      jmax,icit,iclass)
c      print *,'Avg distance(km) and magnitude(Mw): ',rbarr,mbarr
      end

      subroutine gmtscript(gmt,rmax,emmin,emmax,pcmax,rbarr,mbarr,
     1       rmode,mmode,deltar,imax,jmin,jmax,icit,iclass)
c writes a gmt script. For the most part, run it silently. -V is allowed once,
c when the thing is about finished. Diagnostics are supposed to be going to a log file.
      parameter (deltam=0.2, small=0.05)
      parameter (ndmx=100, nlev=80,nmax=10000,nmmax=50,nsig=6)
c Now Jun 23, 2000: using loc(1:iloc) for locations of all files, newmaps no
c longer used.
      logical epson,maybe
      real mbarr,mmode,m3d
      character*10 psa,rate
      character*14 e3dmsg
      character*4 label,labrs,per,retrnt*5,units,shift,unsh
      character*3 dya,xlab*6
c shift will shift labels depending on location of rmode,mmode
c preferred central location, but they can go left or right
      character*8 sitecl(7)/'NEHRP A','NEHRP B','NEHRP BC','NEHRP C','NEHRP CD',
     +'NEHRP D','NEHRP DE'/
      character*10 red/'255/020/20'/,
     1      orange/'255/134/00'/,yellow/'255/255/00'/,
     1      green/'085/255/85'/,blue/'10/244/145'/,violet/'80/150/255'/
      character*10 color(6),metric*20
      character*12 pid
      character*80 loc,fileout
      character*80 pidpr
      character*24 msg,comm*16
c      character*28 sitepr
      common/vs30/vs30
      common/lab/msg,label,per,retrnt,units,comm
      common/pid/pid,ipid
      common/loca/loc,iloc
      common/lab2/psa,rate,xlong,ylat,metric
      common/lowlim/aratemn,ratett
      common/e3d/r3d,m3d,e3dmsg
c ebar0 is the average of the epsilon0s, i.e., epsilons of SA0 or PGA0
c ebarr is the average of the mean epsilons, where the mean epsilon for
c a given source is the integral from eps0 to infinity of epsilon*normal dens
c divided by uppertail from eps0 to infinity of the normal density
      dimension gmt(nsig+2,ndmx,nmmax),eb0(ndmx,nmmax)
      integer iloc,j1
      color(1)=red
      color(2)=orange
      color(3)=yellow
      color(4)=green
      color(5)=blue
      color(6)=violet
      y=5.0
      dowhile(emmax.gt.y+small)
      y=y+deltam
      enddo
      emmax=y
c M5.0 isn't quite low enough in 2002. try 4.9?
      y=6.9
      dowhile(emmin.lt.y)
      y=y-.5
      enddo
      emmin=y
      if(vs30.lt.660.)then
      labrs='soil'
      else
      labrs='rock'
      endif
      xlab='f25a50'
      ystart=100.
            dy=50.
      if(deltar.lt.10.5)then
      dy=10.0
      ystart=50.
      xlab='f05a10'
      elseif(deltar.lt.20.5)then
      dy=20.
      ystart=60.
      xlab='f10a20'
      endif
      y=ystart

      dowhile(rmax.gt.y)
      y=y+dy
      enddo
      rmax=min(y,990.)
       if(rmax.gt.200.and.rmax.le.300..and.dy.lt.24)then
           xlab='f10a20'
       elseif(rmax.gt.300..and.dy.lt.24.)then
       xlab='f10a50'
c this check is to prevent too many labels for long distance axis
       endif
      y=2.0
      dy=2.
      dya='2.0'
      dowhile(pcmax.gt.y)
      y=y+dy
      if(y.eq.12.)then
      dy=4.
      dya='4.0'
      elseif(y.eq.20.)then
      dy=5.
      dya='5.0'
      elseif(y.eq.40.)then
      dy=10.
      dya='10.'
      endif
      enddo
      pcmax=max(10.0,y)
      xsc=8./rmax
      ysc=4./(emmax-emmin)
      zsc=4./pcmax
      xw=0.19
c the below reductions in xw make thinner columns for the situation where
c rmax is large but deltar is small, i.e., many columns which might be trying
c to occupy the same space.
      if(deltar.lt.20..and.rmax.gt.599.)then
      xw=0.08
      elseif(deltar.lt.20..and.rmax.gt.300.)then
      xw=0.12
      endif
      yw=0.5/(emmax-emmin)
      offset=0.062
c 1996 version Offset = 0.125, apply paint to column face with a little force.
c 2002 version: columns are skinnier
c do the label shift calcs:
      if(rmode.gt.0.31*rmax.and.mmode.gt.0.5*emmax)then
      shift='-2.8'
      unsh='2.80'
      elseif(rmode.gt.0.2*rmax.and.mmode.gt.0.8*emmax)then
      shift='1.60'
      unsh='-1.6'
      elseif(rmode.lt..25*rmax.and.mmode.gt.0.3*emmax)then
      shift='1.00'
      unsh='-1.00'
      elseif(rmode.gt.0.8*rmax.and.mmode.lt.0.5*emmax)then
      shift='-1.0'
      unsh='1.00'
      else
      shift='0.00'
      unsh=shift
      endif
      pidpr=loc(1:iloc)//pid(1:ipid)//per
      ipidr=iloc+ipid+4
      fileout=loc(1:iloc)//'gmt.'//per//pid(1:ipid)
      open(unit=1,file=fileout,status='unknown')
      kd=5
      if(rmax.lt.99.)then
      kd=4
      xlab='f05a10'
      endif
      jpc=3
      if(pcmax.gt.99.)jpc=4
c change "hypocentral" to more generic "Source-site" distance in graph labels      
      write(1,12)rmax,emmin,emmax,pcmax,xsc,ysc,zsc,xlab,
     1 'SOURCE-to-SITE ',dya,xw,yw,pidpr(1:ipidr)
12      format('$1psxyz -R0/',e10.5,'/',f4.2,'/',f4.2,'/0/',e10.5,
     1 '  -Jx',f5.3,'/',f5.3,' -Jz',f5.3,' \\',/,' -B',a6,
     2':"',a15,'DISTANCE (km)":/0.5:"MAGNITUDE (Mw)":/f1a',a3,
     3':"% Contribution to Hazard":4 \\',/,' -So',f4.2,'/',f5.3,
     4' -C$2PL/tone.cpt -L -K +$2PL/wwwdefaults -W2 -E136/27', 
     5' << END > ',a,'.ps')
c The above wwwdefaults is the .gmtdefaults file that should be used for the
c gmt-script duration, by the anon www user. This is the idea anyway.
      do i=1,imax
      gmt(2,i,jmax)=min(gmt(2,i,jmax),0.999*emmax)
c this insures that max mag do not disappear because they are smidge too large.
c this issue may need more work.
c mod May 7,2002: color code columns by average eps0 in the bin. New format 12 and
c four columns in the output file below. jmax is puzzling me.
      if(gmt(3,i,jmax).gt.aratemn)write(1,20)gmt(1,i,jmax),
     1 gmt(2,i,jmax),gmt(3,i,jmax),eb0(i,jmax)
20      format(4f12.3)
      enddo
      write(1,30)
30      format('END')
c run through the epsilons (colors)
      maybe=.true.
      k=4
      j=jmax
      yww=0.001
      off=offset
      dowhile(j.ge.jmin)
      dowhile(maybe.and.k.le.8)
      icount=0
      do i=1,imax
      if(gmt(k,i,j).gt.aratemn)icount=icount+1
      enddo
      if(icount.gt.0)then
      if(yww.gt.0.05)then
      write(1,40)xw,yww,pidpr(1:ipidr)
      else
      write(1,444)xw,yww,pidpr(1:ipidr)
      endif
40      format('$1psxyz -R -Jx -Jz -So',f4.2,'/',f5.3,' -C$2PL/tone.cpt', 
     1      ' -L +$2PL/wwwdefaults -K -O -W2 -E136/27  << END >> ',
     2      a,'.ps')
444      format('$1psxyz -R -Jx -Jz -So',f4.2,'/',f5.3,' -G200 ', 
     1      '-L +$2PL/wwwdefaults -K -O -W2 -E136/27  << END >> ',
     2      a,'.ps')
      do i=1,imax
c the below check prevents pile-up of columns due to faults being near
c bin boundaries. This is likely to happen west of SAF, California.
      if(gmt(k,i,j).gt.aratemn)then
      if(i.eq.1)goto 44
      if(gmt(1,i,j).lt.gmt(1,i-1,j)+deltar*.2)
     1      gmt(1,i,j)=gmt(1,i,j)+0.2*deltar
44      write(1,20)gmt(1,i,j),gmt(2,i,j)-off,
     2      gmt(k,i,j),eb0(i,j)
      endif
      enddo
      write(1,30)
      k=k+1
      else
      maybe=.false.
      endif
      yww=0.001
      off=offset
      enddo
      k=3
      yww=yw
      off=0.0
      j=j-1
      maybe=.true.
      enddo
            j1=index(psa,' ')
            if(j1.eq.0)j1=9
            if(epson.and.per.ne.'0.00')then
      write(1,480)shift,pidpr(1:ipidr),sitecl(iclass),
     + labrs,msg(1:icit),abs(xlong),ylat,
     1 per,psa(1:j1),ratett,retrnt,rbarr,mbarr,ebar0,
     1 rmode,mmode,emodeb0,r3d,m3d,e3dmsg,deltar,unsh,pidpr(1:ipidr)
      elseif(epson)then
       write(1,880)shift,pidpr(1:ipidr),sitecl(iclass),
     + labrs,msg(1:icit),abs(xlong),ylat,
     1 psa(1:j1),ratett,retrnt,rbarr,mbarr,ebar0,
     1 rmode,mmode,emodeb0,r3d,m3d,e3dmsg,deltar,unsh,pidpr(1:ipidr)
      elseif(per.eq.'0.00')then
       write(1,88)shift,pidpr(1:ipidr),sitecl(iclass),
     + labrs,msg(1:icit),abs(xlong),ylat,
     1 ratett,retrnt,psa(1:j1),rbarr,mbarr,
     1 rmode,mmode,deltar,unsh,pidpr(1:ipidr)
      else
      write(1,48)shift,pidpr(1:ipidr),sitecl(iclass),
     1 labrs,msg(1:icit),abs(xlong),ylat,
     1 psa(1:j1),per,ratett,retrnt,rbarr,mbarr,
     1 rmode,mmode,deltar,unsh,pidpr(1:ipidr)
      endif
48      format('$1pstext -X',a4,' -Jx1/1  +$2PL/wwwdefaults ',
     1'-R0/11/0/8.5 -K -O << END >> ',
     1a,'.ps',/,
     1 '3.0 6.45 18 0. 4 1 PSH Deaggregation on ',a8,1x,a4,/,
     1 '3.0 6.2 18 0. 4 1 ',a,1x,f7.3,'@+o@+ W. ',f6.3,' N.',
     1 /,'3.0 6.0 14 0. 4 1 Spectral Accel. >=',a,
     1 'g for SA period ',a4,' sec',/,
     1'3.0 5.8 14 0 4 1 Ann. Exceedance Rate ',e7.3,
     1'. Mean Return Time ',a4,' years',/,
     2 /,'3.0 5.6 14 0. 4 1 Mean R ',f5.1,' km, mean M ',f4.2,/,
     3 '3.0 5.4 14 0. 4 1 Modal R ',f5.1,' km, modal M ',f4.2,
     3/,'3.0 5.2 14 0. 4 1 Binning: DeltaR ',f3.0,' km, deltaM 0.2'
     4/,'END',/,
     5'$1pstext +$2PL/wwwdefaults ',
     5'$2PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O >> ',a,'.ps')
88      format('$1pstext -X',a4,' -Jx1/1  +$2PL/wwwdefaults ',
     1'-R0/11/0/8.5 -K -O << END >> ',
     1a,'.ps',/,
     1 '3.0 6.45 18 0. 4 1 PSH Deaggregation on ',a8,1x,a4,/,
     1 '3.0 6.2 18 0. 4 1 ',a,1x,f7.3,'@+o@+ W. ',f6.3,' N.',/,
     1'3.0 6.0 14 0 4 1 Ann. Exceedance Rate ',e8.2,
     1'. Mean Return Time ',a4,' years',/,
     1 '3.0 5.8 14 0. 4 1 Peak Horiz. Ground Accel.>=',a,' g'
     2 ,/,'3.0 5.6 14 0. 4 1 Mean R ',f5.1,' km, mean M ',f4.2,/,
     3 '3.0 5.4 14 0. 4 1 Modal R ',f5.1,' km, modal M ',f4.2,
     3/,'3.0 5.2 14 0. 4 1 Binning: DeltaR ',f3.0,' km, deltaM 0.2'
     4/,'END',/,
     5'$1pstext +$2PL/wwwdefaults ',
     5'$2PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O >> ',a,'.ps')
480      format
     1('$1pstext -X',a4,' -Jx1/1 -R0/11/0/8.5 -K +$2PL/wwwdefaults ',
     1' -O << END >> ',
     1a,'.ps',/,
     1 '3.0 6.45 18 0. 4 1 PSH Deaggregation on ',a8,1x,a4,/,
     1 '3.0 6.2 18 0. 4 1 ',a,1x,f7.3,'@+o@+ W, ',f6.3,' N.',
     1 /,'3.0 6.0 14 0. 4 1 SA period ',a4,' sec. Accel.>=',a,' g'
     1/,'3.0 5.8 14 0 4 1 Ann. Exceedance Rate ',e8.3,
     1'. Mean Return Time ',a5,' yrs',/,
     2 '3.0 5.6 14 0. 4 1 Mean (R,M,@~e@~@-0@-) '
     2,f5.1,' km,',f4.2,
     2', ',f5.2,/,
     3 '3.0 5.4 14 0. 4 1 Modal (R,M,@~e@~@-0@-) ='
     3,f5.1,' km,',f5.2,
     3',',f5.2,' (from peak R,M bin)',/,
     3 '3.0 5.2 14 0. 4 1 Modal (R,M,@~e@~*) =',f5.1,' km,',f5.2,
     3',',a14,' (from peak R,M,@~e@~ bin)',
     3/,'3.0 5.0 14 0. 4 1 Binning: DeltaR=',f3.0,' km, deltaM=0.2'
     3,', Delta@~e@~=1.0',
     4/,'END',/,
     5'$1pstext $2PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O +$2PL/wwwdefaults >> ',a,'.ps')
880      format
     1('$1pstext -X',a4,' -Jx1/1 -R0/11/0/8.5 -K +$2PL/wwwdefaults ',
     1' -O << END >> ',a,'.ps',/,
     1 '3.0 6.45 18 0. 4 1 PSH Deaggregation on ',a8,1x,a4,/,
     1 '3.0 6.2 18 0. 4 1 ',a,1x,f7.3,'@+o@+ W, ',f6.3,' N.',
     1 /,'3.0 6.0 14 0. 4 1 Peak Horiz. Ground Accel.>=',a,' g'
     1/,'3.0 5.8 14 0 4 1 Ann. Exceedance Rate ',e8.3,
     1'. Mean Return Time ',a5,' years',/,
     2 '3.0 5.6 14 0. 4 1 Mean (R,M,@~e@~@-0@-) ',f5.1,' km, ',f4.2,
     2', ',f5.2,/,
     3 '3.0 5.4 14 0. 4 1 Modal (R,M,@~e@~@-0@-) = ',f5.1,' km, ',f4.2,
     3', ',f5.2,' (from peak R,M bin)',/,
     3 '3.0 5.2 14 0. 4 1 Modal (R,M,@~e@~*) =',f5.1,' km,',f5.2,
     3',',a14,' (from peak R,M,@~e@~ bin)'
     3/,'3.0 5.0 14 0. 4 1 Binning: DeltaR ',f3.0,' km, deltaM=0.2'
     3,', Delta@~e@~=1.0',
     4/,'END',/,
     5'$1pstext $2PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O +$2PL/wwwdefaults >> ',a,'.ps')
      write(1,50)labrs,vs30,comm,pidpr(1:ipidr)
50      format('$1psxy ','$2PL/squares_e ',
     2'-U"Distance (R), magnitude (M), epsilon ',
     2'(E0,E) deaggregation for a site on ',a,' with ',
     2'average vs=',f5.0,' m/s top 30 m.', 
     2' USGS CGHT PSHA',a14,' Bins with lt 0.05% contrib. omitted"',
     2' \\',/,
     2' -Jx1/1 -R -O -Ss.25 -C$2PL/tone.cpt +$2PL/wwwdefaults >> ',a,'.ps')
      close(1)
      return
      end      
      

