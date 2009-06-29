c sum_haz_2009.f used on the WWW computer, gldplone . 
c May 2009. multiple GMPE contributions may be reported. Make arg6 1 to get report
c   for each atten model as well as mean hazard report.
c jan 7 2009. minor revs to initialization
c --- this pgm is set up for NSHMP 2008 update; it may need work in these areas:
c --- (1) seismograms for modal event. or mean event based on RV point src.
c --- (2)  seismogram input files are for rock not soil. 
c --- A better solution might be a variety of plausible scenarios with finite source
c ---
c Primary purpose: sums w(i)*haz(i,j), i=1,nfil, j=1,...,nsig
c   in (distance, magnitude, epsilon) [ (R,M,e) ] bins. nfil=number of files
c   of PSHA data, can be up to 20 (or modify parameter nfil). haz() is mean
c   annual rate of ground motion exceedances. w(i) is pre-assigned weight.
c   Nsig is number of epsilon bins, where sig represents aleatory and possibly
c   epistemic GM variability. (epistemic if more than one attenuation model used)
c Second purpose: generates 3 gmt scripts for graphic PSHA display The first is
c   McGuire's (R,M,e) bins, 2nd is geographic deaggregation of hazard (not ready),
c  third script prepares 'grams.
c Third purpose: determines hazard associated with specific faults, fault zones,
c  subduction zones. This is in demand: not just the hazard as a function of distance
c  and magnitude, but also hazard as a function of active structures that people
c  recognize. Where fault density is high, user wants to know which F dominates
c  the hazard. Differentiating hazard by fault is frequently a problem in 
c Los Angeles & other W California sites. 
c
c May 8, 2001: modified to plot 6 seismograms for the modal (M,R) and output
c    an ascii file of these 6 seismograms, from Boore's smsim_td, version 2.10. 
c Apr 6, 2001: 
c    modified to report modal R* M* e* in both ascii and graphic output files
c
c No makefile. use gnu fortran
c PC Linux compile (In 2008, this is gldplone): 
c gfortran ./sum_haz_2009.f -static -O datetime.o  iosubs.o  recipes.o rvtdsubs.o  \
c    smsim_tdww2.o -o ./sum_haz_2009 -ffixed-line-length-none
c Solaris compile: 
c f77 sum_haz_2009.f  -e -fast datetime.o  iosubs.o  recipes.o rvtdsubs.o smsim_tdww2.o -o sum_haz_2009
c  note: -Bstatic with f77 or f95 compilation goes to pieces in a hurry, Solaris System objects to -Bstatic. the
c /usr/lib/libmvec.so is shared. 
c   
c
c Makefile links up several of Boore's subroutine packages & iosubs.o. Use -fast.
c The -static argument invokes static library linking, needed for Web computers.
c
c Usage: sum_wt_haz argj  j=1,...,13
c
c arg1 = site name (no blanks!)
c arg2 = site longitude (E+, dec degrees)
c arg3 = site latitiude (dec degrees)
c arg4 = target Return period (a4) units: years or kiloyears
c arg5 = vs30 m/s
c arg6 = iattn. iattn=1 to get individual gmpe reports new may 2009.
c arg7 = PSA (PGA) period. PGA convention is 0.00 sec.
c arg8 = PSA (PGA) level, units: g
c arg9 = computed annual rate of exceedances (exponential notation please)
c arg10 = pid (process id, for names)
c arg11 = absolute path name of data files  Out
c arg12 = location of input files (some need to distinguish these)
c arg13 = TSA is an estimate of approx. UHS long-period corner.
c arg14 = SA5, the 5-hz spectral accel level.
c arg15 =indicator, 1=> generate and plot 6 seismograms for modal (M,R).
c            0 => don't generate or plot seismograms. (Boore's smsim
c            initially used to generate the accelerograms.)
c Below, read instructions assume that input files have 1 header record 
c Epsilon bin boundaries are assumed at -2, -1, 0, 1, and 2 sigma. The hazard data are
c  presumed to be arranged "all data"; GM<mu+2 sigma;
c  GM < mu + 1 sigma; GM < mu; GM<mu-1sigma, and
c  GM<mu-2 sigma.  Thus, the binned hazard data overlap. (This was for graphics.)
c Mean rates are averaged from srcs in nfil files. Each file has an input weight w(i).
c Input records have a distance, R, and magnitude, M, preceeding the hazard vector.
c these R and M values may be averages within binned distance and magnitude
c intervals, ie, R=Rbar and m=mbar based on hazards determined in those intervals.
c In this program, Mbar,Rbar values are computed and output for the 
c averaged or weighted sum. 
c CEUS bin sizes are usually 25 km for delta-R, and 0.2 units
c for delta-M or dM. WUS delta-R is 10 km or 25 km. 
c dM may be changed by modifying parameter  deltam. dR is changed for different
c geographic regions: higher hazard generally has lower dR, for example, NMSZ.
c 
c When OKed, sum_haz_2009 is for automated web-site PSHA deaggregations.
c Input Magnitudes assumed to be greater than 4.5 (Mw). For lower magnitudes,
c some code modifications are required. The 4.5 limit corresponds to a 5.0(MbLg)
c limit in the 2002 National Seismic Hazard Maps. .
c Two output files: name specified and gmt.period (period is 4 char SA period)
c Steve Harmsen updated 6/2001 to 8/2008
c
      parameter ( deltam=0.2, amagmin=4.5,nsrcmx=35,nattmx=30)
      parameter (ndmx=100, nfmx=80,nmax=10000, nsig=6,nmmax=50)
c nmax=max number of records per file. Be careful here!
      logical avgctr,epson,getwt,header,ispga,aok,bypass
      logical fltok,ice,countme
      real ratecat(40,0:30)
      real, dimension(80,0:30):: rlow
      real d,eps,m,mag,m3d,hmax/0.0/
      real rmax/10./,deltar/25./,h_int,h_intm/0.0/
      character*32 name1,name2,a32
      character*1 ch1
      character*32 srcdesc(35)
     +/'Cascadia M8.0-M8.2 Floating','Cascadia M8.3-M8.7 Floating','Cascadia Megathrust',
     +'Wash-Oreg Cascades-West faults +','Wash-Oreg Extensional faults','Basin & Range M<6.5 flts',
     +'Basin & Range M>6.5 Char','Basin & Range M>6.5 GR','NV and Utah M<6.5 flts',
     +'Nevada M>6.5 Char','Nevada M>6.5 GR','Utah M>6.5 Char no Wasatch','UT M>6.5 GR no Wasatch',
     +'Wasatch Segmented Char','Wasatch Segmented GR','Wasatch Unsegmented M7.4',
     +'WUS Compr crustal gridded','CA-NV Shear Zones','California B-faults Char',
     +'California B-faults GR','California A-faults',
     +'CA Compr. crustal gridded','50-km Deep Intraplate','Puget Lowlands gridded',
     +'SAF Creeping Sec:gridded','Brawley Zone gridded','Mendocino Zone gridded',
     +'Mojave Zone gridded', 'San Gorgonio Zone gridded','Extensional Zones gridded', 
     +'Cheraw CO or Meers OK faults','New Madrid SZ no clustering','CEUS gridded',
     +'Charleston SC M<7.2; 2 zones','Charleston SC M>7.2; 2 zones'/
       character*52 amsg(2)/' (some WUS atten. models use Site Class not Vs30).',
     +' CEUS atten. model site cl BC(firm) or A(hard).'/
      character*24 attnmodel(0:30)
     +/'Mean Hazard w/all GMPEs','Spudich et al. Ext','Toro et al. 1997',
     +'Sadigh et al. 1997','Atkinson-Boore06,140 bar','Attn model 5',
     +'Frankel et al., 1996','Somerville Rifted FinFlt','Attn model 8',
     +'Attn model 9','Campbell CEUS Hybrid','Attn model 11',
     +'Atkinson-Boore03 Casc.','Youngs et al. 1997','Atten Model 14',
     +'Silva 1-corner','Atten Model 16','Atten Model 17','Atkinson-Boore03 Glbl.',
     +'Tavakoli and Pezeshk 05','Atkinson-Boore06,200 bar',   
     +'Boore-Atkinson 2008','Campbell-Bozorgnia 2008',
     +'Chiou-Youngs 2008','Abrahamson-Silva NGA','Atten Model 25',
     +'Crouse 1991','Zhao et al. 2006','Kanno et al. 2006','Gregor et al. 2006',
     +'Atten Model 30'/
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
      common/epsilon/epson,ebar0,emodeb0,eb0
      common/e3d/r3d,m3d,e3dmsg
       character*4 epr(10)
     1/'5.00','4.00','3.00','2.00','1.00','0.50','0.30','0.20','0.10','0.00'/  
      character*80 loc,rec,logdir
      character*32 fltnm(40)
c fltnm will hold short descriptive message (Char or GR; dip of fault)
c about significant faults.
       character*1 cseis,flg*2,flg3*3
      logical maybe,crustal(0:30)
	data crustal/12*.true.,3*.false.,3*.true.,.false.,
     + 7*.true.,5*.false./
      real, dimension(nfmx):: azi
      real, dimension (nfmx,0:30):: amag,rsav,esav
      integer k00/1/,k01/1/,k02/1/
      integer iloc,ireg,kf,kfmin,kfmax,kflim/1/,nocc,nvac
      real, dimension(0:200,0:30) :: flthaz,fm,fd,feps
c flthaz(0) hazard for non-faults, otherwise for a fault having an index
       real tmp
       real, dimension(0:30) :: avgr, avgm,eps1,pr50,rbarr,mbarr,ebar0,emodeb0,
     +  ratelow,ratett,ratesum,rmode,mmode,pcmax,gmpewt
c the array gmt needs to be upped to 4D to be able to plot individual GMPE
c data with the GMT subroutine. This was not done in initial code preparation.
c SHarmsen May 2009.
      real h(nsig),gmt(nsig+2,ndmx,nmmax)
c awesome 4d arrays:
      real, dimension(nsig,ndmx,nmmax,0:30) :: hnew,mnew,rnew
c  mnew, rnew will store avg values in each 3D bin      
      real hsum(nsig,0:30),emmin/10./,emmax/1./
      real, dimension(80,0:30):: rbf,embf,ebf,ratef
      real, dimension(ndmx,nmmax,0:30):: rbar,eb0, mbar,ebar
      data gmpewt/1.0,.1,.2,.1,.1,.1,.1,.2,.1,
     +.1,.1,.1,.25,.5,.1,.1,.1,.1,.25,
     +0.1,.1,.3333,.3333,.3333,.1,.1,
     +0.1,0.5,.1,.1,.1/
c vector initialization
      rbarr=0.0
      ratett=0.0
      ebar0=0.0
      mbarr=0.0
      mmode=0.0
      rmode=0.0
      pcmax=0.0
      label='PSA '
c 2008: all deagg distances hypocentral      or r_cd
      metric='Rrup(Rcd) and Rjb '
      bypass=.true.
      embf=0.0
      rbf=0.0
      ebf=0.0
      ratef=0.0
      lstft=0
c array initialiazation
      rbar=0.
      ebar=0.0
      eb0=0.0
      mbar=0.
      ratecat=0.0
c Also initialze 4D arrays to zero      
      hnew=0.
      mnew=0.
      rnew=0.
c
      hsum=0
      nocc=0
      nvac=0
      imax=1
      jmax=1
      jmin=1
      header=.true.
      if(iargc().ge.9)then
      call getarg(1,msg)
      icit = min(20,index(msg,' ')-1)
      call getarg(2,slong)
      read(slong,'(f9.4)',err=44)xlong
c      if(xlong.lt.-114.95)metric='Rrup and Rjb '
      comm='2008 UPDATE'
      if(xlong.lt.-129.)comm='From OFR 99-36. '
44      call getarg(3,slat)
      read(slat,'(f7.4)',err=46)ylat
      fac=(ylat-38.)/(xlong+121.75)
      if(ylat.lt.12.0)then
      deltar=10.
       comm='Panama PSHA'
      elseif(ylat.lt.38..and. fac .lt. -4./5.75)then
      deltar=10.
       elseif(abs(ylat-36.7).lt.0.3.and.abs(xlong+116.4).lt.0.3)then
       deltar=10.
c ca-nv border region, lots of faults
      elseif(abs(ylat-38.).lt.2..and.abs(xlong+118.).lt.2.)then
       deltar=10.
c Yucca Mountain area (potential national nuclear waste repository site)
      elseif(xlong.le.-119..or.xlong.gt.0.)then
      deltar=10.0
      elseif(xlong.lt.-114.01.and.ylat.lt.42.0)then
c imperial valley region to Yuma; eastern Nv
      deltar=10.0
      elseif(abs(xlong+118.1).lt.0.4.and.abs(ylat-36.5).lt.0.4)then
c Owens Valley area
       deltar=10.0
      elseif(abs(xlong+116.6).lt.1.1.and.abs(ylat-34.6).lt.0.5)then
c landers area
       deltar=10.0
c West COast of USA y
      elseif(abs(ylat-40.7).lt.0.4.and.abs(xlong+111.9).lt.0.3)then
      deltar=10.0
c Salt Lake City and Wasatch front area
       elseif(abs(ylat-36.).lt.1.5.and.abs(xlong+106.0).lt.1.35)then
c interesting place in new mexico
      deltar=10.0
       elseif(abs(ylat-43.5).lt.4.5.and.abs(xlong+111.0).lt.1.35)then
c Yellowstone- Snake River plain region
      deltar=10.0
      
      elseif(abs(ylat-36.0).lt.0.9.and.abs(xlong+90.).lt.0.95)then
c NMSZ
      deltar=10.0
       elseif(abs(ylat-32.8).lt.0.3.and.abs(xlong+80.).lt.0.3)then
       deltar=10.0
c Charleston, SC, vicinity - - interesting effects from "areal source" dilution.
      endif
46      call getarg(4,retrnt)
c new 2008: variable vs30 allowed
      call getarg(5,avs30)
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
      call getarg(6,ch1)
      read(ch1,'(i1)')iattn	!new item 6 in 2009
      call getarg(7,per)
       if(per.eq.'0.00')then
       label='PGA '
       ispga=.true.
       labs='PSHA'
       ipr=10	!pga is 10th per (long to short here)
      per1=0.0
       else
       ispga=.false.
       labs=per(1:3)//'s'
      read(per,'(f4.2)')per1
       ipr=1
       dowhile (per.ne.epr(ipr).and.ipr.lt.10)
       ipr=ipr+1
       enddo
c Might use ipr to gather attenuation coefficients (not in use 2008-2009)
       endif
       aok=per.eq.epr(ipr)
c aok "true" means that program knows what period was input, otherwise, we
c  can't do further attenuation-model epsilon analysis.
      call getarg(8,psa)
       read(psa,'(f11.7)')sa1
       if(sa1.gt.0.15)then
       write(psa(1:6),'(f6.4)')sa1
       psa(7:10)='    '
       elseif(sa1.gt.0.015)then
       write(psa(1:7),'(f7.5)')sa1
       psa(8:10)='    '
       endif
c psa is a char10. We need to get a real value units g. This step can kill!
      call getarg(9,rate)
      read(rate,'(e10.5)')arate
      waterlev=0.01*arate
      if(arate.gt.0.)return=1./arate
c return = mean return time (years) of an exceedance. Should be close to
c the ascii arg retrnt      
      call getarg(10,pid)
      i=1
      iaxmx=0           
      ipid=index(pid,' ')-1
      call getarg(11,loc)
       iloc=index(loc,' ')-1
          call getarg(12,logdir)
       iloq=index(logdir,' ')-1
       call getarg(13,TSA)
       call getarg(14,SA5)
       read(SA5,'(f8.5)')sa5w
       sa5w=sa5w*980.5
c sa5w and sa2 are used to scale grams; sa2 from sa1
       call getarg(15,cseis)
      else
      print *,'Usage sum_haz_2009 site long lat rate vs30 iattn per gmlev rate2'
      endif
c      print 50,' How many files to average : '
5      format(a)
50      format(a,$)
      read *,nfil
      if(nfil.gt.nfmx)stop ' sum_haz_2009: 80 files max. Bump & rerun'
c      print 50,' Uniform weights summing to 1 y/n? '
      read 5,ch1
      getwt= .true.
c      print 50,' Use mbar,rbar to represent bin center y/n? '
      read 5,ch1
      epson=.true.
      avgctr=epson
      ratelow=0.0
      do i=1,nfil
1      continue
c      print 50,' Enter file name : '
      read 5,name1
c      print 5,name1
      ifile=i
      rlow(ifile,0:30)=0.0
       idot=index(name1,' ')-1
c      print 50,' Enter weight for these data: '
      read *,wt,key
c key is a source-category key, added Jan 2003 to combine similar-source hazard.
c in 2008 there are 35 such categories.
      key=min(key,35)       
c      print *,wt
      ratesum=0.0
      open(1,file=logdir(1:iloq)//name1,status='old',err=20)
c kfmin,kfmax fault indexes. We will save summary info for each fault
      kfmin=1000
      kfmax=0
      eps1=0.0
      avgr=0.0
      avgm=0.0
      flthaz=0.0
      fm=0.0
      fd=0.0
      feps=0.0
c eps1(iax) is the avg epsilon for that GMPE or for mean hazard (iax=0)
c      ratepc=0.0
      if(header)then
      read (1,5,end=20)rec
c there are files that are empty: they correspond to 0-hazard cases.
c these files may not even have a header record.
c      print 5,rec
      read(1,5,end=20)rec
      endif
      wtsum=wtsum+wt
      indb=index(name1,per)+3
      name2=name1(1:idot)//'.FDEAG'
      inquire(file=logdir(1:iloq)//name2,exist=fltok)
      if(fltok)then
      open(3,file=logdir(1:iloq)//name2,status='old')
            read(3,5,end=16)rec
            read(3,5,end=16)rec
            read(3,5,end=16)rec
c information indexed to specific faults or f. segments. 
c Add back-azimuth field (averaged over rupture models) Mar 2009.
	countme=.false.
      do ift=1,600
      read(3,27,end=16)rbft,embft,ebft,pr,backaz,a32,iax
27      format(1x,f9.3,1x,f7.3,1x,f5.2,1x,e11.5,1x,f8.2,1x,a32,1x,i2)
      if(iax.eq.0..and.wt*pr.gt.0.02*arate)then
      write(anum,'(i4)')ift
      lstft=lstft+1
      rbf(lstft,iax)=rbft
      embf(lstft,iax)=embft
      ebf(lstft,iax)=ebft
      ratef(lstft,iax)=wt*pr/arate*100.
      fltnm(lstft)=a32
      azi(lstft)=backaz
      countme=.true.
      elseif(iax.eq.0)then
      countme=.false.
      elseif(iax.gt.0.and.countme)then
c iax=0 is used to increment lstft. azi() and fltnm() do not change with GMPE
      rbf(lstft,iax)=rbft
      embf(lstft,iax)=embft
      ebf(lstft,iax)=ebft
      ratef(lstft,iax)=wt*pr/arate*100.
      endif
      enddo
16      close(3)
      endif  !fltok 
      do j=1,nmax
       read(1,*,err=20,end=20)d,m,h,kf,eps,iax      !eps is epsilon0 averged for this bin
	if(iattn.eq.0.and.iax.gt.0)goto 101	!skip these accumulations for 
c individual gmpes if not requested.
      kfmin=min(kf,kfmin)
      kfmax=max(kf,kfmax)
      iaxmx=max(iax,iaxmx)
c removed the clause about working sans eps. Eps just has to be there.
      ibin=1+d/deltar
      jbin=1+(m-amagmin)/deltam
      if(jbin.gt.nmmax)stop 'Sum_haz_2009: a mag exceeds array bound'
      tmp=h(1)*wt
      rbar(ibin,jbin,iax)=rbar(ibin,jbin,iax)+d*tmp
      mbar(ibin,jbin,iax)=mbar(ibin,jbin,iax)+m*tmp
      eb0(ibin,jbin,iax)=eb0(ibin,jbin,iax) + eps*tmp
      ebar0(iax)=ebar0(iax) + eps*tmp
      if(eps.le.0.)then
c ratelow(iax) stores rate of events with mu > SA0      
      ratelow(iax)=ratelow(iax)+tmp
      rlow(ifile,iax)=rlow(ifile,iax)+h(1)
      endif
      avgr(iax)=avgr(iax)+d*tmp
      avgm(iax)=avgm(iax)+m*tmp
      eps1(iax)=eps1(iax)+eps*tmp
      flthaz(kf,iax)=flthaz(kf,iax)+tmp
      fm(kf,iax)=fm(kf,iax)+m*tmp
      fd(kf,iax)=fd(kf,iax)+d*tmp
      feps(kf,iax)=feps(kf,iax)+eps*tmp
c feps() is the epsilon for that site, fault, for that SA(period) or PGA.
      k=1
      dowhile(h(k).gt.0..and.k.lt.6)
      h_int=(h(k)-h(k+1))*wt
      mnew(k,ibin,jbin,iax)=mnew(k,ibin,jbin,iax) + m*h_int
      rnew(k,ibin,jbin,iax)=rnew(k,ibin,jbin,iax) + d*h_int
      k=k+1

      enddo
      rbarr(iax)=rbarr(iax)+d*tmp
      mbarr(iax)=mbarr(iax)+m*tmp
      ratesum(iax)=ratesum(iax)+tmp
      do k=1,nsig
c h(1) is the rate integrated from SA0 to mu+emax*sigma (inf in 1996)      
      hnew(k,ibin,jbin,iax)=hnew(k,ibin,jbin,iax)+wt*h(k)
      hsum(k,iax)=hsum(k,iax)+wt*h(k)
      enddo
      jmax=max(jmax,jbin)
      imax=max(imax,ibin)
      jmin=min(jbin,jmin)
101      continue
      enddo
20      close(1)
	do iax=0,iaxmx
        ratecat(key,iax)=ratecat(key,iax)+ratesum(iax)
        rsav(key,iax)=rsav(key,iax)+avgr(iax)
        amag(key,iax)=amag(key,iax)+avgm(iax)
        esav(key,iax)=esav(key,iax)+eps1(iax)
        ratett(iax)=ratett(iax)+ratesum(iax)
	if(ratesum(iax).gt.0)
     + print *, name1(1:12),' annual rate: ',ratesum(iax), ' cumu ',ratett(iax),iax
        enddo	!iax
c      ratepc=ratesum(0)/arate*100.0
      enddo      !do i=1,nfil
	iax=0	
      read 5,name2
      i=3
      dowhile(name2(i:i).ne.'.')
      i=i+1
      enddo
      name1=name2(1:20)
      open(1,file=logdir(1:iloq)//name1,status='unknown')
      if(slong(1:1).eq.'-')slong=slong(2:10)
c      pr50(iax)=1.-exp(-50.*ratelow(iax))
      if(ratett(iax).gt.0.)return=1./ratett(iax)
c pr50(iax) is the probability of >=1 event with median >= SA0 in a 50-year interval
c pr50(iax) uses different "medians" from the attenuation models, weighed. Also can
c compute a pr50(iax) for each attenuation model (this part should be commented
c out for www runs)
c      4,file=pid(1:ipid)//'.50yr.prob.'//per,status='unknown')
c      write(4,804)msg,slong(1:8),ylat,
c     1 psa(1:9),ratett(iax)
c      do ifile=1,nfil
c      if(rlow(ifile,iax).gt.0.0001)then
c      pr50i(iax)=1.-exp(-50.*rlow(ifile,iax))
c      write(4,904)pr50(iax)i,msav(ifile)
c904    format('#Pr[at least one eq with median motion>=SA0 in 50 yrs]='
c     1,f7.5,' using srcs described as ',a)
       call get_date(dateis)
c       close(4)
      if(xlong.gt.-104..and.ylat.gt.30.)then
      ireg=2
      elseif(xlong.gt.-101.)then
      ireg=2
      else
      ireg=1
      endif
c Remove the negative sign and label west. ratett is the mean rate of exceedance
c computed here. It might differ from the rate on command line.
	do iax=0,iaxmx
	if(ratett(iax).gt.waterlev)then
c write table for each atten model with 1% contribution or more. This will
c be every GMPE except for those with a really token appearance (none such in
c USGS scheme. GMPE minweight is about 0.1).
c pr50 is a conditional probability, based on this GMPE being "true". Or
c that the mean is "true". This explains the 1/gmpewt() factor below.
      pr50(iax)=(1.-exp(-50.*ratelow(iax)))/gmpewt(iax)
       if(ispga)then
       write(1,704)msg(1:icit),abs(xlong),ylat,vs30,amsg(ireg),retrnt,
     1 psa(1:9),ratett(iax),pr50(iax),attnmodel(iax)
       else
      write(1,4)msg(1:icit),abs(xlong),ylat,vs30,amsg(ireg),dateis,retrnt,
     1per,label,psa(1:9),ratett(iax),pr50(iax),attnmodel(iax)
           endif
4      format('PSHA Deaggregation. % Contributions. Site: ',a,
     1' long: ',f7.3,' d W., lat: ',f6.3,' N.',/,'Input Vs30 (m/s) =',f6.1,a,/,
     1'NSHMP 2007-08 update. See USGS OFR 2008-1128. Analysis on DaMoYr:',
     1a10,/,'Mean Return Period: ',a5,' years. ',a4,' s. ',a4,'=',
     1a9,'g. Weight * Computed_Rate_Ex ',e9.3,/,
     1'#Pr[at least one eq with median motion>=PSA in 50 yrs]=',f7.5,/,
     1'#This deaggregation corresponds to ',a,/,
     1'DIST(km) MAG(Mw) ALL_EPS EPSILON>2  1<EPS<2 0<EPS<1 -1<EPS<0 -2<E',
     2'PS<-1 EPS<-2')
704      format('PSHA Deaggregation. %contributions. site: ',a,
     1' long: ',f7.3,' W., lat: ',f6.3,' N.',/,' Vs30(m/s)=',f6.1,a,
     1/,'NSHMP 2007-08  See USGS OFR 2008-1128. dM=0.2 below'
     1/,'Return period: ',a5,' yrs. Exceedance PGA =',
     1a9,'g. Weight * Computed_Rate_Ex ',e9.3,/,
     1'#Pr[at least one eq with median motion>=PGA in 50 yrs]=',f7.5,/,
     1'#This deaggregation corresponds to ',a,/,
     1'DIST(KM) MAG(MW) ALL_EPS EPSILON>2  1<EPS<2 0<EPS<1 -1<EPS<0 -2<E',
     2'PS<-1 EPS<-2')
      aratemn=ratett(iax)*0.00050
c aratemn is 0.05% of the mean or GMPE hazard (1/2-mil). Don't output if less than that.
      do j=1,jmax
      m=amagmin+(j-0.5)*deltam
      m=min(8.,m)
       do i=1,imax
      d=deltar*i-deltar*0.5
      if(hnew(1,i,j,iax).gt.aratemn)then
      dist=rbar(i,j,iax)/hnew(1,i,j,iax)
      epsilon=eb0(i,j,iax)/hnew(1,i,j,iax)
      eb0(i,j,iax)=epsilon
      rmax=max(dist,rmax)
      mag=mbar(i,j,iax)/hnew(1,i,j,iax)
      emmax=max(emmax,mag)
      emmin=min(emmin,mag)
c convert to percent contribution. The percent is with respect to mean hazard. 
c Thus, if Moe, Joe and Larry each yields 8% for a given R,M, the mean hazard for that R,M
c should be 24%.
      do k=1,nsig
      hnew(k,i,j,iax)=100.*hnew(k,i,j,iax)/ratett(0)
      if(iax.eq.0)gmt(k+2,i,j)=hnew(k,i,j,iax)
      enddo
      if(iax.eq.0)then
      gmt(1,i,j)=dist
      gmt(2,i,j)=mag
      endif
      all=hnew(1,i,j,iax)
      if(pcmax(iax).lt.all)then
c mode candidate
      rmode(iax)=dist
      mmode(iax)=mag
      emodeb0(iax)=epsilon
      pcmax(iax)=all
      endif
      do k=1,nsig-1
      hnew(k,i,j,iax)=hnew(k,i,j,iax)-hnew(k+1,i,j,iax)
      if(hnew(k,i,j,iax).gt.h_intm)then
      h_intm=hnew(k,i,j,iax)
      kbc=k
      ibc=i
      jbc=j
c Bazzurro-Cornell mode candidate in three dimensions.
      endif
      enddo
c above changed cumulative to interval.
c Note on (dist,mag):
c Although computed, these values are only optionally output
      if(avgctr)then
      d=dist
      m=mag
      endif
      write(1,10)d,m,all,(hnew(k,i,j,iax),k=1,nsig)
      nocc=nocc+1
      else
      nvac=nvac+1
      endif

10      format(1x,f6.1,4x,f4.2,7(3x,f6.3))
      enddo
      enddo
      if(ratett(iax).gt.1e-12)then
      rbarr(iax)=rbarr(iax)/ratett(iax)
      ebar0(iax)=ebar0(iax)/ratett(iax)
      mbarr(iax)=mbarr(iax)/ratett(iax)
      endif
      write(1,710)labs,label,ratett(iax)/ratett(0)*100.,rbarr(iax),mbarr(iax),ebar0(iax),
     1 rmode(iax),mmode(iax),emodeb0(iax)
           if(h_intm.gt.0.)then
           h_intm= ratett(0)*0.01*h_intm
      if(kbc.eq.1)then
      e3dmsg='> 2 sigma'
      elseif(kbc.eq.2)then
      e3dmsg=' 1 to 2 sigma'
      elseif(kbc.eq.3)then
      e3dmsg=' 0 to 1 sigma'
      elseif(kbc.eq.4)then
      e3dmsg=' -1 to 0 sigma'
      else
      e3dmsg='-2 to -1 sigma'
      endif	!kbc business
      r3d=rnew(kbc,ibc,jbc,iax)/h_intm
      m3d=mnew(kbc,ibc,jbc,iax)/h_intm
c write the percent contribution to total hazard, not conditional based
c on a specific GMPE (or iax). Total hazard has iax=0.
           write(1,910)r3d,m3d,e3dmsg,h_intm*100./ratett(0)
910      format(' MODE R*=',f6.1, 'km; M*=',f5.2,
     1 '; EPS.INTERVAL:',a14,' % CONTRIB.= ',f6.3)
       endif
       if(mmode(iax).gt.8.)write(1,715)
       if(iax.eq.0.and.(cseis.eq.'1'.or.cseis.eq.'2'))then
       sa2=sa1*980.5
       read(TSA,'(f6.3)')ts
c The value sa2 is used to normalize grams; this is why sa1 had to be correctly interpreted.
c g to cm/s2
       if(cseis.eq.'1')then

      call td_drvrwww(xlong,ylat,rmode(iax),mmode(iax),sa2,sa5w,ts,per1,
     1      pid,ipid,ice,'MODAL',vs30)
       else
c determine stochastic seismograms from mean R,M pair rather than modal
       call td_drvrwww(xlong,ylat,rbarr(iax),mbarr(iax),sa2,sa5w,ts,per1,
     1      pid,ipid,ice,' MEAN',vs30)
       endif
       endif	!cseis=1 or 2 and iax=0
     
710      format(/,'Summary statistics for above ',2(a4,1x),
     +'deaggregation, R=distance, e=epsilon:',/,
     +'Contribution from this GMPE(%): ',f6.1,/,
     +' Mean src-site R= ',f6.1,' km;',
     1' M= ',f4.2,
     2'; eps0= ',f6.2,'. Mean calculated for all sources.'
     2 /,'Modal src-site R= ',f6.1
     3,' km; M= ',
     3f4.2,'; eps0= ',f6.2,' from peak (R,M) bin')
715      format('Modal-source dmetric: distance to rupture surface 
     1(Rrup or Rcd)')
717      format('Modal source dmetric: closest horizontal dist.',
     1/,'****** Further analysis of modal (M,D) properties: ******',
     2/,'Atten.Model. epsilon(SA0)   Median Motion  Median+1SD Motion')
999      continue

      write(1,775)
775      format(/,'Principal sources (faults, subduction, random seismicity',
     1' having >10% contribution)',/,
     2 'Source Category:',t33,' % contr.  R(km)    M   epsilon0 (mean values).')
      do k=1,nsrcmx
      pc=100.*ratecat(k,iax)/ratett(0)*arate/ratett(0)
      if(pc.ge.5.)then
      write(1,777)srcdesc(k),pc,rsav(k,iax)/ratecat(k,iax),amag(k,iax)/ratecat(k,iax),
     1      esav(k,iax)/ratecat(k,iax)
      endif
      enddo
	if(crustal(iax))then
      write(1,788)
      do k=1,lstft
c recompute rate on fault using most recently found rate for sources
      ratef(k,iax)=ratef(k,iax)*arate/ratett(0)      
      write(1,777)fltnm(k),ratef(k,iax),rbf(k,iax),embf(k,iax),ebf(k,iax),azi(k)
      enddo
      endif	!extra writes for crustal GMPEs only
777      format(a,2x,f6.2,2x,f6.1,2x,f5.2,2x,f6.2,2x,f8.1)
788      format('Individual fault hazard details if its contribution to mean hazard > 2%: ',/,
     1'Fault ID',t33,' % contr.   Rcd(km)  M   epsilon0 Site-to-src azimuth(d)')
      write(1,887)attnmodel(iax)
  887 format('#*********End of deaggregation corresponding to ',a,' *********#',/)
  	endif	!if 10% contributor.
  	enddo	!loop through attn models with 10% contribution.
      close(1)
      aratemn=0.05
	iax=0
c aratemn is now in units percent, so half a mil needed to plot.
      call gmtscript(gmt,rmax,emmin,emmax,pcmax,rbarr,mbarr,rmode,
     1      mmode,deltar,imax,jmin,
     1      jmax,icit,iclass,iax)
      end

      subroutine gmtscript(gmt,rmax,emmin,emmax,pcmax,rbarr,mbarr,
     1       rmode,mmode,deltar,imax,jmin,jmax,icit,iclass,iax)
c writes a gmt script. For the most part, run it silently. -V is allowed once,
c when the thing is about finished. Diagnostics are supposed to be going to a log file.
      parameter (deltam=0.2, amagmin=4.5, small=0.05)
      parameter (ndmx=100, nfmx=80,nmax=10000,nmmax=50,nsig=6)
c Now Jun 23, 2000: using loc(1:iloc) for locations of all files, newmaps no
c longer used.
      logical epson,maybe
      real, dimension(0:30):: mbarr,mmode,rbarr,rmode,pcmax,ratett,ebar0,emodeb0
      real m3d
      character*10 psa,rate
      character*14 e3dmsg
      character*4 label,labrs,per,retrnt*5,units,shift,unsh
      character*3 dya,xlab*6
c shift will shift labels depending on location of rmode(iax),mmode(iax)
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
       common/epsilon/epson,ebar0,emodeb0,eb0
      common/e3d/r3d,m3d,e3dmsg
c ebar0(iax) is the average of the epsilon0s, i.e., epsilons of SA0 or PGA0
      dimension gmt(nsig+2,ndmx,nmmax),eb0(ndmx,nmmax,0:30)
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
      dowhile(pcmax(iax).gt.y)
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
      pcmax(iax)=max(10.0,y)
      xsc=8./rmax
      ysc=4./(emmax-emmin)
      zsc=4./pcmax(iax)
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
      if(rmode(iax).gt.0.31*rmax.and.mmode(iax).gt.0.5*emmax)then
      shift='-2.8'
      unsh='2.80'
      elseif(rmode(iax).gt.0.2*rmax.and.mmode(iax).gt.0.8*emmax)then
      shift='1.60'
      unsh='-1.6'
      elseif(rmode(iax).lt..25*rmax.and.mmode(iax).gt.0.3*emmax)then
      shift='1.00'
      unsh='-1.00'
      elseif(rmode(iax).gt.0.8*rmax.and.mmode(iax).lt.0.5*emmax)then
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
      if(pcmax(iax).gt.99.)jpc=4
c 2008: change  "Source-site" to Closest distance in graph labels      
      write(1,12)rmax,emmin,emmax,pcmax(iax),xsc,ysc,zsc,xlab,
     1 dya,xw,yw,pidpr(1:ipidr)
12      format('$1psxyz -R0/',e10.5,'/',f4.2,'/',f4.2,'/0/',e10.5,
     1 '  -Jx',f5.3,'/',f5.3,' -Jz',f5.3,' \\',/,' -B',a6,
     2':"Closest Distance, Rcd (km)":/0.5:"MAGNITUDE (Mw)":/f1a',a3,
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
     1 gmt(2,i,jmax),gmt(3,i,jmax),eb0(i,jmax,iax)
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
     2      gmt(k,i,j),eb0(i,j,iax)
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
     1 per,psa(1:j1),ratett(iax),retrnt,rbarr(iax),mbarr(iax),ebar0(iax),
     1 rmode(iax),mmode(iax),emodeb0(iax),r3d,m3d,e3dmsg,deltar,unsh,pidpr(1:ipidr)
      elseif(epson)then
       write(1,880)shift,pidpr(1:ipidr),sitecl(iclass),
     + labrs,msg(1:icit),abs(xlong),ylat,
     1 psa(1:j1),ratett(iax),retrnt,rbarr(iax),mbarr(iax),ebar0(iax),
     1 rmode(iax),mmode(iax),emodeb0(iax),r3d,m3d,e3dmsg,deltar,unsh,pidpr(1:ipidr)
      elseif(per.eq.'0.00')then
       write(1,88)shift,pidpr(1:ipidr),sitecl(iclass),
     + labrs,msg(1:icit),abs(xlong),ylat,
     1 ratett(iax),retrnt,psa(1:j1),rbarr(iax),mbarr(iax),
     1 rmode(iax),mmode(iax),deltar,unsh,pidpr(1:ipidr)
      else
      write(1,48)shift,pidpr(1:ipidr),sitecl(iclass),
     1 labrs,msg(1:icit),abs(xlong),ylat,
     1 psa(1:j1),per,ratett(iax),retrnt,rbarr(iax),mbarr(iax),
     1 rmode(iax),mmode(iax),deltar,unsh,pidpr(1:ipidr)
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
      end  subroutine gmtscript    
      

c * --------------- BEGIN TD_DRVRWWW ---------------------------------
      subroutine td_drvrwww(xlong,ylat,rin,amagin,pshasa,sa5,
     1       ts,perin,pid,ipid,icus,msgsrc,vs30)
c 2008: we really need to get past the point source model!
c     www version: mods of apr 27 2001. S Harmsen. Any requests for
c      code should be forwarded to David Boore. Changes made here
c      result in output modifications that are specific to product
c      desired for deagg. hazard seismograms.
c *  Obtains input parameters and calls time domain simulation
c 2008: input vs30. This helps select the input file to Boore's smsim. 
c set up for hard rock in CEUS only as of Aug 20 2008.
c * Notes regarding modification of driver:
c *  1. None of "include" statements (here and in the other routines) contain
c *     the path "\smsim\"; the user must modify this appropriately, if the
c *     SMSIM program are placed in a directory with a different name and/or
c *     in a different drive than the drive from which the program is being
c *     compiled.
c *  2. The logical variable "new_mr" must be explicitly set (= .true. for a
c *     new value of magnitude and distance); its value is passed through the
c *     block common /misc/ (see SMSIM.FI).
c *  3. If values of stress other than those in the input parameter file are
c *     to be used, set "stressc" equal to the new value, rather than "stress"
c *     (and also set the other parameters if the stress changes with moment).

c * Dates: 06/01/95 - Written by D. M. Boore;
c *             Renamed and slightly modified PSV_DRVR
c *       06/09/95 - Pass r, amag through common rather than parameter lists
c *       06/12/95 - Added optional writing of acc, vel time series to a file
c *       08/11/95 - Added frequency column to psv output
c *       08/11/95 - Changed places where parameter-file and column-file names
c *             are requested
c *       08/18/95 - Added call to Write_Params
c *      1      0/17/95 - Added flag to get_params and write_params that tells 
c *             whether are dealing with time domain or rv.
c *      1      1/14/95 - Modified output slightly
c *      1     2      /08/95 - Added switch to use 91 standard periods
c *      1     2      /14/95 - Print out frequency and duration of excitation
c *      1     2      /17/95 - Changed location of asking for col file name and
c *             reorder some output
c *      1     2      /28/95 - Used new_mr to correct bug related to use of loop over
c *             amag,r
c *      1     2      /31/95 - Modified some of the *.sum output
c *       01/05/96 - Added column file stem name to column headings
c *       01/22/96 - Added total duration, npts to output
c *       02/06/96 - Minor formatting improvements
c *       04/12/96 - changed names psv, psa to psv, psa; changed
c *             format of time in accvel file to fixed format.
c *       02/28/97 - added output of kappa (because I now allow it to be
c *             magnitude dependent)
c *       04/20/97 - changed dimension of per from 91 to 400
c *       09/01/98 - added calculation of standard deviations
c *       09/11/98 - changed array dimensions from 16400 to 33000
c *       01/13/99 - added calculation of Arias intensity
c *       01/23/99 - modified output to include Sd, and changed
c *             "psv" and "psa" to "psv" and "psa"
c *       02/09/99 - Added computation of displacement 
c *       02/10/99 - Include option to print pga, pgd, pgd and similarly
c *             pga/w, pga/w*w, etc in a column file to provide peak
c *             motion levels in plots of response spectra.  The file
c *             has so many columns that I decided to write the
c *             information to a separate file rather than include it
c *             in the response spectra file.  
c *       02/13/99 - Changed "accvel" to "avd"
c *       03/05/99 - Added include statements at end and write smsim version
c *       03/12/99 - Change from 'a\' to 'a' in format statements.  This is
c *             less pleasing aesthetically, but it allows
c *             for redirection of input/output, which did not
c *             work before (a test program showed that redirection
c *             works as long as the total of the sum of the characters
c *             in the query strings written to the screen do not exceed 
c *             some undetermined number). Also, change "call banner(5)" 
c *             to "call banner(6)" and remove the "cr to quit" in the 
c *             "dist (cr to quit)" query
c *       03/12/99 - Added check for nper = 1 in loop over computing log-spaced
c *             periods (this is "fool-proofing"; if only one period is
c *             desired the user should reply "y" to the individual
c *             periods query and the computations with nper-1 in the 
c *             denominator should never be reached).
c *       02/05/00 - Used trim_c for f_out, f_col and corrected error in 
c *             making column heads for f_out stem less than 4 characters
c *       06/22/00 - Print out standard error of the mean rather than std
c *       01/28/01 - Increased array size from 33000 to 66000 (should use
c *             allocatable arrays)
c *       02/12/01 - Use dynamically allocatable arrays; renamed "npts" to "npw2" (SUNS NO)
c *         Apr 2003 - consider using AB95 for CEUS param file. However, how???
      parameter (nsavmx=60)
c nsavmx = max number of seismograms to save; this is an array dimension
       character*80 f_params, f_params2,f_out,fpar_e,loc,message
      character  msg4*4/'PSA='/,date_begin*10,pid*12,msgsrc*6
      character*4 label,perx,retrnt*5,units,comm*16,msg*24
      common/lab/msg,label,perx,retrnt,units,comm
      logical tdflag,icus
c icus=.true. if > 45% exceedances from gridded CEUS sources  
      integer indx(nsavmx)
      real xlong,ylat,pshasa,perin,sa5
       real acc_save(66000,nsavmx)
      real targpsa(10)
      real psvsim_m1(nsavmx)
      REAL PER(10)
      common/loca/loc,iloc
      common/lab2/psa,rate,xlon,yla,metric
      character*10 psa,rate,metric*20

       include 'smsim.fi'

      message =
     1      'PSHA INTERACTIVE DEAGG RVT TIME SERIES FOR '//
     1      msgsrc//' (R,M)'
c      read(*, '(a)') message
    
333   continue
      amag=amagin
      r=rin
c 'E_BCK01L.DAT' is Boore-Frankel CEUS BC parameter file for smsim_td. 
50      format(a,$)
	if(vs30.lt.1300.)then
	fpar_e='SM/E_BCK01L.DAT'
	else
	fpar_e='SM/4ARTK006.DAT'
	endif
      if(xlong.gt.-105.5.or.amag.lt.5.0.or.
     1      (icus.and.amag.lt.6.5))then
c amag < 5 means that MbLg in the 5. to 5.3 range converted to Mw<5. Low M
c sources that determine the mode must have been handled using the Boore-Frankel
c CEUS attenuation model in Map96. 
c For some cases it is appropriate to use the Boore Frankel. model here, too. The
c low mag occurs as far west as the IMSB. Also if icus=true and amag<6.5,
c the Boore Frankel model is appropriate. The constraint amag<6.5 is intended
c to eliminate possibility that fault hazard using WUS attenuation is the
c controlling event. 
c 2002: the AB95 file with BC siteamp might be appropriate for many sources to
c get the 2-corner source spectrum effect. How to include?
c The longer period SA will be low compared to probabilistic SA if using AB95.
C Could scale grams amplitude up as is done for other inputs, but this will cause
c high freq. part to scale up too, by too much. Needs a consensus approach.
c Maybe should not scale up, but go with low amps of AB95 for longer period. This was
c their intention anyway. Herrmann suggests, use larger M or closer R when AB95 in use.
c Mhat = 7.4 has about same mu for FEA as M=8.0 for AB95, 1s SA.
c There is  a mixture of WUS/CEUS contribs.
c for WUS sites when icus=.true., but the WUS portion should be less than CEUS.
      f_params = fpar_e
c      f_params2 = 'SM/AB95.BC.DAT'
      elseif(xlong.gt.-110.and.ylat.gt.40.)then
       f_params = fpar_e
      elseif(amag.gt.8.1.or.
     1      (ylat.gt.41..and.xlong.lt.-120.))then
c Atkinson-Boore '97 designed for moderate mag (M<6.5) PNW seismicity
c We use it quite arbitrarily for big AK & Cascadia subduction as well. 
      f_params = 'SM/AB97-SU.DAT'
      elseif(r.gt.100..or.
     1 amag.gt.7.75)then
      f_params = 'SM/AS97-CA.DAT'
      elseif(amag.le.7.75)then
      f_params = 'SM/WR032496.DAT'
      else
      call get_lun(nout)
      f_out=loc(1:iloc)//'smsim.acc.'//pid(1:ipid)
c All WWW I/O must be to locations with full path name specified.
       open(nout,file=f_out,status='unknown')
      write(nout,950)rin,amagin,ylat,xlong
      close(nout)
c Also a default postscript file needs to be published. Otherwise breakdown
c of communications with Nancy Dickman's driver. Copy western hemisphere,
c west.hem.ps, to the default name.
            f_out=loc(1:iloc)//pid(1:ipid)//'grams.ps'
      call system('cp '//loc(1:iloc)//'PL/west.hem.ps '//f_out)
      return
950      format('RVT time series not available for this set of params:'
     1 /,' R ',f6.1,' M ',f6.2,' site loc: ',2(1x,f8.3))
      endif
c      f_out = ' '
c      write(*, '(a)') 
c     :    ' Enter name of summary file: '
c      read(*, '(a)') f_out
      f_out=loc(1:iloc)//'smsim.acc.'//pid(1:ipid)
      call get_lun(nout)
      open(nout,file=f_out,status='unknown')

      call banner(nout)

      write(nout, '(a)') message
      write(nout,'("Site location: ",f6.3,1x,f8.3)')ylat,xlong

c      write(nout, '(a)') 
c     :    ' *** Results computed using TD_DRVRWWW ***'

      call get_date(date_begin)
      write(nout, '(2a)') ' Date: ', date_begin
c      call get_time(time_begin)
c      geotechnical new.
	write(nout,'(a,1x,f6.1)')' Vs30(m/s) = ',vs30
      write(nout, '(2a)') ' file with parameters: ',
     :   f_params(1:20)

      damp=0.05
      nper=1
      per(1)=perin
      targpsa(1)=pshasa      
      tdflag = .true.  ! controls extent of the input file read and written
      f_params = loc(1:iloc)//f_params
      call get_params( f_params, tdflag )
      if(perin.ne.0.)then
      seed=seed+10./perin
c      print *,'seed ',seed
      endif
      call write_params(nout, tdflag)
      
100   continue
      new_mr = .true.
       nacc_save=min(nsavmx,nsims)
c      

      call Get_Npts(nstart, nstop, npw2, te)  ! included in smsim_td
      ndimen = 66000
      if (npw2.gt. ndimen) then
       write(*,'(2x,a,i5,a,i5,a)') 
     :       ' *** FATAL ERROR: NPW2 (',
     :                      npw2, 
     :       ') > NDIMEN (', 
     :                       ndimen, 
     :       '); QUITTING!'
       write(*, '(2x,a)') ' TRY INCREASING DT'
       stop
      end if
c * ALLOCATABLE: end of group
       if(per(1).eq.0.)then
       nacc_tmp=nacc_save
       else
       nacc_tmp=60
c we will compute 60 psrvs choose the six with lowest "l1" variability
       endif
       do isim=1,nacc_tmp
       indx(isim)=isim
       psvsim_m1(isim)=0.0
       enddo
c smsim_tdww2 is a more recent version of smsim_tdwww. Some question about avgpsa
c 8/2008.
      call smsim_tdww2(per,nper,targpsa,sa5,ts,psvsim_m1, 
     :   pgasim, pgasim_std, avgpsa,
     :   avgscafac, nacc_tmp, acc_save)

c      if (save_avd .eq. 'Y') then
c       call get_lun(nuacc)
c       unit=nuacc,file=f_avd,status='unknown')
      if(perin.gt.0.)then
       call sort2(nacc_tmp,psvsim_m1,indx)
       write(nout,'(2x,a,f6.3,a,f6.3) ')
     :' Fractional oscillator damping = ',
     :damp,' period= ',perin
      write(nout,'(2x,a,f7.4)')'**Agrams scaled to PSHA SA
     1 level. Avg scale factor = ',avgscafac
      else
      write(nout,404)avgscafac
404      format('Agrams scaled to PSHA HPGA. avg factor=',f7.4)
      endif
      write(nout, '(a7,1x,a6,a, 1p3(1x,e10.3))')
     : '**PSHA ',msgsrc,' r(km), amag, PSHA sa or pga(cm/s2)= ',
     : r, amag, pshasa
      write(nout,'(a,2(1x,a4))')
     1 '**PSHA Ground motion mean return time ',retrnt,units 
      if(perin.gt.0.)then
      write(nout,'(a,f8.3,a,/,a,f5.3,1x,f5.3,a)')'** PSHA approximate  
     1UHS target ordinate at 0.2 s is ',sa5,' cm/s/s',
     2'** PSHA approximate UHS short- and long-period corners ',
     30.2*ts,ts,' seconds.'
      endif
      write(nout, '(2x,a,f10.6)') ' kappa= ', kappa_f(amag)
      write(nout, '(2x,a,1x,1pe10.3)') ' const= ', const
      write(nout, '(2x,a,f6.3,1pe9.2,2e10.3,e9.2)') 
     :  ' amag, stress, fa, fb, durex= ', 
     :    amag, stress, fa, fb, durex
      write(nout, '(2x,a,1p2(1x,e10.3))') ' am0, am0b_m0fa= ', 
     :                      am0, am0b_m0
      write(nout, '(2x,a,i6, f8.5, f6.1)')
     :  ' npw2, dt, total duration = ', npw2, dt, npw2*dt
c      write(nout,'(t5,a, t16,a, t25,a, t38,a, t46,a, t60,a)') 
c     :    'pgd(cm)',    'sem/pgd', 
c     :    'pgv(cm/s)',  'sem/pgv', 
c     :    'pga/(cm/s2)', 'sem/pga'
c      write(nout,'(1p,6(1x,e10.2))') 
c     :     pgdsim, (pgdsim_std/sqrt(float(nsims)))/pgdsim, 
c     :     pgvsim, (pgvsim_std/sqrt(float(nsims)))/pgvsim, 
c     :     pgasim, (pgasim_std/sqrt(float(nsims)))/pgasim
c      write(nout,'( t2,a, t25,a)') 
c     :    'Arias intensity(cm/s)', 'sem/Arias'
c      write(nout,'(1p, 11x, 2(1x,e10.2))') 
c     :     arias, (arias_std/sqrt(float(nsims)))/arias

c      write(*,*)

c      if (psvcalc .eq. 'Y') then
c       write(nout,'( t5,a, t17,a, t22,a, t35,a, t43,a, t57,a)')
c     :         'per(s)', 'freq', 'psv(cm/s)', 
c     :         'sem/psv', 'psa(cm/s2)', 'psd(cm)'
c       do i = 1, nper
c       write(nout,
c     :      '( t4,f7.3, t13,f8.3, 1p, t22,e9.2, t32,e10.2,
c     :       t43,e10.2, 1x,e10.3)')
c     :      per(i), 1.0/per(i), psvsim(i), 
c     :      (psvsim_std(i)/sqrt(float(nsims)))/psvsim(i),
c     :      (twopi/per(i))*psvsim(i), psvsim(i)/(twopi/per(i))
c       write(ncol,'(1x,f6.3,1x,f8.3,1p,t19,e10.3,t32,e10.3,
c     :            t45,e10.3, t62,e10.3)')
c     :      per(i), 1.0/per(i), 
c     :      psvsim(i)/(twopi/per(i)), psvsim(i), 
c     :      (twopi/per(i))*psvsim(i), 
c     :      (psvsim_std(i)/sqrt(float(nsims)))/psvsim(i)
c       end do
c       close(unit=ncol)
c       if (write_avd .eq. 'Y') then
c       do i = 1, 2
c        write(npavd, '(3(0p1x,f7.3,1p3(1x,e10.3)))')
c     :       t_pgd(i),pgdsim,pgdsim*(twopi/t_pgd(i)),
c     :           pgdsim*(twopi/t_pgd(i))**2,
c     :       t_pgv(i),pgvsim/(twopi/t_pgv(i)),pgvsim,
c     :           pgvsim*(twopi/t_pgv(i)),
c     :       t_pga(i),pgasim/(twopi/t_pga(i))**2,
c     :           pgasim/(twopi/t_pga(i)),pgasim
c       end do
c       close(unit=npavd)
c       end if
c      end if

      write(nout, *)       
c      call get_time(time_stop)
      write(nout, '(2x,a)') '*** Begin Scaled Accelerogram Data ***'
      if(amag.le.6.)then
      kpts=min(int(  50./dt)+1,npw2)
      elseif(amag.le.7.1)then
      kpts=min(int(  64./dt)+1,npw2)
      elseif(amag.lt.7.5)then
      kpts=min(int(  90./dt)+1,npw2)
      elseif(amag.le.8.3)then
      kpts=min(int(  160./dt)+1,npw2)
      else
      kpts=min(int(  300./dt)+1,npw2)
      endif
      accmax=0.0
      write(nout, '(t7,a,t13,a,t32,a,t43,a,t55,a,t67,a,t79,a)')
     :    'T', 'A1(cm/s2)', 'A2', 'A3','A4','A5','A6'
      nsav=min(nsavmx,nacc_save)
      tmax=kpts*dt*0.95
      i2=tmax/dt +1
      if(amag.gt.6.4.and.amag.lt.6.7)then
      i1=nint(20./dt)
      t_s=20.0
      elseif(amag.lt.7..and.dt.gt.0.007)then
      i1=nint(25./dt)
      t_s=25.0
      else
      i1=nint(15./dt)
      t_s=15.0
      endif
       do i = 1, kpts
       tnow=float(i-1)*dt
c on maxwell the below 1p6 replaces a 1p<nsav> that is used on SUNs
c I dont dknow why the <nsims> is a fatal error on maxwell.       
       write(nout, '(f8.4,1x,1p6(1x,e11.4))') 
     :       tnow, (acc_save(i,indx(j)),j=1,nsav)
      if(i.gt.i1.and.i.lt.i2)
     1      accmax=max(accmax,abs(acc_save(i,indx(1))),
     2      abs(acc_save(i,indx(nsav))))
       end do
c       close(unit=nuacc)
c      end if 
       
c * ALLOCATABLE: comment out following lines to disable dynamic allocation
c      deallocate (acc_save, vel_save, dis_save )
c * ALLOCATABLE: end of group

      write(nout, '(a)') ' *********** END OF AGRAM DATA***********'
      ipgamin=int(accmax)
      ipgamax=int(11.*accmax)
      f_out=loc(1:iloc)//'smsim.gmt.'//pid(1:ipid)
c      print *,'about to open ',f_out
       close(unit=nout)
       call get_lun(nout)
      open(nout,file=f_out,status='unknown')
      f_out=loc(1:iloc)//pid(1:ipid)//'grams.ps'
      do j=1,nsav
      if(j.eq.1.and.tmax.gt.160..and.accmax.lt.500.)then
      write(nout,200)int(t_s),int(tmax),ipgamin,ipgamax,
     1      ' -JX6/9 -Ba50f5/a100f10WSen -K -W2
     1       -P +$2PL/wwwdefaults << END > '
     2      ,f_out
      elseif(j.eq.1.and.tmax.gt.160.)then
      write(nout,200)int(t_s),int(tmax),ipgamin,ipgamax,
     1      ' -JX6/9 -Ba50f5/a500f50WSen -K -W2
     1       -P +$2PL/wwwdefaults << END > '
     2      ,f_out
      
      elseif(j.eq.1.and.accmax.lt.500.)then
      write(nout,200)int(t_s),int(tmax),ipgamin,ipgamax,
     1      ' -JX6/9 -Ba10f1/a100f10WSen -K -W2
     1       -P +$2PL/wwwdefaults << END > '
     2      ,f_out
      elseif(j.eq.1.and.accmax.lt.2000.)then
      write(nout,200)int(t_s),int(tmax),ipgamin,ipgamax,
     1      ' -JX6/9 -Ba10f1/a500f100WSen -K -W2
     1       -P +$2PL/wwwdefaults << END > '
     2      ,f_out
      elseif(j.eq.1)then
      write(nout,200)int(t_s),int(tmax),ipgamin,ipgamax,
     1      ' -JX6/9 -Ba10f2/a5000f500WSen -K -W2
     1       -P +$2PL/wwwdefaults << END > '
     2      ,f_out
      else

      write(nout,200)int(t_s),int(tmax),ipgamin,ipgamax,
     1      ' -JX -B -K -O -W +$2PL/wwwdefaults << END >> ',f_out
      endif
      t=t_s
      off=2.*(j-1)*accmax
      do i=i1,i2
      write(nout,'(f8.4,1x,e11.5)')t,acc_save(i,indx(j))+off
      t=t+dt
      enddo
      write(nout,'("END")')
      enddo
c negative sign outside of i4.4 format which is designed for non-neg. integers only.
200	format('$1psxy +PL/wwwdefaults -R',
     +i3.3,'/',i4.4,'/-',i4.4,'/',i7.7,a,a)
      if(perin.eq.0.)msg4='PGA='
      write(nout,300)f_out,msgsrc,ylat,xlong,rate,perin,msg4,psa,r
     1      ,amag,kappa_f(amag)
300      format('$1pstext -R0/8/0/10 -Jx1/1  -O -S +$2PL/wwwdefaults << END >> ',a,/,
     1      '1.4 9.30 12 0 0 1 Scaled Accelerograms for ',a,' Source',/,
     1      '0.1 9.1 10 0 0 1 Loc: ',f6.3,1x,f8.3,' rate ',a,' per=',f4.2,
     2      ' ',a4,a,' g. From SMSIM_TD',/
     2      '0.25 8.75 10 0 0 1 (R,M)=',f6.1,' km,',f5.2,/,
     2      '4.75 8.75 10 0 0 1 kappa=',f6.4,' sec',/,
     2      '0.25 0.8 10 0 0 1 A1',/,
     2      '0.25 8.35 10 0 0 1 A6',/,
     2      '0.25 2.6 10 90. 0 1 cm/s/s',/,
     3      'END',/,'exit')
      close(unit=nout)
c      write(*,*)
c      write(*, '(a)') 
c     : ' Compute results for another r and M (y/n;cr=quit)? '
c      read(*, '(a)') buf
c      if (buf(1:4) .eq. '    ') go to 999
c      if (buf(1:1) .eq. 'n' .or. buf(1:1) .eq. 'N') go to 999
c
c      goto 100

999   continue      
      return
      end
c * --------------- END TD_DRVRWWW ---------------------------------

      function len2b(name,lmax)
c finds string length up to first encounter of blank-blank or lmax.
      character*32 name
      i2=4
      dowhile(name(i2:i2+1).ne.'  '.and.i2.lt.lmax)
      i2=i2+1
      enddo
      len2b=i2
      return
      end
