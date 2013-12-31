c sum_haz_2013.f used on the WWW computer. Compile with f77 if SUN. 
c dec 2 2013: change aratemn to 1 mil (0.1 %) from 1/2 mil. SH.
c Last update Nov 26 2013 GMPE indexes and weights for 2014 NSHMP model.
c May 2013. multiple GMPE contributions may be reported. Make arg6 1 to get report
c   for each atten model as well as mean hazard report.
c sept 24 2013: kbc index needed to be a vector. "e3dpl" from "e3dmsg" for mean haz
c jan 7 2013. minor revs to initialization
c --- this pgm is set up for NSHMP 2013 update; it may need work in these areas:
c --- (1) Removed seismograms for modal event or mean event based on RV point src.
c --- (2)  NMSZ descriptions
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
c plot modal R* M* e* in both ascii and graphic output files
c
c No makefile. use gnu fortran
c PC Linux compile (In 2008, this is gldplone): 
c gfortran ./sum_haz_2013.f -static -O \
c     -o ./sum_haz_2013 -ffixed-line-length-none
c Solaris compile: 
c f77 sum_haz_2013.f  -e -fast -o ~harmsen/bin/sum_haz_2013
c  note: -Bstatic with f77 or f95 compilation goes to pieces in a hurry, Solaris System objects to -Bstatic. the
c /usr/lib/libmvec.so is shared. 
c  Mystery: compiling with f95 gives us unplottable gmt script but f77 works. Why? 
c
c Makefile links up several of Boores subroutine packages & iosubs.o. Use -fast.
c The -static argument invokes static library linking, needed for Web computers.
c
c Usage: sum_wt_haz argj  j=1,...,13
c
c arg1 = site name (no blanks!)
c arg2 = site longitude (E+, dec degrees)
c arg3 = site latitiude (dec degrees)
c arg4 = target Return period (a4) units: years or kiloyears
c arg5 = vs30 m/s
c arg6 = iattn. iattn=1 to get individual gmpe reports new may 2013.
c arg7 = PSA (PGA) period. PGA convention is 0.00 sec.
c arg8 = PSA (PGA) level, units: g
c arg9 = computed annual rate of exceedances (exponential notation please)
c arg10 = pid (process id, for names)
c arg11 = absolute path name of data files  Out
c arg12 = location of input files (some need to distinguish these)
c arg13 = revised weight for clustered.. work in progress nov 8 2013
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
c When OKed, sum_haz_2013 is for automated web-site PSHA deaggregations.
c Input Magnitudes assumed to be greater than 4.5 (Mw). For lower magnitudes,
c some code modifications are required. The 4.5 limit corresponds to a 5.0(MbLg)
c limit in the 2002 National Seismic Hazard Maps. .
c Two output files: name specified and gmt.period (period is 4 char SA period)
c Steve Harmsen updated 6/2001 to 8/2008
c
      parameter ( deltam=0.2, amagmin=4.5,nsrcmx=43,nattmx=30)
      parameter (ndmx=100, nfmx=180,nmax=10000, nsig=6,nmmax=50)
c nmax=max number of records per file. Be careful here!
      logical avgctr,epson,getwt,header,ispga,aok,bypass
      logical fltok,ice,countme
      real ratecat(0:43,0:43)
      real, dimension(0:100,0:43):: rlow
      real d,eps,m,mag,m3d,m3dpl,hmax/0.0/
      real rmax/10./,deltar/25./,h_int,h_intm/0.0/
      character*64 name1,name2,a32*32
      character*1 ch1
      character*32 srcdesc(0:43)/'sub0(full) Cascadia Megathrust',
     +'Cascadia M8.0-M8.7 Floating','S. Cascadia M8.12-M8.52 Char.','sub1(south) Cascadia Megathrust',
     +'Ch. faults Geodetic rate Bird','G-R faults Geodetic rate Bird',
     +'Ch. faults Geodetic rate Zeng','G-R faults Geodetic rate Zeng',
     +'Ch. faults Geologic rate WC-SRL','G-R faults Geologic rate WC-SRL',
c going into index 10
     +'Extensional gridded ad_Sm Mx7.45','Extensional gridded ad_Sm Mx7.95',
     +'Extensional gridded fixSm Mx7.45',
     +'Extensional gridded fixSm Mx7.95','Wasatch Segmented GR','Wasatch Unsegmented M7.4',
c going into index 16
     +'Compr gridded adSm Mx7.45','Compr  gridded adSm Mx7.95','Compr gridded fixSm Mx7.45',
     +'Compr gridded fixSm Mx7.95','CA-NV Shear Zones',
c going into index 21
     +'New Madrid SZ USGS clustered','Deep Intraplate w/Stairs','Puget Lowlands gridded',
c going into index 24
     +'Wasatch SLC clustered','Wasatch SLC Pechmann','sub4(north) Cascadia Megathrust','Charl. SC Regional Zone',
     +'Compressional gridded ad_Sm Mx7.45', 'Compressional gridded ad_Sm Mx7.95','Compressional gridded fixSm Mx7.45', 
c going into index 31
     +'Cheraw CO or Meers OK faults','New Madrid SZ SSC no cluster','CEUS gridded',
     +'Charleston SC Local zone','Charleston SC Narrow Zone','SSC-Marianna RLME',
c going into index 37
     +'SSC-Wabash RLME','SSC-Commerce RLME','SSC ERM-N and S RLME',
c going into index 40
     +'SSC Charlevoix RLME','CA UCERF2 a-faults','CA UCERF2 b-faults','NMSZ USGS'/
       character*52 amsg(2)/' (some WUS atten. models use Site Class not Vs30).',
     +' CEUS atten. model site cl BC(firm) or A(hard).'/
      character*24 attnmodel(0:38)
     +/'Mean Hazard w/all GMPEs','Spudich et al. Ext','Toro et al. 1997',
     +'Sadigh et al. 1997','Atkinson-Boore06,140 bar','Attn model 5',
     +'Frankel et al., 1996','Somerville Rifted FinFlt','Attn model 8',
     +'Attn model 9','Campbell CEUS Hybrid','Attn model 11',
c index 12 to 14
     +'Atkinson-Boore03 Casc.','Youngs et al. 1997','Atten Model 14',
c index 15 to 18
     +'Silva 1-corner','Atten Model 16','AB03 Global Intraplate','Atkinson-Boore03 Glbl.',
     +'Tavakoli and Pezeshk 05','Atkinson-Boore06,200 bar',   
     +'Boore-Atkinson 2008','Campbell-Bozorgnia 2008',
     +'Chiou-Youngs 2008','Abrahamson-Silva 2008','Pezeshk 2011',
c index 26 to 29:
     +'Crouse 1991','Zhao et al. 2006','BCHydro Subd 2012','Atkinson-Macias Subd',
c index 30 to 33
     +'Atkinson 08prime','BCHydro Deep 2012','AB06prime','BSSA NGA2013',
     +'Campbell-Bozorg NGA2013',
     +'Chiou-Youngs NGA2013','Abr-Silva-Kamai 2013','I.M. Idriss NGA2013',
     +'Graizer-Kalkan 2013'/
c the gmpe list needs updating 2013
      character*8 date,time*10,zone*5
      character*10 avs30
      character*7 faults/'faults/'/
      character*14 e3dmsg,e3dpl
      character*16 comm/'2013 edition. '/
      character*12 pid
      character*10 psa,rate, slong,slat,metric*20
      character*4 labs,label,per,retrnt*5,units/'yrs.'/,anum      !temp hold record fault number
      character*24 msg
      character*8 dum8
      common/lab/msg,label,per,retrnt,units,comm
      common/vs30/vs30
      common/pid/pid,ipid
      common/loca/loc,iloc
      common/lab2/psa,rate,xlong,ylat,metric
       common/lowlim/aratemn,ratett
      common/epsilon/epson,ebar0,emodeb0,eb0
      common/e3d/r3dpl,m3dpl,e3dpl
       character*4 epr(10)
     1/'5.00','4.00','3.00','2.00','1.00','0.50','0.30','0.20','0.10','0.00'/  
      character*80 loc,rec,logdir
      character*32 fltnm(43)
c fltnm will hold short descriptive message (Char or GR; dip of fault)
c about significant faults.
       character*1 cseis,flg*2,flg3*3
      logical maybe,crustal(0:43)
c crustal: faults and subduction "true" gridded "false"
c indexes 0 to 18
	data crustal/10*.true.,4*.false.,4*.true.,.false.,
c 19 to 43
     + 7*.true.,4*.false.,.true.,.false.,12*.true./
      real, dimension(nfmx):: azi
      real, dimension (0:nfmx,0:43):: amag,rsav,esav
      integer k00/1/,k01/1/,k02/1/,ival(8)
      integer iloc,ireg,kf,kfmin,kfmax,kflim/1/,nocc,nvac
c several scalars had to be converted to vectors. kbc is an index for
c a "Bazzuro-Cornell" modal triple, in magnitude, distance, and epsilon.
c kbc(0) corresponds to the mode from the mean hazard. kbc(j) for j>0 correspond
c to the mode for specific attenuation models. Plot is just for the mean hazard.
      integer, dimension(0:43) :: kbc
      real, dimension(0:200,0:43) :: flthaz,fm,fd,feps
c flthaz(0) hazard for non-faults, otherwise for a fault having an index
       real tmp
       real, dimension(0:43) :: avgr, avgm,eps1,pr50,rbarr,mbarr,ebar0,emodeb0,
     +  ratelow,ratett,ratesum,rmode,mmode,pcmax,gmpewt
c the array gmt needs to be upped to 4D to be able to plot individual GMPE
c data with the GMT subroutine. This was not done in initial code preparation.
c SHarmsen May 2013.
      real h(nsig),gmt(nsig+2,ndmx,nmmax)
c awesome 4d arrays:
      real, dimension(nsig,ndmx,nmmax,0:43) :: hnew,mnew,rnew
c  mnew, rnew will store avg values in each 3D bin      
      real hsum(nsig,0:43),emmin/10./,emmax/1./
      real, dimension(80,0:43):: rbf,embf,ebf,ratef
      real, dimension(ndmx,nmmax,0:43):: rbar,eb0, mbar,ebar
c gmpewt needs work nov 2013: looks like 2008 vintage
c indexes 0 to 9. Example: 7 corresponds to Somerville Rifted wt 0.1
      data gmpewt/1.,.1,.11,.2,.1,.1,.06,0.1,.2,.1,
c 10 to 19. 17 is AB03 global intraslab. 18 is AB03 global interface
     +.11,.1,.1665,.25,.5,.06,.1,0.1665,.25,.11,
c 20 to 26
     +0.1,0.1,0.3333,0.3333,0.3333,0.15,0.1,
c 27 to 32
     +0.3333,0.3333,0.1,0.08,0.3333,0.22,
c 33 to 37 These indexes are NGA-W(2) BSSA to IM Idriss wts
     +4*0.22,0.12,
c 38 to 43
     +0.10,0.,0.,0.,0.,0./
c vector initialization
      rbarr=0.0
      ratett=0.0
      ebar0=0.0
      mbarr=0.0
      mmode=0.0
      rmode=0.0
      pcmax=0.0
      rsav=0.
      amag=0.
      esav=0.
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
      comm=' Exceedances.'
      if(xlong.lt.-129.)comm='From OFR 99-36. '
44      call getarg(3,slat)
      read(slat,'(f7.4)',err=46)ylat
      fac=(ylat-38.)/(xlong+121.75)
      if(ylat.lt.41.0.and.ylat.gt.40.and.xlong.gt.-75.)then
      deltar=10.
       comm='New York area'
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
       elseif(abs(ylat-32.8).lt.0.9.and.abs(xlong+80.).lt.0.9)then
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
      read(ch1,'(i1)')iattn	!new item 6 in 2013
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
c Might use ipr to gather attenuation coefficients (not in use 2008-2013)
       endif
       aok=per.eq.epr(ipr)
c aok "true" means that program knows what period was input, otherwise, we
c  cant do further attenuation-model epsilon analysis.
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
      if(arate.gt.0.)ireturn=10*nint(0.1/arate)
      print *,arate,ireturn
      if(ireturn.lt.10000)then
      write(retrnt,1857)ireturn
      else
      write(retrnt,1859)ireturn
      endif
1857	format(i4)
1859	format(i5)
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
c arg 13 is a revised weight for clustered event which depends on another model. Testing Nov 8
c as an alternative to rewriting deaggFLTH.2013.
       call getarg(13,dum8)
       read(dum8,'(f8.5)')wt_rev
      else
      print *,'Usage sum_haz_2013 site long lat rate vs30 iattn per gmlev rate2'
      endif
c      print 50,' How many files to average : '
5      format(a)
50      format(a,$)
      read *,nfil
      if(nfil.gt.nfmx)stop ' sum_haz_2013: 180 files max. Bump & rerun'
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
c      print *,name1,wt,key
c key is a source-category key, added Jan 2003 to combine similar-source hazard.
c in 2013 there are 44 such categories starting with key=0 (full cascadia megathr).
      key=min(key,43) 
      if(key.eq.21 .or. key.eq.24)then
      print *,wt,wt_rev
      wt = wt_rev	!from arg 13. Reweigh clustered NMSZ  or clustered SLC-Wasatch     
	endif
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
c Add back-azimuth field (averaged over rupture models) Mar 2013.
	countme=.false.
      do ift=1,600
1357	continue
      read(3,27,end=16,err=1357)rbft,embft,ebft,pr,backaz,a32,iax
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
      xmlim=max(m,4.5)	!put anything smaller into 4.5 bin new nov 8 2013
      jbin=1+(xmlim-amagmin)/deltam
      if(jbin.gt.nmmax)stop 'Sum_haz_2013: a mag exceeds array bound'
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
     + print *, name1(1:32),' annual rate: ',ratesum(iax), ' cumu ',ratett(iax),iax
        enddo	!iax
c      ratepc=ratesum(0)/arate*100.0
      enddo      !do i=1,nfil
	iax=0	
      read 5,name2
      i=3
      dowhile(name2(i:i).ne.'.')
      i=i+1
      enddo
      name1=name2(1:28)
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
      call date_and_time(date,time,zone,ival)
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
      write(1,4)msg(1:icit),abs(xlong),ylat,vs30,amsg(ireg),date,retrnt,
     1per,label,psa(1:9),ratett(iax),pr50(iax),attnmodel(iax)
           endif
4      format('PSHA Deaggregation. % Contributions. Site: ',a,
     1' long: ',f7.3,' d W., lat: ',f6.3,' N.',/,'Input Vs30 (m/s) =',f6.1,a,/,
     1'NSHMP 2014 update. See URL earthquake.usgs.gov/hazards/2014prelim. Analysis on DaMoYr:',
     1a10,/,'Mean Return Period: ',a5,' years. ',a4,' s. ',a4,'=',
     1a9,'g. Weight * Computed_Rate_Ex ',e9.3,/,
     1'#Pr[at least one eq with median motion>=PSA in 50 yrs]=',f7.5,/,
     1'#This deaggregation corresponds to ',a,/,
     1'DIST(km) MAG(Mw) ALL_EPS EPSILON>2  1<EPS<2 0<EPS<1 -1<EPS<0 -2<E',
     2'PS<-1 EPS<-2')
704      format('PSHA Deaggregation. %contributions. site: ',a,
     1' long: ',f7.3,' W., lat: ',f6.3,' N.',/,' Vs30(m/s)=',f6.1,a,
     1/,'NSHMP URL earthquake.usgs.gov/hazards/2014prelim. dM=0.2 below'
     1/,'Return period: ',a5,' yrs. Exceedance PGA =',
     1a9,'g. Weight * Computed_Rate_Ex ',e9.3,/,
     1'#Pr[at least one eq with median motion>=PGA in 50 yrs]=',f7.5,/,
     1'#This deaggregation corresponds to ',a,/,
     1'DIST(KM) MAG(MW) ALL_EPS EPSILON>2  1<EPS<2 0<EPS<1 -1<EPS<0 -2<E',
     2'PS<-1 EPS<-2')
      aratemn=ratett(iax)*0.0010
c aratemn is 0.1% of the mean or GMPE hazard (1-mil). Dont output if less than that.
      do j=1,jmax
      m=amagmin+(j-0.5)*deltam
      m=min(8.,m)
       do i=1,imax
      d=deltar*i-deltar*0.5
      if(hnew(1,i,j,iax).gt.aratemn)then
      dist=rbar(i,j,iax)/hnew(1,i,j,iax)
      epsilon=eb0(i,j,iax)/hnew(1,i,j,iax)
      eb0(i,j,iax)=epsilon
      if(iax.eq.0)rmax=max(dist,rmax)
      mag=mbar(i,j,iax)/hnew(1,i,j,iax)
      emmax=max(emmax,mag)
      emmin=min(emmin,mag)
c      if(rmax.eq.dist.and.rmax.gt.250)print *,rmax,mag,hnew(1,i,j,iax)
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
      kbc(iax)=k
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
      if(kbc(0).eq.1)then
      e3dmsg='> 2 sigma'
      elseif(kbc(0).eq.2)then
      e3dmsg=' 1 to 2 sigma'
      elseif(kbc(0).eq.3)then
      e3dmsg=' 0 to 1 sigma'
      elseif(kbc(0).eq.4)then
      e3dmsg=' -1 to 0 sigma'
      else
      e3dmsg='-2 to -1 sigma'
      endif	!kbc(iax) business
      r3d=rnew(kbc(iax),ibc,jbc,iax)/h_intm
      m3d=mnew(kbc(iax),ibc,jbc,iax)/h_intm
      if(iax.eq.0)then
      r3dpl=r3d; m3dpl=m3d
      e3dpl=e3dmsg
      endif
c write the percent contribution to total hazard, not conditional based
c on a specific GMPE (or iax). Total hazard has iax=0.
           write(1,910)r3d,m3d,e3dmsg,h_intm*100./ratett(0)
910      format(' MODE R*=',f6.1, 'km; M*=',f5.2,
     1 '; EPS.INTERVAL:',a14,' % CONTRIB.= ',f6.3)
       endif
       if(mmode(iax).gt.8.)write(1,715)
       if(iax.eq.0.and.(cseis.eq.'1'.or.cseis.eq.'2'))then
       sa2=sa1*980.5
c The value sa2 is used to normalize grams; this is why sa1 had to be correctly interpreted.
c g to cm/s2
c       if(cseis.eq.'1')then

c      call td_drvrwww(xlong,ylat,rmode(iax),mmode(iax),sa2,sa5w,ts,per1,
c     1      pid,ipid,ice,'MODAL',vs30)
c       else
c determine stochastic seismograms from mean R,M pair rather than modal
c       call td_drvrwww(xlong,ylat,rbarr(iax),mbarr(iax),sa2,sa5w,ts,per1,
c     1      pid,ipid,ice,' MEAN',vs30)
c       endif
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
     1' having > 3% contribution)',/,
     2 'Source Category:',t33,' % contr.  R(km)    M   epsilon0 (mean values).')
      do k=0,nsrcmx
      pc=100.*ratecat(k,iax)/ratett(0)*arate/ratett(0)
      if(pc.ge.3.)then
      write(1,777)srcdesc(k),pc,rsav(k,iax)/ratecat(k,iax),amag(k,iax)/ratecat(k,iax),
     1      esav(k,iax)/ratecat(k,iax)
      endif
      enddo
	if(crustal(iax))then
      write(1,788)
      do k=1,lstft
c Recompute rate on fault using most recently found rate for sources. Rate can be zero
c if site is mixing CEUS and WUS faults because diff GMMs are associated with these.
	if(ratef(k,iax).gt.0.)then
         ratef(k,iax)=ratef(k,iax)*arate/ratett(0)      
         write(1,777)fltnm(k),ratef(k,iax),rbf(k,iax),embf(k,iax),ebf(k,iax),azi(k)
        endif	!rate>0?
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
c      aratemn=0.05
	aratemn=0.1	!changed dec 2 2013.
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
c Rmax is the X-axis maximum distance (units km)
c Emmax is the Y-axis maximum magnitude (moment mag)
      parameter (deltam=0.2, amagmin=4.5, small=0.05)
      parameter (ndmx=100, nfmx=180,nmax=10000,nmmax=50,nsig=6)
c Now Jun 23, 2000: using loc(1:iloc) for locations of all files, newmaps no
c longer used.
      logical epson,maybe
      real, dimension(0:43):: mbarr,mmode,rbarr,rmode,pcmax,ratett,ebar0,emodeb0
      real m3d,m3dpl
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
      common/e3d/r3dpl,m3dpl,e3dmsg
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
c M5.0 isnt quite low enough in 2002. try 4.9?
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
c header
      write(1,17)
      write(1,18)
17    format('#!/bin/csh')
18    format('#plotting script')
      write(1,*)' '
c 2008: change  "Source-site" to Closest distance in graph labels      
      write(1,12)rmax,emmin,emmax,pcmax(iax),xsc,ysc,zsc,xlab,
     1 dya,xw,yw,pidpr(1:ipidr)
12      format('psxyz -R0/',e10.5,'/',f4.2,'/',f4.2,'/0/',e10.5,
     1 '  -Jx',f5.3,'/',f5.3,' -Jz',f5.3,' -B',a6,
     2':"Closest Distance, Rcd (km)":/0.5:"MAGNITUDE (Mw)":/f1a',a3,
     3':"% Contribution to Hazard":wEnSZ -So',f4.2,'/',f5.3,
     4' -CPL/tone.cpt -L -K  -W2 -E136/27', 
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
40      format('psxyz -R -Jx -Jz -So',f4.2,'/',f5.3,' -CPL/tone.cpt', 
     1      ' -L  -K -O -W2 -E136/27  << END >> ',
     2      a,'.ps')
444      format('psxyz -R -Jx -Jz -So',f4.2,'/',f5.3,' -G200 ', 
     1      '-L  -K -O -W2 -E136/27  << END >> ',
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
     1 rmode(iax),mmode(iax),emodeb0(iax),r3dpl,m3dpl,e3dmsg,deltar,unsh,pidpr(1:ipidr)
      elseif(epson)then
       write(1,880)shift,pidpr(1:ipidr),sitecl(iclass),
     + labrs,msg(1:icit),abs(xlong),ylat,
     1 psa(1:j1),ratett(iax),retrnt,rbarr(iax),mbarr(iax),ebar0(iax),
     1 rmode(iax),mmode(iax),emodeb0(iax),r3dpl,m3dpl,e3dmsg,deltar,unsh,pidpr(1:ipidr)
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
48      format('pstext -X',a4,' -Jx1/1   ',
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
     5'pstext  ',
     5'PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O >> ',a,'.ps')
88      format('pstext -X',a4,' -Jx1/1   ',
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
     5'pstext  ',
     5'PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O >> ',a,'.ps')
480      format('pstext -X',a4,' -Jx1/1 -R0/11/0/8.5 -K  ',
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
     5'pstext PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O  >> ',a,'.ps')
880      format
     1('pstext -X',a4,' -Jx1/1 -R0/11/0/8.5 -K  ',
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
     5'pstext PL/epsp.txt -X',a4,' -Jx1/1 -R -K -O  >> ',a,'.ps')
      write(1,50)labrs,vs30,comm,pidpr(1:ipidr)
50      format('psxy ','PL/squares_e ',
     2'-U"Distance (R), magnitude (M), epsilon ',
     2'(E0,E) deaggregation for a site on ',a,' with ',
     2'average vs=',f5.0,' m/s top 30 m.', 
     2' NSHMP-2013 PSHA',a14,' Bins with lt 0.1% contrib. omitted"',
     2' -Jx1/1 -R -O -Ss.25 -CPL/tone.cpt  >> ',a,'.ps')
      close(1)
      return
      end  subroutine gmtscript    
      


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
