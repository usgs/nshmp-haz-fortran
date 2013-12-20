c--- program deaggSUBD.2013.f from deaggSUBD.2011.f for web 2013 interactive deagg tool.
c 11/15/2013: put in some anelastic atten. into AB03 for T>1s and R>400km. This
c is a proposal and has not been approved by staff. Furthermore, this change
c has not been imported into hazSUBXnga yet.
c
c 6/25/2013: add CMS computations to NAAsub and atkinsub.
c                These need to be revised because the crustal event correlation
c                not appropriate to subduction events Meera.
c 6/10/2013: add BCHydro (getNAAsub) GMPE Index 20; CMS not ready. When combining,
c                global index for BCHydro is 28.
c 6/10/2013: add Atkinson Macias GMPE. Index 21. CMS not ready
c 6/06/2011: add a vector to CMS output that specifies whether SA at each period T can
c be computed for the given GMPE. The code "combine_cms.f" was also modified to
c examine this vector. 
c 4/08/2011: minor changes to  getGeom 
c 12/27/2010-1/04/2011: add conditional mean spectra determination. 
c 3/01/2011: corrected ss in zhao subroutine. Corrected ss raises interface medians.
c Periods where several GMPEs defined are too restricted so addnl coeffs have been interpolated.
c   For meanspec, There is ambiguity between 0.01 and 0.00 s SA. 
c In this code: Define SA(0.01s)=PGA to keep things simple.
c NEW OUTPUT FILES *.SPCTRA (binned cms) and *.MEANSPC (mean of cms, one for each attn model
c  if using the deagg_attn_model option, or just one meanspec if not.
c June, 2010 version: add geographic deagg as an option. Dont want a separate program just
c to do the geographic deagg calcs. Goal: Fewer progs to debug and maintain.
c Modified mar 16 2010: getGeom has the petersen-hartzell corr. if iatten=-2
c revised Mar 23 2010: getGeom and getAB03 have the harmsen corr. if iatten=-2 or -5, resp.
c modified july 14 2009: getGeom now has nonlinear siteamp based on the ba-nga formulation,
c idea from Art Frankel. The other subduction GMPE routines are not modified. 
c Apr 30 2009: output file can include deaggs associated with each attenuation
c model separately tabulated. The mean is still output
c and its data are associated with index 0. Positive indexes corrrespond to specific GMPEs. A run can
c have up to 3 GMPEs. Thus, the output file size is potentially 3 times larger than that of 
c  the earlier code deaggSUBD.2008.f
c
c Modified Mar 2009. Mar 3 09: Output backazimuth binned by M and irup?
c Originally from hazSUBXv3.f, a Frankel-Wesson code.
c This code deaggregates subduction hazard at one spectral period (or pga),
c one return time (SA is input corresponding to that return) and at one specific site
c    whose Vs30 is specified (NEHRP site classes are used for site response). 1/3 of the
c    2008 deagg program set called by cshell script "run_2008" on gldplone.
c
c  Technical Contact: Steve Harmsen USGS. 303 273 8567. harmsen@usgs.gov
c
c  NGA in program name because some "modern" attenuation models of 2006 are present here, even
c  though these are not the product of the famous PEER NGA teams.
c Deaggregates probabilistic hazard by R,M,eps for subduction events. These can be
c characteristic (slab filling) or GR (floating) ruptures. If GR, you specify the a and b values
c where a is the interval rate of M0+-0.05 events and b is slope parameter:
c             rate(M)=10**(a-bM),
c as well as the magnitude interval over which this is valid [Mmin, Mmax].
c The code will consider M from Mmin to Mmax (no shifting of M to interior).
c  If "cascadia" is in the input file name, M8.8 and up are slab filling even though they may
c enter as itype 2 (GR or floating ruptures). This is a 2007 special feature.
c 
c 11/21: attenuation models are read in, because SP set can be different from LP  set.
c 10/20/2008: ABsub modify 2hz and 3.33 hz as well. add 1.33 hz (0.75s) These frequencies' coeffs
c are interpolated from data which include the modified 2.5 and 5hz coefs.
c
c 10/10/08: modify ABsub to get 5hz from a mix of 5hz and 2.5hz (AB BSSA oct 2008)
c 2/13/2008: trim edge of final floating source to stay within defined slab
c 2/01/2008: add more period coeffs to Geomatrix.
c 10/29/2007: deaggregation option available. Just add site Lat Long nper sa(j),j=1,...,nper
c These fixed sa values must be in same order as the period set in input file.
c
c 5/21/2007: add Zhao atten for 11 periods and 3 site classes: B, C, and D.
c 5/22/2007: add Kanno atten for 11 periods andc
c 6/25/2007: change some parameters to real*8 variables. PC issue with erf
c Added Gregor et al. & ABsub subroutine for Cascadia subd, See list of indexes below.
c Sadigh and Geomatrix are now subroutines. 10/18/2006. Sadigh has a 520-to-760 site corr.
c based on BJF97 siteamp model. Affects longer periods the most, 2sc reduction 0.2486 in logspace.
c Steve Harmsen. modified some of the header data and used f95 "type" for compatibility
c with gfortran and ifort (intel) on PC Cluster. This code uses a 1-d p() array for general
c probability of exceedance. It pulls elements from p() as they are needed.
c Limit on slab length: about 3000 km.
cc
c Depth to slab H fixed at 20 km in Geomatrix. Probably should program variable H(?),
c but small variation has very little effect.
c Important mod: THe GR relation does not recenter the lower and upper bound
c magnitudes. They are as stated in input file. For example, if you input
c M8.0 and 8.7 as the bounds, then in the computations Mmin is 8.0 and Mmax is 8.7. 
c This conforming to input is in distinction
c to hazFXnga and other FX codes, where code resets Mmin to 8.05 and Mmax to 8.65 in this example.
c  It is important to keep this distinction in mind for the 2007 Cascadia models
c that we consider.
c  This code can have up to 12 fixed-M megathrust mags (in 2007, there are 3:8.8, 9, 9.2) described
c    on separate lines.
c 
c----Frankel comments:
c --- added distance taper of weighting
c---- Sadigh et al. for M9+-; set it to M8.5
c--- fixed M9 rupture to not float
c--- modified to specify fault length for subduction zone EQ
c--- also increased various array dimensions
c---  hazSUBX.f
c--- calculates mean annual frequency of exceedance for dipping planar subduction
c--- Most of geometry-related code written by R. Wesson, USGS. Original code by A Frankel,
c --- USGS. Deagg additions and several other revisions by Steve Harmsen, USGS.
c--- this version handles one spectral period and several gnd motion relations. The
c --- GMPEs are listed in the input file. For long period analysis, a different set
c --- of GMPEs is used than for shorter-period analysis. The corresponding input files
c --- have extension "LP" to indicate that they should be used in that instance.
c--- to compile on SUN Solaris system: cc -c iosubs.c 
c                 f95 -o deaggSUBD.2013 deaggSUBD.2013.f  -fast -e -ftrap=%none
c
c  The last ftrap flag seems to be necessary for many runs to avoid fatal errors. 
c
c --- To compile with gnu fortran (often available on PCs with Windows OS):
c          gfortran deaggSUBD.2013.f -o deaggSUBD.2013.exe -O -ffixed-line-length-none -ffpe-trap= -static
c Use the latest available gfortran e.g. March 2008 for best results.
c --- To compile with Intel fortran (ifort):
c      ifort deaggSUBD.2013.f  -o deaggSUBD.2013 -O -132  (perhaps other flags)
c
c Note: this code does not use any of the iosubs.o routines. Do  not link to iosubs.o.
c
c Deaggregation example: 
c
c  deaggSUBD.2013 cascadia.bot.8387.LP cabo 44.4 -123.4 5. .06 760. ./ ./ 1
c
c The first arg after pgm name is the input_file name, 2nd is output file; 3rd is lat, 4th is long, 5th (5) spectral_per,
c and 5  SA, Then Vs30 and directory info
c   Deagg output: one file per period, 2 header lines and one to several data lines. 
c Each Data line contains:
c        Rbin Mbin Haz(R,M) EPS<2 EPS<1 EPS<0 EPS<...  EPS0 SALev(g)
c where Rbin is the Rbar in the bin, Mbin is Mbar in the bin, Haz(R,M) is the hazard
c associated with this bin, and EPSILON is further partitioned.
c Deagg bins: dM is 0.1. dR is variable: 10 km inner 100 km, 50 km to 300 km, 
c       then two catchall bins, 300-500 km and  more than 500 km from the site, resp.
c            Epsilon bins are 1 sigma in size with edges at k*sigma, k=0, +-1,+-2
c            Epsilon binning might be considered coarse. Follows 2002 1996 web sites.
c  Command-line info follows the interactive deagg web site precedents of 2002 and 1996.
c-----  your choice of attenuation relations: -, Sadigh,
c Geomatrix, Gregor et al (2002&6), ABsub-2003, Crouse. 
c This set is reduced from earlier code,
c just including subduction-related models in deaggSUBD.2013 (Sadigh grandfathered in).
c--- output of this program is input into hazallXL.v2 to get probabilistic motions
c iatten=1 unused, would correspond to a model that uses Rjb if there was one.
c iatten=2 Geomatrix subduction (Youngs et al, SRL, 1997). 3 states:
c            rock for vs30>700, firmsoil(avg rock&soil) for vs>444, or soil
c iatten= -2. Geomatrix, with Petersen distance decay for pga, 1s and 0.2s SA only
c iatten=3 Sadigh et al. Rock 
c iatten=-3 Sadigh Soil vers. This subroutine is called with fractional wt for Soil C
c    SadighS corrected 4/16/2007. 1hz Median appears too strong compared to others
c iatten=4 Atkinson and Boore, Cascadia, BSSA, August 2003. Has NEHRP site class factors
c iatten=5 Atkinson and Boore, Global, ditto. Has NEHRP site class factors. added Nov11 06 
c            Note: Global>Cascadia median generally. Global is smidge lower for 1s though
c iatten=6 Crouse
c iatten = 7 Zhao et al(BSSA,  June 2006). rock, C soil, and D soil. with mag-cor term.
c            Zhao has intraplate source class and this model is active when islab = 1 in calling list.
c 5/22/2007: add Kanno atten for 11 periods and variable site. Kanno predicts RMS accel, not
c                geometric mean. Be careful about mixing Kanno with other GMPEs
c iatten=15 Gregor et al. (BSSA,2002) rock and soil. Do not use. Replace with 16.
c iatten=16 Gregor et al. (SRL, 2006) w/variable siteamp. 
c    Ground motions are damped response spectra (5%) and PGA. PGV, PGD n/a.
c--- ground motion levels should be in units of g
c--- period=0 indicates PGA
      parameter (npmx=1)
      parameter ( pi=3.14159265,pi2=pi/2.,fourpisq=39.4784175,tiny=1.e-16)
      dimension dmin2(50,800,6)
      real, dimension(2):: ztop
      real, dimension(2,npmx) :: sumwt
      real, dimension (0:37) :: perka
      integer, dimension(10):: iperg,iperk,ipers
c type replaces structure & requires gfortran or f95 to compile. It is included
c so that this code has compatibility with Linux cluster computer compilers. For one-sta
c Subduction deagg, we dont use iosubs or "type"  variables.
c Also compatible with Windows computers and gfortran compiler.
c      type header_type
c        character*30 :: name(6)
c        real*4 :: period
c        integer*4 :: nlev
c        real*4 :: xlev(20)
c        real*4 :: extra(10)
      common/ipindx/iperba,ipgeom,lvs
        common/sdi/sdi,dy_sdi,fac_sde
      logical isbig,isclose,p_corr,pcorab,mendo,sdi/.false./
      logical l_ms,lvs,isok      !lvs=true if nonlinear site response for getGeom.
      logical verbose, on_top_of_rup, dea_attn,lmsp(0:8,21),ltmp(21)
c dea_attn is true if separate deaggs for each atten model are requested.
      real, dimension(500):: magmin,magmax,ar
c  ar=aspect ratio for floating ruptures.
      real, dimension(2):: slat,slong
      real st(21)        !store components of CMS
      character*4 sname(1)
      character*8 msg(2)
c storage for geographic deagg.
      real, dimension (24,2):: rstar,mstar,hstar
      real, dimension (24,2,2):: azbar
      real*8 dp,pr0,emax/3.0/, sqrt2/1.4142135623/
      real vgeo(3),vs30,pr00,prl,prlr,prr,rtot,eps
      integer readn,nrec,i,ival(8),ivs,immax/1/,npd,nper
c ivs indicator variable for NEHRP class, using Zhao subroutine, added May 21 2007
      character*80 nameout,name,adum,bindir,outdir
      character*12, dimension(-5:16) :: att_typ   
      common/mech/ss,rev,normal,obl
      common/ms/ia,npnga,jms,l_ms,lmsp,imsp,meanspec,sigspec
      dimension p(25010),perg(25)
c new v30 array: can have different vs30 at each grid point. Or fixed vs30       
      logical ss,rev,normal,obl,pgacalc,geog/.false./
c add geog deagg capability june 2010.
c icode() is  used in this code. controls delCi branch in NAAsub
      integer, dimension(5) ::  m_att,iatten,icode
        real, dimension(21) :: ms_p
      real, dimension(16):: pdAM        !Atkinson Macias
      real perba(23),perx(8),perz(26),safix
      integer  nattn
c npnga = number of periods available in all NGA-W GMPEs, 21 or more.
c Goal to define mean spectrum at 21 periods. This is not even remotely possible with current 
c interface GMPE T-set.
c
        integer imsp(8,21),npnga/21/,kmsp(0:8),ipzh
c kmsp(ia) is the index of the maximum period available for att model ia.
       real, dimension(12):: a,b,dmag,cmag,crate,wtm
c new variable wtm stores the branch weight for the a,b,dmag or cmag,crate model.
      real, dimension (17,33,10,21,0:9) :: meanspecs        !save avg meanspec in bins.
      real, dimension (21,0:9):: mean2spec        !mean of the mean specs
c mean2haz, etc mean of the means, averaged over the r- m- and epsilon- bins, one for each
c attn model and one for the avg of all attn models (in slot 0)
      real, dimension (0:9):: mean2haz,mean2r,mean2m,mean2e0
c      real, dimension (17,33,10,0:9) :: hazmsp,epsmsp,rmsp,mmsp
c meanspec and sigmaspec: for each source. accumulate in meanspecs. 
      real, dimension(21) :: meanspec,sigspec,rho_jb        !for each request. Fills inside each GMPE subroutine.
       integer, dimension(12):: nmag0
       integer itype,iftype,iperab(19)
       integer, dimension(10):: delCi
       real, dimension(22):: NAAper
      dimension npts(5)
      real, dimension (19):: perab
      real,dimension(6):: ptail
      real, dimension(20):: pergeo
       dimension rate(50,12),ratem(50)
      dimension nrup(50,12)
      dimension tarray(2)
      dimension wtdist(6),wt(6,2)
      real, dimension (16,16,10,0:4):: rbar,mbar,ebar,haz
      real prob5(16,16,10,0:4,5),magminn
      real period,piriod,xlev     !piriod=gregor approximation march 14 2007
       logical v30a,cascadia,abb03(7,7), rock,soil,go_gail/.false./,go_norm/.false./
c initial set of predefined periods T at which mean spectrum is to be determined. Oct 28 2010.
c Note, ms_p(1)=0.00, but on output ms_p(1) is 0.01 s. This is because most relations do not distinguish
c these.
       ztop = (/20.,50./)
       ms_p = (/0.00, 0.020, 0.030, 0.050, 0.075, 0.100,
     + 0.150, 0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000,
     + 1.500, 2.000, 3.000, 4.000, 5.0, 7.5, 10.0/)
c The longer periods are the ones where response is going to be extrapolated...Or what?
c  Seismic Hazard deagg, july 2009... add a few pergeo T to round out mean spectrum (linear interp in logT)
        vgeo=(/760.,300.,475./)
         pergeo= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,0.075,0.4,0.75,1.5,3.,4.,5.,
     + 0.02,0.03,0.05,0.15,0.25,3.33/)      
       perab = (/0.,0.04,0.2,1.0,0.1,0.3,0.4,0.5,0.75,2.0,3.0,0.02,0.03,0.05,0.075,0.15,
     +0.25,1.50,3.33/)
c NAAper = period set for Norm's BC Hydro subduction model. 0.00 for 0 s (PGA) to 0.02 s SA
        NAAper =(/0.00, 0.050, 0.075,0.100, 0.150,0.200,0.250,0.300,0.400,
     + 0.500, 0.600, 0.750, 1.000, 1.500, 2.000, 2.500, 3.000, 4.000, 5.000, 6.0, 7.500, 10.0/)
c PERba Boore Atkinson NGA. from ba_02apr07_usnr.xls coef file
      perba= (/-1.000, 0.000, 0.010, 0.020, 0.030, 0.050, 0.075, 0.100,
     + 0.150, 0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000,
     + 1.500, 2.000, 3.000, 4.000, 5.000, 7.500,10.0/)
       perx = (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
c gregor period set, perg
c pdAM = period set for Atkinson-Macias subduction GMPE(BSSA, 2009)
      pdAM = (/0.0,0.05,0.1,0.2,0.3,0.4,0.5,1.,2.,2.5,3.,4.,5.,7.7,10.,-1./)
c the 2006 version of perg below. 0.0 replaces 0.01-s for consistency with others
      perg = (/ 0.0, 0.02, 0.025, 0.03, 0.04, 0.05, 0.056, 0.063, 0.071428571,
     1               0.083, 0.1, 0.125, 0.143, 0.167, 0.2, 0.25, 0.333, 0.4, 
     2               0.5, 0.769, 1.0, 1.667, 2.0, 2.5, 5.0 /)
c perka = period set for Kanno et al., BSSA 2006. 0 = pga.
       perka=(/0.0,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.15,0.17,0.20,
     +0.22,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,
     +1.30,1.50,1.70,2.00,2.20,2.50,3.00,3.50,4.00,4.50,5.00/)
c zhao period set, perz. added 3 more 1/04/2011. SH.
      perz = (/0.0,0.02,0.03,0.075,0.05,0.10,0.15,0.20,0.25,0.30,0.40,
     &     0.50,0.60,0.70,0.75,0.80,0.90,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0,10./)
       ptail = (/0.9772499,.8413447,0.5,.1586553,.022750139,
     + 0.0013499/)
      att_typ =  (/ 'AB03 w/R200C','Space Unoccu','Sadigh Soil ',        !-5,-4,-3
     + 'Geom w/R200C','Space Unoccu','WtMean(GMPE)',                	!-2,-1,0,
     + 'Space Unoccu','Geomatrix Yn','Sadigh rock ',        !1,2,3
     + 'AB03 Cascadi','AB03 WorldCo','Crouse 1997 ',
     + 'Zhao et al. ','Kanno et al.','Space Unoccu',        !7,8,9
     +'Space Unoccu','Space Unoccu','Space Unoccu',
     +'Space Unoccu','Space Unoccu','Gregor et al',        !13,14,15
     +'Gregor et al'/)
c new 2010: p_corr controls whether to use Petersen-Hartzell distance decay effect in
c Geomatrix relation.
        
        p_corr = .false.
        pcorab = .false.        !new Mar 23 2010. Addn'l anelastic atten for AB03 at select fr.
      call getarg(2,nameout)
      call getarg(3,adum)
      read(adum,'(f7.4)')rlatd
      call getarg(4,adum)
      read(adum,'(f9.4)')rlond
c rlatd rlond are the coordinates of the deagg site. These override the stuff
c in the input file.
      call getarg(5,adum)
        npd=1
      read(adum,'(f5.3)')per
      call getarg(6,adum)
      read(adum,'(f9.6)')safix
      call getarg(7,adum)
      read(adum,'(f6.1)')vs30
c are we really using arg8 and 9?
      call getarg(8,bindir)
      call getarg(9,outdir)
        call getarg(10,adum)
        read(adum,'(i1)')i
        dea_attn = i .gt.0
      if(iargc().gt.10)then
      call getarg(11,adum)
        read(adum,'(i1)')i
        geog = i.eq.1
        l_ms = i.eq.2
        msg(1)='floating'
        msg(2)='megathrs'
      rstar=0.; hstar=0.; mstar=0.
      azbar=0.        !azimuth info.
        endif
      isize=index(outdir,' ')-1
      nameout=outdir(1:isize)//nameout
      call getarg(1,name)
c if running on PC, may need the following call, commented out on Solaris
c      call initialize()
      abb03=.false.      !initialize 
      rtot=0.0
      prob=0.0
      sumwt= 0.0      !initialize. used to check atten model weights
      haz=0.
      rbar=0.
      ebar=0.
      mbar=0.
      prob5=0.
      magminn=100.
c      headr%name(1)= name
      cascadia=index(name,'cascadia').gt.0
c      print *,cascadia
c if cascadia, all M>=8.8 fill the subduction surface with slip. AF mod added here.
      dum=dtime(tarray)
c      write(*,*)tarray
      inquire(file=name,exist=isok)
      if(isok)then
      open(unit=1,file=name,status='old')
      write(6,*)'deaggSUBD.2013 (6/10/2013);  input file ',name
      else
      write(*,*)'deaggSUBD.2013: File not found: ',name
      write(6,*)'deaggSUBD.2013: File not found: ',name
      stop 'Put in working dir and retry'
      endif
c Initialize truncated normal array Pex, store in p(). Arg to erf( ) has a
c real*8 flavor as of 6/25/2007.
        prl=0.5*(erf(emax/sqrt2) + 1.0)
        prlr=1.0/prl
        prr=1.0-prl
        plim=-emax/sqrt2
        pr0=plim
        dp=0.00004*(3.3-plim)
        dp2=1.0/dp
        p(1)=1.e-9
        do i=2,25001
        pr0=pr0+dp
        p(i)=((erf(pr0)+1.)*0.5-prr)*prlr
        enddo
        p(25002)=1.0
c End initializing the truncated normal PEx() array=p()
c The indep. variable is a real*8 to reduce discretization error.
        nrec=1
580      format(a,$)      
       i=1
        slat(1)=rlatd; slong(1)=rlond; sname(1)='DEAG'
c      write(6,*)slat(i),slong(i),sname(i)
      if(slat(i).lt.-88.)stop 'invalid station location. No Antarctica'
      nx=1
      ny=1
      v30a = .false.
      d2p5=2.      !not used, but depth of rock with Vs of 2.5 km/s
c      write(6,*)' Vs30 m/s is ',vs30
       nper=1
         period=per
        if(sdi)fac_sde = alog(980.*max(0.01,period)**2/fourpisq)
       isz=index(nameout,' ')-1
      open(21,file=nameout,status='unknown')
      write(21,907)safix,period,name,rlatd,rlond,vs30
        mendo=rlatd.lt.40.5.and.rlond.lt.-100..and.rlond.gt.-130.        
c Patch because backazimuths are known to fail for sites S of 40.5.
        mendo = mendo .and. rlatd.gt.26.5        !do not extend S of USA
c this is patched in an unsatisfactory manner here
 907  format('#deaggSUBD.2013 @fixed GM=',f6.3,' g; T=',f6.3,' s, infi ',a,
     + /,'#site lat long = ',f8.4,1x,f9.4,' Vs30 ',f6.1)      
      if(geog)then
      open(41,file=nameout(1:isz)//'.GEOG',status='unknown')
      write(41,907)safix,period,name,rlatd,rlond,vs30
      elseif(l_ms)then
c      open(31,file=nameout(1:isz)//'.SPCTRA',status='unknown')
c unit 33 is binary, used when combining.
      open(33,file=nameout(1:isz)//'.bin',status='unknown',form='unformatted')
c      write(31,907)safix,period,name,rlatd,rlond,vs30
      write(33)safix,period,rlatd,rlond,vs30
c      open(32,file=nameout(1:isz)//'.MEANSPC',status='unknown')
c      write(32,907)safix,period,name,rlatd,rlond,vs30
c      hazmsp=0.        !initialize if needed
c      rmsp=0.; mmsp=0.        !rbar,mbar for the source spectra in bins.
c      epsmsp=0.
      meanspecs=0.
      mean2spec=0.
      mean2haz=0.
      mean2r=0.
      mean2m=0.
      mean2e0=0.
      kmsp(0)=npnga        !initialize at the max numbr of periods. We know LP not
c available.
c intialize rho_jb from the article in Baker Jayaram, Eq Spectra Mar 2008.
        t1=period
        do i=1,npnga
        t_min=min(t1,ms_p(i))
        t_max=max(t1,ms_p(i))
        c1=(1.-cos(pi2-alog(t_max/max(t_min, 0.109))*0.366))
        if(t_min .eq. t_max)then
        rho_jb(i)=1.0
        jms=i
c        print *,'index in spectrum of input period ',jms
        else
        if(t_max.lt.0.2)
     +   c2=1.0-0.105*(1.-1./(1.+exp(100*t_max-5.)))*(t_max-t_min)/(t_max-0.0099)
        if (t_max .lt.0.109)then
        c3=c2
        else
        c3=c1
        endif
        c4=c1+0.5*(sqrt(c3)-c3)*(1.+cos(pi*t_min/0.109))
        if (t_max .le. 0.109)then
        rho_jb(i) = c2
        elseif(t_min .gt. 0.109)then
        rho_jb(i) = c1
        elseif(t_max .lt. 0.2)then
        rho_jb(i) = min(c2,c4)
        else
        rho_jb(i) = c4
        endif
        endif
c        print *,ms_p(i),rho_jb(i)
        enddo
      endif

c      endif
c         write(6,580) "Number of atten. relations: "
c different attenuation model sets used for SP and LP. New 11/2008.
        read (1,*)nattn
        do k=1,nattn
        read (1,*) iatten(k),wt(k,1)
        wt(k,2)=wt(k,1)
        enddo
         wtdist=2000.
        ip=1
c--- loop through atten relations for that period
         do 20 ia=1,nattn
              sumwt(1,ip)=sumwt(1,ip)+wt(ia,1)
              sumwt(2,ip)=sumwt(2,ip)+wt(ia,2)
c--------  Geomatrix interface atten.
         if(abs(iatten(ia)).eq.2) then
         if(period.le.0.01)then
         ipgeom=1
         else
         j=2
         dowhile(abs(period-pergeo(j)).gt.0.01)
         j=j+1
         if(j.eq.20)stop 'period doesnt match Geomatrix subd. set'
         enddo
         ipgeom=j      !period index for Geomatrix relations.
         endif        !period <=0.01?
        m_att(ia)=13        !global index for Geomatrix is 13
         if(iatten(ia).lt.0)then
         p_corr=.true.
c         write(6,*)'Petersen correction will apply to 1hz 5hz and pga'
         iatten(ia)=2
         endif
c may also need the companion BA period index 
cc       write(6,*)'Period map for Geomatrix ',j
c         headr%name(6)='Geomatrix attn'      !period-indep assignment.
cc       write(6,*)'Geomatrix is subroutine, no coeff. readin Oct 2006'
         j=1
         dowhile(period.ne.perba(j))
         j=j+1
         if(j.eq.24)stop 'period doesnt match BA siteamp period set for geom'
         enddo
         iperba=j
        lvs = abs(vs30-vgeo(1)).gt.3. .and. abs(vs30-vgeo(2)).gt.3.
       if(l_ms)then
c fill imsp and lmsp arrays with can do information.
        kmsp(ia)=19        !5s
        kmsp(0)=min(kmsp(0),kmsp(ia))
        lmsp(ia,1)=.true.
        imsp(ia,1)=1
        do j=2,kmsp(ia)
        ka=2
        dowhile(abs(ms_p(j)-pergeo(ka)).gt.0.01 .and.ka.lt.20)
        ka=ka+1
        enddo
        if (ka.gt.19)then
        lmsp(ia,j)=.false.
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
c        print *,ms_p(j),' geomatrix can do this T'
        endif
        enddo        !j=1,npnga  
        endif        !mean spectrum details 
c modify siteamp iff vs30 is at least 1%  different from all reference velocities.
         elseif(iatten(ia).eq.7) then
         j=1
         dowhile(period.ne.perz(j))
         j=j+1
         if(j.gt.25)stop 'period doesnt match available Zhao set 11-2008'
         enddo
         ipzh=j      !period index for Geomatrix relations.
         m_att(ia)=27        !27 is global Zhao index.
cc       write(6,*)'Period map for Zhao et al.',j
       if(l_ms)then
c fill imsp and lmsp arrays with can do information.
        do j=1,npnga
        ka=1
        dowhile(ms_p(j).ne.perz(ka).and.ka.lt.26)
        ka=ka+1
        enddo
        if (ka.gt.25)then
        lmsp(ia,j)=.false.
c        print *,ms_p(j),' Zhao cannot do this period'
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
c        print *,ms_p(j),' Zhao can do this period'
        endif
        enddo        !j=1,npnga  
        kmsp(ia)=19        !5s
        kmsp(0)=min(kmsp(0),kmsp(ia))
        endif        !mean spectrum details 
c         headr%name(6)='Zhao et al attn'
c        zhao et al, new may 2007.
c------- Sadigh type
       elseif(iatten(ia).eq.20)then
      ka=1
c New BC Hydro Subd GMPE Jan 2011. For spectral period <=0.02 s, use the pga model.
      if(period.le.0.02)then
      piriod=0.0
      else
      piriod=period
      endif
      dowhile(abs(piriod-NAAper(ka)).gt.0.001)
      ka=ka+1
      if(ka.gt.22)then
      write(6,*)'Period ',period
      stop 'For this period please remove the Abrahamson BC hydro relation from input file'
      endif
      enddo
       ipernaa=ka     
             kbackarc=0        !don't use backarc in PacNW runs.
       izpar=1
       m_att(ia)= 28        !now 28 for BChydro
       write(6,*)period,ipernaa,' BC_Hydro ip map'
       go_norm=.true.
c for NAA relation, use epistemic branch stored in "icode". If icode is invalid, use central branch
c This is a modification May 9 2013. SH.
       if(icode(ia).gt.0 .and. icode(ia).lt.4)then
       delCi(ia)=icode(ia)
       else
       delCi(ia)=2        !a default value (central branch)
       endif
       write(6,*)ia,delCi(ia),' ia and  BC_Hydro GM uncert. branch'
       if(l_ms)then
       m_att(ia) = 30        !associate NAASub with global GMPE index 30 in combine_cms.
c fill imsp and lmsp arrays with can do information.
        imsp(ia,1)=1
        lmsp(ia,1)=.true.        !can do pga or 0.01s SA equivalent
        do j=2,npnga
        ka=1
        dowhile(ms_p(j).ne.NAAper(ka).and.ka.lt.22)
        ka=ka+1
        enddo
        if (ka.eq.22.and.ms_p(j).ne.NAAper(ka))then
        lmsp(ia,j)=.false.
c        print *,ms_p(j),' NAA cannot do this period'
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
c        print *,ms_p(j),' NAA (BCHydro) can do this period'
        endif
        enddo        !j=1,npnga  
        kmsp(ia)=21        !10s
        kmsp(0)=min(kmsp(0),kmsp(ia))
        print *,imsp(ia,1:21),'NAA', imsp(ia,4),lmsp(ia,4)
        print *,lmsp(ia,1:21),'NAA logicals' 
        endif        !mean spectrum details 
       elseif(iatten(ia).eq.21)then
      ka=1
c New atkinson-macias (BSSA,2009) GMPE nov 2012. For spectral period <=0.02 s, use the pga model.
      if(period.gt.0.02)then
        ka=2
      dowhile(abs(period-pdAM(ka)).gt.0.001)
      ka=ka+1
      if(ka.gt.15)then
      write(6,*)'Period ',period
      stop 'For this period please remove the Atkinson Macias relation from input file'
      endif
      enddo
      endif
       ipergail=ka     
             kbackarc=0        !don't use backarc in PacNW runs.
         m_att(ia)=29        !use global index 29 for A-M subduction.
       write(6,*)period,ipergail,' Atkinson Macias ip map'
       go_gail=.true.
c may also need the companion BA period index for soil amp function "basiteamp"
         j=1
         dowhile(period.ne.perba(j))
         j=j+1
         if(j.eq.24)stop 'period doesnt match BA siteamp period set for including with AtMac09'
         enddo
         iperba=j
       if(l_ms)then
c fill imsp and lmsp arrays with can do information.
        imsp(ia,1)=1
        lmsp(ia,1)=.true.
        do j=2,npnga
        ka=1
        dowhile(abs(ms_p(j)-pdAM(ka)).gt.0.001.and.ka.lt.16)
        ka=ka+1
        enddo
        if (ka.gt.15)then
        lmsp(ia,j)=.false.
c        print *,ms_p(j),' Atkin&Macias sub cannot do this period'
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
c        print *,ms_p(j),'  Atkinson & Macias sub can do this period'
        endif
        enddo        !j=1,npnga 
        print *,imsp(ia,1:21),'A&M'
        print *,lmsp(ia,1:21),'A&M logicals' 
        kmsp(ia)=21        !10s
        kmsp(0)=min(kmsp(0),kmsp(ia))
        endif        !mean spectrum details 
            elseif(abs(iatten(ia)).eq.3) then
         m_att(ia)=3
            if(iatten(ia).eq.-3)then
c  c       write(6,*)'Sadigh model with soil coeffs.'
            else
c  c       write(6,*)'Sadigh model with rock coeffs.'
            endif
         j=1
         dowhile(period.ne.perx(j).and.j.lt.npmx)
         j=j+1
         enddo
         if(j.gt.7)stop 'period doesnt match standard set'
         ipers(ip)=j      !period index for Sadigh and perhaps other relations.
       if(l_ms)then
c fill imsp and lmsp arrays with can do information.
        do j=1,npnga
        ka=1
        dowhile(ms_p(j).ne.perx(ka).and.ka.lt.8)
        ka=ka+1
        enddo
        if (ka.gt.7)then
        lmsp(ia,j)=.false.
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
        endif
        enddo        !j=1,npnga  
        endif        !mean spectrum details 
cc       write(6,*)'Period map for Sadigh ',j
       elseif(iatten(ia).eq. 8)then
      ka=0
      dowhile(period.ne.perka(ka))
      ka=ka+1
      if(ka.gt.37)then
      write(6,*)'Period ',per
      stop 'For this period please remove the Kanno relation from input file'
      endif
      enddo
       iperk(ip)=ka
       m_att(ia)=28        !kanno global index 48.     
c      write(6,*)ip,iperk(ip),' Kanno 10/06 ip map'
       if(l_ms)then
c fill imsp and lmsp arrays with can do information.
        do j=1,npnga
        ka=1
        dowhile(ms_p(j).ne.perka(ka).and.ka.lt.38)
        ka=ka+1
        enddo
        if (ka.gt.37)then
        lmsp(ia,j)=.false.
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
        endif
        enddo        !j=1,npnga  
        endif        !mean spectrum details 
c vsfac used for site-amp in  the Kanno et al. relation.
        vsfac=alog10(vs30)          
         elseif(iatten(ia).eq.15)then
c Gregor 2002 added oct 2006. Gregor Email says this has been updated.
cc       write(6,*)'Gregor 2002 subduction model out of date'
         stop 'please use atten model with index 16 for Gregor.'
         elseif(iatten(ia).eq.16)then         
         j=1
c gregor : treat .333 T as 0.3 to combine with other modelers' output.
c alternate might be to interpolate coeffs (not done as of mar14 20007)
         if(period.gt.0.299.and.period.lt.0.345)then
         piriod=0.333
cc       write(6,*)'Note: Using gregor 0.333 T for 0.3 s analysis'
         else
         piriod=period
         endif
         dowhile(abs(piriod-perg(j)).gt.0.01.and.j.lt.26)
         j=j+1
         enddo
c Gregor period set needs to be interpolated to the nga T values.
         if(j.gt.25)stop 'period doesnt match standard Gregor set'
         iperg(ip)=j      !period index for Gregor subduction
         m_att(ia)=29
c updated Gregor logic Jan 2007. Steve Harmsen         
       if(l_ms)then
c fill imsp and lmsp arrays with can do information.
        do j=1,npnga
        ka=1
        dowhile(ms_p(j).ne.perg(ka).and.ka.lt.26)
        ka=ka+1
        enddo
        if (ka.gt.25)then
        lmsp(ia,j)=.false.
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
        endif
        enddo        !j=1,npnga  
        kmsp(ia)=19        !5s
        kmsp(0)=min(kmsp(0),kmsp(ia))
        endif        !mean spectrum details 
c------- AB03 type
            elseif(iatten(ia).eq.4.or.abs(iatten(ia)).eq.5) then
cc       write(6,*)'Atkinson Boore 2003 subduction model'
        abb03(ip,ia)=.true.
         j=1
         dowhile(period.ne.perab(j).and.j.le.18)
         j=j+1
         enddo
         if(j.gt.18)stop 'period doesnt match standard AB subduction set'
         iperab(ip)=j      !period index for Sadigh and perhaps other relations.
         if(iatten(ia).eq.4)then
cc       write(6,*)'Using ABsub with Cascadia coef, period map is ',j
c      headr%name(6)='ABsub Cascadia'
c iclass=1 for PacNW, iclass=3 for worldwide...
        iclass=1
        m_att(ia)=12
         else
cc       write(6,*)'Using ABsub with global coef, period map is ',j
         iclass=3
        m_att(ia)=18
      if(iatten(ia).lt.0)then
      pcorab=.true.
      iatten(ia)=5
c      write(6,*)'AB03: Including the Harmsen correction to R-decay when R>200km'
      endif
         endif
       if(l_ms)then
c fill imsp and lmsp arrays with can do information. Period(10)= 3s. Max available for AB03.
        do j=1,npnga
        ka=1
        dowhile(ms_p(j).ne.perab(ka).and.ka.le.18)
        ka=ka+1
        enddo
        if (ka.gt.18)then
        lmsp(ia,j)=.false.
        else
        imsp(ia,j)=ka
        lmsp(ia,j)=.true.
c        print *,'Can do AB sub @ T= ',ms_p(j),'s? ',lmsp(ia,j)
        endif
        enddo        !j=1,npnga  
        kmsp(ia)=17        !3s
        kmsp(0)=min(kmsp(0),kmsp(ia))
        endif        !mean spectrum details 

c------ Crouse attenuation form
            elseif(iatten(ia).eq.6) then
        m_att(ia)=26
            stop 'Crouse model not supported at this time'
c   c  c       write(6,*) "enter p1,p2,p3,p4,p5,p6,p7,sigma,depth"
c               read(1,*) p1(ip),p2(ip),p3(ip),p4(ip),p5(ip),p6(ip), p7
c     +         (ip),psigma(ip),ph(ip)
            else
c  c       write(6,*)' Surprise atten index. ',iatten(ia)
c  c       write(6,*)' Code does not know how to handle this one;'
c  c       write(6,*)'Please review your input file. Valid -3, 2 to 8 or 16'
            stop 'deaggSUBD.2013.f: did not finish. Check log file.'
            endif
c----- Boore look-up table: removed 10-2006. SH
   20    continue
         nlev=1
        ip=1
   30    xlev= alog(safix)
   40 continue
c Add june 6 2011: what spectral periods can be computed for mean spectrum?
c Ans: those that are in common to all invoked GMPEs. Below logic puts this
c information in lmsp(0,*), which is output later on.
        if(l_ms)then
        lmsp(0,1:21)=lmsp(1,1:21)
        do ia=2,nattn
        lmsp(0,1:21)=lmsp(0,1:21) .and. lmsp(ia,1:21)
        enddo        ! Potential cms periods
        endif        !compute cms?
cccccccccccccccccccccccccccccccccccc      endif
c      write(6,580) "Increment dlen (km) for source segment "
      read(1,*) dlen
c      write(6,*)dlen
c      write(6,580) "Distance increment, dmax (km): "
      read(1,*) di,dmax
c      write(6,*) di,dmax
ccc-------read fault data. there is only one.
         ift=1
c         write(6,580)
c     +   "Enter 1 for char. 2 for G-R; 1 for SS, 2 for reverse: "
          read(1,*,end=60) itype,iftype,nmodel
        rev=iftype.eq.2
         ss=iftype.eq.1
         normal=iftype.eq.3
         if(normal)stop' Attenuation models not prepared for normal slip'
         do 50 imod=1,nmodel
         if(itype.eq.1) then
         write(6,580) "Char. M and rate and weight: "
4442      FORMAT('Rx    Ry  Mchar  ',f6.2,' Rate_Mchar ',e10.5,' filein ',a20)
            read(1,*) cmag(imod),crate(imod),wtm(imod)
          write(6,*)cmag(imod),crate(imod),wtm(imod)
          magmax(imod)=cmag(imod); magmin(imod)=cmag(imod); dmag(imod)=0.1
          nmag0(imod)=1
             magminn=min(magminn,magmin(imod))
         elseif(itype.eq.2) then
c            write(6,580) "enter a,b,min M,max M,dmag,wt for freq mag curve"
            read (1,*) a(imod),b(imod),magmin(imod),magmax(imod),dmag(imod),wtm(imod)
c  c       write(6,*) a(imod),b(imod),magmin(imod),magmax(imod),dmag(imod),wtm(imod)
c DO NOT Reset the GR Magmin & Magmax as in fltrate, mod of Oct 27 2006. SH.
c Unmodded feb 28 after discussion with Art. Be sure incoming rate is appropriate.
c
c            if(magmax(imod).ge.magmin(imod)+dmag(imod))then
c            magmax(imod)=magmax(imod)-0.5*dmag(imod)
cc            magmin(imod)=magmin(imod)+0.5*dmag(imod)
c            endif
            nmag0(imod)= nint((magmax(imod)-magmin(imod))/dmag(imod)) + 1
            if(nmodel.gt.1.and.magmin(imod).ne.magmax(imod))then
          write(6,*)magmin(imod),magmax(imod),' must be equal for this test code.'
             stop ' rupturesetup will need work for the case of unequal mag'
             endif
             magminn=min(magminn,magmin(imod))
4444      FORMAT('Rx    Ry  Total_Rate  Rate_M(j) from ',f6.2,' to ',f6.2,' n ',i2,
     +' filein ',a20)
            if(magmin(imod).lt.6.5)stop'deaggSUBD.2013: magmin below 6.5'
            else
            stop 'Valid itypes = 1 or 2'
         endif        !itype=2
50     continue
c******************************************************************
c      Replace code to read and setup model of subduction zone
c******************************************************************
c     Insert new code to read and setup model
c***************************************************************
         call readmodel(1,tlen) 
c****************************************************************
c     End new code to read and setup model
c****************************************************************
cc       write(6,*) tlen, ' total length (km)'
         coef= pi/180.
   60 nft= 1
c      write(6,*) "nft=",nft
ccccccccccccccccccccccccccccc
c---- write header --------------
c Discussion at Software meeting Dec 16 2005: header record needs to be
c lengthened to include more input-file information.
      ndist= dmax/di +1.5

c no longer building a table of p(). Just get the mean and sd when you need 'em
c mod oct 19 2006. SH
ccccccccccccccccccccccccccccccccccccccccccc

c      dum=dtime(tarray)
c      write(*,*) tarray


ccccccccccccccccccccccccccccccccccc
      icnt=1
       prob=1.e-21      !initialize array
ccccccccccccccccccccccccccccccccccc
c---Assume only one fault
      ift =1
c---Set up parameters of rupture zones
      do 285 imod=1,nmodel
      do 185 m=1,nmag0(imod)
          xmag = magmin(imod)+(m-1)*dmag(imod)
c          print *,xmag,imod,' xmag,imod'
c------find rupture length from Geomatrix interface eqs relation
          ruplen= (xmag-4.94)/1.39
          ruplen= 10**ruplen
          if((xmag.ge.8.8.and.cascadia).or.itype.eq.1) then
c above mod feb 28 to be consistent with AF models for 2007 PSHA update.
                nrup(m,imod)=1
                ruplen= tlen
         write(6,*) xmag, ruplen, 'Characteristic rupture M & ruplen'
c                read(5,*) idum
          else
          nrup(m,imod)= max(1,nint((tlen-ruplen)/dlen) +1)
c max of 1 and the other expression because if tlen < ruplen, will get 0 or less than 0         
c      write(6,*) "nrup(m,imod),ruplen,tlen",nrup(m,imod),ruplen,tlen
      if(nrup(m,imod).gt.799)stop' nrup(m,imod) exceeds 799. Up dim 2 of rupture()'
          if(ruplen.ge.tlen) then
            nrup(m,imod)=1
                ruplen=tlen
          endif
          endif
          if(itype.eq.1)then
          rate(m,imod)=crate(imod)*wtm(imod)
          else
           rate(m,imod)= a(imod) - b(imod)*xmag
c apply epistemic model weight here corresponding to branch imod.
           rate(m,imod)= wtm(imod)*10.**rate(m,imod)
           endif
       rtot=rtot+rate(m,imod)
           rate(m,imod)= rate(m,imod)/float(nrup(m,imod))
c       write(6,*)'Rate of a given scenario at this M ',rate(m,imod)
           if(imod.eq.1)then
c rupturesetup is not designed for multiple "imods" but if there is just one rupture it works
c for all imods (the Cascadia megathrust case). Code needs work for the GR case
          do 182 irup=1,nrup(m,imod)
               xmove= dlen*(irup-1)
               s1 = xmove
               s2 = xmove + ruplen
c known problem if s2 extends even 1 km beyond end of fault with length tlen. next line protects.
               s2=min(s2,tlen-0.1)      !new feb 13 2008
               call rupturesetup(m,irup,s1,s2)
 182      continue
c      write(6,*)'For this mag the final distance s2 is ',s2,' km'
       endif        !if imod = 1
 185   continue
 285        continue
c
c call rupture setup for a fictitious event that covers the whole fault surface
c store this in m= 24 bin.
      irup=1
      call rupturesetup(24,irup,0.,tlen)
c      print *,' megathrust setup ok '
c      endif
      write(6,*)' Cumulative annual rate of subd events: ',rtot
c      write(6,*)' This is equivalent to recurr. interval of ',1./rtot,' yrs'
c---Here's the guts
c
c---One receiver for deagg
      i = 0
      nloop = 1
 
c     loop through floating rupture zones
c     find minimum distance dmin2 to receiver for each rupture zone
c     nrup(m,imod) rupture zones per magnitude bin m
c     move rupture zones by dlen horizontally along fault
 
c*****calculate distances  from all sources*********
c  One station for deagg analysis.        
                  icnt = 1
c                       i is the index of the receiver
                  i = 1
        rx=slong(i)
        ry=slat(i)
c       write(6,*)'Computing seismic hazard for ',sname(i)
        vsfac=alog10(vs30)      !for kanno
        if(vs30 > 599.)then
        rock = .true.
        soil = .false.      !this is used for getSadigh for example
        ivs=1            !ivs is used in getzhao for site response variability. 1=bc rock
        else
        soil = .true.
        rock = .false.
        if(vs30 > 299.)then
        ivs=2      !zhao soilC
        else
        ivs=3      !zhao soilD
        endif
        endif
          taper= 15.
       taperfac=1./30.  
c make a determination if site is more than Rmax from fault. If so do not
c bother with any more calcs in the do 300 below.
              call qkdist(24,1,rx,ry,d_min,fltdepmin, d_minxy,bazminxy)
c              print *,d_min,' km=min dist to site'
c  d_min is minimum distance from fault to receiver
c d_minxy is the minimum distance for megathrust that ruptures the entire slab 
c ibaz definition same as in deaggSUBv3g.f. Will it work?
       if(d_min.gt.dmax)goto 305   
c Only report the rupture that is closest to site, whose distance should about equal
c distance to megathrust. That is the hope anyway.
        do 300 imod=1,nmodel
         dm_inv=1./dmag(imod)
          do 300 m=1,nmag0(imod)
          xmag = magmin(imod)+(m-1)*dmag(imod)
c dM bins have width dm_inv here 
                  im=nint((xmag-magminn)*dm_inv)+1
                  im=min(16,im)      !current dimension 16
                  immax=max(immax,im)
                  if(xmag.gt.8.7)then        !coarse mag bin for geog deagg.
                  im2=2
                  else
                  im2=1
                  endif        
              do 290 irup=1,nrup(m,imod) 
c-----------------------------------------------------------------------
c  dmin is minimum distance from fault to receiver
               if(imod.eq.1)then
                  call qkdist(m,irup,rx,ry,distmin,fltdepmin, distminxy,bazminxy)
c compute index for distance in deagg matrices
                if(distmin.le.100.)then
                ir=0.1*distmin + 1
                elseif(distmin.le.300.)then
                ir=11+(distmin-100.)/50.
                elseif(distmin.le.500.)then
                ir=15
                else
                ir=16
                endif
c ibaz  defn. taken from 2002 version of deagg, deaggSUBv3g.f. Will it still work now?
                ibaz=4+int(4.*bazminxy/3.1415926)
                    if(distminxy.lt.40.)then
                              inow=ibaz
                           elseif(distminxy.lt.100.)then
                              inow=ibaz+8
                           else
                              inow=ibaz+16
                    endif
                endif         !if imod=1
c                print *,bazminxy,distmin,inow,irup
                  dmin2(m,irup,1) = distminxy
c                  if(ry.eq.49.and.m.eq.1)
c     + write(*,*) "m,irup,rx,ry,distmin",m,irup,rx,ry,distmin,distminxy
c         dmin2(m,irup,2) = fltdepmin
c                  dmin2(m,irup,2) = sqrt(distminxy *distminxy +
c     +            400.)
                  dmin2(m,irup,2) = distmin
                  dmin2(m,irup,3) = distmin
                  dmin2(m,irup,4) = distmin
c By delaying calculation of p() at receiver R, we can use variable vs30 if requested.
c  Mod Oct 2006 SH 
              ip=1
             do  ia=1,nattn
c always use R_cd or Rrup as the distance measure. No subduct. relations use R_jb
                           dist= dmin2(m,irup,3)
                           rcd=distmin
                              ii= dist/di +1.5
                              if(ii.gt.ndist) go to 270
         weight= wt(ia,1)
         if(dist.gt.wtdist(ia)) weight=wt(ia,2)
         wt2= wt(ia,2)
         wt1= wt(ia,1)
         dist0= wtdist(ia)
         if((dist.ge.dist0-taper).and.(dist.le.dist0+taper)) 
     &   weight= wt1 + (wt2-wt1)*(dist-dist0+taper)*taperfac
        
       if(weight.lt.0.0001)goto 270      !some relations downweighed with distance
       meanspec=0.0
      if(abb03(ip,ia))then
      call getABsub(iperab(ip),iclass,0,vs30,xmag,dist,gnd,sigmaf,pcorab)
      elseif(iatten(ia).eq.2) then
c the 2nd arg, 0, refers to subduction or interface. 0 for subduction. Added oct 26
c
            call getGeom(ip,0,vs30,xmag,rcd,gnd,sigmaf,p_corr)
      elseif(iatten(ia).eq.7) then
c the 2nd arg, 0, refers to subduction or interface. 0 for subduction. Added oct 26
c
            hslab = 50.      !not used for subduction. hi is used (param)
             call zhao(ipzh,xmag,rcd,gnd,sigmaf,ivs, 0,hslab)
c Zhao: the second-to-last argument, 0, indicates compute for subduction, not for in-slab source.
      elseif(iatten(ia).eq.8) then
      call kanno(iperk(ip),gnd,sigmaf,xmag,rcd,vsfac)
c--- following are Sadigh forms. Might allow a near-source sigma someday.
         elseif(iatten(ia).eq.3)then
               call getSadighR(ipers(ip),xmag,rcd,gnd,sigmaf,0.,0.)
         elseif(iatten(ia).eq.-3)then
               call getSadighS(ipers(ip),xmag,rcd,gnd,sigmaf,0.,0.)
       elseif(iatten(ia).eq.16)then
             call Gregor06(iperg(ip),xmag,rcd,vs30,gnd,sigmaf)
c gregor 2006 only works for subduction, is therefore outputting "gnd2"
c add BCHydro Nov 2012. from hazSUBXnga.f code.
      elseif(iatten(ia).eq.20)then
c Norm Abrahamson BC Hydro program  
      kevnt=0;
      zH = ztop(izpar)
c      print *,zH,' going in to getNAAsub hypo depth (km)'
      
      call getNAAsub(ipernaa, kevnt, kbackarc, xmag, rcd, zH, sigmaf, gnd, delCi(ia), vs30)
c kevnt=0 interface, not for in-slab source. kbackarc=1 is back-arc. This setting needs
c to be geometry dependent in final version.
        elseif(iatten(ia).eq.21)then
        call getAtkinsub(ipergail,rcd,xmag,vs30,kbackarc,gnd,sigmaf)
c        print *,exp(gnd),sigmaf,'getatkinsub'
c Atkinson Macias without backarc effect (workshop suggestion 12-15-2012)
         endif
      wtrate = weight*rate(m,imod)
          k=1
       pr=(gnd - xlev)*sigmaf
       if(pr.gt.3.3)then
       ipr=25002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 270      !transfer out if ground motion above mu+3sigma
       endif
       cfac=wtrate*p(ipr)
               prob= prob+cfac
c k is supposed to equal 1 if deagg. So we dont have a k index below, rbar etc
c you could deaggregate at each input spectral sa level. then you would need a k
c index on rbar, mbar, ebar, haz, prob5.
               eps= -pr*sqrt2
c               print *,eps,gnd,xlev,cfac,meanspec(jms),xmag
            if(eps.lt.emax)then
            ceps = cfac*eps
            cma = cfac*xmag
            cfr=cfac*rcd
            ieps=int(max(1.,min((eps+2.)*2.,10.)))
            rbar(ir,im,ieps,0)=rbar(ir,im,ieps,0)+cfr
            mbar(ir,im,ieps,0)=mbar(ir,im,ieps,0)+cma
            ebar(ir,im,ieps,0)=ebar(ir,im,ieps,0)+ceps
            haz(ir,im,ieps,0)=haz(ir,im,ieps,0)+cfac
c add geog option june 2010 SHarmsen.
            if(geog)then
             rstar(inow,im2)=rstar(inow,im2) + cfr
             mstar(inow,im2)=mstar(inow,im2) + cma
             hstar(inow,im2)=hstar(inow,im2) + cfac
              azbar(inow,im2,1)=azbar(inow,im2,1)+100.*cfac*cos(bazminxy)
              azbar(inow,im2,2)=azbar(inow,im2,2)+100.*cfac*sin(bazminxy)
             endif
          if(l_ms)then
c Add epsilon*rho*sigma. This is for a specific attenuation model, with index ia.
        meanspec=meanspec + eps*rho_jb*sigspec
c        print *,'Mean SA & S.D. for input period ',exp(meanspec(jms)),sigspec(jms)
c     + ,sigmaf, rho_jb(jms),xmag,eps
          do jp=1,npnga
          if(lmsp(ia,jp))then
          meanspecs(ir,im,ieps,jp,0)= meanspecs(ir,im,ieps,jp,0)+ meanspec(jp)*cfac
          mean2spec(jp,0)= mean2spec(jp,0)+ meanspec(jp)*cfac
          endif
          enddo
c mean2spec(jp,0)=big average spectrum and its big average distance, M, e0 over all variables.
          mean2haz(0)=mean2haz(0)+cfac
          mean2r(0)=mean2r(0)+cfr
          mean2m(0)=mean2m(0)+cma
          mean2e0(0)=mean2e0(0)+ceps
c wttmp all weight factors except probability of exceedance
          endif
            if(dea_attn)then
              rbar(ir,im,ieps,ia)=rbar(ir,im,ieps,ia)+cfr
              mbar(ir,im,ieps,ia)=mbar(ir,im,ieps,ia)+cma
              ebar(ir,im,ieps,ia)=ebar(ir,im,ieps,ia)+ceps
              haz(ir,im,ieps,ia)=haz(ir,im,ieps,ia)+cfac
          if(l_ms)then
          do jp=1,npnga
          if(lmsp(ia,jp))then
          meanspecs(ir,im,ieps,jp,ia)= meanspecs(ir,im,ieps,jp,ia)+ meanspec(jp)*cfac
          mean2spec(jp,ia)= mean2spec(jp,ia)+ meanspec(jp)*cfac
          endif
          enddo
c mean2spec(jp,ia)=big average spectrum and its big average distance, M, e0 for attn model ia.
          mean2haz(ia)=mean2haz(ia)+cfac
          mean2r(ia)=mean2r(ia)+cfr
          mean2m(ia)=mean2m(ia)+cma
          mean2e0(ia)=mean2e0(ia)+ceps
          endif
            endif
            if(eps.ge.2.)goto 266
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
            ka=5
            temp=p(ipr)
            prx= temp-ptail(ka)
            dowhile(prx.gt.1.e-9)
            prob5(ir,im,ieps,0,ka)=prob5(ir,im,ieps,0,ka)+
     +       wtrate*prx      
            if(dea_attn)prob5(ir,im,ieps,ia,ka)=prob5(ir,im,ieps,ia,ka)+
     +       wtrate*prx      
            ka=ka-1
            if(ka.eq.0)goto 266
            prx= temp-ptail(ka)
            enddo
            endif      !eps <emax
  266            continue
  270                continue      
        enddo      !ia index
  280                continue
  290                continue
  300             continue
            ratem=0.
344      format(f7.2,1x,f7.2,6(1x,e10.5))
  305            continue      !we go here if site is > dmax from anyplace on
c the subducting slab interface
  310          continue
c
  320          continue
 
67      format(/,'#Station ',f7.4,1x,f10.4,1x,a,' vs30 ',f6.1)
        ip=1
       do ieps=1,10
       ie = ieps
       do im=1,immax
       do ir=1,16
       prx0=haz(ir,im,ieps,0)
       if(prx0.gt.tiny)then
       e00=max(-9.99,ebar(ir,im,ieps,0)/prx0)
c Future: the 0 below could be  ift, which rupture is contributing the
c mode at that ir, im, ieps location?  
        rcd0=  rbar(ir,im,ieps,0)/prx0  
        xm0=  mbar(ir,im,ieps,0)/prx0
       write(21,25)rcd0,xm0,prx0,(prob5(ir,im,ieps,0,k),k=5,1,-1),
     + 0,e00,0
       if(dea_attn)then
       do ia=1,nattn
              prx=haz(ir,im,ieps,ia)
       if(prx.gt.tiny)then
         e0=max(-9.99,ebar(ir,im,ieps,ia)/prx)
c Future: the 0 below could be  ift, which rupture is contributing the
c mode at that ir, im, ieps location?        
         write(21,25)rbar(ir,im,ieps,ia)/prx,
     +   mbar(ir,im,ieps,ia)/prx,prx,(prob5(ir,im,ieps,ia,k),k=5,1,-1),
     +   0,e0,m_att(ia)
        endif        !inner prx>tiny
        enddo        !IA INDEX FOR INDIV. ATTN MODELS
      endif        !deagg individual attn model contribs? ADD apr 2009.
      endif
25      format(f9.3,1x,f7.3,6(1x,e11.5),
     + 1x,i2,1x,f5.2,1x,i2)
       if(l_ms )then
       ms_p(1)=0.01        !T > 0 for output.
c write  mean spectrum for each binned src.
       if(prx0.gt.tiny)then
c        write(31,255)rcd0,xm0,prx0,e00,att_typ(0)
                  imout=nint((xm0-5.049)*10)+1
       write(33)0,ir,imout,ie,rcd0,xm0,prx0,e00
255      format(/,'#T(s)  Mean_SA(g) for Rcdbar(km), Mbar, haz, eps0, GMPE: ',f6.1,1x,f5.2,1x,e11.5,
     + 1x,f6.2,1x,a)
c the mean spectra are weighed by event frequency branch wt, etc, but not by prob. of exceed.
      do jp=1,kmsp(0)
      st(jp)=meanspecs(ir,im,ie,jp,0)
      ltmp(jp)=lmsp(0,jp)
c      write(31,257)ms_p(jp),exp(st(jp)/prx0)
257        format(f7.3,1x,1pe11.4)
       enddo
      write(33)kmsp(0),st,ltmp                !temp storage used to output CMS
       if(dea_attn)then
       do ia=1,nattn
        prx=haz(ir,im,ie,ia)
        if(ia.eq.2)print *,'gail ',prx,m_att(ia)
        if(ia.eq.1)print *,'NAA ',prx,m_att(ia)
       if(prx.gt.tiny)then
         e0=max(-9.99,ebar(ir,im,ieps,ia)/prx)
       rcd=rbar(ir,im,ie,ia)/prx
       xm=mbar(ir,im,ie,ia)/prx
c        write(31,255)rcd,xm,prx,e0,att_typ(iatten(ia))
                  imout=nint((xm-5.049)*10)+1
      write(33)m_att(ia),ir,imout,ie,rcd,xm,prx,e0
      do jp=1,kmsp(ia)
      st(jp)=meanspecs(ir,im,ie,jp,ia)
      ltmp(jp)=lmsp(ia,jp)
c      write(31,257)ms_p(jp),exp(st(jp)/prx)
       enddo
       write(33)kmsp(ia),st,ltmp                !temp storage used to output CMS
       endif        !latest prx>tiny
       enddo        !attn models
       endif        !deaggregate the individ attn models dea_attn?
       endif        !2nd latest prx>tiny
       endif        !l_ms
       enddo
       enddo
       enddo
   360       continue
       close(21)
       if(geog)then
c Write the data that are accumulated in rstar,mstar, etc.
        do im=1,2        !just 2 bins for size, M<8.8 and all else
        do inow=1,24        !combines R and azimuth
        prx=hstar(inow,im)
        if(prx.gt.1e-15)then
c baz as defined generally works but a few anomalies are known.
        baz=atan2(azbar(inow,im,2),azbar(inow,im,1))/coef
c the below kludge attempts to correct some known bad backazimths
        if(mendo.and.rlond.gt.-120.)then
        baz=-85.
        elseif(mendo.and.baz.lt.-80.)then
        baz=max(-5.-5*(rlond+123.),baz+90.)
        endif
        write(41,47)im,prx,rstar(inow,im)/prx,mstar(inow,im)/prx,baz,msg(im)
 47        format(i3,1x,e11.5,1x,f7.1,1x,f5.2,1x,f7.2,1x,a)
         endif
         enddo
         enddo
         close(41)
         elseif(l_ms)then
c         close(31)
         close(33)        !binary version of 31, used for combining
c         print *,'Mean conditional spectrum for binned sources in *.SPCTRA'
c         print *,'Mean conditional spectrum for mean of all sources in *.MEANSPC'
        prx=mean2haz(0)
c       if(prx.gt.tiny)then
c       eps00=mean2e0(0)/prx
c       rcd=mean2r(0)/prx
c       xm=mean2m(0)/prx
c        write(32,255)rcd,xm,prx,eps00,att_typ(0)
c the mean spectra are weighed by prob. of exceed.
c      do jp=1,npnga
c      if (mean2spec(jp,0) .ne. 0.)write(32,257)ms_p(jp),exp(mean2spec(jp,0)/prx)
c       enddo
c       if(dea_attn)then
c       do ia=1,nattn
c        prx=mean2haz(ia)
c       if(prx.gt.tiny)then
c       eps00=mean2e0(ia)/prx
c       rcd=mean2r(ia)/prx
c       xm=mean2m(ia)/prx
c        write(32,255)rcd,xm,prx,eps00,att_typ(iatten(ia))
c      do jp=1,npnga
c      if(mean2spec(jp,ia).ne.0.)write(32,257)ms_p(jp),exp(mean2spec(jp,ia)/prx)
c       enddo
c       endif        !latest prx>tiny
c       enddo        !attn models
c       endif        !deaggregate the individ attn models dea_attn?
c       endif        !2nd latest prx>tiny
        endif        !geog or l_ms deagg data
            end program

      subroutine readmodel(iflt,tlen)
c     yf,xf,zf(iflt,l,i) are the coordinates (latitude, longitude
c         and depth) of the i-th point on the l-th level
c         (top, hinge, bottom) of the iflt-th fault
      common /blk0/ xf(1,5,1000),yf(1,5,1000),zf(1,5,1000), npts
     +(10,5),dist(5),xfhld(5),yfhld(5),zfhld(5), shold(5), v0x,v0y, vnx,
     +vny, t0x(10,2),t0y(10,2),t0z(10,2), tnx(10,2),tny(10,2),tnz(10,2),
     +t0bx(10),t0by(10),t0bz(10), tnbx(10),tnby(10),tnbz(10)
 
      common /blk1/ rx(1,5,7000),ry(1,5,7000),rz(1,5,7000), rsx
     +(1,5,7000),rsy(1,5,7000),rsz(1,5,7000), send(10,5), vx
     +(1,3,7000),vy(1,3,7000), vz(1,3,7000)
 
c there is only one model imod or iflt 
      real tlen
      parameter (coef=3.14159/180.,nlevels=2)
      external fltdst
      tol = 1.e-3
      do 20 l=1,nlevels
         read(1,*) npts(iflt,l)
c         write(*,*) npts(iflt,l)
         do 10 i=1,npts(iflt,l)
            read(1,*) yf(iflt,l,i),xf(iflt,l,i),zf(iflt,l,i)
c             write(*,*) yf(iflt,l,i),xf(iflt,l,i),zf(iflt,l,i)
   10    continue
 
c         set up and resample
 
c         write(*,*) 'resample'
         call resample(iflt,l,npts,xf,yf,zf)
   20 continue
 
c      write(*,*) 'normvec'
 
      call normvec2(nlevels,tol)
c      tlen(iflt) = send(iflt,1)
       tlen= send(iflt,1)
 
      return
      end
 
 
 
      subroutine rupturesetup(imag,jrup,s1,s2)
c
      type ruptparam
          sequence
          real s1
          real s2
          integer ist1
          integer ist2
          real xt1
          real yt1
          real zt1
          real xt2
          real yt2
          real zt2
          real sb1
          real sb2
          integer isb1
          integer isb2
          real xb1
          real yb1
          real zb1
          real xb2
          real yb2
          real zb2
          real v0x
          real v0y
          real vnx
          real vny
      end type ruptparam
      type (ruptparam) :: rupture(50,800)

c
      common /blk0/ xf(1,5,1000),yf(1,5,1000),zf(1,5,1000), npts
     +(10,5),dist(5),xfhld(5),yfhld(5),zfhld(5), shold(5), v0x,v0y, vnx,
     +vny, t0x(10,2),t0y(10,2),t0z(10,2), tnx(10,2),tny(10,2),tnz(10,2),
     +t0bx(10),t0by(10),t0bz(10), tnbx(10),tnby(10),tnbz(10)
 
      common /blk1/ rx(1,5,7000),ry(1,5,7000),rz(1,5,7000), rsx
     +(1,5,7000),rsy(1,5,7000),rsz(1,5,7000), send(10,5), vx
     +(1,3,7000),vy(1,3,7000), vz(1,3,7000)
 
 
 
      common /blk2/ x,y,z,xflt,yflt,zflt,iflt,l
 
      common /blk4/ rupture
 
      parameter (coef=3.14159/180.,nlevels=2)
 
      logical above
 
      external fltdst
 
c      write(*,*) 'rupturesetup'
      tol = 1.e-3
 
 
c     Get x,y,z coordinates of corners of rupture area corresponding to
c         s1,s2
c     s1,s2 are the distances along the top level of the rupture
      ist1 = int(s1) + 2
      ds1 = s1 - int(s1)
      xt1 = rx(iflt,1,ist1) + rsx(iflt,1,ist1)*ds1
      yt1 = ry(iflt,1,ist1) + rsy(iflt,1,ist1)*ds1
      zt1 = rz(iflt,1,ist1) + rsz(iflt,1,ist1)*ds1
      ist2 = int(s2) + 2
      ds2 = s2 - int(s2)
      xt2 = rx(iflt,1,ist2) + rsx(iflt,1,ist2)*ds2
      yt2 = ry(iflt,1,ist2) + rsy(iflt,1,ist2)*ds2
      zt2 = rz(iflt,1,ist2) + rsz(iflt,1,ist2)*ds2
 
c     sb1,sb2 are the distances along the bottom level of the rupture
      sb1 = s1*send(iflt,2)/send(iflt,1)
      sb2 = s2*send(iflt,2)/send(iflt,1)
      isb1 = int(sb1) + 2
      dsb1 = sb1 - int (sb1)
      xb1 = rx(iflt,2,isb1) + rsx(iflt,2,isb1)*dsb1
      yb1 = ry(iflt,2,isb1) + rsy(iflt,2,isb1)*dsb1
      zb1 = rz(iflt,2,isb1) + rsz(iflt,2,isb1)*dsb1
      isb2 = int(sb2) + 2
      dsb2 = sb2 - int(sb2)
      xb2 = rx(iflt,2,isb2) + rsx(iflt,2,isb2)*dsb2
      yb2 = ry(iflt,2,isb2) + rsy(iflt,2,isb2)*dsb2
      zb2 = rz(iflt,2,isb2) + rsz(iflt,2,isb2)*dsb2
 
 
c      write(*,*) xt1,yt1,zt1,xb1,yb1,zb1,xm1,ym1,zm1
c      write(*,*) xt2,yt2,zt2,xb2,yb2,zb2,xm2,ym2,zm2
 
c     Calculate inward pointing vectors perpendicular to end segments
 
 
      cosfac = cos(((yt1+yb1)/2.)*coef)
      v0x =  yb1-yt1
      v0y = -(xb1-xt1)*cosfac
      signfac = 1.
      signfac = sign(signfac, v0x*(xt2-xt1)*cosfac + v0y*(yt2-yt1) )
 
      vmag = sqrt(v0x**2 + v0y**2)
      v0x = v0x * signfac/vmag
      v0y = v0y * signfac/vmag
c      write(*,*) 'signfac,v0x,v0y=',
c     1            signfac,v0x,v0y
 
      cosfac = cos(((yt2+yb2)/2.)*coef)
      vnx =  yb2-yt2
      vny = -(xb2-xt2)*cosfac
      signfac = 1.
      signfac = sign(signfac, vnx*(xt1-xt2)*cosfac + vny*(yt1-yt2) )
 
      vmag = sqrt(vnx**2 + vny**2)
      vnx = vnx * signfac/vmag
      vny = vny * signfac/vmag
c      write(*,*) 'signfac,vnx,vny=',
c     1            signfac,vnx,vny
c what are all of these variables and more important which ones really matter? 
      rupture(imag,jrup)%s1 = s1
      rupture(imag,jrup)%s2 = s2
      rupture(imag,jrup)%ist1 = ist1
      rupture(imag,jrup)%ist2 = ist2
      rupture(imag,jrup)%xt1 = xt1
      rupture(imag,jrup)%yt1 = yt1
      rupture(imag,jrup)%zt1 = zt1
      rupture(imag,jrup)%xt2 = xt2
      rupture(imag,jrup)%yt2 = yt2
      rupture(imag,jrup)%zt2 = zt2
      rupture(imag,jrup)%sb1 = sb1
      rupture(imag,jrup)%sb2 = sb2
      rupture(imag,jrup)%isb1 = isb1
      rupture(imag,jrup)%isb2 = isb2
      rupture(imag,jrup)%xb1 = xb1
      rupture(imag,jrup)%yb1 = yb1
      rupture(imag,jrup)%zb1 = zb1
      rupture(imag,jrup)%xb2 = xb2
      rupture(imag,jrup)%yb2 = yb2
      rupture(imag,jrup)%zb2 = zb2
      rupture(imag,jrup)%v0x = v0x
      rupture(imag,jrup)%v0y = v0y
      rupture(imag,jrup)%vnx = vnx
      rupture(imag,jrup)%vny = vny

      return
      end


      subroutine qkdist(imag,jrup,xpt,ypt,distmin,fltdepmin,distminxy,baz)
c
      type ruptparam
          sequence
          real s1
          real s2
          integer ist1
          integer ist2
          real xt1
          real yt1
          real zt1
          real xt2
          real yt2
          real zt2
          real sb1
          real sb2
          integer isb1
          integer isb2
          real xb1
          real yb1
          real zb1
          real xb2
          real yb2
          real zb2
          real v0x
          real v0y
          real vnx
          real vny
      end type ruptparam
      type (ruptparam) :: rupture(50,800)
        real baz
      common /blk0/ xf(1,5,1000),yf(1,5,1000),zf(1,5,1000), npts
     +(10,5),dist(5),xfhld(5),yfhld(5),zfhld(5), shold(5), v0x,v0y, vnx,
     +vny, t0x(10,2),t0y(10,2),t0z(10,2), tnx(10,2),tny(10,2),tnz(10,2),
     +t0bx(10),t0by(10),t0bz(10), tnbx(10),tnby(10),tnbz(10)
  
      common /blk1/ rx(1,5,7000),ry(1,5,7000),rz(1,5,7000), rsx
     +(1,5,7000),rsy(1,5,7000),rsz(1,5,7000), send(10,5), vx
     +(1,3,7000),vy(1,3,7000), vz(1,3,7000)
  
      common /blk2/ x,y,z,xflt,yflt,zflt,iflt,l
      common /blk4/ rupture
      parameter (pi=3.1415926,coef=pi/180.,nlevels=2)
      logical above, offend
      external fltdst
c      write(*,*) 'qkdist'
      tol = 1.e-3
      x = xpt
      y = ypt
      z = 0.
 
      s1 = rupture(imag,jrup)%s1 
      s2 = rupture(imag,jrup)%s2 
      ist1 = rupture(imag,jrup)%ist1  
      ist2 = rupture(imag,jrup)%ist2  
      xt1 = rupture(imag,jrup)%xt1 
      yt1 = rupture(imag,jrup)%yt1 
      zt1 = rupture(imag,jrup)%zt1  
      xt2 = rupture(imag,jrup)%xt2 
      yt2 = rupture(imag,jrup)%yt2 
      zt2 = rupture(imag,jrup)%zt2  
      sb1 = rupture(imag,jrup)%sb1  
      sb2 = rupture(imag,jrup)%sb2  
      isb1 = rupture(imag,jrup)%isb1 
      isb2 = rupture(imag,jrup)%isb2 
      xb1 = rupture(imag,jrup)%xb1  
      yb1 = rupture(imag,jrup)%yb1  
      zb1 = rupture(imag,jrup)%zb1  
      xb2 = rupture(imag,jrup)%xb2  
      yb2 = rupture(imag,jrup)%yb2  
      zb2 = rupture(imag,jrup)%zb2  
      v0x = rupture(imag,jrup)%v0x  
      v0y = rupture(imag,jrup)%v0y  
      vnx = rupture(imag,jrup)%vnx 
      vny = rupture(imag,jrup)%vny

c     Begin the analysis of where x,y,z is relative to the fault
 
c      write(*,*) 'x,y,s12,s2',x,y,s1,s2
 
c     The logical variable "above" will indicate whether the observation point
c     is directly above the the rupture surface, initially set to true
      above =.true.
      offend = .false.
 
c     Find the closest distance and coordinates of the closest points
c            of the fault level boundaries
 
      fltmindist = 1.e10
      do 20 l=1,nlevels 
 
         if(l.eq.1.) then
            as = s1
            cs = s2
         elseif(l.eq.2) then
            as = sb1
            cs = sb2
         else
            as = 0.
            cs = send(iflt,l)
         endif
         call mnbrak(as,bs,cs,fa,fb,fc,fltdst,tol)
c         write(*,*) 'as,bs,cs,fa,fb,fc  ',as,bs,cs,fa,fb,fc
         if(bs.eq.cs) then
            smin=cs
            dist2=fltdst(smin)
         else
            dist2=golden(as,bs,cs,fltdst,tol,smin)
c            write(*,*) 'as,bs,cs,smin  ',as,bs,cs,smin
c             write(*,*) 'yflt,xflt  ', yflt,xflt
         endif
         dist(l)=sqrt(dist2)
c         write(*,*) 'l,smin,dist(l) = ',l,smin,dist(l),yflt,xflt
         shold(l) = smin
         xfhld(l) = xflt
         yfhld(l) = yflt
         zfhld(l) = zflt
         if(dist(l).lt.fltmindist) then
            fltmindist = dist(l)
            ifltmin = l
         endif
c         write(*,*) 'l,yfhld(l)',l,yfhld(l)
   20 continue
 
c     Where is x,y,z relative to the top fault boundary?
      cosfac=cos((yfhld(1) + y)*coef/2.)
      xdiff1 = (x-xfhld(1))*cosfac*111.195
      ydiff1 = (y-yfhld(1))*111.195
      zdiff1 =  z-zfhld(1)
      is = int(shold(1))+1
      if(shold(1).eq.send(iflt,1)) is=is+1
      dotprod1 = xdiff1 * vx(iflt,1,is) + ydiff1 * vy(iflt,1,is) +
     +zdiff1 * vz(iflt,1,is)
 
c      write(*,*) 'dotprod1=',dotprod1
      dotprod01 = xdiff1 *vx(iflt,1,is) + ydiff1 * vy(iflt,1,is)
 
      if(dotprod01.le.0.) above = .false.
      if(dotprod1.lt.0) then
c        point on "updip" side of top fault boundary
         ireg=1
         distmin=dist(1)
         ifltdepmin = zfhld(1)
            baz=-pi*0.5-atan2(ydiff1,xdiff1)
c above baz is just a fillin until the master speaks            
 
      else
c     Where is x,y,z relative to the bottom fault boundary?
         cosfac=cos((yfhld(2) + y)*coef/2.)
         xdiff2 = (x - xfhld(2))*cosfac*111.195
         ydiff2 = (y - yfhld(2))*111.195
         zdiff2 = z - zfhld(2)
         is = int(shold(2))+1
         if(shold(2).eq.send(iflt,2)) is=is+1
         dotprod2 = xdiff2 * vx(iflt,2,is) + ydiff2 * vy(iflt,2,is) +
     +   zdiff2 * vz(iflt,2,is)
 
c          write(*,*) 'dotprod2=', dotprod2
         dotprod02 = xdiff2 * vx(iflt,2,is) + ydiff2 * vy(iflt,2,is)
 
         if(dotprod02.ge.0.) above = .false.
         if(dotprod2.gt.0.) then
c             point on "downdip" side of bottom fault boundary
            ireg=2
            distmin=dist(2)
            fltdepmin=zfhld(2)
            baz=-pi*0.5-atan2(ydiff2,xdiff2)
         else
            cosfac=cos((yt1 + y)*coef/2.)
            xdiff4 = (x-xt1)*cosfac*111.195
            ydiff4 = (y-yt1)*111.195
 
            dotprod4 = v0x*xdiff4 + v0y*ydiff4
            cosfac=cos((yt2 + y)*coef/2.)
            xdiff5 = (x-xt2)*cosfac*111.195
            ydiff5 = (y-yt2)*111.195
 
            dotprod5 = vnx*xdiff5 + vny*ydiff5
 
            if(dotprod5.lt.0.) above = .false.
            if(dotprod4.lt.0.) above = .false.
 
               if(dotprod4.le.0.) then
                  ireg = 4
                  baz=-pi*0.5-atan2(ydiff4,xdiff4)
               else
                  if(dotprod5.le.0.) then
                     ireg = 5
                     baz=-pi*0.5-atan2(ydiff5,xdiff5)
                  else
                     ireg = 3
                     cosfac=cos((yfhld(1) + yfhld(2))*coef/2.)
                     dxlev1 = (xfhld(2) - xfhld(1))*cosfac*111.195
                     dylev1 = (yfhld(2) - yfhld(1))*111.195
                     dzlev1 = zfhld(2) - zfhld(1)
                     isplus = int(shold(2)) + 5
                     if(isplus.gt.send(iflt,2)) isplus = isplus -10
                     dxlev2 = (rx(iflt,2,isplus) - xfhld(1))*cosfac
     +               *111.195
                     dylev2 = (ry(iflt,2,isplus) - yfhld(1))*111.195
                     dzlev2 = rz(iflt,2,isplus) - zfhld(1)
c                        write(*,*) x,y,z,xfhld(2),xfhld(3),
c     1                                   yfhld(2),yfhld(3)
 
c                        write(*,*) dxlev1,dylev1,dzlev1,
c     1                             dxlev2,dylev2,dzlev2,
c     2                             xdiff2,ydiff2,zdiff2
                     dxcros = dylev1*dzlev2-dzlev1*dylev2
                     dycros = -(dxlev1*dzlev2-dzlev1*dxlev2)
                     dzcros = dxlev1*dylev2-dylev1*dxlev2
c                        write(*,*) xdiff2,ydiff2,zdiff2,dxcros,dycros,
c     1                             dzcros
c                        write(*,*) dxlev1,dylev1,dzlev1,
c     1                             dxlev2,dylev2,dzlev2
                     crosmag = sqrt(dxcros**2+dycros**2+dzcros**2)
                     distmin = abs ( (dxcros*xdiff2 + dycros*ydiff2 +
     +               dzcros*zdiff2)) / crosmag
                     baz=-pi*0.5-atan2(dylev1,dxlev1)
c a temp fillin until Wesson or ? provides a better estimate.           
 
                     fltdepmin = distmin * dzcros / crosmag
 
                  endif
               endif

 
c     Calculate distance for off end points
            if(ireg.eq.4) then
               offend = .true.
               cosfac = cos((y + yt1)*coef/2.)
               ax = (x - xt1)*cosfac*111.195
               ay = (y - yt1)*111.195
               az = - zt1
               cosfac = cos((yb1 + yt1)*coef/2.)
               bx = (xb1 - xt1)*cosfac*111.195
               by = (yb1 - yt1)*111.195
               bz = zb1 - zt1
               bmag = sqrt(bx*bx + by*by + bz*bz)
               cx = ay*bz - az*by
               cy = -(ax*bz - az*bx)
               cz = ax*by - ay*bx
               distmin = sqrt(cx*cx + cy*cy + cz*cz)/bmag
               fltdepmin = zt1 + ( ax*bx+ay*by+az*bz)*(zb1-zt1)
     +         /(bmag*bmag)
               distminxy = abs(dotprod4) 
               baz=-pi*0.5-atan2(by,bx)
c a  guess. use bottom of fault because most sites closer to bottom...                
 
 
            elseif(ireg.eq.5) then
               offend = .true.
               cosfac = cos((y + yt2)*coef/2.)
               ax = (x - xt2)*cosfac*111.195
               ay = (y - yt2)*111.195
               az = - zt2
               cosfac = cos((yb2 + yt2)*coef/2.)
               bx = (xb2 - xt2)*cosfac*111.195
               by = (yb2 - yt2)*111.195
               bz = zb2 - zt2
               bmag = sqrt(bx*bx + by*by + bz*bz)
               cx = ay*bz - az*by
               cy = -(ax*bz - az*bx)
               cz = ax*by - ay*bx
c                  write(*,*) ax,ay,az,bx,by,bz,icx,cy,cz,bmag
               distmin = sqrt(cx*cx + cy*cy + cz*cz)/bmag
               fltdepmin = zt2 + ( ax*bx+ay*by+az*bz)*(zb2-zt2)
     +         /(bmag*bmag)
               distminxy = abs(dotprod5)
               baz=-pi*0.5-atan2(ay,ax)
c temp off to the west...                
 
  
            endif
 
 
         endif
 
 
 
      endif
c     Final checks
c     Check for wierdness off ends
      if(fltmindist.lt.distmin) then
         distmin = fltmindist
         if(ifltmin.eq.1) ireg=1
         if(ifltmin.eq.2) ireg=2
      endif
 
      if(above.and.(ireg.eq.1.or.ireg.eq.2)) then
         cosfac=cos((yb1 + y)*coef/2.)
         xdiff4 = (x-xb1)*cosfac*111.195
         ydiff4 = (y-yb1)*111.195
 
         dotprod4 = v0x*xdiff4 + v0y*ydiff4
         cosfac=cos((yb2 + y)*coef/2.)
         xdiff5 = (x-xb2)*cosfac*111.195
         ydiff5 = (y-yb2)*111.195
 
         dotprod5 = vnx*xdiff5 + vny*ydiff5
 
         if(dotprod5.lt.0.) above = .false.
         if(dotprod4.lt.0.) above = .false.
      endif
      if(above) then
         distminxy = 0.
         baz=0.        !no calculation needed.
      elseif(.not.offend)then
         xc=(x-xfhld(ifltmin))*cosfac
         yc=y-yfhld(ifltmin)
         distminxy = sqrt(xc** 2+ yc**2) * 111.195
         baz=-pi*0.5-atan2(yc,xc)
      endif
   30 continue
   40 return
 
      end
 
 
 
 
      subroutine resample(iflt,l,npts,xf,yf,zf)
 
      dimension xf(1,5,1000),yf(1,5,1000),zf(1,5,1000)
c     yf,xf,zf(iflt,l,i) are the coordinates (latitude, longitude
c         and depth) of the i-th point on the l-th level
c         (top, hinge, bottom) of the iflt-th fault
 
      dimension npts(10,5)
 
      dimension s(7000),xs(7000),ys(7000),zs(7000)
 
      common /blk1/ rx(1,5,7000),ry(1,5,7000),rz(1,5,7000), rsx
     +(1,5,7000),rsy(1,5,7000),rsz(1,5,7000), send(10,5), vx
     +(1,3,7000),vy(1,3,7000), vz(1,3,7000)
 
 
 
 
      parameter (coef=3.14159/180.)
 
      s(1)= 0.
      do 10 i=2,npts(iflt,l)
         ydiff = yf(iflt,l,i)-yf(iflt,l,i-1)
         cosfac=cos((yf(iflt,l,i-1)+.5*ydiff)*coef)
 
         s(i) = s(i-1) + 111.195* sqrt((cosfac*(xf(iflt,l,i)-xf
     +   (iflt,l,i-1)))** 2+ ydiff**2)
 
         xs(i) = (xf(iflt,l,i)-xf(iflt,l,i-1))/ (s(i)-s(i-1))
 
         ys(i) = (yf(iflt,l,i)-yf(iflt,l,i-1))/ (s(i)-s(i-1))
 
         zs(i) = (zf(iflt,l,i)-zf(iflt,l,i-1))/ (s(i)-s(i-1))
 
c          write(*,*) i,xf(iflt,l,i),yf(iflt,l,i),s(i),xs(i),ys(i)
 
   10 continue
      send(iflt,l) = s(npts(iflt,l))
 
c     Resample faults at points 1 km along trace
      rx(iflt,l,1)= xf(iflt,l,1)
      ry(iflt,l,1)= yf(iflt,l,1)
      rz(iflt,l,1)= zf(iflt,l,1)
      intv=2
      nptsr = int(s(npts(iflt,l))) + 2
c      write(*,*) 's(npts(iflt,l)),nptsr= ',
c     1            s(npts(iflt,l)),nptsr
      do 40 i=2,nptsr-1
         s1=i-1
         do 20 j=intv,npts(iflt,l)
            if(s1.le.s(j)) then
               intv=j
               goto 30
            endif
   20    continue
   30    ds= s1-s(intv-1)
         rx(iflt,l,i)= xf(iflt,l,intv-1) + ds*xs(intv)
         ry(iflt,l,i)= yf(iflt,l,intv-1) + ds*ys(intv)
         rz(iflt,l,i)= zf(iflt,l,intv-1) + ds*zs(intv)
         rsx(iflt,l,i)= rx(iflt,l,i) - rx(iflt,l,i-1)
         rsy(iflt,l,i)= ry(iflt,l,i) - ry(iflt,l,i-1)
         rsz(iflt,l,i)= rz(iflt,l,i) - rz(iflt,l,i-1)
c        write(*,*) s1,intv,xf(iflt,l,intv-1),xs(intv),
c     1             yf(iflt,l,intv-1),
c     1          ys(intv),rx(iflt,l,i),ry(iflt,l,i)
   40 continue
 
      rx(iflt,l,nptsr)=xf(iflt,l,npts(iflt,l))
      ry(iflt,l,nptsr)=yf(iflt,l,npts(iflt,l))
      rz(iflt,l,nptsr)=zf(iflt,l,npts(iflt,l))
      ds= s(npts(iflt,l)) - aint(s(npts(iflt,l)))
      rsx(iflt,l,nptsr)= (rx(iflt,l,nptsr)- rx(iflt,l,nptsr-1))/ds
 
      rsy(iflt,l,nptsr)= (ry(iflt,l,nptsr)- ry(iflt,l,nptsr-1))/ds
 
      rsz(iflt,l,nptsr)= (rz(iflt,l,nptsr)- rz(iflt,l,nptsr-1))/ds
 
      return
      end
 
 
 
 
      function fltdst(s)
c     On return fltdst is the squared distance from x,y,z to
c     a point on the fault, xflt(s),yflt(s),zflt(s)
c     s is arc length along the fault
 
      common /blk1/ rx(1,5,7000),ry(1,5,7000),rz(1,5,7000), rsx
     +(1,5,7000),rsy(1,5,7000),rsz(1,5,7000), send(10,5), vx
     +(1,3,7000),vy(1,3,7000), vz(1,3,7000)
 
 
 
      common /blk2/ x,y,z,xflt,yflt,zflt,iflt,l
      parameter (coef=3.14159/180.,distsq=111.195**2)
 
 
c      write(*,*) 's= ',s
      intv=int(s)+2
   10 ds=s-aint(s)
c      write(*,*) intv,ds
      xflt = rx (iflt,l,intv-1) + ds * rsx(iflt,l,intv)
      yflt = ry (iflt,l,intv-1) + ds * rsy(iflt,l,intv)
      zflt = rz (iflt,l,intv-1) + ds * rsz(iflt,l,intv)
 
      ydiff = y-yflt
      cosfac = cos((yflt+.5*ydiff)*coef)
      fltdst = (((x-xflt)*cosfac)**2 + ydiff**2 )*distsq + (z-zflt)**2
c      write(*,*) 'xflt,yflt,zflt,flddst= ', xflt,yflt,zflt,fltdst
      return
      end
 
 
      subroutine normvec2(nlevel,tol)
c     Find unit vectors normal to top of fault pointing in down
c          dip direction
c     At ends of fault, make vector parallel to end segment
      common /blk1/ rx(1,5,7000),ry(1,5,7000),rz(1,5,7000), rsx
     +(1,5,7000),rsy(1,5,7000),rsz(1,5,7000), send(10,5), vx
     +(1,3,7000),vy(1,3,7000), vz(1,3,7000)
 
 
      common /blk2/ x,y,z,xflt,yflt,zflt,iflt,lf
 
      parameter (coef=3.14159/180.)
 
      external fltdst
 
      iflt = 1
 
c     lv will refer to the level along which you seek the normal vectors
c     lf will refer to other level used to define the dip (i.e.the level
c          below if lv=1 or 2, above if lv =3)
 
 
      do 20 lv=1,nlevel
         if(lv.lt.nlevel) then
            lf = lv + 1
            signfac = 1.
            ls = 1
         else
            lf = lv - 1
            signfac = -1.
            ls = 2
         endif
c       write(*,*) iflt,int(send(iflt,lv))+2
         do 10 i=1,int(send(iflt,lv))+2
 
            x = rx(iflt,lv,i)
            y = ry(iflt,lv,i)
            z = rz(iflt,lv,i)
 
            as = 0.
            cs = send(iflt,lf)
 
c         write(*,*) lv,lf
c         write(*,*) 'x,y,z', x,y,z
            call mnbrak(as,bs,cs,fa,fb,fc,fltdst,tol)
c         write(*,*) 'as,bs,cs,fa,fb,fc  ',as,bs,cs,fa,fb,fc
            if(bs.eq.cs) then
               smin=cs
               dist2=fltdst(smin)
            else
               dist2=golden(as,bs,cs,fltdst,tol,smin)
c            write(*,*) 'as,bs,cs,smin  ',as,bs,cs,smin
c            write(*,*) 'yflt,xflt  ', yflt,xflt
            endif
c         write(*,*) 'x,y,z,xflt,yflt,zflt',x,y,z,xflt,yflt,zflt
            cosfac = cos(ry(iflt,lv,i)*coef)
 
            vx(iflt,lv,i) = signfac * (xflt - x) *cosfac * 111.195
            vy(iflt,lv,i) = signfac * (yflt - y) * 111.195
            vz(iflt,lv,i) = signfac * (zflt - z)
            vmag =sqrt(vx(iflt,lv,i)** 2+ vy(iflt,lv,i)** 2+ vz
     +      (iflt,lv,i)**2)
c     1                  iflt,l,i,rx(iflt,l,i),rx(iflt,l,i-1),
c     2                  ry(iflt,l,i),ry(iflt,l,i-1),
c     3                  vx(iflt,l,i),vy(iflt,l,i)
            vx(iflt,lv,i) = vx(iflt,lv,i)/vmag
            vy(iflt,lv,i) = vy(iflt,lv,i)/vmag
            vz(iflt,lv,i) = vz(iflt,lv,i)/vmag
 
 
 
c           write(*,*) 'vx(iflt,lv,i),vy(iflt,lv,i)',
c     1                vx(iflt,lv,i),vy(iflt,lv,i)
 
 
   10    continue
   20 continue
      return
      end
 
 
 
      function zbrent(func,x1,x2,tol)
      parameter (itmax=100,eps=3.e-8)
c      write(*,*) "zbrent"
      a=x1
      b=x2
c      write(*,*) "a,b", a,b
      fa=func(a)
      fb=func(b)
c      write(*,*) "fa,fb", fa,fb
      if(fb*fa.gt.0.) write(6,*) 'Root must be bracketed for zbrent.'
      fc=fb
      do 10 iter=1,itmax
         if(fb*fc.gt.0.) then
            c=a
            fc=fa
            d=b-a
            e=d
         endif
         if(abs(fc).lt.abs(fb)) then
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
         endif
         tol1=2.*eps*abs(b)+0.5*tol
         xm=.5*(c-b)
         if(abs(xm).le.tol1 .or. fb.eq.0.) then
            zbrent=b
            return
         endif
         if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
            s=fb/fa
            if(a.eq.c) then
               p=2.*xm*s
               q=1.-s
            else
               q=fa/fc
               r=fb/fc
               p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
               q=(q-1.)*(r-1.)*(s-1.)
            endif
            if(p.gt.0.) q=-q
            p=abs(p)
            if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
               e=d
               d=p/q
            else
               d=xm
               e=d
            endif
         else
            d=xm
            e=d
         endif
         a=b
         fa=fb
         if(abs(d) .gt. tol1) then
            b=b+d
         else
            b=b+sign(tol1,xm)
         endif
         fb=func(b)
   10 continue
      write(6,*) 'zbrent exceeding maximum iterations.'
      zbrent=b
      return
      end
 
  
      SUBROUTINE MNBRAK(AX,BX,CX,FA,FB,FC,FUNC,tol)
c     Modified from Numerical Recipes
c     Want to establish whether minimum is within interval [ax,cx],
c     or at end point.  If minimum is within interval, then ax,bx,cx
c     will bracket minimum
      PARAMETER (GOLD=0.618034, GLIMIT=100., TINY=1.E-20)
      FA=FUNC(AX)
      FC=FUNC(CX)
      IF(FC.GT.FA)THEN
         DUM=AX
         AX=CX
         CX=DUM
         DUM=FC
         FC=FA
         FA=DUM
      ENDIF
   10 BX=AX+GOLD*(CX-AX)
      if(abs(cx-bx).le.tol) then
         bx=cx
         go to 30
      endif
      FB=FUNC(BX)
   20 IF(FB.LE.FC)THEN
         go to 30
c        Minimum is bracketed by ax,bx,cx
      else
c        Replace ax with bx, and choose new bx
         ax=bx
         fa=fb
         go to 10
      endif
   30 return
      end
 
 
 
 
      FUNCTION GOLDEN(AX,BX,CX,F,TOL,XMIN)
      PARAMETER (R=.61803399,C=.38196602)
      X0=AX
      X3=CX
      IF(ABS(CX-BX).GT.ABS(BX-AX))THEN
         X1=BX
         X2=BX+C*(CX-BX)
      ELSE
         X2=BX
         X1=BX-C*(BX-AX)
      ENDIF
      F1=F(X1)
      F2=F(X2)
   10 IF(ABS(X3-X0).GT.TOL*(ABS(X1)+ABS(X2)))THEN
         IF(F2.LT.F1)THEN
            X0=X1
            X1=X2
            X2=R*X1+C*X3
            F0=F1
            F1=F2
            F2=F(X2)
         ELSE
            X3=X2
            X2=X1
            X1=R*X2+C*X0
            F3=F2
            F2=F1
            F1=F(X1)
         ENDIF
         GOTO 10
      ENDIF
      IF(F1.LT.F2)THEN
         GOLDEN=F1
         XMIN=X1
      ELSE
         GOLDEN=F2
         XMIN=X2
      ENDIF
      RETURN
      END
 
 
 
 
      FUNCTION BRENT(AX,BX,CX,F,TOL,XMIN)
c
c     Numerical Recipes, 1st Edition
c
      PARAMETER (ITMAX=100,CGOLD=.3819660,ZEPS=1.0E-10)
      A=MIN(AX,CX)
      B=MAX(AX,CX)
      V=BX
      W=V
      X=V
      E=0.
      FX=F(X)
      FV=FX
      FW=FX
      DO 30 ITER=1,ITMAX
         XM=0.5*(A+B)
         TOL1=TOL*ABS(X)+ZEPS
         TOL2=2.*TOL1
         IF(ABS(X-XM).LE.(TOL2-.5*(B-A))) GOTO 40
         IF(ABS(E).GT.TOL1) THEN
            R=(X-W)*(FX-FV)
            Q=(X-V)*(FX-FW)
            P=(X-V)*Q-(X-W)*R
            Q=2.*(Q-R)
            IF(Q.GT.0.) P=-P
            Q=ABS(Q)
            ETEMP=E
            E=D
            IF(ABS(P).GE.ABS(.5*Q*ETEMP).OR.P.LE.Q*(A-X).OR. P.GE.Q
     +      *(B-X)) GOTO 10
            D=P/Q
            U=X+D
            IF(U-A.LT.TOL2 .OR. B-U.LT.TOL2) D=SIGN(TOL1,XM-X)
            GOTO 20
         ENDIF
   10    IF(X.GE.XM) THEN
            E=A-X
         ELSE
            E=B-X
         ENDIF
         D=CGOLD*E
   20    IF(ABS(D).GE.TOL1) THEN
            U=X+D
         ELSE
            U=X+SIGN(TOL1,D)
         ENDIF
         FU=F(U)
         IF(FU.LE.FX) THEN
            IF(U.GE.X) THEN
               A=X
            ELSE
               B=X
            ENDIF
            V=W
            FV=FW
            W=X
            FW=FX
            X=U
            FX=FU
         ELSE
            IF(U.LT.X) THEN
               A=U
            ELSE
               B=U
            ENDIF
            IF(FU.LE.FW .OR. W.EQ.X) THEN
               V=W
               FV=FW
               W=U
               FW=FU
            ELSE IF(FU.LE.FV .OR. V.EQ.X .OR. V.EQ.W) THEN
               V=U
               FV=FU
            ENDIF
         ENDIF
   30 CONTINUE
      write(6,*) 'deaggSUBD.2013: Brent exceed maximum iterations.'
   40 XMIN=X
      BRENT=FX
      RETURN
      END
 
 
 

      subroutine Gregor(ip,xmag,rcd,vs30,gnd,sigmaf)
c
c Gregor ea 2002. predicted ground motion for Cascadia subduction megathrust
c compile this subroutine with f95 and with -e or equivalent (extended line-length)
c Programmed Oct 17 2006.
c
c In: ip = period index in per() array below, period T goes to 5 s for tall bldgs
c       xmag=moment mag
c       rcd= closest distance to slab (km)
c       vs30= top30 m Vs, (m/s)
c Out: gnd = ln(median), sigmaf = 1/sigma/sqrt2
c
c Oct 2006 version of Gregor et al, BSSA v92, pp 1923-1932
c
      parameter (pi=3.14159265,sqrt2=1.414213562,vref=760.,np=10)
      real gnd, sigmaf,mfac
      real, dimension(np):: per, cr1,cr2,cr3,cr4,cr5,cr6,sigmart
      real, dimension(np):: cs1,cs2,cs3,cs4,cs5,cs6,sigmast      
      per = (/0.0,0.1,0.2,0.25,0.333,0.5,1.,2.,2.5,5./)
c rock coeffs, table2. avg vs30=363 m/s
      cr1 = (/21.0686,30.005,39.345,37.69,34.787,29.159,6.528,8.657,6.637,8.013/)
      cr2 = (/-1.7712,-2.349,-3.087,-2.96,-2.899,-2.424,-0.406,-0.851,-0.651,-0.943/)
      cr3= (/-5.0631,-6.3862,-7.6002,-7.379,-6.7855,-6.2114,-3.1991,-2.7398,-2.3124,-2.4087/)
      cr4 =(/0.4153,0.5009,0.5972,0.5842,0.5616,0.5216,0.2589,0.2339,0.1879,0.2154/)
      cr5 =(/4.2,4.7,5.1,5.1,4.9,4.7,3.2,2.8,2.8,2.3/)
      cr6 =(/0.0017,-0.0019,0.006,-0.0023,0.0256,0.0161,-0.0225,0.037,0.0364,0.0647/)
      sigmart=(/0.724,0.7954,0.8679,0.8444,0.8776,0.8039,0.7567,0.6305,0.6657,0.773/)
c soil coeffs, table3 avg vs30=182 m/s
      cs1 =(/23.8613,29.9693,75.8218,100.3357,71.7967,56.0088,17.233,17.9124,16.1666,7.4856/)
      cs2 =(/-2.2742,-2.7254,-6.8396,-9.0324,-6.499,-5.1176,-1.5506,-1.7505,-1.5091,-0.836/)
      cs3=(/-4.8803,-5.8054,-12.0687,-15.3511,-11.6056,-9.5083,-4.3287,-3.815,-3.7101,-2.0627/)
      cs4=(/0.4399,0.5098,1.0753,1.3731,1.0415,0.8632,0.393,0.3574,0.3344,0.1779/)
      cs5=(/4.7,5.2,6.3,6.6,6.2,5.9,4.2,4.1,4.1,-0.2/)      !air jordan for 5s
      cs6=(/0.0366,0.0226,0.0096,0-0.0043,0.0102,0.0164,0.0133,0.0583,0.0473,0.0821/)
      sigmast=(/0.5436,0.5926,0.6618,0.6371,0.6431,0.6139,0.6606,0.6276,0.6676,0.8207/)
c for 760 rock need to correct. The authors were emailed 10/17/2006 asking for recommendations on
c how to accomplish this shift to NSHMP generic rock at NEHRP BC boundary. SH.
c
c Another item: depth to slab does not have a specific term unlike Geomatrix relation
      if(vs30.gt.360.)then
      c1=cr1(ip); c2=cr2(ip); c3=cr3(ip); c4=cr4(ip); c5=cr5(ip); c6=cr6(ip)
      sig=sigmart(ip)
      else
      c1=cs1(ip); c2=cs2(ip); c3=cs3(ip); c4=cs4(ip); c5=cs5(ip); c6=cs6(ip)
      sig=sigmast(ip)
      endif
      mfac=(xmag-10.)**3
      fac1=alog(rcd+exp(c5))
      gnd=c1+c2*xmag+(c3+c4*xmag)*fac1+c6*mfac
      sigmaf=1./sig/sqrt2
      return
      end subroutine Gregor

cccccccccccccccccc
      subroutine getSadighR(iq,xmag,rrup,gnd,sigmaf,
     1 sigmanf,distnf)
c modified to perform like nga functions. SH June 2 2006 . 
c added to hazSUBX with mag0 capped at 8.499 oct 2006. sc1 and sd1 corrected feb 27 2007
c Inputs:
c iq is the index of perx() associated with period being run with index ip in main
c rrup or rcd is the distance to fault. no additional distance fac needed.
c  sigmanf, distnf = near source sigma information
c we are not using sigmanf or distnf at this time...
        parameter (sqrt2=1.414213562,rockcoef=0.240336)
c Rock-site coefs only, written below for 7 periods. How fast is that rock?
c walt silva says about 500 m/s (NGA discussion hot topic).
c For soil site condition, use getSadighS instead of this subr.
c Replaced data statements with array constructors oct 2006.
      real sc1(7),sc2(7),sc3(7),sc4(7),sc5(7),sc6(7),perx(7)
      real sd1(7),sd2(7),sd3(7),sd4(7),sd5(7),sd6(7),sd7(7)      
      real ssigma1(7),ssigma2(7),ssigmacoef(7),smagsig(7)
        common/mech/ss,rev,normal,obl
        common/dipinf/dipang,cdipsq
        logical ss,rev,normal,obl
c prepared for 7 periods jun 22 2006. gnd is modified
c for sense-of-slip  from common/mech/ 
      real sthrust(7)/7*.1823/
c      sc1 = (/-0.624,.153,-1.705,0.2750,-0.057,-0.5880,-2.945/)
c sc1 has been reduced somewhat for 10,5,3.33,2, 1 hz. Based on Frankel mods. also pga
      sc1 = (/-0.765,0.048,-1.9699,0.1945,-0.2092,-0.7979,-3.1936/)
      sc2 = (/1.0,1.,1.,1.,1.,1.,1./)
      sc3 = (/0.0,-0.004,-0.055,0.006,-0.017,-0.04,-0.07/)
      sc4 = (/-2.1,-2.08,-1.8,-2.148,-2.028,-1.945,-1.67/)
      sc5 = (/1.29649,1.29649,1.29649,1.29649,1.29649,1.29649,1.29649/)
      sc6 = (/0.25,.25,.25,.25,.25,.25,.25/)
      perx = (/0.0,0.2,1.0,0.1,0.3,0.5,2.0/)
c Sadigh97 sd1 has been modified for pga 1 2 3.3 and 5 10 hz only. 
c Lowered to get a 760 rock from a presumed 520 rock. USE BJF97 siteamp fcn. feb 28
c frankel did not lower the pga. but to be consistent I think it also requires the treatment.
      sd1 = (/-1.4148,-.602,-2.643,-0.4555,-0.8592,-1.4479,-3.8436/)
      sd2 = (/1.1,1.1,1.1,1.1,1.1,1.1,1.1/)
      sd3 = (/0.,-.004,-0.055,0.006,-0.071,-0.040,-0.07/)
      sd4 = (/-2.1,-2.08,-1.8,-2.148,-2.028,-1.945,-1.67/)
      sd5 = (/-.48451,-.48451, -.48451,-.48451,-.48451, -.48451,-.48451/)
      sd6 = (/.524,.524,.524,.524,0.524,.524,.524/)
      sd7 = (/0.,0.,0.,-0.041,0.,0.,0./)
      ssigma1 = (/1.39,1.43,1.53,1.41,1.45,1.5,1.53/)
      ssigmacoef = (/0.14,0.14,0.14,0.14,0.14,0.14,0.14/)
      ssigma2 = (/0.38,0.42,0.52,0.4,0.44,0.49,0.52/)
      smagsig = (/7.21,7.21,7.21,7.21,7.21,7.21,7.21/)
          if(xmag.lt.smagsig(iq))then
           sig= ssigma1(iq)- ssigmacoef(iq)*xmag
          else
           sig= ssigma2(iq)
           endif
          if(rrup.lt.distnf) sig= sig+ sigmanf
          sigmaf= 1./(sig*sqrt2)
          dist= rrup      !ch. april 11 2007
          if(xmag.lt.6.5) then
      gnd= sc1(iq)+sc2(iq)*xmag+sc3(iq)*((8.5-xmag)**2.5)
      gnd= gnd+ sc4(iq)*alog(dist+ exp(sc5(iq)+sc6(iq)*xmag))
           else
        xmag0=min(xmag,8.4999)
      gnd= sd1(iq)+sd2(iq)*xmag0 + sd3(iq)*((8.5-xmag0)**2.5)
      gnd= gnd+ sd4(iq)*alog(dist+ exp(sd5(iq)+sd6(iq)*xmag0))
      if(rev)gnd=gnd + sthrust(iq)
            endif
          if(sd7(iq).ne.0.)gnd= gnd+sd7(iq)*alog(dist0+2.)
      return
      end subroutine getSadighR

      subroutine getSadighS(iq,xmag,rrup,gnd,sigmaf,
     1 sigmanf,distnf)
      parameter (sqrt2=1.414213562,np=7)
c added to hazSUBX with mag0 capped at 8.499. corrected c7(1hz) apr 16 07
c Sadigh et al (SRL 1997), soil site. Separate subroutine for rock, above.
c Modified to look like the new nga subroutines. Original code from T. Cao
c ip is input file period index. 
c The distance  corresponds to Rjb. 
c Coeffs have been written for pga, 0.2, .1,.3, 0.5, 1.0, and 2.0 s. SH)
c iq is index of atten model period corresponding to ip
      real ssigma1(7),ssigma2(7),ssigmacoef(7),smagsig(7)
c this routine is for soil. according to Cao hazFXv3-s.f. With minor changes
c for large subduction eqs.
        common/mech/ss,rev,normal,obl
        common/dipinf/dipang,cdipsq
        logical ss,rev,normal,obl
       real C1SS/-2.170/,C1RV/-1.920/,C2/1./,C3/1.70/
       real C4M1/2.1863/,C4M2/0.3825/, C5M1/0.32/,C5M2/0.5882/
       real C6SS(np),C6RV(np),C7(np),perx(8)
      perx = (/0.0,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
c strike slip term included but is not used for subduction sources.
      C6SS = (/0.,0.9187,0.5665,.6395,.9547,.8494,0.1001/)
      C6RV= (/0.,0.9187,0.5075,.6395,.9547,.8285,-.0526/)
      C7= (/0.,-.004,-0.065,0.005,-.014,-0.033,-0.108/)
      ssigma1= (/1.52,1.565,1.66,1.54,1.58,1.61,1.71/)
      ssigmacoef= (/0.16,0.16,0.16,0.16,0.16,0.16,0.16/)
c      if(vs30.gt.510.)then
c      write(6,*)'getSadighS has been called with probable rock conditions'
c      write(6,*)'Consider calling getSadighR instead; vs30=',vs30
c      endif
c First term is indep of M and R. Subduction is expected to be
c      reverse-slip for this calculation. However, if user makes it SS, so be it.
      if(rev)then
      c1= C1RV+C6RV(iq)
      else
      c1= C1SS + C6SS(iq)
      endif
c-- magnitude limited to 8.5. 
      xmag0=min(xmag,8.5)
      dist = rrup
           sigp= ssigma1(iq)- ssigmacoef(iq)*7.0
          gndm= c1 +C2*xmag0 +C7(iq)*(8.5-xmag0)**2.5
c subduction M>6.5 always. Use the C4M2 and C5M2 versions of C4 and C5
          facm = C4M2*EXP(C5M2*xmag0)
          sigmaf= 1./(sigp*sqrt2)
      gnd= gndm -C3*alog(dist +facm)
      return
      end subroutine getSadighS

      subroutine getGeom(ip,islab,vs30,xmag,dist0,gnd,sigmaf,p_corr)
c update nov 19 2008: add 4 & 5 s.
c  Geomatrix (Youngs et al. subduction). modified to NGA style. option: with
c        petersen-hartzell dist. correction, for dist0 > 200 km.
c        ip = outer loop period index. will need iq = per index in pergeo
c islab=0 for subduction, islab=1 for intraplate.
c        xmag = moment magnitude of source
c        dist0 = slant distance (km)
c        p_corr = logical variable, .true. if petersen correction is used
c                 when R>200 km. Effect is only programmed for 1hz 5hz and pga.
c                Mar 2010.
c modified April 2007. Steve Harmsen. FOR rock or firm or deep soil.
c added coeffs for 0.075, 0.4 0.75 1.5 and 3s feb 1 2008 SH.
c vs30 >= 700  means rock site 
c vs30 > 444 Csoil (avg of rock and soil)
c otherwise, soil site.  rock & soil-coef tri-chotomy is weak. 
c should have continuous variation. 
c returns median and 1/sigma/sqrt2.
c Jan 2011: All coeffs corresponding to T in NGA-W are interpolated to T=5s
      parameter (np=19,sqrt2=1.41421356)
      common/ipindx/iperba,ipgeom,lvs
      common/ms/ia,npnga,jms,l_ms,lmsp,imsp,meanspec,sigspec
        logical p_corr,l_ms,lmsp(0:8,21),lvs
      integer imsp(8,21)
      real, dimension (21) :: meanspec,sigspec
      real gc0/0.2418/,gcs0/-0.6687/,ci/0.3846/,cis/0.3643/
      real gch/0.00607/,gchs/0.00648/,gmr/1.414/,gms/1.438/
      real period,gnd0,gndz,gz,g3,g4,gndm
       real, dimension(np):: gc1,gc2,gc3,gc4,gc5,pergeo,gc1s,gc2s,gc3s,gp
      real vgeo(3)
c array constructors oct 2006. add some periods to build mean spectum dec 2010.
C vgeo is a reference vs30 for Geomatrix, 760 m/s rock, 300 m/s soil.
c Additional siteamp will be wrt these values from Frankel discussion july 7.
c the 475 value isn't currently used. A coeff. set that represents very stiff soil (NEHRP C)
c with something like 475 to 500 m/s could go with this. Currently there is a discontinuity
c when switching between rock and soil coeffs, currently at 520 m/s Vs30.
        vgeo=(/760.,300.,475./)
         pergeo= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,0.075,0.4,0.75,1.5,3.,4.,5.,    
     + 0.02,0.03,0.05,0.15,0.25/)      
       gc1= (/0.,0.722,-1.736,1.1880,0.246,-0.4,-3.3280,1.275,
     + -0.115,-1.149,-2.634,-4.511,-5.350,-6.025,
     +0.4386,0.6952,1.0184,0.9154,0.4600 /)
       gc2= (/0.,-0.0027,-0.0064,-0.0011,-0.0036,-0.0048,-0.0080,0.0,
     + -0.0042,-.0057,-.0073,-.0089,-.0096,-.0101,0.,0.,0.,-0.002,-0.0032/)
       gc1s= (/0.,1.549,-2.87,2.516,0.793,-.438,-6.4330,2.40,0.144,-1.704,
     + -5.101,-6.672,-7.618,-8.352,
     + 0.8256,1.3086,1.9170,1.9503,1.1329/)
       gc2s= (/0.,-0.0019,-0.0066,-0.0019,-0.002,-0.0035,-0.0164,-0.0019,
     + -.002,-.0048,-0.0114,-0.0221,-0.0235,-0.0246,
     + -0.0007,-0.0010,-0.0015,-0.0019,-0.002/)
       gc3= (/-2.556,-2.528,-2.234,-2.6550,-2.454,-2.36,-2.107,-2.707,
     + -2.401,-2.286,-2.16,-2.033,-1.98,-1.939,
     + -2.6079,-2.6383,-2.6766,-2.5807,-2.4873/)
       gc3s= (/-2.329,-2.464,-1.785,-2.697,-2.327,-2.140,-1.29,-2.697,
     + -2.23,-1.95,-1.47,-1.347,-1.272,-1.214,
     + -2.4556,-2.5296,-2.6229,-2.5607,-2.3886/)
       gc4= (/1.45,1.45,1.45,1.45,1.45,1.45,1.55,1.45,
     + 1.45,1.45,1.50,1.65,1.65,1.65,1.45,1.45,1.45,1.45,1.45/)
       gc5= (/-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
     + -0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1/)
c       gp=(/-0.0038,-0.00459,-0.00286,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)
c petersen coeffs were revised during a 2nd regression of the GDSN data. Results generally similar
        gp=(/-.003646,-.004724,-.003208,0.,-.004549,0.,0.,0.,0.,0.,
     + 0.,0.-.001070,0.,0.,0.,0.,0.,0.,0./)
c Also, added 0.3s and 3.0s coeffs because data were also available at these Spectral Ts.
c Always define a reference PGA, determine nonlinear response based on this value
c islab=0 for subduction. skip the irrelevant 0 multiply.
          iq=ipgeom
          gnd0p=gc0                !+ci*islab
c          period=pergeo(iq)
      if(vs30.gt.520.)then
c Use rock coeffs. No interface term ci or cis for subduction
        gnd0=gc0         !+ci*islab
        gz=gch
        g1=gc1(iq)
        g2=gc2(iq)
          g3=gc3(iq)
          g4=1.7818
          ge=0.554
          gm=gmr
          ir=1
      else
c Use soil coeffs
        gnd0=gcs0         !+ cis*islab 
        gz=gchs
        g1=gc1s(iq)
        g2=gc2s(iq)
           g3=gc3s(iq)
           g4=1.097
           ge=0.617
           gm=gms
         ir=2
           endif
c gz term  used for subduction, and for benioff 
          sig= gc4(iq)+gc5(iq)*min(8.,xmag)
          sigmaf= 1./(sig*sqrt2)
c same sigma for soil and rock.
c Thru 2008 we use a constant depth of 20 km but this depth more properly depends on subduction model.
c Some slabs may descend at a steeper or shallower angle affecting this bottom.
          gndm= gnd0 +g1 +gm*xmag +g2*((10.-xmag)**3)  +gz*20.
c add Petersen corrrection if asked for (iatten=-2)
           if(p_corr.and.dist0.gt.200.)gndm=gndm+gp(iq)*(dist0-200.)
          arg= exp(ge*xmag)
c Distance could be hypocentral or distance to top-of-Benioff zone.
          gnd=gndm +g3*alog(dist0+g4*arg)
c additional code appended to getGeom...
      if(lvs )then
c frankel mods for nonlin siteamp July 2009. Use g4p (rock g4) SH.
c Gndzp for nonlin site response.
        gndzp= gnd0p+20.*gch +gc1(1)
        xmagfac=(10.-xmag)**3
          gndmp= gndzp+gmr*xmag+ gc2(1)*xmagfac
          argp= exp(0.554*xmag)
              gndp=gndmp+gc3(1)*alog(dist0+1.7818*argp)
              pganl= exp(gndp)
              gnd= gnd + basiteamp(pganl,vs30,vgeo(ir),jms)
      endif        !vs30 .ne. vgeo(ir), the reference Vs for soil or rock, geom.
              gnds = gnd; sigmas=sig
c play it again Bob Youngs for the mean spectrum.
        if (l_ms)then
        meanspec(jms)=gnds
        sigspec(jms)=sigmas
        do jp=1,npnga-2        !cant go out to 10 s. 5s max.
c        print *,jms,meanspec(jms),sigspec(jms),pergeo(iq)
        if(lmsp(ia,jp).and.jp.ne.jms)then
c compute mean g.m. for the other nga periods.
        j=imsp(ia,jp)        !j = the period index for the geomatrix coeffs.
      if(vs30.gt.520.)then
c Use rock coeffs. No interface term ci or cis for subduction
        g1=gc1(j)
        g2=gc2(j)
          g3=gc3(j)
      else
c Use soil coeffs
        g1=gc1s(j)
        g2=gc2s(j)
           g3=gc3s(j)
           endif
c gz term  used for subduction, and for benioff 
          sig= gc4(j)+gc5(j)*min(8.,xmag)
c same sigma for soil and rock.
c Thru 2008 we use a constant depth of 20 km but this depth more properly depends on subduction model.
c Some slabs may descend at a steeper or shallower angle affecting this bottom.
          gndm= gnd0 +g1 +gm*xmag +g2*xmagfac  +gz*20.
c add Petersen corrrection if asked for (iatten=-2)
           if(p_corr.and.dist0.gt.200.)gndm=gndm+gp(j)*(dist0-200.)
c          arg= exp(ge*xmag)
c Distance could be hypocentral or distance to top-of-Benioff zone.
          gndx=gndm +g3*alog(dist0+g4*arg)
c additional code appended to getGeom...
      if(lvs )then
c frankel mods for nonlin siteamp July 2009. Use g4p (rock g4) SH.
c Gndzp for nonlin site response. these commented-out lines have already been computed above.
c        gndzp= gnd0p+20.*gch +gc1(1)
c          gndmp= gndzp+gmr*xmag+ gc2(1)*xmagfac
c          argp= exp(0.554*xmag)
c              gndp=gndmp+gc3(1)*alog(dist0+1.7818*argp)
c              pganl= exp(gndp)
              gndx= gndx + basiteamp(pganl,vs30,vgeo(ir),jp)
      endif        !vs30 .ne. vgeo(ir), the reference Vs for soil or rock, geom.
      meanspec(jp)=gndx; sigspec(jp)=sig
      endif        !can be caluclated.
      enddo        !meanspec periods.
        endif !if (l_ms)
      return
      end subroutine getGeom

      subroutine getABsub(iq,ir,islab,vs30,xmag,dist0,gnds,sigmaf,p_corr)
c
c from At&Boore BSSA aug 2003, p1715.  C-site and D-site look good apr 11 2007. SH.
c Modify 5hz and add 2.5 hz coeffs., 10/10/2008. AB BSSA erratum Oct 2008. 
c Affects Global but not cascadia. 10/17 also modify 2 hz and 3.33 hz coeffs for global.
c
c dec 27  2010:: add conditional mean spectrum calcs
c Apr 11 2007: added CD and DE boundary response, like the BC boundary
c This modification yields a smoother transition between site classes. SH.
c
c Inputs:
c xmag = moment magnitude, should be > 7.5 according to doc.
c dist0 is R_cd, see p 1706, A&B BSSA aug 2003
c vs30 is top 30 m shear-wave velocity, for determining site class
c iq = period index for coefficients
c 
c ir=1 PNW Cascadia coeff, BC to E site class
c ir=3 Worldwide, BC to E site class
c ir=5 Japan, not programmed. ir=2 or 4 not in use at this time
c islab=0 interface; islab=1 intraslab. have totally different coef sets.
c  Atkinson and Boore subduction zone intraslab or interface.
c p_corr = .true. to use anelastic decay for R>200 from GDSN regr. Mar 2010.
c
c Outputs:
c gnd = ln(median sa or pga)
c sigmaf = sigma factor = 1.0/sigma/sqrt2/aln10
c            where sigma is the log_10 aleatory sigma.
c  Coeffs for 18 spectral pds. also 3s available for interface
c modified for gfortran, f95 Oct 2006. added several periods for mean spec. 12/2010
      parameter (np=18,sqrt2=1.4142136,gfac=2.9912261,aln10=2.30258509)
c gfac = log10(980)
      common/ms/ia,npnga,jms,l_ms,lmsp,imsp,meanspec,sigspec
        logical p_corr,l_ms,lmsp(0:8,21)
      integer imsp(8,21)
      real, dimension (21) :: meanspec,sigspec
c these are two frequencies where special pleading was introduced 10-2008 BSSA.
      real rc1/-.25/,rs1/2.991/
      real rc2/0.6909/,rs2/0.03525/
      real rc3/0.01130/,rs3/0.00759/
      real rc4/-0.00202/,rs4/-0.00206/
c  Dtor allows a distribution if 
c you are uncertain about what that depth is. not using this variable dtor here.
c Global modified coeffs ab 10-2008 affected at 2, 2.5, 3.33, 5 hz.
      real, dimension(np) :: pcor, s1g,s2g,s3g,s4g,s5g,s6g,s7g
c the regressions were done with resids. from the Global model.
      real c1w(np)
c interface, names start with s for "Subduction". s1g global, s1 cascadia only.
c
      real,dimension(np):: s1,s2,s3,s4,s5,s6,s7,ssig
      real perab(18),period
c array constructors oct 2006. 2.5 hz is element number 6 of 10. 1.33 hz is el. 8
      perab= (/0.,0.04,0.2,1.0,0.1,0.3,0.4,0.5,0.75,2.0,3.0,
     + 0.02,0.03,0.05,0.075,0.15,0.25,1.50/)      
c  Inslab coeffs not used here but provided for ref.
c      c1w=(/-0.04713,0.51589,-1.02133,0.43928,0.26067,0.-0.16568,-2.39234,-3.70012/)
c      c1= (/ -0.25,0.40,-0.98,0.160,0.195,-0.172,-2.250,-3.64/)
c      c2= (/0.6909,0.69186,0.8789,0.66675,0.73228,0.7904,0.99640,1.1169/)
c       c3= (/0.01130,0.00572,0.00130,0.0108,0.00372,0.00166,0.00364,.00615/)
c       c4= (/-0.00202,-0.00192,-0.00173,-0.00219,-0.00185,-0.00177,-0.00118,-.00045/)
c      c5= (/0.19,0.15,0.10,.15,0.140,0.125,0.100,0.1 /)
c       c6= (/0.24,0.27,0.30,0.27,0.26,0.23,0.28,0.25/)      !check these
c      sig= (/0.27,0.28,0.29,.28,0.280,0.282,0.300,0.3/)      !BASE 10 SIGMA
c For subduction events, recommendation is M7.5 up and R<300 km.
c coefficients for subduction, Table 1 p 1715 & T3, p 1726. 
c Definitions: s1g global, s1 Cascadia.
c 10/17/2008: Cubic Splines were used for several 3.33 and 2 hz Global-estimation coefs. 
c Interested in the details? See intAB03.f for src code uses Numerical recipes. 
c      s1g=(/2.991,2.6638,2.1442,2.7789,2.5816,2.5249,2.3857867,2.1907,2.301/) prior
c s1g has been recomputed for 2.0, 2.5, 3.33 and 5hz. Ditto s2g ,...
      s1g=(/2.991,2.8753,2.5711536,2.1442,2.7789,2.6168785,2.6175463,2.536019,2.288355,2.1907,2.301,
     +  2.93315,2.89931,2.85182 ,2.80917,2.65738,2.59632,2.17140/)
      s1=(/2.79,2.6,2.54,2.18,2.5,2.516,2.50,2.418,2.241635,2.33,2.36,
     + 2.69500,2.63943,2.57565,2.53140,2.52340,2.52679,2.26774/)
      s2=(/.03525,.07052,.12386,.1345,.09841,.1373,0.1477,.1444,0.14504924,.07148,.02237,
     + 0.05288,0.06320,0.07731,0.08965,0.11330,0.13126,0.09764/)
      s2g=(/.03525,.07052,0.13976128,.1345,.09841,0.13694176,0.13179871,0.1324168,0.13728201,.07148,.02237,
     + 0.05288,0.06320,0.07731,0.08965,0.12260,0.13821,0.09764/)
      s3=(/.00759,.01004,.00884,.00521,.00974,.00789,0.00728,.00671,5.9208343E-3,.00224,.00012,
     + 0.00881,0.00953,0.00997,0.00983,0.00921,0.00832,0.00347/)
      s3g=(/.00759,.01004,7.79948E-3,.00521,.00974,8.12251E-3,8.32052E-3,7.951397E-3,6.429092E-3,.00224,.00012,
     +0.00881,0.00953,0.00997,0.00983,0.00860,0.00798,0.00347/)
      s4=(/-0.00206,-0.00278,-.0028,-0.0011,-.00287,-.00252,-0.00235,-.00195,-1.5903986E-3,0.,0.,
     + -0.00242,-0.00263,-0.00280,-0.00284,-0.00283,-0.00265,-0.00046/)
      s4g=(/-0.00206,-0.00278,-2.49985E-3,-0.0011,-.00287,-2.6353553E-3,-2.6501498E-3,-2.4560375E-3 ,-1.7370114E-3,0.,0.,
     + -0.00242,-0.00263,-0.00280,-0.00284,-0.00265,-0.00257,-0.00046/)
      s5=(/0.19,0.15,0.15,0.10,0.15,0.14,0.13,0.12,0.106354214,0.10,0.10,
     + 0.170,0.1583,0.1500,0.1500, 0.150, 0.1445,0.10/)
      s5g=(/0.19,0.15,0.13666,0.10,0.15,0.1429013,0.14334,0.13579784,0.11287035,0.10,0.10,
     + 0.170,0.1583,0.1500,0.1500,0.1422,0.14009,0.10/)
      s6= (/0.24,0.20,0.23,0.30,0.23,0.253,0.37,0.277,0.34101113,.25,0.25,
     +0.22,0.2083,0.20731,0.22058,0.230,0.24266,0.27075/)
      s6g= (/0.24,0.20,0.32338,0.30,0.23,0.30005046,0.27662,0.32512766,0.2953982,.25,0.25,
     +0.22,0.2083,0.20731,0.22058,0.28462,0.31054,0.27075 /)
c s7 corresponds to site class E 
      s7=(/0.29,0.20,0.25,0.55,0.2,0.319,0.38,0.416,0.53479433,0.4,0.36,
     + 0.245,0.21868,0.200,0.200,0.22925,0.28797,0.46226/)
      s7g=(/0.29,0.20,0.33671,0.55,0.2,0.29974141,0.29329,0.34473022,0.49243947,0.4,0.36,
     + 0.245,0.21868,0.200,0.200,0.27997, 0.3163, 0.46226/)
      ssig=(/0.23,0.26,0.28,0.34,0.27,0.286,0.29,0.31,0.34,0.34,0.36,
     + 0.2450,0.25377,.26244,.26686,.27585,.28833,.340/)
c corrections from GDSN at pga, .2,.3,1,and 3 s. mar 2010.
        pcor= (/-0.00298,0.,-0.00290,-0.00536,0.,-0.00225,0.,0.,0.,0.,-0.0052,
     + 0.,0.,0.,0.,0.,0.,0./)
      period = perab(iq)
c Determine gnd as fcn of dist, M, period, iclass, Vs30
      if(period.ne.0.)then
       freq= 1./period
      else
      freq= 100.
      endif
c
      if(islab.eq.0)then
      depth=20.      !needed. discuss with colleagues.
      dist=dist0
c      five=iq.eq.2
c      twop5=iq.eq.6
      xmag0=min(xmag,8.5)      !p 1706, limit to 8.5
c subduction coeffs rock?
      if(ir.gt.2)then
       gnd0=s1g(iq)
       else
c cascadia case: do not change.
       gnd0=s1(iq)
      endif
c g-term differs for inslab and interface.
          g= 10.**(1.2-0.18*xmag0)
          if(ir.gt.2)then
c modified Global coeffs for f between 2 and 5 hz
          co2=s2g(iq);co3=s3g(iq);co4=s4g(iq)
          co5=s5g(iq);co6=s6g(iq);co7=s7g(iq)
          else
c  Cascadia prediction not modified 10-2008.
          co2=s2(iq);co3=s3(iq);co4=s4(iq)
          co5=s5(iq);co6=s6(iq);co7=s7(iq)
          endif
      r1=rs1;r2=rs2;r3=rs3;r4=rs4
c sigma not discussed 10-2008.
      sigma=ssig(iq)*aln10
      else
c intraslab. do not use 10-10-2008. not set up...
      stop ' 10/13: inslab not usable in SUBDUCTION CODE'
      endif
          sigmaf= 1.0/sigma/sqrt2
          delta= 0.00724*(10.**(0.507*xmag0))
          gnd=gnd0+co2*xmag0
          dist2= sqrt(dist*dist + delta*delta)
c the below if..endif is the modification of Nov 15 2013.
        if (dist.gt.400.and.period.gt.1.)then
        dist2a = sqrt((dist-400.)**2)
        gnd = gnd-0.001*dist2a        !-0.0011 is the 1-s anelastic. kicks in for LP
        endif
c End of modification Nov 15 2013.
          gnd= gnd+co3*depth+co4*dist2
     &         -g*alog10(dist2)
c--- calculate rock PGA for C site amp
c--- Cascadia subduction. rpga units are cm/s/s.
       rpga= r1+ r2*xmag0 + r3*depth+ r4*dist2- g*alog10(dist2)
       rpga= 10.**rpga
       if((rpga.le.100.).or.(freq.le.1.))then
        sl=1.
       elseif((rpga.gt.100.).and.(rpga.lt.500.).and.(freq.gt.1.).and.
     &   (freq.lt.2.)) then
       sl= 1.-(freq-1.)*(rpga-100.)/400.
       elseif((rpga.ge.500.).and.(freq.gt.1.).and.(freq.lt.2.)) then
        sl= 1.-(freq-1.)
       elseif((rpga.gt.100.).and.(rpga.lt.500.).and.(freq.ge.2.)) then
        sl= 1.-(rpga-100.)/400.
      else
cc for f>2 & pga>0.5 g, all ground condition siteamps are null compared to B.
      sl=0.
      endif
c-----
c---   Site Amp for NEHRP classes, AB03 style.
      if(vs30.gt.900.)then
      gnd=gnd-gfac
      elseif(vs30.gt.720)then
c---   take log ave of B (rock) and C site
          gnd= gnd + (sl*co5)*0.5 - gfac
          elseif(vs30.ge.380.)then
c === C site class coeff in c6
          gnd= gnd + sl*co5 -gfac
          elseif(vs30.ge.350.)then
          gnd = gnd + 0.5*sl*(co5+co6) - gfac
c --- CD boundary          
      elseif (vs30.ge.190.)then
c === D site class coeff in c6
          gnd= gnd + sl*co6 - gfac
          elseif(vs30.ge.170.)then
          gnd = gnd +0.5*sl*(co6+co7) -gfac
c--- DE boundary.
          else
          gnd= gnd + sl*co7 - gfac
      endif          
c log base 10 to base e
          gnd= gnd * aln10
          gnds=gnd; sigmas = sigma
c If requested, put the extra anelastic effect  where base e conversion has occurred
        if(p_corr.and.dist0.gt.200.)gnd=gnd+pcor(iq)*(dist0-200.0)
c      write(15,*) period, xmag, dist, exp(gnd), rpga, sl,xmag0,iq,gnd0
        if (l_ms)then
c Play it again, Sam, to compute the mean spectrum. New Oct 29 2010.
        do jp=1,npnga
c The below conditional says, compute if coeffs are available at the period with index jp. This
c will generally be true for the NGA-W relations but less so for older GMPEs with fewer periods.
c  
        j=imsp(ia,jp)
        if ( jp .eq. jms)then
        meanspec(jms)=gnds
        sigspec(jms)=sigmas
c        print *,jms,meanspec(jms),sigspec(jms)
        elseif(lmsp(ia,jp))then
      period = perab(j)
c Determine gnd as fcn of dist, M, period, iclass, Vs30
      if(period.ne.0.)then
       freq= 1./period
      else
      freq= 100.
      endif
c
      if(islab.eq.0)then
      depth=20.      !needed. discuss with colleagues.
      dist=dist0
c      five=iq.eq.2
c      twop5=iq.eq.6
      xmag0=min(xmag,8.5)      !p 1706, limit to 8.5
c subduction coeffs rock?
      if(ir.gt.2)then
       gnd0=s1g(j)
       else
c cascadia case: do not change.
       gnd0=s1(j)
      endif
c g-term differs for inslab and interface. Already calculated above.
c          g= 10.**(1.2-0.18*xmag0)
          if(ir.gt.2)then
c modified Global coeffs for f between 2 and 5 hz
          co2=s2g(j);co3=s3g(j);co4=s4g(j)
          co5=s5g(j);co6=s6g(j);co7=s7g(j)
          else
c  Cascadia prediction not modified 10-2008.
          co2=s2(j);co3=s3(j);co4=s4(j)
          co5=s5(j);co6=s6(j);co7=s7(j)
          endif
c      r1=rs1;r2=rs2;r3=rs3;r4=rs4
c sigma not discussed 10-2008.
      sigma=ssig(j)*aln10
      else
c intraslab. do not use 10-10-2008. not set up...
      stop ' 10/08: inslab needs additional coeff. for 2.5 hz'
c          g= 10.**(0.301-0.01*xmag0)
c      xmag0=min(xmag,8.0)      !p 1706, limit to 8      
c      depth=50.      !vary depth ?
c      dsq=depth*depth
c      dist= sqrt(dist0*dist0+ dsq)
c      if(ir.eq.1)then
c      gnd0=c1(j)
c      elseif(ir.gt.2)then
cc      gnd0=c1w(j)
c      endif
c      co2=c2(j);co3=c3(j);co4=c4(j);
c      co5=c5(j);co6=c6(j)
c      r1=rc1;r2=rc2;r3=rc3;r4=rc4
c      sigma=sig(j)
      endif
          gnd=gnd0+co2*xmag0
          gnd= gnd+co3*depth+co4*dist2
     &         -g*alog10(dist2)
c--- Already calculated rock PGA for C site amp
c--- Cascadia subduction. rpga units are cm/s/s.
c sl factor has already been calculated. Not redone here.
c---   Site Amp for NEHRP classes, AB03 style.
      if(vs30.gt.900.)then
      gnd=gnd-gfac
      elseif(vs30.gt.720)then
c---   take log ave of B (rock) and C site
          gnd= gnd + (sl*co5)*0.5 - gfac
          elseif(vs30.ge.380.)then
c === C site class coeff in c6
          gnd= gnd + sl*co5 -gfac
          elseif(vs30.ge.350.)then
          gnd = gnd + 0.5*sl*(co5+co6) - gfac
c --- CD boundary          
      elseif (vs30.ge.190.)then
c === D site class coeff in c6
          gnd= gnd + sl*co6 - gfac
          elseif(vs30.ge.170.)then
          gnd = gnd +0.5*sl*(co6+co7) -gfac
c--- DE boundary.
          else
          gnd= gnd + sl*co7 - gfac
      endif          
c log base 10 to base e
          gnd= gnd * aln10
          meanspec(jp)=gnd
          sigspec(jp)=sigma
          endif        !if (lmsp)
          enddo        !npnga periods.
          endif        !if(l_ms)
        
      return
      end subroutine getABsub


      subroutine Gregor06(ip,mag,rupDist,vs30,gnd,sigmaf)
c from Gregor's g06.F emailed jan 2 2007.
c Note: Gregor says that an implicit 25-km depth-of-slab is embedded in this
c model. How to adjust to 20 km open question for now. Jan 8 2007. SHarmsen.
C     This subroutine will compute the median ground motion and sigma
C     values from the Gregor et al. (2006) Cascadia Subduction Zone
C     Megathrust earthquake attenuation relationship. This
C     relationship is an update to the Gregor et al. (2002) relationship.
C     This update allows the user to input a given Vs30m value
C     and the ground motion will be computed using the site response
C     model of Walling and Abrahamson (2005). 
C
c     Input Parameters:
c      ip = period index, from period set in Period1
C     mag= moment  Magnitude
C     rupDist =rcd =  Rupture Distance (km)
C     vs30=  Vs30m (m/sec)
c      Output :
C     gnd = log(median motion) (g)
C      sigmaf = sigmafactor
      parameter (np=25,sqrt2=1.4142136,exp2p8=16.44464677)
      real mag, rupDist, gnd, pgarock, sigmaf, vs30                                  
      real n/1.18/, c/1.88/  
      real, dimension (25) :: c1, c2, c3, c4, c5, c6,
     1 period1, sig, b_soil, theta10, vref, sampadj

      Period1= (/ 0.01, 0.02, 0.025, 0.03, 0.04, 0.05, 0.056, 0.063, 0.071428571,
     1               0.083, 0.1, 0.125, 0.143, 0.167, 0.2, 0.25, 0.333, 0.4, 
     2               0.5, 0.769, 1.0, 1.667, 2.0, 2.5, 5.0 /)

      c1  =   (/ 9.30979, 9.41101, 9.55386, 9.93963, 9.89636, 10.43056, 
     1              10.96241, 11.18537, 11.6955, 11.97711, 13.45255, 13.31905,
     2              14.11414, 14.24239, 13.70461, 14.01509, 12.17016, 11.26958,
     3              10.7543, 7.75352, 6.00913, 4.08439, 5.54081, 4.79477, 5.13609/)

      c2  =   (/ -0.72025, -0.72744, -0.73843, -0.76664, -0.75422,
     1               -0.77911, -0.80698, -0.81486, -0.83129, -0.82809,
     2               -0.93278, -0.89503, -0.95795, -0.94895, -0.91753,
     3               -0.9636,  -0.81309, -0.73763, -0.73886, -0.51697,
     4               -0.38538, -0.30129, -0.50319, -0.46478, -0.61463 /)

      c5  =   (/ 2.8, 2.8, 2.8, 2.8, 2.7, 2.7, 2.8, 2.8, 2.9, 3.0,
     1               3.2, 3.3, 3.4, 3.5, 3.5, 3.5, 3.5, 3.5, 3.4, 3.2,
     2               2.9, 2.7, 2.6, 2.4, 2.1 /)

      c3  =   (/ -3.10553, -3.12269, -3.14878, -3.21488, -3.21742, 
     1               -3.23785, -3.31652, -3.34308, -3.44266, -3.47499,
     2               -3.61826, -3.64879, -3.64193, -3.70331, -3.65089,
     3               -3.6063,  -3.53255, -3.51567, -3.31553, -3.12098,
     4               -2.75873, -2.49346, -2.47444, -2.16295, -1.95111/)

      c4  =   (/ 0.23705, 0.2383, 0.2404, 0.2453, 0.24439, 0.24094,
     1               0.24457, 0.24497, 0.24982, 0.24837, 0.25549, 0.2567,
     2               0.25266, 0.2571, 0.257, 0.25674, 0.2595, 0.26323,
     3               0.25162, 0.25014, 0.22003, 0.20614, 0.20765,
     4               0.177, 0.16358 /)

      c6  =  (/ 0.03739, 0.03731, 0.03715, 0.03691, 0.03598, 0.03894, 
     1               0.03911, 0.03914, 0.03752, 0.03751, 0.04228, 0.03879,
     2               0.04621, 0.04379, 0.04239, 0.04843, 0.03491, 0.02456,
     3               0.03481, 0.01789, 0.02998, 0.02732, 0.0438, 0.05213, 0.07373/)

      sig  =  (/ 0.7028, 0.7062, 0.7135, 0.7221, 0.7290, 0.7422, 
     1               0.7516, 0.7610, 0.7655, 0.8060, 0.8240, 0.8690,
     2               0.8536, 0.8374, 0.8240, 0.8093, 0.7914, 0.7962,
     3               0.7237, 0.7871, 0.7095, 0.6570, 0.5959, 0.6552, 0.7897/)
      b_soil= (/ -1.186, -1.219, -1.248718345, -1.273, -1.308, -1.346,
     1               -1.380937866, -1.417248955, -1.455958581, -1.534878688,
     2               -1.624, -1.792954292, -1.894815052, -2.026908441, -2.188,
     3               -2.378282827, -2.568423866, -2.657, -2.669, -2.362214398,
     4               -1.955, -0.758605507, -0.299, -0.134448426, 0.0 /)

      Vref =  (/ 865.1, 865.1, 888.5995058, 907.8, 994.5, 1053.5, 
     1              1062.499994, 1071.853731, 1081.825331, 1058.230721, 
     2              1032.5, 947.2523784, 895.8574867, 829.3099131, 748.2,
     3               655.3244676, 556.5917476, 503.0, 456.6, 409.5868861,
     4               400.0, 400.0, 400.0, 400.0, 400.0 /)

      Theta10= (/ 0.9255, 0.9647, 0.997445213, 1.0242, 1.057, 1.1022,
     1               1.147199972, 1.193968654, 1.243826653, 1.335426432,
     2               1.4285, 1.62722767, 1.747038486, 1.897524986,
     3               2.0788, 2.278848596, 2.463524828, 2.5514, 2.528,
     4               2.075643639, 1.506, -0.00763144, -0.5703, 
     5              -0.742391228, -0.7993 /)
      sampadj= (/ -0.118149083, -0.118084273, -0.105904522, -0.09610953,
     1               -0.05344769, -0.025393821, -0.021093876, -0.016729322, 
     2               -0.01219193, -0.022721662, -0.035306933, -0.077444948,
     3               -0.10477571, -0.144077725, -0.198421599, -0.277997736, 
     4               -0.391539256, -0.462139859, -0.552012084, -0.709610012, 
     5               -0.817439183, -0.921429142, -0.942183118, -0.919647528, 
     6               -0.815806142 /)

      i = ip
         gnd = c1(i) + c2(i)*mag + (c3(i) + c4(i)*mag)*alog(rupdist+exp(c5(i))) +
     1            c6(i)*(mag-10.0)**3

C     Now apply Walling and Abrahamson site amplification model.
c     Site response
c     Compute PGARock (i.e., Vs=1100 m/sec)
         pgaRock =  9.28996 -0.71918*mag + (-3.10146 + 0.23683*mag)*alog(rupdist+exp2p8) +
     1              0.03741*(mag-10.0)**3
         pgarock = exp(pgarock)
          sigmaf= 1./(sig(i)*sqrt2)

c below line was commented out because not used.
c         sigrock = 0.7038
         if (vs30 .lt. vref(i)) then
             soilamp = theta10(i)*alog(vs30/vref(i)) - 1.0*b_soil(i)*alog(c+pgaRock) 
     1                 + b_soil(i)*alog(pgaRock+c*((vs30/vref(i))**n) )
         else
             soilamp = (theta10(i) + b_soil(i)*n) * alog(vs30/vref(i))
         endif

C     Now correct for reference Vs30m = 1100 m/sec
         soilamp = soilamp - sampadj(i)
         gnd = gnd + soilamp
      return
      end  subroutine Gregor06                                                                     
                                                   
      subroutine zhao(ip,xmag,dist,gnds,sigmaf,ivs, islab,hslab)
c predicted interface g.m. followed Advisory Panel meeting & suggs. Apr 2007
c compute median and sigmaf for the Zhao model with Somerville correction
c
c      input: xmag = source M. Considerably lower than Gregor. May 21 Steve Harmsen.
c            dist = distance to slab (r_cd, km)
c            ivs = Vs30 indicator, 1 for B or vs30 > 600 m/s, 2 for C, 3 for D
c            ip = period indicator, corresponding to per below. ip=1 is per(1)=0.0=PGA, etc.
c            islab=0 for interface events, =1 for intraplate events. Islab = -1 for crustal (option not used 2007)
c            hslab= depth (km) of slab events (50 km for PacNW). (only applies to interface
c            or DEEP events, with islab =1 (not used here))
c            If islab =0, an interface term si lowers median for T >= 0.5 s.
c
c added 3 periods jan 4 2011 to remove gaps when computing the mean spectrum 
c (put in slots 2 3 and 4 respectively)
      parameter (nper=25,hi=20.,sqrt2=1.41421356,xmagc=6.3,gcor=6.88806)
      common/ms/ia,npnga,jms,l_ms,lmsp,imsp,meanspec,sigspec
        logical p_corr,l_ms,lmsp(0:8,21),lvs
      integer imsp(8,21)
      real, dimension (21) :: meanspec,sigspec
      real, dimension (nper):: a,b,c,d,e,c1,c2,c3,qi,wi,ps,
     & qs,ws,si,sr,ss,ssl,sig,tau,tauI,tauS,per,sigt
      real afac,dist,hfac,site,sigma,sigmaf
c coefficients from Zhao source code not from BSSA tables: more precision in his software.      
c added 0.75s nov 2008. 0.75s (1.33 Hz) is an in-demand T.
        per = (/0.0,0.02,0.03,0.075,0.05,0.10,0.15,0.20,0.25,0.30,0.40,
     &     0.50,0.60,0.70,0.75,0.80,0.90,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0/)
        a=(/1.101,1.0902331, 1.0839348, 1.1005684,1.076,1.118,1.134,1.147,1.149,1.163,
     &       1.200,1.250,1.293,1.336,1.360956,1.386,1.433,1.479,1.551,1.621,
     &       1.694,1.748,1.759,1.826,1.825/)
        b=(/-0.00564,-6.100824E-3, -6.3703884E-3, -7.3885563E-3,-0.00671,-0.00787,-0.00722,-0.00659,
     & -0.00590,-0.00520,-0.00422,-0.00338,-0.00282,-0.00258,-2.50E-3,-0.00242,
     & -0.00232,-0.00220,-0.00207,-0.00224,-0.00201,-0.00187,-0.00147,
     & -0.00195,-0.00237/)
        c=(/0.0055,6.361353E-3, 6.865212E-3, 8.377444E-3,0.0075,0.0090,0.0100,0.0120,0.0140,
     &     0.0150,0.0060,0.0060,0.0030,0.0025,2.3797892E-3,0.0022,0.0020,0.0020,
     &     0.0020,0.0020,0.0025,0.0028,0.0032,0.0040,0.0050/)
        d=(/1.07967,1.0711297, 1.0661339, 1.0732356,1.05984,1.08274,1.05292,1.01360,
     &   0.96638,0.93427,0.95880,1.00779,1.08773,1.08384,1.0826032,1.08849,
     &   1.10920,1.11474,1.08295,1.09117,1.05492,1.05191,1.02452,
     &   1.04356,1.06518/)
        e=(/0.01412,0.014339645, 0.01446813, 0.014396015,0.01463,0.01423,0.01509,0.01462,
     &   0.01459,0.01458,0.01257,0.01114,0.01019,0.00979,9.560258E-3,0.00944,
     &   0.00972,0.01005,0.01003,0.00928,0.00833,0.00776,0.00644,
     &   0.00590,0.00510/)
        sr=(/0.2509,0.2510723, 0.25117304, 0.24486542,0.2513,0.2403,0.2506,0.2601,0.2690,
     &      0.2590,0.2479,0.2470,0.2326,0.2200,2.9499817E-3,0.2321,0.2196,0.2107,
     &      0.2510,0.2483,0.2631,0.2620,0.3066,0.3529,0.2485/)
        si=(/0.0000,0.,0.,0.,0.0000,0.0000,0.0000,0.0000,0.0000,
     &       0.0000,-0.0412,-0.0528,-0.1034,-0.1460,-0.15377375,-0.1638,-0.2062,
     &      -0.2393,-0.2557,-0.3065,-0.3214,-0.3366,-0.3306,-0.3898,
     &      -0.4978/)
c ss() corrected Mar 1, 2011
        ss=(/2.607,2.6746163,2.7141692,2.4083425,
     + 2.764,2.156,2.161,1.901,1.814,2.181,2.432,2.629,2.702,2.654,2.55,
     + 2.48,2.332,2.233,2.029,1.589,0.966,0.789,1.037,0.561,0.225/)
c sigt = Table 5 total sigma, from intra- and inter-event sigmas sqrt(SS) without the mag-cor term
c       sigt = (/0.723,0.849,0.811,0.77,0.76,0.775,0.779,0.787,0.776,0.751,0.745/)
        sig=(/0.6039,0.6194044, 0.6284739, 0.6713125,0.6399,0.6936,0.7017,0.6917,
     &  0.6823,0.6696,0.6589,0.6530,0.6527,0.6516,0.6483102,0.6467,0.6525,
     &  0.6570,0.6601,0.6640,0.6694,0.6706,0.6671,0.6468,0.6431/)
c interevnt sig tau from table 5,  Somerville says to use the mag-cor term, but not use the lower sigma_total.
        tau=(/0.3976,0.41745418, 0.42906814, 0.47095924,0.4437,0.4903,0.4603,0.4233,0.3908,
     &       0.3790,0.3897,0.3890,0.4014,0.4079,0.41473114,0.4183,0.4106,0.4101,
     &       0.4021,0.4076,0.4138,0.4108,0.3961,0.3821,0.3766/)
c intraevent tauI & tauS with the mag-cor term from table 6. see comment p 910,
c col II, Zhao et al.
        tauI=(/0.308,0.3230737, 0.3318912, 0.37809774,0.343,0.403,0.367,0.328,0.289,0.280,
     &    0.271,0.277,0.296,0.313,0.32331293,0.329,0.324,0.328,0.339,0.352,0.360,
     &    0.356,0.338,0.307,0.272/)
        tauS=(/0.321,0.34554857, 0.35990855, 0.4025684,0.378,0.420,0.372,0.324,0.294,0.284,
     &    0.278,0.272,0.285,0.290,0.29615661,0.299,0.289,0.286,0.277,0.282,0.300,
     &    0.292,0.274,0.281,0.296/)
c c1 = rock site term, for NEHRP A and B, vs30 > 600 m/s. SCI from Zhao table 2 p 901 bssa june 2006
c c2 = stiff soil term for NEHRP C (SCII)
c c3 = soft soil term for NEHRP D  (SCIII)
        c1=(/1.1111,1.3580499, 1.5025064, 1.9046799,1.6845,2.0609,1.9165,1.6688,1.4683,
     &         1.1720,0.6548,0.0713,-0.4288,-0.8656,-1.099344,-1.3250,-1.7322,
     &      -2.1522,-2.9226,-3.5476,-4.4102,-5.0492,-5.4307,-6.1813,
     &      -6.3471/)
        c2=(/1.3440,1.5373738, 1.6504902, 1.9928231,1.7930,2.1346,2.1680,2.0854,1.9416,
     &   1.6829,1.1271,0.5149,-0.0027,-0.4493,-0.69264763,-0.9284,-1.3490,-1.7757,
     &  -2.5422,-3.1689,-4.0387,-4.6979,-5.0890,-5.8821,-6.0512/)
        c3=(/1.3548,1.5238836, 1.6227913, 1.9133539,1.7474,2.0311,2.0518,2.0007,1.9407,
     &   1.8083,1.4825,0.9339,0.3936,-0.1109,-0.37189444,-0.6200,-1.0665,-1.5228,
     &  -2.3272,-2.9789,-3.8714,-4.4963,-4.8932,-5.6981,-5.8733/)
c slab event term not used for subduction
        ssl=(/-0.5284,-0.5380041, -0.54362213, -0.4743039,-0.5507,-0.4201,-0.4315,-0.3715,
     &     -0.3601,-0.4505,-0.5061,-0.5538,-0.5746,-0.5721,-0.5563261,-0.5397,
     &     -0.5216,-0.5094,-0.4692,-0.3787,-0.2484,-0.2215,-0.2625,
     &     -0.1689,-0.1201/)
c qi and wi are for the magnitude square term correction
        qi =(/0.0,0.0,0.,0.,0.,0.0,-0.0138,-0.0256,-0.0348,-0.0423,
     &        -0.0541,-0.0632,-0.0707,-0.0771,-0.07988431,-0.0825,-0.0874,-0.0917,
     &        -0.1009,-0.1083,-0.1202,-0.1293,-0.1368,-0.1486,-0.1578/)
        wi =(/0.0,0.,0.,0.,0.0,0.0,0.0286,0.0352,0.0403,0.0445,
     &    0.0511,0.0562,0.0604,0.0639,0.065498955,0.0670,0.0697,0.0721,0.0772,
     &    0.0814,0.0880,0.0931,0.0972,0.1038,0.1090/)
c Ps, Qs, and Ws for deep intraplate eqs; see last several columns of Table 6
        ps= (/0.1392,0.14970851, 0.1558556, 0.16675879,0.1636,0.1690,0.1669,0.1631,0.1588,
     &      0.1544,0.1460,0.1381,0.1307,0.1239,0.12070625,0.1176,0.1116,0.1060,
     &      0.0933,0.0821,0.0628,0.0465,0.0322,0.0083,-0.0117/)
        qs= (/0.1584,0.17338754, 0.1821547, 0.20051204,0.1932,0.2057,0.1984,0.1856,0.1714,
     &      0.1573,0.1309,0.1078,0.0878,0.0705,0.06278851,0.0556,0.0426,0.0314,
     &      0.0093,-0.0062,-0.0235,-0.0287,-0.0261,-0.0065,0.0246/)
        ws= (/-0.0529,-0.06633711, -0.074197314, -0.08620587,-0.0841,-0.0877,-0.0773,-0.0644,
     &       -0.0515,-0.0395,-0.0183,-0.0008,0.0136,0.0254,0.030538963,0.0352,
     &        0.0432,0.0498,0.0612,0.0674,0.0692,0.0622,0.0496,0.0150,
     &       -0.0268/)
      if(ivs.eq.1)then
      site=c1(ip)
      elseif(ivs.eq.2)then
      site=c2(ip)
      else
      site=c3(ip)
      endif
      if(islab.gt.0)then
      h=min(hslab,125.)
      hfac = (h - 15.)      !is this OK?
      sigma=sqrt(sig(ip)**2+tauS(ip)**2)
      afac= ssl(ip)*alog(dist) + ss(ip)
        xmcor= ps(ip)*(xmag-6.5) + qs(ip)*(xmag-6.5)**2 + ws(ip)
        elseif(islab.eq.0)then
c Interface events use si term. No others include si. Campbell noticed
c that we had called this s1 until late feb 2008.
        afac=si(ip)
        hfac=(hi-15.)
        xmcor = qi(ip)*(xmag-xmagc)**2 + wi(ip)
        h = hi
      sigma=sqrt(sig(ip)**2+tau(ip)**2)
c Frankel email may 22 2007: use sigt from table 5. Not the reduced-tau sigma
c associated with mag correction seen in table 6. Zhao says "truth" is somewhere
c in between.
      else
      h= 5.0    !depth of crustal source
      afac=0.0
      hfac=0.      !dont use the correction term for crustal events.
      sigma=sigt(ip)
c crustal events: haven't coded up the xmcor term. use interface for now.
        xmcor = qi(ip)*(xmag-xmagc)**2 + wi(ip)
      endif
      r= dist+ c(ip)*exp(d(ip)*xmag)
      gnd= a(ip)*xmag + b(ip)*dist - alog(r) + site
      gnd= gnd + e(ip)*hfac + afac
c--- magnitude square term correction
      gnd= gnd + xmcor
c cm/s_sq to g
      gnds= gnd - gcor
      sigmaf = 1./sigma/sqrt2
      sigmas=sigma
      if(l_ms)then
c play it again Sam for the other periods that define the mean spectrum.
        do jp=1,npnga
c The below conditional says, compute if coeffs are available at the period with index jp. This
c will generally be true for the NGA-W relations but less so for older GMPEs with fewer periods.
c  
        j=imsp(ia,jp)
        if ( jp .eq. jms)then
        meanspec(jms)=gnds
        sigspec(jms)=sigmas
c        print *,jms,meanspec(jms),sigspec(jms)
        elseif(lmsp(ia,jp))then
      if(ivs.eq.1)then
      site=c1(j)
      elseif(ivs.eq.2)then
      site=c2(j)
      else
      site=c3(j)
      endif
      if(islab.gt.0)then
      h=min(hslab,125.)
      hfac = (h - 15.)      !is this OK?
      sigma=sqrt(sig(j)**2+tauS(j)**2)
      afac= ssl(j)*alog(dist) + ss(j)
        xmcor= ps(j)*(xmag-6.5) + qs(j)*(xmag-6.5)**2 + ws(j)
        elseif(islab.eq.0)then
c Interface events use si term. No others include si. Campbell noticed
c that we had called this s1 until late feb 2008.
        afac=si(j)
c        hfac=(hi-15.)
        xmcor = qi(j)*(xmag-xmagc)**2 + wi(j)
        h = hi
      sigma=sqrt(sig(j)**2+tau(j)**2)
c Frankel email may 22 2007: use sigt from table 5. Not the reduced-tau sigma
c associated with mag correction seen in table 6. Zhao says "truth" is somewhere
c in between.
      else
      h= 5.0    !depth of crustal source
      afac=0.0
      hfac=0.      !dont use the correction term for crustal events.
      sigma=sigt(j)
c crustal events: haven't coded up the xmcor term. use interface for now.
        xmcor = qi(j)*(xmag-xmagc)**2 + wi(j)
      endif
      r= dist+ c(j)*exp(d(j)*xmag)
      gnd= a(j)*xmag + b(j)*dist - alog(r) + site
      gnd= gnd + e(j)*hfac + afac
c--- magnitude square term correction
      gnd= gnd + xmcor
c cm/s_sq to g
      meanspec(jp)= gnd - gcor
        sigspec(jp)=sigma
        endif        !a doable period
        enddo        !loop thru the periods.
        endif        !play it again
      return
      end subroutine zhao

      subroutine kanno(ip,gnd,sigmaf,amag,rrup,vsfac)
c....................................................................
c  Kanno et al. regression relation, 2006 for shallow earthquakes
c  written by Yuehua Zeng, USGS. Mods Steve Harmsen Nov 7 2006.
c Kanno has a PGV model but this is not included here.
c
c  input: ip  : period index, use 0 for pga.
c         amag : magnitude
c         rrup : closest fault distance
c         vsfac   : log 10 of site S-velocity in m/s
c
c  output: gnd   : ground motion spectral values (ln) (g). to be consistent
c
c          sigmaf : 1/sigma/sqrt2 errors in ln. the sigma factor.
c....................................................................
      parameter (gfac=6.8875526,sfac=2.3025851,sqrt2=1.414213562)
c conversion factors gfac cm/s/s to g, sfac log10 to ln.
      real, dimension (37) :: T,a1,b1,c1,d1,e1,p,q
        real sigma,sigmaf,gnd,vsfac
      T =(/0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.15,0.17,0.20,0.22,
     +0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,
     +1.30,1.50,1.70,2.00,2.20,2.50,3.00,3.50,4.00,4.50,5.00/)
      a1=(/0.54,0.54,0.53,0.52,0.52,0.52,0.50,0.51,0.51,0.52,0.53,0.54,0.54,
     +0.54,0.56,0.56,0.58,0.59,0.59,0.62,0.63,0.65,0.68,0.71,0.72,0.73,
     +0.74,0.77,0.79,0.80,0.82,0.84,0.86,0.90,0.92,0.94,0.92/)
      b1=(/-0.0035,-0.0037,-0.0039,-0.0040,-0.0041,-0.0041,-0.0040,-0.0040,
     +-0.0039,-0.0038,-0.0037,-0.0034,-0.0032,-0.0029,-0.0026,-0.0024,
     +-0.0021,-0.0019,-0.0016,-0.0014,-0.0012,-0.0011,-0.0009,-0.0009,
     +-0.0007,-0.0006,-0.0006,-0.0005,-0.0005,-0.0004,-0.0004,-0.0003,
     +-0.0002,-0.0003,-0.0005,-0.0007,-0.0004/)
      c1=(/0.48,0.57,0.67,0.75,0.80,0.85,0.96,0.93,0.91,0.89,0.84,0.76,0.73,
     +0.66,0.51,0.42,0.26,0.13,0.04,-0.22,-0.37,-0.54,-0.80,-1.04,-1.19,
     +-1.32,-1.44,-1.70,-1.89,-2.08,-2.24,-2.46,-2.72,-2.99,-3.21,-3.39,
     +-3.35/)
      d1=(/0.0061,0.0065,0.0066,0.0069,0.0071,0.0073,0.0061,0.0062,0.0062,
     +0.0060,0.0056,0.0053,0.0048,0.0044,0.0039,0.0036,0.0033,0.0030,
     +0.0022,0.0025,0.0022,0.0020,0.0019,0.0021,0.0018,0.0014,0.0014,
     +0.0017,0.0019,0.0020,0.0022,0.0023,0.0021,0.0032,0.0045,0.0064,
     +0.0030/)
      e1=(/0.37,0.38,0.38,0.39,0.40,0.40,0.40,0.40,0.40,0.41,0.41,0.40,0.40,
     +0.40,0.39,0.40,0.40,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,0.41,
     +0.41,0.40,0.39,0.39,0.38,0.38,0.38,0.37,0.38,0.38,0.38/)
      p= (/-0.32,-0.26,-0.24,-0.26,-0.29,-0.32,-0.35,-0.39,-0.43,-0.53,-0.61,
     + -0.68,-0.72,-0.75,-0.80,-0.85,-0.87,-0.89,-0.91,-0.92,-0.96,-0.98,
     + -0.97,-0.93,-0.92,-0.91,-0.88,-0.85,-0.83,-0.78,-0.76,-0.72,-0.68,
     + -0.66,-0.62,-0.60,-0.59/)
      q= (/0.80,0.65,0.60,0.64,0.72,0.78,0.84,0.94,1.04,1.28,1.47,1.65,1.74,
     + 1.82,1.96,2.09,2.13,2.18,2.25,2.30,2.41,2.46,2.44,2.32,2.30,2.26,
     + 2.20,2.12,2.06,1.92,1.88,1.80,1.70,1.64,1.54,1.50,1.46/)
c
c use d1 for "shallow event" category. eqn 5 of Kanno et al.
c All of Kanno regression's shallow-event data outside of Japan are crustal. No subduction.
c site term is p*vsfac+q. The factor (10*term) is generally < 1 for Vs30 > 300 m/s. 
c Vs30 = 300 m/s seems to be
c a Japan-data mean value. This site function is suspect. What about nonlin soil response?
      if(ip.eq.0)then
       gnd=0.56*amag-0.0031*rrup-alog10(rrup+0.0055*10**(0.5*amag))
     +        +0.26-0.55*vsfac+1.35
        sigma=0.37*sfac
        else
       gnd=a1(ip)*amag+b1(ip)*rrup-alog10(rrup+d1(ip)*10**(0.5*amag))
     +        +c1(ip)+p(ip)*vsfac+q(ip)
       sigma=e1(ip)*sfac
       endif
       gnd=sfac*gnd-gfac
       sigmaf=1.0/sigma/sqrt2
      return
      end subroutine kanno

      real function basiteamp(pganl,vs30,vs30r,j)
c inputs: pganl (g)
c          j : index of period in per below is j+1
c        vs30: vs30 of site
c        vs30r: a reference vs30, usually one value for soil and another for rock
c implicit input:
c         period: spectral period (scalar)  Equals per(j+1)
c Output: basiteamp = log(AMP at vs30)-log(AMP at vs30r)
        parameter (np=23)        !23 periods apr 07. include 0.01 to 10 s 
        parameter (pi=3.14159265,sqrt2=1.414213562,vref=760.)
        parameter (dx=1.098612289,dxsq=1.206948961,dxcube=1.325968960,plfac=-0.510825624)
c dx = ln(a2/a1), made a param. used in a smoothed nonlin calculation sept 2006.
c plfac = ln(pga_low/0.1)        This never changes so shouldnt be calculated.
      common/ipindx/iq,ipgeom,lvs
       logical lvs, compute/.true./
        real e1nl/-0.53804/,e2nl/-0.50350/,e3nl/-0.75472/,e4nl/-0.50970/
        real e5nl/0.28805/,e6nl/-0.10164/,e7nl/0.0/
        real c1nl/-0.66050/,c2nl/0.11970/,c3nl/-0.011510/,hnl/1.35/,b1nl/0./
      real b2nl/0./,pga_low/0.06/,mhnl/6.75/,mrefnl/4.5/,rrefnl/1.0/
      real  pganl, pganlm,pganlmec,site,siter/0.0/, site0, site0r,c,cr,d,dr
      real  a1/ 0.030/,a2/ 0.090/,a2fac/0.405465108/
        real, dimension(np):: per,e1,e2,e3,e4,e5,e6,e7,e8
     + ,mh,c1,c2,c3,c4,mref,rref,h,blin,b1,b2,v1,v2
c several quantities may be computed once, and used many times. These are saved below.
       save bnl, bnlr, c, cr, d, dr, dy, dyr, site0, site0r
        save compute
c e2,e3,e4 are the mech-dependent set. e1 is a mech-unspecified value.
c Notation change from subroutine version of 12/05 and earlier.
c array constructors for f95 Linux
c from ba_02apr07_usnr.xls coef file
      per= (/-1.000, 0.000, 0.010, 0.020, 0.030, 0.050, 0.075, 0.100,
     + 0.150, 0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000,
     + 1.500, 2.000, 3.000, 4.000, 5.000, 7.500,10.0/)
       e1= (/ 5.00121,-0.53804,-0.52883,-0.52192,-0.45285,-0.28476,
     1  0.00767, 0.20109, 0.46128, 0.57180, 0.51884, 0.43825, 0.39220, 0.18957,-0.21338,
     1 -0.46896,-0.86271,-1.22652,-1.82979,-2.24656,-1.28408,-1.43145,-2.15446/)
       e2= (/ 5.04727,-0.50350,-0.49429,-0.48508,-0.41831,-0.25022,
     1  0.04912, 0.23102, 0.48661, 0.59253, 0.53496, 0.44516, 0.40602, 0.19878,-0.19496,
     1 -0.43443,-0.79593,-1.15514,-1.74690,-2.15906,-1.21270,-1.31632,-2.16137/)
c Editorial comment I used the e3(10s)=e3(7.5)+e2(10s)-e2(7.5s) for 10s normal, because BA
c report e3(10s) as 0.0 and this gives a large motion compared to others in nhbd.
c also reported this to Dave Boore in email for his advice. SHarmsen Oct 3 2007. Gail suggests
c that ratio for normal might equal ratio for unspecified or ratio for SS. No consensus. Using
c Gail's sugg. This choice of e3(10s), -2.66, is very low. Normal median is probably too low.
       e3= (/ 4.63188,-0.75472,-0.74551,-0.73906,-0.66722,-0.48462,
     1 -0.20578, 0.03058, 0.30185, 0.40860, 0.33880, 0.25356, 0.21398, 0.00967,-0.49176,
     1 -0.78465,-1.20902,-1.57697,-2.22584,-2.58228,-1.50904,-1.81022, -2.66/)
       e4= (/ 5.08210,-0.50970,-0.49966,-0.48895,-0.42229,-0.26092,
     1  0.02706, 0.22193, 0.49328, 0.61472, 0.57747, 0.51990, 0.46080, 0.26337,-0.10813,
     1 -0.39330,-0.88085,-1.27669,-1.91814,-2.38168,-1.41093,-1.59217,-2.14635/)
       e5= (/ 0.18322, 0.28805, 0.28897, 0.25144, 0.17976, 0.06369,
     1  0.01170, 0.04697, 0.17990, 0.52729, 0.60880, 0.64472, 0.78610, 0.76837, 0.75179,
     1  0.67880, 0.70689, 0.77989, 0.77966, 1.24961, 0.14271, 0.52407, 0.40387/)
       e6= (/-0.12736,-0.10164,-0.10019,-0.11006,-0.12858,-0.15752,
     1 -0.17051,-0.15948,-0.14539,-0.12964,-0.13843,-0.15694,-0.07843,-0.09054,-0.14053,
     1 -0.18257,-0.25950,-0.29657,-0.45384,-0.35874,-0.39006,-0.37578,-0.48492/)
       e7= (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
     1  0.00000, 0.00000, 0.00000, 0.00102, 0.08607, 0.10601, 0.02262, 0.00000, 0.10302,
     1  0.05393, 0.19082, 0.29888, 0.67466, 0.79508, 0.00000, 0.00000, 0.00000/)
       e8= (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
     1  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
     1  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/)
       mh= (/ 8.50000, 6.75000, 6.75000, 6.75000, 6.75000, 6.75000,
     1  6.75000, 6.75000, 6.75000, 6.75000, 6.75000, 6.75000, 6.75000, 6.75000, 6.75000,
     1  6.75000, 6.75000, 6.75000, 6.75000, 6.75000, 8.50000, 8.50000, 8.50000/)
       c1= (/-0.873700,-0.660500,-0.662200,-0.666000,-0.690100,-0.717000,
     1 -0.720500,-0.708100,-0.696100,-0.583000,-0.572600,-0.554300,-0.644300,-0.691400,-0.740800,
     1 -0.818300,-0.830300,-0.828500,-0.784400,-0.685400,-0.509600,-0.372400,-0.098240/)
       c2= (/ 0.100600, 0.119700, 0.120000, 0.122800, 0.128300, 0.131700,
     1  0.123700, 0.111700, 0.098840, 0.042730, 0.029770, 0.019550, 0.043940, 0.060800, 0.075180,
     1  0.102700, 0.097930, 0.094320, 0.072820, 0.037580,-0.023910,-0.065680,-0.138000/)
       c3= (/-0.003340,-0.011510,-0.011510,-0.011510,-0.011510,-0.011510,
     1 -0.011510,-0.011510,-0.011130,-0.009520,-0.008370,-0.007500,-0.006260,-0.005400,-0.004090,
     1 -0.003340,-0.002550,-0.002170,-0.001910,-0.001910,-0.001910,-0.001910,-0.001910/)
       c4= (/ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
     1  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
     1  0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/)
       mref= (/4.500,4.500,4.500,4.500,4.500,4.500,
     1 4.500,4.500,4.500,4.500,4.500,4.500,4.500,4.500,4.500,
     1 4.500,4.500,4.500,4.500,4.500,4.500,4.500,4.500/)
       rref= (/1.000,1.000,1.000,1.000,1.000,1.000,
     1 1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000,
     1 1.000,1.000,1.000,1.000,1.000,1.000,1.000,1.000/)
          h= (/2.540,1.350,1.350,1.350,1.350,1.350,
     1 1.550,1.680,1.860,1.980,2.070,2.140,2.240,2.320,2.460,
     1 2.540,2.660,2.730,2.830,2.890,2.930,3.000,3.040/)
      blin = (/-0.600,-0.360,-0.360,-0.340,-0.330,-0.290,
     1 -0.230,-0.250,-0.280,-0.310,-0.390,-0.440,-0.500,-0.600,-0.690,
     1 -0.700,-0.720,-0.730,-0.740,-0.750,-0.750,-0.692,-0.650/)
       b1= (/-0.500,-0.640,-0.640,-0.630,-0.620,-0.640,
     1 -0.640,-0.600,-0.530,-0.520,-0.520,-0.520,-0.510,-0.500,-0.470,
     1 -0.440,-0.400,-0.380,-0.340,-0.310,-0.291,-0.247,-0.215/)
       b2= (/-0.060,-0.140,-0.140,-0.120,-0.110,-0.110,
     1 -0.110,-0.130,-0.180,-0.190,-0.160,-0.140,-0.100,-0.060, 0.000,
     1  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
         v1= (/ 180., 180., 180., 180., 180., 180.,
     1  180., 180., 180., 180., 180., 180., 180., 180., 180.,
     1  180., 180., 180., 180., 180., 180., 180., 180./)
         v2= (/ 300., 300., 300., 300., 300., 300.,
     1  300., 300., 300., 300., 300., 300., 300., 300., 300.,
     1  300., 300., 300., 300., 300., 300., 300., 300./)
c end of April 07 coeff. updates. Remove sigma inits. not used here.
c----  Input vs30.  
c  ip index corresponds to period. (check per correspondence with main perx)
c--- 
        iq = j+1     ! period index is an argument: 1/2011.
        if(compute)then
           if(v1(iq).lt.vs30.and.vs30.le.v2(iq))then
c some site term precalcs that are not M or d dependent
        bnl=(b1(iq)-b2(iq))*
     + alog(vs30/v2(iq))/alog(v1(iq)/v2(iq)) + b2(iq)
        elseif(v2(iq).lt.vs30.and.vs30.le.vref)then
        bnl=b2(iq)*alog(vs30/vref)/alog(v2(iq)/vref)
        elseif(vs30.le.v1(iq))then
        bnl=b1(iq)
        else
        bnl=0.0
        endif
        if(v1(iq).lt.vs30r.and.vs30r.le.v2(iq))then
c repeat site term precalcs that are not M or d dependent @ ref. vs.
        bnlr=(b1(iq)-b2(iq))*
     + alog(vs30r/v2(iq))/alog(v1(iq)/v2(iq)) + b2(iq)
        elseif(v2(iq).lt.vs30r.and.vs30r.le.vref)then
        bnlr=b2(iq)*alog(vs30r/vref)/alog(v2(iq)/vref)
        elseif(vs30r.le.v1(iq))then
        bnlr=b1(iq)
        else
        bnlr=0.0
        endif
c----- ADF added next line
        dy= bnl*a2fac 
        dyr= bnlr*a2fac 
        c=(3.*dy - bnl*dx)/dxsq
        d=(bnl*dx-2.*dy)/dxcube
        cr=(3.*dyr -bnlr*dx)/dxsq
        dr=(bnlr*dx-2.*dyr)/dxcube
        site0 = blin(iq)*alog(vs30/vref)
        site0r = blin(iq)*alog(vs30r/vref)
        compute = .false.
c        print *, bnl, bnlr, dy,dyr,site0,site0r
        endif        !first-pass compute
c Second part, nonlinear siteamp reductions below.
        if(pganl.le.a1)then
        site=site0+bnl*plfac
        siter=site0r+bnlr*plfac
        pgafac=0.
        elseif(pganl.le.a2)then
c extra lines smooth a kink in siteamp, pp 9-11 of boore sept report.
c c and d from p 10 of boore sept report. Smoothing introduces extra calcs
c in the range a1 < pganl < a2. Otherwise nonlin term same as in june-july.
c many of these terms are fixed and are defined in data or parameter statements
c Of course, if a1 and a2 change from their sept 06 values the parameters will
c also have to be redefined. (a1,a2) represents a siteamp smoothing range (units g)
        pgafac=alog(pganl/a1)
        psq=pgafac*pgafac
        site=site0+bnl*plfac + (c + d*pgafac)*psq
        siter= site0r +bnlr*plfac + (cr + dr*pgafac)*psq
        else
        pgafac=alog(pganl/0.1)
        site=site0+bnl*pgafac
        siter=site0r+bnlr*pgafac
        endif
        basiteamp = site-siter
        return
      end 

      subroutine getNAAsub(iq, Fevnt, Ffaba, xmag, R, zH, sigmaf, gm, delCi, Vs30)
c This is Abrahamson 2010 BCHydro Subduction model that was developed using the dataset of 
c acceleration response spectra of the geometrical mean 5% spectral damping. The coefficients
c  have been written for periods: pga, 0.050, 0.075,0.100, 0.150,0.200,0.250,0.300,0.400,
c 0.500, 0.600, 0.750, 1.000, 1.500, 2.000, 2.500, 3.000, 4.000, 5.000, 6.000, 7.500, 10.000.
c
c
c For values of T=0.02sec, use PGA values
c Recommended values for DelC1 is -0.5, 0.0, 0.5
c PGA1000 = Median PGA value for a Vs30 1,000 m/sec
c***************************************************************************************************
c The parameters are defined as:
c ip = period  index, global
c iq        -        Period index, internal
c Fevnt -        Indicator variable for an intraslab & interface earthquake. Intraslab = 1, interface = 0
c Ffaba -         Indicator variable for a forearc or unknown sites & backarc sites. Forearc or unknown
c                sites = 0, backarc sites = 1.
c xmag -         Magnitdue in Mw
c R        -         Rupture distance for interface events or hypocentral distance for intraslab events.  
c zH -          Hypocentral Distance (km)
c sigmaf -         sigmaf = 1/sigma/sqrt2
c gm  -         logged SA ground motion (units g)
c delCi =         Integer indicator for low middle hi branch, 1 2 or 3 resp.
c***************************************************************************************************
c Note delCi changes the median by about 57%. It affects a magnitude dep. as well.
        implicit none
        real sqrt2
        real, dimension(21):: meanspec,sigspec
        logical l_ms,lmsp(0:8,21)
        integer imsp(8,21),ia,npnga,jms
        real, dimension(22):: NAAper,b,theta1, theta2, theta6, theta7, theta8, vlin
        real theta10(22), theta11(22), theta12(22), theta13(22), theta14(22)
        real theta15(22), theta16(22)
        real delC(3), delC1, dy_sdi, rhat,sdisd,sde,sdi_ratio
        real :: fac_sde
        real xmag, R, zH, fterm
        integer Fevnt, Ffaba, delCi, c4/10/, ip, iq, j,jp
        real c1/7.8/, c/1.88/, n/1.18/
        real theta3/0.1/, theta4/0.9/, theta5/0.0/, theta9/0.4/
        real fMag,PGArock,gm, Rmax
        real fDepth,  fSite, Vs30, VsStar
        real sigma/0.772/,sigmaf
        logical sdi
        parameter (sqrt2=1.414213562)
        common/sdi/sdi,dy_sdi,fac_sde
        common/ms/ia,npnga,jms,l_ms,lmsp,imsp,meanspec,sigspec
c NAAper = period set for Norms BC Hydro subduction model. 0.00 for 0 s (PGA) to 0.02 s SA
        NAAper =(/0.00, 0.050, 0.075,0.100, 0.150,0.200,0.250,0.300,0.400,
     + 0.500, 0.600, 0.750, 1.000, 1.500, 2.000, 2.500, 3.000, 4.000, 5.000, 6.0, 7.500, 10.0/)
        vlin = (/865.1,1053.5,1085.7,1032.5,877.6,748.2,654.3,587.1,503.0,456.6,430.3,410.5,400.0,
     1            400.0, 400.0, 400.0, 400.0, 400.00, 400.0, 400.0, 400.0, 400.0/)
        b = (/-1.186, -1.346, -1.471,-1.624, -1.931,-2.188, -2.381, -2.518, -2.657, -2.669,
     1       -2.599, -2.401, -1.955, -1.025, -0.299, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/)
        theta1  = (/4.2203, 4.5371, 5.0733, 5.2892, 5.4563,5.2684, 5.0594, 4.7945, 4.4644, 4.0181, 
     1                 3.6055, 3.2174, 2.7981, 2.0123,1.4128, 0.9976, 0.6443, 0.0657, -0.4624, -0.9809,
     2                 -1.6017, -2.2937/)
        theta2  = (/-1.350, -1.400, -1.450, -1.450, -1.450, -1.400, -1.350, -1.280, -1.180, -1.080, 
     1                 -0.990, -0.910, -0.850, -0.770, -0.710, -0.670, -0.640, -0.580, -0.540, -0.500,
     2                  -0.460, -0.400/)
        theta6  = (/-0.0012, -0.0012, -0.0012, -0.0012, -0.0014, -0.0018, -0.0023, -0.0027,
     1                 -0.0035, -0.0044, -0.0050, -0.0058, -0.0062, -0.0064 , -0.0064, -0.0064,
     2                 -0.0064, -0.0064, -0.0064, -0.0064 , -0.0064, -0.0064/)
        theta7  = (/1.0988, 1.2536, 1.4175, 1.3997, 1.3582, 1.1648, 0.9940, 0.8821, 0.7046,
     1                 0.5799, 0.5021, 0.3687, 0.1746, -0.0820, -0.2821, -0.4108, -0.4466, -0.4344,
     2               -0.4368, -0.4586, -0.4433, -0.4828/)
        theta8  = (/-1.42, -1.65, -1.80, -1.80, -1.69, -1.49, -1.30, -1.18, -0.98, -0.82, -0.70, 
     1                 -0.54, -0.34, -0.05, 0.12, 0.25, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30/)
        theta10 = (/3.12,3.37,3.37,3.33, 3.25, 3.03, 2.80, 2.59, 2.20, 1.92, 1.70, 1.42, 1.10, 
     1                  0.70, 0.70, 0.70, 0.70, 0.70, 0.70,0.70, 0.70, 0.70/)
        theta11 = (/0.0130, 0.0130,0.0130,0.0130,0.0130,0.0129,0.0129,0.0128,0.0127,0.0125,
     1                    0.0124,0.0120,0.0114,0.0100,0.0085,0.0069,0.0054, 0.0027,0.0005, -0.0013,
     2                   -0.0033, -0.0060/)
        theta12 = (/0.980, 1.288, 1.483, 1.613, 1.882, 2.076, 2.248, 2.348, 2.427, 2.399, 2.273,
     1                   1.993, 1.470, 0.408, -0.401, -0.723, -0.673, -0.627, -0.596, -0.566,
     2                    -0.528, -0.504/)
        theta13 = (/-0.0135, -0.0138,-0.0142,-0.0145,-0.0153,-0.0162,-0.0172, -0.0183,-0.0206,
     1                   -0.0231, -0.0256, -0.0296, -0.0363, -0.0493, -0.0610, -0.0711, -0.0798, 
     2                   -0.0935, -0.0980 , -0.0980, -0.0980, -0.0980/) 
        theta14 = (/-0.40,-0.40,-0.40,-0.40,-0.40,-0.35,-0.31, -0.28, -0.23, -0.19, -0.16, -0.12,
     1                   -0.07, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00/)
        theta15 = (/0.9996, 1.1030, 1.2732, 1.3042, 1.2600, 1.2230, 1.1600, 1.0500, 0.800,
     1                    0.6620, 0.5800, 0.4800, 0.3300, 0.3100, 0.300 , 0.300, 0.300, 0.300, 0.300,
     2                    0.300, 0.300, 0.300/)
        theta16 = (/-1.00, -1.18, -1.36, -1.36, -1.30, -1.25, -1.17, -1.06, -0.78, -0.62, -0.50,
     1                  -0.34, -0.14, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00/) 
        delC = (/-0.5, 0.0,0.5/)
     
C        Define DelC1
        if(delCi.le.3.and.delCi.gt.0) then
                delC1 = delC(delCi)
        else 
                print *, 'Error in delC1'
        endif
C                
C        Calculate Magnitude Scaling Factor fMag(M)
        if(xmag.le.c1+delC1) then 
                fMag = theta4*(xmag-(c1+delC1))+theta13(iq)*(10-xmag)**2
        else !its greater
                fMag = theta5*(xmag-(c1+delC1))+theta13(iq)*(10-xmag)**2
        endif
C        
C         Calculate Depth Scaling Factor fDepth(zH)
        if(Fevnt > 1.or.Fevnt < 0) then 
                print *,'Error in Fevnt term'
        else        
                fDepth = theta11(iq)*(zH-60)*Fevnt
        endif
C       
C        Calculate Forearc/Backarc Scaling fFaba
        if(Ffaba .gt. 1 .or. Ffaba .lt. 0) then 
                print *,fFaba
                print *,'Error in Ffaba term'
        endif
        if(Fevnt.eq.1 .and. Ffaba.eq.1) then ! its a backarc site loc and intraplate src.
                Rmax = max(R, 85.0)
                fterm = (theta7(iq)+theta8(iq)*alog(Rmax/40))	!*Ffaba
        elseif (Ffaba.eq.1)then !its a backarc  and interface src . 
                Rmax = max(R, 100.0)
                fterm = (theta15(iq)+theta16(iq)*alog(Rmax/40))	!*Ffaba
        else        !forearc site or the site is unknown. No effect.
                fterm=0.0
        endif
C
C        Calculate Site Scaling Factor        fSite
c     compute VsStar
        if(Vs30.gt.1000.) then
                VsStar = 1000.
        else
                VsStar = Vs30
        endif
C
C Calculate PGArock
        PGArock = theta12(iq)*alog(1000./vlin(iq))+b(iq)*n*alog(1000./vlin(iq))
C        
        if(Vs30.lt.Vlin(iq)) then ! non-linear site response
                fSite = theta12(iq)*alog(VsStar/vlin(iq)) - b(iq)*alog(PGArock+c)+b(iq)*alog(PGArock+c*(VsStar/vlin(iq))**n)
        else ! linear site repsonse
                fSite = theta12(iq)*alog(VsStar/vlin(iq))+b(iq)*n*alog(VsStar/vlin(iq))
        endif        
C

C Calulate ground motion of Base Function
        gm = theta1(iq)+theta4*delC1+(theta2(iq)+theta14(iq)*Fevnt+theta3*(xmag-7.8))*
     1        alog(R+c4*exp((xmag-6)*theta9))+theta6(iq)*R+theta10(iq)*Fevnt+fMag+fDepth+
     2       fterm +fSite
        if(sdi)then
       sde=gm+fac_sde        !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)        !10 is an upper bound for rhat.
       gm = sdi_ratio(NAAper(iq),xmag,rhat,sigma,sdisd) + sde
       sigmaf= 1./sdisd/sqrt2
       else
        
       sigmaf = 1./sigma/sqrt2
       endif
       if(l_ms)then
c Play it again, Norm, for the CMS.
        do jp=1,npnga
c The below conditional says, compute if coeffs are available at the period with index jp. This
c will generally be true for the NGA-W relations but less so for older GMPEs with fewer periods.
c  
        j=imsp(ia,jp)
        if ( jp .eq. jms)then
        meanspec(jp)=gm
        sigspec(jp)=sigma
c        print *,jms,meanspec(jms),sigspec(jms)
        elseif(lmsp(ia,jp))then
c        print *,jp,j,' inside NAA '
C        Calculate Magnitude Scaling Factor fMag(M)
        if(xmag.le.c1+delC1) then 
                fMag = theta4*(xmag-(c1+delC1))+theta13(j)*(10-xmag)**2
        else !its greater
                fMag = theta5*(xmag-(c1+delC1))+theta13(j)*(10-xmag)**2
        endif
C        
C         Calculate Depth Scaling Factor fDepth(zH)
                fDepth = theta11(j)*(zH-60)*Fevnt
        if(Fevnt.eq.1 .and. Ffaba.eq.1) then ! its a backarc site loc and intraplate src.
                Rmax = max(R, 85.0)
                fterm = (theta7(j)+theta8(j)*alog(Rmax/40))	!*Ffaba
        elseif (Ffaba.eq.1)then !its a backarc  and interface src . 
                Rmax = max(R, 100.0)
                fterm = (theta15(j)+theta16(j)*alog(Rmax/40))	!*Ffaba
        else        !forearc site or the site is unknown. No effect.
                fterm=0.0
        endif
C
C        Calculate Site Scaling Factor        fSite
C
C Calculate PGArock. why is PGArock frequency dependent?
        PGArock = theta12(j)*alog(1000./vlin(j))+b(j)*n*alog(1000./vlin(j))
C        
        if(Vs30.lt.Vlin(j)) then ! non-linear site response
                fSite = theta12(j)*alog(VsStar/vlin(j)) - b(j)*alog(PGArock+c)+b(j)*alog(PGArock+c*(VsStar/vlin(j))**n)
        else ! linear site repsonse
                fSite = theta12(j)*alog(VsStar/vlin(j))+b(j)*n*alog(VsStar/vlin(j))
        endif        
C

C Calulate ground motion of Base Function
        meanspec(jp) = theta1(j)+theta4*delC1+(theta2(j)+theta14(j)*Fevnt+theta3*(xmag-7.8))*
     1        alog(R+c4*exp((xmag-6)*theta9))+theta6(j)*R+theta10(j)*Fevnt+fMag+fDepth+
     2       fterm +fSite
        sigspec(jp) = sigma        !period independent in BCHydro as we received it.
        endif
        enddo        !cms periods
        endif        !do you want cms?
       return
       end subroutine getNAAsub
        
        
        subroutine getAtkinsub(ip,rcd,M,vs30,backarc,gnds,sigmaf)
        real rcd,M,vs30,gnd,sigma,sfac,gfac,gnd0,gndp,gnd0p,pganl
c Atkinson and Macias (BSSA, 2009) for great cascadia subduction eqs M>=7.5
c        input variables:
c        ip = index corresponding to  spectral period (s)
c        iperba = index of period in the BAsiteamp fcn. In common /ipindex/
c        rcd = closest distance to rupture (km)
c        M = moment magnitude
c        vs30 = shear wave veloc top 30 m
c        forearc = logical variable .true. if site is in the forearc region
c        backarc= 1 (integer) if site is in the backarc.
c        Varying response w.r.t. this situation is not in the original formulation but may be helpful.
c c Always define a reference PGA, determine nonlinear response based on this value
c
c        Output:
c        gnds = ln(median SA or PGA) units g
c        sigmaf = 1/sqrt2/sig_lnY based on Campbell&Bozorgnia
c        coded Nov 7 2012 S Harmsen USGS harmsen@usgs.gov
c        Add CMS calcs Jun 28 2013. looks lower than others 
c
c        ipg = global ip index used with site amplification calculation=1.
        integer np
        parameter (np=15)
        parameter (gfac=6.8875526,sfac=2.3025851,sqrt2=1.4142136)
        common/sdi/sdi,dy_sdi,fac_sde
        common/ipindx/iperba,ipgeom,lvs
        common/ms/ia,npnga,jms,l_ms,lmsp,imsp,meanspec,sigspec
        real, dimension(21):: meanspec,sigspec
        logical l_ms,lmsp(0:8,21)
        integer imsp(8,21),ia,npnga,jms
        real dy_sdi,rhat, sigmaf,sdisd
        real fac_sde, xm, xmsq
c base10 to natural logs. cm/s/s to g
        logical sdi,lvs
        integer backarc,ip,ipq
        real, dimension(np):: fr,pd,c0,c1,c2,c2b,c3,c4,sig
        sig =(/0.24,0.26,0.27,0.27,0.27,0.27,0.27,0.29,0.30,0.30,
     + 0.30,0.3,0.32,0.35,0.38/)
             c0=(/5.006,5.843,5.490,4.746,4.303,4.167,3.999,3.621,3.241,3.104,2.978,
     + 2.814,2.671,2.489,2.338/)
             c1=(/-1.5573,-1.9391,-1.6257,-1.1691,-0.9322,-0.8854,
     + -0.8211,-0.7376,-0.6741,-.6585,-0.6431,-0.6108,-0.5942,-0.6412,-0.6311/)
       c2=(/-0.00034,0.0,-0.00115,-0.00212,-0.00231,-0.00211,-0.00195,
     + -0.00128,-0.00081,-0.00063,-0.00057,-0.00046,-0.00040,-0.00003,0.0/)
       c2b=(/-0.0015,-.0015,-0.00225,-0.00332,-.00331,-0.00281,-0.00225,
     + -0.00158,-0.00083,-0.00065,-0.00059,-0.00048,-0.00041,-0.000032,0.0/)
c c2b first attempt at decay in backarc region. hybridizing two Atkinson papers       
      c3=(/0.1774,0.1813,0.1736,0.1593,0.1713,0.1802,0.1870,0.2116,0.2696,0.2990,
     + 0.3258,0.3490,0.3822,0.4760,0.5357/)
      c4 = (/0.0827,0.0199,0.0261,0.0432,0.0270,0.0258,0.0271,0.0328,
     + -0.0064,-0.0074,-0.0103,-0.0299,-0.0417,-0.0629,-0.0737/)
      fr = (/99.,20.,10.,5.,3.16,2.5,2.,1.0,0.5,0.4,0.32,0.25,0.2,0.13,0.10/)
      pd = (/0.01,0.05,0.1,0.2,0.3,0.4,0.5,1.,2.,2.5,3.,4.,5.,7.7,10./)
c h term from M
        h=M**2-3.1*M-14.55        !eqn 6
        xm=M-8.0; xmsq = xm*xm
        gnd0p=c0(1)+c3(1)*xm+c4(1)*xmsq
        gnd0=c0(ip)+c3(ip)*xm+c4(ip)*xmsq
        r=sqrt(rcd**2+h**2)
        if(backarc.eq.0)then
c original GMPE form
        gnd=gnd0 + c1(ip)*alog10(r)+c2(ip)*r
        gndp=gnd0p+c1(1)*alog10(r)+c2(1)*r
        else
c backarc, steeper falloff Qterm. from Gofrani and Atkinson. Trial balloon: not using 2013. 
        gnd=gnd0 + c1(ip)*alog10(r)+c2b(ip)*r
        gndp=gnd0p + c1(1)*alog10(r)+c2b(1)*r
        endif        
        gnd = gnd*sfac - gfac
        if(vs30.ne.760.)then
        pganl = exp(gndp*sfac - gfac)
c siteamp: provisional. the above gnd corresponds to rock at B/C boundary (article, p.1569)
c she suggests BA08 siteamp. pganl is the reference rock pga 
        gnd=gnd+basiteamp(pganl,vs30,760.,iperba)
        endif
        gnds=gnd
        sigma=sig(ip)*sfac
        if(sdi)then
       sde=gnd+fac_sde        !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)        !10 is an upper bound for rhat.
       gnd = sdi_ratio(pd(ip),M,rhat,sigma,sdisd) + sde
       sigmaf= 1./sdisd/sqrt2
       else
        sigmaf=1./sigma/sqrt2
        endif
       if(l_ms)then
c Play it again, Gail, for the CMS.
        do jp=1,npnga
c The below conditional says, compute if coeffs are available at the period with index jp. This
c will generally be true for the NGA-W relations but less so for older GMPEs with fewer periods.
c  
        j=imsp(ia,jp)
        if ( jp .eq. jms)then
        meanspec(jp)=gnds
        sigspec(jp)=sigma
c        print *,jms,meanspec(jms),sigspec(jms),' a-m sub'
        elseif(lmsp(ia,jp))then
        gnd0=c0(j)+c3(j)*xm+c4(j)*xmsq
        if(backarc.eq.0)then
c original GMPE form
        gnd=gnd0 + c1(j)*alog10(r)+c2(j)*r
        else
c backarc, steeper falloff Qterm. from Gofrani and Atkinson. Trial balloon: not using 2013. 
        gnd=gnd0 + c1(j)*alog10(r)+c2b(j)*r
        endif
        gnd = gnd*sfac - gfac
        if(abs(vs30-760.).gt.12.)then
c siteamp: provisional. the above gnd corresponds to rock at B/C boundary (article, p.1569)
c she suggests BA08 siteamp. pganl is the reference rock pga 
c the jp index below is suspect but is used in getGeom. 
c May need  further investigation. SH June 25 2013.
        gnd=gnd+basiteamp(pganl,vs30,760.,jp)
        endif
        meanspec(jp)=gnd
        sigspec(jp)=sig(j)*sfac
        endif        !computable?
        enddo
c        print *,'getAtkinsub: '
c        print *,sigspec
        endif        !do you want CMS?
        return
        end subroutine getAtkinsub

        real function sdi_ratio(per,M,rhat,sige,sdi)
        real M,rhat,td,dt,sde,sdi
c inputs:
c per = spectral period per (s)
c rhat= elastic spectral displacement Sd_e/dY. (Sd_e is not input)
c M = moment magnitude
c elastic ln(S_d dispersion), sige
c 
c Outputs:
C compute tothung-cornell spectral displacement ratio, ln(S_di/S_d)
c Also compute the inelastic log st. deviation, sdi
c steve harmsen feb 14 2013.
c Use regression coeffs from table 2.
c The regression coeffs were designed for a specific GMPE Abrahmson Silva(which?)
c but we will use these with any crustal-event GMPE at least WUS...
c Possibly different coeffs may be found for subduction event GMPEs...
        real, dimension(18):: beta1,beta2,beta3,beta4,beta5,beta6, 
     + c1p,c2p,c3p,c4p,ap,bp,T
c special functions g1, g2
        g1(r,xm,b1,b2,b3,b4,b5)=(b1+b2*xm)*r+
     + (b3+b4*xm)*r*alog(r) +b5*r**2.5
         g2(r,b6)=0.37*b6*(r-0.3)
         T=(/ 0.30, 0.40, 0.50, 0.60, 0.70, 0.750, 0.80, 
     + 0.850, 0.90, 1.0, 1.100, 1.250, 1.50, 2.0, 2.50, 3.0, 4.0, 5.0/)
        beta1=(/-0.10370, 0.12260, 0.04490, 0.11160,-0.07590,-0.13070,-0.32330,-0.37980,-0.37910,
     + -0.38000,-0.33240,-0.35820,-0.42440,-0.13040,-0.22750, 0.06930,-0.16920,-0.41700/)
        beta2=(/ 0.01380,-0.01970,-0.00730,-0.02110, 0.00660, 0.01570, 0.04340, 0.05040, 
     + 0.04890, 0.05150, 0.04440, 0.04630, 0.05610, 0.01300, 0.02640,-0.01950, 0.01610, 0.05170/)
        beta3=(/-0.00590,-0.06960,-0.04890,-0.09530,-0.03500,-0.01560, 0.05680, 0.09710, 
     + 0.09130, 0.06600, 0.03780, 0.07960, 0.11250,-0.00420, 0.05370,-0.06030, 0.03480, 0.07880/)
        beta4=(/ 0.01050, 0.01600, 0.01060, 0.01780, 0.00760, 0.00380,-0.00690,-0.01170,
     + -0.01020,-0.00840,-0.00510,-0.01120,-0.01570, 0.00090,-0.00720, 0.01200,-0.00380,-0.01080/)
        beta5=(/-0.00220,-0.00080,-0.00070,-0.00030,-0.00030, 0.00000, 0.00000,-0.00040,
     + -0.00050,-0.00050, 0.00000, 0.00060, 0.00050, 0.00100, 0.00090, 0.00010, 0.00040, 0.00000/)
        beta6=(/-0.07400, 0.01300, 0.01700, 0.03700, 0.05000, 0.03000,-0.02500,-0.05000,
     + -0.05000,-0.06100,-0.05000,-0.02000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000/)
        c1p=(/ 0.00200, 0.00800, 0.00600, 0.00800, 0.00700, 0.00600, 
     + 0.00900, 0.01100, 0.01400, 0.01300, 0.011, 0.012, 0.008, 0.008, 0.008, 0.013, 0.01800, 0.03500/)
        c2p=(/-0.01000,-0.04000,-0.02800,-0.04000,-0.03300,-0.03000,-0.04700,
     + -0.05600,-0.07200,-0.06400,-0.05600,-0.06200,-0.04000,-0.04100,-0.03800,-0.06700,-0.08800,-0.1770/)
        c3p=(/ 0.09600, 0.09500, 0.07500, 0.06500, 0.05900, 0.05600, 0.04800, 
     + 0.03700, 0.04600, 0.07300, 0.07800, 0.03900, 0.02000, 0.04500, 0.05300,
     + 0.08500, 0.09100, 0.16500/)
        c4p=(/-0.07700,-0.04500,-0.03600,-0.01000,-0.01300,-0.01200, 0.01600, 
     + 0.03500, 0.04100, 0.00700,-0.00800, 0.03700, 0.03800, 0.01400,-0.00600,
     + -0.00900,-0.00500, 0.01700/)
         ap=(/0.9,0.7,0.9,1.0,1.8,1.8,1.2,1.2,1.2,1.4,1.8,1.4,1.0,1.4,1.6,
     + 1.2,1.0,0.8/)
         bp=(/3.8,3.8,3.0,2.9,5.0,4.5,2.0,1.8,1.4,3.4,5.0,2.0,2.2,3.5,5.0,
     + 2.0,4.5,5.0/)
c
c find the index corresponding to the requested period 
       i=1
       if(per.le.T(1))then
c first try: do not extrapolate, just use 0.3s coeffs for shorter periods.
         b1=beta1(1); b2=beta2(1); b3=beta3(1)
         b4=beta4(1); b5=beta5(1); b6=beta6(1)
         a=ap(1); b=bp(1); c1 = c1p(1); c2= c2p(1)
         c4=c4p(1); c3=c3p(1)
         goto 3
       endif
         ip=i+1
         dowhile(per.gt.T(ip).and.ip.lt.18)
         ip=ip+1
         enddo
       if(abs(per-T(ip)).lt.0.01)then
c no need to interpolate
        b1=beta1(ip); b2=beta2(ip); b3=beta3(ip)
        b4=beta4(ip); b5=beta5(ip); b6=beta6(ip)
        a=ap(ip); b=bp(ip); c1 = c1p(ip); c2= c2p(ip)
         c4=c4p(ip); c3=c3p(ip)
        elseif(T(ip).lt.per)then
c extrapolate: Need a rule. for now just use 5s coeffs.
        b1=beta1(ip); b2=beta2(ip); b3=beta3(ip)
        b4=beta4(ip); b5=beta5(ip); b6=beta6(ip)
        a=ap(ip); b=bp(ip); c1 = c1p(ip); c2= c2p(ip)
         c4=c4p(ip); c3=c3p(ip)
         else
c interpolate. From meeting Feb 12 2013 try linear interpolation.
        i=ip-1
        dt=(per-T(i))/(T(ip)-T(i))
        td=1.0-dt
        b1=beta1(ip)*dt+beta1(i)*td
        b2=beta2(ip)*dt+beta2(i)*td
        b3=beta3(ip)*dt+beta3(i)*td
        b4=beta4(ip)*dt+beta4(i)*td
        b5=beta5(ip)*dt+beta5(i)*td
        b6=beta6(ip)*dt+beta6(i)*td
        a = ap(ip)*dt+ap(i)*td
        b=  bp(ip)*dt+bp(i)*td
        c1=  c1p(ip)*dt+c1p(i)*td
        c2=  c2p(ip)*dt+c2p(i)*td
        c3=  c3p(ip)*dt+c3p(i)*td
        c4=  c4p(ip)*dt+c4p(i)*td
        endif
3        continue
c        print *,'b1,b2,b3,b4,b5,b6,a,b:'
c        print *,b1,b2,b3,b4,b5,b6,a,b
       if(rhat.lt.0.2)then
        sdi_ratio=0.0
        sdi = sige
        elseif(rhat.lt.0.3)then
        sdi_ratio=g1(rhat,M,b1,b2,b3,b4,b5) - g1(0.2,M,b1,b2,b3,b4,b5)
c because in this interval g2 is zero we dont see it ablove
        elseif(0.3.le. rhat .and. rhat .lt. 3.0)then
        sdi_ratio=g1(rhat,M,b1,b2,b3,b4,b5) - g1(0.2,M,b1,b2,b3,b4,b5)+
     +  g2(rhat,b6)*(M-6.5)
        elseif(3.0.le.rhat.and. rhat.le.10.)then
        sdi_ratio=g1(rhat,M,b1,b2,b3,b4,b5) - g1(0.2,M,b1,b2,b3,b4,b5)+
     +  b6*(M-6.5)
        else
        stop' no rule when rhat >100'
        endif
c compute sdi for rhat.ge. 0.2
        if(rhat.ge.0.2)then
        sdi = sige+c1+c2*rhat
c        print *,c3,c4
        if(rhat.gt.a)sdi= sdi + c3*(rhat-a)        
        if(rhat.gt.b)sdi= sdi + c4*(rhat-b)
        endif        
        return
        end     function sdi_ratio   
c------------------------------------------------------------------------------
