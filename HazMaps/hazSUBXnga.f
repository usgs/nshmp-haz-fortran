c--- program hazSUBXnga.f 
c last mod Mar 4 2009 Crouse atten for soil. Originally from hazSUBXv3.f, a Frankel-Wesson code.
c  Technical Contact: Steve Harmsen USGS. 303 273 8567. harmsen@usgs.gov
c
c  NGA in program name because some "modern" attenuation models of 2006 are present here, even
c  though these are not the product of the famous PEER NGA teams.
c Computes probabilistic hazard curves for subduction events. These can be
c characteristic (slab filling) or GR (floating) ruptures. If GR, you specify the a and b values
c where a is the interval rate of M0+-0.05 events and b is slope parameter:
c             rate(M)=10**(a-bM),
c as well as the magnitude interval over which this is valid [Mmin, Mmax].
c The code will consider M from Mmin to Mmax (no shifting of M to interior).
c  If "cascadia" is in the input file name, M8.8 and up are slab filling even though they may
c enter as itype 2 (GR or floating ruptures). This is a 2007 special feature.
c 
c  Mar 4 2009 Crouse91 atten for soil. 
c 11/19/2008: add 4 and 5 s to getGeom. These LP are not available for ABsub.
c 11/18/2008: add 0.75 s or 1.33 hz to zhao.
c 10/21/2008: add 1.33 hz to ABsub (also cublic spline). this freq. is in demand.
c 10/17/2008: 2 hz and 3.33 hz coeffs also modified in getABsub. I used cubic splines.
c 10/10/2008: ABsub add 2.5hz to both interface and inslab. Atkinson BSSA erratum Oct 2008 says
c     the 5hz model will be improved by weighing in the 2.5hz . Increase dim of perab.
c 9/16/2008: ABsub add 3s to the inslab set (which is not used). Increase npmx to 8.
c 3/25/2008: cosmetic changes, spaces replace tabs.
c 2/25/2008: SI term in getzhao is corrected. Previously, due to a typo,
c            this attenuation model was using a term with S1 (rather than SI)
c
c 2/13/2008: trim edge of final floating source to stay within defined slab
c 2/01/2008: add more period coeffs to Geomatrix.
c 10/29/2007: deaggregation option available. Just add site Lat Long nper sa(j),j=1,...,nper
c These fixed sa values must be in same order as the period set in input file.
c
c 8/07/2007: Deterministic output median and sigma added. just consider closest
c rupture. To  get this need distance to entire slab, and try to match it to within 1 km
c among the considered mini-ruptures. First one to match is reported
c
c 4/16/2007: 2 getSadigh subroutines, getSadighR foor rock, getSadighS for soil.
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
c
c There is no magmin or maxmax per se in the p() array any more.
c
c Depth to slab fixed at 20 km in Geomatrix. Probably should program variable H(?),
c but small variation has very little effect.
c Important mod: THe GR relation does not recenter the lower and upper bound
c magnitudes. They are as stated in input file. For example, if you input
c M8.0 and 8.7 as the bounds, then in the computations Mmin is 8.0 and Mmax is 8.7. 
c This conforming to input is in distinction
c to hazFXnga and other FX codes, where code recenters Mmin to 8.05 and Mmax to 8.65.
c  It is important to keep this distinction in mind for the 2007 Cascadia models
c that we consider.
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
c --- USGS. 2006 & 2007 revisions by Steve Harmsen, USGS.
c--- this version handles multiple periods and gnd motion relations
c--- to compile on SUN Solaris: cc -c iosubs.c 
c                 f95 -o hazSUBXnga hazSUBXnga.f iosubs.o -fast -e -ftrap=%none
c
c To see source listing (*.lst) use -Xlist:
c f95 -o hazSUBXnga hazSUBXnga.f iosubs.o -e -ftrap=%none -Xlist
c  The last ftrap flag seems to be necessary for many runs to avoid fatal errors. 
c
c --- To compile with gnu fortran (often available on PCs with Windows OS):
c          gfortran hazSUBXnga.f iosubs.o -o hazSUBXnga.exe -O -ffixed-line-length-none -ffpe-trap= -static
c Use the latest available gfortran e.g. March 2008 for best results.
c --- To compile with Intel fortran (ifort):
c      ifort hazSUBXnga.f iosubs.o -o hazSUBXnga -O -132  (perhaps other flags)
c
c
c--- to run:             hazSUBXnga inputfile
c -- or                  hazSUBXnga inputfile > output.log
c Use 2nd version if you want to save some diagnostic information in the log file
c Deaggregation example:: 
c
c  hazSUBXnga cascadia.bot.8387x.in 44.4 -123.4 3 .12 .14 .34
c
c The first arg after pgm name is the file name, 2nd is lat, third is long, 4th (3) nper,
c and 5 thru 7 are PGA,1s SA, and 0.2s SA, resp. 
c PGA, 1s and 0.2s SA are the order of T in cascadia.bot.8087x.2s.in in this example.
c The command-line information replaces (supercedes) the corresponding info in
c the input file. For example, input file might have 20 sa levels, but there is
c always just one sa level per period for deaggregation, which is input on command line.
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
c just including subduction-related models in hazSUBXnga (Sadigh grandfathered in).
c--- output of this program is input into hazallXL.v2 to get probabilistic motions
c iatten=1 unused, would correspond to a model that uses Rjb if there was one.
c iatten=2 Geomatrix subduction (Youngs et al, SRL, 1997). 3 states:
c            rock for vs30>700, firmsoil(avg rock&soil) for vs>444, or soil
c iatten=3 Sadigh et al. Rock 
c iatten=-3 Sadigh Soil vers. This subroutine is called with fractional wt for Soil C
c    SadighS corrected 4/16/2007. 1hz Median appears too strong compared to others
c iatten=4 Atkinson and Boore, Cascadia, BSSA, August 2003. Has NEHRP site class factors
c iatten=5 Atkinson and Boore, Global, ditto. Has NEHRP site class factors. added Nov11 06 
c            Note: Global>Cascadia median generally. Global is smidge lower for 1s though
c iatten=6 Crouse
c  Crouse has a soil model for subduction sources (Eq Spectra,1991, p201). 
c added this one march 2009. checked 1s, 4s and pga against eq spectra graphs.

c iatten = 7 Zhao et al(BSSA,  June 2006). rock, C soil, and D soil. with mag-cor term.
c            Zhao has intraplate source class and this model is active when islab = 1 in calling list.
c iatten = 8 Kanno (BSSA,  June 2006). rock,c
c iatten=15 Gregor et al. (BSSA,2002) rock and soil. Do not use. Replace with 16.
c iatten=16 Gregor et al. (SRL, 2006) w/variable siteamp. 
c    Ground motions are damped response spectra (5%) and PGA. PGV, PGD n/a.
c--- ground motion levels should be in units of g
c--- period=0 indicates PGA
      parameter (npmx=8)
      parameter ( pi=3.14159265,tiny=1.e-16)
      dimension dmin2(50,800,6)
      integer, dimension (5):: nmag0
      real, dimension(2,npmx) :: sumwt
      real, dimension (0:37) :: perka
      integer, dimension(10):: iperg,iperk,ipers,ipgeom,iperab,ipzh,ipercbc
c type replaces structure & requires gfortran or f95 to compile. It is included
c so that this code has compatibility with Linux cluster computer compilers.
c Also compatible with Windows computers and gfortran compiler.
      type header_type
        character*30 :: name(6)
        real*4 :: period
        integer*4 :: nlev
        real*4 :: xlev(20)
        real*4 :: extra(10)
      end type header_type
      type(header_type) :: headr,hdv30
      logical verbose, on_top_of_rup
      real, dimension(500):: magmin,magmax,ar
c  ar=aspect ratio for floating ruptures.
      real, dimension(33):: slat,slong
      character*4 sname(33),date*8,time*10,zone*5
      real*8 dp,pr0,emax/3.0/, sqrt2/1.4142135623/
      real vs30,prl,prlr,prr,rtot,eps
      integer readn,nrec,i,ival(8),ivs,immax/1/,npd,nper
c ivs indicator variable for NEHRP class, using Zhao subroutine, added May 21 2007
      character nameout*78,name*80,adum*80
      common/mech/ss,rev,normal,obl
      dimension p(25010),v30(128000),perg(25)
c new v30 array: can have different vs30 at each grid point. Or fixed vs30       
      logical ss,rev,normal,obl,determ
      logical norpt(8,8),isbig,isclose,deagg/.false./
c determ = .true. if deterministic output desired. To get this, make nrec < 0. -2 = 2 sta,
c         nrec = -100 = grid of stations with deterministic output. added aug 7.
c A report will be generated if source is big (isbig .true.) and close (isclose .true.)      
c These are relative words, and should have the right meaning in below context. The
c purpose of these logical controls is to try to keep the length of the report down to a
c managable size by eliminating not-so-big and/or not-so-close candidates.
      logical grid,isok      !grid=true if stations form a regular grid.
c iconv() is not used in this code. However, it is read in.
      dimension prob(250,20,npmx),xlev(20,npmx),out(800000), nlev(npmx),
     + iatten(npmx,5),iconv(npmx,5)
      dimension period(npmx),perx(8),pergeo(13),perz(22),safix(npmx)
      integer, dimension(npmx) :: nattn,ifp
      real,dimension(500):: tlen,a,b,dmag
      dimension itype(500),iftype(500),cmag(500),crate(500)
      real, dimension (10):: perab,percbc
      real,dimension(6):: ptail
      dimension rate(50),ratem(50)
      dimension nrup(50)
      dimension tarray(2)
      dimension wtdist(npmx,6),wt(npmx,6,2)
      real, dimension (16,16,10,npmx):: rbar,mbar,ebar,haz
      real prob5(16,16,10,npmx,5)
      real piriod      !gregor approximation march 14 2007
       logical v30a,cascadia,abb03(8,7), rock,soil,yakatag
c Diff relations have different spectral periods available. 
       perab = (/0.,0.2,1.0,0.1,0.3,0.4,0.5,0.75,2.0,3.0/)	
cPERab: Atkinson-Boore 2003. add 0.4 oct 10 2008. add 0.75 s Oct 21
c percbc is spectral period set for Crouse91.
	percbc = (/0.0,0.1,0.2,0.3,0.5,0.75,1.0,2.0,3.,4.0/)
       perx = (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
c gregor period set, perg
c the 2006 version of perg below. 0.0 replaces 0.01-s for consistency with others
      perg = (/ 0.0, 0.02, 0.025, 0.03, 0.04, 0.05, 0.056, 0.063, 0.071428571,
     1               0.083, 0.1, 0.125, 0.143, 0.167, 0.2, 0.25, 0.333, 0.4, 
     2               0.5, 0.769, 1.0, 1.667, 2.0, 2.5, 5.0 /)
c perka = period set for Kanno et al., BSSA 2006. 0 = pga.
       perka=(/0.0,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.15,0.17,0.20,
     +0.22,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,
     +1.30,1.50,1.70,2.00,2.20,2.50,3.00,3.50,4.00,4.50,5.00/)
c Geomatrix period set, pergeo. Addd 4 and 5 s. Nov 2008
         pergeo= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,.4,0.75,1.5,3.,4.,5./)    
c zhao period set, perz
      perz = (/0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.40,
     &     0.50,0.60,0.70,0.75,0.80,0.90,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0/)
       ptail = (/0.9772499,.8413447,0.5,.1586553,.022750139,
     + 0.0013499/)
c On solaris systems, at least, software can keep tract of user
      call getlog(adum)
      headr%name(5)='Run by '//adum(1:24)
      call getarg(0,adum)      !keep tract of program name
      headr%name(3)= adum      !new 11/05. prog name (needs to be shortpath)
      if(iargc().lt.1)stop'Usage hazSUBXnga input.file > log.file'
      if(iargc().eq.2)then
      open(44,file='subduct.rates.tmp',status='unknown')
      write(6,*)'An output file subduct.rates.tmp will be written'
      verbose=.true.
      write(6,*)'this file contains rates of various mag slab events'
      write(6,*)'under each site. Testing feb 4 2008'
      endif
      if(iargc().gt.2)then
c deaggregation information
      deagg=.true.
      safix=0.1
      call getarg(2,adum)
      read(adum,'(f7.4)')rlatd
      call getarg(3,adum)
      read(adum,'(f9.4)')rlond
c rlatd rlond are the coordinates of the deagg site. These override the stuff
c in the input file.
      call getarg(4,adum)
      read(adum,'(i1)')npd
      do i=5,npd +4
      call getarg(i,adum)
      read(adum,'(f6.4)')safix(i-4)
      write(6,*)'#command line sa ',safix(i-4)
      enddo
      endif
        call date_and_time(date,time,zone,ival)
      call getarg(1,name)
c if running on PC, may need the following call, commented out on Solaris
c      call initialize()
      abb03=.false.      !initialize 
      rtot=0.0
      sumwt= 0.0      !initialize. used to check atten model weights
      haz=0.
      rbar=0.
      ebar=0.
      mbar=0.
      prob5=0.
      headr%name(1)= name
      cascadia=index(name,'cascadia').gt.0
	yakatag=index(name,'AK.yak').gt.0
c if cascadia, all M>=8.8 fill the subduction surface with slip. AF mod added here.
c if yakatag, all M>=8.1 fill the slab
      dum=dtime(tarray)
c      write(*,*)tarray
      inquire(file=name,exist=isok)
      if(isok)then
      open(unit=1,file=name,status='old')
      write(6,*)'HazSUBXnga updated Nov 18 2008;  input file ',name
      write(6,*)'Date of run ',date,' at ',time
      else
      write(6,*)'hazSUBXnga: File not found: ',name
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
cc could use single precision for erf call. However, could be a d.p. routine
c        pr00=pr0
        p(i)=((erf(pr0)+1.)*0.5-prr)*prlr
        enddo
        determ = .false.
        p(25002)=1.0
c End initializing the truncated normal PEx() array=p()
c The indep. variable is a real*8 to reduce discretization error.
c      write(6,*) "enter name of output file"
      write(6,580)'Enter a zero for grid of sites 1 to 30 for list: '
      read(1,*)nrec
      write(6,*)nrec
      if(nrec.eq.0.or.nrec.lt.-40)then
      if(nrec.lt.-40)then
      determ=.true.
      nrec=0
      else
      determ=.false.
      endif
      grid=.true.
      write(6,580) "For sites: enter min lat, max lat, dlat: "
      read(1,*) ymin, ymax, dy
      write(6,*)ymin,ymax,dy
      if(ymin.gt.ymax)then
      tmp=ymax
      ymax=ymin
      ymin=tmp
      write(6,*)' input file error was corrected ymin, ymax'
      endif
      write(6,580) "for sites: enter min lon, max lon, dlon: "
      read(1,*) xmin, xmax, dx
      if(xmin.gt.xmax)then
      tmp=xmax
      xmax=xmin
      xmin=tmp
      write(6,*)' input file error was corrected xmin, xmax'
      endif
      write(6,*)xmin,xmax,dx
      nx= nint((xmax-xmin)/dx) +1
      ny= nint((ymax-ymin)/dy) +1      !nint will make this calc more robust,
c portable among different computers
      write(6,*) nx,ny,' nx ny for discrete grid of sites'
      nrec= nx*ny
      elseif(nrec.lt.33.and.nrec.gt.-33)then
      grid=.false.
      if(nrec.lt.0)then
      nrec=-nrec
      determ=.true.
      deagg=.false.
      endif
      do i=1,nrec
      write(6,580)'Station lat,long (dec deg) and name(1 word): '
580      format(a,$)      
      read(1,*)slat(i),slong(i),sname(i)
      write(6,*)slat(i),slong(i),sname(i)
      if(slat(i).lt.-88.)stop 'invalid station location. Version 4?'
      enddo
      else
      write(6,*)'hazSUBXnga expects first line of input to be nrec'
      write(6,*)'Valid nrec are 0, 1, ..., 29,30. Just read in ',nrec
      stop 'Please correct input file.'
      endif            !if grid
      if(deagg)then
c replace whatever was in the file with the command line location
      nrec=1
      grid=.false.
      slat(1)=rlatd
      slong(1)=rlond
      sname(1)='DEAG1'
      nx=1
      ny=1
      endif
c **** Enter soil Vs30 condition  ******NEW*******
c      write(6,*)"For sites, enter Vs30(m/s) "
      v30a = .false.
      read(1,*)vs30
c currently, code is reading in Vs30 but not depth to hardrock. Contrast with
c hazFXnga7, etc
c
      d2p5=2.      !not used, but depth of rock with Vs of 2.5 km/s
      write(6,*)' Vs30 m/s is ',vs30
      if(vs30.lt.90..and.vs30.gt.0.)then
      write(6,*)'Vs30 = ',vs30,'. This looks unreasonable.'
      stop'please check input file just before distance incr,dmax'
      elseif(vs30.eq.0)then
      write(6,580)'Enter the name of the binary vs30 array: '
      read (1,909)name
      inquire(file=name,exist=isok)
      if(.not.isok)stop'Vs30 array not in working directory'
      call openr(name)
      write(6,580)'Enter ymin, ymax, dy (lat, d) of Vs30 array: '
      read (1,*)vymin,vymax,dvy
      write(6,580)'Enter xmin, xmax, dx (long, d) of Vs30 array: '
      read (1,*)vxmin,vxmax,dvx
      write(6,580)'Enter a default vs30 in case site not inbounds: '
      read (1,*)vs30dflt
      if(vymin.gt.vymax)then
      y0=vymin
      vymin=vymax
      vymax=y0
c obvious repair if user has transposed. common input error.
      endif
      nvx= nint((vxmax-vxmin)/dvx) +1
      nvy= nint((vymax-vymin)/dvy) +1      
c Using the function nint will make this calc more robust,
c portable among different computers
      v30a=.true.
      nv30=nvx*nvy
      nhd=308
      call gethead(hdv30,nhd,nrd)
      call getbuf2(v30,nv30,nvread)
      write(6,*)'For vs30 expect ',nv30,' got ',nvread
      write(6,*)'header xmin,xmax,dx=',hdv30%extra(2),hdv30%extra(3),hdv30%extra(4)
      endif
c      read(5,*) idum
cccccccccccccccccccccccccccccccccc
c      write(6,*) "enter number of periods"
      read(1,*) nper
      write(6,*)'Number of periods ',nper
      if(nper.gt.npmx)stop 'this number exceeds npmx. Software limit'
c---loop through periods
      do 40 ip=1,nper
         read(1,*) period(ip)
         if(period(ip).lt.0.)stop'hazSUBXnga: PGV PGD are not available in subduction code'
c         write(6,*) "enter name of output file for this period"
         read(1,909) nameout
 909     format(a)
       if(grid)then
      call openwx(ifp(ip),nameout)
      write(6,*) ifp(ip),nameout
      else
      open(9+ip,file=nameout,status='unknown')      !ascii
       write(6,895)nameout
895      format('Ascii output file: ',a)
      if(deagg)then
c cannot have both deaggregation and determinstic output. Safe to use same unit number.
      if(nper.ne.npd)stop'Deaggregation command-line np does not equal input nper'
       isz=index(nameout,' ')-1
      open(20+ip,file=nameout(1:isz)//'.DEAG',status='unknown')
      write(20+ip,907)period(ip),date,time,name,rlatd,rlond
 907      format('#hazSUBXnga deagg @fixed GM, T='f5.2,' s, run on ',a,' at ',a,' fi ',a,
     +/,'#sta lat long = ',f7.4,1x,f9.4)      
      endif      !if deagg. new oct 26
      endif      !if grid
      if(determ)then
       isz=index(nameout,' ')-1
      open(20+ip,file=nameout(1:isz)//'.DET')
      write(20+ip,677)period(ip),date,time,name
      endif
677      format('#hazSUBXnga sources. Sp_Per= ',f5.2,' s. Run on ',a,' at ',a,/,
     +'# *** Input control file: ',a,/,
     +'#MedianSA(g)     sd   iattn  irup    M  Rjb(km) Rcd(km)  Wt*Rate',
     +' Dtor(km)')
         write(6,580) "Number of atten. relations for this period: "
         read(1,*) nattn(ip)
         write(6,*)nattn(ip)
c--- loop through atten relations for that period
         do 20 ia=1,nattn(ip)
            write(6,580)
     +      "Type of atten. relation, weight, mb to M conv. "
            read(1,*) iatten(ip,ia),wt(ip,ia,1),wtdist(ip,ia),
     +      wt(ip,ia,2),iconv(ip,ia)
            write(6,*) iatten(ip,ia),wt(ip,ia,1),wtdist(ip,ia),
     +      wt(ip,ia,2),iconv(ip,ia)
              sumwt(1,ip)=sumwt(1,ip)+wt(ip,ia,1)
              sumwt(2,ip)=sumwt(1,ip)+wt(ip,ia,2)
c--------  Geomatrix interface atten.
         if(iatten(ip,ia).eq.2) then
         j=1
         dowhile(period(ip).ne.pergeo(j).and.j.lt.13)
         j=j+1
         enddo
         if(j.gt.13)stop 'period doesnt match Geomatrix subd. set.'
         ipgeom(ip)=j      !period index for Geomatrix relations.
         write(6,*)'Period map for Geomatrix ',j
         headr%name(6)='Geomatrix attn'      !period-indep assignment.
c        write(6,*)'Geomatrix is subroutine, no coeff. readin Oct 2006'
         elseif(iatten(ip,ia).eq.7) then
         j=1
         dowhile(period(ip).ne.perz(j).and.j.le.22)
         j=j+1
         enddo
         if(j.gt.22)stop 'period doesnt match available Zhao set 11-2008'
         ipzh(ip)=j      !period index for Geomatrix relations.
         write(6,*)'Period map for Zhao et al.',j
         headr%name(6)='Zhao et al attn'
c        zhao et al, new may 2007.
c------- Sadigh type
            elseif(abs(iatten(ip,ia)).eq.3) then
            if(iatten(ip,ia).eq.-3)then
            write(6,*)'Sadigh model with soil coeffs.'
            else
            write(6,*)'Sadigh model with rock coeffs.'
            endif
         j=1
         dowhile(period(ip).ne.perx(j).and.j.lt.npmx)
         j=j+1
         enddo
         if(j.gt.7)stop 'period doesnt match standard set'
         ipers(ip)=j      !period index for Sadigh and perhaps other relations.
         write(6,*)'Period map for Sadigh ',j
       elseif(iatten(ip,ia).eq. 8)then
      ka=0
      dowhile(period(ip).ne.perka(ka))
      ka=ka+1
      if(ka.gt.37)then
      write(6,*)'Period ',period(ip)
      stop 'For this period please remove the Kanno relation from input file'
      endif
      enddo
       iperk(ip)=ka     
       write(6,*)ip,iperk(ip),' Kanno 10/06 ip map'
c vsfac used for site-amp in  the Kanno et al. relation.
        vsfac=alog10(vs30)          
         elseif(iatten(ip,ia).eq.15)then
c Gregor 2002 added oct 2006. Gregor Email says this has been updated.
         write(6,*)'Gregor 2002 subduction model out of date'
         stop 'please use atten model with index 16 for Gregor.'
         elseif(iatten(ip,ia).eq.16)then         
         j=1
c gregor : treat .333 T as 0.3 to combine with other modelers' output.
c alternate might be to interpolate coeffs (not done as of mar14 2007)
         if(period(ip).gt.0.299.and.period(ip).lt.0.345)then
         piriod=0.333
         write(6,*)'Note: Using gregor 0.333 T for 0.3 s analysis'
         else
         piriod=period(ip)
         endif
         dowhile(abs(piriod-perg(j)).gt.0.01.and.j.lt.26)
         j=j+1
         enddo
         if(j.gt.25)stop 'period doesnt match standard Gregor set'
         iperg(ip)=j      !period index for Gregor subduction
         write(6,*)'Period map for Gregor ',j
c updated Gregor logic Jan 2007. Steve Harmsen         
c------- AB03 type
            elseif(iatten(ip,ia).eq.4.or.iatten(ip,ia).eq.5) then
         write(6,*)'Atkinson Boore 2003 subduction model'
        abb03(ip,ia)=.true.
         j=1
         dowhile(period(ip).ne.perab(j).and.j.lt.10)
         j=j+1
         enddo
         if(j.gt.10)stop 'period doesnt match standard AB subduction set'
         iperab(ip)=j      !period index for ABsub and perhaps other relations.
         if(iatten(ip,ia).eq.4)then
         write(6,*)'Using ABsub with Cascadia coef, period map is ',j
      headr%name(6)='ABsub Cascadia'
c iclass=1 for PacNW, iclass=3 for worldwide...
        iclass=1
         else
         write(6,*)'Using ABsub (10/17/2008) with global coef, period map is ',j
         iclass=3
      headr%name(6)='ABsub Global'
         endif
c------ Crouse attenuation form
            elseif(iatten(ip,ia).eq.6) then
            if(vs30.gt.550.)stop'Crouse model was designed for soil site-condition only'
         j=1
         dowhile(period(ip).ne.percbc(j).and.j.le.10)
         j=j+1
         enddo
         if(j.gt.10)stop 'period doesnt match CB Crouse subduction set'
         ipercbc(ip)=j      !period index for ABsub and perhaps other relations.
            print *,j, 'Crouse index for soil available (2009)'
            else
            write(6,*)' Surprise atten index. ',iatten(ip,ia)
            write(6,*)' Code does not know how to handle this one;'
            write(6,*)'Please review your input file. Valid -3, 2 to 8 or 16'
            stop 'hazSUBXnga.f: did not finish. Check log file.'
            endif
c----- Boore look-up table: removed 10-2006. SH
   20    continue
         write(6,580) "Number of ground motion levels: "
         read(1,*) nlev(ip)
         write(6,*)nlev(ip)
         if(nlev(ip).gt.20)stop 'number exceeds limit of 20'
         write(6,*) "Ground motion levels "
         read(1,*) (xlev(k,ip),k=1,nlev(ip))
         write(6,*) (xlev(k,ip),k=1,nlev(ip))
         if(deagg)then
         nlev(ip)=1
         xlev(1,ip)=safix(ip)
         endif
         do 30 k=1,nlev(ip)
   30    xlev(k,ip)= alog(xlev(k,ip))
   40 continue
ccccccccccccccccccccccccccccccccccc
c warn the user if atten model weights don't sum to 1.0.
      do ip=1,nper
      if(sumwt(1,ip).lt.0.99.or.sumwt(1,ip).gt.1.01)then
      write(6,*)'Unusual sum of att. model wts for inner annulus ',sumwt(1,ip)
      write(6,*)'period index ',ip
      endif
      if(sumwt(2,ip).lt.0.99.or.sumwt(2,ip).gt.1.01)then
      write(6,*)'Unusual sum of att. model wts for outer annulus ',sumwt(2,ip)
      write(6,*)'period index ',ip
      endif
      enddo
      write(6,580) "Increment dlen (km) for source segment "
      read(1,*) dlen
      write(6,*)dlen
      write(6,580) "Distance increment, dmax (km): "
      read(1,*) di,dmax
      write(6,*) di,dmax
ccc-------read fault data
      do 50 ift=1,1000
         write(6,*) ift
         write(6,580)
     +   "Enter 1 for char. 2 for G-R; 1 for SS, 2 for reverse: "
         read(1,*,end=60) itype(ift),iftype(ift)
          write(6,*) itype(ift),iftype(ift)
        rev=iftype(ift).eq.2
         ss=iftype(ift).eq.1
         normal=iftype(ift).eq.3
         if(normal)stop' Attenuation models not prepared for normal slip'
         if(itype(ift).eq.1) then
            write(6,*) "Char. M and rate: "
       if(verbose)write(44,4442)cmag(ift),crate(ift),name
4442      FORMAT('Rx    Ry  Mchar  ',f6.2,' Rate_Mchar ',e10.5,' filein ',a20)
            read(1,*) cmag(ift),crate(ift)
            write(6,*)cmag(ift),crate(ift)
         endif
         if(itype(ift).eq.2) then
            write(6,580) "enter a,b,min M,max M,dmag for freq mag curve"
            read (1,*) a(ift),b(ift),magmin(ift),magmax(ift),dmag(ift)
            write(6,*) a(ift),b(ift),magmin(ift),magmax(ift),dmag(ift)
c DO NOT Reset the GR Magmin & Magmax as in fltrate, mod of Oct 27 2006. SH.
c Unmodded feb 28 after discussion with Art. Be sure incoming rate is appropriate.
c
c            if(magmax(ift).ge.magmin(ift)+dmag(ift))then
c            magmax(ift)=magmax(ift)-0.5*dmag(ift)
cc            magmin(ift)=magmin(ift)+0.5*dmag(ift)
c            endif
            nmag0(ift)= nint((magmax(ift)-magmin(ift))/dmag(ift)) + 1
      if(verbose)write(44,4444)magmin(ift),magmax(ift),nmag0(ift),name
4444      FORMAT('Rx    Ry  Total_Rate  Rate_M(j) from ',f6.2,' to ',f6.2,' n ',i2,
     +' filein ',a20)
            if(magmin(ift).lt.6.5)stop'hazSUBXnga: magmin below 6.5'
         endif
c******************************************************************
c      Replace code to read and setup model of subduction zone
c******************************************************************
c     Insert new code to read and setup model
c***************************************************************
         call readmodel(ift,tlen) 
c****************************************************************
c     End new code to read and setup model
c****************************************************************
         write(6,*) tlen(ift), ' total length (km)'
         coef= pi/180.
c      dip1(ift) =dip1(ift)*coef
c      dip2(ift) =dip2(ift)*coef
   50 continue
   60 nft= ift-1
      write(6,*) "nft=",nft
ccccccccccccccccccccccccccccc
c---- write header --------------
c Discussion at Software meeting Dec 16 2005: header record needs to be
c lengthened to include more input-file information.
      headr%name(4)=date//' at '//time
      do 80 ip=1,nper
      headr%period= period(ip)
      headr%nlev= nlev(ip)
      do 702 k=1,nlev(ip)
 702  headr%xlev(k)= exp(xlev(k,ip))
      ndata= 308
c New 10/05: store solution grid info.
      headr%extra(2)=xmin
      headr%extra(3)=xmax
      headr%extra(4)=dx
      headr%extra(5)=ymin
      headr%extra(6)=ymax
      headr%extra(7)=dy
      headr%extra(8)=float(nrec)
c store other run-specific info. extra(9) and (10)      
      headr%extra(9)=vs30
      headr%extra(10)=d2p5      !depth of fast rock (not used here)
         ndata= 308
      if(grid)then
         call puthead(ifp(ip),headr,ndata,readn)
         write(6,*) 'Finished writing header period index ',ip
      endif
   80 continue
ccccccccccccccccccccccccccccccccccc
      ndist= dmax/di +1.5

c no longer building a table of p(). Just get the mean and sd when you need 'em
c mod oct 19 2006. SH
ccccccccccccccccccccccccccccccccccccccccccc

      dum=dtime(tarray)
c      write(*,*) tarray


ccccccccccccccccccccccccccccccccccc
      icnt=1
       prob=1.e-21      !initialize array
ccccccccccccccccccccccccccccccccccc
c---Assume only one fault
      ift =1
c---Set up parameters of rupture zones
      do 185 m=1,nmag0(ift)
          xmag = magmin(ift)+(m-1)*dmag(ift)
c------find rupture length from Geomatrix interface eqs relation
          ruplen= (xmag-4.94)/1.39
          ruplen= 10**ruplen
          nrup(m)= max(1,nint((tlen(ift)-ruplen)/dlen) +1)
c max of 1 and the other expression because if tlen < ruplen, will get 0 or less than 0         
       write(6,*) "nrup(m),ruplen,tlen",nrup(m),ruplen,tlen(ift)
      if(nrup(m).gt.799)stop' nrup(m) exceeds 799. Up dim 2 of rupture()'
          if(xmag.ge.8.8.and.cascadia) then
c above mod feb 28 to be consistent with AF models for 2007 PSHA update.
                nrup(m)=1
                ruplen= tlen(ift)
                write(6,*) xmag, ruplen, tlen(ift)
c                read(5,*) idum
	elseif(xmag.ge.8.1.and.yakatag)then
c above mod mar 10 2009 for Yakataga largest M slab filling
		nrup(m)=1
		ruplen= tlen(ift)
		write(6,*)'Yakatag filling zone for M8.1+',ruplen
		endif
          if(ruplen.ge.tlen(ift)) then
            nrup(m)=1
                ruplen=tlen(ift)
          endif
          rate(m)= a(ift) - b(ift)*xmag
          rate(m)= 10.**rate(m)
      rtot=rtot+rate(m)
      write(6,*) xmag,rate(m),' M ann rate resp'
c       read(5,*) idum
          rate(m)= rate(m)/float(nrup(m))
          write(6,*)'Rate of a given scenario at this M ',rate(m)
          do 182 irup=1,nrup(m)
               xmove= dlen*(irup-1)
               s1 = xmove
               s2 = xmove + ruplen
c known problem if s2 extends even 1 km beyond end of fault with length tlen. next line protects.
               s2=min(s2,tlen(ift)-0.1)      !new feb 13 2008
               call rupturesetup(m,irup,s1,s2)
 182      continue
       write(6,*)'For this mag the final distance s2 is ',s2,' km'
 185   continue
c
c       if(determ)then
c call rupture setup for a fictitious event that covers the whole fault surface
c store this in m= 24 bin.
      irup=1
      call rupturesetup(24,irup,0.,tlen(ift))
      print *,' megathrust setup ok '
c      endif
      write(6,*)' Cumulative annual rate of subd events: ',rtot
      write(6,*)' This is equivalent to recurr. interval of ',1./rtot,' yrs'
c---Here's the guts
c
c---loop through receiver sites 250 at a time. Why 250 at a time? big mystery.
      i = 0
      nloop = (nrec-1) / 250+ 1
      do 370 iloop=1,nloop
 
c     loop through floating rupture zones
c     find minimum distance dmin2 to receiver for each rupture zone
c     nrup(m) rupture zones per magnitude bin m
c     move rupture zones by dlen horizontally along fault
 
c*****calculate distances for 250 receivers from all sources*********
        
         mloop = 250

       if(.not.grid)mloop =nrec      
c !max 30 stations if they occur in a list
         do 310 jloop=1,mloop
                  icnt = jloop
c                       i is the index of the receiver
                  i = (iloop-1)*250 + jloop
                  if(i.gt.nrec) then
                     i = i - 1
                     icnt = jloop -1          
                     go to 320
                  endif
c the below logic worked for hazFXnga5. With all the loops above, doubts creep in. 
      ratem=0.
        if(grid)then
      iy= (i-1)/nx 
      ix= i-1-iy*nx
      rx= xmin + float(ix)*dx
      ry= ymax - float(iy)*dy
        else
        rx=slong(i)
        ry=slat(i)
        write(6,*)'Computing seismic hazard for ',sname(i)
        endif
        if(v30a)then
c find nearest gridpoint in vs30 array to the site with coords rx,ry
        if(rx.gt.vxmax.or.rx.lt.vxmin.or.ry.gt.vymax.or.ry.lt.vymin)then
        vs30=vs30dflt
        else
        ivx=1+nint((rx-vxmin)/dvx)
        ivy=nint((vymax-ry)/dvy)        !same organization as we are used to.
        vs30=v30(nvx*ivy+ivx)
        endif   !inbounds
        endif   !array rather than scalar vs30
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
       if (determ)then
      do ip=1,nper
      write(20+ip,6684)rx,ry,i,vs30
6684      format('! ',f8.3,1x,f6.3,' site ',i7,' vs30 ',f6.2,' m/s')
      enddo
                  norpt = .true.      !start with no report
      endif      !this if clause added aug 7 for comparing ruptures with mega.
c make a determination if site is more than Rmax from fault. If so do not
c bother with any more calcs in the do 300 below.
              call qkdist(24,1,rx,ry,d_min,fltdepmin, d_minxy)
c  d_min is minimum distance from fault to receiver
c d_minxy is the minimum distance for megathrust that ruptures the entire slab 
       if(d_min.gt.dmax)goto 305   
c Only report the rupture that is closest to site, whose distance should about equal
c distance to megathrust. That is the hope anyway.
          do 300 m=1,nmag0(ift)
          xmag = magmin(ift)+(m-1)*dmag(ift)
          isbig=m.eq.nmag0(ift)
          if(i.eq.1)write(6,*)' Mag ',xmag,' nrup is ',nrup(m)
c dM bins have width 0.1 here 
                  if(deagg)then
                  im=nint((xmag-magmin(ift))*10.)+1
                  im=min(16,im)      !current dimension 16
                  immax=max(immax,im)
                  endif
              do 290 irup=1,nrup(m) 
c-----------------------------------------------------------------------
c  dmin is minimum distance from fault to receiver
                  call qkdist(m,irup,rx,ry,distmin,fltdepmin, distminxy)
                  on_top_of_rup=distminxy.lt.1..and.distmin.lt.20..and.verbose
c                   on_top_of_rup=distmin.lt.20..and.verbose
                  dmin2(m,irup,1) = distminxy
c                  if(ry.eq.49.and.m.eq.1)
c     + write(*,*) "m,irup,rx,ry,distmin",m,irup,rx,ry,distmin,distminxy
c         dmin2(m,irup,2) = fltdepmin
c                  dmin2(m,irup,2) = sqrt(distminxy *distminxy +
c     +            400.)
                  dmin2(m,irup,2) = distmin
                  dmin2(m,irup,3) = distmin
                  dmin2(m,irup,4) = distmin
                  isclose= abs(distmin-d_min).lt.1.5      
      if(deagg)then
        if(distmin.le.100.)then
        ir=0.1*distmin + 1
        elseif(distmin.le.300.)then
        ir=11+(distmin-100.)/50.
        elseif(distmin.le.500.)then
        ir=15
        else
        ir=16
        endif 
       endif            !if deagg
c*****calculate hazard for 250 receivers from all sources*********
c By delaying calculation of p() at receiver R, we can use variable vs30 if requested.
c  Mod Oct 2006 SH 
              do 280 ip=1,nper
              ipdet=ip+20
                           do  ia=1,nattn(ip)
c always use R_cd or Rrup as the distance measure. No subduct. relations use R_jb
                           dist= dmin2(m,irup,3)
                           rcd=distmin
                              ii= dist/di +1.5
                              if(ii.gt.ndist) go to 270
         weight= wt(ip,ia,1)
         if(dist.gt.wtdist(ip,ia)) weight=wt(ip,ia,2)
         wt2= wt(ip,ia,2)
         wt1= wt(ip,ia,1)
         dist0= wtdist(ip,ia)
         if((dist.ge.dist0-taper).and.(dist.le.dist0+taper)) 
     &   weight= wt1 + (wt2-wt1)*(dist-dist0+taper)*taperfac
        
       if(weight.lt.0.0001)goto 270      !some relations downweighed with distance
      if(abb03(ip,ia))then
      call getABsub(iperab(ip),iclass,0,vs30,xmag,dist,gnd,sigmaf)
      elseif(iatten(ip,ia).eq.2) then
c the 2nd arg, 0, refers to subduction or interface. 0 for subduction. Added oct 26
c
            call getGeom(ipgeom(ip),0,vs30,xmag,rcd,gnd,sigmaf)
      elseif(iatten(ip,ia).eq.7) then
c the 2nd arg, 0, refers to subduction or interface. 0 for subduction. Added oct 26
c
            hslab = 50.      !not used for subduction. hi is used (param)
             call zhao(ipzh(ip),xmag,rcd,gnd,sigmaf,ivs, 0,hslab)
c Zhao: the second-to-last argument, 0, indicates compute for subduction, not for in-slab source.
      elseif(iatten(ip,ia).eq.8) then
      call kanno(iperk(ip),gnd,sigmaf,xmag,rcd,vsfac)
c--- following are Sadigh forms. Might allow a near-source sigma someday.
         elseif(iatten(ip,ia).eq.3)then
               call getSadighR(ipers(ip),xmag,rcd,gnd,sigmaf,0.,0.)
         elseif(iatten(ip,ia).eq.-3)then
               call getSadighS(ipers(ip),xmag,rcd,gnd,sigmaf,0.,0.)
      elseif(iatten(ip,ia).eq.6) then
c C.B. Crouse relation Added Mar 4 2009. The 20 is depth of interface.
c
            call Crouse91(ipercbc(ip),xmag,rcd,vs30,20.,gnd,sigmaf)
c       elseif(iatten(ip,ia).eq.15)then
c             call Gregor02(iperg(ip),xmag,rcd,vs30,gnd,sigmaf)
c gregor 2002 only works for subduction, is therefore outputting "gnd2"
       elseif(iatten(ip,ia).eq.16)then
             call Gregor06(iperg(ip),xmag,rcd,vs30,gnd,sigmaf)
c gregor 2006 only works for subduction, is therefore outputting "gnd2"
         endif
      wtrate = weight*rate(m)
c keep books on events that are directly under site (For Ned Field) (verbose)
      if(on_top_of_rup)then
      ratem(m)=ratem(m)+rate(m)
      ra_tot=ra_tot+rate(m)
      endif
         if(isbig.and.isclose.and.determ.and.norpt(ip,ia))then
         sigma=1./sigmaf/sqrt2
         write(ipdet,679)exp(gnd),sigma,iatten(ip,ia),irup,xmag,distminxy,rcd,wtrate
         norpt(ip,ia)=.false.
679     format(e11.5,1x,e11.5,1x,i2,1x,i3,1x,f6.2,1x,f7.1,1x,f7.1,1x,e11.5,
     + f5.1,1x,i2)
         endif
         do k=1,nlev(ip)
       pr=(gnd - xlev(k,ip))*sigmaf
       if(pr.gt.3.3)then
       ipr=25002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 270      !transfer out if ground motion above mu+3sigma
       endif
       cfac=wtrate*p(ipr)
               prob(icnt,k,ip)= prob(icnt,k,ip)+cfac
               if(deagg)then
c k is supposed to equal 1 if deagg. So we dont have a k index below, rbar etc
c you could deaggregate at each input spectral sa level. then you would need a k
c index on rbar, mbar, ebar, haz, prob5.
               eps= -pr*sqrt2
            if(eps.lt.emax)then
            ieps=max(1.,min((eps+2.)*2.,10.))
            rbar(ir,im,ieps,ip)=rbar(ir,im,ieps,ip)+cfac*rcd
            mbar(ir,im,ieps,ip)=mbar(ir,im,ieps,ip)+cfac*xmag
            ebar(ir,im,ieps,ip)=ebar(ir,im,ieps,ip)+cfac*eps
            haz(ir,im,ieps,ip)=haz(ir,im,ieps,ip)+cfac
            if(eps.ge.2.)goto 266
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
            ka=5
            temp=p(ipr)
            prx= temp-ptail(ka)
            dowhile(prx.gt.1.e-9)
            prob5(ir,im,ieps,ip,ka)=prob5(ir,im,ieps,ip,ka)+
     +       wtrate*prx      ! *prlr       dont know about this
            ka=ka-1
            if(ka.eq.0)goto 266
            prx= temp-ptail(ka)
            enddo
            endif      !eps <emax
  266            continue
              endif      !if deagg

            enddo      !k index
  270                continue      
        enddo      !ia index
  280                continue
  290                continue
  300             continue
            if(verbose)write(44,344)rx,ry,ra_tot,(ratem(m),m=1,nmag0(ift))
            ra_tot=0.
            ratem=0.
344      format(f7.2,1x,f7.2,6(1x,e10.5))
  305            continue      !we go here if site is > dmax from anyplace on
c the subducting interface
  310          continue
c
  320          continue
 
c*****write output for 250 receivers************************
      if(grid)then
                  write(6,*) i,prob(icnt,1,1)
                  iend= 250
                  if(i.eq.nrec) iend=icnt
                  do 340 ip=1,nper
                     ndata= iend*nlev(ip)
                     do 330 ii=1,ndata
                        i2= ii-1
                        irec= i2/nlev(ip)
                        ilev= i2-irec*nlev(ip)
                        irec= irec+1
                        ilev= ilev+1
                        out(ii)= prob(irec,ilev,ip)
c                       write(*,*) ip,irec,ilev,out(ii)
  330                continue
c                     write(*,*) "hi"
                     call putbufx(ifp(ip),out,ndata,readn)
c                      write(*,*) "bye"
                     readn= readn/4
c                     write(6,*) ndata,readn
 340            continue
            prob = 1.e-21
            else
c For list of stations, formatted output in place of binary
67      format(/,'#Station ',f7.4,1x,f10.4,1x,a,' vs30 ',f6.1)
      do 615 i=1,nrec       
      do 615 ip=1,nper
       write(9+ip,67)slat(i),slong(i),sname(i),headr%extra(9)
      write(9+ip,69) period(ip), nlev(ip)

69      format('#SP(s)',f9.2,' Nlev ',i2)
      do 614 k=1,nlev(ip)

      write(9+ip,612) exp(xlev(k,ip)), prob(i,k,ip)
 612      format(f10.6,1x,e11.5) 
 614  continue
 615  continue
c       write(10,68)sname(i)
c 68      format('#End receiver ',a,/)
      if(deagg)then
       do ip=1,nper
       do ieps=1,10
       do im=1,immax
       do ir=1,16
      prx=haz(ir,im,ieps,ip)
       if(prx.gt.tiny)then
       e0=max(-9.99,ebar(ir,im,ieps,ip)/prx)
c Future: the 0 below could be  ift, which rupture is contributing the
c mode at that ir, im, ieps location?        
       write(ip+20,25)rbar(ir,im,ieps,ip)/prx,
     + mbar(ir,im,ieps,ip)/prx,prx,(prob5(ir,im,ieps,ip,k),k=5,1,-1),
     + 0,e0,exp(xlev(1,ip))
      endif
25      format(f9.3,1x,f7.3,6(1x,e11.5),
     + 1x,i2,1x,f5.2,1x,f8.5)
       enddo
       enddo
       enddo
       close(ip+20)
       enddo
      endif      !if deagg option is requested.
        endif      !if grid or not
  370       continue
       if(.not.grid)then
       do ip=1,nper
       close(9+ip)
       enddo
       write(6,*)'An ascii file was written for each period in input'
       endif
            dum=dtime(tarray)
            write(*,*) tarray(1),' s = time to complete hazSUBXnga'
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
 
 
      dimension tlen(500)
 
 
      parameter (coef=3.14159/180.,nlevels=2)
 
      tol = 1.e-3
 
 
      do 20 l=1,nlevels
         read(1,*) npts(iflt,l)
         write(*,*) npts(iflt,l)
         do 10 i=1,npts(iflt,l)
            read(1,*) yf(iflt,l,i),xf(iflt,l,i),zf(iflt,l,i)
c             write(*,*) yf(iflt,l,i),xf(iflt,l,i),zf(iflt,l,i)
   10    continue
 
c         set up and resample
 
         write(*,*) 'resample'
         call resample(iflt,l,npts,xf,yf,zf)
   20 continue
 
      write(*,*) 'normvec'
 
      call normvec2(nlevels,tol)
c      tlen(iflt) = send(iflt,1)
       tlen(iflt)= send(iflt,1)
 
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


      subroutine qkdist(imag,jrup,xpt,ypt,distmin,fltdepmin,distminxy)
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
 
      logical above, offend
 
      external fltdst
 
c      write(*,*) 'qkdist'
      tol = 1.e-3
      x = xpt
      y = ypt
      z = 0.
 
      s1 = rupture(imag,jrup)%s1 
      s2 = rupture(imag,jrup)%s2 
c      ist1 = rupture(imag,jrup)%ist1  
c      ist2 = rupture(imag,jrup)%ist2  
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
c         ifltdepmin = zfhld(1)
 
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
               else
                  if(dotprod5.le.0.) then
                     ireg = 5
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
      else
         if(.not.offend)
     +   distminxy = sqrt(((x-xfhld(ifltmin))*cosfac)** 2+ (y-yfhld
     +   (ifltmin))**2) * 111.195
      endif
c      write(*,*) 'ireg,distmin,distminxy =', ireg,distmin,distminxy
c      write(9,*) y,x,ireg,distmin
      return
 
      end
 
 
 
 
      subroutine resample(iflt,l,npts,xf,yf,zf)
 
      dimension xf(1,5,1000),yf(1,5,1000),zf(1,5,1000)
c     yf,xf,zf(iflt,l,i) are the coordinates (latitude, longitude
c         and depth) of the i-th point on the l-th level
c         (top, hinge, bottom) of the iflt-th fault
 
      dimension npts(10,5)
 
      dimension s(5000),xs(5000),ys(5000),zs(5000)
 
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
      write(6,*) 'hazSUBXnga: Brent exceed maximum iterations.'
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


      subroutine getGeom(iq,islab,vs30,xmag,dist0,gnd,sigmaf)
c 
c  Geomatrix (Youngs et al. subduction). modified to NGA style.
c Inputs:
c	iq = index in pergeo
c 	islab=0 for subduction, islab=1 for intraplate.
c	vs30: site's avg Vs30 (m/s), with limited variation.
c vs30 >= 700  means rock site 
c vs30 > 444 Csoil (avg of rock and soil)
c otherwise, soil site.  rock & soil-coef tri-chotomy is weak. 
c should have continuous variation. 
c	xmag = moment magnitude of source
c	dist0 = slant distance (km)
c Outputs: gnd and sigmaf.
c
c modified April 2007. Steve Harmsen. FOR rock or firm or deep soil.
c added coeffs for 0.4 0.75 1.5 and 3s feb 1 2008 SH.
c added coeffs foor 4 and 5 s SA from Youngs' spreadsheet Nov 19 2008.
c vs30 >= 700  means rock site 
c vs30 > 444 Csoil (avg of rock and soil)
c otherwise, soil site.  rock & soil-coef tri-chotomy is weak. 
c should have continuous variation. 
c returns median and 1/sigma/sqrt2
      parameter (np=13,sqrt2=1.4142136)
      real gc0/0.2418/,gcs0/-0.6687/,ci/0.3846/,cis/0.3643/
      real gch/0.00607/,gchs/0.00648/,gmr/1.414/,gms/1.438/
c gch is the factor applied to interface depth, 20 km for Cascadia, for rock sites.
c gchs is the corresponding factor for soil sites.
      real period,gnd0,gndz,gz,g3,g4,gndm
       real, dimension(np):: gc1,gc2,gc3,gc4,gc5,pergeo
       real gc1s(np),gc2s(np),gc3s(np)
c array constructors oct 2006. gc1 is Youngs C8 in Nov08 spreadsheet.
         pergeo= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,.4,0.75,1.5,3.,4.,5./)    
       gc1= (/0.,0.722,-1.736,1.1880,0.246,-0.4,-3.3280,
     + -0.115,-1.149,-2.634,-4.511,-5.350,-6.025/)
       gc2= (/0.,-0.0027,-0.0064,-0.0011,-0.0036,-0.0048,-0.0080,
     + -0.0042,-.0057,-.0073,-.0089,-.0096,-.0101/)
       gc1s= (/0.,1.549,-2.87,2.516,0.793,-.438,-6.4330,0.144,-1.704,
     + -5.101,-6.672,-7.618,-8.352/)
       gc2s= (/0.,-0.0019,-0.0066,-0.0019,-0.002,-0.0035,-0.0164,
     + -.002,-.0048,-0.0114,-0.0221,-0.0235,-0.0246/)
       gc3= (/-2.556,-2.528,-2.234,-2.6550,-2.454,-2.36,-2.107,
     + -2.401,-2.286,-2.16,-2.033,-1.98,-1.939/)
       gc3s= (/-2.329,-2.464,-1.785,-2.697,-2.327,-2.140,-1.29,
     + -2.23,-1.95,-1.47,-1.347,-1.272,-1.214/)
       gc4= (/1.45,1.45,1.45,1.45,1.45,1.45,1.55,
     + 1.45,1.45,1.50,1.65,1.65,1.65/)
       gc5= (/-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,-0.1,
     + -0.1,-0.1,-0.1,-0.1,-0.1,-0.1/)
      if(vs30.gt.699.)then
c Use rock coeffs. No interface term ci or cis for subduction
        gnd0=gc0 +ci*islab
        gz=gch
        g1=gc1(iq)
        g2=gc2(iq)
          g3=gc3(iq)
          g4=1.7818
          ge=0.554
          gm=gmr
          elseif(vs30.gt.444.)then
c C soil zone: average rock and soil
        gnd0=0.5*(gc0+gcs0 +(ci+cis)*islab)
        gz=0.5*(gch+gchs)
        g1=0.5*(gc1(iq)+gc1s(iq))
        g2=0.5*(gc2(iq)+gc2s(iq))
          g3=0.5*(gc3(iq)+gc3s(iq))
          g4=0.5*(1.7818+1.097)
          ge=0.5*(0.554+0.617)
          gm=0.5*(gmr+gms)
          
      else
c Use soil coeffs
        gnd0=gcs0 + cis*islab 
        gz=gchs
        g1=gc1s(iq)
        g2=gc2s(iq)
           g3=gc3s(iq)
           g4=1.097
           ge=0.617
           gm=gms
           endif
c gz term  used for subduction, and for benioff (see hazgridXnga2)
          sig= gc4(iq)+gc5(iq)*min(8.,xmag)
          sigmaf= 1./(sig*sqrt2)
c same sigma for soil and rock.
c Thru 2008 we use a constant depth of 20 km but this depth more properly depends on subduction model.
c Some slabs may descend at a steeper or shallower angle affecting this bottom.
          gndm= gnd0 +g1 +gm*xmag +g2*((10.-xmag)**3)  +gz*20.
          arg= exp(ge*xmag)
c Distance could be hypocentral or distance to top-of-Benioff zone.
          gnd=gndm +g3*alog(dist0+g4*arg)
      return
      end subroutine getGeom

      subroutine getABsub(iq,ir,islab,vs30,xmag,dist0,gnd,sigmaf)
c from At&Boore BSSA aug 2003, p1715.  C-site and D-site look good apr 11 2007. SH.
c Modify 5hz and add 2.5 hz coeffs., 10/10/2008. AB BSSA erratum Oct 2008. 
c Affects Global but not cascadia. 10/17 also modify 2 hz and 3.33 hz coeffs for global.
c
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
c Outputs:
c gnd = ln(median sa or pga)
c sigmaf = sigma factor = 1.0/sigma/sqrt2/aln10
c            where sigma is the log_10 aleatory sigma.
c  Coeffs for 8 spectral pds. also 3s available for interface
c modified for gfortran, f95 Oct 2006
      parameter (np=10,sqrt2=1.4142136,gfac=2.9912261,aln10=2.30258509)
c gfac = log10(980)
      real rc1/-.25/,rs1/2.991/
      real rc2/0.6909/,rs2/0.03525/
      real rc3/0.01130/,rs3/0.00759/
      real rc4/-0.00202/,rs4/-0.00206/
c  Dtor allows a distribution if 
c you are uncertain about what that depth is. not using this variable dtor here.
c Global modified coeffs ab 10-2008 affected at 2, 2.5, 3.33, 5 hz.
      real, dimension(np) :: s2g,s3g,s4g,s5g,s6g,s7g
      real c1w(np)
c interface, names start with s for "Subduction". s1g global, s1 cascadia only.
c
      real,dimension(np):: s1,s1g,s2,s3,s4,s5,s6,s7,ssig
      real perab(10),period
c array constructors oct 2006. 2.5 hz is element number 6 of 10. 1.33 hz is el. 8
      perab= (/0.,0.2,1.0,0.1,0.3,0.4,0.5,0.75,2.0,3.0/)      
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
c Cannot Add 4 and 5 s. Not available. Nov 19 2008.
      s1g=(/2.991,2.5711536,2.1442,2.7789,2.6168785,2.6175463,2.536019,2.288355,2.1907,2.301
     + /)
      s1=(/2.79,2.54,2.18,2.5,2.516,2.50,2.418,2.241635,2.33,2.36/)
      s2=(/.03525,.12386,.1345,.09841,.1373,0.1477,.1444,0.14504924,.07148,.02237/)
      s2g=(/.03525,0.13976128,.1345,.09841,0.13694176,0.13179871,0.1324168,0.13728201,.07148,.02237/)
      s3=(/.00759,.00884,.00521,.00974,.00789,0.00728,.00671,5.9208343E-3,.00224,.00012/)
      s3g=(/.00759,7.79948E-3,.00521,.00974,8.12251E-3,8.32052E-3,7.951397E-3,6.429092E-3,.00224,.00012/)
      s4=(/-0.00206,-.0028,-0.0011,-.00287,-.00252,-0.00235,-.00195,-1.5903986E-3,0.,0./)
      s4g=(/-0.00206,-2.49985E-3,-0.0011,-.00287,-2.6353553E-3,-2.6501498E-3,-2.4560375E-3 ,-1.7370114E-3,0.,0./)
      s5=(/0.19,0.15,0.10,0.15,0.14,0.13,0.12,0.106354214,0.10,0.10/)
      s5g=(/0.19,0.13666,0.10,0.15,0.1429013,0.14334,0.13579784,0.11287035,0.10,0.10/)
      s6= (/0.24,0.23,0.30,0.23,0.253,0.37,0.277,0.34101113,.25,0.25/)
      s6g= (/0.24,0.32338,0.30,0.23,0.30005046,0.27662,0.32512766,0.2953982,.25,0.25/)
c s7 corresponds to site class E 
      s7=(/0.29,0.25,0.55,0.2,0.319,0.38,0.416,0.53479433,0.4,0.36/)
      s7g=(/0.29,0.33671,0.55,0.2,0.29974141,0.29329,0.34473022,0.49243947,0.4,0.36/)
      ssig=(/0.23,0.28,0.34,0.27,0.286,0.29,0.31,0.34,0.34,0.36/)
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
      sigma=ssig(iq)
      else
c intraslab. do not use 10-10-2008. not set up...
      stop ' 10/08: inslab needs additional coeff. for 2.5 hz'
c          g= 10.**(0.301-0.01*xmag0)
c      xmag0=min(xmag,8.0)      !p 1706, limit to 8      
c      depth=50.      !vary depth ?
c      dsq=depth*depth
c      dist= sqrt(dist0*dist0+ dsq)
c      if(ir.eq.1)then
c      gnd0=c1(iq)
c      elseif(ir.gt.2)then
cc      gnd0=c1w(iq)
c      endif
c      co2=c2(iq);co3=c3(iq);co4=c4(iq);
c      co5=c5(iq);co6=c6(iq)
c      r1=rc1;r2=rc2;r3=rc3;r4=rc4
c      sigma=sig(iq)
      endif
          sigmaf= 1.0/sigma/sqrt2/aln10
          delta= 0.00724*(10.**(0.507*xmag0))
          gnd=gnd0+co2*xmag0
          dist2= sqrt(dist*dist + delta*delta)
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
          gnd = gnd + 0.5*sl*(co5+c06) - gfac
c --- CD boundary          
      elseif (vs30.ge.190.)then
c === D site class coeff in c6
          gnd= gnd + sl*co6 - gfac
          elseif(vs30.ge.170.)then
          gnd = gnd +0.5*sl*(c06+c07) -gfac
c--- DE boundary.
          else
          gnd= gnd + sl*co7 - gfac
      endif          
c log base 10 to base e
          gnd= gnd * aln10
c          if(dist0.lt.30.)period,dist,exp(gnd),xmag0,rpga,sl
c      write(15,*) period, xmag, dist, exp(gnd), rpga, sl,xmag0,iq,gnd0
      return
      end subroutine getABsub
      
      subroutine Crouse91(ip,xmag,R,vs30,h,gnd,sigmaf)
c  Crouse earthquake spectra 1991 (p 191-227)
c  his model output is psv (cm/s) and pga (gal= cm/s/s). We convert units to
c g here with coefficient p0 (p0 does not appear in his original model).
c psv -> psa (g) factor is: 2 pi f /980. Additional periods are included by linear
c interpolation of his table 2 values, p. 225, using the logarithm of period (T).
c interpolated values below are evident from number of sig figures.
c input :
c       xmag = moment mag
c      vs30 top 30-m shearwave veloc. Not currently used. But, Crouse91 is a
c soil-only relation. Calling program freezes if Vs30 is too rock-like. Further
c exclusions "soft" soil. Engineers typically equate soft soil to NEHRP E or 
c F, but soft may be for D soil near the DE boundary (180 m/s). NEHRP C could also
c be firm soil. In 1990 these finer soil distinctions were not made.
c	R=rupture distance "distance to center of energy release" p 209.
c we have an uncertainty here, for future events along a long slab. Rrup could
c be close but Rcen*** could be far away.
c	h = depth of top of slab (km). Cascadia models use 20 km. Crouse calls
c	h the focal depth, again a source of uncert because 5 to 20 km are all
c	focal depths that might make sense.
c Output:
c	gnd = log of median motion (Units g)
c 	sigmaf = 1./sigma/sqrt2
c allowed periods are listed in the array period. 0=pga
c
	parameter (p3=0.0,sqrt2=1.4142136)	!p3=the M**2 coeff. dont use this term
	real, dimension(10):: period, p0,p1,p2,p4,p5,p6,p7,sigma
	period = (/0.0,0.1,0.2,0.3,0.5,0.75,1.0,2.0,3.,4.0/)
	p0=(/-6.8876,-2.747,-3.4402,-3.8458,-4.3565,-4.762,-5.0497,-5.7428,
     + -6.1484,-6.43597/)
        p1=(/11.5,3.26,4.44,3.61523,2.93644,2.05331,
     + 1.43,-0.987,-1.67,-2.20/)
	p2=(/0.657, 1.12,1.09,1.14264,1.30657,1.47981,
     + 1.56,1.50,1.59,1.67/)
	p4=(/-2.09,-1.93,-1.92,-1.782,-1.82208,-1.85244,
     + -1.83,-1.38,-1.41,-1.46/)
	p5=(/63.7,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58,1.58/)
	p6=(/0.128,.608,.608,.608,.608,.608,.608,.608,.608,.608/)
	p7=(/0.00397,0.00566,0.00531,4.29220E-03,3.01968E-03,2.24423E-03,
     + 0.00114,-0.0022,-.00367,-.00439/)
	sigma=(/0.6309,0.7376,0.6745,0.65268,0.666606, 0.738604,
     + 0.7483,0.7190,0.8044,0.8099/)
       fac = p4(ip)*alog(R+p5(ip)*exp(p6(ip)*xmag))
       gnd=p0(ip)+p1(ip)+p2(ip)*xmag+fac+p7(ip)*h
       sigmaf=1.0/sigma(ip)/sqrt2
	return
	end subroutine Crouse91
		
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
                                                   
      subroutine zhao(ip,xmag,dist,gnd,sigmaf,ivs, islab,hslab)
c predicted interface g.m. followed Advisory Panel meeting & suggs. Apr 2007
c compute median and sigmaf for the Zhao model with Somerville correction.
c Add 0.75s Nov 18 2008. SH.
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
      parameter (nper=22,hi=20.,sqrt2=1.41421356,xmagc=6.3,gcor=6.88806)
      real, dimension (nper):: a,b,c,d,e,c1,c2,c3,qi,wi,ps,
     & qs,ws,si,sr,ss,ssl,sig,tau,tauI,tauS,per,sigt
      real afac,dist,hfac,site,sigma,sigmaf
c coefficients from Zhao source code not from BSSA tables: more precision in his software.      
c added 0.75s nov 2008. 0.75s (1.33 Hz) is an in-demand T.
	per = (/0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.40,
     &     0.50,0.60,0.70,0.75,0.80,0.90,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0/)
	a=(/1.101,1.076,1.118,1.134,1.147,1.149,1.163,
     &       1.200,1.250,1.293,1.336,1.360956,1.386,1.433,1.479,1.551,1.621,
     &       1.694,1.748,1.759,1.826,1.825/)
	b=(/-0.00564,-0.00671,-0.00787,-0.00722,-0.00659,
     & -0.00590,-0.00520,-0.00422,-0.00338,-0.00282,-0.00258,-2.50E-3,-0.00242,
     & -0.00232,-0.00220,-0.00207,-0.00224,-0.00201,-0.00187,-0.00147,
     & -0.00195,-0.00237/)
	c=(/0.0055,0.0075,0.0090,0.0100,0.0120,0.0140,
     &     0.0150,0.0060,0.0060,0.0030,0.0025,2.3797892E-3,0.0022,0.0020,0.0020,
     &     0.0020,0.0020,0.0025,0.0028,0.0032,0.0040,0.0050/)
	d=(/1.07967,1.05984,1.08274,1.05292,1.01360,
     &   0.96638,0.93427,0.95880,1.00779,1.08773,1.08384,1.0826032,1.08849,
     &   1.10920,1.11474,1.08295,1.09117,1.05492,1.05191,1.02452,
     &   1.04356,1.06518/)
	e=(/0.01412,0.01463,0.01423,0.01509,0.01462,
     &   0.01459,0.01458,0.01257,0.01114,0.01019,0.00979,9.560258E-3,0.00944,
     &   0.00972,0.01005,0.01003,0.00928,0.00833,0.00776,0.00644,
     &   0.00590,0.00510/)
	sr=(/0.2509,0.2513,0.2403,0.2506,0.2601,0.2690,
     &      0.2590,0.2479,0.2470,0.2326,0.2200,2.9499817E-3,0.2321,0.2196,0.2107,
     &      0.2510,0.2483,0.2631,0.2620,0.3066,0.3529,0.2485/)
        si=(/0.0000,0.0000,0.0000,0.0000,0.0000,0.0000,
     &       0.0000,-0.0412,-0.0528,-0.1034,-0.1460,-0.15377375,-0.1638,-0.2062,
     &      -0.2393,-0.2557,-0.3065,-0.3214,-0.3366,-0.3306,-0.3898,
     &      -0.4978/)
	ss=(/0.0557,0.1047,0.1276,0.0780,0.1074,0.0753,
     &        0.0058,-0.0120,-0.0448,-0.0727,-0.1082,-0.11372811,-0.1257,-0.1859,
     &       -0.2268,-0.2370,-0.2400,-0.2332,-0.2804,-0.2305,-0.2548,
     &       -0.3551/)
c sigt = Table 5 total sigma, from intra- and inter-event sigmas sqrt(SS) without the mag-cor term
c	sigt = (/0.723,0.849,0.811,0.77,0.76,0.775,0.779,0.787,0.776,0.751,0.745/)
	sig=(/0.6039,0.6399,0.6936,0.7017,0.6917,
     &  0.6823,0.6696,0.6589,0.6530,0.6527,0.6516,0.6483102,0.6467,0.6525,
     &  0.6570,0.6601,0.6640,0.6694,0.6706,0.6671,0.6468,0.6431/)
c interevnt sig tau from table 5,  Somerville says to use the mag-cor term, but not use the lower sigma_total.	
	tau=(/0.3976,0.4437,0.4903,0.4603,0.4233,0.3908,
     &       0.3790,0.3897,0.3890,0.4014,0.4079,0.41473,0.4183,0.4106,0.4101,
     &       0.4021,0.4076,0.4138,0.4108,0.3961,0.3821,0.3766/)
c intraevent tauI & tauS with the mag-cor term from table 6. see comment p 910,
c col II, Zhao et al.
	tauI=(/0.308,0.343,0.403,0.367,0.328,0.289,0.280,
     &    0.271,0.277,0.296,0.313,0.32331293,0.329,0.324,0.328,0.339,0.352,0.360,
     &    0.356,0.338,0.307,0.272/)
	tauS=(/0.321,0.378,0.420,0.372,0.324,0.294,0.284,
     &    0.278,0.272,0.285,0.290,0.29615661,0.299,0.289,0.286,0.277,0.282,0.300,
     &    0.292,0.274,0.281,0.296/)
c c1 = rock site term, for NEHRP A and B, vs > 600 m/s. SCI from Zhao table 2 p 901 bssa june 2006
c c2 = stiff soil term for NEHRP C (SCII)
c c3 = soft soil term for NEHRP D  (SCIII)
	c1=(/1.1111,1.6845,2.0609,1.9165,1.6688,1.4683,
     &         1.1720,0.6548,0.0713,-0.4288,-0.8656,-1.099344,-1.3250,-1.7322,
     &      -2.1522,-2.9226,-3.5476,-4.4102,-5.0492,-5.4307,-6.1813,
     &      -6.3471/)
	c2=(/1.3440,1.7930,2.1346,2.1680,2.0854,1.9416,
     &   1.6829,1.1271,0.5149,-0.0027,-0.4493,-0.69264763,-0.9284,-1.3490,-1.7757,
     &  -2.5422,-3.1689,-4.0387,-4.6979,-5.0890,-5.8821,-6.0512/)
	c3=(/1.3548,1.7474,2.0311,2.0518,2.0007,1.9407,
     &   1.8083,1.4825,0.9339,0.3936,-0.1109,-0.37189444,-0.6200,-1.0665,-1.5228,
     &  -2.3272,-2.9789,-3.8714,-4.4963,-4.8932,-5.6981,-5.8733/)
c slab event term not used for subduction
	ssl=(/-0.5284,-0.5507,-0.4201,-0.4315,-0.3715,
     &     -0.3601,-0.4505,-0.5061,-0.5538,-0.5746,-0.5721,-0.5563261,-0.5397,
     &     -0.5216,-0.5094,-0.4692,-0.3787,-0.2484,-0.2215,-0.2625,
     &     -0.1689,-0.1201/)
c qi and wi are for the magnitude square term correction
	qi =(/0.0,0.0,0.0,-0.0138,-0.0256,-0.0348,-0.0423,
     &        -0.0541,-0.0632,-0.0707,-0.0771,-0.07988431,-0.0825,-0.0874,-0.0917,
     &        -0.1009,-0.1083,-0.1202,-0.1293,-0.1368,-0.1486,-0.1578/)
	wi =(/0.0,0.0,0.0,0.0286,0.0352,0.0403,0.0445,
     &    0.0511,0.0562,0.0604,0.0639,0.065498955,0.0670,0.0697,0.0721,0.0772,
     &    0.0814,0.0880,0.0931,0.0972,0.1038,0.1090/)
c Ps, Qs, and Ws for deep intraplate eqs; see last several columns of Table 6
        ps= (/0.1392,0.1636,0.1690,0.1669,0.1631,0.1588,
     &      0.1544,0.1460,0.1381,0.1307,0.1239,0.12070625,0.1176,0.1116,0.1060,
     &      0.0933,0.0821,0.0628,0.0465,0.0322,0.0083,-0.0117/)
        qs= (/0.1584,0.1932,0.2057,0.1984,0.1856,0.1714,
     &      0.1573,0.1309,0.1078,0.0878,0.0705,0.06278851,0.0556,0.0426,0.0314,
     &      0.0093,-0.0062,-0.0235,-0.0287,-0.0261,-0.0065,0.0246/)
        ws= (/-0.0529,-0.0841,-0.0877,-0.0773,-0.0644,
     &       -0.0515,-0.0395,-0.0183,-0.0008,0.0136,0.0254,0.030538963,0.0352,
     &        0.0432,0.0498,0.0612,0.0674,0.0692,0.0622,0.0496,0.0150,
     &       -0.0268/)
       sigt=sig
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
      gnd= gnd - gcor
      sigmaf = 1./sigma/sqrt2
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

