c--- program  hazFXnga13l.f; 03/20/2013; Use  with NGA relations, or others.
c 3/19/2013 : increase nfltmx to 1000.
c 3/13/2013: revise AS13 to v10. The subroutine with v10 is according to Norm closer to AS08 than the v5 version.
c 3/14/2013: revise AS08 to use our Rx. 3/15:Use vs30_class=1 (measured). Index 16
C 3/04/2013: working on BSSA13. Mar 11. Running. Simplified the coding some. Perhaps more needed (SH).
c 3/06/2013: working on CY2013. Dirctivity, cDPP, is turned off for initial exercise.
c       A long document discusses some suggested ways to implement...
c 3/01/2013: update A&S GMPE. 3/04/2013: The A&S subr. is expecting z1 to come in with units km/s not m/s.
c
c 2/25/2013: Update CB-NGA2 GMPE Index is 34
c 2/19/2013: Add option to output inelastic spectral displacement, SD_i. 
c This option is invoked if the second argument is "dy" or "dY" and the third is the
c value of dy to use (units: cm). Expected use, with WUS crustal-src GMPEs.
c 1/07/2013: can use a fixed aleatory sigma for all NGA-W gmpes.
c To invoke this option:
c      hazFXnga13l.test input.file fixsigma 0.6
c That is make arg(2) 'fixsigma' and make arg(3) the value you want (natural log sigma)
c 12/17/2012: update the CB2012 input file to "v2". Update affects the tau part of sigma_aleatory.
c 12/11/2012: add and repair the AS2012 GMPE for NGA-W. Use index 36. Implementation may have the Ry0 term.
c     To get the Ry0 version, set icode(ip,ia)=1 in input line
c 12/04/2012: add the latest Graizer Kalkan gmpe. Includes basin effect. See dgkbasin. Index 38
c
c 11/28/2012: add the BSSA 2012 GMPE: index 33
c 12/03/2012: CB12 with v3 update. Previously CB12 had v2 update. Index 34. (removed 2/25/2013)
c 12/03/2012: Idriss 2012 relation added index 37.
c 11/14/2012: the Grazier GMPE is included for use in tectonically active areas. Index 39.
c 11/20/2012: the CY2012 GMPE is included for use in tectonically active areas. Index 35. Use the
c   sigma model which assumes Vs30 is known. From Chiou email 11/20/2012
c 3/21/2011: add 3 2011 GMPES for CENA from G. Atkinson email. indexes 25, 26, 27.
c These 2011 CENA GMPEs are set up for hard rock and BC(firm rock) only. Don't use them for the soil site condition.
c 10/30/2012: use Pezeshk11 sigma from article. not Gail's "stopgap". 
c 11/01/2012: repair index 
c      	problem when accessing sigPez11. 
c 10/30/2012: Add Silva single-corner, m-dependent stress, index 28. BC or A
c FOr variable Vs30 do not use this code. Instead use hazFXnga7vs.f.
c from hazFXnga7x.f. This version has event-clustering option. If this option is
c selected, all events in input file must be part of the cluster model. Initial
c work on this was for NMSZ scenarios. These scenarios dictate several limits
c on the available model. Do not mix unclustered sources (faults) with clustered
c in any given input file,
c but you can  use this code on either clustered or independent events. Clustered events
c are set up for Characteristic without uncertainty. Other rup. models need more work. 
c--- Solaris compile: cc -c iosubs.c 
c                 f95 -o hazFXnga13l hazFXnga13l.f iosubs_128.o -fast -e -ftrap=%none
c for debug:
c  f95 -o hazFXnga13l hazFXnga13l.f iosubs_128.o -C  -e -ftrap=%none
c for typical runs:
c  f95 -o hazFXnga13l hazFXnga13l.f iosubs_128.o -O  -e -ftrap=%none
c -e to extend linelength. Several statements go past column 72. Flag -fast no longer works due to some
c excess stuff in the BSSA initialization. Best might be -O 
c The ftrap flag is required to prevent crashing when NaNs such as log(0) are
c attempted. Code has been modified to avoid some of these where known.
c for detailed source listing -Xlist is the argument. -Xlistf list no compile
c
c Linux compile, with Intel fortran compiler:
c  ifort -o hazFXnga13l hazFXnga13l.f iosubs.o -O -132
c Sometimes it is necessary to use a static link flag:
c  ifort -o hazFXnga13l hazFXnga13l.f iosubs.o -O -132 -static
c
c Linux OS and older versions of gfortran (thru 2006) are a bad combination 
c   when trying to run this pgm. 
c  Note: on cluster, if using Intel version must use flag -V:  qsub -V job
c
c PC, windows, if you use gnu fortran:
c   gfortran -static -O -ffixed-line-length-none -ffpe-trap= -o hazFXnga13l hazFXnga13l.f iosubs.o
c
c 8/28/2007: add option GR branches with b=0 wherever there is currently a branch with some other b. 
c This new set of branches is weighed equally with the old ones. For example, if a fault
c GR branch on input has b=0.8 and weight 0.25, 
c on output the b=.8 branch has weight 0.125 and the b=0 branch
c also has weight 0.125. This option is invoked by making itype(ift)=-2, previously an illegal value.
c The effect is to slightly decrease the rate of M6.5 eqs in the model and slightly increase the
c rate of M7.0 eqs in the model. The M6.5 rate model bulge is slightly reduced by trying this.
c
c 10/30/2007: deaggregation option available. Just add site Lat Long nper sa(j),j=1,...,nper
c These fixed sa values should be in same order as the period set in input file. This
c deagg work was not completed for clustered events or for CEUS. Should work for WUS fault files.
c More work is needed for the CEUS or anyplace with clamped motions or  clustered-event option
c
c 10/13/2009: correct TP05 c1(0.3 s) for BC rock.
c 11/12/2008: initial ymax for set of stations. Thanks R. Dobbs.
c july 1 2008: add 2.5 and 25 hz to getSilva
c 5/05/2008: in resample_fault deltax is now defined as deltaz*cot(dip). This mod affects GR scenarios for some
c steep but non-vertical dipping faults such as Imperial flt, southern CA-Mexico.
c 5/01/2008: fix the clustered-event output for 1 station, every curve for a given SA outputs to 1 file.
c 3/20/2008: update BA-NGA relation to compute pga_nl using the final PGA coeffs, which
c include mechanism dependency, for example. Follows a new pdf file from Dave Boore, 3/19/2008.
c This modification does not affect the BC rock (reference site condition) calculations. 
c 3/06/2008: add some logic to allow input fault dips to be > 90 or < -90 degrees (K.Campbell sugg);
c      also check for VS30 > 1500 and stop if using Campbell NGA relation.
c 2/22/2008: GR source model epistemicM uncertainty is no longer removed if sdal is 0.
c Note: For the single magnitude GR case you must include nonzero sdal if dM_epistemic is non-zero
c
c 2/21/2008: Only changes are to comment lines (updated) and cosmetic format changes
c 1/08/2008: uncommented out a required statement,  cs=1.88 in getCBNGA1107
c 12/07/2007: add median clamps to getSilva and getTP05 (inadvertent omission prior to)
c 12/05/2007: the irab() array is now 2-dimensional (was 1-d before). to get proper AB06, need 2-D.
c 11/26/2007: AB2006 (CEUS) nonlinear S term is modified for Vs<760. Does not affect BC rock
c 11/19/2007: Campbell sigma terms are modified slightly for nonlinear soil. These mods
c can affect pga when Vs30 = 760 because Campbell k1 is 865 m/s, and nonlinearity kicks in
c when Vs < k1. However, the effect is quite minute for PGA of rock at BC boundary.
c 11/07: Chiou makes another small correction to his tau-term. Increases tau very slightly.
c 11/06: Chiou makes a small correction to his tau-term. Reduces tau very slightly.
c 11/02: use very latest Chiou-youngs model with a new dlt1 effect in the sigma calc.
c 11/01: Use Chiou M-dependent sigma 
c 10/25/2007: Add the new CY2007 NGA model (B. Chiou emailed this), 
c               which uses new R_x signed distance, now computed
c      in the revised mindist1 (Y. Zeng)
c 10/17/2007: Invoke the downdip rupture option only if top of fault within 1 km of surface
c      of Earth. Blind thrusts and so on will not have additional rupture tops.
c 10.12.2007: improve log file output info for b=0 branch. Show crossover M.
c 10/03/2007: Update Boore Atkinson atten with Apr 02, 2007 coeffs.
c
c 8/28/2007: make the dMa on aleatory a function of the sigma. dMa previously was fixed at 0.05, now
c       dM is 2sigma_ale / 5 = 0.4 sigma_ale The Mmin=5.8 low limit has been removed.
c
c 8/13/2007: add aleatory dM with sigma<0 to keep event rate fixed at input value.
c      if sigma>0, moment rate is fixed.
c july 20 2007: add deterministic source GND output option: invoked if nrec<0
c august 3 2007: modify output of determinstic GR to only show the nearest src
c      i.e., remove all but one of the aleatory rupture location scenarios
c May 16: add scenario and geographic index to the output 
c
c April 23 2007: add polygon feature. Only include a fault source if at least
c      one point on fault when projected to earth surface is 
c      interior/on a polygon surface. The polygon-checking
c      feature is optional, and is included if the polygon file name
c      is included as arg 2 in the command line input stream.
c      May 25: this check was refined to include check whether individual
c       segments are inbounds or not for GR segments or floating ruptures. In
c      other words a fault can be inbounds but several floating segments on that
c      fault may be out-of-bounds. When determining event rates, this distinction
c      can be important for many relatively long WUS faults.
c Clustered source model has been checked for characteristic-without-uncertainty-in-M
c types of ruptures only. Future work, add aleatory uncert to clustered source model.
c
c late March: use Zeng's mindist1 algorithm for all rjb, rcd dist calculations.
c 8/12/2009: add getDahle95 for Panama hazard assessments.
c 3/19/2007: use Frankel HR->FR factors for TP05 at 7 periods (this still needs work at other
c      periods). Somerville HR revised slightly.
c 3/13/2007: revise BA NGA relation, now has 21 spectral periods. see feb 14 doc.
c 1/27/2007: add/subtr. additional epistemic uncertainty to all median relations with the uncert
c in a 3x3xnp matrix gnd_ep(3,3,npmx) if l_gnd_ep is true. Advisory panel
c has recommended smoothing this gnd_ep matrix in some way. 
c 
c Mod Nov 7 to report rates of events when the "list of stations" option is
c invoked. For each station or receiver, a file is written, eqrate.sta_name.
c  NGA7: This code has potential downdip ruptures with tops at 2 and 4 km depth
c when GR-distributed sources and 6.5<=M<=7.0. Testing is nearly complete. 
c  Added GR with 6 < M < 6.5. For these, rupture tops are equally distributed 
c               at ztor+0,2,4,and 6 km. Revised so that AspectRatio>=1 for all rups,
c            even those associated with M6.0 sources.
c What's in a name?
c   Some of the downdip rupture branching is triggered or not triggered by
c a match to part of a filename. This feature definitely needs a more failsafe
c logic approach. For example, filename containing 'aFault_unseg' is a trigger to
c not include downdip rupture scenarios even though these are type 2 faults (GR)
c Because different rules are being used in different regions, it might
c be possible to make rule sets depending on state boundaries or other
c geopolitical features. This has not been done as of Feb 2008. SHarmsen.
c Jan 2007: increase or decrease in median is possible with amount of increase
c contained in a 3 x 3 x nper matrix, gnd_ep. If you do not want to modify
c gnd, make the indicator variable wind 0 for that period. All attenuation models will
c be affected. We are planning to use this only with NGA relations, because there
c is already enough variety of ideas about the median in CEUS, intraplate, and
c so on. If this option is on, three curves per site are computed and output,
c one for median, one for median+augment, and one for median - augment. They
c go to different files. This option can be turned on for some spectral periods,
c off for others. If on, and not NGA, the .p and .m files will be zero (10**-21)
c
c--- to run: hazFXnga13l inputfile > report.log &
c ---  Alternate:  hazFXnga13l inputfile polygon.file > report.log &
c --- 2nd Alternate: hazFXnga13l inputfile dy 0.5 > report.log &
c The second alternate computes inelastic response spectrum with dy=0.5 cm.
c Deaggregation example (new Oct 2007): 
c
c  hazFXnga13l orwa.char  44.4 -123.4 3 .12 .34 .14 
c A companion deagg code now exists, called deaggFLTH.2011.f. Used at web site
c for public deagg. requests. 
c
c The first arg after pgm name is the file name, 2nd is lat, third is long, 4th (3) nper,
c and 5 thru 7 are PGA, 0.2s SA, and 1s SA, resp. 
c PGA,  0.2s  and 1s SA are the order of T in orwa.char in this example.
c The command-line information replaces (supercedes) the corresponding info in
c the input file. For example, input file might have 20 sa levels, but there is
c always just one sa level per period for deaggregation, which is input on command line.
c Deagg bins: dM is 0.1, dR is variable, 10 km inner 100, 25 km for awhile, then a catchall.
c      Epsilon bins are 1 sigma in size with edges at k*sigma, k=0, +-1,+-2
c      Epsilon binning might be considered coarse. Follows 2002 1996 web sites.
c  Command-line info follows the interactive deagg web site precedents of 2002 and 1996.
c This code should be compatible with certain Linux Cluster as well as Solaris
c compilers. Solaris version: -e in gfortran is fixed-line-length-none. 
c Linux: gfortran -o hazFXnga13l hazFXnga13l.f iosubs.o -ffixed-line-length-none
c   linux compiler is accepting this code Oct 13 2006. 
c "array constructors" are modern replacement for "data" statements. This is the
c main modification to make code gfortran-compatible. SHarmsen
c
c new 6/06:  sources now can have up to 12 magnitude branches per fault.
c For each mag there is a specific cmag and crate. If no epistemic uncert on
c magnitude-frequency is modeled for a given fault, ift line must
c have a 1 for number of mags per fault (this means all 2002-era input files must be modified)
c Has similar mods for GR-distribution of sources. This required dimension increase on
c several arrays including xmo2
c  
c Modified Mar 8 2006 to improve identification of HW sites.
c Also, transfer out of labor-intensive downdip calcs if appropriate. This change
c may bypass the Frankel abrahamson-hw part of code which is not used with
c NGA. If this code is to work with old subroutines, need to revisit...
c
c Some older relations are in the nga-style subroutine format with coeffs 
c      explicitly included as statements. These are:
c iatten=0 Truth. This one hasn't been programmed. It is the oldest but least accessible
c iatten=1 Spudich ea, 2000 Extensional-tectonics regions, with BJF97 site amp
c iatten=2 Toro ea; CEUS 7 periods and BC rock geotech. M has to be Mw. finite fault
c iatten=-2 Toro ea; CEUS 7 periods and A rock geotech. New 11/06: finite fault
c iatten=3 Sadigh 7 periods rock
c iatten=-3 Somerville IMW. 5 periods available. No site response var. add jan19 2007.
c iatten=4 or -4. AB06, -4 for hardrock (is overridden if Vs30<1500 m/s). Mod Mar 14 2008. SH.
c iatten=6 Fea ceus 7 periods BC rock
c iatten=-6 Fea ceus 7 periods A rock
c iatten=5 AB94 BC ceus (New index, previously had same 6 index as fea. Too confusing)
c iatten=-5 AB94 A ceus (New, hard rock only)
c iatten=7 Somerville finite, BC rock. ceus
c iatten=-7 Somerville finite, A rock. Will get -7 if Vs30>= 1500. ceus.
c    Somerville distance modified to rjb Apr 7 2009. K. Campbell noticed previous
c      use of rrup in the getSomer calls.
c iatten=8 Abrahamson & Silva 7 periods rock. Using the logical variable "hwsite" to
c       note for A&S: increase ground motion when hwsite is true. THis is a quick fix and does not include
c       the 22 degree wings.
c iatten= 9 Campbell and Bozorgnia 2003. 7 periods BCrock. july 27 2006
c   added iatten=10 Campbell CEUS. bc rock. july 29
c   added iatten=-10 Campbell CEUS.  hardrock. 
c iatten=19 Tabakoli and Pezeshk can be firm rock or hard rock.
c  BACK TO WUS atten models:
c iatten=11 BJF97, any vs30. 7 periods. july 27.
c Important note, if any of the above require a different than moment mag, must
c do the conversion outside of the subroutine. icode=1 do something before calling...
c   iatten=12 Motazetti and Atkinson, for PR/VI. Limited period set.
c
c ---  New relations: 
c --- Include the latest NGA developer-team updates here asap. 
c --- Iatten=13 for B&A 4/07; to 10-s SA maximum T, has PGV model.
c For all relations, PGV output units: cm/s.
c B&A April update includes 7.5 and 10-s SA predictions; coeffs change from Feb
c for some T, no change for 1s or 0.2-s SA and the reference 760 m/s site class.
c --- Iatten=14 C&B, 11/10/06. Includes PGV&PGD models. Sigma for random H-comp.
c              C&B hanging-wall term (f4) applies to both Normal and Reverse.
c --- Iatten=15 C&Y. This is model of Oct 2007. No PGV for C&Y, but 105 T-vals  (rev Oct 2007)
c --- Iatten=35 C&Y. This is model of Nov 2012. No PGV for C&Y, just 3 pds CYpd
c --- 
c --- 
c --- Iatten=16 for A&S 3/08; added Mar 13 2013. With Sanaz mods HW taper. 108 per 
c ********* 5
c --- Iatten=17 Idriss. Added nov 2005. Note Sept 2006: Idriss has coeffs for several
c --- iatten=18 kanno
c --- iatten=22 Dahle 1995 added aug 12 2009. Has  a rock or soil flag, 2 state response.
c --- iatten=20. Silva, Gregor, Darragh 2002. Hard rock and BC rock for CEUS. New jan 30 2007. 8 periods
c --- iatten=21 BA06, CEUS, 200 bar . should work with any reasonable vs30. (S term may or may
c     not be present. Hardrock coeffs used  if Vs30>=2000 m/s, in which case S is 0.0)
c --- iatten=-21 BA06, CEUS, 200 bar. Specifically request hardrock coeffs even though Vs30 is in
c      the 1500 to 2000 m/s range.
c New Mar 18 2011
c      iatten=25 Atkinson 08 CENA distance attn. based on Rjb
c      iatten=26 Atkinson 06 CENA distance attn. based on Rcd
c      iatten=27 Pezeshk11 CENA distance attn. based on Rcd
c       iatten=28 Silva (variable stress, with M dependence). uses Rjb. Large-amp median. be careful.
c       iatten=38 Graizer and Kalkan, WUS attenuation 2012 update. Distance Rcd
c ***   iatten=39 Graizer and Kalkan, WUS attenuation 2009 model. Not using basin
c term. 
c    
c ***
c
c *** New Oct 27 2005: Enter vs30 after the site description. Also enter depth of
c *** basin 
c *** dbasin= 0 for surficial hardrock; Depth to 2.5 km/s VS in
c *** Campbell&Bozorgnia (Nov 2007 NGAmodels) *** Also used in Graizer&Kalkan 2012
c *** what about firm rock basement? Z1 value. New default Z1=.025 * dbasin Oct 2007, converted to m.
c ***
c *** Towards parallel processing: Y. Zeng has developed a separate code for
c *** parallel processing. It is now used on gldlabm, a PC Cluster Linux OS

c       iosubs.o:       ELF 64-bit MSB relocatable SPARCV9 Version 1
c
c --- Nov 2005: New Feature corresponding to feature of hazgridXnga2:
c --- code works for a large grid of sites or a small set (<=30) of sites.
c ---   Input file first line is number of sites, nrec.
c --- Make nrec 0 to perform the calculations for a grid of sites.
c --- Make nrec 1 to 30 to specify set of sites <lat(i),long(i),i=1,...,nrec>
c If you use small station set, output includes eq rates*model weight files,
c one for each station.
c --- Output file names are specified in input file, one per period.
c --- This version calls atten. subroutine each time the mean,sd are needed.
c--- uses NGA  subroutines for attenuation relations 
c--- Important: if sdal is zero, the epistemic dM branches will be ignored and
c --- only the central branch will be used (dM=0). This is a hard to remember
c --- feature of this code. 
c----  HazFX originally written by Art Frankel, with many modifications & updates 
c        by Steve Harmsen, USGS, Golden,CO
c---- averages strike of fault segment. dip direction is from an avg strike
c---- Frankel revised hanging wall test for AS97, including 22.5 degree wings.
c ---   the 22.5 degree wings are no longer supported. SH 2007.
c-- Can use up to 8 attenuation relations for each period. Atten model weights should sum to 1.
c --- This vers. does not save atten-model pex info in arrays. SH 
c---- added getMota 10/16/02
c---- revised Campbell2002
c---- added hanging wall term for Campbell 2002
c---- hanging wall term, separate reverse from thrust faulting
c---- fixed hanging wall for floating rupture
c---- changed to dmove increment for floating rupture 9/30/02
c---- added getBJF97 9/30/02
c--- changed length to WC all 7/1/02
c--- changed to only epistemic uncertainty for GR faults
c----- changed to Mmin= 5.8  nmag= 55
c---- added clamp for CEUS   max gnd motions
c---- added Mchar epistemic uncertainty 10/16/01
c---- added Mchar aleatory uncertainty 10/16/01
c---- checks to see if distribution goes below M6.5 for GR
c---- or below 6.0 for char.
c---- if so, ignores uncertainties
cc----- truncates at 3 sigma
c--- with clipping for Toro and Boore new table
c--- calculates mean annual frequency of exceedance for dipping faults
c--- Orig version was written by Art Frankel 4/95, modified 10/95
c--- Many updates since then for the 2002 and 2007 NSHMP updates.
c-----  
c--- Output of this program is input into hazallXL.v2 to get probabilistic
c    ground motions
c--- PGA & SA ground motion levels should be in units of g
c--- period=0 indicates PGA
c --- PGV ground motions are in units cm/s (BA and CB as of feb07)
c --- PGD units cm (CB only, oct 2006)
      parameter ( pi=3.14159265,fourpisq=39.4784175,tiny=1e-12)
c Prob. calcs associated with fault rates less than tiny will be omitted.
c Some of the Aug 2007 AFault segmodels have rates of 0 and of 1e-20 for example.
      parameter (nfltmx=1000, nlvmx=20, npmx=8)  
      REAL, PARAMETER :: vref = 760.0
      dimension xlev(nlvmx,npmx),ylev(nlvmx,npmx),nlev(npmx),icode(npmx,8)
      dimension p(250005)      !table of complementary normal probab.
c You can raise p dim & make some minor code changes to improve accuracy. Currently
c p() has about 4 or 5 decimal place accuracy which is likely good enough.
      real, dimension(nfltmx):: xaz,xlen
      real, dimension(800,nfltmx):: x,y
      real prob_s(3,8,5,20,8)
      real, dimension(8):: dmin,rate_cl,pdSilva2
      real, dimension(11) :: pdSilva
c wind used to determine if additional epistemic uncert will be added to
c or subtracted from log(median) of each relation for that period
      real, dimension (6,6) :: dminr
c      real, dimension (8,8) :: nga13l
      real, dimension (3,3,8) :: gnd_ep
      real, dimension(3) :: gnd,gmwt
c gnd(1) contains the median; gnd(2) contains median+gnd_ep(); and gnd(3)
c      contains median - gnd_ep(). Ultimately, curves go to different files,
c One goes to output name given, two, to output//.p and three, to output//.m
      real, dimension(2) :: mcut,dcut,tarray
      real, dimension(3) :: CYpd
        real, dimension(32) :: slat,slong 
        real, dimension(100) :: ale, ale2
c  array of station coordinates if using list option
      integer, dimension(nfltmx):: igroup,jseg,nmagf,itest,mx_index
      integer, dimension (50,nfltmx) :: nrupd,nrups
c add some arrays for storing information related to downdip ruptures
      real rupln(50,nfltmx),dmin2(20,1000,6,6),wtscene(8,5)
      integer nmag0(nfltmx,12),immax/1/,indtmp(npmx),Vs30_class

      real ratenew(nfltmx,12,6),xmo2(nfltmx,12),relwt(nfltmx,12)
      real aperiod(22),camper(25)      !abrhamson, campbell-boz period sets resp. Code pgv=-1
       real perb(23),abper(26),tpper(16),xdiff,ydiff
c abper atkinson boore ceus2006 tp=tabakoli & pezeshk
c Deagg array storage: rbar,mbar,ebar,haz indexed by R,M,epsilon,GMuncert, and spectral period
      real, dimension (16,30,10,5,npmx):: rbar,mbar,ebar,haz
c Fault deagg array storage, indexed by fault number (ift), GMuncert, and spectral period
      real, dimension (300,5,npmx):: frbar,fmbar,febar,fhaz
      real prob5(16,30,10,npmx,5,5)      !30 possible M bins M5.8 to M8.7 by 0.1 dM. 4megabyte storage
      integer ip,nft_out,nrupdd,ii,igpmax/0/,isoild
c isoild is a soil indicator for the Dahle95 relation, 1 if Vs30<550 m/s.
c           
      character*12 sname(32)    !station names might be useful.
      character*48 polygon
      character*12, dimension(-10:39) :: att_typ     
      real ar(400,6)      !aspect ratio for floating ruptures. Can have 5 downdip
      logical ss,rev,normal,obl,pgacalc,nearsrc,noas97,l_gnd_ep(8)      
      logical lb0,cluster,grid,isok,isbig,isclose,norpt(8,8),sdi/.false./
      logical lxyin,segin,poly,slist,openme/.true./,override_vs
c     grid = true if stations form a regular grid
c firstf is a logical variable to instruct code to perform certain fault calcs
c the first time the fault is worked on. When many stations are present, it makes
c sense to economize on these repeated fault calcs. New Nov 15 2006.
c Some of these laborious calcs only necessary for nearby sources. Logical
c variable nearsrc is true if nearest distance to current fault < 30 km.
c Variable poly is true if arg2 is polygon  file.
c fault sources will be included only if f_is_in(ift)=.true. 
c norpt(ip,ia) is true if no report has been filed for fault(ift) (deterministic
c      calcs only)
       logical, dimension(nfltmx):: firstf,f_is_in,nodowndip
c type replaces structure & requires gfortran or f95 to compile. It is included
c so that this code has compatibility with Linux cluster-computer compilers.
      type header_type
        character*128 :: name(6)
        real*4 :: period
        integer*4 :: nlev
        real*4 :: xlev(20)
        real*4 :: extra(10)
      end type header_type
      type(header_type) :: headr, hdvs30
c hdvs30 might be used if variable vs30 is read in.
      real magmin(nfltmx,12),magmax(nfltmx,12), mmax,clamp(8),wts
c clamp can limit the probabilistic motion in CEUS. THis is applied in main
c in hazFXnga13l (Clamp applied in CEUS only. Set clamp(i) to 0 to skip this constraint).
c
      real z1,z1km	!depth to sustained 1 km/s Vs . Z1 is defaulted to a CY2013 function of Vs30.
      real v30(100000)      !possible array of site-specific vs30 new mar 2006.
      real dbasin,vs30,pr,rx,ry,dy_sdi
      real, dimension (5) :: prdsom
      integer nmag,readn,ix,ix0,ix1,iy,iy0,jsegmin/10/,jsegmax/0/
      character nameout*78,name*60,adum*80,date*8,time*10,zone*5
      character*12 pithy(3),g_name(5),cmt
c cluster group names, currently g1,g2,g3,g4,g5. Could be more geographic such as w1, w2, c0, e1, e2.
c Could be input
      character*3 gname(5)      !grouped cluster file-name extensions, .g1, .g2, ... 
       common/gail3/freq(13)
       common/fix_sigma/fix_sigma,sigma_fx
       common/sdi/sdi,dy_sdi,fac_sde
       common/perris/sigma_down,sigt_fac
      common/mech/ss,rev,normal,obl
      common/soils/vs30
      common/dipinf/dipang,cosDELTA,cdipsq,cyhwfac,cbhwfac
c pgacalc is a logical variable, true if first period in set is calculated
c hardrock is true if Vs30 >=2000 m/s or greater (why not 1500 m/s)
       common/ceus_sig/lceus_sigma,ceus_sigma
      common/hardrock/hardrock
      common/pgac/pgacalc
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/fault/x,y,u,v
      common/cb13p/Percb13
      real*8 dp,dp2,pr0,prl,prlr,prr,plim
      real*8 emax/3.0/, sqrt2/1.4142135623730951/
       real, dimension(8) :: fac_sde
      real dum(5)      !for bssa nov 2012
      real a10,xmob0,xmorate,sum,alet/0.0/,sigt_fac,lnSa
      logical fix_sigma/.false./,byeca,byeext,byeoreg,hwsite,v30a,byenv,cal_fl,determ,bfault,callme(3),
     1 deagg/.false./, perris/.false./,sigma_down,lceus_sigma/.false./,readbssa/.true./,
     2 readcb13/.true./,useRy0(10)
c readbssa logical to read in bssa coeffs for NGA-W 2012 revision
c readcb13  logical to read in CB12 coeffs for NGA-W 2012 revision
c determ = true if deterministic model: make nrec -100 for gridded smaple or -n,
c for n=1,30 for station set.
c cal_fl: true if  California floating eq. If true use length based on xmag(area)
c for california faults Hanks? or Ellsworth?
c
c hwsite: once ascertained that a site is on hanging wall many calcs can be bypassed.
c for footwall, same is true but code not altered for this case.
      real, dimension (990,4,nfltmx) :: u,v      
      real, dimension (10) :: dmbranch,wtbranch,grp_rate
c not using testhw3(1000,50)
      real prd(106),prob(1000,20,npmx,3,5),out(800000),wt(8,8,2),wtdist(8,8)
      real, dimension(80):: px,py 
      real, dimension(0:35) :: pdgk      !graizer and kalkan
      integer ival(8),icc(8,5)
      integer, dimension(8,8):: iatten,irab
c ipertp = map to tabakoli-pezeshsk;iperb = map from input file to boore set of per
c irab is an index for 4 possible varieties of the AB06 model. 2 for 140 bar & 2 for 200 bar stress
c irab became 2-D on December 5, 2007. S Harmsen
      integer, dimension(npmx) :: ipercy,iper,nattn,ipera,iperb,icy13,ipcb12,
     + iperab,ipercb,iperk,ipertp,isilva,isilva2,isomer,ngroup,nfi,iperdahl,
     + ia06,ia08,ip11,ipgk,indbssa,idriss,ipas13 !new models Mar 2011 to dec 2012
c      real, dimension (22):: Percb13,PerIMIdriss
      real, dimension (23):: Percb13
      real, dimension (22):: PerIMIdriss
      real, dimension (23):: perAS13
      integer ifn,isz,ipia,nscene
      integer, dimension(npmx,5,3) :: ifp
      integer, dimension (nfltmx) :: itype,npts,npts1,iftype,ibtype
      logical, dimension(npmx,8) :: nga,wus02,ceus02,ceus11     
c logical variables for subsets of attenuation models should help narrow
c the search more efficiently.   CEUS11 added mar 18 2011: Gail Atkinson's 3 new
C CENA models, with indexes 25, 26, and 27. 
c Benioff or deep-seismicity relations n/a here: see gridded hazard code hazgridXnga2.f      
c new 6/06: potentially variable a and b values for up to 12 branches for each fault
c Provides extra flexibility to model epistemic uncertainty of size of eq. 
c dmag also needs a branch-specific
c value. The customary 0.1 dM wont work when M precision is carried to 2 dec. places
      real, dimension(13) :: a11fr    !,a11per        Atkinson 2011. frequency set 99=pga, 89=pgv
      real, dimension (0:37) :: perka
      real, dimension (24) :: percy13
      real, dimension (108) :: peras08
      real, dimension (23):: perbssa13
      real, dimension (6) :: ptail,perdahl
      real, dimension(npmx):: period,perx,wtsum,safix
      real, dimension(nfltmx,12) :: a,b,dmag
      real, dimension(nfltmx,12):: cmag,crate
      real, dimension(nfltmx):: dip,dip0,width,depth0,tlen,cDPP
c      integer iargc,numarg,hwflag,vs30_class
      integer iargc,numarg,hwflag
c Below, the array constructor business. No repeat (8*.01) unlike data statement.
c gfortran required replacement of "data" statements in our Linux PC system. 10/06. SH.
c abper is the set of spectral periods for the AB2006 CEUS model. -1 => PGV
      abper = (/5.0000, 4.0000, 3.1250, 2.5000, 2.0000, 1.5873, 1.2500, 1.0000,
     1 0.7937, 0.6289, 0.5000, 0.3968, 0.3, 0.2506, 0.2000, 0.1580,
     1 0.1255, 0.1, 0.0791, 0.0629, 0.0499, 0.0396, 0.0315, 0.0250,
     1 0.0000, -1.000/)
      perbssa13=(/-1.0, 0.0, 0.010000, 0.020000, 0.030000, 0.040000, 0.050000
     +, 0.075000, 0.10, 0.150000, 0.20, 0.250000, 0.30, 0.40, 0.50
     +, 0.750000, 1.0, 1.50, 2.0, 3.0, 4.0, 5.0,10.0/)
       a11fr =(/0.20,    0.33, 0.50, 1.00, 2.00, 3.33, 5.00,10.00,20.,33.00,50.00,99.00,89.00/)
c Somerville IMW period set (5 of 'em). pga is 0.0
       prdsom=(/0.000,0.200,0.300,1.000,5.000/)
c Abrahamson Silva 2012 model 22 perios
      data perAS13 /0, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3,
     1             0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6., 7.5, 10., -1. /
c Silva period set.
      pdSilva=(/0.,0.04,0.05,0.1,0.2,0.3,0.4,0.5,1.,2.,5./)
       pdSilva2=(/0.,0.1,0.2,0.3,0.5,1.,2.,5./)
c Tavakoli periods 0 = pga. added 0.04 & 0.4 s july 8 2008 (spline interpolation)
      tpper = (/0.0,.04, 0.05, 0.08, 1.00e-01,1.50e-01,2.00e-01,
     1 0.3, 0.40, 0.5, 0.75, 1.0, 1.50, 2.0, 3.0, 4.0/)
c perka = period set for Kanno et al., BSSA 2006. 0 = pga.
        perdahl=(/0.0,0.1,0.2,0.5,1.,2.0/)
       perka=(/0.0,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.15,0.17,0.20,
     +0.22,0.25,0.30,0.35,0.40,0.45,0.50,0.60,0.70,0.80,0.90,1.00,1.10,1.20,
     +1.30,1.50,1.70,2.00,2.20,2.50,3.00,3.50,4.00,4.50,5.00/)
c perb = Boore-Atkinson NGA 4/2007 period set, -1 = pgv. N(T)= 23. 10 s is longest.
      perb= (/-1.000, 0.000, 0.010, 0.020, 0.030, 0.050, 0.075, 0.100,
     + 0.150, 0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000,
     + 1.500, 2.000, 3.000, 4.000, 5.0, 7.5, 10.0/)
c perx = a default set of spectral periods used in 2002 PSHA maps. PGV=-1 was not
c available then and should not be tried for those models
       perx = (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
       percy13=     (/ -1.,     -2.,      0.0100, 0.0200, 0.0300,
     1                0.0400, 0.0500, 0.0750, 0.1000, 0.1500,
     1                0.2000, 0.2500, 0.3000, 0.4000, 0.5000,
     1                0.7500, 1.0000, 1.5000, 2.0000, 3.0000,
     1                4.0000, 5.0000, 7.5000,10.0000/)
       gmwt = (/0.63, 0.185, 0.185/)      !weights for gm uncert branches
       clamp = (/3.,6.,3.,6.,6.,6.,3.,300./)
      pdgk= (/0.,0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.12,0.14,
     &         0.16,0.18,0.20,0.22,0.24,0.27,0.30,0.33,
     &         0.36,0.4,0.46,0.5,0.6,0.75,0.85,1.0,1.5,
     &         2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/)
       aperiod = (/ 0.0,0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,
     1  0.2,0.25,0.3,0.4,0.5,0.75,1.,1.5,2.,3.,4.,5.,-1.0/)
c camper=Campbell-Bozorgnia NGA spectral period set, from 11-06 final(?) report. 
c 0.0=pga here, equiv to 0.01 in their report. -1=pgv, -2=pgd set. -3  = CAV added Oct 2012
      camper =(/0.01,0.020,0.030,0.050,0.075,0.100,0.150,0.200,0.250,0.300,0.400,
     1 0.500,0.750,1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5,10.0,0.,-1.0,-2.,-3.0/)
        PerIMIdriss  = (/0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.,1.5,2.,3.,4.,5.,7.5,10./)
       peras08= (/ 0.0, -1.0, -2.0, 0.01, 0.02, 0.022, 0.025, 0.029,
     1              0.03, 0.032, 0.035, 0.036, 0.04, 0.042, 0.044, 
     2              0.045, 0.046, 0.048, 0.05, 0.055, 0.06, 0.065, 
     3              0.067, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 
     4              0.11, 0.12, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 
     5              0.18, 0.19, 0.2, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 
     6              0.3, 0.32, 0.34, 0.35, 0.36, 0.38, 0.4, 0.42, 0.44,
     7              0.45, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 0.667, 0.7,
     8              0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 
     9              1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 
     1              3.0, 3.2, 3.4, 3.5, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 
     1              5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 /)
c prd is the C&Y period set, 106 of 'em, jan 2009. Same, Oct 2007. PGA=0.0s in our code. 
      prd= (/0.0,0.020,0.022,0.025,0.029,0.030,0.032,0.035,0.036,0.040,0.042,0.044,0.045,0.046,
     10.048,0.050,0.055,0.060,0.065,0.067,0.070,0.075,0.080,0.085,0.090,0.095,0.100,0.110,0.120,
     10.130,0.133,0.140,0.150,0.160,0.170,0.180,0.190,0.200,0.220,0.240,0.250,0.260,0.280,0.290,
     10.300,0.320,0.340,0.350,0.360,0.380,0.400,0.420,0.440,0.450,0.460,0.480,0.500,0.550,0.600,
     10.650,0.667,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,
     11.700,1.800,1.900,2.000,2.200,2.400,2.500,2.600,2.800,3.000,3.200,3.400,3.500,3.600,3.800,
     14.000,4.200,4.400,4.600,4.800,5.000,5.500,6.000,6.500,7.000,7.500,8.000,8.500,9.000,9.500,
     110.0,-1.0/)
      wtsum = (/0.,0.,0.,0.,0.,0.,0.,0./)
       ptail = (/0.9772499,.8413447,0.5,.1586553,.022750139, 0.0013499/)
      pithy = (/'Using Median','Median+EpUnc','Median-EpUnc'/)
      g_name =(/'Cluster GRP1','Cluster GRP2','Cluster GRP3','Cluster GRP4','Cluster GRP5'/)
      att_typ =  (/ 'CampH CEUS A','Space Unoccu','Space Unoccu',
     1 'Somervill HR','FeaHard CEUS','AB-hard CEUS','A&BHR CEUS06',
     1 'Somervil IMW','TORO-CEUS HR','Space Unoccu','Truth -  N/A',
     1 'Spudich 2000','TORO-CEUS BC','Sadigh ea 97',
     1 'A&B06BC CEUS','AB95-CEUS BC','Frankel eaBC','Somervill BC',
     1 'Abr-Silva 97','Campbell2003',
     2 'Camp CEUS BC','BJF     1997','Mota ea PRVI',
     2 'B&A 03/08NGA','C&B03/08 NGA',
     2 'ChiouY 03/08','Ab-Sil NGA08','Idriss   PGA','Kanno ea2006',
     2 'T&P CEUS  06','Silva CEUS02','AB06s200CEUS',
     3 'Space Unoccu','Space Unoccu','Space Unoccu',
     3 'AB08p Prelim','AB06p Prelim','Pez11 Prelim','Silva M-depS',
     3 'Space Unoccu','Space Unoccu','Space Unoccu',
     3 'Space Unoccu','BSSA 03/2013','C&B02/13 NGA',
     3 'ChiouY 03/13','Ab-Silva2013','Idriss 2012a','GraizKalkn12',
     4 'GrazKalkan09'/)
       gname=(/'.g1','.g2','.g3','.g4','.g5'/)
      do ia=1,8
      do ip=1,npmx
      nga(ip,ia)=.false.
      wus02(ip,ia)=.false.
      ceus02(ip,ia)=.false.
       ceus11(ip,ia)=.false.
c      benioff02(ip,ia)=.false.
      enddo
      enddo      !ia
      noas97 = .true.
      mcut(1)=6.0
      mcut(2)=7.0
      dcut(1)=10.
      dcut(2)=30.
      rate_cl=0.      !initialize clustered event rate to 0
      nft_out=0      !keeps tract of number of faults outside a specified perimeter
      gnd_ep=0.      !initialize epistemic variation of gnd to 0.
      poly = .false.
c Initialize truncated normal array Pex, store in p().
        prl=0.5*(erf(emax/sqrt2) + 1.0)
        prlr=1.0/prl
        prr=1.0-prl
        plim=-emax/sqrt2
      v30a = .false.
        pr0=plim
        dp=0.000004*(3.3-plim)
        dp2=1.0/dp
        p(1)=1.e-9
        do i=2,250001
        pr0=pr0+dp
        p(i)=((erf(pr0)+1.)*0.5-prr)*prlr
        enddo
        ymax=-90.0
        ymin=90.
        p(250002)=1.0
        grp_rate=0.0      !may accumulate wt*rate for clustered event groups
c End initializing the truncated normal PEx() array=p()
c The indep. variable is a real*8 to reduce discretization error.
      call date_and_time(date,time,zone,ival)
c      instructions from command-line file, argument 1
 900  format(a)
            override_vs=.false.
            callme=.true.
      call getarg(0,adum)      !keep tract of program name
      headr%name(3)= adum      !new 11/05. prog name (needs to be shortpath)
      headr%name(4)=date//' at '//time
c user name might help. Who to blame if all goes to H in a handbasket?
      call getlog(adum)
      headr%name(6)=adum
c In 2007 file names are sometimes > 30 chars long.
      call getarg(1,name)
      headr%name(1)= name(1:30)
      dum=dtime(tarray)
c      write(*,*)tarray
      inquire(file=name,exist=isok)
      if(isok)then
      open(unit=1,file=name,status='old')
      else
      write(6,*)'hazFXnga13l: File not found: ',name
      stop 'Put in working dir and retry'
      endif
      numarg=iargc()
c iargc() has been reported to be unstable with gfortran and Linux. By explicitly
c declaring iargc and numarg as integers, current code attempts a repair.
      if(numarg.ge.2.and.numarg.le.3)then
            call getarg(2,polygon)
            if(polygon.eq.'lowsigma')then
            lceus_sigma=.true.
            ceus_sigma=0.4
c lceus_sigma lower all CEUS GMPE sigmas to 0.4 (natural log sigma) for sens. study Mar 2011
            perris=.false.
            poly=.false.
            elseif(polygon.eq.'fixsigma'.or.polygon.eq.'FIXSIGMA')then
            fix_sigma=.true.
            call getarg(3,adum)
            read(adum,'(f6.3)')sigma_fx
            print *,'The code will use a fixed aleatory sigma of ',sigma_fx,' all periods all NGA-W gmpes'
            perris=.false.
            poly=.false.
            elseif(polygon.eq.'dy'.or.polygon.eq.'dY'.or.polygon.eq.'DY')then
            poly=.false.; perris=.false.
            sdi=.true.
            call getarg(3,adum)
            read(adum,'(f6.3)')dy_sdi
            print *,'The code will compute inelastic displ spectra with dy(cm) ',dy_sdi
            else
             inquire(file=name,exist=isok)
            if(isok)then
              open(unit=2,file=polygon,status='old')
              do i=1,80
              read(2,*,end=53)px(i),py(i)
53            enddo
            npmax=i-1
             else
             write(6,*)'hazFXnga13l: Polygon file not found: ',polygon
       	     stop 'Put in working dir and retry'
            endif
            poly=.true.
            endif      
      elseif(numarg.gt.3)then
c deaggregation information
            deagg=.true.
            safix=0.1
            haz=0.
            rbar=0.
            ebar=0.
            mbar=0.
            fmbar=0.; frbar=0.; febar=0.; fhaz=0.
            prob5=0.
            call getarg(2,adum)
            read(adum,'(f7.4)')rlatd
            call getarg(3,adum)
            read(adum,'(f9.4)')rlond
c rlatd rlond are the coordinates of the deagg-analysis site. This location overrides the stuff
c in the input file.
            call getarg(4,adum)
            read(adum,'(i1)')npd
            do i=5,npd +4
            call getarg(i,adum)
c adum could be sa(g) or pgv (cm/s). need flexible format
            read(adum,'(f8.4)')safix(i-4)
            write(6,*)'#command line sa ',safix(i-4)
            enddo
            if(numarg.gt.npd+4)then
            call getarg(npd+5,adum)
            read(adum,'(f6.1)')vs30d
            override_vs=.true.
            write(6,*)' vs30 will be ',vs30d,' m/s for this run. Command line'
            else
            override_vs=.false.
            endif
      endif
      write (6,61)date,time,name
61      format('# *** hazFXnga13l 03/04/2013 log file. Pgm run on ',a,' at ',a,/,
     +'# *** Input control file: ',a)
      if(poly)write(6,*)'hazFXnga13l: Polygon file &npts: ',polygon,npmax
c Below bypasses are based on file name. Bypass wont work if file names change
      byeext=index(name,'brange').gt.0
      byenv = index(name,'nv.').gt.0
c Calif fault file names 6/2007: aFault bFault
      byeca=index(name,'ca-a').gt.0.or.index(name,'aFault_').gt.0
     1  .or.index(name,'ca-wg').gt.0.or.index(name,'bFault_').gt.0
      byeoreg=index(name,'orwa_').gt.0.or.index(name,'orwa.').gt.0
      cal_fl = index(name,'unsegA').gt.0 .or. index(name,'aFault_unseg').gt.0
     + .or. index(name,'creepflt').gt.0
     + .or. (index(name,'aF01').gt.0..and.index(name,'_unseg').gt.0)
      bfault = index(name,'bFault_').gt.0
     + .or.index(name,'bF01').gt.0
c above name is used to determine if CAL floater: need a better way.
c creeping section is also host to california floaters. Parkfield did not rupture to
c surface. After short discussion with Wesson, I did not include creeping sec
c set of calif. floaters that always rupture to surface. 
c      write(6,*)'Enter a zero for grid of sites 1 to 30 for list: '
      read(1,*)nrec
      if(nrec.lt.0)then
      determ=.true.
      deagg=.false.      ! determ trumps deagg
      if(nrec.lt.-35)then
      nrec=0
      else
      nrec=-nrec
      endif      !gridded or specified stations option
      else
      determ=.false.
      endif
      if(nrec.eq.0)then
      grid=.true.
      slist=.false.
      write(6,580) "For sites: enter min lat, max lat, dlat: "
      read(1,*) ymin, ymax, dy
      write(6,*)ymin,ymax,dy
      write(6,580) "for sites: enter min lon, max lon, dlon: "
      read(1,*) xmin, xmax, dx
      write(6,*)xmin,xmax,dx
      nx= nint((xmax-xmin)/dx) +1
      ny= nint((ymax-ymin)/dy) +1      !nint will make this calc more robust,
c portable among different computers
      write(6,*) nx,ny,' nx ny for discrete grid of sites'
      nrec= nx*ny
      elseif(nrec.lt.33)then
      write(6,*)'*** Station-list option has been selected.'
      grid=.false.
      slist=.true.
      do i=1,nrec
      write(6,*)'Enter station lat,long (dec deg) and name(1 word): '
580      format(a,$)      
      read(1,*)slat(i),slong(i),sname(i)
      write(6,*)slat(i),slong(i),sname(i)
c Richard Dobbs noticed that ymax needs to be initialized 11/08.
      ymax= max(ymax,slat(i))
      ymin= min(ymin,slat(i))
      if(.not.deagg)open(30+i,file='eqrate.'//sname(i),status='unknown')
      if(poly)then
      write(30+i,3890)headr%name(1),polygon
      elseif(.not.deagg)then
      write(30+i,389)headr%name(1)
      endif
 389      format('#Pgm hazFXnga13l.f file of rates*weights for fault ruptures',/,
     + '#Dist(km)  Mw   Rate*wt, ift,Ztor(km),idis Infile: ',a)     
 3890      format('#Pgm hazFXnga13l.f file of rates*weights for fault ruptures',/,
     + '#Dist(km)  Mw   Rate*wt, ift,Ztor(km),idis Infile: ',a,' polygon ',a)     
c write eqrates to these files New, nov 7 2006.      
      if(slat(i).lt.-88.)stop 'invalid station location. Version 7c?'
      enddo
      else
      write(6,*)'hazFXnga13l expects first line of input to be nrec'
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
      slist=.true.
      endif
c **** Enter soil Vs30 condition  ******NEW*******
c      write(6,*)"For sites, enter Vs30(m/s) and basin depth(km)"
      read(1,*)vs30,dbasin
      dgkbasin=max(1.0,dbasin*0.5)      !Vladimir says his basin is 1.5 km/s isosurface. Campbell is
c       				2.5 km/s isosurface.
      if(override_vs)vs30=vs30d
c the below Z1 is defined so we do not need to rewrite input files with a Z1 value.
c Units m. Z1cal was modified to equal the CY report eqn 2. Mar 11 2013.
        Z1cal = exp(-7.15/4 * log((VS30**4 + 570.94**4)/(1360.**4 + 570.94**4)))
c     Norm Abrahamson's CA z1 reference (eq 18)
       z1_ref = exp ( -7.67/4. * alog( (Vs30**4 + 610.**4)/(1360.**4+610.**4) ) ) / 1000.
      z1=z1cal	!CY2013 function used until we know better for wus...
      z1km=Z1*0.001	!for AS need units km
c this Z1 is  the recommended value of CY. May need to be changed (40 m for 760 m/s vs30)
      write(6,*)' Vs30 (m/s), Z1 (m) and depth of basin (km): ',vs30,Z1,dbasin
      if(vs30.lt.90..and.vs30.gt.0.)then
      write(6,*)'Vs30 = ',vs30,'. This looks unreasonable.'
      stop'please check input file just before distance incr,dmax'
      elseif(vs30.eq.0)then
      write(6,580)'Enter the name of the binary vs30 array: '
      read (1,900)name
      stop'Please use code hazFXnga7vs.f for variable Vs30 case.'
      inquire(file=name,exist=isok)
      if(.not.isok)stop'Vs30 array not in working directory'
      call openr(name)
      write(6,580)'Enter ymin, ymax, dy (lat, d) of Vs30 array: '
      read (1,*)vymin,vymax,vdy
      write(6,580)'Enter xmin, xmax, dx (long, d) of Vs30 array: '
      read (1,*)vxmin,vxmax,vdx
      write(6,580)'Enter a default vs30 in case site not inbounds: '
      read (1,*)vs30dflt
      if(vymin.gt.vymax)then
      y0=vymin
      vymin=vymax
      vymax=y0
c obvious repair if user has transposed. common input error.
      endif
      nvx= nint((vxmax-vxmin)/vdx) +1
      nvy= nint((vymax-vymin)/vdy) +1      
c Using the function nint will make this calc more robust,
c portable among different computers
      v30a=.true.
      nv30=nvx*nvy
c 30 character name
c      nhd=308
c 128 character name
      nhd=896
      call gethead(hdvs30,nhd,nrd)
      if(abs(hdvs30%extra(2)-vxmin).gt.0.001.or.abs(hdvs30%extra(3)-vxmax).gt.0.001)then
      print *,'Warning longitude info for  Vs30 does not match header. USe header'
      vxmin=hdvs30%extra(2); vxmax=hdvs30%extra(3)
      endif
      if(abs(hdvs30%extra(5)-vymin).gt.0.001.or.abs(hdvs30%extra(6)-vymax).gt.0.001)then
      print *,'Warning latitude info for  Vs30 does not match header'
      vymin=hdvs30%extra(5); vymax=hdvs30%extra(6)
      endif 
      if(abs(hdvs30%extra(4)-vdx).gt.0.00001.or. abs(hdvs30%extra(7)-vdy).gt.0.00001)then
      print *,'Warning, sampling interval vdx or vdy not equal to header.Using header'
      vdx=  hdvs30%extra(4); vdy= hdvs30%extra(7);
      endif
      nvx= nint((vxmax-vxmin)/vdx) +1
      nvy= nint((vymax-vymin)/vdy) +1      
c revise the data to read to agree with header,
c portable among different computers
      nv30=nvx*nvy
      call getbuf2(v30,nv30,nvread)
      write(6,*)'For vs30 expect ',nv30,' got ',nvread
      write(6,*)'If using CY2007 model, need to redefine Z1 to correspond to Vs30'
      if(vxmax.le.-100.)print *,'Please be sure all vs30 are less than 1500 m/s'
c something like v30=min(v30,1500.) works for the west. with NGA relations. SH Mar 6 2008. We
c would like to let user know that code is changing his /her model if redefining v30 at some nodes. 
      endif
      srcSiteA=90.      !this angle should vary with site location. See delaz calls
c      write(6,*) "enter distance increment, dmax"
      read(1,*) di,dmax
      write(6,*)'Distance increment and dmax (km): ',di,dmax
cccccccccccccccccccccccccccccccccc
c      write(6,*) "enter number of periods"
      read(1,*) nper
c hazFXnga13l: this version  IS NOT defining the 22 degree wings of AS97.
c      angtest= 22.5*3.14159/180.
c      angtest= tan(angtest)
c---loop through periods
      do 700 ip=1,nper
      read(1,*) period(ip), wind
      per= period(ip)
      if(sdi)fac_sde(ip) = alog(980.*max(0.01,per)**2/fourpisq)
      k=1
      dowhile(per.ne.perx(k))
      k=k+1
      if(k.gt.8)then
      write(6,*) '*** Input period not in 2002 group ***', per
      write(6,*)'Warning: some atten relations may not be ready.'
      write(6,*)'Warning 2: clamp is only prepared for 2002 group.'
      k=6
      goto 801      !not really a solution
      endif
      enddo
801      write(6,*)'Period ',per,' underway'
      iper(ip)=k
      if(wind.gt.0.)then
      cluster=.false.
      l_gnd_ep(ip)=.true.
      do im=1,3
      if(im.eq.1)then
      write(6,505)mcut(1)
      elseif(im.eq.2)then
      write(6,506)mcut(1),mcut(2)
      else
      write(6,507)mcut(2)
      endif
      write(6,509)dcut(1),dcut(1),dcut(2),dcut(2)
509      format('3 DeltaGnd, for d < ',f4.1,', ',f4.1,' <=d < ',f4.1,
     + '& d >=',f4.1,' km: ',$) 
      read(1,*)gnd_ep(1,im,ip),gnd_ep(2,im,ip),gnd_ep(3,im,ip)
c increase or decrease gm by equal amounts in the atten subroutines
      write(6,*)gnd_ep(1,im,ip),gnd_ep(2,im,ip),gnd_ep(3,im,ip)
505      format('Additional epistemic gnd, for M < ',f4.1)
506      format('Next, for ',f4.1,'<=M < ',f4.1)
507      format('Finally, for M >= ',f4.1)
      enddo      !im loop
      nfi(ip)=3
      ngroup(ip)=1
      icc(ip,1)=9+ip
      elseif(wind.eq.0.)then
      nfi(ip)=1
      ifn=1      !always use index 1 for this case.
      ngroup(ip)=1
      cluster=.false.
      icc(ip,1)=9+ip
      else
c current indicator to perform clustering is wind < 0. You can make up to 5 independent clustered-event curves
c set wind=-1 for one, wind=-2 for 2, ..., wind = -5 for 5. These could be the 5 virtual NMSZ faults in 2002 hz model
c The use of this variable for different things depending on its sign might be considered clumsy & contemptible.
c Should be revised.
      nfi(ip)=1
      ifn=1
      cluster=.true.
      ngroup(ip)= iabs(int(wind))
      if(ip.gt.1.and.ngroup(ip).ne.ngroup(1))stop'number of cluster groups cannot be period dep.'
      endif      !if additional epistemic sigma is read in or model clustering is asked for: one or the other
c but not both at least initially
      nameout= '                                             '
c      write(6,*) "enter name of output file for this period"
      read(1,909) nameout
 909  format(a)
       isz=index(nameout,' ')-1
      if(determ)then
      open(66+ip,file=nameout(1:isz)//'.DET')
      write(66+ip,677)per,date,time,name
      endif
677      format('#hazFXnga13l sources. Sp_Per= ',f5.2,' s. Run on ',a,' at ',a,/,
     +'# *** Input control file: ',a,/,
     +'#MedianSA(g)     sd   iattn  iflt    M  Rjb(km) Rcd(km)  Wt*Rate',
     +' Dtor(km) SlipCode EPS_AT')
       do kg=1,ngroup(ip)
       if(cluster)then
       nameout=nameout(1:isz)//gname(kg)
       icc(ip,kg)=9+ip+9+kg
       endif
       if(grid)then
      call openwx(ifp(ip,kg,1),nameout)
      write(6,*) ifp(ip,kg,1),nameout
      if(nfi(ip).eq.3)then
      write(6,*)'Additional files for epistemic gm branches:'
c i have seen garbage appended after the .p and .m. This needs to be corrected.
      adum='                                               '
      adum=nameout(1:isz)//'.p'
      call openwx(ifp(ip,1,2),adum)
      write(6,*) ifp(ip,1,2),adum
      adum=nameout(1:isz)//'.m'
      call openwx(ifp(ip,1,3),adum)
      write(6,*) ifp(ip,1,3),adum
      endif      !nfi = 3. Open  2 extra files for epistemic branch curves
      else
c if ascii put each curve in same output file. Slight repair May 1 2008. SH.
      if(kg.eq.1)then
      nameout=nameout(1:isz)//'.asc'
      open(icc(ip,1),file=nameout,status='unknown')      !ascii
       write(6,895)nameout
       write(icc(ip,1),402)date,headr%name(1)
402      format('#hazFXnga13l(1/29/2009) run on ',a9,' with file ',a)
895      format('Ascii output file: ',a)
      endif
      if(deagg.and.kg.eq.1)then
c cannot have both deaggregation and determinstic output. Safe to use same unit number.
c First try: deagg. epistemic gm uncert or CEUS geographic flt loc uncert will get combined into
c same rbar, mbar, ebar, haz array. These could be separated if nec (add a "kg" dimension).
      if(nper.ne.npd)stop'Deaggregation command-line np does not equal input nper'
       isz=index(nameout,' ')-1
c open ascii deagg output. Index arrays by distance, magnitude, epsilon, gmuncert, period.
      open(ip+20,file=nameout(1:isz)//'.DEAG',status='unknown')
      write(20+ip,907)safix(ip),period(ip),date,time,name,rlatd,rlond,vs30
c Separate file for individual fault hazard. Does SAF dominate S Hayward at site S, etc?
      open(40+ip,file=nameout(1:isz)//'.FDEAG',status='unknown')
      write(40+ip,929)safix(ip),period(ip),date,time,name,rlatd,rlond,vs30
 907      format('#hazFXnga13l deagg @SA=',f5.3,' g. T='f5.2,
     +' s, run on ',a,' at ',a,' fi ',a,
     +/,'#sta lat long = ',f7.4,1x,f9.4,' Vs30 (m/s) is ',f6.1)      
 929      format('#hazFXnga13l deagg @SA=',f5.3,' g. T='f5.2,
     +' s, run on ',a,' at ',a,' fi ',a,
     +/,'#sta lat long = ',f7.4,1x,f9.4,' Vs30 (m/s) is ',f6.1,/,
     +'#ift Rbar(km)  Mbar(Mw) E0bar Wt*Rate_Exc gm_uncert_indx')      
      endif      !if deagg. new oct 30
      endif      !if grid
      enddo      !kg=1,ngroup      added may 16 2007
c      write(6,*) "enter number of ground motion levels"
      read(1,*) nlev(ip)
      if(nlev(ip).gt.nlvmx)stop'number of gm levels exceeds max'
c      write(6,*) "enter SA,PGA ground motion levels"
      read(1,*) (xlev(k,ip),k=1,nlev(ip))
         if(deagg)then
         nlev(ip)=1
         xlev(1,ip)=safix(ip)
         elseif(sdi)then
c convert sa to sd
        do k=1,nlev(ip)
        xlev(k,ip)=exp(fac_sde(ip))*xlev(k,ip)
        enddo
         endif
      write(6,*)'Nlev ',nlev(ip),' min max ',
     1 xlev(1,ip),xlev(nlev(ip),ip)
         if(.not.deagg)write(6,*) (xlev(k,ip),k=1,nlev(ip))

      do 414 k=1,nlev(ip)
      ylev(k,ip)=xlev(k,ip)
      if(k.gt.1.and.xlev(k,ip).lt.1.05*xlev(k-1,ip))then
      write(6,*)'*** The following ground motion progression is not good:'
      write(6,*)xlev(k-1,ip),xlev(k,ip)
      stop'*** Please be sure ground motion increases at next sample ***'
      endif
 414  xlev(k,ip)= alog(xlev(k,ip))
c-----------
c      write(6,*) "enter number of atten. relations for this period"
      read(1,*) nattn(ip)
      write(6,*)"number of atten. relations for this period ",nattn(ip)
c--- loop through atten relations for that period
      do 701 ia=1,nattn(ip)
c      write(6,*) "enter type of atten. relation, weight1, wtdist, 
c     &  weight2 , mb to M conv."
      read(1,*) iatten(ip,ia),wt(ip,ia,1),wtdist(ip,ia),wt(ip,ia,2),
     &   icode(ip,ia)
      ipia=iatten(ip,ia)
       if(ipia.eq.8)noas97=.false.
c      Certain calcs can be bypassed if not using Abrahamson-Silva 97 attn model
c hazFXnga13l: there is no longer a precalc of a pr() array. Note how 
c for valid iatten() numbers the calls are commentd out.
c add Somerville  imw to WUS02 set (even though it is more recent than 2002.
c This relation will need more periods and site response if  it is to be 
c routinely used. Added for special studies Jan 19 2007. SHarmsen.
       wus02(ip,ia)=ipia.eq.1 .or.ipia.eq.3.or.
     1 ipia.eq.8.or.ipia.eq.9.or.ipia.eq.11.or.ipia.eq.-3
       ceus02(ip,ia)=abs(ipia).eq.2 .or.ipia.eq.19.or.ipia.eq.20.or.
     1 abs(ipia).eq.4.or.abs(ipia).eq.5.or.ipia.eq.21.or.ipia.eq.28.or.
     1 abs(ipia).eq.6.or.abs(ipia).eq.7.or.abs(ipia).eq.10
       ceus11(ip,ia)=ipia.gt.24.and.ipia.lt.28  !new mar 2011.
      nga(ip,ia)=(ipia.gt.12.and.ipia.lt.19).or.ipia.gt.30
c kanno et. al. is included with NGA even though it's not. But is modern.
c prepare look-up tables for certain CEUS relations.
        if(ceus11(ip,ia))then
        kf=1
        if(per.gt.0.01)then
        fr=1.0/per
        elseif(per.ge.0.0)then
        fr=99.
        elseif(per.eq.-1)then
        fr=89.
        else
        print *,'period requested ',per
        stop' CEUS11 models close encounter with an unknown period'
        endif
        dowhile(fr.gt.a11fr(kf)+.01)
        kf=kf+1
        if(kf.gt.13)stop' period not in A06,A08,P11 set'
        enddo
        freq(ip)=a11fr(kf)
        if(ipia.eq.25)then
        if(callme(1))then
        print *,'calling GailTable A08revA_Rjb fr,per,kf=',fr,per,kf
        name ='GR/A08revA_Rjb.dat'
        name=trim(name)
             open(3,file=name,status='old',err=202)
        call GailTable(1)
        callme(1)=.false.
        endif 
        ia08(ip)=kf  
        elseif(ipia.eq.26)then
         if(callme(2))then
       name ='GR/AB06revA_Rcd.dat'
        print *,'calling GailTable AB06revA_Rcd fr,per=',fr,per,kf
             open(3,file=name,status='old',err=202)
       call GailTable(2)
       callme(2)=.false.
       endif
        ia06(ip)=kf 
        elseif(ipia.eq.27)then
        if(per.lt.0.0)stop' Pezeshk 2011 does not have PGV coeffs'
        if(callme(3))then
        name= 'GR/P11A_Rcd.dat'
             open(3,file=name,status='old',err=202)
       call GailTable(3) 
        print *," just read Gail's Pezeshk 2011 A-table",fr,per,kf
       callme(3)=.false.
       endif   !read table 3?
        ip11(ip)=kf
       endif       !if ipia = 25, 26, or 27
       elseif(ipia.eq.6) then
       call getFEAtab(iper(ip),1)
       write(6,*)'FEA BC table OK for period ',perx(iper(ip))
       elseif(ipia.eq.5) then
       call getFEAtab(iper(ip),2)       
       write(6,*)'ABceus BC table OK for period ',perx(iper(ip))
       elseif(ipia.eq.-6)then
       call getFEAtab(iper(ip),3)
       write(6,*)'FEA hardrock table OK for period ',perx(iper(ip))
      if(vs30.lt.800.)write(6,*)'Warning: vs30 ',vs30,' not consistent with HR'
       elseif(ipia.eq.-5)then
       call getFEAtab(iper(ip),4)
       write(6,*)'ABceus hardrock table OK for period ',perx(iper(ip))       
      if(vs30.lt.800.)write(6,*)'Warning: vs30 ',vs30,' not consistent with HR'
       endif
             wtsum(ip)=wtsum(ip)+wt(ip,ia,1)
             if(wtsum(ip).gt.1.01)then
             write(6,*)' For period ',per,' sum of atten-model wts = ',wtsum(ip)
             stop'This is unreasonable. Please check input file.'
             endif
c
      write(6,807)att_typ(ipia),ia
        if(abs(ipia).eq.4.or.ipia.eq.21)then
          ka=1
          dowhile(per.ne.abper(ka))
            ka=ka+1
            if(ka.gt.26)then
              write(6,*) 'As of 6/06 input period doesnt correspond to A&B06 set'
              stop 'Please remove this relation from input file'
            endif
          enddo
          iperab(ip)=ka     
c You can call the hardrock table when Vs30 >= 1500 m/s. But not clear that AB would approve
          if(vs30.lt.1500..and.ipia.lt.0)write(6,808)vs30
          if(vs30.ge.2000.and.ipia.eq.21)then
          irab(ip,ia)=4
          iatten(ip,ia)=4      !redefine
          write(6,*)'A (no S) rock & 200-bar CEUS stress parameter assumed for AB06'
          elseif(vs30.lt.2000..and.ipia.eq.21)then
          irab(ip,ia)=3
          iatten(ip,ia)=4      !redefine
          write(6,*)'Site-corrected (S) rock & 200bar CEUS stress parameter assumed for AB06'
          elseif(vs30.ge.1500..and.ipia.eq.-21)then
          irab(ip,ia)=4
          iatten(ip,ia)=4
          write(6,*)'A (no S) rock & 200-bar CEUS stress parameter assumed for AB06'          
          elseif(ipia.eq.-4)then
          write(6,*)'A (no S) rock & 140-bar CEUS stress parameter assumed for AB06'
          iatten(ip,ia)=4
          irab(ip,ia)=2
          else
          write(6,*)'Site-corrected (S) rock & 140-bar CEUS stress parameter assumed for AB06'
          irab(ip,ia)=1
          endif
          write(6,*)ip,ka,' A&B 9/06 ip map; irab(ip,ia)= ',irab(ip,ia)
  808      format('You have called hardrock version of AB06 even though vs30 is ',f6.1)   
        endif
        if(ipia.eq.19)then
          ka=1
          dowhile(per.ne.tpper(ka))
            ka=ka+1
            if(ka.gt.16)then
              write(6,*) 'As of 7/2008 input period doesnt correspond to T&P05 set'
              stop 'Please remove this relation from input file'
            endif
          enddo
          ipertp(ip)=ka
          write(6,*)ip,ka,' T&P 7/08 ip map'
      if(vs30.ge.1500.)then
      irtb=-1
      write(6,*)'TP05 relation is called with hardrock coeff'
      else
      irtb=1
      endif
        elseif(ipia.eq.20)then
          ka=1
          dowhile(per.ne.pdSilva(ka))
            ka=ka+1
            if(ka.gt.11)then
              write(6,*) 'As of 7/2008 input period doesnt correspond to Silva2002 set'
              stop 'Suggestion: remove this relation from input file for this period'
            endif
          enddo
          isilva(ip)=ka
          write(6,*)ip,ka,' Silva 7/2008 ip map'
          if(vs30.ge.1500.)then
          irsilva=-1
      write(6,*)'Silva relation is called with hardrock coeff'
          else
          irsilva=1
          endif
         elseif(ipia.eq.28)then
           ka=1
           dowhile(per.ne.pdSilva2(ka))
             ka=ka+1
             if(ka.gt.8)then
               write(6,*) 'As of 7/2008 input period doesnt correspond to Silva Variable stress set'
               stop 'Suggestion: remove this relation from input file for this period'
             endif
           enddo
           isilva2(ip)=ka
           write(6,*)ip,ka,' SilvaV 7/2008 ip map'
           if(vs30.ge.1500.)then
           irsilva=-1
       write(6,*)'Silva M-dep. stress relation is called with hardrock coeff'
           else
           irsilva=1
           endif
        endif
        if(ipia.eq.-3)then
c Somerville IMW, 5 periods
      ka=1
      dowhile(per.ne.prdsom(ka))
      ka=ka+1
      if(ka.gt.5)stop 'input period doesnt correspond to Somerville IMW set'
      enddo
       isomer(ip)=ka
      if(vs30.lt.600..or.vs30.gt.900.)then
      write(6,*)'********** WARNING Somerville IMW called with non-BC rock *******'
      write(6,*)'** Please consider modifying Somerville for this site condition**'
      endif
c For NGA relations, is the requested spectral period available?
      elseif(ipia.eq.16)then
c Below, find index from abrahamson 2008 period set
      ka=1
      dowhile(per.ne.peras08(ka))
      ka=ka+1
      if(ka.gt.108)stop 'input period doesnt correspond to abram.silva08 set'
      enddo
       ipera(ip)=ka
       vs30_class=1      !Vs30 measured not estimated (SR recommends email 3/15/13)
       elseif( ipia.eq.14 )then   
c Below, campbell-bozorgnia 11-07  period set
      ka=1
      dowhile(per.ne.camper(ka))
      ka=ka+1
      if(ka.gt.25)stop 'Input period doesnt correspond to Campbell-Bozorgnia 11-07 set'
      enddo
       ipercb(ip)=ka  
       write(6,*)ip,ipercb(ip),' C&B 11/07 ip map'
        if(ka.eq.25)print *,' This run computes probabilistic CAV'
c below added mar 6 2008 from email comment by Ken Campbell
       if(vs30.gt.1500.)stop 'Vs30 >1500 m/s and CB NGA relation does not permit this.'   
cc Below, boore/atkinson period set. 4/07 version of BA model has 23 periods. 
c  3 s SA corresponds to index 7.
      elseif(ipia.eq.13)then
      ka=1
      dowhile(per.ne.perb(ka))
      ka=ka+1
c BA revision of march 2007: 23 spectral periods incl pga and pgv.
      if(ka.gt.23)then
      write(6,*)'Period ',per
      stop 'For this period please remove the BA relation from input file'
      endif
      enddo
       iperb(ip)=ka     
       write(6,*)ip,iperb(ip),' BA 4/2007 ip map'   
c Below, Chiou/Youngs 11/2007 period set
      elseif(ipia.eq.15)then
      ka=1
      dowhile(per.ne.prd(ka))
      ka=ka+1
      if(ka.gt.106)then
      write(6,*) 'As of 01/2009 input period doesnt correspond to chiou&y.set'
      stop 'Please remove this relation from input file'
      endif
      enddo
       ipercy(ip)=ka     
       write(6,*)ip,ka,' CY-NGA 03/2008 ip map'   
        call CY2007I(ip,ka, vs30, Z1)
c initialize some site-related terms for the new CY relation for each spectral period.
      elseif(ipia.eq.35)then
      deltaZ1=0.0	!dont know use 0. from guidance in CY doc.
      if(per.eq.-1.)then
       icy13(ip)=1
      write(6,*)'Calling CY2013 NGA-W with period index 1: PGV'
      elseif(per.eq.-2.)then
       icy13(ip)=2
      write(6,*)'Calling CY2013 NGA-W with period index 2: PGD'
      elseif (per.ge.0..and. per.le.0.01)then
      icy13(ip)=3	!PGA index is 3 in Mar 2013 update
      write(6,*)'Calling CY2013 NGA-W with period index 3: PGA'
      else
      k=4
      dowhile(percy13(k).ne.per.and.k.lt.24)
      k=k+1
      enddo
      if(abs(percy13(k)-per).gt.0.002)stop'period not available 
     * for CY2013 GMPE'
      write(6,*)'Calling CY2013 NGA-W revision with deltaz1=',deltaz1
      write(6,*)'Calling CY2013 NGA-W with period index ',k
      icy13(ip) = k
      endif
c!pga or other?
      elseif(ipia.eq.36)then
      if (per.le.0.01)then
      ipas13(ip)=1
      else
      k=2
      dowhile(per.ne.Peras13(k).and.k.lt.23)
      k=k+1
      enddo
      if(per.ne.Peras13(k))then
      print *,' period ',per,' not available for AS13(V10) model'
      stop 'please remove this GMPE from your input file'
      endif
      ipas13(ip)=k
      	endif
      write(6,*)'Calling AS2013 NGA-W revision with z1=',z1km,' period index ',ipas13(ip)
c 12/17L commandeer icode to use as a decision on Ry0 usage. 1 = use it.
      if(icode(ip,ia).eq.1)then
      useRy0(ip)=.true.
      else
      useRy0(ip)=.false.
      endif
      print *,'Using Ry0 metric for this spectral period?',useRy0(ip)
        vs30_rock = 1100.
        z10_rock = 0.006
c above coefficients are needed by the AS2013 model      
      elseif(ipia.eq.34)then
      if(readcb13)then
      if(per.gt.0.01)stop'CB13 must be called with PGA as first g.m.'
      call cb13_nga_spec_in
      readcb13=.false.
      k=22
      else
      k=2
      dowhile (Percb13(k).ne.per.and.k.lt.22)
      k=k+1
      enddo
      if(k.eq.22)stop' Period not found for CB13 relation.'
      endif	!pga or other spectral accel?
      ipcb12(ip)=k
      print *,'CB13 relation period index ',k,' for period ',per
      elseif(ipia.eq.37)then
      if(per.le.0.01)then
      k=1
      else
      k=2
      dowhile(PerIMIdriss(k).ne.per.and.k.lt.22)
      k=k+1
      enddo
      if(k.eq.22.and.per.ne.10.)stop' Period not found for IMIdriss relation.'
      endif	!pga or other spectral accel?
      idriss(ip)=k
      print *,'Idriss relation period index ',k,' for period ',per
      
      elseif(ipia.eq.33)then
c 2013 bssa just include the coeffs in the subr. do not read in. too much detail
c
      indx_pga=2;indx_pgv=1
      if(per.eq.-1.)then
      k=1
      print *,'BSSA 2013 relation called for PGV '
      elseif(per.eq.0.0)then
      k=2
      print *,'BSSA 2013 relation called for PGA '
      print *,' BSSA index for pga ',indx_pga
      else
      k=3
      dowhile(Perbssa13(k).ne.per.and.k.lt.23)
      k=k+1
      enddo
      if(k.eq.23.and.per.ne.10.)stop' Period not found for BSSA relation.'
c      if(fix_sigma)sigt_gmpe=sigma_fx      !override table with fixed sigma jan 7 2012.
       nper_gmpe = 23
c      print *,nper_gmpe,' number of periods having coeffs BSSA'
      print *,per,Perbssa13(k),' BSSA period match? Index is ',k
      endif
      indbssa(ip)=k
      elseif(ipia.eq.39.or.ipia.eq.38)then
      if (per.lt.0.01)then
      ipgk(ip)=0
      else
      k=1
      dowhile(pdgk(k).ne.per)
      k=k+1
      if(k.gt.35)stop' Graizer Kalkan model does not have this pd'
      enddo
      ipgk(ip)=k
      endif
      if(ipia.eq.39)then
      write(6,*)'Calling Grazier Kalkan09 model. 0 for basin. Pd index ',k
      else
      if(per.le.0.01)ipgk(ip)=1	!use the sa1 for new graizer model
      write(6,*)'Calling Grazier Kalkan12 model.  Pd index ',ipgk(ip)
      write(6,*)'Graizer12 basin depth (km) set to ',dgkbasin
      Q_CA=435
      write(6,*)'Graizer12 Quality factor = ',Q_CA
      endif
        elseif(ipia.eq.18)then
      ka=0
      dowhile(per.ne.perka(ka))
      ka=ka+1
      if(ka.gt.37)then
      write(6,*)'Period ',per
      stop 'For this period please remove the Kanno relation from input file'
      endif
      enddo
       iperk(ip)=ka     
       write(6,*)ip,iperk(ip),' Kanno 10/06 ip map '   
c vsfac used for site-amp in  the Kanno et al. relation.
        vsfac=alog10(vs30)          
c Below,Dahle 1995 period set
      elseif(ipia.eq.22)then
      ka=1
      dowhile(per.ne.perdahl(ka))
      ka=ka+1
      if(ka.gt.6)then
      write(6,*) 'As of 08/2009 input period doesnt correspond to Dahle.set'
      stop 'Please remove this relation from input file'
      endif
      enddo
       iperdahl(ip)=ka     
       write(6,*)ip,ka,' Dahle ip map'
       isoild=0
       if(vs30.lt.550.)isoild=1   
      endif      !spectral period indexes for NGA relations
      if(ia.eq.1)headr%name(2)=att_typ(ipia)
807      format(a12,1x,'attenuation model assoc. with index ',i1)
 701  continue
 700  continue
ccccccccccccccccccccccccccccccccccc
c      write(6,*) "enter increment dlen and dmove for floating rup"
      read(1,*) dlen, dmove
      write(6,*)'dlen, dmove (km)=',dlen,dmove
c      write(6,*) "enter number of branches for Mchar logic tree"
      read(1,*) nbranch
      write(6,*)'Number of epistemic M-branches: ',nbranch
      if(cluster.and.nbranch.gt.1)then
      write(*,*)'Clustering is not currently compatible with Mchar branching.'
      write(*,*)'Clustering requires fixed event rates, branching induces variable rate'
      stop'hazFXnga13l: Please use no mag variability when clustering eq events.'
      endif
c      write(6,*) "enter dM of each branch"
      read(1,*) (dmbranch(i),i=1,nbranch)
c      write(6,*) "enter weight of each branch"
      read(1,*) (wtbranch(i),i=1,nbranch)
      write(6,*)'epistemic dM= ',(dmbranch(i),i=1,nbranch)
c      write(6,*) "enter std dev and width for Mchar aleatory uncer"
      wts=0.
      do i=1,nbranch
      wts=wts+wtbranch(i)
      enddo
      if(abs(wts-1.).gt.0.001)stop 'epistemic Mbranch wts do not sum to 1.'
      read(1,*) sdal,mwid
      dma=0.4*sdal      !so that 5 dma  = 2 sd_aleatory
      write(6,*)sdal,' sd_aleatory, step is ',dma
      if(sdal.gt.0.) then
         denom=0.
         alet=0.0
c sum of moments to equal initial moment (slip? constraint)
         do 670 i= -mwid, mwid
         pwr=float(i)*dma
         ale(i+mwid+1)= exp(-pwr*pwr/(2.*sdal*sdal))
         denom= denom+ ale(i+mwid+1)*10.**(1.5*pwr)
c         write(6,*) ale(i+mwid+1)
 670     continue
         do i=1,2*mwid+1
       ale2(i)= ale(i)/denom
       alet=alet+ale2(i)
       enddo      
         write(6,*) 'aletot ',alet,' denom ',denom
c         read(5,*) idum
      elseif(sdal.lt.0.) then
c sum of rates to equal initial rate (other geologic constraint)
c added aug 13 2007. SH.
         denom=0.
         alet=0.0
         sdfac=2.*sdal**2
         do 680 i= -mwid, mwid
         pwr=float(i)*dma
         ale(i+mwid+1)= exp(-pwr*pwr/sdfac)
         denom= denom+ ale(i+mwid+1)
c         write(6,*) ale(i+mwid+1)
 680     continue
         do i=1,2*mwid+1
       ale2(i)= ale(i)/denom
       alet=alet+ale2(i)
       enddo      
         write(6,*) 'Rate-fixed aletot ',alet,' denom ',denom
c         read(5,*) idum
         endif            !sdal .ne. 0
ccc-------read fault data
      do 1 ift=1,nfltmx
      write(6,*) ift
c      write(6,*) "enter 1 for char. 2 for G-R; 1 for SS, 2 for reverse; nmagep"
      read(1,900,end=999) adum
      write (6,900)adum
c blank lines at end of file are ok. should not cause hiccough below.
c nmagf (ift) is a new required field 6-06. Number of mag branches. Used for
c several SF bay area faults (WG99), for NMSZ and others where the epistemic
c model needs some flexibility. 
c Added relwt to each fault-mag line. This relwt is relative to all faults
c in the current input file, not just this particular fault. For example, if
c combining 25-wt and 50-wt California faults into one file, and the relwt for the
c 25-wt items is X, then the relwt for the 50-wt items is 2X. This was discussed
c at a software meeting late Oct 2006. S Harmsen.
c add moment rate to log file output 11/2010
      if(cluster)then
      read(adum,*,err=999,end=999)itype(ift),iftype(ift),nmagf(ift),igroup(ift),jseg(ift)
c May 16 2007: igroup is a new index, for the geographic or other group index for clustering. 
c Idea is : Do not mix scenarios from different groups.
c All scenarios within a group must have same recurrence (rate).
c jseg is the segmentation model, should be 1, 2, or 3 for north, central, or south segment of NMSZ fault system
      else
c code should be able to work with earlier file format if clustering is not requested.
      read(adum,*,err=999,end=999)itype(ift),iftype(ift),nmagf(ift)
      endif
            if(iftype(ift).eq.1)then
      ibtype(ift)=1
      cDPP(ift)=0.5
      elseif(iftype(ift).eq.2)then
      ibtype(ift)=3	!boore 2012 swithces
      cDPP(ift)=0.25
      else
      ibtype(ift)=2
      cDPP(ift)=0.125	!temporary setting for amount of directivity effect
      endif

      if(itype(ift).gt.0)then
      nmg= nmagf(ift)
      lb0=.false.
      elseif(itype(ift).lt.-1)then
c -2 is a trick number, it means n branches are specified, and n others with bvalue=0 are derived here.
c mod of aug 28 2007 to study sensitivity of rate to b=0 (half relative wt)
        lb0= .true.
      nmg=2*nmagf(ift)
      print *,ift,' modified nmagf ',nmagf(ift),' to ',nmg
      else
      print *,adum
      stop 'invalid itype 0'
      endif
        if(lb0)itype(ift)=2
      if(nmg.lt.1.or.nmg.gt.12)then
      write(6,*)'Number of magnitude branches outside of acceptable range (1 to 12).'
      write(6,*)'Please check input file near ',adum
      stop'hazFXnga13l: fatal error.'
      endif
      if(cluster)then
      jsegmin=min(jsegmin,jseg(ift))
      jsegmax=max(jsegmax,jseg(ift))
c checks associated with clustering
      igpmax=max(igpmax,igroup(ift))
      if(igroup(ift).lt.1.or.igpmax.gt.5)then
      write(6,*)'igroup outside of acceptable range (1 to 5).'
      write(6,*)'Please check input file near ',adum
      stop'hazFXnga13l: fatal error.'
      endif      !igroup bound check
      if(ift.gt.1.and.nmagf(ift).ne.nmagf(ift-1))then
      write(6,*)'Number of M-recurrence scenarios must be same on all fault segs'
      write(6,*)'For flt ',ift,' nmagf is ',nmagf(ift),' but for previous nmagf was ',nmagf(ift-1)
      stop'hazFXnga13l: fatal error.'
      elseif(ift.eq.1)then
      nscene=nmagf(1)
      endif      !incompatible nmagf
      endif      !if cluster
      if(itype(ift).eq.1) then
      cmnow=3.0      !initialize a characteristic magnitude at a small value
      do imag=1,nmagf(ift)      !new: combine all fault mags
c        write(6,*) "enter char. M, rate, and relative weight"
c Added relwt to each fault-mag line. This relwt is relative to all faults
c in the current input file, not just this particular fault. For example, if
c combining 25-wt and 50-wt California faults into one file, and the relwt for the
c 25-wt items is X, then the relwt for the 50-wt items is 2X. This was discussed
c at a software meeting late Oct 2006. S Harmsen. This relwt should be an epistemic
c uncert weight. For example, related to fault location for NMSZ, or to freq
c of recurrence, to fault dip, or other issues.
        read(1,*) cmag(ift,imag),crate(ift,imag),relwt(ift,imag)
        if(cmag(ift,imag).gt.cmnow)then
        cmnow=cmag(ift,imag)
        mx_index(ift)=imag
        endif
        if(cluster)then
        wtscene(imag,igroup(ift))=relwt(ift,imag)
        if(rate_cl(igroup(ift)).eq.0.)then
        rate_cl(igroup(ift))=crate(ift,imag)
        elseif(crate(ift,imag).ne.rate_cl(igroup(ift)))then
        write(6,*)'Rates of clustered events within a group must all be equal. Trouble:'
        write(6,*) rate_cl(igroup(ift)),crate(ift,imag),' for ift ',ift
        stop 'please correct input file'
        endif !check that event rates are same
        endif      !clustered events
c for clustered events all faults in a group should have the same weight. This is new may 17.        
         xmorate= 10.**(1.5*cmag(ift,imag)+9.05)*crate(ift,imag)
      write(6,*) 'Cmag, rate, morate, rel_wt: ',cmag(ift,imag),crate(ift,imag),xmorate,relwt(ift,imag)
      if(crate(ift,imag).le.tiny)write(6,*)'**** Warning, eq rate <= 1e-12 for this fault. ****'
c mod 0ct 31, 2006 from Oliver Boyd discovery
      if(imag.gt.1)then
      if(abs(cmag(ift,imag)-cmag(ift,imag-1)).gt.0.8)then
      write(6,*)' characteristic magnitudes varying too rapidly.'
      write(6,*)' This seems like an error. Please check input file near ',adum
      stop'hazFXnga13l: fatal error.'
      endif
      endif
        itest(ift)=0      !initially lets not worry about this test for multi-mag
        test= cmag(ift,imag)+dmbranch(1)-mwid*0.05
       if(test.lt.5.8) then
c          itest(ift)=1
           write(6,*) "ale mag < 5.8, but not doing anything about it ",test,itest(ift)
c
c           read(5,*) idum
           endif
        enddo      !imag loop, multi mags.
        endif
      if(itype(ift).eq.2) then
      nmg=nmagf(ift)
      write(6,*) 'Zeng algorithm for floating ruptures w/variable dip direction'
      write(6,*) '**************************#####************************************'
      do imag=1,nmagf(ift)      !new: combine all fault mags
c         write(6,*) "enter a,b,min M,max M,dmag for freq mag curve"
c Added relwt to each fault-mag line. This relwt is relative to all faults
c in the current input file, not just this particular fault. For example, if
c combining 25-wt and 50-wt California faults into one file, and the relwt for the
c 25-wt items is X, then the relwt for the 50-wt items is 2X. This was discussed
c at a software meeting late Oct 2006. S Harmsen.
c Note, a(ift,imag) is an interval rate in this code (i.e., rate of M0+-dmag)
        read (1,*) a(ift,imag),b(ift,imag),magmin(ift,imag),magmax(ift,imag),
     + dmag(ift,imag),relwt(ift,imag)
      write(6,*)'a b magmin magmax relwt',
     + a(ift,imag),b(ift,imag),magmin(ift,imag),magmax(ift,imag),relwt(ift,imag)
c---- center magnitude bins. Follow earlier hazFX code, but generalize
c on dmag. Dmag was 0.1, now dmag is a flexible quantity.
c---- But first, check that dmag(ift,imag) is valid. If not valid, fix it.
       if(dmag(ift,imag).le.0.004)then
       write(6,*)dmag(ift,imag),' invalid dmag for ',adum
       write(6,*)"Code reset this fault's dmag to 0.1"
       dmag(ift,imag)=0.1
       endif
        if(magmin(ift,imag).eq.magmax(ift,imag))nmag0(ift,imag)=1
        if(magmin(ift,imag).ne.magmax(ift,imag)) then
          magmin(ift,imag)= magmin(ift,imag)+dmag(ift,imag)/2.
          magmax(ift,imag)= magmax(ift,imag)-dmag(ift,imag)/2.+.0001      !small delta may be
        if(lb0)then
        magmin(ift,imag+nmg)=magmin(ift,imag)
        magmax(ift,imag+nmg)=magmax(ift,imag)
        dmag(ift,imag+nmg)=dmag(ift,imag)
        relwt(ift,imag)=relwt(ift,imag)*0.5
        relwt(ift,imag+nmg)=relwt(ift,imag)
        print *,magmin(ift,imag+nmg),magmax(ift,imag+nmg),dmag(ift,imag+nmg),relwt(ift,imag+nmg)
        if(sdal.ne.0)then
        print *,'hazFXnga13l does not work with extra GR branches w/b=0 and sdal >0'
        print *,'Please revise input file to have sdal = 0; or revise this code.'
        stop 'hazFXnga13l: unsuccessful input combination'
        endif
        endif
c needed to insure that magmax(ift,imag) > magmin(ift,imag)          
          endif      !magmin.ne.magmax
        nmag0(ift,imag)= int((magmax(ift,imag)-magmin(ift,imag))/dmag(ift,imag) + 1.4)
      if(lb0)nmag0(ift,imag+nmg)=nmag0(ift,imag)
        itest(ift)=0
        test= magmax(ift,imag)+dmbranch(1)
        if((test.lt.6.5).and.(nmag0(ift,imag).gt.1)) then
          itest(ift)=1
c if one of the mag branches yields test < 6.5, then all branches will have itest=1
c this is going to give a different answer than if each mag branch is run separately.
c May want to take another look at this. SHarmsen. Feb 5 2007.
          write(6,*) "test.lt.6.5", ift, itest(ift)
          endif
        if(nmag0(ift,imag).eq.1) then
           test= magmax(ift,imag)+dmbranch(1)-mwid*dma
           if(test.lt.6.5) then
              itest(ift)=1
              write(6,*) "test.lt.6.5 for fault #", ift, ' itest ',itest(ift)
           endif
           endif
c---- calculate moment rate
c initialize to nonzero because a log will be taken
        xmorate=1.0e-10
        do 600 i=1,nmag0(ift,imag)
        xmag= magmin(ift,imag)+ dmag(ift,imag)*(i-1)
        xmorate= xmorate+10.**(a(ift,imag)-b(ift,imag)*xmag+1.5*xmag+9.05)
 600    continue
       if(lb0)then
       xmob0=1.e-20
       do 620 i=1,nmag0(ift,imag+nmg)
        xmag= magmin(ift,imag+nmg)+ dmag(ift,imag+nmg)*(i-1)
620      xmob0=xmob0+10.**(1.5*xmag+ 9.05) 
      afo=xmorate/xmob0
      a(ift,imag+nmg)      = alog10(afo)
      b(ift,imag+nmg) = 0.0
      if(nmag0(ift,imag).eq.1)write(6,*)'xmag, Rate(M0|b0) a ',magmin(ift,imag+nmg),afo,a(ift,imag+nmg)
      xmo2(ift,imag+nmg) = xmorate
      endif       !lb0 true
c xmo2 can be different for each "imag" branch. Needs dimension corresponding
c to this potential variation. SH Nov 2006.
        xmo2(ift,imag)= xmorate
        write(6,*) xmorate,' moment rate J/yr'
c----- find rates for epistemic uncertainty
       if(nmag0(ift,imag).gt.1) then
         do 603 ilt= 1, nbranch
         mmax= magmax(ift,imag)+ dmbranch(ilt)
         nmagg= int((mmax-magmin(ift,imag))/dmag(ift,imag) + 1.4)
         nmagg = max(1,nmagg)      !added Mar 4 2008. 
         sum= 1.e-20
         do 602 m=1,nmagg
         xmag= magmin(ift,imag)+ (m-1)*dmag(ift,imag)
  602    sum= sum + 10.**(-b(ift,imag)*xmag+1.5*xmag+9.05)
         a10 = xmorate/sum
         al10=alog10(a10)
         write(6,*) 'xmag, rate(M0),a ',xmag,a10,al10
         ratenew(ift,imag,ilt)=al10
        
c         write(6,*) "sum2=",sum2
c         read(5,*) idum
  603    continue
         endif      !nmag0(ift,imag).gt.1
c added below aug 28 2007
c----- find rates for epistemic uncertainty, b=0 branch
       if(lb0.and.nmag0(ift,imag+nmg).gt.1) then
         do 623 ilt= 1, nbranch
         mmax= magmax(ift,imag+nmg)+ dmbranch(ilt)
         nmagg= int((mmax-magmin(ift,imag+nmg))/dmag(ift,imag+nmg) + 1.4)
         nmagg=max(1,nmagg)      	!new mar 4 2008 
         sum= 1.e-20
         do 622 m=1,nmagg
         xmag= magmin(ift,imag+nmg)+ (m-1)*dmag(ift,imag+nmg)
  622    sum= sum + 10.**(1.5*xmag+9.05)
         a10= xmorate/sum
         write(6,*) 'xmag, rate(M0|b0),a ',xmag,a10,a(ift,imag+nmg)
         a10=alog10(a10)
c display crossover M, where rate|b0 equals rate|input b. for M>Mcross, the b=0 model 
c predicts more events, whereas for M<Mcross, the b=.8 model (or whatever) predicts more.
         ratenew(ift,imag+nmg,ilt) =a10
      xmag= (ratenew(ift,imag,ilt)-a10)/b(ift,imag)
      write(6,*)'Crossover M beyond which b0 branch rate exceeds other one:',xmag         
c         write(6,*) "sum2=",sum2
c         read(5,*) idum
  623    continue
         endif      !nmag0(ift,imag).gt.1
         enddo      !imag
         if(lb0)then
         nmagf(ift)=2*nmagf(ift)
         print *,'ift ',ift,' revised nmagf(ift)=',nmagf(ift)
         endif
        endif      !itype=2 
c the below line was changed feb 22 2008 Old version is commented out. New vers
c requires both epistemic and aleatory Muncert to be zero to set itest() to 1
      if(sdal.eq.0. .and. itype(ift).eq.1) itest(ift)=1
      if(sdal.eq.0. .and. itype(ift).eq.2 .and. nbranch.eq.1) itest(ift)=1      !line added feb28 2008
c      if(sdal.eq.0.) itest(ift)=1      	!this line seems to remove epistemic branching
c  
      write(6,*)ift,itest(ift),nmag0(ift,1),' fault number, itest,nmag0(ift,1)'
       write(6,*) "enter dip, downdip width, depth0 of fault plane"
      read(1,*) dip(ift), width(ift), depth0(ift)
      write(6,*)'Dip width depth0 ',dip(ift), width(ift), depth0(ift)
      if(abs(dip(ift)).gt.180.)stop'hazFXnga13l: Allowed fault dip in the range -180 to 180 degrees'
      if(depth0(ift).gt.8.)print 2458,depth0(ift),ift
2458      format('***********************************************************',/,
     + '************ Warning: Implausible depth to top of fault ****',
     +/,'************ Depth ',f5.1,' km for fault number ',i3,' ***********')
c      write(6,*) "enter number of segment points"
      read(1,*) npts(ift)
c      write(6,*) "enter lat,lon for each point"
      tlen(ift)= 0.
c nodowndip new Oct 17 2007. Include downdip rupture scenarios only if original
c top of fault is at or near Earth surface. Deep blind thrusts don't need this.
      nodowndip(ift)=depth0(ift).gt.1.      !km. 
c Non-daylighting faults will not have additional downdip tops, just 1 at depth0.
c Below is new apr 3 2007.
c Make dip positive. Distance algorithm is known to give us problems if dip<0.
c Below enforces a standard: look right of strike to find dip direction. 
      if(dip(ift).lt.0..and.dip(ift).ge.-90.)then
      dip(ift)=-dip(ift)
      j1=npts(ift)
      j2=1
      j3=-1
      elseif(dip(ift).gt.90.)then
c Ken Campbell suggests that unconventional fault descriptions should be handled. email mar 6 2008.
      dip(ift)= 180. - dip(ift)
      j1=npts(ift)
      j2=1
      j3=-1
      elseif(dip(ift).lt.-90.)then
c Ken Campbell ditto.
      dip(ift)= 180. - dip(ift)
      j1=1
      j2=npts(ift)
      j3=1
      else
      j1=1
      j2=npts(ift)
      j3=1
      endif 
      do j=j1,j2,j3
      read(1,*)y(j,ift),x(j,ift)
      enddo 
c end new apr 3.    
      do 2 j=1,npts(ift)
c      read(1,*) y(j,ift),x(j,ift)
      if(j.gt.1) then
        sy= y(j-1,ift)
        sx= x(j-1,ift)
        ry= y(j,ift)
        rx= x(j,ift)
        if((sy.eq.ry).and.(sx.eq.rx)) then
             delta=0.
             az=0.
             go to 1020
        else
           call delaz(sy,sx,ry,rx,delta,az,baz)
        endif
 1020   if(j.gt.2) xlen(j-1)= delta +xlen(j-2)
        if(j.eq.2) xlen(1)= delta
        xaz(j-1)= az
        tlen(ift)= tlen(ift)+ delta
        endif
c      write(6,*) xlen(j-1),delta,xaz(j-1)
   2  continue
c----- next line added 10/18 based on Harmsen comment
       xaz(npts(ift))= xaz(npts(ift)-1)
      write(6,*) npts(ift),tlen(ift)
      if(tlen(ift).le.0.) then
          write (6,*) "tlen le 0",ift
          read(5,*) idum
          elseif(tlen(ift).gt.970.*dmove)then
          write(6,*)"Fault length > 970*dmove km. Please increase u,v dim 1"
          stop 'hazFXnga13l: underdimensioned arrays u and v'
          endif
      coef= 3.14159/180.
      dip0(ift)= dip(ift)
c dip(ift) in radians henceforth
      dip(ift) =dip(ift)*coef
      if(poly)then
      zwide = 15.-depth0(ift)
c zwide is the vertical extent of fault=  15 km brittle crust - top depth
c Here is the test to determine if fault is at least partly inside a polygon
      f_is_in(ift)=.false.      !presumed outside initially
      if(ift.eq.1.and..not.deagg)open(29,file='resample.flt',status='unknown')
c the following resampling computes points along bottom of fault to determine if
c downdip part of fault is inside the polygon
      call resample(ift,npts(ift),xlen,xaz,dip(ift),depth0(ift),
     + 2,tlen(ift),dmove,zwide,npts1(ift))
      if(.not.deagg)write(29,55)'# ',ift,npts1(ift)
55      format(a2,i4,1x,i3)      
      do j=1,2
      imax=npts1(ift)
      idel=imax-1
c short version check endpoints only
      do i=1,imax,idel
      f_is_in(ift)= f_is_in(ift).or. LXYIN(u(i,j,ift),v(i,j,ift),PX,PY,npmax)
      if(.not.deagg)write(29,*)u(i,j,ift),v(i,j,ift),depth0(ift)+(j-1)*zwide,f_is_in(ift)
      enddo      !i or along-strike loop 
c      if(.not.deagg)write(29,*)
c      if(.not.deagg)write(29,900)'#'      !gmt marker
      enddo      !j or depth loop
c Note that just because some part of the fault is inside
c the perimeter, this is not the last word. For GR events with partial rupture,
c some of these might be entirely outside the polygon. More checking
c is required for these cases.
      if(.not.f_is_in(ift))then
      write(6,57)adum
57      format('#Fault outside perimeter: ',a,/,'# This source will not be included')
      nft_out = nft_out+1
      endif
      else
      f_is_in(ift)=.true.      !standard presumption when not testing.
       endif      !polygon test
c added March 28 2007: need uniform sampling along the fault for floaters
      if(itype(ift).eq.2)then
      dlen2 = 2.      ! 2 km vertical separation between downdip ruptures
      nleny = 4      ! 4 contours, 1 at depth0, 1 at d+2, 1 at d+4km, 1 at d+6,
c the last is for lowmag faults such as creeping section saf, with M6.2
      call resample(ift,npts(ift),xlen,xaz,dip(ift),depth0(ift),
     + nleny,tlen(ift),dmove,dlen2,npts1(ift))
       write(6,*)'Resampling at ',dmove,' km for fault number ',ift,npts1(ift)
      if(openme)then
      open(49,file='resample.fault',status='unknown')
      openme=.false.
      endif
      write(49,*)'# ',ift,npts1(ift)
      do j=1,3
      do i=1,npts1(ift)
      write(49,*)u(i,j,ift),v(i,j,ift),i,j
      enddo      
      write(49,*)
      enddo
       endif      !iftype=2
   1  continue
 999  nft= ift-1
      write(6,*) "# nft=",nft
      close(49)
      write(6,*)'Fault points w/resampled coords in resample.fault'
      if(poly.and..not.deagg)then
      close(29)
      write(6,*)'fault w/resampled bottom coords in resample.flt'
      endif
      if(nft_out.gt.0)write(6,*)"# Number of excluded faults: ",nft_out
      if(nft.eq.nfltmx)write(6,*)'WARNING: maximum fault index reached. Is this all?'
       do j=1,nft
      firstf(j)=.true.
      enddo
ccccccccccccccccccccccccccccc
c---- write header --------------
c Discussion at Software meeting Dec 16 2005: header record needs to be
c lengthened to include more input-file information.
      do 201 ip=1,nper
      headr%period= period(ip)
      headr%nlev= nlev(ip)
      do 702 k=1,nlev(ip)
 702  headr%xlev(k)= exp(xlev(k,ip))
c 30 character name
c      ndata= 308
c 128 character name
      ndata= 896
c New 10/05: store solution grid info.
      if(sdi)then
      headr%extra(1)=dy_sdi
      else
      headr%extra(1)=-1.
      endif
      headr%extra(2)=xmin
      headr%extra(3)=xmax
      headr%extra(4)=dx
      headr%extra(5)=ymin
      headr%extra(6)=ymax
      headr%extra(7)=dy
      headr%extra(8)=float(nrec)
c store other run-specific info. extra(9) and (10)      
      headr%extra(9)=vs30
      headr%extra(10)=dbasin      !bookkeeping
      if(grid.and.cluster)then
      do ifn=1,ngroup(ip)
      headr%name(5)=g_name(ifn)
      call puthead(ifp(ip,ifn,1),headr,ndata,readn)
      enddo
      elseif(grid)then
      do ifn=1,nfi(ip)
      headr%name(5)=pithy(ifn)
      call puthead(ifp(ip,1,ifn),headr,ndata,readn)
c      print *, "PP: puthead 2 ",ifp(ip,1,ifn)
      enddo
      endif      !if gridded output
 201  continue
ccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccc
      icnt=1
      prob = 0.1e-20      !initialze the array to negligible
      if(cluster)then
      do i=1,ngroup(1)
      write(6,885)i,rate_cl(i)
885      format('group # ',i1,' mean rate = ',e10.5)
      enddo
      write(6,*)'jsegmin jsegmax ',jsegmin,jsegmax
      endif      !write out group weight * rate. added diagnostic june 5 2007
c---Here's the guts
c
c---loop through receiver sites
      do 100 i=1,nrec
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
      ivx=1+nint((rx-vxmin)/vdx)      
      ivy=nint((vymax-ry)/vdy)      !same organization as we are used to.
      vs30=v30(nvx*ivy+ivx)
      endif      !inbounds
      endif      !array rather than scalar vs30      
      if(determ)then
      xmagdet=5.0
      rjbdet=dmax
      do ip=1,nper
      write(66+ip,6684)rx,ry,i,vs30
6684      format('! ',f8.3,1x,f7.3,' site ',i7,' vs30 ',f7.1,' m/s')
      enddo
      endif
c following lines omit several site checks for ca and/or extensional faults.
      if(byeca.and.ry.gt.43.9) goto 860
c      if(byewg99.and.rx.gt.-118.)goto 860      !wg faults 120.5+
      if(byenv.and.(rx.gt.-111.1.or.rx.lt.-122.5))goto 860
      if(byenv.and.(ry.lt.32.1.or.ry.gt.44.1))goto 860
      if(byeext.and.ry.lt.38..and.rx.lt.-121.)goto 860
      if(byeext.and.ry.lt.40..and.rx.lt.-122.)goto 860
      if(byeoreg.and.(ry.lt.40..or.rx.gt.-113.))goto 860
c--- loop through faults
      prob_s = 0.
      do 850 ift=1,nft
c new apr 23 2007: exit if fault not inside polygon
      if(.not.f_is_in(ift))goto 850
      dipang=dip(ift)
      sinedip=sin(abs(dipang))
      if(dip0(ift).le.70..or.dip0(ift).ge.110..or.dip0(ift).lt.0.)then
      cbhwfac=1.
      else
      cbhwfac=abs(90.-dip0(ift))/20.      !changed july 12 2006 SH
      endif
      cosDELTA=cos(dipang)
      cdipsq=cosDELTA**2
      dtor=depth0(ift)
c Heads up: the below factor changes when top of rupture is moved downdip
        cyhwfac=atan(width(ift)*0.5*cosDELTA/(dtor+1.0))/(pi*0.5)
      if(iftype(ift).eq.1)then
      ss=.true.
      rev=.false.
      normal=.false.
      rake=0.
      F_RV=0.0
      F_NM=0.0
      elseif(iftype(ift).eq.2)then
      rev=.true.
      normal=.false.
      ss=.false.
      rake=90.
      F_RV=1.0
      F_NM=0.0
      elseif(iftype(ift).eq.3)then
      normal=.true.
      rev=.false.
      ss=.false.
      F_RV=0.0
      F_NM=1.0
      rake=-90.
      else
c oblique thrust currently assumed. If transtensional, consider -135.
      obl=.true.
      rev=.true.      !check this Oblique in DB?
      ss=.false.
      normal=.false.
      rake=135.
      F_RV=0.5      !intuitive choice. May require colleagues OK.
      F_NM=0.0
      endif
      xdiff=1.e6
      ydiff=1.e6
      do ii=1,npts(ift)
      xdiff= min(xdiff,abs(rx-x(ii,ift)))
      ydiff= min(ydiff,abs(ry-y(ii,ift)))
      enddo
      xdmax= (dmax+200.)/(111.11*cos(ymax*coef))
      ydmax= (dmax+200.)/111.11
      if(xdiff.gt.xdmax) go to 850
      if(ydiff.gt.ydmax) go to 850
      do 35 iatype=1,5
 35   dmin(iatype)=10000.
c Currently, all char events rupture to top of fault. Only gr includes downdip top possibilitty      
c Steve Harmsen, Nov 16 2006. This needs to be discussed. Downdip only applies to
c faults with tops at Earth surface. Oct 2007.
c----- new for AS hanging wall test 
      aspectratio=tlen(ift)/width(ift)      !if rupture occupies entire fault.
c--- calculate distance to site for points (ix,iy) on fault plane
c---loop through  segments 
      hwsite=.false.
c added below dec 14 2006 to perform checks on dipping faults that might have rupture top below
c top-of-fault. If rupture top is below fault top, more dmin checks must be performed
c      if(magmin(ift,1).ge.7.0)then
c      iymx=1
c      elseif(magmin(ift,1).ge.6.5)then
c      iymx=5
c      else
c      iymx=9
c      endif
c      write(6,*)ift,iymx,sinedip,nseg,nlenx,' ift iymx sin(dip) nseg,nlenx'
c do 101 loop entirely eliminated from hazFXnga13l. Just trying char event ruptures (entire fault surface) initially
c with Yuehua Zeng code mindist1. 
c azm and azb are azimuths to a seg boundary not to nearest place on the segment, therefore ,in general, not the best
c choice (when deaggregating).
c Zeng code mindist1 uses a continuous fault surface rather than series of rectangles
c added R_x, distance to top of fault, oct 22 2007
      call mindist1(rcd,azm,rjb,azb,R_x,ry,rx,y(1,ift),x(1,ift),npts(ift),
     +                   dip0(ift),width(ift),depth0(ift))
c       if(rcd.gt.100.)print *,rcd,azm,ift,dip0(ift),width(ift),depth0(ift),rx,ry,
c     + x(1,ift),y(1,ift)
        dmin(1)=rjb
        dmin(2)=rcd; dmin(3) = rcd
        dmin(4)= R_x      !used in CY NGA 11/07
        dmin(5)= R_x*tan(azb*coef)
        if(dmin(5).lt.0.)dmin(5)=1000.      !This step is supposed to nullify footwall. Stopgap dec 10 2012.
c this is the end of the previous do 101 loop. Initially just itype()=1 is programmed.
c        write(6,*) ry,rx,testhw,dhy,angtest,dmin(5)
c  dmin is minimum distance from fault to receiver 
      nearsrc=dmin(1).lt.30.0.and.depth0(ift).lt.1.5      !certain calcs can be bypassed for distant  flts
      if(dmin(1).gt.dmax) go to 850 
      if(itype(ift).eq.2) then
c----  loop through floating rupture zones along fault
c--- find minimum distance dmin2 to receiver for each rupture zone
c--- nrups(m,ift) rupture zones per magnitude bin m
c-- move rupture zones by dlen horizontally along fault
      ntmp=0
      do imag=1,nmagf(ift)
      ntmp=max(ntmp,nmag0(ift,imag))
      enddo
        do 222 m=1,ntmp      !+2 no need for overkill.
        xmag= magmin(ift,1)+(m-1)*dmag(ift,1)      
c magmin: can change slightly with imag. Using first branch value here
c------find rupture length from Wells Coppersmith relation 
c        if(iftype(ift).eq.1) ruplen=  (xmag-5.16)/1.12
c        if(iftype(ift).eq.2) ruplen=  (xmag-5.00)/1.22
c        if(iftype(ift).eq.3) ruplen=  (xmag-4.86)/1.32
c        if(iftype(ift).eq.1) ruplen=  -3.55+0.74*xmag
c        if(iftype(ift).eq.2) ruplen=  -2.86+0.63*xmag
c        if(iftype(ift).eq.3) ruplen=  -3.22+0.69*xmag
c---- changed 5/29/02
c        if(iftype(ift).eq.3) ruplen=  -2.01+0.50*xmag
c here is one place where certain calcs are done only once. New nov 2006
      if(firstf(ift))then
c---- In general, ruplen based on W&C all flt slip types. changed 7/1/02
c--- New Apr 2007: ruplen for CA floaters based on H&B, if input file is CA floater
      if(cal_fl.or.bfault)then
      area1 = 10.**(xmag-4.2)
      area2 = 10.**((xmag-4.07)/0.98)
      if(area1.ge. 500.)then
      ruplen = min(area1/width(ift),tlen(ift))
      else
      ruplen = area2/width(ift)
      endif
c      write (6,*) xmag,ruplen,area1,area2,ift      !temp check.
      else
c general case: ruplen is W&C length based on Mchar      
        ruplen= -3.22+ 0.69*xmag
        ruplen= 10**ruplen
        endif
c down-dip rupture possibilities for M<= 7 . This information
c  xxx not ready for any but top-of-fault to bottom-of-fault rupture. 
c Mods of mar 26 are for entire downdip rupture. Active fault is assumed to be there below these
c downdip tops. Some 2007 California fault models are shallow: Hosgri 0 to 7 km, Bartlett Spr. 0 to 7 km.
c However, these downdip limits do not appear to have been thought out very carefully in some cases.
      if(cal_fl.or.nodowndip(ift))then
      nrupdd =1      !floaters are characteristic events. Occupy full fault width
      elseif(xmag.lt.6.5)then
      nrupdd=4      !smaller events can have deeper tops. new dec 2006. 
      print *, ift,'xmag < 6.5 trying 4 top-of-rupture with Zeng code mindist1'
      elseif(xmag.le.6.75)then
      nrupdd=3
      elseif(xmag.le.7.0)then
      nrupdd=2
      else
      nrupdd=1
      endif
      deltaw = dlen2/sin(dip(ift))      !for 2 km loss in Z, what is reduction in flt width?
      do idn=1,nrupdd
c below guards against the long strip downdip in case M is less than 6.2 or so.
      widthrup=max(1.0, min(width(ift)-deltaw*float(idn-1),ruplen))
        ar(m,idn)=ruplen/widthrup
        enddo      !variable aspect ratio for floating rupture downdip.
c--- changed to dmove move for floating rupture
        nrups(m,ift)= int((tlen(ift)-ruplen)/dmove) +1
        nrupd(m,ift)=nrupdd
      rupln(m,ift)=min(ruplen,tlen(ift))
        if(ruplen.ge.tlen(ift))nrups(m,ift)=1
c the above calcs are done once per fault, when firstf(ift) is .true.      
      endif
c do jrup below is the new downdip rupture loop: if M<7 then nrupd>1.
c nrupd is 4 for M<6.5. Maximum top of rupture 8 km deep in that case.
c nrup is nrups * nrupd (number along strike * number downdip)
c simple rectangle is assumed. Does this have problems for complicated geom?
c
        iy0 = 1
c each station -source pair has calculation branch that is based on "nearsrc"
        if(nearsrc)then
        nrupdd=nrupd(m,ift)
        else
        nrupdd=1
c nearsource sensitivity. far away ASSUME no sensitivity to depth of rupture top.
c Note: this has been tested and found correct for CA GR sources.
      endif
      nxxp = rupln(m,ift)/dmove +1
      deltaw = dlen2/sin(dip(ift))      !for 2 km loss in Z, what is reduction in flt width?
        do jrup=1,nrupdd
        ztop = depth0(ift)+dlen2*(jrup-1)
c rupture width decreases as top of rupture goes deeper. Bottom is fixed at Z=15km
c or other fault bottom except for skinny faults (CALIF). Widthp is width of
c rupture at this time. Minimum width 2 km.
        widthp = max(2.0,width(ift)-deltaw*(jrup-1))
        do 302 irup=1,nrups(m,ift)
        if(jrup.eq.1)segin=.false.
c        xmove= dmove*(irup-1)
c The below process will only work for faults sampled at dmove km increments.
        ix0= irup
        ix1=min(nxxp+ix0,npts1(ift))      
      nx1=1+ix1-ix0
c ix1 is supposed to be index of last point on rup seg. nx1 is the number of points
c in the segment.
c Next: Another check if polygon option is on. Added May 25 2007 SH.
c This test checks whether the top of flt segment is inside the polygon
      if(poly.and.jrup.eq.1)then
      do ix=ix0,ix1
      segin=segin.or. LXYIN(u(ix,1,ift),v(ix,1,ift),PX,PY,npmax)
      enddo
      endif
c Code ready for downdip ruptures, Oct 22 2007
        call mindist1(rcd,azm,rjb,azb,R_x,ry,rx,v(ix0,jrup,ift),u(ix0,jrup,ift),nx1,
     +                   dip0(ift),widthp,ztop)
c       if(rcd.gt.100.)print *,rcd,azm,ift,dip0(ift),width(ift),depth0(ift),rx,ry,
c     + v(ix0,jrup,ift),u(ix0,jrup,ift)
      if(poly.and..not.segin)then
c how far away is the moon?
      dminr(1,jrup)=200000.
      dminr(2,jrup)=200000.
      dminr(3,jrup)=200000.
      else
        dminr(1,jrup)=rjb
        dminr(2,jrup)=rcd
        dminr(3,jrup)=rcd
        dminr(4,jrup)= R_x      !used in CY NGA 11/07
        dminr(5,jrup)=R_x*tan(azb*coef)         !used in AS12. Dec 2012.   
c index 5: wings were clipped, winged metric no longer available. SH Mar 2007.
      endif
c----- Hanging wall test for AS97. Skip it if as97 is not invoked. Skip it.
c        if(dip0(ift).eq.90..or.noas97) go to 304
c        if(nyy2.le.(nyy+2)) go to 304
c        do 3031 iy=nyy+2, nyy2
c        dist= fdist(ix,iy,5)
c        if(dist.lt.dminr(5,jrup)) then
c         dminr(5,jrup)= dist
c         dhy= (iy-(nyy+2))*dlen*cos(abs(dip(ift)))
c         testhw3(irup,m)= dhy*angtest+ 2.
c         endif
c 3031   continue
c-----------------------
  304   continue
        do 39 iatype=1,4      !type 4 added oct 23 2007
 39     dmin2(m,irup,jrup,iatype)=dminr(iatype,jrup)
  302   continue
      iy0=iy0+2
      enddo      !downdip index
  222   continue
      firstf(ift)=.false.
c end of down-dip related changes
      endif
c 
c
c-----------------------------------------
c----  GR with nmag0>1 with Mmax uncertainties
      if((itype(ift).eq.2).and.(nmag0(ift,1).gt.1).and.
     &   (itest(ift).eq.0)) then
      do imag=1,nmagf(ift)      !can have two or more distributions often
c with diffrent Mmax.
      dtor= depth0(ift)
      if(cluster)then
      ks=imag
      kg=igroup(ift)
      else
      ks=1
      kg=1
      endif
       do 403 ilt=1,nbranch
        mmax= magmax(ift,imag)+ dmbranch(ilt)
        totmo= 0.
        nmag= (mmax- magmin(ift,imag))/dmag(ift,imag) + 1.4
c-- loop through magnitudes
c       if(abs(rx+101.).lt.0.1.and.abs(ry-38.).lt.0.1)then
c       if(ilt.eq.1)write(6,*)rx,ry,' gr with uncrt ift',ift
c       write(6,*) "nmag=",nmag,"nmag0(ift,imag)",nmag0(ift,imag),nrups(m,ift)
c       endif
        do 285 m=1,nmag
      isbig=nmag.eq.m
        xmag= magmin(ift,imag)+(m-1)*dmag(ift,imag)
        rate= ratenew(ift,imag,ilt) - b(ift,imag)*xmag
c rate is now a relative rate, because of new relwt factor Nov 2 2006. SH.
        rate= relwt(ift,imag) * 10.**rate
        if(rate.lt.tiny)goto 403      !bypass irrelevant models aug 23 2007
        totmo= totmo+ rate*(10.**(1.5*xmag+ 9.05))
c        torate(mm)= torate(mm)+ rate*ale(idis+mwid+1)
c cluster model new may 16 2007. Base cluster on mag and geographic branch for Char. For GR
c no cluster model has been specifically discussed
      if(cluster)then
      ks=1
      kg=1
      else
      ks=1
      kg=1
      endif
c dM bins have width 0.1 here 
                  if(deagg)then
                  im=nint((xmag-5.7)*10.)+1
                  im=min(30,im)      !current dimension 30
                  im=max(im,1)
                  immax=max(immax,im)
                  endif
c new, possible downdip top surface of rupture. nrupd can be 1, 2 or 3 
      if(nearsrc)then
      nrupdd=nrupd(m,ift)
      else
      nrupdd=1
      endif 
        rate= rate/float(nrups(m,ift)*nrupdd)
          if(xmag.lt.mcut(1))then
          ime=1
          elseif(xmag.lt.mcut(2))then
          ime=2
          else
          ime=3
          endif
      if(slist)write(30+i,8899)dmin(1),nearsrc,nrupd(m,ift)
8899      format('#',t30,f6.1,'km; nearsrc=',l4,1x,i4,'=nrupd(m,ift)')
      do 284 jrup=1,nrupdd      !new loop downdip rup.
       norpt=.true.      !can have a report for rupture at each downdip level
        aspectratio=ar(m,jrup)
        dz=dlen2*float(jrup-1)      !dlen2 is probably 2 km
        dtor1=dtor + dz      !all important depth to top of rupture.
        W_rup=width(ift)-dz
        cyhwfac=atan(W_rup*0.5*cosDELTA/(dtor1+1.0))/(pi*0.5)
c--- loop through floating rupture zones
        do 284 irup=1,nrups(m,ift)
c report rate of eqs * weight applied to logic-tree branch, R,M,rate'
        rjb=dmin2(m,irup,jrup,1)
        isclose=abs(rjb-dmin(1)).lt.0.95
          if(rjb.lt.dcut(1))then
          ide=1
          elseif(rjb.lt.dcut(2))then
          ide=2
          else
          ide=3
          endif
c rjb for segments not inside the polygon has been set to 200000 km. Code has
c not yet checked for r < dmax so this has to be done here. New May 25.
      if(slist.and.rjb.lt.10000.)
     +  write(30+i,399)rjb,xmag,rate*wtbranch(ilt),ift,dtor1,irup,ilt
399      format(f7.2,1x,f6.2,1x,e12.5,1x,i3,1x,f5.1,3(1x,i3))
        do  ip=1,nper
        iq=iper(ip)      !std periods in 2002. 
        idet=ip+66
        R_x=dmin2(m,irup,jrup,4)      !type 4 for the new R_x
        do 282 ia=1,nattn(ip)
        rrup= dmin2(m,irup,jrup,3)      !type 3 supposed to be rrup
        Ry0= dmin2(m,irup,jrup,5)      !type 5 new dec 2012.
        weight= wt(ip,ia,1)
        if(rjb.gt.wtdist(ip,ia)) weight=wt(ip,ia,2)
        if(rjb.gt.dmax) go to 282
      if(deagg)then
        if(rrup.le.100.)then
        ir=0.1*rrup + 1
        elseif(rrup.le.300.)then
        ir=11+(rrup-100.)/50.
        elseif(rrup.le.500.)then
        ir=15
        else
        ir=16
        endif 
       endif      !if deagg
      ipia=iatten(ip,ia)
c nga relations are collected into the first set of models to check
      if(nga(ip,ia)) then
        if(ipia.eq.13) then
      call getBooreNGA07m(ip,iperb(ip),xmag,rjb,vs30,gnd,sigmaf)
         elseif(ipia.eq.14) then
      pgacalc=ip.eq.1.or.period(ip).eq.0.
        call getCBNGA1107 (ip,ipercb(ip),xmag,rrup,rjb,dtor1,
     1 vs30,dbasin,gnd,sigmaf)
      elseif(ipia.eq.15)then
c 10/2007 implementation of getCYNGA 
c        call CY2006H(ip,ipercy(ip), xmag, width(ift), rrup, rjb, vs30,  
c     1      dtor1, srcSiteA, F_RV, F_NM,gnd,sigmaf)
      call CY2007H(ip,ipercy(ip), xmag, rrup, rjb, R_x, vs30,z1, dtor1,
     1                 F_RV, F_NM,  gnd, sigmaf)
      elseif(ipia.eq.33)then
      indx1=indbssa(ip)
      mech = ibtype(ift)
      iregion=1	!CA-Taiwan version
      call bssa2013drv(ip,indx_pga,indx1,xmag,rjb,vs30,z1km,mech,iregion,alny,sigmaf)
c skip the 2nd period for now assume the period of interest was tabulated.
      gnd(1)=alny
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
          if(ift.eq.-1)then
          print *,ip,exp(alny),rjb,xmag,sigmaf
          print *, per1, xmag, rjb, r1,' per1, xmag, rjb, r1,sigmaf'
      endif
      elseif(ipia.eq.37)then
      call  getIdriss2012(idriss(ip),ip,xmag,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.38)then
      call gksa12(ipgk(ip),ip,xmag,rrup,gnd,sigmaf,vs30,iftype(ift),dgkbasin,Q_CA)
      elseif(ipia.eq.39)then
      ib=0	!basin term Graizer & kalkan model
      call  gksa(ipgk(ip),ip,xmag,rrup,gnd,sigmaf,vs30,iftype(ift),IB)
      elseif(ipia.eq.35)then
c from BChiou email use fmeasured=1. dtor1 is top of floating rupture may be downdip
      F_Measured=1.
      F_Inferred=0.0	!could change these. Reference rock.
c cDpp new predictor variable 3/2013: 
c      call  CY2012_NGA(ip,icy12(ip), xmag, rrup, rjb, R_x, vs30,
c     1                   F_Measured, F_Inferred, z1, dip0(ift), dtor1,
c     1                   F_RV, F_NM,
c     1                   psa_ref, gnd, tau, sigma_NL0, sigmaf)
       call CY2013_NGA(ip,icy13(ip), xmag, rrup, rjb, R_x, vs30,
     1                   F_Measured, F_Inferred, deltaZ1, dip0(ift), dtor1,
     1                   F_RV, F_NM, cDPP(ift), gnd, tau, sigma_NL0, sigmaf)
c      print *,gnd(1),ift
      elseif(ipia.eq.36)then
c     Compute SA1100 for AS2013. New Mar 1 2013. @W_rup is the rupture width (km). Downdip ruptures.
      SA_rock = 0.
c hanging-wall flag for as12 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
      call AS2013_v10_model ( ipas13(ip),xmag,dip0(ift), W_rup, dtor1, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      Sa1100 = exp(lnSa)
c 
c     Compute Sa at spectral period for given Vs30
      call AS2013_v10_model ( ipas13(ip),xmag, dip0(ift), W_rup, dtor1, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30, SA1100, z1km, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      gnd(1)=lnSa
c      print *,exp(lnSa),phi,tau,width(ift)
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
        sigmaf = 1./sqrt( phi**2 + tau**2 )/sqrt2

        elseif (ipia.eq.16) then
c hanging-wall flag for as08 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif

c Although iflag makes a token appearance in SR code, it isn't used. So I dropped iflag (SH).
       call AS_072007 ( ip,ipera(ip),xmag, dip0(ift), F_NM,F_RV, W_rup, rRup, rjb,R_x,
     1                     vs30, hwflag, gnd(1), sigma1, dtor1,  vs30_class,
     3                     z1km )
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
c gnd is lnSA in AS code. need to fill in some other input. vs30_class ?
c period1 is output by AS_modelXX. period1 should equal per(ip)
       sigmaf=1./sqrt2/sigma1
      elseif(ipia.eq.34)then
      SJ=0.0	!Are we in Japan? 1 = yes.
      Hhyp=min(12.,dtor1+8.)	!Take a stab at hypocenter depth.
      call CB13_NGA_SPEC  (ip,ipcb12(ip), xmag,Rrup,R_x,Rjb,F_RV, F_NM,dtor1,Hhyp,W_rup,dbasin,Vs30,
     +Dip0(ift),SJ,gnd,sigmaf)
c      print *,exp(gnd(1)),ift,xmag,Rrup,period(ip)
      elseif(ipia.eq.17)then
      call getIdriss(ip,iper(ip),xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.18)then
c new nov 7,2006. Kanno
      call kanno(ip,iperk(ip),gnd,sigmaf,xmag,rrup,vsfac)
      endif
      elseif(wus02(ip,ia))then
c WUS pre-nga (2002 atten models)
      if(ipia.eq.1)then
      call getSpu2000(ip,iq,xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.3)then
      call getsadigh(ip,iq,xmag,rrup,gnd,sigmaf)
      elseif(ipia.eq.8)then
      call getAS97(ip,iq,xmag,rrup,hwsite,gnd,sigmaf)
       elseif(ipia.eq.9)then
       call getCamp2003(ip,iq,xmag,rrup,rjb,gnd,sigmaf)
      elseif(ipia.eq.11)then
      call getBJF97(ip,iq,xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.-3)then
        call SomerIMW(ip,isomer(ip), xmag,  rrup, rjb,vs30,gnd,sigmaf)
      endif
      elseif(ceus02(ip,ia).or. ceus11(ip,ia))then
c CEUS pre-nga (2002 atten models) added 3 new CENA models defined by tables. Mar 2011
       rjbp=max(rjb,0.11)
       if(ipia.eq.25)then
      rkm=rjbp
      jf=ia08(ip)
             ka=1
             sigma = 0.3*2.303
       gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
       sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.26)then
       rkm=rrup
       jf=ia06(ip)
       ka=2
             sigma = 0.3*2.303
       gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
       sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.27)then
       rkm=rrup
       jf=ip11(ip)
c 2nd variable in sigPez11 is the frequency index according to the GailA set 2011.
      sigma = sigPez11(xmag,jf)	!new 10/30/2012. SH.
       ka=3
       gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
       sigmaf = 1./sqrt2/sigma
c the rest of these are pre-2011 models
       elseif(ipia.eq.2)then
      call getToro(ip,iq,1,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-2)then
      call getToro(ip,iq,2,xmag,rjb,gnd,sigma,sigmaf)
        elseif(ipia.eq.4)then
        call getAB06(ip,iperab(ip),irab(ip,ia),xmag,rrup,gnd,sigma,sigmaf,vs30)
        elseif(ipia.eq.19)then
      call getTP05(ip,ipertp(ip),irtb,xmag,rrup,gnd,sigma,sigmaf)
        elseif(ipia.eq.20)then
      call getSilva(ip,isilva(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
         elseif(ipia.eq.28)then
       call getSilvaV(ip,isilva2(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.6)then
      call getFEA(ip,iq,1,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-6)then      
      call getFEA(ip,iq,3,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.5)then      !to distinguish ab95 from fea
      call getFEA(ip,iq,2,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-5)then      !HR for atkinson-boore
      call getFEA(ip,iq,4,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.7) then
        call getSomer(ip,iq,1,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-7) then
        call getSomer(ip,iq,2,xmag,rjb,gnd,sigma,sigmaf)
             elseif(ipia.eq.10)then
        call getCampCEUS(ip,iq,1,xmag,rrup,gnd,sigma,sigmaf)
             elseif(ipia.eq.-10)then
        call getCampCEUS(ip,iq,2,xmag,rrup,gnd,sigma,sigmaf)
      endif
c CEUS only: This is a good place to apply the extra clamp constraint which
c in 2002 was used only for CEUS sources. The precomputed p() array does
c not work for this case because truncation is no longer at mu+3sig
      test0=gnd(1)+3.*sigma
      test= exp(test0)
      if(clamp(iq).lt.test .and. clamp(iq).gt.0.) then
      clamp2= alog(clamp(iq))
      sigmasq=1./sigma/sqrt2
      tempgt3= (gnd(1) - clamp2)*sigmasq
      probgt3= (erf(tempgt3)+1.)*0.5
      do  k=1,nlev(ip)
      temp= (gnd(1) - xlev(k,ip))*sigmasq
      temp1= (erf(temp)+1.)*0.5
      temp1= (temp1-probgt3)/(1.-probgt3)
      if(temp1.lt.0.) goto 282      !no more calcs once p<0
        prob(icnt,k,ip,1,1)= prob(icnt,k,ip,1,1)+ 
     &         wtbranch(ilt)*weight*rate*temp1
        enddo
        goto 282
        endif
      elseif(ipia.eq.12)then
      call getMota(ip,iper(ip),xmag,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.22)then
      call getDahle95(iperdahl(ip),isoild,xmag,rjb,gnd,sigmaf)
      endif
      wttmp=wtbranch(ilt)*weight*rate
        do ifn=1,nfi(ip)
c for deterministic calcs write the largest M cases from GR distributions.
        if(determ.and.isbig.and.isclose.and.norpt(ip,ia))then
        write(idet,679)exp(gnd(ifn)),1./sigmaf/sqrt2,ipia,ift,xmag,rjb,rrup,
     +     wttmp,dtor1,iftype(ift), ifn

c there is now a report dont tell it again, Sam. (you could have more than one
c that is big and close but these are all the same w.r.t. saved params)
      if(ifn .eq. nfi(ip)) norpt(ip,ia)=.false.
      endif
679      format(e11.5,1x,e11.5,1x,i2,1x,i3,1x,f6.2,1x,f7.1,1x,f7.1,1x,e11.5,
     + f5.1,1x,i2,1x,i1)
      do  k=1,nlev(ip)
       pr=(gnd(ifn) - xlev(k,ip))*sigmaf
       if(pr.gt.3.3)then
       ipr=250002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 283      !transfer out if ground motion above mu+3sigma
       endif
        cfac= wttmp*p(ipr)
c kg index is present: 2nd to last frontier.
        prob(icnt,k,ip,ifn,kg)= prob(icnt,k,ip,ifn,kg)+ cfac
        if(deagg)then
c Deagg is set up for one ground motion per period. Thus no k index here.
c for individual fault
       eps= -pr*sqrt2
      frbar(ift,ifn,ip)=frbar(ift,ifn,ip)+cfac*rrup
      fmbar(ift,ifn,ip)=fmbar(ift,ifn,ip)+cfac*xmag
      febar(ift,ifn,ip)=febar(ift,ifn,ip)+cfac*eps
      fhaz(ift,ifn,ip)=fhaz(ift,ifn,ip)+cfac
      if(eps.lt.emax)then
      ieps=max(1,min(int((eps+2.)*2.),10))
      rbar(ir,im,ieps,ifn,ip)=rbar(ir,im,ieps,ifn,ip)+cfac*rrup
      mbar(ir,im,ieps,ifn,ip)=mbar(ir,im,ieps,ifn,ip)+cfac*xmag
      ebar(ir,im,ieps,ifn,ip)=ebar(ir,im,ieps,ifn,ip)+cfac*eps
      haz(ir,im,ieps,ifn,ip)=haz(ir,im,ieps,ifn,ip)+cfac
      if(eps.ge.2.)goto 283
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
      ka=5
      temp=p(ipr)
      prx= temp-ptail(ka)
      dowhile(prx.gt.1.e-9)
      prob5(ir,im,ieps,ip,ifn,ka)=prob5(ir,im,ieps,ip,ifn,ka)+
     +       wttmp*prx      ! *prlr 	dont know about this
      ka=ka-1
      if(ka.eq.0)goto 283
      prx= temp-ptail(ka)
      enddo
      endif      !eps <emax
      endif      !deagg
         enddo      !levels
283      continue
      enddo      !files
 282    continue
       enddo      ! this used to be a do 283,...283 continue loop
 284    continue
 285    continue

 403    continue
       enddo      !imag loop, diff. distributions.
        go to 850
        endif
ccccccccccccccccccccccccccccc
ccc---- GR with nmag0=1, one magnitude floated, with uncertainties,
c       used for some faults outside of CA and in 2002 Maacama northern CA
ccc--- for very long faults when Mmmax set to 7.5. In 2007 Maacama char has itype 1 
      if((itype(ift).eq.2).and.(nmag0(ift,1).eq.1).and.
     &  (itest(ift).eq.0)) then
      do imag=1,nmagf(ift)      !can have two or more distributions often
c with different Mmax.
        totmo=0.
        aspectratio=ar(1,1)      
c        write(6,*)ift,aspectratio,' one mag floated'
        do 1406 ilt=1,nbranch
        isbig=ilt.gt.nbranch/2      !smaller M-branch omit from determ report
        xmag= magmax(ift,imag)+ dmbranch(ilt)
c cluster model new may 16 2007. Base cluster on mag and geographic branch for Char. For GR
c no cluster model has been specifically discussed
      if(cluster)then
      ks=1
      kg=1
      else
      ks=1
      kg=1
      endif
        xmo= 10.**(1.5*xmag+9.05)
c relative rate, nov 2006
        rate= relwt(ift,imag) * xmo2(ift,imag)/xmo
        if(rate.lt.tiny)goto 1406      !bypass irrelevant models aug 23 2007
        rate= rate/float(nrups(1,ift))
        do 1403 idis= -mwid, mwid
        xmag2= xmag+ idis*dma      !note change from 0.05 to the variable dma
c      print *,idis,xmag2,rate
          if(xmag2.lt.mcut(1))then
          ime=1
          elseif(xmag2.lt.mcut(2))then
          ime=2
          else
          ime=3
          endif
          if(deagg)then
          im = nint((xmag2-5.7)*10.) +1
          im= max(1,im)
          im= min(30,im)
          endif
        wt1wt2=ale2(idis+mwid+1)*wtbranch(ilt)
c--- loop through floating rupture zones
      norpt=.true.      !Array ASSIGNMENT (ip, ia)
        do 1284 irup=1,nrups(1,ift)
c        totmo= totmo+ ale2(idis+mwid+1)*rate*10.**(1.5*xmag2+ 9.05)
C      write rate data if station list option is in effect.
c for big long faults assume all ruptures to surface... Discuss with colleagues?
c this case can happen for short faults such as creeping section. could have down-
c  dip rupture on these smaller sources.
        rjb=dmin2(1,irup,1,1)
      xdif=abs(rjb-dmin(1))
      isclose=xdif.lt.0.95.or.(rjb.gt.xdif.and.xdif.lt.2.5)
          if(rjb.lt.dcut(1))then
          ide=1
          elseif(rjb.lt.dcut(2))then
          ide=2
          else
          ide=3
          endif
      if(slist.and.rjb.lt.10000.)
     + write(30+i,399)rjb,xmag2,rate*wt1wt2,ift,depth0(ift),irup,ilt,idis
        do  ip=1,nper
        idet=ip+66
        iq=iper(ip)
        R_x = dmin2(1,irup,1,4)      !new signed distance oct 2007
        do 1282 ia=1,nattn(ip)
        rrup= dmin2(1,irup,1,3)
        Ry0= dmin2(1,irup,1,5)      !type 5 new dec 2012.
        weight= wt(ip,ia,1)
        if(rjb.gt.wtdist(ip,ia)) weight=wt(ip,ia,2)
        if(rjb.gt.dmax) go to 1282
      if(deagg)then
        if(rrup.le.100.)then
        ir=0.1*rrup + 1
        elseif(rrup.le.300.)then
        ir=11+(rrup-100.)/50.
        elseif(rrup.le.500.)then
        ir=15
        else
        ir=16
        endif 
       endif      !if deagg
      ipia=iatten(ip,ia)
c
c nga relations are collected into the first set of models to check
      if(nga(ip,ia)) then
        if(ipia.eq.13) then
      call getBooreNGA07m(ip,iperb(ip),xmag2,rjb,vs30,gnd,sigmaf)
         elseif(ipia.eq.14) then
      pgacalc=ip.eq.1.or.period(ip).eq.0.     !campbell pga calc
        call getCBNGA1107 (ip,ipercb(ip),xmag2,rrup,rjb,dtor,
     1 vs30,dbasin,gnd,sigmaf)
      elseif(ipia.eq.15)then
c 3/2008 implementation of getCYNGA 
c        call CY2006H(ip,ipercy(ip), xmag2, width(ift), rrup, rjb, vs30,  
c     1      dtor, srcSiteA, F_RV, F_NM,gnd,sigmaf)
      call CY2007H(ip,ipercy(ip), xmag2, rrup, rjb, R_x, vs30,z1, dtor,
     1                 F_RV, F_NM,  gnd, sigmaf)
c these particular floating ruptures always rupture to top of fault, dtor.
      elseif(ipia.eq.33)then
      indx1=indbssa(ip)
      mech=ibtype(ift)
      iregion=1	!CA-Taiwan version
      call bssa2013drv(ip,indx_pga,indx1,xmag2,rjb,vs30,z1km,mech,iregion,alny,sigmaf)
c skip the 2nd period for now assume the period of interest was tabulated.
      gnd(1)=alny
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
          if(ift.eq.-1)then
          print *,ip,y1,rjb,xmag2,sigmaf
          print *, per1, xmag2, rjb, r1,' per1, xmag, rjb, r1'
        print *,  mech,vs30,y1,expsiglny,' mech,vs30,y1,expsiglny'
      endif
      elseif(ipia.eq.37)then
      call  getIdriss2012(idriss(ip),ip,xmag2,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.38)then
      call gksa12(ipgk(ip),ip,xmag2,rrup,gnd,sigmaf,vs30,iftype(ift),dgkbasin,Q_CA)
      elseif(ipia.eq.39)then
      ib=0	!basin term Graizer & kalkan model
      call  gksa(ipgk(ip),ip,xmag2,rrup,gnd,sigmaf,vs30,iftype(ift),IB)
c      if(ift.eq.1)print *,ip,exp(gnd(1)),rrup,xmag,1/sigmaf/sqrt2
      elseif(ipia.eq.35)then
c from BChiou email use fmeasured=1.
      F_Measured=1.
      F_Inferred=0.0	!could change these. Reference rock.
c      call  CY2012_NGA(ip,icy13(ip), xmag2, rrup, rjb, R_x, vs30,
c     1                   F_Measured, F_Inferred, z1, dip0(ift), dtor,
c     1                   F_RV, F_NM,
c     1                   psa_ref, gnd, tau, sigma_NL0, sigmaf)
       call CY2013_NGA(ip,icy13(ip), xmag2, rrup, rjb, R_x, vs30,
     1                   F_Measured, F_Inferred, deltaZ1, dip0(ift), dtor,
     1                   F_RV, F_NM, cDPP(ift), gnd, tau, sigma_NL0, sigmaf)

c       print *,gnd(1),ift
      elseif(ipia.eq.36)then
c     Compute SA1100 for AS2013. New Mar 1 2013.W_rup
      SA_rock = 0.
c hanging-wall flag for as12 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
      call AS2013_v10_model ( ipas13(ip),xmag2,dip0(ift), Width(ift), dtor, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      Sa1100 = exp(lnSa)
c 
c     Compute Sa at spectral period for given Vs30
      call AS2013_v10_model ( ipas13(ip),xmag2, dip0(ift), Width(ift), dtor, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30, SA1100, z1km, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      gnd(1)=lnSa
c      print *,exp(lnSa),phi,tau,width(ift),xmag2,Width(ift)
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
        sigmaf = 1./sqrt( phi**2 + tau**2 )/sqrt2
      elseif(ipia.eq.34)then
      SJ=0.0	!Are we in Japan? 1 = yes.
      Hhyp=min(12.,dtor+8.)	!Take a stab at hypocenter depth.
      call CB13_NGA_SPEC  (ip,ipcb12(ip), xmag2,Rrup,R_x,Rjb,F_RV, F_NM,dtor,Hhyp,Width(ift),dbasin,Vs30,
     +Dip0(ift),SJ,gnd,sigmaf)
c      print *,exp(gnd(1)),ift,xmag2,Rrup
       elseif (ipia.eq.16) then
c hanging-wall flag for as08 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif

c Although iflag makes a token appearance in SR code, it isn't used. So I dropped iflag (SH).
       call AS_072007 ( ip,ipera(ip),xmag2, dip0(ift), F_NM,F_RV, Width(ift), rRup, rjb,R_x,
     1                     vs30, hwflag, gnd(1), sigma1, dtor,  vs30_class,
     3                     z1km )
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
c gnd is lnSA in AS code. need to fill in some other input. vs30_class ?
c period1 is output by AS_modelXX. period1 should equal per(ip)
       sigmaf=1./sqrt2/sigma1
c gnd is lnSA in AS code. need to fill in some other input. esp. aspectratio.
c period1 is output by AS_modelXX. period1 should equal per(ip)
      sigmaf=1./sqrt2/sigma1
      elseif(ipia.eq.17)then
      call getIdriss(ip,iper(ip),xmag2,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.18)then
c new nov 7,2006. Kanno
      call kanno(ip,iperk(ip),gnd,sigmaf,xmag2,rrup,vsfac)
      endif
      elseif(wus02(ip,ia))then
c WUS pre-nga (2002 atten models)
      if(ipia.eq.1)then
      call getSpu2000(ip,iq,xmag2,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.3)then
      call getsadigh(ip,iq,xmag2,rrup,gnd,sigmaf)
      elseif(ipia.eq.8)then
      call getAS97(ip,iq,xmag2,rrup,hwsite,gnd,sigmaf)
       elseif(ipia.eq.9)then
       call getCamp2003(ip,iq,xmag2,rrup,rjb,gnd,sigmaf)
      elseif(ipia.eq.11)then
      call getBJF97(ip,iq,xmag2,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.-3)then
        call SomerIMW(ip,isomer(ip), xmag2,  rrup, rjb,vs30,gnd,sigmaf)
      endif
      elseif(ceus02(ip,ia).or.ceus11(ip,ia))then
c CEUS  (2002 and 2006 atten models)
c
c added 3 new CENA models defined by tables. Mar 2011. 25&26 are Gail Atkinson. 27 is Pezeshk
       rjbp=max(rjb,0.11)
       if(ipia.eq.25)then
      rkm=rjbp
       jf=ia08(ip)
       ka=1
             sigma = 0.3*2.303
       gnd(1) = amean11(xmag2,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
      sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.26)then
       rkm=rrup
       jf=ia06(ip)
       ka=2
             sigma = 0.3*2.303
       gnd(1) = amean11(xmag2,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
      sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.27)then
       rkm=rrup
       jf=ip11(ip)
      sigma = sigPez11(xmag2,jf)
       ka=3
       gnd(1) = amean11(xmag2,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
       sigmaf = 1./sqrt2/sigma

      elseif(ipia.eq.2)then
      call getToro(ip,iq,1,xmag2,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-2)then
      call getToro(ip,iq,2,xmag2,rjb,gnd,sigma,sigmaf)
        elseif(ipia.eq.4)then
        call getAB06(ip,iperab(ip),irab(ip,ia),xmag2,rrup,gnd,sigma,sigmaf,vs30)
        elseif(ipia.eq.19)then
      call getTP05(ip,ipertp(ip),irtb,xmag2,rrup,gnd,sigma,sigmaf)
        elseif(ipia.eq.20)then
      call getSilva(ip,isilva(ip),irsilva,xmag2,rjb,gnd,sigma,sigmaf)
         elseif(ipia.eq.28)then
       call getSilvaV(ip,isilva2(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.6)then
      call getFEA(ip,iq,1,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-6)then      !to distinguish ab95 from fea
      call getFEA(ip,iq,3,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.5)then      !HR
c hard rock possibilities. 
      call getFEA(ip,iq,2,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-5)then      !HR for atkinson-boore
      call getFEA(ip,iq,4,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.7) then
        call getSomer(ip,iq,1,xmag2,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-7) then
        call getSomer(ip,iq,2,xmag2,rjb,gnd,sigma,sigmaf)
             elseif(ipia.eq.10)then
        call getCampCEUS(ip,iq,1,xmag2,rrup,gnd,sigma,sigmaf)
             elseif(ipia.eq.-10)then
        call getCampCEUS(ip,iq,2,xmag2,rrup,gnd,sigma,sigmaf)
      endif
c This is a good place to apply the extra clamp constraint which
c in 2002 was used only for CEUS sources. The precomputed p() matrix does
c not work for this case because truncation is no longer at mu+3sig
      test0=gnd(1)+3.*sigma
      test= exp(test0)
      if(clamp(iq).lt.test .and. clamp(iq).gt.0.) then
      clamp2= alog(clamp(iq))
      sigmasq=1./sigma/sqrt2
      tempgt3= (gnd(1) - clamp2)*sigmasq
      probgt3= (erf(tempgt3)+1.)*0.5
      prr= 1./(1.-probgt3)
c GR-region is not ready for clustered event. KG index (last dim of prob) is 1.
      do  k=1,nlev(ip)
      temp= (gnd(1) - xlev(k,ip))*sigmasq
      temp1= (erf(temp)+1.)*0.5
      temp1= (temp1-probgt3)*prr
      if(temp1.lt.0.) goto 1282      !no more calcs once p<0
        prob(icnt,k,ip,1,1)= prob(icnt,k,ip,1,1)+ wt1wt2*weight*rate*temp1
c deaggregation for the special case. not prepared oct 30 2007.
         enddo      !do k...
         goto 1282
        endif
      elseif(ipia.eq.12)then
      call getMota(ip,iper(ip),xmag2,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.22)then
      call getDahle95(iperdahl(ip),isoild,xmag2,rjb,gnd,sigmaf)
      endif
      wttmp=wt1wt2*weight*rate
        do ifn=1,nfi(ip)
        if(determ.and.isbig.and.isclose.and.norpt(ip,ia))then
       write(idet,679)exp(gnd(ifn)),1./sigmaf/sqrt2,ipia,ift,xmag2,rjb,rrup,
     + wttmp,depth0(ift),iftype(ift),ifn
      if(ifn .eq. nfi(ip))norpt(ip,ia)=.false.
      endif
        do  k=1,nlev(ip)
       pr=(gnd(ifn) - xlev(k,ip))*sigmaf
       if(pr.gt.3.3)then
       ipr=250002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 1283      !transfer out if ground motion above mu+3sigma
       endif
      cfac=wttmp*p(ipr)
        prob(icnt,k,ip,ifn,kg)= prob(icnt,k,ip,ifn,kg)+cfac
       if(deagg)then
c for individual fault
      frbar(ift,ifn,ip)=frbar(ift,ifn,ip)+cfac*rrup
      fmbar(ift,ifn,ip)=fmbar(ift,ifn,ip)+cfac*xmag
               eps= -pr*sqrt2
      febar(ift,ifn,ip)=febar(ift,ifn,ip)+cfac*eps
      fhaz(ift,ifn,ip) =fhaz(ift,ifn,ip)+cfac
c ifn index is present: 2nd to last frontier. ieps is  an explicit dimension, recomm. by Bazzuro
c Some attn. models will say a given M,R is a low eps0 combination, others will say a higher eps0.
c 
      if(eps.lt.emax)then
      ieps=max(1,min(int((eps+2.)*2.),10))
      rbar(ir,im,ieps,ifn,ip)=rbar(ir,im,ieps,ifn,ip)+cfac*rrup
      mbar(ir,im,ieps,ifn,ip)=mbar(ir,im,ieps,ifn,ip)+cfac*xmag2
      ebar(ir,im,ieps,ifn,ip)=ebar(ir,im,ieps,ifn,ip)+cfac*eps
      haz(ir,im,ieps,ifn,ip) =haz(ir,im,ieps,ifn,ip)+cfac
      if(eps.ge.2.)goto 1283
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
      ka=5
      temp=p(ipr)
      prx= temp-ptail(ka)
      dowhile(prx.gt.1.e-9)
      prob5(ir,im,ieps,ip,ifn,ka)=prob5(ir,im,ieps,ip,ifn,ka)+wttmp*prx
c  *prlr       this factor is not present now. Why?
      ka=ka-1
      if(ka.eq.0)goto 1283
      prx= temp-ptail(ka)
      enddo
      endif      !eps <emax
      endif      !if deagg
           enddo
1283      continue
      enddo      !nlev
1282    continue
      enddo        !ifn    
1284    continue
1285    continue
1403    continue
1406      continue      !added outer loop 8/2007, for goto if rate<tiny
      enddo      !imag new index 6/2006. Multiple mmax rate per fault.
c        write(6,*) totmo,(xmo2(ift,j),j=1,nmagf(ift))
c        read(5,*) idum
        go to 850
        endif
c---------------------
c------ all GR without uncertainties, with possible downdip ruptures, z=dtor1
      if((itype(ift).eq.2).and.(itest(ift).eq.1)) then
c-- loop through magnitudes
c        iftype2= iftype(ift)
        dtor = depth0(ift)      !depth to top of fault
c        print *,rx,ry
      do imag=1,nmagf(ift)      !can have two or more distributions often
c with diffrent Mmax.
        do 2285 m=1,nmag0(ift,imag)
        isbig= m .eq. nmag0(ift,imag)
        xmag= magmin(ift,imag)+(m-1)*dmag(ift,imag)
        rate= a(ift,imag) - b(ift,imag)*xmag
c relative rate
        rate= relwt(ift,imag) * 10.**rate
        if(rate.lt.tiny)goto 2285      !bypass irrelevant models aug 23 2007
      if(deagg)then
       im = nint((xmag-5.7)*10.)+1
       im=max(im, 1)
       im=min(im,30)
                  immax=max(immax,im)
       endif
          if(xmag.lt.mcut(1))then
          ime=1
          elseif(xmag.lt.mcut(2))then
          ime=2
          else
          ime=3
          endif
      if(nearsrc)then
      nrupdd=nrupd(m,ift)
      else
      nrupdd=1
      endif 
        rate= rate/float(nrups(m,ift)*nrupdd)
        kg=1      !no clustering in GR region june 2007
c--- down-dip loop through floating rupture zones
      do 2284 jrup=1,nrupdd
      norpt=.true.      !report deterministic only once vector assignment
        aspectratio=ar(m,jrup)
        dz=dlen2*float(jrup-1)
        dtor1=dtor + dz      !all important depth to top of rupture.
        w_rup=width(ift)-dz
        cyhwfac=atan(w_rup*0.5*cosDELTA/(dtor1+1.0))/(pi*0.5)
c variable factor above for CY relation when downdip rupture occurs        
c--- along-strike loop through floating rupture zones
        do 2284 irup=1,nrups(m,ift)
        rjb  = dmin2(m,irup,jrup,1)
c        print *,irup,rjb,dmin(1)
      xdif=abs(rjb-dmin(1))
c this complicated logic is trying to save data near fault endpoint properly.
c it is easy to stumble near fault endpoints in deterministic calcs. SHarmsen feb 27 2009
      isclose=xdif.lt.0.95.or.(rjb.gt.xdif.and.xdif.lt.2.5)
          if(rjb.lt.dcut(1))then
          ide=1
          elseif(rjb.lt.dcut(2))then
          ide=2
          else
          ide=3
          endif
      if(slist.and.rjb.lt.10000.)
     + write(30+i,399)rjb,xmag,rate,ift,dtor1,imag,irup
        do  ip=1,nper
        iq=iper(ip)      !std period index
        idet=ip+66      !possibly writing deterministic src data
        R_x = dmin2(m,irup,jrup,4)      !new R_x signed distance 10/07
        Ry0= dmin2(m,irup,jrup,5)      !type 5 new dec 2012.
        do 2282 ia=1,nattn(ip)
c         iftype2= iftype(ift)
         rrup = dmin2(m,irup,jrup,3)
        weight= wt(ip,ia,1)
        if(rjb.gt.wtdist(ip,ia)) weight=wt(ip,ia,2)
        if(rjb.gt.dmax) go to 2282
      if(deagg)then
        if(rrup.le.100.)then
        ir=0.1*rrup + 1
        elseif(rrup.le.300.)then
        ir=11+(rrup-100.)/50.
        elseif(rrup.le.500.)then
        ir=15
        else
        ir=16
        endif 
       endif      !if deagg
      ipia=iatten(ip,ia)
c nga relations are collected into the first set of models to check
      if(nga(ip,ia)) then
        if(ipia.eq.13) then
      call getBooreNGA07m(ip,iperb(ip),xmag,rjb,vs30,gnd,sigmaf)
         elseif(ipia.eq.14) then
      pgacalc=ip.eq.1.or.period(ip).eq.0.     !campbell pga calc 
        call getCBNGA1107 (ip,ipercb(ip),xmag,rrup,rjb,dtor1,
     1 vs30,dbasin,gnd,sigmaf)
      elseif(ipia.eq.15)then
c 10/2007 implementation of getCYNGA 
c        call CY2006H(ip,ipercy(ip), xmag, width(ift), rrup, rjb, vs30,  
c     1      dtor1, srcSiteA, F_RV, F_NM,gnd,sigmaf)
      call CY2007H(ip,ipercy(ip), xmag, rrup, rjb, R_x, vs30,z1, dtor1,
     1                 F_RV, F_NM,  gnd, sigmaf)
      elseif(ipia.eq.33)then
      indx1=indbssa(ip)
      mech=ibtype(ift)
      iregion=1	!CA-Taiwan version
      call bssa2013drv(ip,indx_pga,indx1,xmag,rjb,vs30,z1km,mech,iregion,alny,sigmaf)
c skip the 2nd period for now assume the period of interest was tabulated.
      gnd(1)=alny
         if(l_gnd_ep(ip))then      !additional epistemic uncert.
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
          if(ift.eq.-1)then
          print *,ip,y1,rjb,xmag,sigmaf
        print *,  mech,vs30,y1,expsiglny,' mech,vs30,y1,expsiglny'
      endif
      elseif(ipia.eq.37)then
      call  getIdriss2012(idriss(ip),ip,xmag,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.38)then
      call gksa12(ipgk(ip),ip,xmag,rrup,gnd,sigmaf,vs30,iftype(ift),dgkbasin,Q_CA)
      elseif(ipia.eq.39)then
      ib=0	!basin term Graizer & kalkan model
      call  gksa(ipgk(ip),ip,xmag,rrup,gnd,sigmaf,vs30,iftype(ift),IB)
c      if(ift.eq.1)print *,ip,exp(gnd(1)),rrup,xmag,1/sigmaf/sqrt2
      elseif(ipia.eq.35)then
c from BChiou email use fmeasured=1.
      F_Measured=1.
      F_Inferred=0.0	!could change these. Reference rock.
c      call  CY2012_NGA(ip,icy12(ip), xmag, rrup, rjb, R_x, vs30,
c     1                   F_Measured, F_Inferred, z1, dip0(ift), dtor1,
c     1                   F_RV, F_NM,
c     1                   psa_ref, gnd, tau, sigma_NL0, sigmaf)
       call CY2013_NGA(ip,icy13(ip), xmag, rrup, rjb, R_x, vs30,
     1                   F_Measured, F_Inferred, deltaZ1, dip0(ift), dtor1,
     1                   F_RV, F_NM, cDPP(ift), gnd, tau, sigma_NL0, sigmaf)

c      print *,exp(gnd(1)),ift,xmag,Rrup
      elseif(ipia.eq.36)then
c     Compute SA1100 for AS2012. New Mar 1 2013.
c hanging-wall flag for as12 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
      SA_rock = 0.
      call AS2013_v10_model ( ipas13(ip),xmag,dip0(ift), W_rup, dtor1, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      Sa1100 = exp(lnSa)
c 
c     Compute Sa at spectral period for given Vs30
      call AS2013_v10_model ( ipas13(ip),xmag, dip0(ift), W_rup, dtor1, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30, SA1100, z1km, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
c      print *,exp(lnSa),phi,tau,w_rup
      gnd(1)=lnSa
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
        sigmaf = 1./sqrt( phi**2 + tau**2 )/sqrt2
      elseif(ipia.eq.34)then
      SJ=0.0	!Are we in Japan? 1 = yes.
      Hhyp=min(12.,dtor1+8.)	!Take a stab at hypocenter depth.
      call CB13_NGA_SPEC  (ip,ipcb12(ip), xmag,Rrup,R_x,Rjb,F_RV, F_NM,dtor1,Hhyp,w_rup,dbasin,Vs30,
     +Dip0(ift),SJ,gnd,sigmaf)
c      print *,exp(gnd(1)),ift,xmag,Rrup
        elseif (ipia.eq.16) then
c hanging-wall flag for as08 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
c Although iflag makes a token appearance in SR code, it isn't used. So I dropped iflag (SH).
       call AS_072007 ( ip,ipera(ip),xmag, dip0(ift), F_NM,F_RV, W_rup, rRup, rjb,R_x,
     1                     vs30, hwflag, gnd(1), sigma1, dtor1,  vs30_class,
     3                     z1km )
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
c gnd is lnSA in AS code. need to fill in some other input. esp. aspectratio.
c period1 is output by AS_modelXX. period1 should equal per(ip)
      sigmaf=1./sqrt2/sigma1
      elseif(ipia.eq.17)then
      call getIdriss(ip,iper(ip),xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.18)then
c new nov 7,2006. Kanno
      call kanno(ip,iperk(ip),gnd,sigmaf,xmag,rrup,vsfac)
      endif
      elseif(wus02(ip,ia))then
c WUS pre-nga (2002 atten models)
      if(ipia.eq.1)then
      call getSpu2000(ip,iq,xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.3)then
      call getsadigh(ip,iq,xmag,rrup,gnd,sigmaf)
      elseif(ipia.eq.8)then
      call getAS97(ip,iq,xmag,rrup,hwsite,gnd,sigmaf)
       elseif(ipia.eq.9)then
       call getCamp2003(ip,iq,xmag,rrup,rjb,gnd,sigmaf)
      elseif(ipia.eq.11)then
      call getBJF97(ip,iq,xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.-3)then
        call SomerIMW(ip,isomer(ip), xmag,  rrup, rjb,vs30,gnd,sigmaf)
      endif
      elseif(ceus02(ip,ia).or.ceus11(ip,ia))then
c
c CEUS pre-nga (2002 atten models). These motions can be "clamped"
c
       elseif(ceus02(ip,ia).or.ceus11(ip,ia))then
c added 3 new CENA models defined by tables. Mar 2011
       rjbp=max(rjb,0.11)
       if(ipia.eq.25)then
      rkm=rjbp
      jf=ia08(ip)
       ka=1
            sigma = 0.3*2.303
      gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
       sigmaf = 1./sqrt2/sigma
      elseif(ipia.eq.26)then
       rkm=rrup
       jf=ia06(ip)
       ka=2
             sigma = 0.3*2.303
       gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
       sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.27) then
       rkm=rrup
       jf=ip11(ip)
      sigma = sigPez11(xmag,jf)
       ka=3
       gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
      sigmaf = 1./sqrt2/sigma
      elseif(ipia.eq.2)then
      call getToro(ip,iq,1,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-2)then
      call getToro(ip,iq,2,xmag,rjb,gnd,sigma,sigmaf)
        elseif(ipia.eq.4)then
        call getAB06(ip,iperab(ip),irab(ip,ia),xmag,rrup,gnd,sigma,sigmaf,vs30)
        elseif(ipia.eq.19)then
      call getTP05(ip,ipertp(ip),irtb,xmag,rrup,gnd,sigma,sigmaf)
        elseif(ipia.eq.20)then
      call getSilva(ip,isilva(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
         elseif(ipia.eq.28)then
       call getSilvaV(ip,isilva2(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.6)then
      call getFEA(ip,iq,1,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-6)then      !HR
      call getFEA(ip,iq,3,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.5)then      !Rto distinguish ab95 from fea
      call getFEA(ip,iq,2,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-5)then      !HR for atkinson-boore
      call getFEA(ip,iq,4,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.7) then
        call getSomer(ip,iq,1,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-7) then
        call getSomer(ip,iq,2,xmag,rjb,gnd,sigma,sigmaf)
             elseif(ipia.eq.10)then
        call getCampCEUS(ip,iq,1,xmag,rrup,gnd,sigma,sigmaf)
             elseif(ipia.eq.-10)then
        call getCampCEUS(ip,iq,2,xmag,rrup,gnd,sigma,sigmaf)
      endif
c This is a good place to apply the extra clamp constraint which
c in 2002 was used only for CEUS sources. The precomputed p() matrix does
c not work for this case because truncation is no longer at mu+3sig
      test0=gnd(1)+3.*sigma
      test= exp(test0)
      if(clamp(iq).lt.test .and. clamp(iq).gt.0.) then
      clamp2= alog(clamp(iq))
      sigmasq=1./sigma/sqrt2
      tempgt3= (gnd(1) - clamp2)*sigmasq
      probgt3= (erf(tempgt3)+1.)*0.5
      prr=1./(1.-probgt3)
        if(determ.and.isbig.and.isclose.and.norpt(ip,ia))then
c write deterministic median if it's big and close and you havent written but should
      write(idet,679)exp(gnd(ifn)),1./sigmaf/sqrt2,ipia,
     + ift,xmag,rjb,rrup,wttmp,
     + dtor1,iftype(ift),ifn
c depth to top of rupture = dtor1. can have report for different dtors, same 
c fault segment.
       norpt(ip,ia)=.false.      !there is now a report
       endif
      do  k=1,nlev(ip)
      temp= (gnd(1) - xlev(k,ip))*sigmasq
      temp1= (erf(temp)+1.)*0.5
      temp1= (temp1-probgt3)*prr
      if(temp1.lt.0.) goto 2282      !no more calcs once p<0
        prob(icnt,k,ip,1,kg)= prob(icnt,k,ip,1,kg)+ 
     &         weight*rate*temp1
        enddo      !do k=1,nlev
        goto 2282
        endif
      elseif(ipia.eq.12)then
      call getMota(ip,iper(ip),xmag,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.22)then
      call getDahle95(iperdahl(ip),isoild,xmag,rjb,gnd,sigmaf)
      endif
      wttmp=weight*rate
        do ifn=1,nfi(ip)
        if(determ.and.isbig.and.isclose.and.norpt(ip,ia))then
c write home to mom if it's big and close and you havent written but should
      write(idet,679)exp(gnd(ifn)),1./sigmaf/sqrt2,ipia,
     + ift,xmag,rjb,rrup,wttmp,
     + dtor1,iftype(ift),ifn
c depth to top of rupture = dtor1. can have report for different dtors, same 
c fault segment.
       if(ifn .eq. nfi(ip)) norpt(ip,ia)=.false.      !there is now a report
       endif
        do  k=1,nlev(ip)
       pr=(gnd(ifn) - xlev(k,ip))*sigmaf
       if(pr.gt.3.3)then
       ipr=250002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 2283      !transfer out if ground motion above mu+3sigma
       endif
      cfac =  wttmp*p(ipr)
        prob(icnt,k,ip,ifn,kg)= prob(icnt,k,ip,ifn,kg)+cfac 
               if(deagg)then
c for individual fault
      frbar(ift,ifn,ip)=frbar(ift,ifn,ip)+cfac*rrup
      fmbar(ift,ifn,ip)=fmbar(ift,ifn,ip)+cfac*xmag
               eps= -pr*sqrt2
      febar(ift,ifn,ip)=febar(ift,ifn,ip)+cfac*eps
      fhaz(ift,ifn,ip)=fhaz(ift,ifn,ip)+cfac
c ifn index is present: 2nd to last frontier. ieps is  an explicit dimension, recomm. by Bazzuro
c Some attn. models will say a given M,R is a low eps0 combination, others will say a higher eps0.
c store separate ifn in different records. Why? to give 'em diff. weights when combining.
      if(eps.lt.emax)then
      ieps=max(1,min(int((eps+2.)*2.),10))
      rbar(ir,im,ieps,ifn,ip)=rbar(ir,im,ieps,ifn,ip)+cfac*rrup
      mbar(ir,im,ieps,ifn,ip)=mbar(ir,im,ieps,ifn,ip)+cfac*xmag
      ebar(ir,im,ieps,ifn,ip)=ebar(ir,im,ieps,ifn,ip)+cfac*eps
      haz(ir,im,ieps,ifn,ip)=haz(ir,im,ieps,ifn,ip)+cfac
      if(eps.ge.2.)goto 2283
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
      ka=5
      temp=p(ipr)
      prx= temp-ptail(ka)
      dowhile(prx.gt.1.e-9)
      prob5(ir,im,ieps,ip,ifn,ka)=prob5(ir,im,ieps,ip,ifn,ka)+wttmp*prx
c  *prlr       this factor is not present now. Why?
      ka=ka-1
      if(ka.eq.0)goto 2283
      prx= temp-ptail(ka)
      enddo
      endif      !eps <emax
      endif      !if deagg
      enddo      !levels
2283      continue
      enddo      !files
2282    continue
      enddo
2284    continue
2285    continue
      enddo      !imag loop, different GR distributions.
        go to 850
        endif
c---------------------------
c--- for characteristic event, modified to include multiple mags per fault
c---- characteristic with uncertainties
c To some extent diff M(RA) uncertainties take care of mag. variation. We seem to
c be repeating some mag uncertainty here.
c characteristic events fill the fault by definition. no variation in top of rup.
c discuss with colleagues? nov 15 2006
      if((itype(ift).eq.1).and.(itest(ift).eq.0)) then
      do imag=1,nmagf(ift)      !new 6/06 consolidate mag variation for each fault
        xmag= cmag(ift,imag)
        xmo= 1.5*xmag +9.05
        xmo= 10.**xmo
c relative rate nov 2006. SH. But not for cluster model. May 2007.
      if(cluster)then
      rate= 1      !weights will be treated differently, using wtscene(*,*)
      else
        rate= relwt(ift,imag) * crate(ift,imag)
        endif
        if(rate.lt.tiny)goto 4203      !bypass irrelevant models aug 23 2007
        xmorate= xmo *rate
c cluster model new may 16 2007. Base cluster on mag and geographic branch for Char. For GR
c no cluster model has been specifically discussed
      if(cluster)then
      ks=imag
      kg=igroup(ift)
      else
      ks=1
      kg=1
      endif
        do 203 ilt=1,nbranch
        xmag= cmag(ift,imag) + dmbranch(ilt)
        xmo= 1.5*xmag+9.05
        xmo= 10.**xmo
        rate= xmorate/xmo
        totmo =0.
        do 244 idis= -mwid, mwid
        xmag2= xmag + idis*dma
          if(xmag2.lt.mcut(1))then
          ime=1
          elseif(xmag2.lt.mcut(2))then
          ime=2
          else
          ime=3
          endif
        rjb = dmin(1)
          if(rjb.lt.dcut(1))then
          ide=1
          elseif(rjb.lt.dcut(2))then
          ide=2
          else
          ide=3
          endif
        wt1wt2=wtbranch(ilt)*ale2(idis+mwid+1)
       if(deagg)then
       im = nint((xmag2-5.7)*10.)+1
       im=max(1,im)
       im=min(30,im)
                  immax=max(immax,im)
       endif
      if(slist)write(30+i,399)rjb,xmag2,rate*wt1wt2,ift,depth0(ift),ilt,idis
        do  ip=1,nper
        iq = iper(ip)
        idet=ip+66
        R_x = dmin(4)      !new signed distance 10/07
        do 204 ia=1,nattn(ip)
c        iftype2= iftype(ift)
        rrup = dmin(3)      !type3 r_cd or rrup
        Ry0= dmin(5)      !type 5 new dec 2012.
c        write(6,*)rjb,rrup,ip,' rjb rrup ip'
        if(rjb.gt.dmax)goto 204      !sail on out if too far.
      ipia=iatten(ip,ia)
c       if(imag.eq.1.and.ia.eq.1.and.ip.eq.1.and.ilt.eq.1)write(6,*)rx,ry,R_x,rjb,ift,crate(ift,imag)
        weight= wt(ip,ia,1)
        if(rjb.gt.wtdist(ip,ia)) weight=wt(ip,ia,2)
      if(deagg)then
        if(rrup.le.100.)then
        ir=0.1*rrup + 1
        elseif(rrup.le.300.)then
        ir=11+(rrup-100.)/50.
        elseif(rrup.le.500.)then
        ir=15
        else
        ir=16
        endif 
       endif      !if deagg
c nga relations are collected into the first set of models to check
      if(nga(ip,ia)) then
        if(ipia.eq.13) then
      call getBooreNGA07m(ip,iperb(ip),xmag2,rjb,vs30,gnd,sigmaf)
c      if(xmag2.ge.8.)print *,rjb,exp(gnd(1)),xmag2,1./sigmaf/sqrt2
         elseif(ipia.eq.14) then
      pgacalc=ip.eq.1.or.period(ip).eq.0.     !campbell pga calc d
        call getCBNGA1107 (ip,ipercb(ip),xmag2,rrup,rjb,dtor,
     1 vs30,dbasin,gnd,sigmaf)
      elseif(ipia.eq.15)then
c 3/2008 implementation of getCYNGA 
c        call CY2006H(ip,ipercy(ip), xmag2, width(ift), rrup, rjb, vs30,  
c     1      dtor, srcSiteA, F_RV, F_NM,gnd,sigmaf)
      call CY2007H(ip,ipercy(ip), xmag2, rrup, rjb, R_x, vs30,z1, dtor,
     1                 F_RV, F_NM,  gnd, sigmaf)
      elseif(ipia.eq.33)then
      indx1=indbssa(ip)
      mech=ibtype(ift)
      iregion=1	!CA-Taiwan version
      call bssa2013drv(ip,indx_pga,indx1,xmag2,rjb,vs30,z1km,mech,iregion,alny,sigmaf)
c skip the 2nd period for now assume the period of interest was tabulated.
      gnd(1)=alny
c      if(xmag2.ge.8)print *,rjb,exp(alny),xmag2,1./sigmaf/sqrt2
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
c          if(ift.eq.-1)then
c          if(ip.eq.2)print *,rjb,xmag2,exp(alny),1./sigmaf/sqrt2
c          print *, per1,  rjb, r1,' per1, rjb, r1'
c      endif
      elseif(ipia.eq.37)then
      call  getIdriss2012(idriss(ip),ip,xmag2,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.38)then
      call gksa12(ipgk(ip),ip,xmag2,rrup,gnd,sigmaf,vs30,iftype(ift),dgkbasin,Q_CA)
c      print *,xmag2,rrup,gnd(1),sigmaf
      elseif(ipia.eq.39)then
      ib=0	!basin term Graizer & kalkan model
      call  gksa(ipgk(ip),ip,xmag2,rrup,gnd,sigmaf,vs30,iftype(ift),IB)
c      if(ift.eq.1)print *,ip,exp(gnd(1)),rrup,xmag,1/sigmaf/sqrt2
      elseif(ipia.eq.35)then
c from BChiou email use fmeasured=1.
      F_Measured=1.
      F_Inferred=0.0	!could change these. Reference rock.
c      call  CY2012_NGA(ip,icy12(ip), xmag2, rrup, rjb, R_x, vs30,
c     1                   F_Measured, F_Inferred, z1, dip0(ift), dtor,
c     1                   F_RV, F_NM,
c     1                   psa_ref, gnd, tau, sigma_NL0, sigmaf)
       call CY2013_NGA(ip,icy13(ip), xmag2, rrup, rjb, R_x, vs30,
     1                   F_Measured, F_Inferred, deltaZ1, dip0(ift), dtor,
     1                   F_RV, F_NM, cDPP(ift), gnd, tau, sigma_NL0, sigmaf)

c      print *,exp(gnd(1)),ift,xmag2,Rrup
      elseif(ipia.eq.36)then
c     Compute SA1100 for AS2012. New Mar 1 2013.
      SA_rock = 0.
c hanging-wall flag for as12 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
      dippy=max(dip0(ift),65.)
      Ry0=dmin(5)      !new dec 10
      call AS2013_v10_model ( ipas13(ip),xmag2,dip0(ift), Width(ift), dtor, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      Sa1100 = exp(lnSa)
c 
c     Compute Sa at spectral period for given Vs30
      call AS2013_v10_model ( ipas13(ip),xmag2, dip0(ift), Width(ift), dtor, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30, SA1100, z1km, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      gnd(1)=lnSa
      sig=sqrt( phi**2 + tau**2 )
c      print *,exp(lnSa),xmag2,rrup,rjb,r_x,Sa1100,sig,period(ip)
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
        sigmaf = 1./sig/sqrt2
      elseif(ipia.eq.34)then
      SJ=0.0	!Are we in Japan? 1 = yes.
      Hhyp=min(12.,dtor+8.)	!Take a stab at hypocenter depth.
      call CB13_NGA_SPEC  (ip,ipcb12(ip), xmag2,Rrup,R_x,Rjb,F_RV, F_NM,dtor,Hhyp,width(ift),dbasin,Vs30,
     +Dip0(ift),SJ,gnd,sigmaf)
c      print *,exp(gnd(1)),ift,xmag2,Rrup
        elseif (ipia.eq.16) then
c hanging-wall flag for as08 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
c Although iflag makes a token appearance in SR code, it isn't used. So I dropped iflag (SH).
c the following call should include Rx rather than trying  to compute Rx.
       call AS_072007 ( ip,ipera(ip),xmag2, dip0(ift), F_NM,F_RV, Width(ift), rRup, rjb,R_x,
     1                     vs30, hwflag, gnd(1), sigma1, dtor,  vs30_class,
     3                     z1km )
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
c gnd is lnSA in AS code. need to fill in some other input. esp. aspectratio.
c period1 is output by AS_modelXX. period1 should equal per(ip)
      sigmaf=1./sqrt2/sigma1
      elseif(ipia.eq.17)then
      call getIdriss(ip,iper(ip),xmag2,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.18)then
c new nov 7,2006. Kanno
      call kanno(ip,iperk(ip),gnd,sigmaf,xmag2,rrup,vsfac)
      endif
      elseif(wus02(ip,ia))then
c WUS pre-nga (2002 atten models)
      if(ipia.eq.1)then
      call getSpu2000(ip,iq,xmag2,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.3)then
      call getsadigh(ip,iq,xmag2,rrup,gnd,sigmaf)
      elseif(ipia.eq.8)then
      call getAS97(ip,iq,xmag2,rrup,hwsite,gnd,sigmaf)
       elseif(ipia.eq.9)then
       call getCamp2003(ip,iq,xmag2,rrup,rjb,gnd,sigmaf)
      elseif(ipia.eq.11)then
      call getBJF97(ip,iq,xmag2,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.-3)then
        call SomerIMW(ip,isomer(ip), xmag2,  rrup, rjb,vs30,gnd,sigmaf)
      endif
      elseif(ceus02(ip,ia).or.ceus11(ip,ia))then
c CEUS pre-nga (2002 atten models)
       rjbp=max(rjb,0.11)
       if(ipia.eq.25)then
             rkm=rjbp
             jf=ia08(ip)
             ka=1
             sigma = 0.3*2.303
             gnd(1) = amean11(xmag2,rkm,rjbp,vs30,ip,jf,ka)
              if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
            sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.26)then
             rkm=rrup
             jf=ia06(ip)
             ka=2
               sigma = 0.3*2.303
            gnd(1) = amean11(xmag2,rkm,rjbp,vs30,ip,jf,ka)
           if(lceus_sigma)sigma=ceus_sigma     !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
            sigmaf = 1./sqrt2/sigma
       elseif(ipia.eq.27)then
            rkm=rrup
             jf=ip11(ip)
      sigma = sigPez11(xmag2,jf)
            ka=3
             gnd(1) = amean11(xmag2,rkm,rjbp,vs30,ip,jf,ka)
               if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
c     print *,gnd,rkm,ip,jf,ka
             sigmaf = 1./sqrt2/sigma
      elseif(ipia.eq.2)then
c        write(6,*)'iq clamp(iq) gnd rrup ip Toro ',iq,clamp(iq),gnd,rrup,ip
      call getToro(ip,iq,1,xmag2,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-2)then
      call getToro(ip,iq,2,xmag2,rjb,gnd,sigma,sigmaf)
        elseif(ipia.eq.4)then
        call getAB06(ip,iperab(ip),irab(ip,ia),xmag2,rrup,gnd,sigma,sigmaf,vs30)
        elseif(ipia.eq.19)then
      call getTP05(ip,ipertp(ip),irtb,xmag2,rrup,gnd,sigma,sigmaf)
        elseif(ipia.eq.20)then
      call getSilva(ip,isilva(ip),irsilva,xmag2,rjb,gnd,sigma,sigmaf)
         elseif(ipia.eq.28)then
      call getSilvaV(ip,isilva2(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.6)then
      call getFEA(ip,iq,1,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-6)then      !FEA- HR
      call getFEA(ip,iq,3,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.5)then      !to distinguish ab95 from fea
c hard rock possibilities. These are  available for CEUS atten models as
c of July or August 2006.
      call getFEA(ip,iq,2,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-5)then      !HR for atkinson-boore
      call getFEA(ip,iq,4,xmag2,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.7) then
        call getSomer(ip,iq,1,xmag2,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-7) then
        call getSomer(ip,iq,2,xmag2,rjb,gnd,sigma,sigmaf)
             elseif(ipia.eq.10)then
        call getCampCEUS(ip,iq,1,xmag2,rrup,gnd,sigma,sigmaf)
             elseif(ipia.eq.-10)then
        call getCampCEUS(ip,iq,2,xmag2,rrup,gnd,sigma,sigmaf)
      endif
c This is a good place to apply the extra clamp constraint which
c in 2002 was used only for CEUS sources. The precomputed p() matrix does
c not work for this case because truncation is no longer at mu+3sig
      test0=gnd(1)+3.*sigma
      test= exp(test0)
      if(clamp(iq).lt.test .and. clamp(iq).gt.0.) then
      clamp2= alog(clamp(iq))
      sigmasq=1./sigma/sqrt2
      tempgt3= (gnd(1) - clamp2)*sigmasq
      probgt3= (erf(tempgt3)+1.)*0.5
      prr=1./(1.-probgt3)
c in the CEUS, only one epistemic branch on gnd. This could change.
c CEUS deagg not attmpted Oct 30 2007.
c Deterministic output has to occur whether special clamping is in effect or not...
        if(determ)then
        isbig=.false.
        isclose=.false.
        if(rjb.le.rjbdet)then
        isclose=.true.
        rjbdet=rjb
        endif
        if(idis.ge.0.and.xmag2.gt.xmagdet)then
        xmagdet=xmag2
        isbig=.true.
        endif
c mx_index is  the index of the maximum magnitude in the set 1,...,nmagf(ift)
       if(mx_index(ift).eq.imag.and.(isbig.or.isclose))
     + write(idet,679)exp(gnd(1)),rrup,xmag,1./sigmaf/sqrt2,ipia,ift,xmag2,rjb,
     + rrup,wttmp,dtor1,iftype(ift),1
      endif      !if determ
      do  k=1,nlev(ip)
      temp= (gnd(1) - xlev(k,ip))*sigmasq
      temp1= (erf(temp)+1.)*0.5
      temp1= (temp1-probgt3)*prr
      if(temp1.lt.0.) goto 204      !no more calcs once p<0
      if(cluster)then
            prob_s(jseg(ift),ks,kg,k,ip)=prob_s(jseg(ift),ks,kg,k,ip)+ wt1wt2*weight*rate*temp1
      else
        prob(icnt,k,ip,1,kg)= prob(icnt,k,ip,1,kg)+ wt1wt2*weight*rate*temp1
        endif
         enddo      !do k...
         goto 204
        endif
      elseif(ipia.eq.12)then
      call getMota(ip,iper(ip),xmag2,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.22)then
      call getDahle95(iperdahl(ip),isoild,xmag2,rjb,gnd,sigmaf)
      endif
      wttmp=wt1wt2*weight*rate
        do ifn = 1,nfi(ip)      !can have 3 files one per epistemic branch
        if(determ)then
        isbig=.false.
        isclose=.false.
        if(rjb.le.rjbdet)then
        isclose=.true.
        rjbdet=rjb
        endif
        if(idis.ge.0.and.xmag2.gt.xmagdet)then
        xmagdet=xmag2
        isbig=.true.
        endif
c mx_index is  the index of the maximum magnitude in the set 1,...,nmagf(ift)
       if(mx_index(ift).eq.imag.and.(isbig.or.isclose))
     + write(idet,679)exp(gnd(ifn)),1./sigmaf/sqrt2,ipia,ift,xmag2,rjb,
     + rrup,wttmp,dtor1,iftype(ift),ifn
      endif      !if determ
        do  k=1,nlev(ip)
       pr=(gnd(ifn) - xlev(k,ip))*sigmaf
       if(pr.gt.3.3)then
       ipr=250002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 206      !transfer out if ground motion above mu+3sigma
       endif
      if(cluster)then
            prob_s(jseg(ift),ks,kg,k,ip)=prob_s(jseg(ift),ks,kg,k,ip)+ wt1wt2*weight*rate*p(ipr)
      else
       cfac= wttmp*p(ipr)
        prob(icnt,k,ip,ifn,kg)= prob(icnt,k,ip,ifn,kg)+ cfac
               if(deagg)then
c deagg not ready for clustered events oct 2007
c some of the standard arrays have been added but this is not enough. 
c for individual fault
      frbar(ift,ifn,ip)=frbar(ift,ifn,ip)+cfac*rrup
      fmbar(ift,ifn,ip)=fmbar(ift,ifn,ip)+cfac*xmag
               eps= -pr*sqrt2
      febar(ift,ifn,ip)=febar(ift,ifn,ip)+cfac*eps
      fhaz(ift,ifn,ip)=fhaz(ift,ifn,ip)+cfac
c
c ieps is  an explicit dimension, recomm. by P. Bazzuro
c Some attn. models will say a given M,R is a low eps0 combination, others will say a higher eps0.
c 
      if(eps.lt.emax)then
      ieps=max(1,min(int((eps+2.)*2.),10))
      rbar(ir,im,ieps,ifn,ip)=rbar(ir,im,ieps,ifn,ip)+cfac*rrup
      mbar(ir,im,ieps,ifn,ip)=mbar(ir,im,ieps,ifn,ip)+cfac*xmag2
      ebar(ir,im,ieps,ifn,ip)=ebar(ir,im,ieps,ifn,ip)+cfac*eps
      haz(ir,im,ieps,ifn,ip)=haz(ir,im,ieps,ifn,ip)+cfac
      if(eps.ge.2.)goto 206
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
      ka=5
      temp=p(ipr)
      prx= temp-ptail(ka)
      dowhile(prx.gt.1.e-9)
      prob5(ir,im,ieps,ip,ifn,ka)=prob5(ir,im,ieps,ip,ifn,ka)+wttmp*prx
c  *prlr       this factor is not present now. Why?
      ka=ka-1
      if(ka.eq.0)goto 206
      prx= temp-ptail(ka)
      enddo
      endif      !eps <emax
      endif      !if deagg
        endif
        enddo      !k=1,nlev
206      continue
      enddo       !ifn files
 
ccccccccc
 204    continue
      enddo      
 244    continue
c        write(6,*) "totmo, xmorate",ift,ilt,totmo, xmorate
c        read(5,*) idum
 203    continue
4203      continue
      enddo      !imag index
        endif
ccccccccccccccccc
cc-- characteristic without uncertainties
      if((itype(ift).eq.1).and.(itest(ift).eq.1)) then
      do imag=1,nmagf(ift)      !new 6/06 consolidate mag variation for given fault
        xmag= cmag(ift,imag)
        xmag0=xmag
c cluster model new may 16 2007. Base cluster on mag and geographic branch for Char. For GR
c no cluster model has been specifically discussed
      if(cluster)then
      ks=imag
      kg=igroup(ift)
      rate= 1      !weights will be treated differently, using wtscene(*,*)
      else
c relative rate nov 2006. SH. But not for cluster model. May 2007.
        rate= relwt(ift,imag) * crate(ift,imag)
      ks=1
      kg=1
      endif
      if(rate.lt.tiny)goto 3303
      if(deagg)then
      im=nint((xmag-5.7)*10.) +1
      im=max(im,1)
      im=min(im,30)
                  immax=max(immax,im)
      endif
        rjb = dmin(1)
c        iftype2= iftype(ift)
      if(l_gnd_ep(1))then
          if(xmag.lt.mcut(1))then
          ime=1
          elseif(xmag.lt.mcut(2))then
          ime=2
          else
          ime=3
          endif
          if(rjb.lt.dcut(1))then
          ide=1
          elseif(rjb.lt.dcut(2))then
          ide=2
          else
          ide=3
          endif
          endif       !(gnd_ep calcs to be performed)
      if(slist)write(30+i,399)rjb,xmag,rate,ift,depth0(ift),imag
        do  ip=1,nper
      idet=ip+66
      R_x = dmin(4)      !new signed distance, 10/2007
        Ry0= dmin(5)      !type 5 new dec 2012.
        do 3203 ia=1,nattn(ip)
        rrup= dmin(3)
        if(rjb.gt.dmax)goto 3203
c       if(imag.eq.1.and.ia.eq.1.and.ip.eq.1)write(6,*)rx,ry,R_x,rjb,ift,crate(ift,imag)
      if(deagg)then
        if(rrup.le.100.)then
        ir=0.1*rrup + 1
        elseif(rrup.le.300.)then
        ir=11+(rrup-100.)/50.
        elseif(rrup.le.500.)then
        ir=15
        else
        ir=16
        endif 
       endif      !if deagg
        iq=iper(ip)
      ipia=iatten(ip,ia)
        weight= wt(ip,ia,1)
        if(rjb.gt.wtdist(ip,ia)) weight=wt(ip,ia,2)
c nga relations are collected into the first set of models to check
      if(nga(ip,ia)) then
        if(ipia.eq.13) then
      call getBooreNGA07m(ip,iperb(ip),xmag,rjb,vs30,gnd,sigmaf)
         elseif(ipia.eq.14) then
      pgacalc=ip.eq.1
        call getCBNGA1107 (ip,ipercb(ip),xmag,rrup,rjb,dtor,
     1 vs30,dbasin,gnd,sigmaf)
      elseif(ipia.eq.15)then
c 10/2007 implementation of getCYNGA 
c        call CY2006H(ip,ipercy(ip), xmag, width(ift), rrup, rjb, vs30,  
c     1      dtor, srcSiteA, F_RV, F_NM,gnd,sigmaf)
      call CY2007H(ip,ipercy(ip), xmag, rrup, rjb, R_x, vs30,z1, dtor,
     1                 F_RV, F_NM,  gnd, sigmaf)
      elseif(ipia.eq.33)then
      indx1=indbssa(ip)
      mech=ibtype(ift)
      iregion=1	!CA-Taiwan version
      call bssa2013drv(ip,indx_pga,indx1,xmag,rjb,vs30,z1km,mech,iregion,alny,sigmaf)
c skip the 2nd period for now assume the period of interest was tabulated.
      gnd(1)=alny
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
          if(ift.eq.-1)then
          print *,ip,alny,rjb,xmag,sigmaf
          print *, per1, xmag, rjb, r1,' per1, xmag, rjb, r1'
      endif
      elseif(ipia.eq.37)then
      call  getIdriss2012(idriss(ip),ip,xmag,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.38)then
      call gksa12(ipgk(ip),ip,xmag,rrup,gnd,sigmaf,vs30,iftype(ift),dgkbasin,Q_CA)
c      print *,xmag,rrup,gnd(1),sigmaf
      elseif(ipia.eq.39)then
      ib=0	!basin term Graizer & kalkan model
      call  gksa(ipgk(ip),ip,xmag,rrup,gnd,sigmaf,vs30,iftype(ift),IB)
c      if(ift.eq.1)print *,ip,exp(gnd(1)),rrup,xmag,1/sigmaf/sqrt2
      elseif(ipia.eq.35)then
c from BChiou email use fmeasured=1.
      F_Measured=1.
      F_Inferred=0.0	!could change these. Reference rock.
c      call  CY2012_NGA(ip,icy12(ip), xmag, rrup, rjb, R_x, vs30,
c     1                   F_Measured, F_Inferred, z1, dip0(ift), dtor,
c     1                   F_RV, F_NM,
c     1                   psa_ref, gnd, tau, sigma_NL0, sigmaf)
       call CY2013_NGA(ip,icy13(ip), xmag, rrup, rjb, R_x, vs30,
     1                   F_Measured, F_Inferred, deltaZ1, dip0(ift), dtor,
     1                   F_RV, F_NM, cDPP(ift), gnd, tau, sigma_NL0, sigmaf)


c      print *,gnd(1),ift
      elseif(ipia.eq.36)then
c     Compute SA1100 for AS2012. New Mar 1 2013.
      SA_rock = 0.
c hanging-wall flag for as12 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
      call AS2013_v10_model ( ipas13(ip),xmag,dip0(ift), Width(ift), dtor, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30_rock, SA_rock, Z10_rock, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      Sa1100 = exp(lnSa)
c      print *,rx,ry,xmag,rRup,r_x,Sa1100,period(ip),' AS12 '
c 
c     Compute Sa at spectral period for given Vs30
      call AS2013_v10_model ( ipas13(ip),xmag, dip0(ift), Width(ift), dtor, F_rv, F_NM, rRup, rjb, r_x, Ry0, 
     1                     vs30, SA1100, z1km, z1_ref, hwflag, vs30_class,lnSa, phi, tau, useRy0(ip) )
      gnd(1)=lnSa
      sig=sqrt( phi**2 + tau**2 )
      print *,exp(lnSa),sig,Sa1100,period(ip),xmag,rRup
               if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
        sigmaf = 1./sig/sqrt2
      elseif(ipia.eq.34)then
      SJ=0.0	!Are we in Japan? 1 = yes.
      Hhyp=min(12.,dtor+8.)	!Take a stab at hypocenter depth.
      call CB13_NGA_SPEC  (ip,ipcb12(ip), xmag,Rrup,R_x,Rjb,F_RV, F_NM,dtor,Hhyp,width(ift),dbasin,Vs30,
     +Dip0(ift),SJ,gnd,sigmaf)
c      print *,exp(gnd(1)),ift,xmag,Rrup,Period(ip)
        elseif (ipia.eq.16) then
c hanging-wall flag for as08 model added dec 7 2012.
      if(dip0(ift).eq.90.)then
      hwflag=0
      elseif(rjb.le.0.01.and.r_x.ge.0.)then
      hwflag=1
      else
      hwflag=0
      endif
c Although iflag makes a token appearance in SR code, it isn't used. So I dropped iflag (SH).
       call AS_072007 ( ip,ipera(ip),xmag, dip0(ift), F_NM,F_RV, Width(ift), rRup, rjb,R_x,
     1                     vs30, hwflag, gnd(1), sigma1, dtor,  vs30_class,
     3                     z1km )
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
c gnd is lnSA in AS code. need to fill in some other input. vs30_class.
c 
      sigmaf=1./sqrt2/sigma1
      elseif(ipia.eq.17)then
      call getIdriss(ip,iper(ip),xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.18)then
c new nov 7,2006. Kanno
      call kanno(ip,iperk(ip),gnd,sigmaf,xmag,rrup,vsfac)
      endif
      elseif(wus02(ip,ia))then
c WUS pre-nga (2002 atten models)
      if(ipia.eq.1)then
      call getSpu2000(ip,iq,xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.3)then
      call getsadigh(ip,iq,xmag,rrup,gnd,sigmaf)
      elseif(ipia.eq.8)then
      call getAS97(ip,iq,xmag,rrup,hwsite,gnd,sigmaf)
       elseif(ipia.eq.9)then
       call getCamp2003(ip,iq,xmag,rrup,rjb,gnd,sigmaf)
      elseif(ipia.eq.11)then
      call getBJF97(ip,iq,xmag,rjb,vs30,gnd,sigmaf)
      elseif(ipia.eq.-3)then
        call SomerIMW(ip,isomer(ip), xmag,  rrup, rjb,vs30,gnd,sigmaf)
      endif
      elseif(ceus02(ip,ia).or. ceus11(ip,ia))then
      ifn=1
c CEUS pre-nga (2002 atten models)
c
c CEUS pre-nga (2002 atten models) and Gail Atkinson's latest 2011 models
c added 3 new CENA models defined by tables. Mar 2011
       rjbp=max(rjb,0.11)
       if(ipia.eq.25)then
      rkm=rjbp
      jf=ia08(ip)
       ka=1
            sigma = 0.3*2.303
      gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
       if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
       sigmaf = 1./sqrt2/sigma
      elseif(ipia.eq.26)then
      rkm=rrup
      jf=ia06(ip)
      ka=2
            sigma = 0.3*2.303
      gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
      if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
      sigmaf = 1./sqrt2/sigma
      elseif(ipia.eq.27)then
      rkm=rrup
      jf=ip11(ip)
      sigma = sigPez11(xmag,jf)
       ka=3
      gnd(1) = amean11(xmag,rkm,rjbp,vs30,ip,jf,ka)
      if(lceus_sigma)sigma=ceus_sigma !for special study 3/2011
       sigmaf = 1./sqrt2/sigma
      elseif(ipia.eq.2)then
      call getToro(ip,iq,1,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-2)then
      call getToro(ip,iq,2,xmag,rjb,gnd,sigma,sigmaf )
        elseif(ipia.eq.4)then
        call getAB06(ip,iperab(ip),irab(ip,ia),xmag,rrup,gnd,sigma,sigmaf,vs30)
        elseif(ipia.eq.19)then
      call getTP05(ip,ipertp(ip),irtb,xmag,rrup,gnd,sigma,sigmaf)
c        write(6,*)'iq clamp(iq) gnd rrup ip TP ',iq,clamp(iq),gnd,rrup,ip
        elseif(ipia.eq.20)then
      call getSilva(ip,isilva(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
        elseif(ipia.eq.28)then
      call getSilvaV(ip,isilva2(ip),irsilva,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.5)then      !to distinguish ab95 from fea
      call getFEA(ip,iq,2,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-5)then      !HR for atkinson-boore
      call getFEA(ip,iq,4,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.6)then
      call getFEA(ip,iq,1,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.-6)then      !hr
      call getFEA(ip,iq,3,xmag,rrup,gnd,sigma,sigmaf)
      elseif(ipia.eq.7) then
        call getSomer(ip,iq,1,xmag,rjb,gnd,sigma,sigmaf)
      elseif(ipia.eq.-7) then
        call getSomer(ip,iq,2,xmag,rjb,gnd,sigma,sigmaf)
             elseif(ipia.eq.10)then
        call getCampCEUS(ip,iq,1,xmag,rrup,gnd,sigma,sigmaf)
             elseif(ipia.eq.-10)then
        call getCampCEUS(ip,iq,2,xmag,rrup,gnd,sigma,sigmaf)
      endif
c This is a good place to apply the extra clamp constraint which
c in 2002 was used only for CEUS sources. The precomputed p() matrix does
c not work for this case because truncation is no longer at mu+3sig
      test0=gnd(1)+3.*sigma
      test= exp(test0)
c clamp is relying on the iq index. Some CEUS models now have different set than the 7 std
      if(clamp(iq).lt.test .and. clamp(iq).gt.0.) then
      clamp2= alog(clamp(iq))
      sigmasq=1./sigma/sqrt2
      tempgt3= (gnd(1) - clamp2)*sigmasq
      probgt3= (erf(tempgt3)+1.)*0.5
c Deterministic output must occur whether clamping is in effect or not.
        if(determ)then
        isbig=.false.
        isclose=.false.
c  Nico suggestion apr 3 2008. Can increase file size a lot. But can add reasonable candidates.
        if(rjb.le.rjbdet.or.xmag.ge.xmagdet)then
        isclose=.true.
        rjbdet=rjb
        xmagdet=xmag
c        isbig=.true.
        endif
c mx_index is  the index of the maximum magnitude in the set 1,...,nmagf(ift)
       if(imag.eq.mx_index(ift).and.isclose)write(idet,679)exp(gnd(1)),rrup,xmag,1./sigmaf/sqrt2,ipia,ift,xmag,rjb,rrup,
     + wttmp,dtor1,iftype(ift),1
      endif      !if determ
      do k=1,nlev(ip)
      temp= (gnd(1) - xlev(k,ip))*sigmasq
      temp1= (erf(temp)+1.)*0.5
      temp1= (temp1-probgt3)/(1.-probgt3)
      if(temp1.lt.0.) goto 3203      !no more calcs once p<0
      if(cluster)then
            prob_s(jseg(ift),ks,kg,k,ip)=prob_s(jseg(ift),ks,kg,k,ip)+ weight*temp1
      else
        prob(icnt,k,ip,1,kg)= prob(icnt,k,ip,1,kg)+ 
     &         weight*rate*temp1
        endif      !if cluster
        enddo      !k=1,nlev
        goto 3203
        endif
      elseif(ipia.eq.12)then
      call getMota(ip,iper(ip),xmag,rrup,vs30,gnd,sigmaf)
      elseif(ipia.eq.22)then
      call getDahle95(iperdahl(ip),isoild,xmag,rjb,gnd,sigmaf)
      endif
      wttmp=weight*rate
        do ifn=1,nfi(ip)
        if(determ)then
        isbig=.false.
        isclose=.false.
c  Nico suggestion apr 3 2008
        if(rjb.le.rjbdet.or.xmag.ge.xmagdet)then
        isclose=.true.
        rjbdet=rjb
        xmagdet=xmag
c        isbig=.true.
        endif
c mx_index is  the index of the maximum magnitude in the set 1,...,nmagf(ift)
       if(imag.eq.mx_index(ift).and.isclose)write(idet,679)exp(gnd(ifn)),1./sigmaf/sqrt2,ipia,ift,xmag,rjb,rrup,
     + wttmp,dtor1,iftype(ift),ifn
      endif      !if determ
        do  k=1,nlev(ip)
       pr=(gnd(ifn) - xlev(k,ip))*sigmaf
       if(pr.gt.3.3)then
       ipr=250002
       elseif(pr.gt.plim)then
       ipr= 1+nint(dp2*(pr-plim))      !3sigma cutoff n'(mu,sig)
       else
       goto 3201      !transfer out if ground motion above mu+3sigma
       endif
      if(cluster)then
            prob_s(jseg(ift),ks,kg,k,ip)=prob_s(jseg(ift),ks,kg,k,ip)+ weight*p(ipr)
      else
      cfac = wttmp*p(ipr)
        prob(icnt,k,ip,ifn,kg)= prob(icnt,k,ip,ifn,kg)+ cfac
               if(deagg)then
c for individual fault
      frbar(ift,ifn,ip)=frbar(ift,ifn,ip)+cfac*rrup
      fmbar(ift,ifn,ip)=fmbar(ift,ifn,ip)+cfac*xmag
               eps= -pr*sqrt2
      febar(ift,ifn,ip)=febar(ift,ifn,ip)+cfac*eps
      fhaz(ift,ifn,ip)=fhaz(ift,ifn,ip)+cfac
c ifn index is present: 2nd to last frontier. ieps is  an explicit dimension, recomm. by Bazzuro
c Some attn. models may claim a given M,R is a low eps0 combination, others will say a higher eps0.
c 
      if(eps.lt.emax)then
      ieps=max(1,min(int((eps+2.)*2.),10))
      rbar(ir,im,ieps,ifn,ip)=rbar(ir,im,ieps,ifn,ip)+cfac*rrup
      mbar(ir,im,ieps,ifn,ip)=mbar(ir,im,ieps,ifn,ip)+cfac*xmag
      ebar(ir,im,ieps,ifn,ip)=ebar(ir,im,ieps,ifn,ip)+cfac*eps
      haz(ir,im,ieps,ifn,ip)=haz(ir,im,ieps,ifn,ip)+cfac
      if(eps.ge.2.)goto 3201
c Below keeps the books on epsilon distribution for e0<2. if e0>2 dont bother.
      ka=5
      temp=p(ipr)
      prx= temp-ptail(ka)
      dowhile(prx.gt.1.e-9)
      prob5(ir,im,ieps,ip,ifn,ka)=prob5(ir,im,ieps,ip,ifn,ka)+wttmp*prx
c  *prlr       this factor is not present now. Why?
      ka=ka-1
      if(ka.eq.0)goto 3201
      prx= temp-ptail(ka)
      enddo
      endif      !eps <emax
      endif      !if deagg
        endif      !if cluster
      enddo      !levels
3201      continue
      enddo      !files
 3203   continue
      enddo      !period set
3303      continue      !if char event rate is tiny skip above calcs and go here
      enddo      !imag set
        endif
ccccccccccccccccccc
 850  continue
c all faults have been dealt with; if clustered events we next convert segment rates to mean rates in prob() array
c New May 17 2007. 
       if(cluster)then
       call cluster_me(prob,prob_s,wtscene,nper,nlev,nscene,
     +  ngroup,jsegmin,jsegmax,icnt)
       endif
c--- write output files every 1000 sites
 860    continue
        if(grid.and.((icnt.eq.1000).or.(i.eq.nrec))) then
      write(6,*) i,prob(icnt,1,1,1,1)
        iend= 1000
        if(i.eq.nrec) iend=icnt
      if(cluster)then
        do 315 kg=1,ngroup(1)
c the above do 315 is to separately reorganize data for each epistemic  branch location NMSZ. 
c has to be dominant loop here. For Clustered events, clustered-event rate is
c applied now (this is an important distinction from non-clustered models).
        do 315 ip=1,nper
        ndata= iend*nlev(ip)
        do 350 ii=1,ndata
        i2= ii-1
        irec= i2/nlev(ip)
        ilev= i2-irec*nlev(ip)
        irec= irec+1
        ilev= ilev+1
        out(ii)= rate_cl(kg)*prob(irec,ilev,ip,1,kg)
 350    continue
        call putbufx(ifp(ip,kg,1),out,ndata,readn)
        readn= readn/4
c        write(6,*) ndata,readn
 315    continue
       icnt=0
      prob=1.e-21
      else
        do 115 ip=1,nper
        ndata= iend*nlev(ip)
        do 115 ifn = 1,nfi(ip)
c the above do 115 is to separately reorganize data for each epistemic gm branch. 
c has to be dominant loop here.
        do 50 ii=1,ndata
        i2= ii-1
        irec= i2/nlev(ip)
        ilev= i2-irec*nlev(ip)
        irec= irec+1
        ilev= ilev+1
        out(ii)= prob(irec,ilev,ip,ifn,1)
   50   continue
        call putbufx(ifp(ip,1,ifn),out,ndata,readn)
        readn= readn/4
c        write(6,*) ndata,readn
 115    continue
       endif
        icnt= 0
      prob=1.e-21
       elseif(.not.grid)then
67      format(/,'#Station ',f7.4,1x,f10.4,1x,a,' vs30 ',f6.1)       
      fac=1.0
      do 615 ip=1,nper
      do 615 kg=1,ngroup(ip)
       if(kg.eq.1)write(icc(ip,kg),67)slat(i),slong(i),sname(i),vs30
       if(cluster)fac=rate_cl(kg)
      do 614 ifn=1,nfi(ip)
      if(.not.cluster)then
      cmt=pithy(ifn)
      else
      cmt=g_name(kg)
      endif
c write everything to icc(ip,1) new may 1 2008. Separate branch curves
c with white space
      write(icc(ip,1),69) period(ip), nlev(ip),cmt

69      format(/,'#',f9.2,1x,i2,1x,a12)
      do 614 k=1,nlev(ip)

      write(icc(ip,1),612) ylev(k,ip), fac*prob(i,k,ip,ifn,kg)
 612      format(f10.6,1x,e11.5) 
 614  continue
 615  continue
c       write(10,68)sname(i)
c 68      format('#End hazcurves for ',a,/)      
        endif
  100 icnt= icnt+1
      dum=dtime(tarray)
       if(slist)then
       do ip=1,nper
       do kg=1,ngroup(ip)
       close(icc(ip,kg))
       enddo
       enddo
       write(6,*)'An ascii file was written for each period in input'
       do i=1,nrec
       close(30+i)
       enddo
       write(6,*)'A file of eqrates was written for each station'
      if(deagg)then
       do ip=1,nper
      print *,ip,immax,' pd index and immax deagg output'
       do ifn=1,nfi(ip)
       do ieps=1,10
       do im=1,immax
       do ir=1,16
      prx=haz(ir,im,ieps,ifn,ip)
       if(prx.gt.tiny)then
       e0=max(-9.99,ebar(ir,im,ieps,ifn,ip)/prx)
c Need to keep tract of important faults. May need another array to write out
c ifn below lets downstream programs know which gm uncert branch this item came from       
       write(ip+20,25)rbar(ir,im,ieps,ifn,ip)/prx,
     + mbar(ir,im,ieps,ifn,ip)/prx,prx,(prob5(ir,im,ieps,ip,ifn,k),k=5,1,-1),
     + ifn,e0
      endif
25      format(f9.3,1x,f7.3,6(1x,e11.5),
     + 1x,i2,1x,f5.2,1x,i4)
       enddo      !ir 
       enddo      !im (I code, therefore I R)
       enddo      !ieps
       enddo      !ifn (typically NGA gm median uncertainty)
c write some mean deagg info for each fault with more than tiny contribution
       do iflt=1,nft
       prx=0.0
       faebar=0.0
       fambar=0.0
       farbar=0.0      !fault rbar
       wtf=1.0
       do ifn=1,nfi(ip)
       if(nfi(ip).eq.3)then
       wtf=gmwt(ifn)
       elseif(nfi(ip).gt.3)then
       stop'notready for deagg with nfi not 1 or 3'
       endif
       prx=prx+fhaz(iflt,ifn,ip)*wtf
       faebar=faebar+febar(iflt,ifn,ip)*wtf
       farbar=farbar+frbar(iflt,ifn,ip)*wtf
       fambar=fambar+fmbar(iflt,ifn,ip)*wtf
       enddo      !ifn
       if(prx.gt.tiny)write(40+ip,27)iflt,farbar/prx,
     + fambar/prx,faebar/prx,prx,nfi(ip)
27      format(i3,1x,f9.3,1x,f7.3,1x,f5.2,1x,e11.5,4x,i2)
      enddo	!iflt
       close(ip+20)
       close(40+ip)
       enddo      !ip (.true.)
      endif      !if deagg, true if one receiver with loc on command line.
       endif
       write(*,*)tarray(1),' s = time to run hazFXnga13l.'
       stop
202   print *,fname,' not found'
       stop 'put in W.D.'
2014      print *,'GR/BSSAcoef.dat not found. coef file for BSSA model needed'
      stop 'put in the folder indicated'
      end
c
      subroutine delaz(sorlat,sorlon,stnlat,stnlon,delta,az,baz)
      parameter (coef= 3.1415926/180.)
      xlat= sorlat*coef
      xlon= sorlon*coef
      st0= cos(xlat)
      ct0= sin(xlat) 
      phi0= xlon
      xlat= stnlat*coef
      xlon= stnlon*coef
      ct1= sin(xlat)
      st1= cos(xlat)
      sdlon= sin(xlon-phi0)
      cdlon= cos(xlon-phi0)
      cdelt= st0*st1*cdlon+ct0*ct1
      x= st0*ct1-st1*ct0*cdlon
      y= st1*sdlon
      sdelt= sqrt(x*x+y*y)
      delta= atan2(sdelt,cdelt)
      delta= delta/coef
      az= atan2(y,x)
      az= az/coef
      x= st1*ct0-st0*ct1*cdlon
      y= -sdlon*st0
      baz= atan2(y,x)
      baz= baz/coef
      delta= delta*111.2
      return
      end         

      subroutine back(delta0,az,sorlat,sorlon,stnlat,stnlon)
      parameter (coef= 3.1415926/180.)
      delta= delta0/111.11
      delta= delta*coef
      sdelt= sin(delta)
      cdelt= cos(delta)
      xlat= sorlat*coef
      xlon= sorlon*coef
      az2= az*coef
      st0= cos(xlat)
      ct0= sin(xlat)
      phi0= xlon
      cz0= cos(az2)
      ct1= st0*sdelt*cz0+ct0*cdelt
      x= st0*cdelt-ct0*sdelt*cz0
      y= sdelt*sin(az2)
      st1= sqrt(x*x+y*y)
      dlon= atan2(y,x)
      stnlat= atan2(ct1,st1)
      stnlat= stnlat/coef
      stnlon= phi0+dlon
      stnlon= stnlon/coef
      return
      end

      subroutine gettab(gma,f,xmag,rlog)
c      gets corresponding log PSA values (cm/s**2) from interpolation
c      of file ab94.tab.
      dimension xmag(20),gma(22,30,20), rlog(30), f(20)
      character*60 header
c      Written by G.M. Atkinson, Jan. 1994
C----  modified by A. Frankel 10/95
c
c     Read AB94 simulated grd. motion table.
      open(3,file='ab94.tab',status='old')
 900  format(i2)
      do 20 j = 1, 5
        read(3,10)header
        read(3,900) idum
10      format(a)
20    continue
      nf = 11
      nmag = 14
      nd = 18
      jdl = 0
      do 30 jf = 1, nf-2
        f(jf) = 10.**(-0.3 + float(jf-1)*0.2)
30    continue
      f(nf) = 89.
      f(nf-1) = 99.
      do 50 jm = 1, nmag
        read(3,900) idum
        read(3,900) idum
        read(3,35)xmag(jm)

35      format(13x,f5.2)
        read(3,900) idum

        read(3,45)
        read(3,900) idum
        read(3,45)
        read(3,900) idum
        do 40 jd = 1, nd
          read(3,45)rlog(jd),(gma(jf,jd,jm), jf = 1,nf)
          read(3,900) idum
45        format(f6.2,2x,11f6.2)
40      continue
50    continue
      close(3)
c      write(6,*) gma(1,1,1)
      return
      end subroutine gettab

c***************************************************************
c*****************************************************************

      subroutine table(gma,amag,R,freq,amean,f,xmag,rlog)
      dimension xmag(20),gma(22,30,20), rlog(30), f(20)
      dimension avg(20)
c     Use interpolation to get amean for given freq(jfreq), amag, R.
c
      nf= 11
      nmag= 14
      nd= 18

      rl = alog10(R)

      jfl = 0

c     Find bounding jf values.

      do 10 j = 1,nf-3

        if (freq .ge. f(j) .and. freq .lt. f(j+1))jfl = j

10    continue

      if (abs(freq-f(1)) .lt. 0.01) jfl = 1

      if (abs(freq-f(nf)) .lt. 0.01) jfl = nf

      if (abs(freq-f(nf-1)) .lt. 0.01) jfl = nf-1

      if (abs(freq-20.) .lt. 0.01)  jfl = nf-2

      if (jfl .eq. 0) write(6,*) ' ERROR. FREQ. OUTSIDE RANGE.'

      if (freq .gt. 20. .and. freq .lt. 89.)

     *   write(6,*)' ERROR. FREQ. OUTSIDE RANGE'

      jfu = jfl + 1

      fracf = 0.

      if (freq .ge. 20.) go to 15

      fracf = (alog10 (freq) - alog10 (f(jfl)) ) /

     *        (alog10 (f(jfu)) - alog10 (f(jfl)) )

15    continue

c     Find bounding jm values.

      jml = 0

      do 20 j = 1, nmag-1

         if (amag .ge. xmag(j) .and. amag .lt. xmag(j+1)) jml = j

20    continue
      if (amag .eq. xmag(nmag)) jml = nmag
      if (jml .eq. 0) write(6,*)' ERROR. MAGNITUDE OUTSIDE RANGE.'
      jmu = jml + 1
      fracm = (amag - xmag(jml)) / (xmag(jmu)-xmag(jml))
c     find bounding distance values.
      jdl= 0
      do 30 j = 1, nd-1
        if (rl .ge. rlog(j) .and. rl .lt. rlog(j+1)) jdl= j
30    continue
      if (rl .eq. rlog(nd)) jdl = nd
      if (jdl .eq. 0) write(6,*)' ERROR. DISTANCE OUTSIDE RANGE.'
      jdu = jdl + 1
      fracd = (rl - rlog(jdl)) / (rlog(jdu)-rlog(jdl))
      do 40 j = jfl, jfu
        arl = gma(j,jdl,jml) + fracm * (gma(j,jdl,jmu)-gma(j,jdl,jml))
        aru = gma(j,jdu,jml) + fracm * (gma(j,jdu,jmu)-gma(j,jdu,jml))
c       arl is the amplitude for lower dist. bound, so arl>aru
        avg(j) = arl - fracd*(arl-aru)
40    continue
      amean = avg(jfl) + fracf * (avg(jfu) - avg(jfl))
      return
      end

cccccccccccccccccccccccccccccc
      subroutine getSpu2000(ip,iq,xmag,dist0,vs30,gndout,sigmaf)
c commandeered by Spudich ea. from BJF93 of "hazFXv3" etc.
c To distinguish from BJF97. 7 periods.
c Inputs:
c iq is the index of the current period, see perx(iq)
c dist0 is Rjb
c xmag = M or moment magnitude
c vs30 = site Vs top 30 m
c Outputs gnd = ln(median motion,g)
c        sigmaf = 1/sigma/sqrt2  where sigma is dispersion natural log
c
c adapted to Fortran 95 july 26 2006.SH
      parameter (sqrt2=1.414213562, pi=3.141592654,np=7,aln10=2.30258509)
      parameter (vref=760.)
c spudich coeffs...BSSA Oct 1999
      real magmin,perx(8),gndout(3)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8)
      logical l_gnd_ep(8)

c spudich coeffs...
c site amp: for now use that of BJF97. THis should be reviewed July26 2006.
      real bv(np),va(np)
        real b1(np),b2(np),b3(np),b4(np),b5(np),h(np),sigma(np)
c array constructors OCT 2006. No coeffs corresponding to PGV, perx(8)
        perx=(/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
        bv = (/-.371,-.292,-.698,-.212,-.401,-.553,-.655 /)
        va = (/1396.,2118.,1406.,1112.,2133.,1816.,1795. /)
        b1 = (/0.299,2.224,2.276,2.144,2.263, 2.292, 2.168/)           
        b2 = (/ 0.229, 0.309, 0.450,0.327, 0.334,0.384,0.471/)
        b3 = (/0., -0.090,-.014,-0.098,-0.070,-0.039, -0.037/)
        b4 = (/ 0., 0.,0.,0.0,0.,0.,0.0/)
        b5 =(/-1.052,-1.047,-1.083,-1.250,-1.020,-1.038,-1.049/)
        h = (/7.27,8.63, 6.01,9.99,7.72,6.70, 6.71/)
         sigma= (/0.225,0.262,0.301,0.295,0.263 ,0.275,0.341/)
         sig=sigma(iq)*aln10
          sigmaf = 1./(sig*sqrt2)
c site amp from BJF97. THis will need review. Placeholder july 27 2006.
      period=perx(iq)
       dist= sqrt(dist0*dist0+h(iq)**2)
          gnd=  b1(iq) +b2(iq)*(xmag-6.) +b3(iq)*(xmag-6.)**2
          gnd= gnd + b5(iq)*alog10(dist)      !+ b4(iq)*dist nuisance is 0
c--- convert from psv in cm/sec to psa in g
          if(period.ne.0.) gnd=gnd+alog10(2.*pi/(980.*period))
c bse 10 to base e
c site amp is last term because it is a natural log term. 
c Not too useful when soil siteamp has known nonlinear response.
          gnd= gnd *aln10 +bv(iq)*alog(vs30/vref)
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
      return
      end subroutine getSpu2000

cccccccccccccccc
      subroutine getToro(iper,ip,ir,xmag0,dist0,gndout,sigma,sigmaf)
c midcontinent w/moment mag. Adapted to NGA code SH June 2006. 7 periods used in2002.
c ip = index in perx() array, 
c ir=1 use BC rock; ir=2 use hardrock model.
c Hard-rock in tc1h, otherwise same regression model & coeffs.
c dist0 is rjb (km). Eqn 4 SRL converts to "R_M".
c median motions are capped. Max motions are not really constrained here.
c Warning: Clamp business will have to be done in main.
c replaced data statements with array constructors oct 2006.
        parameter (sqrt2=1.414213562, pi=3.141592654)
c the gnd_ep branching will not be done for CEUS relations. 
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/soils/vs30
      real gnd_ep(3,3,8),gndout(3)
       logical lceus_sigma,l_gnd_ep(8)
      common/ceus_sig/lceus_sigma,ceus_sigma

             real, dimension(7):: tc1,tc2,tc3,tc4,tc5,tc6
      real, dimension(7):: tc1h,th,tsigma,clamp
           real   perx(8) , pganl
           save pganl
c array constructors
       perx = (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
c Below tc coeffs correspond to midcontinent, equations using moment mag.
      tc1 = (/2.619,2.295,0.383,2.924,1.8823,1.288,-0.558/)
      tc1h = (/2.20,1.73,0.09,2.37,1.34,0.8306,-0.74/)      
      tc2 = (/0.81,0.84,1.42,0.81,0.964,1.14,1.86/)
      tc3 = (/0.,0.0,-0.2,0.,-0.059,-0.1244,-0.31/)
      tc4 = (/1.27,0.98,0.90,1.1,0.951,0.9227,0.92/)
      tc5 = (/1.16,0.66,0.49,1.02,0.601,0.5429,0.46/)
      tc6 = (/0.0021,0.0042,0.0023,0.004,0.00367,0.00306,0.0017/)
      
           th = (/9.3,7.5,6.8,8.3,7.26,7.027,6.9/)
c write sigma in nat log units. Saves a divide
           tsigma = (/0.7506,0.7506,0.799,.7506,.7506,.7506,0.799/)      
           clamp = (/3.,6.,4.,6.,6.,6.,6./)
        xmag= xmag0
c correct mag will be used in calling program. never need to convert here.
      if(lceus_sigma)then
      sigma=ceus_sigma	!special study 3/2011
      else
          sigma= tsigma(ip)
          endif
          sigmaf= 1./sigma/sqrt2
        cor = exp(-1.25 + 0.227*xmag)
c        dist= sqrt(dist0*dist0+th(ip)*th(ip))
c Modified Toro for finite-fault Oliver. Nov 2006.The empirical model eqn 4-3 Toro 2002
          dist= sqrt(dist0*dist0+th(ip)*th(ip)*cor*cor)
      if(ir.eq.1)then
c firm rock const coeff.
          gnd= tc1(ip)
          else
c hard rock. Could have other possibilities as well. For now default to h.r.
          gnd=tc1h(ip)
          endif
          gnd=gnd+tc2(ip)*(xmag-6.)+ tc3(ip)*((xmag-6.)**2)
          gnd= gnd-tc4(ip)*alog(dist)-tc6(ip)*dist
          factor= alog(dist/100.)
          if(factor.gt.0.) gnd= gnd-(tc5(ip)-tc4(ip))*factor
c---following is for clipping gnd motions: 1.5g PGA, 3.0g for rest 
          if(ip.eq.1.or.ip.eq.7) then
c for 2s I also put the PGA 1.5g median limit (g). SH June 2006.
           gnd=min(gnd,0.405)
          else
           gnd=min(1.099,gnd)
          endif
          if(ip.eq.1)pganl=gnd      !pga for ip=1.
c          if(lsoil)gnd = gnd+basiteamp(pganl,vs30,760.)
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
      return
      end subroutine getToro

cccccccccccccccccc
      subroutine getsadigh(ip,iq,xmag,rrup,gndout,sigmaf)
c modified to perform like nga functions. SH June 2 2006 workin progress...
c iq is the index of perx() associated with period being run with index ip in main
c we are not using sigmanf or distnf at this time in hazFXnga13l
        parameter (sqrt2=1.414213562)
c Rock-site coefs only, written below for 7 periods. How fast is that rock?
c Walt Silva says SM sites avg about 500 m/s (NGA discussion hot topic). 
c This code does
c not vary the median response with rock Vs30. You get one. That is it.
c replaced data statements with array constructors oct 2006.
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

      real sc1(7),sc2(7),sc3(7),sc4(7),sc5(7),sc6(7),perx(7)
      real sd1(7),sd2(7),sd3(7),sd4(7),sd5(7),sd6(7),sd7(7)      
      real ssigma1(7),ssigma2(7),ssigmacoef(7),smagsig(7)
        common/mech/ss,rev,normal,obl
c        common/dipinf/dipang,cdipsq
        logical ss,rev,normal,obl
c prepared for 7 periods jun 22 2006. gnd is modified
c for sense-of-slip  from common/mech/ 
      real sh(7)/7*5./
      real sthrust(7)/7*.1823/
c sc1 and sd1 are from the SRL table and might correspond to 600 m/s? rock
      sc1 = (/-0.624,.153,-1.705,0.2750,-0.057,-0.5880,-2.945/)
      sc2 = (/1.0,1.,1.,1.,1.,1.,1./)
      sc3 = (/0.0,-0.004,-0.055,0.006,-0.017,-0.04,-0.07/)
      sc4 = (/-2.1,-2.08,-1.8,-2.148,-2.028,-1.945,-1.67/)
      sc5 = (/1.29649,1.29649,1.29649,1.29649,1.29649,1.29649,1.29649/)
      sc6 = (/0.25,.25,.25,.25,.25,.25,.25/)
      perx = (/0.0,0.2,1.0,0.1,0.3,0.5,2.0/)
      sd1 = (/-1.274,-.497,-2.355,-0.375,-0.707,-1.238,-3.595/)
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
        xmag0=xmag
c        if(icode(ip,ia).eq.1) xmag= 2.715 -0.277*xmag0+0.127*xmag0*xmag0
          if(xmag.lt.smagsig(iq))then
           sig= ssigma1(iq)- ssigmacoef(iq)*xmag
          else
           sig= ssigma2(iq)
           endif
          sigmaf= 1./(sig*sqrt2)
          dist= sqrt(rrup*rrup + sh(iq)**2)
          if(xmag.lt.6.5) then
      gnd= sc1(iq)+sc2(iq)*xmag+sc3(iq)*((8.5-xmag)**2.5)
      gnd= gnd+ sc4(iq)*alog(dist+ exp(sc5(iq)+sc6(iq)*xmag))
           else
      gnd= sd1(iq)+sd2(iq)*xmag+sd3(iq)*((8.5-xmag)**2.5)
      gnd= gnd+ sd4(iq)*alog(dist+ exp(sd5(iq)+sd6(iq)*xmag))
      if(rev)gnd=gnd + sthrust(iq)
            endif
          if(sd7(iq).ne.0.)gnd= gnd+sd7(iq)*alog(rrup+2.)
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
      return
      end subroutine getsadigh

c getCamp removed: no longer used. Replace with getCamp2003
c getAB95 removed. Seems to be too archaic.

      subroutine getFEAtab(ip,iq)
c  hardrock and   BC rock tables can be read in here.
c 
c replaced data statements with array constructors oct 2006.
      parameter (np=7)
      real gma(20,21,4*np),tabdist(21,4)
      logical isok
      common/gmat/gma,tabdist
      character*30 nametab(np),nameab(np),hardtab(np),hardab(np)
      character*4 adum
c assumes these files are in working directory or GR/ sub-directory.
      nametab=(/'pgak01l.tbl ','t0p2k01l.tbl','t1p0k01l.tbl',
     1 't0p1k01l.tbl','t0p3k01l.tbl','t0p5k01l.tbl','t2p0k01l.tbl'/)
      nameab=(/'Abbc_pga.tbl','Abbc0p20.tbl','Abbc1p00.tbl',
     1 'Abbc0p10.tbl','Abbc0p30.tbl','Abbc0p50.tbl','Abbc2p00.tbl'/)
       hardtab=(/'pgak006.tbl ','t0p2k006.tbl','t1p0k006.tbl',
     1 't0p1k006.tbl','t0p3k006.tbl','t0p5k006.tbl','t2p0k006.tbl'/)
      hardab=(/'ABHR_PGA.TBL','ABHR0P20.TBL','ABHR1P00.TBL',
     1'ABHR0P10.TBL','ABHR0P30.TBL','ABHR0P50.TBL','ABHR2P00.TBL'/)
      if(iq.eq.1)then            !fea BCtable
      inquire(file=nametab(ip),exist=isok)
      if(.not.isok)then
c look in GR subdirectory if not in working dir
      nametab(ip)='GR/'//nametab(ip)
      inquire(file=nametab(ip),exist=isok)
      if(.not.isok)goto 20
      endif
         open(unit=15,file=nametab(ip),status='old')
         ipq=ip
         elseif(iq.eq.2)then      !ab BC table
      inquire(file=nameab(ip),exist=isok)
      if(.not.isok)then
c look in GR subdirectory if not in working dir
      nameab(ip)='GR/'//nameab(ip)
      inquire(file=nameab(ip),exist=isok)
      if(.not.isok)goto 21
      endif
         open(unit=15,file=nameab(ip),status='old')
         ipq=ip+np
         elseif(iq.eq.3)then      !fea A table
      inquire(file=hardtab(ip),exist=isok)
      if(.not.isok)then
c look in GR subdirectory if not in working dir
      hardtab(ip)='GR/'//hardtab(ip)
      inquire(file=hardtab(ip),exist=isok)
      if(.not.isok)goto 23
      endif
         open(unit=15,file=hardtab(ip),status='old')
         ipq=ip+np+np
         elseif(iq.eq.4)then      !ab A table
      inquire(file=hardab(ip),exist=isok)
      if(.not.isok)then
c look in GR subdirectory if not in working dir
      hardab(ip)='GR/'//hardab(ip)
      inquire(file=hardab(ip),exist=isok)
      if(.not.isok)goto 24
      endif
         open(unit=15,file=hardab(ip),status='old')
         ipq=ip+np*3
         endif
         read(15,900) adum
900      format(a)
         do 80 idist=1,21
   80    read(15,*) tabdist(idist,iq),(gma(imag,idist,ipq),imag=1,20)
         close(15)
      return
20      write(6,*)'getFEAtab open failed for ',nametab(ip)
      stop 'hazFXnga13l: put in working dir or rewrite program'
21      write(6,*)'getFEAtab  open failed for ',nameab(ip)
      stop 'hazFXnga13l: put in working dir or rewrite program'
23      write(6,*)'getFEAtab  open failed for ',hardtab(ip)
      stop 'hazFXnga13l: put in working dir or rewrite program'
24      write(6,*)'getFEAtab  open failed for ',hardab(ip)
      stop 'hazFXnga13l: put in working dir or rewrite program'
      end subroutine getFEAtab

ccccccccccccccccccccccc
      subroutine getFEA(iper,ip,iq,xmag,dist0,gndout,sigma,sigmaf)
c table info is available, from common/gmat/. See getFEAtab.
c this routine should work for finite CEUS faults. What if 0 depth suchas meers?
c four tables are set up and good to go as of july 28 2006.
c if iq = 1 use feaBC. if  iq=2, use ab95 BC table. 
c if iq = 3 fea A rock, and iq=4, ab95 A rock table.
c primary outputs:
c gnd=log median motion
c sigma = logarthmic s.d.
c sigmaf = 1/sigma/sqrt(2)
c inputs:
c ip=period index in perx() below. This controls which period coeffs are used
c iq=table index used with ip for table lookup
c xmag= moment mag (only)
c dist0= hypo distance used in getfea
c sigmanf, distnf = nearfield sigma apply when distance < distnf (not used 2002)
c replaced data statements with array constructors oct 2006.
c
        parameter (np=7,sqrt2i=0.707106781)
      common/gmat/gma,tabdist
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/ceus_sig/lceus_sigma,ceus_sigma
      real gnd_ep(3,3,8),gndout(3)
      logical lceus_sigma,l_gnd_ep(8),sp

      real bdepth(7),bsigma(7),clamp(7),xlogfac(7)
      real gma(20,21,4*np),tabdist(21,4)
             real perx(8)      !a reference set including pgv which wasn't used in 2002
        perx= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./ )
           bdepth= (/5.0,5.0,5.0,5.0,5.0,5.0,5.0/)
           bsigma= (/0.7506,0.7506,0.799,0.7506,0.7506,.7506,0.799/)      !converted to nat log
c---- following for  FEA or Atkinson-Boore look-up table
      ipq=ip+np*(iq-1)
c        if(icode(ip,ia).eq.1) xmag= 2.715 -0.277*xmag0+0.127*xmag0*xmag0
          sp = perx(ip).gt.0.02 .and.perx(ip).lt.0.5
          if(lceus_sigma)then
          sigma=ceus_sigma      !for special study 3/2011
          else
          sigma= bsigma(ip)      !in correct log base already
          endif
          sigmaf= sqrt2i/sigma
          dist= dist0
ccccccccc if fault depth0 ne 0. then this depth term not needed 
c          dist= sqrt(dist0*dist0+bdepth*bdepth)
          if(dist.lt.10.) dist=10.
          imag= (xmag-4.4)/0.2 +1
          rdist= alog10(dist)
          idist= (rdist-tabdist(1,iq))/0.1 + 1
          xm1= (imag-1)*0.2 + 4.4
          fracm= (xmag-xm1)/0.2
          xd1= (idist-1)*0.1 +tabdist(1,iq)
          fracd= (rdist-xd1)/0.1
          idist1= idist+1
          if(idist.gt.21) idist=21
          if(idist1.gt.21) idist1=21
          gm1= gma(imag,idist,ipq)
          gm2= gma(imag+1,idist,ipq)
          gm3= gma(imag,idist1,ipq)
          gm4= gma(imag+1,idist1,ipq)
          arl= gm1 + fracm*(gm2-gm1)
          aru= gm3 + fracm*(gm4-gm3)
          gnd= arl +fracd*(aru-arl)
c          gnd= gnd+ xlogfac      !nuisance. is 0
          gnd= gnd/0.4343
c--- following is for clipping 
          if(ip.eq.1)then
c pga has index ip=1
           gnd=min(0.405,gnd)
           elseif(sp)then
c No Clipping for LP: 1s SA haz index 3. 2s has period index 7  Noo clipping for2hz either.
           gnd=min(gnd,1.099)
           endif
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
  103 continue
      return
      end subroutine getFEA

ccccccccccccccccccccccccc
      subroutine getSomer(iper,ip,ir,xmag0,dist0,gndout,sig,sigmaf)
c---- Somerville et al. (2001) for CEUS
c input dist0 (rjb ) (km) and xmag0 (moment mag) 
c Input ir to control rock conditions:
c ip = index in period array that is currently of interest, ip=1 => PGA for example
c ir=1 BC or firm rock
c ir=2 hard rock (2800 m/s?)
c output gnd=log(median) and sigmaf=1/sigma/sqrt2
c replaced data statements with array constructors oct 2006.
c
      parameter(dist1=50.358713,sqrt2=1.41421356)
c      dist1= sqrt(50.*50.+ 6.*6.) in above parameter statement
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/ceus_sig/lceus_sigma,ceus_sigma
      real gnd_ep(3,3,8),gndout(3)
      logical lceus_sigma,l_gnd_ep(8)

      real a1(7),a2(7),a3(7),a4(7),a5(7),a6(7),a7(7),sig0(7)
c enter statements with coeff values.
             real perx(8),a1h(7)      !a reference set including pgv which wasn't used in 2002
      perx = (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
      a1 = (/0.658,1.358,-0.0143,1.442,1.2353,.8532,-0.9497/)
      a1h = (/0.239,0.793,-0.307,0.888,0.6930,0.3958,-1.132/)
      a2 = (/0.805,0.805,0.805,0.805,0.805,0.805,0.805/)
      a3 = (/-0.679,-.679,-.696,-.679,-.67023,-.671792,-0.728/)
      a4 = (/0.0861,0.0861,.0861,.0861,0.0861,.0861,.0861/)
      a5 = (/-0.00498,-.00498,-0.00362,-.00498,-.0048045,-.00442189,-0.00221/)
      a6 = (/-0.477,-.477,-0.755,-.477,-.523792,-.605213,-.946/)
      a7 = (/0.,0.,-0.102,0.,-.030298,-.0640237,-.140/)
      sig0 = (/0.587,0.611,0.693,0.595,.6057,.6242,0.824/)
c      clamp = (/3.,6.,3.,6.,6.,6.,6./)
c compute SOmerville median and dispersion estimates.
      if(lceus_sigma)then
      sig=ceus_sigma	!special study 3/2011
      else
        sig= sig0(ip)
        endif
      sigmaf= 1/sig/sqrt2
      dist= sqrt(dist0*dist0+ 36.)
      if(ir.eq.1)then
      gnd=a1(ip)
      else
      gnd=a1h(ip)
c hard rock/firm rock variation in Somerville only affects a1 coef.
      endif      
      xmag= xmag0
c      if(icode(ip,ia).eq.1) xmag= 2.715 -0.277*xmag0+0.127*xmag0*xmag0
      if(dist0.lt.50.) then
      gnd= gnd + a2(ip)*(xmag-6.4) + a3(ip)*alog(dist)
     1 +a4(ip)*(xmag-6.4)*alog(dist)
     2  + a5(ip)*dist0 + a7(ip)*(8.5-xmag)*(8.5-xmag)
        else
      gnd= gnd + a2(ip)*(xmag-6.4) + a3(ip)*alog(dist1)
     1 +a4(ip)*(xmag-6.4)*alog(dist)+a5(ip)*dist0
     2 + a6(ip)*(alog(dist)-alog(dist1))+a7(ip)*(8.5-xmag)*(8.5-xmag)
      endif
c      write(16,*) period,dist0,xmag,exp(gnd)
          if(ip.eq.1.or.ip.eq.7) then
c for long period I also put the PGA  limit (g). SH June 2006.
           gnd=min(gnd,0.405)
          else
           gnd=min(1.099,gnd)
          endif
c apply before or after these limits? here, trying after.
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
      return
      end subroutine getSomer

ccccccccccccccccccccccccc
      subroutine getAS97(iper,ip,xmag0,dist0,hw,gndout,sigmaf)
c modified for nga style, 7 periods, rock. SHarmsen June 2006.
c input xmag0=moment mag
c dist0 = r_cd or r_rup (km)
c : near-field terms to modify s.d. 0 if none.
c hw = logical variable, true if site on hanging wall
c output:
c gnd,sigmaf
c replaced data statements with array constructors oct 2006.
        parameter (sqrt2=1.414213562)
      common/mech/ss,rev,normal,obl
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

      logical hw      !.true. if site on hanging wall established in calling pgm
      logical ss,rev,normal,obl
      real, dimension(7):: as1,as2,as3,as4,as5,as6
      real, dimension(7)::  as9,as12,as13,asc1,asc4,b5,b6
             real perx(8)      !a reference set including pgv which wasn't used in 2002
        perx= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
      as1= (/1.64,2.406,0.828,2.16,2.114,1.615,-0.15/)
      as2= (/0.512,0.512,0.512,0.512,0.512,0.512,0.512/)
      as3= (/-1.145,-1.115,-0.8383,-1.145,-1.035,-0.9515,-0.725/)
      as4= (/-0.144,-0.144,-0.144,-0.144,-0.144,-0.144,-0.144/)
      as5= (/0.61,0.610,0.49,0.610,0.610,0.5810,0.40/)
      as6= (/0.26,0.26,0.013,0.26,0.198,0.119,-0.094/)
      as9= (/0.37,0.37,0.281,0.37,0.37,0.37,0.16/)
      as12= (/0.,-0.0138,-0.102,0.0280,-0.036,-0.0635,-0.14/)
      as13= (/0.17,0.17,0.17,0.17,0.17,0.17,0.17/)
      asc1= (/6.4,6.4,6.4,6.4,6.4,6.4,6.4/)
      asc4= (/5.6,5.1,3.7,5.5,4.8,4.3,3.5/)
      b5= (/0.70,0.77,0.83,0.74,0.78,0.8,0.85/)
      b6= (/0.135,0.135,0.118,0.135,0.135,0.13,0.105/)
c from 2003 input files: no correction (20% median redux?) for normal faults. Is this
c the current understanding?
c-- normal faults do not use the HW term in our 2002 model
        xmag= xmag0
        if(xmag.lt.7.0)then
         sig= b5(ip)- b6(ip)*(xmag-5.0)
        else
         sig= b5(ip)- b6(ip)*2.0
         endif
        sigmaf= 1./(sig*sqrt2)
          r= sqrt(dist0*dist0+ asc4(ip)*asc4(ip))
          if(xmag.le.asc1(ip)) then
          gnd= as1(ip)+as2(ip)*(xmag-asc1(ip))+as12(ip)*((8.5-xmag)**2)
     &         + (as3(ip)+ as13(ip)*(xmag-asc1(ip)))*alog(r)
          
          else 
          gnd= as1(ip)+as4(ip)*(xmag-asc1(ip))+as12(ip)*((8.5-xmag)**2)
     &         + (as3(ip)+ as13(ip)*(xmag-asc1(ip)))*alog(r)
          endif
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
          if(.not.rev)return      
c Subsequent AS97 calculations are for reverse-slip sources. 
c There could be a HW term
c for normal slip but this is not considered in USGS 2002 PSHA. And is not
c considered here for consistency.
        if(xmag.le.5.8) then
          fltfac= as5(ip)
        elseif(xmag.lt.asc1(ip))then
           fltfac= as5(ip)
     &        +(as6(ip)-as5(ip))*(xmag-5.8)/(asc1(ip)-5.8)
        else
           fltfac= as6(ip)
        endif
          gnd= gnd+ fltfac
        if(.not.hw)return
c------- for hanging wall sites some additional increments of ground motion
           if(xmag.lt.6.5) then
             fac1= xmag-5.5
           else
              fac1= 1.
           endif
          if(dist0.le.4..or.dist0.gt.25.)return
          if(dist0.lt.8.)then
               fac2= as9(ip)*(dist0-4.)*0.25
          elseif(dist0.lt.18.)then
              fac2= as9(ip)
          else
               fac2= as9(ip)*(1.-(dist0-18.)/7.)
          endif
             gnd= gnd+ fac1*fac2
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
c      write(14,*) period,dist0,xmag,exp(gnd)
      return
      end

ccccccccccccccccccccccccccccc
      subroutine getCamp2003(iper,ip,xmag,dist,rjb,gndout,sigmaf)
c  Campbell and Bozorgnia (2003). modifying for nga style. Input rock coeffs.
c Later, add soil coeffs?
c Hangingwall indicated by relation between dist and rjb (both input)
c Coeffs written but not checked july 27 2006. SH. Checked oct 2007. 
        parameter (sqrt2=1.414213562)
      common/mech/ss,rev,normal,obl
      logical ss,rev,normal,obl      
      common/dipinf/dipang,cosDELTA,cdipsq,cyhwfac,cbhwfac
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

      real, dimension(7):: c1,c2,c3,c4,c5,c6,c7,c8,c9,cmagsig
      real, dimension(7):: csigma1,csigmacoef,csigma2,c10,c11,csite,c15
      real mfac
      c1= (/-4.033,-2.771,-3.867,-2.661,-2.999,-3.556,-4.311/)
      c2= (/0.812,0.812,0.812, 0.812,0.812,0.812,0.812/)
      c3= (/0.036,0.030,-.101,0.060,.007,-0.035,-0.180/)
      c4= (/ -1.061,-1.153,-.964,-1.308,-1.080,-0.964,-0.964/)
      c5= (/0.041,0.098,0.019,0.166,0.059,0.023,0.019/)
      c6= (/-0.005,-0.014,0.,-0.009,-0.007,-0.002,0./)
      c7= (/-.018,-0.038,0.,-0.068,-0.022,-0.004,0./)
      c8= (/0.766,0.704,0.842,0.621,0.752,0.842,0.842/)
      c9= (/0.034,0.026,-.105,0.046,0.007,-0.036,-0.187/)
      c10= (/0.343,0.296,0.329,0.224,0.359,0.406,0.060/)
      c11= (/0.351,0.342,0.338,0.313,0.385,0.479,0.064/)
      c15= (/.370,.370,0.281,0.370,0.370,0.370,0.160/)      !HW effect.
      csite= (/-0.289,-0.331,-0.607,-0.299,-0.453,-0.528,-0.649/)      !is this a vs30-dependent term?
      csigma1= (/0.920,0.981,1.021,0.958,0.984,0.990,1.021/)
      csigmacoef= (/0.07, 0.07,0.07,0.07, 0.07,0.07,0.07/)
      csigma2= (/.402,0.463,0.503,0.44,0.466,0.472,0.503/)
      cmagsig= (/ 7.4,7.4,7.4,7.4,7.4,7.4,7.4/)
      mfac=(8.5-xmag)*(8.5-xmag)
c
c        if(icode(ip,ia).eq.1) xmag= 2.715 -0.277*xmag0+0.127*xmag0*xmag0
      gnd= c1(ip) + c2(ip)*xmag + c3(ip)*mfac
      arg0= c8(ip)*xmag + c9(ip)*mfac
      arg= dist*dist + ((c5(ip) + 0.5*c6(ip) +0.5*c7(ip))*exp(arg0))**2
      arg= sqrt(arg)
      gnd = gnd + c4(ip)*alog(arg) + csite(ip)
      if(.not.rev)goto 2
c----  for thrust faults and reverse faults, exta kick. Not for normal
c though.
      if(cosDELTA.lt.0.76)then
      gnd= gnd+ c11(ip)
c--- Another coef for steeper dipping reverse faults
      else
      gnd= gnd+ c10(ip)
      endif
      if(dist.le.30.0 .and. rjb.lt.5.) then
      fhw= (5.-rjb)/5.
      if(xmag.le.6.5) fhw= fhw*(xmag-5.5)
      if(dist.lt.8.)then
       fhw= fhw*c15(ip)*dist/8.
      else
       fhw= fhw*c15(ip)
       endif
c---- hanging wall term for thrust faults or more steeply dipping reverse faults.
      gnd= gnd + fhw
      endif
ccc   sigma C&B with mag-dependence
2          if(xmag.lt.cmagsig(ip)) csigma= csigma1(ip)- 
     &     csigmacoef(ip)*xmag
          if(xmag.ge.cmagsig(ip)) csigma= csigma2(ip)
      sigmaf= 1./csigma/sqrt2
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
      return
      end

cccccccccccccc
      subroutine getCampCEUS(iper,ip,ir,xmag,dist,gndout,csigma,sigmaf)
        parameter (sqrt2=1.414213562,alg70=4.2484952,alg130=4.8675345)
        parameter (np=7)
      common/mech/ss,rev,normal,obl
      common/ceus_sig/lceus_sigma,ceus_sigma
      logical lceus_sigma,ss,rev,normal,obl      
c      iper = period index in iatten(*,*)
c      ip = period index in coeff arrays below.
c  Campbell CEUS BSSA June 2003. BC or A rock. 7 periods.
c Campbell has exceptionally low sigma for higher M events (NMSZ M7.7 for example) 
c testing june 23 2006. SH. No clamp during
c initial effort. 6-23-06. Clipping applied for range of T, 0.02 to 0.5 s.
c input xmag= moment mag
c input dist = rrup (km) appropriate for Camp Ceus.
c ir=1 BC or firm rock
c ir=2 A or hard rock. Only difference is in constant term (check this)
c output gnd = log(median SA) (geometric mean), sigmaf= sigmafactor, and csigma
c
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8),sp

      real,dimension(np):: c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,
     1 clamp,c1h,csigma1,csigmacoef,csigma2,cmagsig,perx
       perx = (/0.,0.2,1.0,0.1,0.3,0.5,2.0/)
      c1= (/0.4492,.1325,-.3177,.4064,-0.1483,-0.1333,-1.2483/)
        c1h= (/0.0305,-0.4328,-.6104,-0.1475,-.6906,-.5907,-1.4306/)
      c2= (/.633,.617,.451,.613,0.609,0.534,0.459/)
      c3= (/-.0427,-.0586,-.2090,-.0353,-0.0786,-0.1379,-0.2552/)
      c4= (/-1.591,-1.32,-1.158,-1.369,-1.28,-1.216,-1.124/)
      c5= (/.683,.399,.299,0.484,0.349,0.318,.310/)
      c6= (/.416,.493,.503,0.467,0.502,0.503,.499/)
      c7= (/1.140,1.25,1.067,1.096,1.241,1.116,1.015/)
      c8= (/-.873,-.928,-.482,-1.284,-.753,-0.606,-.417/)
      c9= (/-.00428,-.0046,-.00255,-.00454,-.00414,-.00341,-.00187/)
      c10= (/.000483,.000337,.000141,.00046,.000263,.000194,.000103/)
      csigma1= (/1.030,1.077,1.110,1.059,1.081,1.098,1.093 /)
      csigmacoef= (/-.0860,-.0838,-.0793,-.0838,-0.0838,-0.0824,-.0758/)
      csigma2= (/0.414,.478,.543,0.460,0.482,0.508,0.551/)
      cmagsig= (/7.16,7.16,7.16,7.16,7.16,7.16,7.16/)
      sp = perx(ip).gt.0.02.and.perx(ip).lt.0.5
      if(ir.eq.1)then
      gnd=c1(ip)
      else
      gnd= c1h(ip)
      endif
      gnd= gnd + c2(ip)*xmag + c3(ip)*(8.5-xmag)*(8.5-xmag)
      arg= dist*dist + (c5(ip)*exp(c6(ip)*xmag))**2
      arg= sqrt(arg)
      fac=0.
      if(dist.gt.70.) fac= c7(ip)*(alog(dist)- alg70)
      if(dist.gt.130.) fac= fac+ c8(ip)*(alog(dist)-alg130)
      gnd = gnd + c4(ip)*alog(arg) + fac +(c9(ip)+c10(ip)*xmag)*dist
c      write(12,*) period,dist0,xmag,exp(gnd)
c--- following is for clipping 
         if(ip.eq.1) then
          gnd=min(gnd,0.405)
         elseif(sp)then
           gnd=min(gnd,1.099)
          endif
          if(lceus_sigma)then
          csigma = sigma_ceus      !special study 3/2011
          elseif(xmag.lt.cmagsig(ip))then
           csigma= csigma1(ip)+ csigmacoef(ip)*xmag
          else
           csigma= csigma2(ip)
          endif
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
      sigmaf= 1./csigma/sqrt2
c apply clamp in main if needed.
      return
      end subroutine getCampCEUS

ccccccc
      subroutine getBJF97(ip,iq,xmag,dist0,vs30,gndout,sigmaf)
c adapted to nga 7 periods. any vs30.
c dist0 is rjb
c returns gnd, sigmaf
         parameter (sqrt2=1.414213562,np=7)
      common/mech/ss,rev,normal,obl
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

      logical ss,rev,normal,obl      
        real,dimension(np):: b1ss,b2,b3,b4,b5,hsq,sigma,
     1 b1rv,b1all,bv,va
             real,dimension(8):: perx
c      !a reference set including pgv which wasn't used in 2002
c array constructors
      b1all= (/-0.242,1.089,-1.08,1.059,0.70,-.25,-1.743/)
        bv= (/-.371,-.292,-.698,-.212,-.401,-.553,-.655/)
        va= (/1396.,2118.,1406.,1112.,2133.,1816.,1795./)
      b1ss= (/-.313,0.999,-1.133,1.006,0.598,-.268,-1.699/)
      b1rv= (/-.117,1.17,-1.009,1.087,.803,.087,-1.801/)
      b1all= (/-.242,1.089,-1.08,1.059,0.70,-0.025,-1.743/)
      b2= (/0.527,0.711,1.036,0.753,0.769,0.884,1.085/)
      b3= (/0.,-0.207,-0.032,-.226,-.161,-.09,-.085/)
      b4= (/0.,0.,0.,0.,0.,0.,0./)
      b5= (/-0.778,-0.924,-0.798 ,-.934,-.893,-.846,-.812/)   
c      h= (/5.57,7.02,2.90,6.27,5.94,4.13,5.85/)
      hsq= (/31.0249,49.2804,8.41,39.3129,35.2836,17.0569,34.2225/)
      sigma= (/0.520,0.502,0.613,0.479,0.522,0.556,0.672/)
       perx= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)
c what about near-source dispersion?
      sig = sigma(iq)
          sigmaf= 1./sig/sqrt2
c mech-dependent and site Vs30 dependency
      if(ss)then
      gnd0=b1ss(iq)
      elseif(rev)then
      gnd0 =b1rv(iq)
      else
      gnd0=b1all(iq)
      endif
      gnd0=gnd0+bv(iq)*alog(vs30/va(iq))
          dist= sqrt(dist0*dist0+hsq(iq))
          gnd= gnd0 +b2(iq)*(xmag-6.)+b3(iq)*(xmag-6.)**2
          gnd= gnd + b4(iq)*dist + b5(iq)*alog(dist)
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
      return
      end subroutine getBJF97

cccccccc
      subroutine getMota(ip,iq,xmag,dist0,vs30,gndout,sigmaf)
c input xmag = moment mag
c      dist0 = rrup (check with Chuck) (km)
c      vs30  = top30 m Vs. for site amp term. (ditto)
c--  for Motazedian and Atkinson 2003 Puerto Rico relations
c adaptations: predefine some frequently used scalars. put in bjf97 siteamp
       parameter (np=7,alg75=1.8750613,sqrt2=1.41421356)
       parameter (vref=760.0,aln10=2.30258509,alg=2.99122608)
      real bv(np),va(np)
c alg = log10(980)      
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

      real deltasq,fac,perx(8)
      real c1(np),c2(np),c3(np),c4(np),sigma(np)
c array constructors 10/2006. SH
          perx= (/0.,0.2,1.0,0.1,0.3,0.5,2.0,-1./)      !-1 shall be reserved for pgv
        bv= (/-.371,-.292,-.698,-.212,-.401,-.553,-.655/)
        va= (/1396.,2118.,1406.,1112.,2133.,1816.,1795./)
      c1= (/3.87,4.33,3.40,3.,3.,3.,2.86/)
      c2= (/0.39062,0.38815,0.64818,1.,1.,1.,0.77055/)
      c3= (/-0.11289,-0.13977,-0.15222,0.,0.,0.,-0.11963/)
      c4= (/-0.00213,-0.00189,-0.00091,0.,0.,0.,-0.00082/)
      sigma= (/0.28,0.28,0.28,0.28,0.28,0.28,0.28/)
c      h= (/5.,5.,5.,5.,5.,5.,5./)
c site amp first, convert to base 10 for compatibility
c first try at siteamp is the BJF 97 version. no site nonlinearity. Needs work for better compatibility w/nga
      gnd0 = bv(iq)*alog(vs30/vref)/aln10
      sig = sigma(iq)* aln10 !second factor added Mar 11 2011.
          sigmaf= 1./sig/sqrt2
      period = perx(iq)
c gnd0 the site term is added to the magnitude-dependent terms
        deltasq = (-7.333 + 2.333*xmag)**2
        dist= sqrt(dist0*dist0 + deltasq)
        gnd=gnd0 +c1(iq)+c2(iq)*(xmag-6.)+c3(iq)*((xmag-6.)**2)+ c4(iq)*dist
        if(dist.le.75.) then
        fac= (-1.88+0.14*xmag)*alog10(dist)
        elseif(dist.le.100.)then
        fac= (-1.88+0.14*xmag)*alg75
        else 
        fac= (-1.88+0.14*xmag)*alg75 -0.5*alog10(dist/100.)
        endif
c alg serves to convert cm/s/s to g
        gnd= gnd + fac - alg
c base10 to base e
        gnd= gnd*aln10
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
      return
      end subroutine getMota
c

      subroutine getIdriss(iper,ip,xmag,rjb,vs30,gndout,sigmaf)
c Oct 2005: for pga only.
c ip is period index in calling program. 
c iq is period index in this subroutine. But there is no need for period
c subscripting because only pga is available.
      parameter (pi=3.14159265,sqrt2=1.414213562,vref=760.)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

      common/mech/ss,rev,normal,obl
      common/dipinf/dipang,cosDELTA,cdipsq,cyhwfac,cbhwfac
c----  This version assumes v30 in neighborhood of 760 m/sec (no Vs30 dependency)
c----  uses ln coefficients
      real magmin,dmag,xmag
      real a1/2.14/,a2/0.134/,b1/2.8/,b2/-0.197/
      real phi/0.08/,sigma/0.68/
      real a16/6.052/,a26/-.473/,b16/3.256/,b26/-0.273/
      logical ss,rev,normal,obl
      if(ip.ne.1)stop'getIdriss: pga only.'
c coeffs. from oct 5 2005 powerpoint progress report
      sig=sigma
          sigmaf= 1./sig/sqrt2
c-- 
      dist0= rjb
          dist= dist0
          if(xmag.lt.6.0) then
          gnd= a1+a2*xmag-(b1+b2*xmag)*alog(dist+10.)
        else
          gnd= a16+a26*xmag- (b16+b26*xmag)*alog(dist+10.)
          endif
          if(rev)then
          gnd = gnd+ phi
          endif
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
      return
      end subroutine getIdriss

      subroutine getBooreNGA07m(iper,ip,xmag,rjb,vs30,gndout,sigmaf)
c
c In: 
c      iper = period index in the epistemic gm uncert named common.
c      ip = period index in per() array below,
c      xmag=moment mag
c      rjb= Joyner Boore distance (km)
c      vs30= top30 m Vs, (m/s)
c Out: 
c      gndout = ln(median), sigmaf = 1/sigma/sqrt2
c
c April 2007 update of Boore-Atkinson documentation.  Includes 23 periods
c up to 10-s T.
c----  enter vs30 (need not be 760 m/sec)
c *** gnd is sensitive to mech type, based on logical variable ss, rev, normal
c----  non-lin soil terms.  
c --- returns log(median) and sigmaf, or "sigma-factor"=1/sigma/sqrt(2)
c --- testing mar 07. Several coeffs have changed meaning compared to earlier.
      parameter (pi=3.14159265,sqrt2=1.414213562,vref=760.,np=23)
      parameter (dx=1.098612289,dxsq=1.206948961,dxcube=1.325968960,plfac=-0.510825624)
c plfac = ln(pga_low/0.1)      This never changes so shouldnt be calculated.
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)
       common/fix_sigma/fix_sigma,sigma_fx
       common/sdi/sdi,dy_sdi,fac_sde
       real, dimension(8) :: fac_sde

c new feb 19 2013 sdi computations
      common/mech/ss,rev,normal,obl
      logical fix_sigma,sdi,ss,rev,normal,obl      
      real per(np),e1(np),e2(np),e3(np),e4(np),e5(np),e6(np),e7(np),e8(np)
     1 ,mh(np),c1(np),c2(np),c3(np),c4(np),mref(np),rref(np),h(np),
     2blin(np),b1(np),b2(np),v1(np),v2(np),
     2 sig1(np),sig2u(np),sigtu(np),sig2m(np),sigtm(np)
c Boore & Atkinson consider 23 periods in April 2007 rev. 
c Smooths a kink pp 9-11 of boore sept report.
c
c The e_jnl coeffs prior to Mar 20 2008, used for some SoilC and SoilD analysis for TS3 subcomm.
c      real e1nl/-0.03279/,e2nl/-0.03279/,e3nl/-0.03279/,e4nl/-0.03279/
c      real e5nl/0.29795/,e6nl/-0.20341/,e7nl/0.0/
c      real c1nl/-0.55/,c2nl/0./,c3nl/-0.01151/,hnl/1.35/,b1nl/0./
c BA now suggest using the final PGA coefficients for PGAnl in a followup article to the Mar 08 Eq
c Spectra paper. The above commented-out coeffs are the original "nl" coefficients, and the
c below coefficients are the pga coeffs which correspond to element 2 of the below arrays. 
c Update of Mar 20 2008.
      real e1nl/-0.53804/,e2nl/-0.50350/,e3nl/-0.75472/,e4nl/-0.50970/
      real e5nl/0.28805/,e6nl/-0.10164/,e7nl/0.0/
      real c1nl/-0.66050/,c2nl/0.11970/,c3nl/-0.011510/,hnl/1.35/,b1nl/0./
      real b2nl/0./,pga_low/0.06/,mhnl/6.75/,mrefnl/4.5/,rrefnl/1.0/
c a1,a2 are scalars here. Boore has 'em as vectors in the spreadsheet, but
c unchanging with spectral period.
      real  a1/ 0.030/,sigma_fx,sdisd
      real  a2/ 0.090/,a2fac/0.405465108/
      integer i
c dx = ln(a2/a1), made a param. used in a smoothed nonlin calculation. Mod of sept 2006.
c a2fac=ln(a2/pga_low)
c e2,e3,e4 are the mech-dependent set. e1 is a mech-unspecified value. 
c Mech-dep. coeffs are now used for pganl as well
c   as all other periods, mar 20, 2008. From BA pdf file of Mar 19, 2008. SH.
c 
c array constructors for f95 Linux
c from ba_02apr07_usnr.txt or xlx coef file.
c e1 unspecified; e2 ss; e3 normal; e4 thrust/reverse. Table 2 BA Sept 07 doc.
c
      per= (/-1.000, 0.000, 0.010, 0.020, 0.030, 0.050, 0.075, 0.100,
     + 0.150, 0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000,
     + 1.500, 2.000, 3.000, 4.000, 5.000, 7.500,10.0/)
       e1= (/ 5.00121,-0.53804,-0.52883,-0.52192,-0.45285,-0.28476,
     1  0.00767, 0.20109, 0.46128, 0.57180, 0.51884, 0.43825, 0.39220, 0.18957,-0.21338,
     1 -0.46896,-0.86271,-1.22652,-1.82979,-2.24656,-1.28408,-1.43145,-2.15446/)
       e2= (/ 5.04727,-0.50350,-0.49429,-0.48508,-0.41831,-0.25022,
     1  0.04912, 0.23102, 0.48661, 0.59253, 0.53496, 0.44516, 0.40602, 0.19878,-0.19496,
     1 -0.43443,-0.79593,-1.15514,-1.74690,-2.15906,-1.21270,-1.31632,-2.16137/)
c Editorial comment Oct 4 2007:
c I (SH) used the e3(10s)=e3(7.5)+e2(10s)-e2(7.5s) for 10s normal, because BA
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
       sig1= (/0.500,0.502,0.502,0.502,0.507,0.516,
     1 0.513,0.520,0.518,0.523,0.527,0.546,0.541,0.555,0.571,
     1 0.573,0.566,0.580,0.566,0.583,0.601,0.626,0.645/)
      sig2u= (/0.286,0.265,0.267,0.267,0.276,0.286,
     1 0.322,0.313,0.288,0.283,0.267,0.272,0.267,0.265,0.311,
     1 0.318,0.382,0.398,0.410,0.394,0.414,0.465,0.355/)
      sigtu= (/0.576,0.566,0.569,0.569,0.578,0.589,
     1 0.606,0.608,0.592,0.596,0.592,0.608,0.603,0.615,0.649,
     1 0.654,0.684,0.702,0.700,0.702,0.730,0.781,0.735/)
      sig2m= (/0.256,0.260,0.262,0.262,0.274,0.286,
     1 0.320,0.318,0.290,0.288,0.267,0.269,0.267,0.265,0.299,
     1 0.302,0.373,0.389,0.401,0.385,0.437,0.477,0.477/)
      sigtm= (/0.560,0.564,0.566,0.566,0.576,0.589,
     1 0.606,0.608,0.594,0.596,0.592,0.608,0.603,0.615,0.645,
     1 0.647,0.679,0.700,0.695,0.698,0.744,0.787,0.801/)

c end of April 2007 coeff. updates.
c      
c----  Input vs30. Ref v30= 760 m/sec
c
c  ip index corresponds to period. (check per correspondence with main perx)
c--- 
      site=0.0      !site term, nonzero for nonrefernce vs30
c sigma for specific mech type. If unspecified, sigmaf is redefined.
      if(fix_sigma)then
      sig=sigma_fx
      else
      sig=sigtm(ip)
      endif
          sigmaf= 1./sig/sqrt2      ! used with erf( ) in rate/prob calcs.
          dist= sqrt(rjb*rjb+h(ip)**2)
          if(xmag.le.mh(ip)) then
          gnd= e5(ip)*(xmag-mh(ip))
     1  +e6(ip)*((xmag-mh(ip))**2)
          else
           gnd= e7(ip)*(xmag-mh(ip))
c      1+e8(ip)*((xmag-mh(ip))**2)      !commented out because e8 is zero
           endif
          gnd= gnd+ (c1(ip)+ c2(ip)*(xmag-mref(ip)))*alog(dist/rref(ip))
     1  +c3(ip)*(dist-rref(ip))      
c      2 +c4(ip)*(dist-rref(ip))*(xmag-mref(ip))             
c 2nd continuation line commented out:nuisance
c computation because c4 is zero. Could change in the future. Still 0 june 06.     
          if(rev) then
          gnd = gnd + e4(ip)      !e4 can be lower than e2, the ss-mech coeff.
          elseif(normal) then
          gnd = gnd + e3(ip)      !e3 is lower.
          elseif(ss) then
          gnd = gnd + e2(ip)
          else                  !unspecified slip case
          gnd = gnd + e1(ip)
          sig=sigtu(ip)
           sigmaf= 1./sig/sqrt2      ! used with erf( ) in rate/prob calcs.
          endif
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
       if(sdi)then
       sde=gndout(1)+fac_sde(iper)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)      !10 is an upper bound for rhat.
       gndout(1) = sdi_ratio(per(ip),xmag,rhat,sig,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
         if(l_gnd_ep(iper))then
         gndout(2)= gndout(1)+gnd_ep(ide,ime,iper)
         gndout(3)= gndout(1)-gnd_ep(ide,ime,iper)
         endif
         endif      !if(sdi)
          if(vs30.eq.vref)return
c Otherwise correct for non-reference site conditions.
c nonlinear effects, updated March 20 2008. Mech -dependent pga4nl
c  pga estimation for the nonlinearly responding earth
          distnl = sqrt(rjb**2+hnl*hnl)
c Initialize pganl with mech-dep. value. New Mar 20 2008: see BA article following EqSpectra one.
          if(rev) then
          pganl = e4nl      !e4nl slightly lower than e2nl, the ss-mech coeff.
          elseif(normal) then
          pganl = e3nl     !e3 is lower than e4.
          elseif(ss) then
          pganl =  e2nl
          else                  !unspecified slip case
          pganl = e1nl
          endif
c compute magnitude-dependent part of pganl  
c The slip-type dependency was not present in BA_NGA, oct 2007; but is present Mar 2008.
          if(xmag.le.mhnl)then
      pganl = pganl+e5nl*(xmag-mhnl)+e6nl*((xmag-mhnl)**2)
          else
c no increase in Near-Field PGA above mhnl (6.75?), because e7nl is 0.
      pganl = pganl + e7nl*(xmag-mhnl)
           endif
c This pganl is needed when vs30 is not equal to vref. Now compute the
c distance dependency of pganl and add it to current value (log domain). 
          pganl= pganl+ c1nl*alog(distnl/rrefnl)+c3nl*(distnl-rrefnl)
     1 + c2nl*(xmag-mrefnl)*alog(distnl/rrefnl)   
c Note***: c2nl is no longer zero. ***
          pganl=exp(pganl)      !units g
c          
        if(v1(ip).lt.vs30.and.vs30.le.v2(ip))then
        bnl=(b1(ip)-b2(ip))*
     1 alog(vs30/v2(ip))/alog(v1(ip)/v2(ip)) + b2(ip)
        elseif(v2(ip).lt.vs30.and.vs30.le.vref)then
        bnl=b2(ip)*alog(vs30/vref)/alog(v2(ip)/vref)
        elseif(vs30.le.v1(ip))then
        bnl=b1(ip)
        else
        bnl=0.0
        endif
        dy=bnl*a2fac
c Below: include site term. C.f., site with 760 m/s vref.
        site=blin(ip)*alog(vs30/vref)
        if(pganl.le.a1)then
        site=site+bnl*plfac
        pgafac=0.
        elseif(pganl.le.a2)then
c c and d from p 10 of boore sept report. Smoothing introduces extra calcs
c in the range a1 < pganl < a2. Otherwise nonlin term same as in june-july.
c many of these terms are fixed and are defined in or parameter statements
c Of course, if a1 and a2 change from their sept 06 values the parameters will
c also have to be redefined.
        c=(3.*dy-bnl*dx)/dxsq
        d=(bnl*dx-2.*dy)/dxcube
        pgafac=alog(pganl/a1)
        site=site+bnl*plfac + c*pgafac**2 + d*pgafac**3
        else
        site=site+bnl*alog(pganl/0.1)
        pgafac=0.
        endif
        gnd=gnd+site
c modify gndout if gnd is modified by nonlinear soil amp. But, how to modify gnd_ep?
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
c There may be a different sigma for nonlin soil. Not included.
c  There is a sigma associated with pganl. 
c
       if(sdi)then
       sde=gndout(1)+fac_sde(iper)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)
       gndout(1) = sdi_ratio(per(ip),xmag,rhat,sig,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
         if(l_gnd_ep(iper))then
         gndout(2)= gndout(1)+gnd_ep(ide,ime,iper)
         gndout(3)= gndout(1)-gnd_ep(ide,ime,iper)
         endif
         endif      !if(sdi)
      return
      end subroutine getBooreNGA07m

        subroutine CY2007H(ip,iprd, M, R_Rup, R_JB, R_x, V_S30, Z1, Z_TOR,
     1                 F_RV, F_NM,  gndout, sigmaf)

        parameter (mxnprd=106)
      parameter (pi=3.14159265,d2r=0.0174533,sqrt2=1.414213562)
c   Jan 29 2009: add PGV coeffs to CY, from MAR 2008 Earthquake Spectra
c Previously output psa=ln(median) and sigmaf, two scalar quantities needed for Pex calcs.
c Jan 24 2007: psa replaced by gndout, a 3-component vector
c gndout(1)=psa, gndout(2)=psa+gnd_ep, gndout(3)=psa-gnd_ep. Mod was inspired
c by need to increase variability not present in the 3 or 4 NGA relations
c From Nov 2007 Chiou & Youngs NGA subroutine "cy2007".
c  Steve Harmsen 10/24/2007. Fortran efficiencies added.
c Predictor variables
c DELTA = fault dip (degrees) no longer an argument
c dipang (radians) is available in common/dipinf/
c z1 = depth (m) where vs is 1.0 km/s. THis is now used Oct 2007 code.
c vs30 top 30 m avg vs30,
c R_x new signed distance, < 0 on footwall for big extension of fault
c R_rup and R_JB the usual distance metrics.
c iprd the period index in the 105 element arrays below
c M the source moment mag
c Z_TOR depth to top of rupture (km).
c 
c This version returns gndout and sigmaf. It does not return psa_ref.
c psa_ref is the median SA for the "reference" site class, might be useful.
c New, Feb 19 2013: if sdi is .true., returns inelastic SD and its associated
c  sigmaf using Tothing & Cornell module to compute. Units : cm
       common/fix_sigma/fix_sigma,sigma_fx
      common/dipinf/dipang,cosDELTA,cdipsq,cyhwfac,cbhwfac
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/cyinit/a,b,c,rkdepth
       common/sdi/sdi,dy_sdi,fac_sde
       real, dimension(8) :: fac_sde
      real gnd_ep(3,3,8),gndout(3),dy_sdi
      logical l_gnd_ep(8),fix_sigma,sdi

c Predictor variables
        real PERIOD, M, Width, R_rup, R_JB, R_x, V_S30, Z1, 
     1       Z_TOR, F_RV, F_NM,sigma_fx

c Model coefficients. Phi moved to cy2007i, oct 25 2007. SHarmsen.
        real, dimension(mxnprd):: prd,
     1       c1, c1a, c1b, c2,
     1       c3, cn,  cM,  c4,
     1       c4a,cRB, c5,  c6,
     1       cHM,c7,  c9, c9a,
     1       cgamma1, cgamma2, cgamma3,
     1       sigma_t,tau1,tau2,sigma1,sigma2,dlt1
c sigma_t is from the 2006 CY model and is retained here only to make comparisons with current estimate
c of this important quantity. Steve Harmsen.
      real cc,gamma, cosDELTA, cyhwfac, cbhwfac, psa,psa_ref,sdisd
      real dipang,cdipsq,sigmaf, NL,sdi_ratio,fac
      real r1,r2, r3, r4, fw, hw
      real, dimension(8) :: a,b,c,rkdepth
      integer iprd, i
c prd array not used below. Assume SA 0.01 s is PGA. -1 corresponds to PGV.
      prd= (/0.010,0.020,0.022,0.025,0.029,0.030,0.032,0.035,0.036,0.040,0.042,0.044,0.045,0.046,
     10.048,0.050,0.055,0.060,0.065,0.067,0.070,0.075,0.080,0.085,0.090,0.095,0.100,0.110,0.120,
     10.130,0.133,0.140,0.150,0.160,0.170,0.180,0.190,0.200,0.220,0.240,0.250,0.260,0.280,0.290,
     10.300,0.320,0.340,0.350,0.360,0.380,0.400,0.420,0.440,0.450,0.460,0.480,0.500,0.550,0.600,
     10.650,0.667,0.700,0.750,0.800,0.850,0.900,0.950,1.000,1.100,1.200,1.300,1.400,1.500,1.600,
     11.700,1.800,1.900,2.000,2.200,2.400,2.500,2.600,2.800,3.000,3.200,3.400,3.500,3.600,3.800,
     14.000,4.200,4.400,4.600,4.800,5.000,5.500,6.000,6.500,7.000,7.500,8.000,8.500,9.000,9.500,
     110.0,-1.0/)
       c1= (/-1.26869,-1.25148,-1.23809,-1.21643,-1.18418,-1.17437,
     1 -1.15447,-1.12331,-1.11187,-1.06711,-1.04307,-1.01882,-1.00656,-0.99433,-0.97022,
     1 -0.94639,-0.88949,-0.83610,-0.78729,-0.76939,-0.74395,-0.70507,-0.67083,-0.64093,
     1 -0.61510,-0.59310,-0.57470,-0.54528,-0.52757,-0.51999,-0.52032,-0.52270,-0.53087,
     1 -0.54469,-0.56296,-0.58466,-0.60898,-0.63524,-0.68925,-0.74654,-0.77656,-0.80683,
     1 -0.86661,-0.89720,-0.92776,-0.98769,-1.04691,-1.07638,-1.10551,-1.16228,-1.21760,
     1 -1.27136,-1.32348,-1.34896,-1.37399,-1.42262,-1.46945,-1.57899,-1.67858,-1.76913,
     1 -1.79793,-1.85160,-1.92784,-1.99884,-2.06549,-2.12847,-2.18828,-2.24533,-2.35211,
     1 -2.45195,-2.54744,-2.64012,-2.73065,-2.81887,-2.90440,-2.98685,-3.06587,-3.14125,
     1 -3.28154,-3.40957,-3.46975,-3.52774,-3.63778,-3.74126,-3.83896,-3.93140,-3.97588,
     1 -4.01921,-4.10246,-4.18136,-4.25603,-4.32676,-4.39387,-4.45775,-4.51873,-4.66058,
     1 -4.78938,-4.90790,-5.01833,-5.12243,-5.22154,-5.31657,-5.40817,-5.49675,-5.58722,
     + 2.2884/)      !final value for pgv see table 2 p 198
      c1a= (/ 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
     1  0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
     1  0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
     1  0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
     1  0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100, 0.100,
     1  0.09990, 0.09990, 0.09990, 0.09990, 0.09990, 0.09980, 0.09980, 0.09980, 0.09970,
     1  0.09960, 0.09950, 0.09950, 0.09940, 0.09930, 0.09910, 0.09860, 0.09780, 0.09680,
     1  0.09630, 0.09540, 0.09360, 0.09140, 0.08870, 0.08530, 0.08130, 0.07660, 0.06510,
     1  0.05120, 0.03550, 0.01880, 0.00220,-0.01350,-0.02750,-0.03990,-0.05040,-0.05910,
     1 -0.07220,-0.08080,-0.08400,-0.08660,-0.09050,-0.09310,-0.09490,-0.09620,-0.09670,
     1 -0.09710,-0.09780,-0.09820,-0.09860,-0.09890,-0.09910,-0.09930,-0.09940,-0.09960,
     1 -0.09980,-0.09980,-0.09990,-0.09990,-0.09990,-0.100,-0.100,-0.100,-0.10000,
     + 0.1094/)
      c1b= (/-0.25500,-0.25500,-0.25500,-0.25500,-0.25500,-0.25500,
     1 -0.25500,-0.25500,-0.25500,-0.25500,-0.25500,-0.25500,-0.25500,-0.25500,-0.25500,
     1 -0.25500,-0.25500,-0.25500,-0.25470,-0.25400,-0.25400,-0.25400,-0.25400,-0.25400,
     1 -0.25400,-0.25300,-0.25300,-0.25290,-0.25200,-0.25100,-0.25100,-0.25040,-0.25000,
     1 -0.24900,-0.24800,-0.24700,-0.24600,-0.24490,-0.24280,-0.24000,-0.23820,-0.23700,
     1 -0.23430,-0.23280,-0.23130,-0.22750,-0.22470,-0.22260,-0.22140,-0.21800,-0.21460,
     1 -0.21070,-0.20730,-0.20540,-0.20370,-0.20080,-0.19720,-0.18890,-0.18140,-0.17440,
     1 -0.17220,-0.16800,-0.16200,-0.15640,-0.15110,-0.14720,-0.14320,-0.14000,-0.13370,
     1 -0.12820,-0.12460,-0.12140,-0.11840,-0.11660,-0.11400,-0.11250,-0.11110,-0.11000,
     1 -0.10800,-0.10700,-0.10600,-0.10600,-0.10500,-0.10400,-0.10400,-0.10300,-0.10300,
     1 -0.10300,-0.10200,-0.10200,-0.10200,-0.10200,-0.10200,-0.10100,-0.10100,-0.10100,
     1 -0.10100,-0.10100,-0.10100,-0.10100,-0.100,-0.100,-0.100,-0.100,-0.10000,
     + -0.0626/)
        c2= (/1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,
     1 1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.060,1.06/)
        c3= (/3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,
     1 3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.450,3.45/)
       cn= (/ 2.99600, 3.29240, 3.35160, 3.42930, 3.50120, 3.51370,
     1  3.53310, 3.55120, 3.55490, 3.56300, 3.56270, 3.56070, 3.55940, 3.55740, 3.55260,
     1  3.54730, 3.53080, 3.51290, 3.49300, 3.48440, 3.47140, 3.44800, 3.42320, 3.39660,
     1  3.36890, 3.34070, 3.31200, 3.25490, 3.19920, 3.14470, 3.12880, 3.09320, 3.04360,
     1  2.99700, 2.95220, 2.91000, 2.87010, 2.83120, 2.75960, 2.69160, 2.65790, 2.62560,
     1  2.56370, 2.53340, 2.50480, 2.44870, 2.39700, 2.37190, 2.34800, 2.30340, 2.26110,
     1  2.22180, 2.18480, 2.16720, 2.15000, 2.11780, 2.08680, 2.01740, 1.95700, 1.90360,
     1  1.88680, 1.85540, 1.81190, 1.77270, 1.73680, 1.70390, 1.67470, 1.64800, 1.60460,
     1  1.57170, 1.54620, 1.52630, 1.51100, 1.49840, 1.48890, 1.48140, 1.47440, 1.46980,
     1  1.46250, 1.45800, 1.45620, 1.45600, 1.45500, 1.45570, 1.45650, 1.45810, 1.45940,
     1  1.46060, 1.46300, 1.46520, 1.46760, 1.47030, 1.47330, 1.47520, 1.47790, 1.48310,
     1  1.48780, 1.49230, 1.49550, 1.49750, 1.49900, 1.500, 1.50100, 1.50100, 1.50200,
     + 1.648/)
       cM= (/ 4.18400, 4.18790, 4.18280, 4.17340, 4.15930, 4.15560,
     1  4.14850, 4.13820, 4.13510, 4.12260, 4.11740, 4.11230, 4.11040, 4.10840, 4.10380,
     1  4.10110, 4.09400, 4.08920, 4.08670, 4.08600, 4.08600, 4.08600, 4.08730, 4.08990,
     1  4.09380, 4.09850, 4.10300, 4.11440, 4.12770, 4.14160, 4.14590, 4.15650, 4.17170,
     1  4.18710, 4.20230, 4.21720, 4.23230, 4.24760, 4.27590, 4.30420, 4.31840, 4.33200,
     1  4.35840, 4.37120, 4.38440, 4.40860, 4.43230, 4.44410, 4.45570, 4.47680, 4.49790,
     1  4.51720, 4.53610, 4.54520, 4.55450, 4.57120, 4.58810, 4.62730, 4.66320, 4.69590,
     1  4.70710, 4.72760, 4.75710, 4.78510, 4.81140, 4.83620, 4.85970, 4.88200, 4.92450,
     1  4.96410, 5.00130, 5.03670, 5.06970, 5.10190, 5.13250, 5.16230, 5.19050, 5.21730,
     1  5.26910, 5.31730, 5.33930, 5.36100, 5.40130, 5.43850, 5.47370, 5.50690, 5.52290,
     1  5.53820, 5.56870, 5.59770, 5.62520, 5.65180, 5.67760, 5.70270, 5.72760, 5.78550,
     1  5.84040, 5.89240, 5.94220, 5.98910, 6.03390, 6.07700, 6.11720, 6.15610, 6.19300,
     + 4.2979/)
        c4= (/-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,
     1 -2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1,-2.1/)
      c4a= (/-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,
     1 -0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5,-.5/)
      cRB= (/50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,
     1 50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50.0,50./)
       c5= (/ 6.16000, 6.15800, 6.15800, 6.15700, 6.15580, 6.15500,
     1  6.15450, 6.15300, 6.15300, 6.15080, 6.14970, 6.14870, 6.14770, 6.14700, 6.14590,
     1  6.14410, 6.14090, 6.13620, 6.13140, 6.12940, 6.12600, 6.12000, 6.11440, 6.10720,
     1  6.10070, 6.09290, 6.08500, 6.06830, 6.04940, 6.02960, 6.02370, 6.00870, 5.98710,
     1  5.96470, 5.94160, 5.91770, 5.89420, 5.86990, 5.82310, 5.77670, 5.75470, 5.73350,
     1  5.69170, 5.67190, 5.65270, 5.61630, 5.58320, 5.56810, 5.55280, 5.52520, 5.49970,
     1  5.47640, 5.45550, 5.44580, 5.43620, 5.41890, 5.40290, 5.36970, 5.34310, 5.32130,
     1  5.31490, 5.30450, 5.29000, 5.27880, 5.26920, 5.26070, 5.25370, 5.24800, 5.23870,
     1  5.23210, 5.22660, 5.22240, 5.21940, 5.21660, 5.21400, 5.21250, 5.21110, 5.20990,
     1  5.20800, 5.20600, 5.20600, 5.20500, 5.20430, 5.20400, 5.20300, 5.20300, 5.20300,
     1  5.20240, 5.20200, 5.20200, 5.20200, 5.20170, 5.20100, 5.20100, 5.20100, 5.20100,
     1  5.20100, 5.20100, 5.20100, 5.200, 5.200, 5.200, 5.200, 5.200, 5.20000,
     + 5.17/)
       c6= (/ 0.48930, 0.48920, 0.48920, 0.48910, 0.48910, 0.48900,
     1  0.48900, 0.48890, 0.48890, 0.48880, 0.48870, 0.48870, 0.48860, 0.48860, 0.48850,
     1  0.48840, 0.48830, 0.48800, 0.48780, 0.48760, 0.48750, 0.48720, 0.48690, 0.48650,
     1  0.48620, 0.48580, 0.48540, 0.48460, 0.48370, 0.48280, 0.48250, 0.48180, 0.48080,
     1  0.47970, 0.47870, 0.47760, 0.47650, 0.47550, 0.47350, 0.47150, 0.47060, 0.46980,
     1  0.46800, 0.46730, 0.46650, 0.46510, 0.46380, 0.46320, 0.46260, 0.46160, 0.46070,
     1  0.45980, 0.45910, 0.45870, 0.45830, 0.45780, 0.45710, 0.45600, 0.45500, 0.45420,
     1  0.45400, 0.45360, 0.45310, 0.45280, 0.45240, 0.45220, 0.45190, 0.45170, 0.45140,
     1  0.45110, 0.45100, 0.45080, 0.45070, 0.45060, 0.45050, 0.45040, 0.45040, 0.45040,
     1  0.45030, 0.45020, 0.45020, 0.45020, 0.45020, 0.45010, 0.45010, 0.45010, 0.45010,
     1  0.45010, 0.45010, 0.45010, 0.45010, 0.45010, 0.45010, 0.45000, 0.45000, 0.45000,
     1  0.45000, 0.45000, 0.45000, 0.45000, 0.45000, 0.45000, 0.45000, 0.45000, 0.45000,
     + 0.4407/)
      cHM= (/ 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,
     1  3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,3.0/)
       c7= (/ 0.05120, 0.05120, 0.05120, 0.05120, 0.05110, 0.05110,
     1  0.05110, 0.05100, 0.05090, 0.05080, 0.05070, 0.05060, 0.05060, 0.05050, 0.05050,
     1  0.05040, 0.05020, 0.05000, 0.04990, 0.04980, 0.04970, 0.04950, 0.04940, 0.04920,
     1  0.04910, 0.04900, 0.04890, 0.04860, 0.04840, 0.04820, 0.04820, 0.04810, 0.04790,
     1  0.04770, 0.04760, 0.04740, 0.04730, 0.04710, 0.04680, 0.04660, 0.04640, 0.04630,
     1  0.04600, 0.04590, 0.04580, 0.04550, 0.04530, 0.04520, 0.04500, 0.04480, 0.04450,
     1  0.04420, 0.04390, 0.04370, 0.04360, 0.04320, 0.04290, 0.04210, 0.04120, 0.04040,
     1  0.04010, 0.03950, 0.03870, 0.03790, 0.03720, 0.03640, 0.03570, 0.03500, 0.03360,
     1  0.03220, 0.03080, 0.02940, 0.02800, 0.02660, 0.02530, 0.02400, 0.02270, 0.02130,
     1  0.01880, 0.01650, 0.01540, 0.01440, 0.01240, 0.01060, 0.00900, 0.00760, 0.00690,
     1  0.00630, 0.00520, 0.00410, 0.00330, 0.00250, 0.00190, 0.00140, 0.00100, 0.00020,
     1  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.00000,
     1 0.0207/)
c No c7a (for aftershocks). not included.
       c9= (/ 0.79000, 0.81290, 0.81880, 0.82800, 0.84070, 0.84390,
     1  0.85020, 0.85940, 0.86240, 0.87400, 0.87950, 0.88480, 0.88740, 0.88990, 0.89490,
     1  0.89960, 0.91050, 0.92040, 0.92920, 0.93250, 0.93710, 0.94420, 0.95050, 0.95590,
     1  0.96060, 0.96450, 0.96770, 0.97180, 0.97330, 0.97250, 0.97190, 0.96990, 0.96600,
     1  0.96090, 0.95490, 0.94820, 0.94100, 0.93340, 0.91790, 0.90230, 0.89460, 0.88710,
     1  0.87260, 0.86570, 0.85900, 0.84620, 0.83420, 0.82840, 0.82280, 0.81210, 0.80190,
     1  0.79220, 0.78300, 0.77860, 0.77430, 0.76590, 0.75780, 0.73910, 0.72210, 0.70650,
     1  0.70150, 0.69220, 0.67880, 0.66620, 0.65400, 0.64230, 0.63080, 0.61960, 0.59750,
     1  0.57560, 0.55380, 0.53200, 0.51010, 0.48770, 0.46490, 0.44120, 0.41690, 0.39170,
     1  0.33900, 0.28330, 0.25460, 0.22620, 0.17220, 0.12440, 0.08460, 0.05380, 0.04200,
     1  0.03220, 0.01770, 0.00860, 0.00310, 0.00040, 0.000, 0.000, 0.000, 0.000,
     1  0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.00000,
     + 0.3079/)
      c9a= (/ 1.50050, 1.50276, 1.50351, 1.50456, 1.50667, 1.50712,
     1  1.50818, 1.51044, 1.51104, 1.51376, 1.51573, 1.51771, 1.51862, 1.51953, 1.52135,
     1  1.52303, 1.52974, 1.53587, 1.54280, 1.54635, 1.55162, 1.55971, 1.56784, 1.57917,
     1  1.59011, 1.60047, 1.61043, 1.63837, 1.66446, 1.69266, 1.70233, 1.72461, 1.75488,
     1  1.78479, 1.81939, 1.85280, 1.88476, 1.91573, 1.98060, 2.04153, 2.07094, 2.09824,
     1  2.15050, 2.17581, 2.20053, 2.24678, 2.28462, 2.30297, 2.32100, 2.35585, 2.38858,
     1  2.41259, 2.43562, 2.44685, 2.45788, 2.47936, 2.50002, 2.53375, 2.56459, 2.58933,
     1  2.59529, 2.60648, 2.62243, 2.63663, 2.64561, 2.65382, 2.66153, 2.66899, 2.67728,
     1  2.68505, 2.69096, 2.69474, 2.69851, 2.70148, 2.70337, 2.70527, 2.70689, 2.70851,
     1  2.71014, 2.71177, 2.71231, 2.71285, 2.71367, 2.71448, 2.71502, 2.71529, 2.71556,
     1  2.71584, 2.71611, 2.71638, 2.71665, 2.71665, 2.71692, 2.71692, 2.71720, 2.71747,
     1  2.71747, 2.71774, 2.71774, 2.71774, 2.71801, 2.71801, 2.71801, 2.71801, 2.71801,
     + 2.6690/)
      cgamma1= (/-0.008040,-0.008113,-0.008155,-0.008230,-0.008351,-0.008387,
     1 -0.008454,-0.008564,-0.008599,-0.008754,-0.008833,-0.008906,-0.008945,-0.008984,-0.009050,
     1 -0.009121,-0.009288,-0.009433,-0.009561,-0.009603,-0.009656,-0.009733,-0.009782,-0.009808,
     1 -0.009809,-0.009788,-0.009753,-0.009627,-0.009460,-0.009269,-0.009205,-0.009054,-0.008832,
     1 -0.008610,-0.008398,-0.008183,-0.007976,-0.007776,-0.007397,-0.007044,-0.006877,-0.006715,
     1 -0.006408,-0.006263,-0.006123,-0.005859,-0.005614,-0.005497,-0.005386,-0.005175,-0.004980,
     1 -0.004799,-0.004632,-0.004552,-0.004477,-0.004332,-0.004199,-0.003903,-0.003652,-0.003435,
     1 -0.003368,-0.003246,-0.003078,-0.002929,-0.002795,-0.002674,-0.002564,-0.002464,-0.002287,
     1 -0.002138,-0.002010,-0.001898,-0.001802,-0.001717,-0.001643,-0.001578,-0.001519,-0.001467,
     1 -0.001381,-0.001311,-0.001281,-0.001255,-0.001209,-0.001172,-0.001142,-0.001117,-0.001106,
     1 -0.001096,-0.001079,-0.001065,-0.001052,-0.001042,-0.001032,-0.001024,-0.001016,-0.000999,
     1 -0.000987,-0.000977,-0.000968,-0.000961,-0.000955,-0.000950,-0.000945,-0.000941,-0.000937,
     + -0.00275/)
      cgamma2= (/-0.007850,-0.007921,-0.007962,-0.008035,-0.008154,-0.008189,
     1 -0.008255,-0.008362,-0.008396,-0.008547,-0.008624,-0.008696,-0.008734,-0.008771,-0.008836,
     1 -0.008906,-0.009068,-0.009210,-0.009335,-0.009376,-0.009428,-0.009503,-0.009550,-0.009577,
     1 -0.009577,-0.009557,-0.009522,-0.009400,-0.009236,-0.009050,-0.008988,-0.008840,-0.008623,
     1 -0.008406,-0.008199,-0.007989,-0.007787,-0.007592,-0.007222,-0.006878,-0.006714,-0.006556,
     1 -0.006257,-0.006115,-0.005979,-0.005720,-0.005481,-0.005368,-0.005259,-0.005053,-0.004862,
     1 -0.004686,-0.004522,-0.004445,-0.004371,-0.004230,-0.004100,-0.003811,-0.003565,-0.003354,
     1 -0.003288,-0.003169,-0.003005,-0.002860,-0.002729,-0.002611,-0.002504,-0.002406,-0.002233,
     1 -0.002087,-0.001962,-0.001854,-0.001759,-0.001677,-0.001604,-0.001540,-0.001483,-0.001433,
     1 -0.001348,-0.001280,-0.001251,-0.001225,-0.001181,-0.001145,-0.001115,-0.001091,-0.001080,
     1 -0.001071,-0.001054,-0.001040,-0.001027,-0.001017,-0.001007,-0.000999,-0.000992,-0.000976,
     1 -0.000964,-0.000954,-0.000945,-0.000938,-0.000933,-0.000927,-0.000923,-0.000919,-0.000914,
     + -.00625/)
      cgamma3= (/ 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000,
     1  4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.0000, 4.000, 4.0/)
c phi values are in subroutine cy2007i, below.
c the below tau1,2, sigma1,2 dlta are new 5pm nov 2 2007.
      tau1= (/0.3437,0.3471,0.3505,0.3538,0.3571,0.3603,
     1 0.3633,0.3663,0.3691,0.3718,0.3744,0.3768,0.3791,0.3811,0.3831,
     1 0.3848,0.3863,0.3876,0.3877,0.3881,0.3883,0.3878,0.3872,0.3865,
     1 0.3856,0.3846,0.3835,0.3816,0.3795,0.3775,0.3761,0.3742,0.3719,
     1 0.3696,0.3672,0.3649,0.3626,0.3601,0.3572,0.3543,0.3522,0.3500,
     1 0.3474,0.3455,0.3438,0.3417,0.3398,0.3386,0.3375,0.3362,0.3351,
     1 0.3344,0.3339,0.3340,0.3344,0.3346,0.3353,0.3354,0.3360,0.3369,
     1 0.3390,0.3409,0.3429,0.3452,0.3478,0.3508,0.3541,0.3577,0.3608,
     1 0.3644,0.3682,0.3725,0.3769,0.3816,0.3865,0.3916,0.3968,0.4023,
     1 0.4085,0.4149,0.4212,0.4277,0.4341,0.4406,0.4470,0.4534,0.4598,
     1 0.4661,0.4723,0.4784,0.4845,0.4904,0.4962,0.5019,0.5074,0.5128,
     1 0.5181,0.5232,0.5281,0.5328,0.5374,0.5419,0.5461,0.5502,0.5542,
     + 0.2539/)
      tau2= (/0.2637,0.2671,0.2705,0.2738,0.2771,0.2803,
     1 0.2833,0.2863,0.2891,0.2918,0.2944,0.2968,0.2991,0.3011,0.3031,
     1 0.3048,0.3063,0.3076,0.3095,0.3106,0.3118,0.3129,0.3138,0.3145,
     1 0.3149,0.3151,0.3152,0.3154,0.3153,0.3151,0.3143,0.3135,0.3128,
     1 0.3120,0.3110,0.3100,0.3089,0.3076,0.3068,0.3060,0.3047,0.3034,
     1 0.3026,0.3015,0.3005,0.2999,0.2993,0.2988,0.2983,0.2983,0.2984,
     1 0.2987,0.2993,0.2999,0.3008,0.3021,0.3036,0.3060,0.3085,0.3113,
     1 0.3139,0.3169,0.3205,0.3243,0.3283,0.3326,0.3371,0.3419,0.3472,
     1 0.3527,0.3584,0.3643,0.3703,0.3765,0.3828,0.3892,0.3957,0.4023,
     1 0.4085,0.4149,0.4212,0.4277,0.4341,0.4406,0.4470,0.4534,0.4598,
     1 0.4661,0.4723,0.4784,0.4845,0.4904,0.4962,0.5019,0.5074,0.5128,
     1 0.5181,0.5232,0.5281,0.5328,0.5374,0.5419,0.5461,0.5502,0.5542,
     + 0.2381/)
      sigma1= (/0.4458,0.4458,0.4476,0.4500,0.4529,0.4535,
     1 0.4547,0.4564,0.4569,0.4589,0.4598,0.4607,0.4611,0.4615,0.4623,
     1 0.4630,0.4647,0.4663,0.4677,0.4682,0.4690,0.4702,0.4712,0.4722,
     1 0.4731,0.4740,0.4747,0.4761,0.4773,0.4782,0.4785,0.4791,0.4798,
     1 0.4803,0.4808,0.4811,0.4814,0.4816,0.4817,0.4816,0.4815,0.4813,
     1 0.4808,0.4805,0.4801,0.4794,0.4786,0.4781,0.4777,0.4768,0.4758,
     1 0.4748,0.4738,0.4734,0.4729,0.4719,0.4710,0.4688,0.4667,0.4650,
     1 0.4644,0.4634,0.4621,0.4610,0.4600,0.4592,0.4586,0.4581,0.4555,
     1 0.4535,0.4518,0.4505,0.4493,0.4484,0.4476,0.4469,0.4463,0.4459,
     1 0.4451,0.4444,0.4442,0.4440,0.4436,0.4433,0.4430,0.4428,0.4427,
     1 0.4426,0.4425,0.4424,0.4423,0.4422,0.4421,0.4420,0.4420,0.4418,
     1 0.4417,0.4417,0.4416,0.4416,0.4415,0.4415,0.4415,0.4415,0.4414,
     + 0.4496/)
      sigma2= (/0.3459,0.3459,0.3477,0.3502,0.3530,0.3537,
     1 0.3549,0.3566,0.3572,0.3592,0.3602,0.3611,0.3615,0.3619,0.3627,
     1 0.3635,0.3654,0.3670,0.3686,0.3692,0.3700,0.3713,0.3726,0.3738,
     1 0.3749,0.3759,0.3769,0.3787,0.3804,0.3819,0.3824,0.3834,0.3847,
     1 0.3859,0.3871,0.3882,0.3893,0.3902,0.3921,0.3938,0.3946,0.3953,
     1 0.3967,0.3974,0.3981,0.3993,0.4005,0.4010,0.4016,0.4026,0.4036,
     1 0.4046,0.4054,0.4059,0.4063,0.4071,0.4079,0.4098,0.4114,0.4130,
     1 0.4135,0.4144,0.4157,0.4170,0.4181,0.4192,0.4203,0.4213,0.4213,
     1 0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,
     1 0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,
     1 0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,
     1 0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,0.4213,
     +0.3554/)
c dlt1 is used for the case of INFERRED Vs30. dlt1 is sigma3 of Eq Spectra paper
        dlt1= (/0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,
     1 0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,
     1 0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,
     1 0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,0.8000,
     1 0.8000,0.8000,0.8000,0.8000,0.8000,0.7999,0.7999,0.7999,0.7998,
     1 0.7998,0.7997,0.7997,0.7996,0.7994,0.7993,0.7992,0.7990,0.7988,
     1 0.7983,0.7979,0.7976,0.7974,0.7970,0.7966,0.7940,0.7917,0.7884,
     1 0.7867,0.7836,0.7792,0.7747,0.7681,0.7619,0.7560,0.7504,0.7400,
     1 0.7304,0.7230,0.7182,0.7136,0.7097,0.7080,0.7064,0.7049,0.7035,
     1 0.7025,0.7017,0.7012,0.7011,0.7008,0.7006,0.7004,0.7003,0.7003,
     1 0.7002,0.7002,0.7001,0.7001,0.7001,0.7001,0.7001,0.7000,0.7000,
     1 0.7000,0.7000,0.7000,0.7000,0.7000,0.7000,0.7000,0.7000,0.7000,0.7/)
c sigma_t = 2006 not to be used with final CY model
              cc = c5(iprd)* cosh(c6(iprd) * max(M-cHM(iprd),0.))
        gamma = cgamma1(iprd) +
     1          cgamma2(iprd)/cosh(max(M-cgamma3(iprd),0.))
c        cosDELTA = cos(DELTA*d2r)

c Magnitude scaling
        r1 = c1(iprd) + c2(iprd) * (M-6.0) +
     1       (c2(iprd)-c3(iprd))/cn(iprd) *
     1             log(1.0 + exp(-cn(iprd)*(M-cM(iprd))))

c Near-field magnitude and distance scaling
        r2 = c4(iprd) * log(R_Rup + cc)

c Far-field distance
        r3 = (c4a(iprd)-c4(iprd))/2.0 *
     1            log( R_Rup*R_Rup+cRB(iprd)*cRB(iprd) ) +
     1       R_Rup * gamma

c Scaling with other source variables (F_RV, F_NM, and Z_TOR)
        r4 = c1a(iprd)*F_RV +
     1       c1b(iprd)*F_NM +
     1       c7(iprd)*(Z_TOR - 4.)

c HW effect. R_x < 0? A signed distance measure. R_x < 0 occurs at footwall locations
        if (R_x .lt. 0) then
          hw = 0.0
        else
          hw = c9(iprd) * tanh(R_x*cosDELTA**2/c9a(iprd)) *
     1        (1.0 - sqrt(R_JB**2+Z_TOR**2)/(R_Rup + 0.001))
        endif


        psa_ref = r1+r2+r3+r4+hw

c......
c Soil effect: linear response. The commented-out lines are done in CY2007I, below
c        a(ip) = phi1(iprd) * min(log(V_S30/1130.), 0.)

c Soil effect: nonlinear response. 1130-360= 770 
c        b(ip) = phi2(iprd) *
c     1(exp(phi3(iprd)*(min(V_S30,1130.)-360.))-exp(phi3(iprd)*(770.)))
c        c (ip)= phi4(iprd)

c Modificaiton to ln(Vs30) scaling: bedrock depth (Z1)
c NOTE: max(0,Z1-15) is capped at 300 to avoid overflow of function cosh
c        rkdepth(ip) = phi5(iprd) *
c     1        ( 1 - 1.0/cosh(phi6(iprd)*max(0.,Z1-phi7(iprd)))) +
c     1        phi8(iprd)/cosh(0.15*min(max(0., Z1-15.),300.))

c......
c Median PSA prediction for reference condition. Keep psa a logged quantity. SH.
        psa = psa_ref + 
     1 (a(ip) + b(ip) * log((exp(psa_ref)+c(ip))/c(ip))) + rkdepth(ip)
        psa_ref = exp(psa_ref)
c....... Aleatory variablility (to be provided soon)
c Tau
        tau = tau1(iprd) + (tau2(iprd)-tau1(iprd))/2. * (min(max(M,5.),7.)-5.)
c To compute NL, psa_ref is Not Logged
        NL = b(ip) * psa_ref/(psa_ref+c(ip))
c correction by B Chiou, 11/06/2007& 11/08
c      tau = tau * sqrt(1.0+ NL)
      tau = tau * (1.0 + NL)
c.......
c Sigma. This is the sigma associated with "Inferred" Vs30. As distinguished from "measured"
c Vs30. Note added Nov 15 2012. SHarmsen.
      if(fix_sigma)then 
      sig=sigma_fx
      else
        sigma_M = sigma1(iprd) +
     1        (sigma2(iprd)-sigma1(iprd))/2*(min(max(M,5.),7.)-5.)
c dlt1 below is an additional feature of the 5pm nov 2 model from Chiou email.
        sig = sigma_M * sqrt(dlt1(iprd)+ (1.0 + NL)**2)
      sig = sqrt(sig**2 + tau**2)
      endif
c 2007 sigma is returned. 
c   sigma is  heteroscedastic. lower sigma for soil sites and high ref. rock PGA.
       gndout(1) = psa
         if(l_gnd_ep(ip))then
         gndout(2)= psa+gnd_ep(ide,ime,ip)
         gndout(3)= psa-gnd_ep(ide,ime,ip)
         endif
       sigmaf= 1.0/sig/sqrt2
       if(sdi)then
       sde=gndout(1)+fac_sde(ip)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)
       gndout(1) = sdi_ratio(prd(iprd),M,rhat,sig,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
         if(l_gnd_ep(ip))then
         gndout(2)= gndout(1)+gnd_ep(ide,ime,ip)
         gndout(3)= gndout(1)-gnd_ep(ide,ime,ip)
         endif
         endif      !if(sdi)
      return
      end subroutine CY2007H

        subroutine CY2007I(ip,iprd, V_S30, Z1)                  
c this subroutine computes some soil terms (a,b,c,rkdepth) for a fixed
c vs30 and z1 model. These are stored in array at element ip. This initialization
c is done to speed up runs. 
c Add PGV coeffs Jan 29 2009. PGV output units: cm/s.
c Input:
c   ip, period index in global psha run. IPMAX is 8.
c   iprd, period index of the CY coefficients associated with period(ip)
c    V_S30 = vs in top 30 m (m/s)
c    Z1 = depth (m) where Vs is >= 1 km/s. THis quantity should be fixed for all
c    receivers in the run (like the standard Vs30 maps we have produced). See phi7.
c
c 
        parameter (mxnprd=106)
      parameter (pi=3.14159265,d2r=0.0174533,sqrt2=1.414213562)
      common/cyinit/a,b,c,rkdepth
        real V_S30, Z1
      real, dimension(8) :: a,b,c,rkdepth
      real, dimension(mxnprd) :: phi1,phi2,phi3,phi4,phi5,phi6,phi7,phi8
      phi1= (/-.4417,-.4340,-.4313,-.4267,-.4196,-.4177,
     1 -.4139,-.4082,-.4064,-.4000,-.3973,-.3949,-.3939,-.3930,-.3914,
     1 -.3903,-.3892,-.3903,-.3934,-.3951,-.3981,-.4040,-.4108,-.4182,
     1 -.4261,-.4341,-.4423,-.4585,-.4743,-.4892,-.4935,-.5032,-.5162,
     1 -.5283,-.5396,-.5502,-.5602,-.5697,-.5873,-.6034,-.6109,-.6182,
     1 -.6319,-.6383,-.6444,-.6559,-.6665,-.6715,-.6762,-.6850,-.6931,
     1 -.7005,-.7072,-.7104,-.7135,-.7193,-.7246,-.7365,-.7468,-.7557,
     1 -.7585,-.7636,-.7708,-.7773,-.7833,-.7888,-.7941,-.7990,-.8082,
     1 -.8165,-.8243,-.8315,-.8382,-.8445,-.8504,-.8560,-.8613,-.8663,
     1 -.8755,-.8836,-.8874,-.8909,-.8974,-.9032,-.9083,-.9130,-.9151,
     1 -.9170,-.9205,-.9231,-.9249,-.9257,-.9255,-.9243,-.9222,-.9129,
     1 -.8982,-.8791,-.8572,-.8346,-.8126,-.7914,-.7711,-.7517,-.7332,-0.7861/)
      phi2= (/-.1417,-.1364,-.1361,-.1365,-.1392,-.1403,
     1 -.1430,-.1482,-.1502,-.1591,-.1641,-.1694,-.1721,-.1748,-.1804,
     1 -.1862,-.2008,-.2153,-.2291,-.2344,-.2420,-.2538,-.2644,-.2739,
     1 -.2819,-.2887,-.2943,-.3025,-.3077,-.3106,-.3111,-.3118,-.3113,
     1 -.3093,-.3062,-.3022,-.2976,-.2927,-.2823,-.2716,-.2662,-.2609,
     1 -.2505,-.2455,-.2405,-.2310,-.2220,-.2177,-.2135,-.2053,-.1975,
     1 -.1901,-.1830,-.1795,-.1762,-.1696,-.1633,-.1487,-.1353,-.1232,
     1 -.1194,-.1124,-.1028,-.0943,-.0869,-.0805,-.0748,-.0699,-.0617,
     1 -.0552,-.0501,-.0459,-.0425,-.0395,-.0369,-.0346,-.0323,-.0302,
     1 -.0262,-.0225,-.0207,-.0190,-.0159,-.0129,-.0102,-.0077,-.0066,
     1 -.0055,-.0036,-.0016,0.00,0.00,0.00,0.00,0.00,0.00,
     1 0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.0000,-0.0699/)
         phi3= (/-0.007010,-0.007279,-0.007301,-0.007364,-0.007378,-0.007354,
     1 -0.007281,-0.007162,-0.007129,-0.006977,-0.006878,-0.006765,-0.006710,-0.006656,-0.006556,
     1 -0.006467,-0.006279,-0.006117,-0.005970,-0.005914,-0.005835,-0.005734,-0.005670,-0.005632,
     1 -0.005607,-0.005597,-0.005604,-0.005644,-0.005696,-0.005744,-0.005758,-0.005794,-0.005845,
     1 -0.005901,-0.005959,-0.006019,-0.006080,-0.006141,-0.006262,-0.006381,-0.006439,-0.006495,
     1 -0.006604,-0.006655,-0.006704,-0.006795,-0.006882,-0.006923,-0.006965,-0.007047,-0.007125,
     1 -0.007194,-0.007259,-0.007290,-0.007320,-0.007378,-0.007435,-0.007579,-0.007720,-0.007863,
     1 -0.007911,-0.008001,-0.008120,-0.008223,-0.008313,-0.008381,-0.008423,-0.008444,-0.008500,
     1 -0.008478,-0.008307,-0.008042,-0.007707,-0.007317,-0.006862,-0.006265,-0.005541,-0.004792,
     1 -0.003555,-0.002764,-0.002497,-0.002292,-0.002007,-0.001828,-0.001713,-0.001636,-0.001608,
     1 -0.001585,-0.001549,-0.001523,-0.001501,-0.001483,-0.001467,-0.001453,-0.001440,-0.001416,
     1 -0.001397,-0.001384,-0.001375,-0.001369,-0.001364,-0.001362,-0.001360,-0.001360,-0.001361,
     + -.008444/)
         phi4= (/ 0.102151, 0.108360, 0.110372, 0.113710, 0.118600, 0.119888,
     1  0.122493, 0.126540, 0.127926, 0.133641, 0.136572, 0.139596, 0.141112, 0.142659, 0.145774,
     1  0.148927, 0.157001, 0.165249, 0.173635, 0.177001, 0.182082, 0.190596, 0.199129, 0.207505,
     1  0.215628, 0.223398, 0.230662, 0.243315, 0.253169, 0.260175, 0.261767, 0.264504, 0.266468,
     1  0.266468, 0.265060, 0.262501, 0.259163, 0.255253, 0.246252, 0.236525, 0.231541, 0.226570,
     1  0.216796, 0.211993, 0.207277, 0.198077, 0.189304, 0.185074, 0.180920, 0.172976, 0.165464,
     1  0.158358, 0.151662, 0.148466, 0.145352, 0.139415, 0.133828, 0.121226, 0.110339, 0.100842,
     1  0.097891, 0.092504, 0.085153, 0.078622, 0.072788, 0.067563, 0.062850, 0.058595, 0.051206,
     1  0.045054, 0.039879, 0.035504, 0.031787, 0.028613, 0.025890, 0.023537, 0.021496, 0.019716,
     1  0.016771, 0.014434, 0.013436, 0.012534, 0.010962, 0.009643, 0.008521, 0.007561, 0.007130,
     1  0.006730, 0.006008, 0.005379, 0.004830, 0.004349, 0.003925, 0.003553, 0.003223, 0.002551,
     1  0.002047, 0.001662, 0.001366, 0.001134, 0.000952, 0.000806, 0.000689, 0.000593, 0.000515,
     + 5.410/)
         phi5= (/ 0.228900, 0.228900, 0.228900, 0.228900, 0.228900, 0.228900,
     1  0.228900, 0.228900, 0.228900, 0.228900, 0.228900, 0.229000, 0.229000, 0.229000, 0.229000,
     1  0.229000, 0.229000, 0.229000, 0.229100, 0.229100, 0.229100, 0.229200, 0.229300, 0.229400,
     1  0.229500, 0.229600, 0.229700, 0.230200, 0.230500, 0.231100, 0.231300, 0.231900, 0.232600,
     1  0.233400, 0.234800, 0.236100, 0.237400, 0.238600, 0.243300, 0.247700, 0.249700, 0.253300,
     1  0.260600, 0.264100, 0.267400, 0.274600, 0.284700, 0.289500, 0.294200, 0.303200, 0.312000,
     1  0.322700, 0.332900, 0.337800, 0.342700, 0.352000, 0.361000, 0.381000, 0.399300, 0.414200,
     1  0.418000, 0.425200, 0.435300, 0.444400, 0.449400, 0.454200, 0.458700, 0.462900, 0.466800,
     1  0.470300, 0.472900, 0.474300, 0.475600, 0.476700, 0.477200, 0.477700, 0.478100, 0.478500,
     1  0.478900, 0.479200, 0.479300, 0.479400, 0.479500, 0.479600, 0.479700, 0.479800, 0.479800,
     1  0.479800, 0.479800, 0.479900, 0.479900, 0.479900, 0.479900, 0.479900, 0.479900, 0.4800,
     1  0.4800, 0.4800, 0.4800, 0.4800, 0.4800, 0.4800, 0.4800, 0.4800, 0.480000,
     + 0.2899/)
         phi6= (/ 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996,
     1  0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996,
     1  0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996, 0.014996,
     1  0.014996, 0.014996, 0.014996, 0.014994, 0.014994, 0.014993, 0.014993, 0.014991, 0.014988,
     1  0.014987, 0.014981, 0.014975, 0.014970, 0.014964, 0.014928, 0.014895, 0.014881, 0.014832,
     1  0.014731, 0.014684, 0.014639, 0.014513, 0.014239, 0.014110, 0.013985, 0.013747, 0.013493,
     1  0.012938, 0.012429, 0.012190, 0.011962, 0.011532, 0.011133, 0.009769, 0.008660, 0.007829,
     1  0.007620, 0.007244, 0.006739, 0.006325, 0.006163, 0.006014, 0.005876, 0.005749, 0.005678,
     1  0.005613, 0.005573, 0.005558, 0.005544, 0.005533, 0.005529, 0.005527, 0.005524, 0.005521,
     1  0.005519, 0.005518, 0.005518, 0.005518, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517,
     1  0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517,
     1  0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517, 0.005517,
     + 0.006718/)
      phi7= (/ 580.0, 580.0, 580.0, 580.0, 580.0, 580.0,
     1  580.0, 580.0, 580.0, 579.9, 579.9, 579.9, 579.9, 579.9, 579.9,
     1  579.9, 579.8, 579.8, 579.8, 579.7, 579.7, 579.6, 579.6, 579.5,
     1  579.4, 579.3, 579.2, 578.8, 578.6, 578.2, 578.0, 577.7, 577.2,
     1  576.7, 576.0, 575.2, 574.6, 573.9, 571.6, 569.5, 568.5, 566.9,
     1  563.6, 562.0, 560.5, 557.3, 552.6, 550.4, 548.3, 544.1, 540.0,
     1  534.0, 528.4, 525.7, 523.0, 517.8, 512.9, 497.1, 482.7, 468.7,
     1  463.9, 454.8, 441.9, 429.9, 419.5, 409.8, 400.5, 391.8, 379.6,
     1  368.5, 359.8, 353.7, 348.1, 343.1, 340.2, 337.5, 334.9, 332.5,
     1  330.0, 327.8, 326.7, 326.1, 325.1, 324.1, 323.3, 322.9, 322.7,
     1  322.5, 322.1, 321.7, 321.6, 321.4, 321.2, 321.1, 320.9, 320.7,
     1  320.6, 320.4, 320.4, 320.3, 320.2, 320.2, 320.2, 320.1, 320.1,
     + 459.0/)
         phi8= (/ 0.070000, 0.070000, 0.070000, 0.070000, 0.070000, 0.070000,
     1  0.070000, 0.070000, 0.070000, 0.070000, 0.070000, 0.070000, 0.070000, 0.070000, 0.070000,
     1  0.070000, 0.070000, 0.069800, 0.069600, 0.069400, 0.069200, 0.068600, 0.067900, 0.067100,
     1  0.066200, 0.065400, 0.064600, 0.063500, 0.062500, 0.060200, 0.059200, 0.056000, 0.049400,
     1  0.040700, 0.030600, 0.019900, 0.008900,-0.001900,-0.022300,-0.040100,-0.047900,-0.054800,
     1 -0.066500,-0.071300,-0.075600,-0.082500,-0.087500,-0.089500,-0.091200,-0.093900,-0.096000,
     1 -0.097500,-0.098700,-0.099100,-0.099400,-0.099800,-0.099800,-0.098300,-0.094800,-0.089600,
     1 -0.087600,-0.083400,-0.076500,-0.069300,-0.062000,-0.054900,-0.047900,-0.041200,-0.028500,
     1 -0.016700,-0.005700, 0.004500, 0.014000, 0.022900, 0.031300, 0.039300, 0.046900, 0.054400,
     1  0.068700, 0.082600, 0.089500, 0.096300, 0.109800, 0.123200, 0.136600, 0.149800, 0.156200,
     1  0.162500, 0.174600, 0.185900, 0.196400, 0.206000, 0.214700, 0.222500, 0.229500, 0.243500,
     1  0.253200, 0.259500, 0.263500, 0.266000, 0.267500, 0.268300, 0.268600, 0.268500, 0.268200,
     + 0.1138/)
c Above phi? coeffs from a file called cy2007.coe emailed by Brian Chiou, Oct 2007
c Added PGV coeffs Jan 29 2009. SHarmsen.
c Soil effect: linear response
      if(Z1 .le.5.)then
        Z1 = exp(28.5-3.82/8*log(V_S30**8+378.8**8)) 
        write(6,*)'CY uses a default Z1 (m) of ',Z1
        endif
        a(ip) = phi1(iprd) * min(log(V_S30/1130.), 0.)

c Soil effect: nonlinear response. 1130-360= 770 
        b(ip) = phi2(iprd) *
     1(exp(phi3(iprd)*(min(V_S30,1130.)-360.))-exp(phi3(iprd)*(770.)))
        c (ip)= phi4(iprd)

c Modificaiton to ln(Vs30) scaling: bedrock depth (Z1)
c NOTE: max(0,Z1-15) is capped at 300 to avoid overflow of function cosh
        rkdepth(ip) = phi5(iprd) *
     1        ( 1 - 1.0/cosh(phi6(iprd)*max(0.,Z1-phi7(iprd)))) +
     1        phi8(iprd)/cosh(0.15*min(max(0., Z1-15.),300.))

      return
      end subroutine CY2007I

      subroutine getCBNGA1107(iper,ip,amag,rrup,rjb,H,
     1 vs,d,gndout,sigmaf)
c....................................................................
c  Campbell and Bozorgnia NGA regression relation, Nov 15 2007. Update of
c  November 2007 affects aleatory uncert only, add PGD (index 24).
c  from getCampNGA0601 based on code of Yuehua Zeng, USGS
c Steve Harmsen. 
c Note: this version saves pgahr computed first time, Other periods need pgar
c hard-rock value for subsequent nonlin soil calcs.
c Coeffs from Campbell_Bozorgnia_NGA_MAR08_EERI_SPECTRA.txt (xls) November 2007
c  input: 
c        
c        iper : index for period in iatten(ip,ia) etc.
c        ip      :index for period in Pd array below
c         amag : moment magnitude
c         rrup : closest fault distance
c         rjb  : distance to the fault projection on the surface
c         rak  : rake angle (d). Replaced by style of faulting coming in in common, SH.
c         vs   : site S-velocity in m/s (should be Vs30 )
c         H    : depth to top of the fault
c         d    : sediment depth now defined as Z2.5, depth to 2.5 km/s Vs.
c
c  output: gnd   : ln(ground motion spectral value, g)
c          
c          sigmaf= 1/sig_t/sqrt2 = sigma-factor in natural log
c....................................................................
c random component sig_t:
c sig_t = sqrt (tau**2+sigma**2+chisq) see coeff values below.
c but we use a lower-sigma "geometric" oct 2007.
      parameter (sqrt2=1.414213562,rockcoef=0.24033595,expm75=0.47236655)
c rockcoef = alog(1100/k1(1)), expm75=exp(-0.75)         
c C-B NGA Report: november 2007. 
        parameter (np=25)
       common/fix_sigma/fix_sigma,sigma_fx
       common/sdi/sdi,dy_sdi,fac_sde
       real, dimension(8) :: fac_sde
      common/pgac/pgacalc
      common/mech/ss,rev,normal,obl
      common/dipinf/dipang,cosDELTA,cdipsq,cyhwfac,cbhwfac
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      real sln_Ab, sln_yb      !sigma at base of soil column. New additions
c of late Nov 2007. These are in the Ken Campbell & Yusef Bozorgnia paper to 
c Earthquake Spectra March 2008
      logical l_gnd_ep(8),fix_sigma

      real sigsqu,spgasq/0.228484/,sigma_fx,dy_sdi,sdisd
      logical ss,sdi,rev,normal,obl,pgacalc,first/.true./
        real pgar
        save pgar,first
        real,dimension(np):: Pd,c0,c1,c2,c3,c4,c5,c6,c7,c8,
     1 c9,c10,c11,c12,k1,k2,k3,
     2 rhos,slny,tlny,slnAF,sC
c Coeffs from Campbell_Bozorgnia_NGA_MAR08_EERI_SPECTRA.txt  (xls) November 2007
c 25 periods available, -1 is PGV. -2 is PGD. -3 is CAV. 0.00 is PGA.
c SHarmsen Nov 15 2006. Pd=spectral period vector. Added CAV coeffs Oct 30 2012.
      Pd=(/0.010,0.020,0.030,0.050,0.075,0.100,0.150,0.200,0.250,0.300,0.400,0.500,0.750,
     + 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5,10.0, 0.0,-1.0,-2.0,-3.0/)
       c0=(/-1.715,-1.680,-1.552,-1.209,-0.657,-0.314,-0.133,-0.486,-0.890,-1.171,-1.466,-2.569,-4.844,
     + -6.406, -8.692, -9.701,-10.556,-11.212,-11.684,-12.505,-13.087, -1.715,  0.954, -5.270,-4.354/)
       c1=(/ 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.500, 0.656, 0.972,
     + 1.196, 1.513, 1.600, 1.600, 1.600, 1.600, 1.600, 1.600, 0.500, 0.696, 1.600,0.942/)
       c2=(/-0.530,-0.530,-0.530,-0.530,-0.530,-0.530,-0.530,-0.446,-0.362,-0.294,-0.186,-0.304,-0.578,
     +-0.772,-1.046,-0.978,-0.638,-0.316,-0.070,-0.070,-0.070,-0.530,-0.309,-0.070,-0.178/)
       c3=(/-0.262,-0.262,-0.262,-0.267,-0.302,-0.324,-0.339,-0.398,-0.458,-0.511,-0.592,-0.536,-0.406,
     +-0.314,-0.185,-0.236,-0.491,-0.770,-0.986,-0.656,-0.422,-0.262,-0.019, 0.000,-0.346/)
       c4=(/-2.118,-2.123,-2.145,-2.199,-2.277,-2.318,-2.309,-2.220,-2.146,-2.095,-2.066,-2.041,-2.000,
     +-2.000,-2.000,-2.000,-2.000,-2.000,-2.000,-2.000,-2.000,-2.118,-2.016,-2.000,-1.309/)
       c5=(/ 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170,
     + 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170, 0.170,0.087/)
       c6=(/ 5.600, 5.600, 5.600, 5.740, 7.090, 8.050, 8.790, 7.600, 6.580, 6.040, 5.300, 4.730, 4.000,
     + 4.000, 4.000, 4.000, 4.000, 4.000, 4.000, 4.000, 4.000, 5.600, 4.000, 4.000,7.24/)
       c7=(/ 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280, 0.280,
     + 0.255, 0.161, 0.094, 0.000, 0.000, 0.000, 0.000, 0.000, 0.280, 0.245, 0.000,0.111/)
       c8=(/-0.120,-0.120,-0.120,-0.120,-0.120,-0.099,-0.048,-0.012, 0.000, 0.000, 0.000, 0.000, 0.000,
     + 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,-0.120, 0.000, 0.000,-0.108/)
       c9=(/ 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490, 0.490,
     + 0.490, 0.490, 0.371, 0.154, 0.000, 0.000, 0.000, 0.000, 0.490, 0.358, 0.000,0.362/)
      c10=(/1.058,1.102,1.174,1.272,1.438,1.604,1.928, 2.194,2.351,2.460,2.587,2.544,2.133,
     + 1.571, 0.406,-0.456,-0.820,-0.820,-0.820,-0.820,-0.820, 1.058, 1.694,-0.820,2.549/)
      c11=(/ 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.040, 0.077,
     + 0.150, 0.253, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.040, 0.092, 0.300,0.090/)
      c12=(/ 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.610, 0.883, 1.000,
     + 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.610, 1.000, 1.000,1.277/)
       k1=(/ 865., 865., 908.,1054.,1086.,1032., 878., 748., 654., 587., 503., 457., 410.,
     + 400., 400., 400., 400., 400., 400., 400., 400., 865., 400., 400.,400./)
       k2=(/-1.186,-1.219,-1.273,-1.346,-1.471,-1.624,-1.931,-2.188,-2.381,-2.518,-2.657,-2.669,-2.401,
     +-1.955,-1.025,-0.299, 0.000, 0.000, 0.000, 0.000, 0.000,-1.186,-1.955, 0.000,-2.690/)
       k3=(/ 1.839, 1.840, 1.841, 1.843, 1.845, 1.847, 1.852, 1.856, 1.861, 1.865, 1.874, 1.883, 1.906,
     + 1.929, 1.974, 2.019, 2.110, 2.200, 2.291, 2.517, 2.744, 1.839, 1.929, 2.744, 1.0/)
c some revised coeffs. Mar 2008 Eq Spectra 
      slny=(/ 0.478, 0.480, 0.489, 0.510, 0.520, 0.531,
     + 0.532, 0.534, 0.534, 0.544, 0.541, 0.550, 0.568, 0.568, 0.564,
     + 0.571, 0.558, 0.576, 0.601, 0.628, 0.667, 0.478, 0.484, 0.667,0.371/)
        tlny=(/ 0.219, 0.219, 0.235, 0.258, 0.292, 0.286,
     + 0.280, 0.249, 0.240, 0.215, 0.217, 0.214, 0.227, 0.255, 0.296,
     + 0.296, 0.326, 0.297, 0.359, 0.428, 0.485, 0.219, 0.203, 0.485,0.196 /)
      slnAF=(/ 0.300, 0.300, 0.300, 0.300, 0.300, 0.300,
     + 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300,
     + 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300, 0.300,0.300/)
      sC   =(/ 0.166, 0.166, 0.165, 0.162, 0.158, 0.170,
     + 0.180, 0.186, 0.191, 0.198, 0.206, 0.208, 0.221, 0.225, 0.222,
     + 0.226, 0.229, 0.237, 0.237, 0.271, 0.290, 0.166, 0.190, 0.290,0.089/)
          rhos=(/ 1.000, 0.999, 0.989, 0.963, 0.922, 0.898,
     + 0.890, 0.871, 0.852, 0.831, 0.785, 0.735, 0.628, 0.534, 0.411,
     + 0.331, 0.289, 0.261, 0.200, 0.174, 0.174, 1.000, 0.691, 0.174,0.735/)
        cs= 1.88
        cn= 1.18      !these 2 c-values did not change from the 1-06 report.
c
        i=ip      !period index
c return to ground-motion independent sigma, C&B 11/06.
            tausq=tlny(i)**2
            sigsq=slny(i)**2
            sigsqu = sigsq      !add for use in non-linear sigma. Nov 2007
c            chisq=sC(i)**2
c random horizontal component, include the chisq boost compared to geom. mean.
c            sigt=sqrt(tausq+sigsq+chisq)      !campbell has a kchi=1 factor for USGS
         sigt = sqrt(tausq + sigsq)      !campbell note Oct 2007. Be consistent with others,
c who use geometric mean.
c  Originally CB recom. RHC for use by USGS. Other apps might want geometric mean, whose sigt is lower.
c sigmaf is a factor in the normal prob. 
            sigmaf=1.0/sigt/sqrt2
c           f(i)=1./per(i)      !dont need
c f1=median dependence on magnitude. f2= joint M,R dependence
           f1=c0(i)+c1(i)*amag
           if(amag.gt.5.5)f1=f1+c2(i)*(amag-5.5)
           if(amag.gt.6.5)f1=f1+c3(i)*(amag-6.5)
           f2=(c4(i)+c5(i)*amag)*0.5*alog(rrup*rrup+c6(i)*c6(i))
           f3=0.0
           f4=0.0
c f3 dependence on source mechanism; f4 hanging-wall dependence N&R both
c           if(rak.gt.30.0.and.rak.lt.150.0)then
      if(rev)then      !less precise but corresponds to current classification
             if(H.gt.1.0)then
               f3=c7(i)
             else
               f3=c7(i)*H
             endif
c           elseif(rak.gt.-150.0.and.rak.lt.-30.0)then
      elseif(normal)then      !less precise but corresponds to current classification
             f3=c8(i)
         endif      !reverse or normal slip. Hangingwall can be on for normal too.
             if(amag.gt.6.0.and.H.lt.20.0)then
             if(rjb.eq.0.)then
             fhw=1.0      !eqn 7 jan 31 report
             elseif(H.lt.1.)then
c June 2006: new smoothing factor when rupture makes it to the surface or nearly makes it.
             rjbf=sqrt(rjb*rjb+1.0)
c next line is a correction of jan 16 2007. SH.
             fhw= (max(rrup,rjbf)-rjb) /max(rrup,rjbf)
             else
             fhw=(rrup-rjb)/rrup
              endif
             f4=c9(i)*fhw*cbhwfac*(20.0-H)*0.05      
c cbhwfac new july 2006, computed in main
               if(amag.lt.6.5)f4=f4*(amag-6.0)*2.0
             endif
c f6 = Shallow sediment thickness dependence:
c Null depth effect for 1<d<3. 
           f6=0.0
           if(d.lt.1.0)then
             f6=c11(i)*(d-1.0)
           elseif(d.gt.3.0)then
             f6=c12(i)*k3(i)*expm75*(1.0-exp(-0.25*(d-3.)))
           endif
         vscap = min(vs,1100.)
           vsrk=vscap/k1(i)
c Ken Campbell email mar 6 2008: cap Vs at 1100 m/s for these calcs. See Eq Spectra
c report.
           if(vsrk.lt.1.0)then
c Code must save pgar for subsequent ip loop indexes. PGA has to be computed first. SH.
c Also, ip index loop has to be inside any distance and mag index loops. It is in current code.
             if(pgacalc)pgar=exp(f1+f2+f3+f4+f6+(c10(1)
     +                      +k2(1)*cn)*rockcoef)
             csfac=cs*vsrk**cn
             f5=c10(i)*alog(vsrk)
     +            +k2(i)*alog((pgar+csfac)/(pgar+cs))
c additional alpha term CB eqn (17) Sept 1 2006, in sigma computation
c need spgasq=spga**2 & tpgasq=tpga**2. Use the sigma(1) & tau(1). 
c these are conveniently defined in a statement to save calculations.
            alpha=k2(i)*pgar*(1./(pgar+csfac)-1./(pgar+cs))
            alfsq=alpha*alpha 
c tau no dependency on soil nonlinearity in the CB EQ Spectra article, see. eqn (14) & eqn (16)
c sln_Ab sigma of pga at base of site profile. see "doc" Eq Spectra paper, p 13 midway thru.
            sln_Ab= sqrt(spgasq - slnAF(1)**2)
            sln_yb= sqrt(slny(i)**2 - slnAF(i)**2)
c sln_yb is sigma for spectral period j at base of soil column. Campbell corr.
c Nov 30 2007
            sigsq=sigsqu +alfsq*sln_Ab**2+2.*alpha*rhos(i)*sln_Ab*sln_yb
c nonlinear motion-dependent sigma, but no motion-dependent tau.
           sigt = sqrt(tausq + sigsq)      
            sigmaf=1.0/sigt/sqrt2
            else
             alpha=0.
             f5=(c10(i)+k2(i)*cn)*alog(vsrk)
           endif
           if(fix_sigma)sigmaf=1.0/sigma_fx/sqrt2
           gnd=f1+f2+f3+f4+f5+f6
       gndout(1)=gnd
c for QA tests, look at near-source log mean values
c      if(rjb.lt.0.1.and.amag.gt.7.)print *,rjb,rrup,amag,gnd,sigt,Pd(ip),H,pgar,f5
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
c
c      if(first)write(6,*)'CB1107 i,gnd pgar,sigt,sigt_o,amag,rrup,alpha,vs30=',vs
c      first=.false. 
c      write(6,*)i,exp(gnd),sigt,sig_t(i),pgar,amag,rrup,alpha
c sigmaf is a factor in the normal prob
       if(sdi)then
       sde=gndout(1)+fac_sde(iper)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)
       gndout(1) = sdi_ratio(Pd(ip),amag,rhat,sigt,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
         if(l_gnd_ep(iper))then
         gndout(2)= gndout(1)+gnd_ep(ide,ime,iper)
         gndout(3)= gndout(1)-gnd_ep(ide,ime,iper)
         endif
         endif      !if(sdi)
        return
      end subroutine getCBNGA1107

      subroutine getAB06(iper,ip,ir,xmag0,dist0,gndout,sigma,sigmaf,v30)
c Atkinson and Boore BSSA 2006. A CEUS relation
c      input:
c iper=period index in l_gnd_ep; ip = period index in coefficients below.
c ir=1 or 3 use BC model. 3 for 200bar
c ir=2 or 4 use hardrock model. 4 for 200bar.
c if stress value is greater than 140 bars, change the numerator of stressfac to this value.
c if stress value is less than 140 bars, there is a net reduction in gnd
c Stress is not currently an input variable, hardwired to 200 bars.
c modified from version written by Oliver Boyd. Steve Harmsen Nov 13 2006. A Frankel feb 13
c the stress-factor was corrected in a BSSA erratum, June 2007, p 1032. This repair was made
c to this code, June 11, 2007. Steve Harmsen.
c different sets of coeffs for hardrock (ir=2&4) and bc rock (ir<2). 
c Change Mar 17 2008: AB06 doc says to use
c the hardrock coeffs if Vs30 >= 2000. However, the Nehrp A begins at 1520 m/s or so, and A usually
c corresponds to hard rock. You will get the hard rock if Vs30>=2000 in input file. 
c 
c   near-source terms not used here.
c xmag0 = moment M
c dist0 = src-station Rcd (km)
c v30==site vs30 (m/s)
c Output:
c       gndout: a 3-component vector, component1 is central gnd. Others are not used
c       for CEUS PSHA in most cases.
c        sigma = model sigma,
c       sigmaf = sigma factor used in prob calcs
      parameter (np=26)
      parameter (sqrt2=1.414213562,stressfac=0.5146)      !stressfac =alog10(200./140.)/alog10(2.)
      parameter (gfac=6.8875526,sfac=2.3025851,tfac=-0.5108256)      !log(980),ln(10),
c  and ln(60/100), resp.
      parameter (fac70=1.8450980,fac140=2.1461280,facv1=-0.5108256,facv2=-0.9295360)      
c  log10(70),log10(140),ln(v1/v2),ln(v2/vref),resp
        parameter (vref = 760., v1 = 180., v2 = 300.)      !for siteamp calcs
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
       common/ceus_sig/lceus_sigma,ceus_sigma
       common/hardrock/hardrock
      real gnd_ep(3,3,8),gndout(3)
      logical  lceus_sigma,l_gnd_ep(8),hardrock
      real f0,f1,f2,R,M,bnl,S
c del,m1, and mh are arrays associated with variable stress drop model. Add feb16 2007
      real,dimension(np):: Frq,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,bln,b1,b2,del,m1,mh
c the frequencies according to AB2006
      Frq = (/2.00e-01,2.50e-01,3.20e-01,4.00e-01,5.00e-01,6.30e-01,8.00e-01,1.00e+00,
     1       1.26e+00,1.59e+00,2.00e+00,2.52e+00,3.17e+00,3.99e+00,5.03e+00,6.33e+00,
     1       7.97e+00,1.00e+01,1.26e+01,1.59e+01,2.00e+01,2.52e+01,3.18e+01,4.00e+01,
     1       0.00e+00,-1.00e+00/)
c frequency - independent sigma
      if(lceus_sigma)then
      sigma=ceus_sigma	!special study 3/2011
      else
      sigma = 0.3*sfac
      endif
      sigmaf = 1./sqrt2/sigma
        if (ir.eq.2.or.ir.eq.4) then
c hr coefficients
        c1 = (/-5.41e+00,-5.79e+00,-6.04e+00,-6.17e+00,-6.18e+00,-6.04e+00,-5.72e+00,
     1        -5.27e+00,-4.60e+00,-3.92e+00,-3.22e+00,-2.44e+00,-1.72e+00,-1.12e+00,
     1        -6.15e-01,-1.46e-01,2.14e-01,4.80e-01,6.91e-01,9.11e-01,1.11e+00,1.26e+00,
     1        1.44e+00,1.52e+00,9.07e-01,-1.44e+00/)
        c2 = (/1.71e+00,1.92e+00,2.08e+00,2.21e+00,2.30e+00,2.34e+00,2.32e+00,2.26e+00,
     1        2.13e+00,1.99e+00,1.83e+00,1.65e+00,1.48e+00,1.34e+00,1.23e+00,1.12e+00,
     1        1.05e+00,1.02e+00,9.97e-01,9.80e-01,9.72e-01,9.68e-01,9.59e-01,9.60e-01,
     1        9.83e-01,9.91e-01/)
        c3 = (/-9.01e-02,-1.07e-01,-1.22e-01,-1.35e-01,-1.44e-01,-1.50e-01,-1.51e-01,
     1        -1.48e-01,-1.41e-01,-1.31e-01,-1.20e-01,-1.08e-01,-9.74e-02,-8.72e-02,
     1        -7.89e-02,-7.14e-02,-6.66e-02,-6.40e-02,-6.28e-02,-6.21e-02,-6.20e-02,
     1        -6.23e-02,-6.28e-02,-6.35e-02,-6.60e-02,-5.85e-02/)
        c4 = (/-2.54e+00,-2.44e+00,-2.37e+00,-2.30e+00,-2.22e+00,-2.16e+00,-2.10e+00,
     1        -2.07e+00,-2.06e+00,-2.05e+00,-2.02e+00,-2.05e+00,-2.08e+00,-2.08e+00,
     1        -2.09e+00,-2.12e+00,-2.15e+00,-2.20e+00,-2.26e+00,-2.36e+00,-2.47e+00,
     1        -2.58e+00,-2.71e+00,-2.81e+00,-2.70e+00,-2.70e+00/)
        c5 = (/2.27e-01,2.11e-01,2.00e-01,1.90e-01,1.77e-01,1.66e-01,1.57e-01,1.50e-01,
     1        1.47e-01,1.42e-01,1.34e-01,1.36e-01,1.38e-01,1.35e-01,1.31e-01,1.30e-01,
     1        1.30e-01,1.27e-01,1.25e-01,1.26e-01,1.28e-01,1.32e-01,1.40e-01,1.46e-01,
     1        1.59e-01,2.16e-01/)
        c6 = (/-1.27e+00,-1.16e+00,-1.07e+00,-9.86e-01,-9.37e-01,-8.70e-01,-8.20e-01,
     1        -8.13e-01,-7.97e-01,-7.82e-01,-8.13e-01,-8.43e-01,-8.89e-01,-9.71e-01,
     1        -1.12e+00,-1.30e+00,-1.61e+00,-2.01e+00,-2.49e+00,-2.97e+00,-3.39e+00,
     1        -3.64e+00,-3.73e+00,-3.65e+00,-2.80e+00,-2.44e+00/)
        c7 = (/1.16e-01,1.02e-01,8.95e-02,7.86e-02,7.07e-02,6.05e-02,5.19e-02,4.67e-02,
     1        4.35e-02,4.30e-02,4.44e-02,4.48e-02,4.87e-02,5.63e-02,6.79e-02,8.31e-02,
     1        1.05e-01,1.33e-01,1.64e-01,1.91e-01,2.14e-01,2.28e-01,2.34e-01,2.36e-01,
     1        2.12e-01,2.66e-01/)
        c8 = (/9.79e-01,1.01e+00,1.00e+00,9.68e-01,9.52e-01,9.21e-01,8.56e-01,8.26e-01,
     1        7.75e-01,7.88e-01,8.84e-01,7.39e-01,6.10e-01,6.14e-01,6.06e-01,5.62e-01,
     1        4.27e-01,3.37e-01,2.14e-01,1.07e-01,-1.39e-01,-3.51e-01,-5.43e-01,-6.54e-01,
     1        -3.01e-01,8.48e-02/)
        c9 = (/-1.77e-01,-1.82e-01,-1.80e-01,-1.77e-01,-1.77e-01,-1.73e-01,-1.66e-01,
     1        -1.62e-01,-1.56e-01,-1.59e-01,-1.75e-01,-1.56e-01,-1.39e-01,-1.43e-01,
     1        -1.46e-01,-1.44e-01,-1.30e-01,-1.27e-01,-1.21e-01,-1.17e-01,-9.84e-02,
     1        -8.13e-02,-6.45e-02,-5.50e-02,-6.53e-02,-6.93e-02/)
        c10 = (/-1.76e-04,-2.01e-04,-2.31e-04,-2.82e-04,-3.22e-04,-3.75e-04,-4.33e-04,
     1        -4.86e-04,-5.79e-04,-6.95e-04,-7.70e-04,-8.51e-04,-9.54e-04,-1.06e-03,
     1        -1.13e-03,-1.18e-03,-1.15e-03,-1.05e-03,-8.47e-04,-5.79e-04,-3.17e-04,
     1        -1.23e-04,-3.23e-05,-4.85e-05,-4.48e-04,-3.73e-04/)
        else
c bc coefficients from AB06. Frankel corrected these feb 13. SH
        c1 = (/-4.85e+00,-5.26e+00,-5.59e+00,-5.80e+00,-5.85e+00,-5.75e+00,-5.49e+00,
     1         -5.06e+00,-4.45e+00,-3.75e+00,-3.01e+00,-2.28e+00,-1.56e+00,-8.76e-01,
     1         -3.06e-01,1.19e-01,5.36e-01,7.82e-01,9.67e-01,1.11e+00,1.21e+00,1.26e+00,
     1         1.19e+00,1.05e+00,5.23e-01,-1.66e+00/)
        c2 = (/1.58e+00,1.79e+00,1.97e+00,2.13e+00,2.23e+00,2.29e+00,2.29e+00,2.23e+00,
     1        2.12e+00,1.97e+00,1.80e+00,1.63e+00,1.46e+00,1.29e+00,1.16e+00,1.06e+00,
     1        9.65e-01,9.24e-01,9.03e-01,8.88e-01,8.83e-01,8.79e-01,8.88e-01,9.03e-01,
     1        9.69e-01,1.05e+00/)
        c3 = (/-8.07e-02,-9.79e-02,-1.14e-01,-1.28e-01,-1.39e-01,-1.45e-01,-1.48e-01,
     1        -1.45e-01,-1.39e-01,-1.29e-01,-1.18e-01,-1.05e-01,-9.31e-02,-8.19e-02,
     1        -7.21e-02,-6.47e-02,-5.84e-02,-5.56e-02,-5.48e-02,-5.39e-02,-5.44e-02,
     1        -5.52e-02,-5.64e-02,-5.77e-02,-6.20e-02,-6.04e-02/)
        c4 = (/-2.53e+00,-2.44e+00,-2.33e+00,-2.26e+00,-2.20e+00,-2.13e+00,-2.08e+00,
     1        -2.03e+00,-2.01e+00,-2.00e+00,-1.98e+00,-1.97e+00,-1.98e+00,-2.01e+00,
     1        -2.04e+00,-2.05e+00,-2.11e+00,-2.17e+00,-2.25e+00,-2.33e+00,-2.44e+00,
     1        -2.54e+00,-2.58e+00,-2.57e+00,-2.44e+00,-2.50e+00/)
        c5 = (/2.22e-01,2.07e-01,1.91e-01,1.79e-01,1.69e-01,1.58e-01,1.50e-01,1.41e-01,
     1        1.36e-01,1.31e-01,1.27e-01,1.23e-01,1.21e-01,1.23e-01,1.22e-01,1.19e-01,
     1        1.21e-01,1.19e-01,1.22e-01,1.23e-01,1.30e-01,1.39e-01,1.45e-01,1.48e-01,
     1        1.47e-01,1.84e-01/)
        c6 = (/-1.43e+00,-1.31e+00,-1.20e+00,-1.12e+00,-1.04e+00,-9.57e-01,-9.00e-01,
     1        -8.74e-01,-8.58e-01,-8.42e-01,-8.47e-01,-8.88e-01,-9.47e-01,-1.03e+00,
     1        -1.15e+00,-1.36e+00,-1.67e+00,-2.10e+00,-2.53e+00,-2.88e+00,-3.04e+00,
     1        -2.99e+00,-2.84e+00,-2.65e+00,-2.34e+00,-2.30e+00/)
        c7 = (/1.36e-01,1.21e-01,1.10e-01,9.54e-02,8.00e-02,6.76e-02,5.79e-02,5.41e-02,
     1        4.98e-02,4.82e-02,4.70e-02,5.03e-02,5.58e-02,6.34e-02,7.38e-02,9.16e-02,
     1        1.16e-01,1.48e-01,1.78e-01,2.01e-01,2.13e-01,2.16e-01,2.12e-01,2.07e-01,
     1        1.91e-01,2.50e-01/)
        c8 = (/6.34e-01,7.34e-01,8.45e-01,8.91e-01,8.67e-01,8.67e-01,8.21e-01,7.92e-01,
     1        7.08e-01,6.77e-01,6.67e-01,6.84e-01,6.50e-01,5.81e-01,5.08e-01,5.16e-01,
     1        3.43e-01,2.85e-01,1.00e-01,-3.19e-02,-2.10e-01,-3.91e-01,-4.37e-01,-4.08e-01,
     1        -8.70e-02,1.27e-01/)
        c9 = (/-1.41e-01,-1.56e-01,-1.72e-01,-1.80e-01,-1.79e-01,-1.79e-01,-1.72e-01,
     1        -1.70e-01,-1.59e-01,-1.56e-01,-1.55e-01,-1.58e-01,-1.56e-01,-1.49e-01,
     1        -1.43e-01,-1.50e-01,-1.32e-01,-1.32e-01,-1.15e-01,-1.07e-01,-9.00e-02,
     1        -6.75e-02,-5.87e-02,-5.77e-02,-8.29e-02,-8.70e-02/)
        c10 = (/-1.61e-04,-1.96e-04,-2.45e-04,-2.60e-04,-2.86e-04,-3.43e-04,-4.07e-04,
     1        -4.89e-04,-5.75e-04,-6.76e-04,-7.68e-04,-8.59e-04,-9.55e-04,-1.05e-03,
     1        -1.14e-03,-1.18e-03,-1.13e-03,-9.90e-04,-7.72e-04,-5.48e-04,-4.15e-04,
     1        -3.88e-04,-4.33e-04,-5.12e-04,-6.30e-04,-4.27e-04/)
      endif
c Soil amplification
      bln = (/-7.52e-01,-7.45e-01,-7.40e-01,-7.35e-01,-7.30e-01,-7.26e-01,-7.16e-01,
     1        -7.00e-01,-6.90e-01,-6.70e-01,-6.00e-01,-5.00e-01,-4.45e-01,-3.90e-01,
     1        -3.06e-01,-2.80e-01,-2.60e-01,-2.50e-01,-2.32e-01,-2.49e-01,-2.86e-01,
     1        -3.14e-01,-3.22e-01,-3.30e-01,-3.61e-01,-6.00e-01/)
      b1 = (/-3.00e-01,-3.10e-01,-3.30e-01,-3.52e-01,-3.75e-01,-3.95e-01,-3.40e-01,
     1        -4.40e-01,-4.65e-01,-4.80e-01,-4.95e-01,-5.08e-01,-5.13e-01,-5.18e-01,
     1        -5.21e-01,-5.28e-01,-5.60e-01,-5.95e-01,-6.37e-01,-6.42e-01,-6.43e-01,
     1        -6.09e-01,-6.18e-01,-6.24e-01,-6.41e-01,-4.95e-01/)
      b2 = (/0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00,0.00e+00,
     1        -2.00e-03,-3.10e-02,-6.00e-02,-9.50e-02,-1.30e-01,-1.60e-01,-1.85e-01,
     1        -1.85e-01,-1.40e-01,-1.32e-01,-1.17e-01,-1.05e-01,-1.05e-01,-1.05e-01,
     1        -1.08e-01,-1.15e-01,-1.44e-01,-6.00e-02/)
c stress adjustment factors, table 7
       del=(/0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,
     1 0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.15,0.11/)
       m1=(/6.,5.75,5.5,5.25,5.,4.84,4.67,4.5,4.34,4.17,4.,3.65,3.3,2.9,2.5,1.85,1.15,0.5,
     1 0.34,0.17,0.,0.,0.,0.,0.5,2.0/)
       mh=(/8.5,8.37,8.25,8.12,8.,7.7,7.45,7.2,6.95,6.7,6.5,6.37,6.25,6.12,6.0,5.84,
     1 5.67,5.5,5.34,5.17,5.0,5.0,5.0,5.0,5.5,5.5/)
      if(ir.gt.2)then
       diff=max(xmag0-m1(ip),0.)
c sf2 is supposed to be eqn(6) of AB06 paper.
       sf2=stressfac*min( del(ip)+0.05, 0.05+del(ip)*diff/(mh(ip)-m1(ip)))
c       write(6,*)diff,sf2,ir,ip,iper,dist0,' diff sf2 ir ip iper dist0'
       else
       sf2=0.0
       endif
      M = xmag0
c R: I put a lower limit of 5 km SH Nov 2006. For near-surface faults singularity is possible
      R = max(5.0,dist0)
      rfac=alog10(R)
      f0 = max(1.-rfac,0.)
      f1 = min(rfac,fac70)
      f2 = max(rfac-fac140,0.)
      if(v30.gt.0.) then
c compute pga on rock
        gnd = c1(25) + c2(25)*M + c3(25)*M*M + (c4(25)+c5(25)*M)*f1 +
     1        (c6(25)+c7(25)*M)*f2 + (c8(25)+c9(25)*M)*f0 + c10(25)*R + sf2
c apply stress factor before nonlinear adjustments, which occur in eqn (7) of ab paper.
        if(v30.le.v1) then
          bnl = b1(ip);
        elseif(v30.le.v2) then
          bnl = (b1(ip) - b2(ip))*log(v30/v2)/facv1 + b1(ip);
        elseif(v30.le.vref) then
          bnl = b2(ip)*log(v30/vref)/facv2
        else
          bnl = 0.;
        endif
        if(ir.eq.2.or.ir.eq.4)then
         S=0.            !hard rock no site term
        elseif(10**gnd.le.60.) then
          S = bln(ip)*log(v30/vref) + bnl*tfac
        else
          S = bln(ip)*log(v30/vref) + bnl*log((10**gnd)/100.)
        endif
c need to take alog10(exp(S)) according to eqns. 7a and 7b AB2006. p. 2200 bssa
c This correction does not affect rock at BC boundary, but does affect soil calcs.
      S = alog10(exp(S))      	!new nov 26 2007.
      endif
      gnd = c1(ip) + c2(ip)*M + c3(ip)*M*M + (c4(ip)+c5(ip)*M)*f1 +
     1      (c6(ip)+c7(ip)*M)*f2 + (c8(ip)+c9(ip)*M)*f0 + c10(ip)*R +sf2 + S
      if (ip.lt.26) then
        gnd = gnd*sfac - gfac
      else
c pgv?
        gnd = gnd*sfac
      endif
          if(Frq(ip).eq.0.) then
c  limit pga median to 1.5 g; 5hz to 3 g. 2006.
           gnd=min(gnd,0.405)
          elseif(Frq(ip).gt.2.5)then
           gnd=min(1.099,gnd)
          endif
       gndout(1)=gnd
c We don't plan to include extra gnd uncert for CEUS models but this is permitted in AB06 and others.
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
c      if(R.lt.10.)write(6,*) R,exp(gnd),M,sigma,exp(sfac*S),ip
      return
      end subroutine getAB06

      subroutine getTP05(iper,ip,ir,xmag0,dist0,gndout,sigma,sigmaf)
c Added getTP05, Oliver Boyd
c ip = index of period in Pd array and other coeff arrays.
c for what range of vs30 is this one supposed to be valid? SH nov 13 2006
c ir seems to have a hard-rock, firm-rock feature
c clamp on median motion added dec 7 2007
      parameter (np=16,sqrt2=1.414213562)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/ceus_sig/lceus_sigma,ceus_sigma
      real gnd_ep(3,3,8),gndout(3)
      logical lceus_sigma,l_gnd_ep(8)
c median limits are applied in main, not in subroutines.
      real f1,f2,f3,R,Rrup,M,cor
      real,dimension(np):: Pd,c1,c1h,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16
      Pd = (/0.00e+00, 0.04, 0.05,8.00e-02,1.00e-01,1.50e-01,2.00e-01,3.00e-01,
     1 0.40,0.50,
     1       7.50e-01,1.0, 1.50, 2.0, 3.0, 4.0/)
      c1h = (/1.14e+00,2.2217764,1.82e+00,6.83e-01,8.69e-01,2.38e+00,-5.48e-01,-5.13e-01,
     1 0.18295555,
     1       2.40e-01,-6.79e-01,-1.55e+00,-2.30e+00,-2.70e+00,-2.42e+00,-3.69e+00/)
c c1 below is based on a CEUS A->BC conversion from c1h where c1h Vs30 is ? and c1 is ?. 
c Use of earlier models of siteamp for ceus. So for all periods we use Fr. 1996 terms.
c c1 modified at 0.1, 0.3, 0.5, and 2.0 s for Frankel ceus amp. mar 19 2007.
c c1 checked for pga, 1hz and 5hz apr 17 2007. c1(0.4s) added June 30. corrected
c c1(0.3 s) to 0.0293 from K Campbell email Oct 13 2009.
      c1 = (/1.559e+00,2.6819708,2.24e+00,1.10e+00,1.4229,2.87e+00,1.70e-02,0.0293,
     1 0.9333664,
     1       6.974e-01,-3.23e-01,-1.257e+00,-1.94e+00,-2.5177,-2.28,-2.28/)
      c2 = (/6.23e-01,0.42201685,5.33e-01,7.43e-01,6.07e-01,5.01e-01,8.57e-01,6.67e-01,
     1 0.5887256,
     1       6.11e-01,6.66e-01,7.64e-01,7.94e-01,8.05e-01,8.01e-01,8.17e-01/)
      c3 = (/-4.83e-02,-0.05892166,-4.75e-02,-2.93e-02,-4.74e-02,-6.42e-02,-2.62e-02,
     1 -4.43e-02,-0.068179324,
     1       -7.89e-02,-8.30e-02,-8.59e-02,-8.84e-02,-9.29e-02,
     1       -1.08e-01,-1.18e-01/)
c c4(2.5hz) spline
      c4 = (/-1.81e+00,-1.5474958,-1.63e+00,-1.71e+00,-1.52e+00,-1.73e+00,-1.68e+00,
     1 -1.42,-1.4782926,
     1       -1.55e+00,-1.48e+00,-1.49e+00,-1.45e+00,-1.44e+00,
     1       -1.65e+00,-1.46e+00/)
      c5 = (/-6.52e-01,-0.46962538,-5.67e-01,-7.56e-01,-7.04e-01,-9.76e-01,-8.61e-01,
     1 -4.70e-01,-0.6740682,
     1       -8.44e-01,-7.34e-01,-9.41e-01,-8.86e-01,-9.23e-01,
     1       -8.98e-01,-8.45e-01/)
      c6 = (/4.46e-01,0.44941282,4.54e-01,4.60e-01,4.49e-01,4.14e-01,4.33e-01,4.68e-01,
     1 0.4376,0.43636485,
     1       4.35e-01,4.24e-01,4.12e-01,4.08e-01,4.37e-01,4.25e-01/)
      c7 = (/-2.93e-05,9.11331E-3,7.77e-03,-9.68e-04,-6.19e-03,6.60e-03,2.79e-03,1.08e-02,
     1  9.212394E-3,7.89e-03,
     1      9.53e-03,-5.84e-03,8.30e-03,2.06e-02,1.67e-02,1.13e-02/)
c c8 and all other coeffs at 2.5 hz use cubic spline interpolations. See intTP05 code. 
      c8 = (/-4.05e-03,-4.7678155E-3,-4.91e-03,-4.94e-03,-4.70e-03,-4.80e-03,-3.65e-03,
     1 -5.41e-03, -4.63148E-3,
     1       -3.65e-03,-3.37e-03,-2.09e-03,-3.27e-03,-2.14e-03,
     1       -2.03e-03,-1.72e-03/)
      c9 = (/9.46e-03,-1.5620958E-3,-3.14e-03,-5.50e-03,-4.24e-03,3.93e-03,-2.02e-03,6.44e-03,
     1 4.467385E-3,
     1       -2.65e-04,-1.19e-03,3.30e-03,2.51e-03,2.30e-03,3.58e-03,-3.34e-03/)
      c10 = (/1.41e+00,0.90152883,9.80e-01,1.13e+00,1.04e+00,1.51e+00,1.64e+00,1.52e+00,
     1 1.543476,
     1       1.59e+00,1.55e+00,1.52e+00,1.71e+00,1.43e+00,1.93e+00,1.69e+00/)
c c11 c10 and c9 log(f) interp 
      c11 = (/-9.61e-01,-0.95111495,-9.39e-01,-9.16e-01,-9.13e-01,-8.65e-01,-9.25e-01,
     1 -9.15e-01,-0.8872426,
     1       -8.59e-01,-7.84e-01,-7.57e-01,-7.69e-01,-7.55e-01,
     1       -8.18e-01,-7.37e-01/)
      c12 = (/4.32e-04,4.995133E-4,5.12e-04,4.82e-04,4.11e-04,3.64e-04,1.61e-04,4.32e-04,
     1 3.839152E-4,
     1       2.77e-04,2.45e-04,1.17e-04,2.33e-04,2.14e-04,1.16e-04,1.10e-04/)
c 
      c13 = (/1.33e-04,8.605951E-4,9.30e-04,7.33e-04,3.58e-04,6.84e-04,6.43e-04,2.87e-04,
     1 1.2648004E-4,
     1       1.46e-04,5.47e-04,7.59e-04,1.66e-04,3.91e-04,3.98e-04,3.59e-04/)
      c14 = (/1.21e+00,1.221863,1.22e+00,1.22e+00,1.23e+00,1.24e+00,1.24e+00,1.26e+00,
     1 1.2742121,
     1       1.28e+00,1.28e+00,1.28e+00,1.27e+00,1.26e+00,1.26e+00,1.25e+00/)
c for c15, corrected value for the 0.5-s or 2 Hz motion, from email Pezeshk dec 7 2007
      c15 = (/-1.11e-01,-0.10814471,-1.08e-01,-1.08e-01,-1.08e-01,-1.08e-01,-1.08e-01,
     1       -1.09e-01,-0.10839832,
     1 -1.073e-01,-1.05e-01,-1.03e-01,-9.99e-02,-9.78e-02,
     1       -9.52e-02,-9.26e-02/)
      c16 = (/4.09e-01,0.43780863,4.41e-01,4.49e-01,4.56e-01,4.64e-01,4.69e-01,4.79e-01,
     1 0.49330785,5.05e-01,5.22e-01,5.37e-01,5.51e-01,5.62e-01,5.73e-01,5.89e-01/)
      M = xmag0
      Rrup = dist0
      if(lceus_sigma)then
      sigma=ceus_sigma      !special study 3/2011.
      elseif (M.lt.7.2) then
        sigma = c14(ip) + c15(ip)*M 
      else
        sigma = c16(ip)
      endif
      sigmaf = 1./sqrt2/sigma
      if (ir.lt.0) then
        f1 = c1h(ip) + c2(ip)*M + c3(ip)*(8.5 - M)**2.5
      else
        f1 = c1(ip) + c2(ip)*M + c3(ip)*(8.5 - M)**2.5
      endif
      f2 = c9(ip)*log(Rrup + 4.5)
      if (Rrup.gt.70.) f2 = f2 + c10(ip)*log(Rrup/70.)
      if (Rrup.gt.130.) f2 = f2 + c11(ip)*log(Rrup/130.)
      cor = exp(c6(ip)*M + c7(ip)*(8.5 - M)**2.5)
      R = sqrt(Rrup*Rrup + c5(ip)*c5(ip)*cor*cor)
      f3 = (c4(ip) + c13(ip)*M)*log(R) + (c8(ip) + c12(ip)*M)*R
      gnd = f1 + f2 + f3
c---following is for clipping gnd motions: 1.5g PGA, 3.0g for T<0.5 s (2hz)
c this clamp was omitted until dec 7 2007. 
          if(pd(ip).lt.0.02) then
c 1.5 g is the PGA  limit (g). SH June 2006.
           gnd=min(gnd,0.405)
          elseif(pd(ip).lt.0.5)then
           gnd=min(1.099,gnd)
          endif
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
c      write(6,*) gnd,sigmaf,M,R,S,ip
      return
      end subroutine getTP05

      subroutine kanno(iper,ip,gndout,sigmaf,amag,rrup,vsfac)
c....................................................................
c  Kanno et al. regression relation, 2006 for shallow earthquakes
c  written by Yuehua Zeng, USGS. Mods Steve Harmsen Nov 7 2006.
c
c  input: per  : period
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
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
       real gnd_ep(3,3,8),gndout(3)
       logical l_gnd_ep(8)

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
       gndout(1)=gnd
         if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
         endif
       sigmaf=1.0/sigma/sqrt2
      return
      end subroutine kanno

        subroutine SomerIMW(ip,iprd, xmag,R_Rup,rjb,V_S30,gndout,sigmaf)
c Somerville et al prepared this model for USGS from finite-fault gm simulations.
c Coding by Stephen Harmsen. This model was designed for fault hazard in
c intermtn west. Should not be used in compressional regimes.
        integer mxnprd
        parameter (mxnprd=5)
      parameter (pi=3.14159265,sqrt2=1.414213562)
      parameter(Hsq = 42.25)      !H=6.5 in Somerville document
      common/mech/ss,rev,normal,obl
c Outputs gnd=ln(median) and sigmaf, two scalar quantities needed for Pex calcs.
c To DO: Find out if we want to include a site_amp feature with this subroutine
c Predictor variables, xmag = moment mag, R_Rup=nearest distance to the 
c fault plane (rupture surface), and rjb=joyner-boore dist. 
c Model requires no dip, no depth to top information. Hmm.
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)

        real  xmag, Width, R_rup,  V_S30, rjb, hw,gnd
      logical ss,rev,normal,obl
c Model coefficients
        real, dimension(mxnprd) :: prd,c1,c2,c3,c4,c5,c6,c7,c9,sigma_t
c prd array not used below. Assume SA 0.01 s is PGA.
c the coefficients
       prd=(/0.010,0.200,0.300,1.000,5.000/)
       c1=(/6.764,7.305,7.183,6.594,-4.598/)
      c2=(/-0.758,-0.758,-0.758,-0.758,0.633/)
      c3=(/-1.8,-1.597,-1.515,-1.480,-0.827/)
      c4=(/0.1375,0.1202,0.1156,0.0954,0.0005/)
      c5=(/-0.0104,-0.0104,-0.0104,-0.0055,-0.0018/)
      c6=(/-0.236,-0.236,-0.236,-0.351,-0.351/)
      c7=(/0.212,0.212,0.212,0.212,0.085/)
      sigma_t=(/0.6016,0.6073,0.6027,0.6270,0.8241/)
c the terms
      rc=sqrt(R_rup**2+Hsq)
      alogR=alog(rc)
      gnd=c1(iprd)+c2(iprd)*xmag+(c3(iprd)+c4(iprd)*xmag)*alogR
      gnd=gnd+c5(iprd)*R_rup+c6(iprd)*(8.5-xmag)**2
      sigma=sigma_t(iprd)
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
       sigmaf=1.0/sigma/sqrt2
c if strike slip or fw, do not include hanging wall effects.
      if(ss.or.rjb.gt.0.5.or.R_Rup.ge.20.)return
c Otherwise boost median with hanging wall term. may need to taper off edges      
      if(R_Rup.lt.5.)then
      hw=0.2*R_Rup
      elseif(R_Rup.lt.15.)then
      hw=1.0
      else
      hw=1.0-(R_Rup-15.0)/5.
      endif
      gnd=gnd+c7(iprd)*hw*(8.5-xmag)
       gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
      return
      end subroutine SomerIMW

      subroutine getSilva(iper,iq,ir,xmag,dist,gndout,sig,sigmaf)
c Inputs: xmag = moment mag, dist = rjb (km), see Nov 1, 2002 doc page 4
c iper = period index in iatten(ip,ia) and elsewhere.
c iq is index in pd array for current spectral period, iq shuld be in 1 to 8 range.
c ir is hardrock (-1) or BCrock indicator (+1).
c Silva (2002) table 5, p. 13. Single corner, constant stress drop, w/saturation.
c parts of this subroutine are from frankel hazFXv7 code. S Harmsen
c Added 2.5 and 25 hz July 1 2008 (for NRC projects). Add 20 hz, or 0.05 s. july 6
c Note: PGV and several other periods are available in Silva Table 5. 
      parameter (sqrt2=1.414213562,npmx=11)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/ceus_sig/lceus_sigma,ceus_sigma
        real gnd_ep(3,3,8),gndout(3),c
        logical lceus_sigma,l_gnd_ep(8)

      real, dimension(npmx) :: c1,c1hr,c2,c3,c4,c5,c6,c7,c8,c9,c10,pd,sigma
      pd=(/0.,0.04,0.05,0.1,0.2,0.3,0.4,0.5,1.,2.,5./)
      c1hr=(/5.53459,6.81012,6.63937,5.43782,3.71953,2.60689,
     1 1.64228,0.69539,-2.89906,-7.42051,-13.69697/)
c c1 from c1hr using A->BC factors, 1.74 for 0.1s, 1.72 for 0.3s, 1.58 for 0.5s, and 1.20 for 2s
c this from A Frankel advice, Mar 14 2007. For 25 hz use PGA amp. 
c For BC at 2.5 hz use interp between .3 and .5. 1.64116 whose log is 0.4953 
        c1=(/5.9533, 7.2288,7.023,5.9917,4.2848,3.14919,2.13759,
     1 1.15279,-2.60639,-7.23821,-13.39/)
       c2=(/-.11691,-0.13594,-.12193,-0.02059,.12490,.23165,.34751,
     1 .45254,.88116,1.41946,2.03488/)
      c4=(/ 2.9,3.0,3.,2.9,2.8,2.8,2.8,2.8,2.8,2.7,2.5/)     
      c6=(/ -3.42173,-3.48499,-3.45478,-3.25499,-3.04591,-2.96321,-2.87774,
     1 -2.818,-2.58296,-2.26433,-1.91969/)
      c7=(/ .26461,.26220,0.26008,0.24527,.22877,.22112,.21215,
     1 .20613,.18098,.14984,.12052/)
      c10=(/ -.06810,-.06115,-.06201,-0.06853,-.08886,-0.11352,-.13838,
     1 -0.16423,-.25757,-0.33999,-0.35463/)
c note very high sigma for longer period SA:
      sigma= (/.8471,0.8870,0.8753,0.8546,.8338,0.8428,.8386,
     1 0.8484,.8785,1.0142,1.2253/)
        if(ir.ge.1)then
        c=c1(iq)
        else
        c=c1hr(iq)
        endif
c formula (3) in Nov 1, 2002 doc
        gnd= c + c2(iq)*xmag+(c6(iq)+c7(iq)*xmag)*alog(dist+exp(c4(iq)))
     &  + c10(iq)*((xmag-6.0)**2)  
c---following is for clipping gnd motions: 1.5g PGA, 3.0g for T<0.5 s (2hz)
c this clamp was omitted until dec 7 2007. 
          if(pd(iq).lt.0.02) then
c  the PGA  limit =1.5 g. SH June 2006.
           gnd=min(gnd,0.405)
          elseif(pd(iq).lt.0.5)then
c 2.5, 3 and 5 hz 10 Hz median limit 3.0 g
           gnd=min(1.099,gnd)
          endif
       gndout(1)=gnd
       if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
       endif
c--- return sigma and sigmaf for clipping 
      if(lceus_sigma)then
      sig=ceus_sigma	!special study 3/2011
      else 
      sig=sigma(iq)
      endif
      sigmaf= 1.0/(sig*sqrt2)
      return
      end subroutine getSilva

      subroutine getSilvaV(iper,iq,ir,xmag,dist,gndout,sig,sigmaf)
c Inputs: xmag = moment mag, dist = rjb (km), see Nov 1, 2002 doc page 4
c Index Inputs:
c iper = period index in iatten(ip,ia) and elsewhere.
c iq is index in pd array for current spectral period, iq shuld be in 1 to 8 range.
c ir is hardrock (any integer but 1) or BCrock indicator (+1). BC not implemented
c Silva table 3c.  S Harmsen Oct 2012
c  Single corner, variable stress drop as fcn of M, w/saturation.
      parameter (sqrt2=1.414213562,npmx=8)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      common/ceus_sig/lceus_sigma,ceus_sigma
        real gnd_ep(3,3,8),gndout(3),c
        logical lceus_sigma,l_gnd_ep(8)
      real, dimension(npmx) :: c1,c1hr,c2,c3,c4,c5,c6,c7,c8,c9,c10,pd,sigma
      pd=(/0.,0.1,0.2,0.3,0.5,1.,2.,5./)
      c1hr=(/5.19757,4.94207,3.46126,1.72806,-0.51056,-4.46472,-9.12315,-15.20886/)
c A->BC factors: 1.52 for pga 1.74 for 0.1s; 1.76 for 0.2s; 1.58 for 0.5s; 1.34 for 1.0s; 1.20 for 2s, 1.15 for 5-s SA.
      c1=(/5.61628, 5.4959551,4.0265738,2.2703843,-0.0531352,-4.1720504,-8.9408284,-15.0690981/)	!copied from getSilva. C1uses 
c frankel A->BC factors. These may need work in future updates (say, relative to 2008)
      c2=(/0.07129,0.18877,0.33544,0.49782,0.76645,1.24156,1.78482,2.34990/)
      c4=(/ 2.8,2.8,2.8,2.7,2.7,2.7,2.6,2.4/)     
      c6=(/ -3.13247,-2.93847,-2.80186,-2.64087,-2.48971,-2.25138,-1.94059,-1.62679/)
      c7=(/ 0.20958,0.18204,0.16992,0.15740,0.14173,0.11635,0.08684,0.06230/)
      c10=(/ -.07375,-0.07672,-0.10713,-0.14160,-0.20549,-0.30377,-0.37239,-0.35359/)
c Note: very high sigma for longer period SA aleatory:
      sigma= (/.7405,0.7549,0.7289,0.7441,0.7601,0.8032,0.9557,1.1928/)
     
c clamping  to be done in main not in subroutines.
        if(ir.eq.1)then
        c=c1(iq)
c      stop 'now prepard for BC rock. sorry'
        else
        c=c1hr(iq)
        endif
c formula (3) in Nov 1, 2002 doc
        gnd= c + c2(iq)*xmag+(c6(iq)+c7(iq)*xmag)*alog(dist+exp(c4(iq)))
     &  + c10(iq)*((xmag-6.0)**2)  
c---following is for clipping gnd motions: 1.5g PGA, 3.0g for T<0.5 s (2hz)
c this clamp is applied inside the subroutine in the hazFX code. 
          if(pd(iq).lt.0.02) then
c  the PGA  limit =1.5 g. SH June 2006.
           gnd=min(gnd,0.405)
          elseif(pd(iq).lt.0.5)then
c 2.5, 3 and 5 hz 10 Hz median limit 3.0 g
           gnd=min(1.099,gnd)
          endif
       gndout(1)=gnd
       if(l_gnd_ep(iper))then
         gndout(2)= gnd+gnd_ep(ide,ime,iper)
         gndout(3)= gnd-gnd_ep(ide,ime,iper)
       endif
c--- return sigma and sigmaf for clipping 
      sig=sigma(iq)
      sigmaf= 1.0/(sig*sqrt2)
      sigmasq = sig*sqrt2
      return
      end subroutine getSilvaV



      subroutine resample(ift,npts,xlen,xaz,dip,depth0,nleny,
     + tlen,dlen,dlen2,npts1)
c this subroutine resamples fault trace to new 1km(?) discretization. It also
c determines downdip coordintes corresponding to those at top of fault.
c these downdip values will be used for buried ruptures with GR distribution.
c Inputs: ift = fault number
	parameter (nfltmx=1000)
c      npts= number of points on discretized input fault
c      x(*,ift)=x -coord
c      y(*,ift)=y-coord (degrees)
c      xlen = segment length, computed in main
c      xaz = segment azimuth, computed in main.
c      dip = dip angle (radians) pi/2>=dip>0      !note units.
c       dip needs to be positive for proper function.
c      xxx width = fault width (not input, not needed here).
c      depth0 = depth of top of fault (km) (not needed here)
c      nleny = number of downdip tops (4 max) for computing fault contours, including top
c      tlen = flt length (km) already computed. Use this value
c      dlen = delta-length for horizontal sampling; probably 1 km
c      dlen2 = vertical distance between downdip contours; probably 2 km
c      width = fault width (km)
c Outputs:
c      npts1 = number of points on resampled fault's contour line
c      u(ix,iy,ift) = x-coord, resampled
c      v(ix,iy,ift) = y-coord, resampled. iy is depth index (1,2,3, maybe 4)
c      
      parameter (d2r=3.14159265/180.0)
        real delta,depth0,dip,sinedip,cosdip,avaz,avbaz
        real, dimension(800) :: xaz, xlen
      real, dimension (800,nfltmx):: x,y
      real, dimension (990,4,nfltmx) :: u,v      !2nd dim is downdip 1=top, 2=down 2km,
       common/fault/x,y,u,v
c 3=down 4km. 4 = down 6 km (contour 4 first tried July 2007), 3rd dim fault index
      integer  npts,npts1
            nseg= npts - 1
      sinedip = sin(dip)
      cosdip= sqrt(1.-sinedip**2)
      if(abs(dip).gt.0.1)then
      ctgdip=cosdip/sinedip
      else
c Can code handle very shallow dipping faults? not with Plan A.
      ctgdip=3.0      !avoid sing
      endif
      n=npts
       call delaz(y(1,ift),x(1,ift),y(n,ift),x(n,ift),delta,avaz,avbaz)
       if(dip.gt.0.)then
       avdip=avaz+90.
       else
       avdip=avaz-90.
       endif
      nlenx=int((tlen+1.)/dlen)
      npts1=nlenx
      write(6,*)nlenx,nleny,avaz
      do 101 ix=1,nlenx
      delta= float(ix-1)*dlen
      do 500 iseg=1,nseg
      if(delta.lt.xlen(iseg)) go to 501
 500  continue
 501  continue
      if(iseg.ge.2) delta= delta-xlen(iseg-1)
      az= xaz(iseg)
      sx= x(iseg,ift)
      sy= y(iseg,ift)
c      write(6,*) "ix, iseg, az", ix, iseg, az
      call back(delta,az,sy,sx,dfy,dfx)
      u(ix,1,ift)=dfx
      v(ix,1,ift)=dfy
      if(dip.lt.0.) az2= az-90.
      if(dip.gt.0.) az2= az+90.
c---loop through increments iy downdip
      do 111 iy=2,nleny
      deltaz= (iy-1)*dlen2
      depth= depth0+ deltaz
c deltax for a given deltaz
      deltax= deltaz * ctgdip
c---find surface point above downdip portion of fault. We use an average dip
c rather than a local dip direction. Works better for many faults.
c
      if(dip.lt.1.57)then
      call back(deltax,avdip,dfy,dfx,dfy2,dfx2)
      u(ix,iy,ift)=dfx2
      v(ix,iy,ift)=dfy2
      else
c vertical dipping faults.
      u(ix,iy,ift)=u(ix,1,ift)
      v(ix,iy,ift)=v(ix,1,ift)
      endif
111      continue
101      continue
        u(npts1,1,ift)=x(n,ift)
        v(npts1,1,ift)=y(n,ift)      !same flt end point
      return
      end subroutine resample
      
      logical function LXYIN (X,Y,PX,PY,N)
c Is point (X,Y) inside the (unclosed) polygon (PX,PY)?
c See "Application of the winding-number algorithm...",
c  by Godkin and Pulli, BSSA, pp. 1845-1848, 1984.
c LXYIN= .true. if (X,Y) is inside or on polygon, .false. otherwise.
c Written by C. Mueller, USGS.
      dimension PX(N),PY(N)
      LXYIN= .true.
      KSUM= 0
      do 1 I=1,N-1
        K= KPSCR(PX(I)-X,PY(I)-Y,PX(I+1)-X,PY(I+1)-Y)
        if (K.eq.4) return
        KSUM= KSUM+K
1       continue
      K= KPSCR(PX(N)-X,PY(N)-Y,PX(1)-X,PY(1)-Y)
      if (K.eq.4) return
      KSUM= KSUM+K
      if (KSUM.eq.0) LXYIN= .false.
      return
      end function LXYIN

      integer function KPSCR (X1,Y1,X2,Y2)
c Compute the signed crossing number of the segment from (X1,Y1) to (X2,Y2).
c See "Application of the winding-number algorithm...",
c  by Godkin and Pulli, BSSA, pp. 1845-1848, 1984.
c KPSCR= +4 if segment passes through the origin
c        +2 if segment crosses -x axis from below
c        +1 if segment ends on -x axis from below or starts up from -x axis
c         0 if no crossing
c        -1 if segment ends on -x axis from above or starts down from -x axis
c        -2 if segment crosses -x axis from above
c Written by C. Mueller, USGS.
      KPSCR= 0
      if (Y1*Y2.gt.0.) return
      if (X1*Y2.eq.X2*Y1.and.X1*X2.le.0.) then
        KPSCR= 4
      elseif (Y1*Y2.lt.0.) then
        if (Y1.lt.0..and.X1*Y2.lt.X2*Y1) KPSCR= +2
        if (Y1.gt.0..and.X1*Y2.gt.X2*Y1) KPSCR= -2
      elseif (Y1.eq.0..and.X1.lt.0.) then
        if (Y2.lt.0.) KPSCR= -1
        if (Y2.gt.0.) KPSCR= +1
      elseif (Y2.eq.0..and.X2.lt.0.) then
        if (Y1.lt.0.) KPSCR= +1
        if (Y1.gt.0.) KPSCR= -1
      endif
      return
      end function KPSCR

       subroutine cluster_me(prob,prob_s,wtscene,nper,nlev,nscene,ngroup,j1,j2,isite)
c this routine computes prob. of exceedance for clustered events grouped into ngroup groups
c from rate of exceedance for individual sources stored in prob_s in input
c Sources are combined, zero or one from jseg(1), one from jseg(2), and zero or one from jseg(3), 
c and for a variety of M and a given (constant) return time, 
c nsrc has to be 2 or 3 because subr is written for exceedance from one or more of 2 or 3 clustered events.
c We want this to be more general, nsrc = n clustered events
c inputs:
c      nscene = number of scenarios per group. Probably 4 or 8 for NMSZ. 
c      nlev = number of gm levels. A vector, i=1,...,nper
c      isite = index for site (1 to 1000)
c      prob_s = probability of ground motion exceedance from one event in a
c group setting.
c  j1 = minimum segment number, 1 or 2
c  j2 = maximum segment number, 2 or 3. (this needs work, if two segments,  their numbers must be consecutive)
c
c There are nper periods, nlev(ip) ground motion levels per period, ngroup=1 to 5 groups
c This subroutine was written to solve the current model space for NMSZ only. There could easily be
c different sets to consider for other sources, but because of time constraints and limits to my imagination,
c I kept the solution focussed on this source. Steve Harmsen, May 17 2007
      parameter (nlvmx=20, npmx=8)  
      integer, dimension(npmx) :: ngroup,nlev
      integer i_anom/0/      !keep tract of anomalous rates. hope this is zero
      real   prob(1000,20,npmx,3,5)
      real*8 p,q(4)
       real prob_s(3,8,5,20,8), wtscene(8,5)
c prob_s: dim1= segment #, 1 to 3, dim2=ks=scenario, dim3=kg=group, dim4=gm level, dim5=period indx
c convert for groups, for scenarios, for each gm level, for each spectral period
       nseg_c=j2-j1+1
      do ip=1,nper
      do igp=1,ngroup(ip)
      do k=1,nlev(ip)
c ps stores sum of scenario rates
      ps=0.0
      do j=1,nscene
      p=0.0
      i=1
      do isrc=j1,j2      !nsrc=2 or 3
      q(i)=prob_s(isrc,j,igp,k,ip)
      p=p+q(i)
      i=i+1
      enddo
      if(nseg_c.eq.3)then
      p= p +q(1)*q(2)*q(3) -q(1)*q(2) -q(1)*q(3) -q(2)*q(3)
      else 
      p= p-q(1)*q(2)
      endif
c above includes all cross terms for this limited model
      if(p.ge.1.0001)then
c Flow should not arrive here, but if so, assume 10 events a year. That should be bad enough
      ps =10.*wtscene(j,igp) + ps
      i_anom=i_anom+1
      write(*,*)q
      write(*,*)p,isite,j,k,igp,ip
      if(i_anom.gt.10)stop
      else
      ps = ps +wtscene(j,igp)*p
      endif
      enddo      !scenarios 8?
      prob(isite,k,ip,1,igp)=ps
      enddo      !nlev gm levels
      enddo      !ngroup 5 or less
      enddo      !spectral period up to 7
      if(i_anom.gt.0)then
      write(*,*)'Cluster_me had ',i_anom, ' high rates. Should be 0'
      write(*,*)'This outcome should be examined further '
      endif
      return
      end subroutine cluster_me
      
      subroutine mindist(dist,azm,disb,azb,rx,stlat,stlon,solat,solon,ns,
     +                   dip,wid,ftop)
c.....................................................................
      parameter (d2r = 0.0174533)
c  input:
c         stlat - station latitude
c         stlon - station longitude
c         solat - latitude of fault segment corners
c         solon - longitude of fault segment corners
c         ns    - number of fault segment corners
c         dip   - fault dip angle
c         wid   - fault width
c         ftop  - depth of fault top 
c  output:
c         dist  - minimum distance to rutpure plane
c         azm   - azimuth for rupture distance
c         disb  - JB distance to rutpure plane
c         azm   - azimuth for JB distance
c         rx    - distance to the surface projection of top edge of fault
c                                              Yuehua Zeng, 2007
c.....................................................................
      dimension solat(500),solon(500)
c
      dp=dip*d2r
      cosdp=cos(dp)
      sindp=sin(dp)
      dist=1.e10
      disb=1.e10
      call delaz2(solat(1),solon(1),solat(ns),solon(ns),flen,str1)
      call delaz2(solat(1),solon(1),stlat,stlon,dista,az)
      sx=dista*cos(az-str1)
      y=dista*sin(az-str1)
      if(sx.le.0.0.or.sx.ge.flen)then; rx=y
      else; rx=sign(sqrt(dista*dista+flen*(flen-sx-sx)),y); endif
      do i=1,ns-1
         call delaz2(solat(i),solon(i),solat(i+1),solon(i+1),flen,str)
         call delaz2(solat(i),solon(i),stlat,stlon,dista,az)
         sx=dista*cos(az-str)
         if(sx.gt.flen)then
           sx=sx-flen
         elseif(sx.gt.0.0)then
           sx=0.0
         endif
         y=dista*sin(az-str)
         r=sqrt(sx*sx+y*y)
         if(sx.le.0.0.and.r.lt.abs(rx))rx=sign(r,y)
c
c  for JB distance
         sy=y
         if(sy.gt.wid*cosdp)then
           sy=sy-wid*cosdp
         elseif(sy.gt.0.0)then
           sy=0.0
         endif
         r=sqrt(sx*sx+sy*sy)
         if(r.lt.disb)then
           disb=r
           azb=az/d2r
         endif
c
c  for rupture distance
         sy=y*cosdp-ftop*sindp
         if(sy.gt.wid)then
           sy=sy-wid
         elseif(sy.gt.0.0)then
           sy=0.0
         endif
         sz=y*sindp+ftop*cosdp
         r=sqrt(sx*sx+sy*sy+sz*sz)
         if(r.lt.dist)then
           dist=r
           azm=az/d2r
         endif
      end do
c
      return
      end subroutine mindist
c
      subroutine delaz2(sorlat,sorlon,stnlat,stnlon,delta,az)
      parameter (coef= 3.1415926/180.)
      xlat= sorlat*coef
      xlon= sorlon*coef
      st0= cos(xlat)
      ct0= sin(xlat) 
      phi0= xlon
      xlat= stnlat*coef
      xlon= stnlon*coef
      ct1= sin(xlat)
      st1= cos(xlat)
      sdlon= sin(xlon-phi0)
      cdlon= cos(xlon-phi0)
      cdelt= st0*st1*cdlon+ct0*ct1
      x= st0*ct1-st1*ct0*cdlon
      y= st1*sdlon
      sdelt= sqrt(x*x+y*y)
      delta= atan2(sdelt,cdelt)
      delta= delta/coef
      az= atan2(y,x)
c     az= az/coef
      delta= delta*111.2
      return
      end    subroutine delaz2

      subroutine mindist1(dist,azm,disb,azb,rx,stlat,stlon,solat,solon,ns,
     +                   dip,wid,ftop)
c.....................................................................
c  For calculating JB or rupture distance to a curved fault defined 
c  by extending the fault top segment trace down its dip direction
c
c  input:
c         stlat - station latitude
c         stlon - station longitude
c         solat - latitude of fault segment corners
c         solon - longitude of fault segment corners
c         ns    - number of fault segment corners
c         dip   - fault dip angle
c         wid   - fault width along the dip direction. wid>0 is a good constraint
c         ftop  - depth of fault top 
c  output:
c         dist  - minimum distance to rutpure plane
c         azm   - azimuth for rupture distance
c         disb  - JB distance to rutpure plane
c         azb   - azimuth for JB distance
c         rx    - distance to the surface projection of top edge of fault
c                                              Yuehua Zeng, 2007
c.....................................................................
      parameter (d2r = 3.1415926/180.)
      dimension solat(*),solon(*)
c
      dp=dip*d2r
      cosdp=cos(dp)
      sindp=sin(dp)
      dist=1.e10
      disb=1.e10
      call delaz2(solat(1),solon(1),solat(ns),solon(ns),flen,str1)
      call delaz2(solat(1),solon(1),stlat,stlon,dista,az)
      sx=dista*cos(az-str1)
      yy=dista*sin(az-str1)
      if(sx.le.0.0.or.sx.ge.flen)then; rx=yy
      else; rx=sign(sqrt(dista*dista+flen*(flen-sx-sx)),yy); endif
      do i=1,ns-1
         call delaz2(solat(i),solon(i),solat(i+1),solon(i+1),flen,str)
         call delaz2(solat(i),solon(i),stlat,stlon,dista,az)
c
         x=sindp
         y=cosdp*cos(str-str1)
         sindpp=x/sqrt(x*x+y*y)
         cosdpp=sqrt(1.0-sindpp*sindpp)
c
         sx=dista*cos(az-str)
         yy=dista*sin(az-str)
c
         if(sx.ge.0.0.and.sx.le.flen.and.abs(yy).lt.abs(rx))then;rx=yy
         elseif(sx.lt.0.0.and.dista.lt.abs(rx))then;rx=sign(dista,yy);endif
c
c  for JB distance
         sy=yy
         wid1=wid*cosdp
c
         cosfi=sin(str-str1)
         sinfi=sqrt(1.0-cosfi*cosfi)+1.0e-10
         syy=sy/sinfi
         sxx=sx-syy*cosfi
         r=-10.0
c
         if(sxx.lt.0.0.or.sxx.gt.flen)then
           y=sxx*cosfi+syy
           if(sxx.gt.flen)y=(sxx-flen)*cosfi+syy
           if(y.lt.0.0.or.y.gt.wid1)then
             if(y.gt.wid1)syy=syy-wid1
             x=syy*cosfi+sxx
             if(x.gt.flen)then
               sxx=sxx-flen
             elseif(x.ge.0.0)then
               r=abs(syy)*sinfi
               az0=90.0*sign(1.0,syy)
             endif
           else
             if(sxx.gt.flen)sxx=sxx-flen
             r=abs(sxx)*sinfi
             az0=(str1-str)/d2r-90.0+90.0*sign(1.0,sxx)
           endif
         elseif(syy.lt.0.0.or.syy.gt.wid1)then
           if(syy.gt.wid1)syy=syy-wid1
           x=syy*cosfi+sxx
           if(x.gt.flen)then
             sxx=sxx-flen
           elseif(x.ge.0.0)then
             r=abs(syy)*sinfi
             az0=90.0*sign(1.0,syy)
           endif
         else
           r=0.0
           az0=90.0
         endif
c
         if(r.lt.-1.0)then
           sxx=sxx+syy*cosfi
           syy=syy*sinfi
           r=sqrt(sxx*sxx+syy*syy)
           az0=atan2(syy,sxx)/d2r
         endif
         if(r.lt.disb)then
           disb=r
           azb=str/d2r+az0
         endif
c
c  for rupture distance
         sy=yy*cosdpp-ftop*sindpp
         sz=yy*sindpp+ftop*cosdpp
c
         cosfi=sin(str-str1)*cosdp
         sinfi=sqrt(1.0-cosfi*cosfi)+1.0e-10
         syy=sy/sinfi
         sxx=sx-syy*cosfi
         r=-10.0
c
         if(sxx.lt.0.0.or.sxx.gt.flen)then
           y=sxx*cosfi+syy
           if(sxx.gt.flen)y=(sxx-flen)*cosfi+syy
           if(y.lt.0.0.or.y.gt.wid)then
             if(y.gt.wid)syy=syy-wid
             x=syy*cosfi+sxx
             if(x.gt.flen)then
               sxx=sxx-flen
             elseif(x.ge.0.0)then
               r=abs(syy)*sinfi
               az0=90.0*sign(1.0,syy)
             endif
           else
             if(sxx.gt.flen)sxx=sxx-flen
             r=abs(sxx)*sinfi
             az0=(str1-str)/d2r-90.0+90.0*sign(1.0,sxx)
           endif
         elseif(syy.lt.0.0.or.syy.gt.wid)then
           if(syy.gt.wid)syy=syy-wid
           x=syy*cosfi+sxx
           if(x.gt.flen)then
             sxx=sxx-flen
           elseif(x.ge.0.0)then
             r=abs(syy)*sinfi
             az0=90.0*sign(1.0,syy)
           endif
         else
           r=0.0
           az0=90.0
         endif
c
         if(r.gt.-1.0)then
           r=sqrt(r*r+sz*sz)
         else
           sxx=sxx+syy*cosfi
           syy=syy*sinfi
           r=sqrt(sxx*sxx+syy*syy+sz*sz)
           az0=atan2(syy*cosdpp+sz*sindpp,sxx)/d2r
         endif
         if(r.lt.dist)then
           dist=r
           azm=str/d2r+az0
         endif
      end do
c
      return
      end subroutine mindist1

      subroutine mindist1n(dist,azm,disb,azb,rx,stlat,stlon,solat,solon,
     +                   ns,njump,dip,wid,ftop,rb)
c.....................................................................
c  For calculating JB or rupture distance to a curved fault defined 
c  by extending the fault top segment trace down its dip direction
c
c  input:
c         stlat - station latitude
c         stlon - station longitude
c         solat - latitude of fault segment corners
c         solon - longitude of fault segment corners
c         ns    - number of fault segment corners
c         njump - the non-rupture segment order number, if no non-rupture
c                 segment, set njump=0
c         dip   - fault dip angle
c         wid   - fault width along the dip direction. wid>0 is a good constraint
c         ftop  - depth of fault top 
c  output:
c         dist  - minimum distance to rutpure plane
c         azm   - azimuth for rupture distance
c         disb  - JB distance to rutpure plane
c         azb   - azimuth for JB distance
c         rx    - Chiou and Yong's distance to the surface projection 
c                 of top edge of fault
c         rb    - Bartlett's distance to the surface projection of 
c                 top edge of fault
c                                              Yuehua Zeng, 2007
c.....................................................................
      parameter (d2r = 3.1415926/180.)
      dimension solat(*),solon(*)
c
      dp=dip*d2r
      cosdp=cos(dp)
      sindp=sin(dp)
      dist=1.e10
      disb=1.e10
      call delaz2(solat(1),solon(1),solat(ns),solon(ns),flen,str1)
      call delaz2(solat(1),solon(1),stlat,stlon,dista,az)
      sx=dista*cos(az-str1)
      yy=dista*sin(az-str1)
      rb=sqrt(dista*dista+flen*(flen-sx-sx))
      if(sx.le.0.0.or.sx.ge.flen)then; rx=yy
      else; rx=sign(sqrt(dista*dista+flen*(flen-sx-sx)),yy); endif
      do i=1,ns-1
         if(i.ne.njump)then
         call delaz2(solat(i),solon(i),solat(i+1),solon(i+1),flen,str)
         call delaz2(solat(i),solon(i),stlat,stlon,dista,az)
c
         x=sindp
         y=cosdp*cos(str-str1)
         sindpp=x/sqrt(x*x+y*y)
         cosdpp=sqrt(1.0-sindpp*sindpp)
c
         sx=dista*cos(az-str)
         yy=dista*sin(az-str)
c
         if(sx.ge.0.0.and.sx.le.flen.and.abs(yy).lt.abs(rx))then;rx=yy
         elseif(dista.lt.abs(rx))then;rx=sign(dista,yy);endif
         if(sx.ge.0.0.and.sx.le.flen.and.abs(yy).lt.rb)then;rb=abs(yy)
         elseif(dista.lt.rb)then;rb=dista;endif
c
c  for JB distance
         sy=yy
         wid1=wid*cosdp
c
         cosfi=sin(str-str1)
         sinfi=sqrt(1.0-cosfi*cosfi)+1.0e-10
         syy=sy/sinfi
         sxx=sx-syy*cosfi
         r=-10.0
c
         if(sxx.lt.0.0.or.sxx.gt.flen)then
           y=sxx*cosfi+syy
           if(sxx.gt.flen)y=(sxx-flen)*cosfi+syy
           if(y.lt.0.0.or.y.gt.wid1)then
             if(y.gt.wid1)syy=syy-wid1
             x=syy*cosfi+sxx
             if(x.gt.flen)then
               sxx=sxx-flen
             elseif(x.ge.0.0)then
               r=abs(syy)*sinfi
               az0=90.0*sign(1.0,syy)
             endif
           else
             if(sxx.gt.flen)sxx=sxx-flen
             r=abs(sxx)*sinfi
             az0=(str1-str)/d2r-90.0+90.0*sign(1.0,sxx)
           endif
         elseif(syy.lt.0.0.or.syy.gt.wid1)then
           if(syy.gt.wid1)syy=syy-wid1
           x=syy*cosfi+sxx
           if(x.gt.flen)then
             sxx=sxx-flen
           elseif(x.ge.0.0)then
             r=abs(syy)*sinfi
             az0=90.0*sign(1.0,syy)
           endif
         else
           r=0.0
           az0=90.0
         endif
c
         if(r.lt.-1.0)then
           sxx=sxx+syy*cosfi
           syy=syy*sinfi
           r=sqrt(sxx*sxx+syy*syy)
           az0=atan2(syy,sxx)/d2r
         endif
         if(r.lt.disb)then
           disb=r
           azb=str/d2r+az0
         endif
c
c  for rupture distance
         sy=yy*cosdpp-ftop*sindpp
         sz=yy*sindpp+ftop*cosdpp
c
         cosfi=sin(str-str1)*cosdp
         sinfi=sqrt(1.0-cosfi*cosfi)+1.0e-10
         syy=sy/sinfi
         sxx=sx-syy*cosfi
         r=-10.0
c
         if(sxx.lt.0.0.or.sxx.gt.flen)then
           y=sxx*cosfi+syy
           if(sxx.gt.flen)y=(sxx-flen)*cosfi+syy
           if(y.lt.0.0.or.y.gt.wid)then
             if(y.gt.wid)syy=syy-wid
             x=syy*cosfi+sxx
             if(x.gt.flen)then
               sxx=sxx-flen
             elseif(x.ge.0.0)then
               r=abs(syy)*sinfi
               az0=90.0*sign(1.0,syy)
             endif
           else
             if(sxx.gt.flen)sxx=sxx-flen
             r=abs(sxx)*sinfi
             az0=(str1-str)/d2r-90.0+90.0*sign(1.0,sxx)
           endif
         elseif(syy.lt.0.0.or.syy.gt.wid)then
           if(syy.gt.wid)syy=syy-wid
           x=syy*cosfi+sxx
           if(x.gt.flen)then
             sxx=sxx-flen
           elseif(x.ge.0.0)then
             r=abs(syy)*sinfi
             az0=90.0*sign(1.0,syy)
           endif
         else
           r=0.0
           az0=90.0
         endif
c
         if(r.gt.-1.0)then
           r=sqrt(r*r+sz*sz)
         else
           sxx=sxx+syy*cosfi
           syy=syy*sinfi
           r=sqrt(sxx*sxx+syy*syy+sz*sz)
           az0=atan2(syy*cosdpp+sz*sindpp,sxx)/d2r
         endif
         if(r.lt.dist)then
           dist=r
           azm=str/d2r+az0
         endif
         endif
      end do
c
      return
      end subroutine mindist1n

      subroutine getDahle95(iq,isoil,xmag,dist0,gnd3,sigmaf)
c Inputs:
c iq is the period index.
c isoil=1 if soil, 0 if rock. We make the transition at 550 m/s.
c Dist0 is the Rjb distance(km). CHECK THIS.
c xmag is moment mag.
c Outputs:
c gnd = median motion (units g) logged.
c sigmaf = 1/sigma/sqrt2
c This version is for use in hazFXnga codes. 
c Fact to ponder: Dahle Regression had < 100 rock records and < 100
c      soil records. 
      parameter (npmx=6,emax=3.0,sqrt2=1.41421356,twopi=6.2831852)
       real omega,gnd
      real, dimension(3):: gnd3
      real, dimension(npmx):: c1,c2,c3,c4,c5,period,sigma
      period=(/0.0,0.1,0.2,0.5,1.,2.0/)
      c1=(/	-1.579,-4.608,-4.746,-5.717,-6.595,-7.205/)
      c2=(/0.554, 0.486,0.645, 0.920, 1.084,1.131/)    
      c3=(/-0.560, -0.609,-0.674,-0.761, -0.792, -0.762/)  
      c4=(/-0.00302,-0.00198, -0.00155,-0.00106,-0.00075, -0.00051/)  
      c5=(/0.326, 0.381,0.470,0.566,0.588, 0.536/)
      sigma=(/0.73, 0.79, 0.80,  0.81, 0.79,  0.75/)
c--  for Dahle et al. 1995 EERI Internatl Seis. Zonation Conf. relats
c      write(6,*) "enter c1,c2,c3,c4,c5,sigma"
c        if(icode(ip,ia).eq.1) xmag= 2.715 -0.277*xmag0+0.127*xmag0*xmag0
        delta= -7.333 + 2.333*xmag
        dist= sqrt(dist0*dist0 + 36.0)
          gnd= c1(iq)+c2(iq)*xmag+c3(iq)*alog(dist)
          gnd= gnd+ c4(iq)*dist + c5(iq)*isoil
c-- Convert from psv to psa (not for pga though)
      if(period(iq).ne.0.)then
          omega=twopi/period(iq)
           gnd=gnd+alog(omega)
           endif
c-- Convert from m/s*s to g.
          gnd=gnd-alog(9.81)
          gnd3(1)=gnd
          sigmaf=  1./sqrt2/sigma(iq)
      return
      end subroutine getDahle95

      real function basiteamp(pganl,vs30,vs30r)
c this version works with just one period of SA and one Vs30. Used in the deagg codes for 2009.
c inputs: pganl (g)
c      vs30: vs30 of site
c      vs30r: a reference vs30, usually one value for soil and another for rock
c implicit input:
c      ip = outer loop period index is 1 => period2: spectral period (scalar) is known
c common input
c       iq = period index in the below arrays. This has been precomputed
c      ipgeom (not used here)
c      lvs	(not used here)
c
c Output: basiteamp = log(AMP at vs30)-log(AMP at vs30r)
      parameter (np=23)	!23 periods apr 07. include 0.01 to 10 s 
      parameter (pi=3.14159265,sqrt2=1.414213562,vref=760.)
      parameter (dx=1.098612289,dxsq=1.206948961,dxcube=1.325968960,plfac=-0.510825624)
c dx = ln(a2/a1), made a param. used in a smoothed nonlin calculation sept 2006.
c plfac = ln(pga_low/0.1)      This never changes so shouldnt be calculated.
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
c        iq = iperba     !Iq comes in in common/ /
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
        dyr= bnl*a2fac 
        c=(3.*dy-bnl*dx)/dxsq
        d=(bnl*dx-2.*dy)/dxcube
        cr=(3.*dy-bnlr*dx)/dxsq
        dr=(bnlr*dx-2.*dy)/dxcube
        site0 = blin(iq)*alog(vs30/vref)
        site0r = blin(iq)*alog(vs30r/vref)
        compute = .false.
c        print *, bnl, bnlr, dy,dyr,site0,site0r
        endif      !first-pass compute
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
      subroutine GailTable(i)
c the "i" index refers to which model's table is read in. This is a mod
c from the snippet Gail sent, which only allows one model.
c moved open to main. problems passing the string not resolved.
c Currently looking at 3 possible models, AB08, AB06', or Pz11
c wrute ti unit 94.
        parameter (gfac=6.8875526,sfac=2.3025851)	!to convert to base e
      common/gail1/xmag(20,3), gma(20,30,20,3), rlog(30,3), f(20,3),itype(3)
      common/gail2/nf(3), nd(3), nm(3)
      common/gnome/fname
      logical m7
      real a,b,c
      character*80 header,fname
      if(i.gt.3)stop 'please redimension the arrays in gail1 and gail2.'
      read(3,3)header
3     format(a)
      i2=index(fname,'.dat')-13
      i3=i2+19
      read(3,*) nm(i), nd(i), nf(i), itype(i)
      write(*,9) fname, itype(i)
9     format(' GMfile, itype(0=Repi,1=Rhypo,2=Rcd,3=Rjb): ',/, a20,0p,i2)
99      format(a)
      read(3,*) (f(jf,i), jf=1,nf(i))
c
      write(94,9)fname(i2:i3), itype(i)
      write(94,99)'R(km)  1hzSA(g)   5hzSA(g)	  PGA(g) for hard rock'
      do 5000 jm = 1, nm(i)
        read(3,*)xmag(jm,i)
        m7=xmag(jm,i).eq.7.0
        do 4000 jd = 1, nd(i)
          read(3,*)rlog(jd,i),(gma(jf,jd,jm,i), jf = 1,nf(i))
          a=gma(4,jd,jm,i)*sfac-gfac;b=gma(7,jd,jm,i)*sfac-gfac;c=gma(12,jd,jm,i)*sfac-gfac
          if(m7)write(94,10)10**rlog(jd,i),exp(a),exp(b),exp(c)
10      format(f7.1,1p,3(1x,e11.5))
4000    continue
5000  continue
      close(3)
      return
      end subroutine GailTable

      real function amean11(amag,r,rjb,Vs30,ip,jf,ka)
c Inputs: amag: Mw
c       r= distance (rjb or rcd, depending on ka)
c       Vs30 = top 30 m Vs. A or BC is expected.
c      ip = global period index
c      jf = frequency index
c      ka = 1,2, or 3 depending on which table to use AB06 AB08 or P11
c     Gets ground motion value (ln PSA) from table
      real sfac,Vs30,r,rjb,amag
        parameter (gfac=6.8875526,sfac=2.3025851)	!to convert to base e
      common/gail1/xmag(20,3), gma(20,30,20,3), rlog(30,3), f(20,3), itype(3)
      common/gail2/nf(3), nd(3), nm(3)
      common/gail3/freq(13)
      real, dimension(13) :: bcfac
        real, dimension(20) :: avg
c frequencies: 0.2, 0.33, 0.5, 1.0, 2.0, 3.33, 5., 10., 20., 33.0, 50., PGA, PGV
            bcfac = (/0.06, 0.08,0.09, 0.11, 0.14, 0.14, 0.12, 0.03, -0.2, -.045, -.045, -.045,0.09/)
c bcfac for pga can be a function of rjb= -0.3+.15 log(rjb) (Repi in Gail's notes)
c     
c     Use interpolation to get amean for given freq(jfreq), amag, R (hazard looping values).
c     Note that freq(jfreq) is the defined block of freq. for hazard calcs (common block)
c     GMPE tables give values(gma) for freq array f(nf), for distances rlog(nd), magnitudes xmag(nm)
c
c The frequency is already known. It is input. Don't need to search.
      jfl=jf
      rl = alog10(r)
      jdflag = 0
c     Find bounding jf values.
c     check near the edges (in case match of frequencies not exact)
      jfu = jfl + 1
      fracf = 0.
      if (jf.gt.11) go to 15
      fracf = (alog10 (freq(ip)) - alog10 (f(jfl,ka)) ) /
     *        (alog10 (f(jfu,ka)) - alog10 (f(jfl,ka)) )
c**   fracf gives interpolation fraction in frequency
15    continue
c     Find bounding jm values.
      jml = 0
      do 20 j = 1, nm(ka)-1
       if (amag .ge. xmag(j,ka) .and. amag .lt. xmag(j+1,ka)) jml = j
20    continue
      if (abs(amag - xmag(nm(ka),ka)) .lt. 0.01) jml = nm(ka)
      if (jml .eq. 0) write(*,*)' ERROR. MAGNITUDE OUTSIDE RANGE.'
      jmu = jml + 1
      fracm = (amag - xmag(jml,ka)) / (xmag(jmu,ka)-xmag(jml,ka))
c**   fracm give interpolation fraction in magnitude
c     find bounding distance values.
      if (rl .lt. rlog(1,ka)) write(*,*)' ERROR. DIST. < RANGE '
        jdl = 0
      do 30 j = 1, nd(ka)-1
      if (rl .ge. rlog(j,ka) .and. rl .lt. rlog(j+1,ka)) jdl= j
30    continue
      if (rl .eq. rlog(nd(ka),ka)) jdl = nd(ka)
        if (rl .gt. rlog(nd(ka),ka)) then
c*       dist. beyond table maximum.  Attenuate ground motion by 1/R from last point in table.
         rscale= rl - rlog(nd(ka),ka)
         jdflag= 1
         jdl = nd(ka)
        endif
c        
      if (jdl .eq. 0) write(*,*)' ERROR. CANNOT FIND DIST'
      jdu = jdl + 1
      fracd = (rl - rlog(jdl,ka)) / (rlog(jdu,ka)-rlog(jdl,ka))
c**   fracd gives interpolation fraction in distanc
      do 40 j = jfl, jfu
      arl = gma(j,jdl,jml,ka) + fracm * (gma(j,jdl,jmu,ka)-gma(j,jdl,jml,ka))
      aru = gma(j,jdu,jml,ka) + fracm * (gma(j,jdu,jmu,ka)-gma(j,jdu,jml,ka))
c       arl is the amplitude for lower dist. bound, so arl>aru
      avg(j) = arl - fracd*(arl-aru)
40    continue
      amean = avg(jfl) + fracf * (avg(jfu) - avg(jfl))
      if (jdflag .eq. 1) amean = amean - rscale
c modify for BC rock. PGA BC lower than A at close-in sites according to AB11.
c PGA median is > 5hz median for A and even for BC for close-in sites.
c Current limitation : either A or BC. Nothing in between is provided for.
c Also, nothing in 2011 models is available for handing soil Vs.
      if(Vs30.lt.800.)then
      if(jf.eq.12.or.jf.eq.11)then
      amean = amean -0.3 + 0.15*rl
      else
      amean = amean + bcfac(jf)
      endif	!pga or not?
      endif	!BC rock from A?
c change to base e log to fit into standard framework. Convert to units g
      amean11 = amean*sfac -gfac
c apply the median clamp for some frequencies. Gail email, Mar 23, 2011.
      if(freq(jf).gt.2.1 .and. freq(jf) .lt.8.)then
      amean11=min(amean11,1.792)
      elseif(freq(jf).gt.90.)then
      amean11=min(amean11,0.4055)
      endif
      return
      end function  amean11

      real function sigPez11(mag,ip)
      parameter (fac=-6.95e-3,afac=2.3025851 )
      real mag	!Mw
c Input magnitude and period index, 
c Output sigma_lnY
c magnitude dependent sigma_lnY for the Pezeshk and Zandieh (BSSA,2011)
c GMPE. coeffs c12, c13,c14 from table 2 of their article
       common/ceus_sig/lceus_sigma,ceus_sigma
       logical lceus_sigma
      real, dimension(13) :: per, c12,c13,c14
      per= (/5.,3.,2.,1.,0.5,0.3,0.2,0.1,0.05,0.03,0.02,0.00,0.01/)	!s
      c12= (/-6.9e-3,-8.509e-3,-9.443e-3,-1.18e-2,-1.556e-2,-1.837e-2,-2.046e-2,-2.259e-2,
     & -2.244e-2,-2.094e-2,-1.974e-2,-2.105e-2,-1.974e-2/)
      c13= (/3.577e-1,3.54e-1,3.56e-1,3.588e-1,3.722e-1,3.867e-1,3.979e-1,4.102e-1,3.990e-1,
     & 3.817e-1,3.691e-1,3.778e-1,3.688e-1/)
      c14= (/3.58e-1,3.43e-1,3.387e-1,3.249e-1,3.119e-1,3.068e-1,3.033e-1,3.007e-1,2.905e-1,
     & 2.838e-1,2.796e-1,2.791e-1,2.792e-1/)
        if(lceus_sigma)then
        sigPez11=ceus_sigma      !for special studies. Use a low sigma for all GMPE periods, magnitudes, etc.
      elseif(mag.le.7.0) then
      sigPez11=c12(ip)*mag+c13(ip)
      else
      sigPez11=fac*mag+c14(ip)
      endif
c the Pezeshk article is in base10 logarithm units. we work in natural log units.
      sigPez11=sigPez11*afac
      return
      end function sigPez11

c-----------------------------------------------------------------------
        subroutine CY2012_NGA(ip,iprd, M, R_Rup, R_JB, R_x, V_S30,
     1                   F_Measured, F_Inferred, Z1, DELTA, Z_TOR,
     1                   F_RV, F_NM,
     1                   psa_ref, gndout, tau, sigma_NL0, sigmaf)

        implicit none
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
       common/fix_sigma/fix_sigma,sigma_fx
      real gnd_ep(3,3,8),gndout(3),sqrt2, sigmaf
      logical l_gnd_ep(8),fix_sigma
c Predictor variables
        real PERIOD, M, Width, R_rup, R_JB, R_x, V_S30, F_Measured,
     1       F_Inferred, Z1, DELTA, Z_TOR, F_RV, F_NM,sigma_fx

c Model cofficients
c ip is global period index, used to assign additional epistemic gm uncert.
        integer nprd,ip,ide,ime
        parameter (nprd=3,sqrt2=1.4142136)
        real prd(nprd),
     1       c1(nprd), c1a(nprd), c1b(nprd), c2(nprd),
     1       c3(nprd), cn(nprd),  cM(nprd),  c4(nprd),
     1       c4a(nprd),cRB(nprd), c5(nprd),  c6(nprd),
     1       cHM(nprd),c7(nprd),  c9(nprd),  c9a(nprd), c11(nprd),
     1       cgamma1(nprd), cgamma2(nprd), cgamma3(nprd),
     1       phi1(nprd), phi2(nprd), phi3(nprd), phi4(nprd),
     1       phi5(nprd), phi6(nprd), phi7(nprd), phi8(nprd),
     1       tau1(nprd), tau2(nprd),
     1       sigma1(nprd), sigma2(nprd), sigma3(nprd)
      data prd     /  0.0100, 0.2000, 1.0000/

      data c1      / -1.0944608,-0.4072448,-2.0178770/

      data c1a     /  0.1000, 0.1000, 0.0766/

      data c1b     / -0.2550,-0.2449,-0.1400/

      data c2      /  1.06,    1.06,   1.06/

      data c3      /  1.975122, 2.186490, 3.019458/

      data cn      /  2.9960, 2.8312, 1.6480/

      data cM      /  5.699819, 5.428224, 5.459493/

      data c4      / -2.1,   -2.1,   -2.1/

      data c4a     / -0.5,   -0.5,   -0.5/

      data cRB     /  50,     50,     50/

      data c5      /  6.1600, 5.8699, 5.2480/

      data c6      /  0.4893, 0.4755, 0.4517/

      data cHM     /  3, 3, 3/

      data c7      /  0.08281534, 0.05232253, 0.0000000/

      data c9      /  0.7900, 0.9334, 0.6196/

      data c9a     /  1.500502, 1.915732, 2.668990/

      data c11     /  0.0125256, 0.010086769, 0.002604898/

      data cgamma1 / -0.0075,-0.0090,-0.0030/

      data cgamma2 / -0.0065,-0.0040,-0.0025/

      data cgamma3 / 4, 4, 4/

      data phi1    / -0.4417,-0.5697,-0.7990/

      data phi2    / -0.1417,-0.2927,-0.0699/

      data phi3    / -0.007010,-0.006141,-0.008444/

      data phi4    /  0.1021510, 0.255253, 0.058595/

      data phi5    /  0.2289, 0.2386, 0.4629/

      data phi6    /  0.014996, 0.014964, 0.005749/

      data phi7    /  580.0, 573.9, 391.8/

      data phi8    /  0.0700,-0.0019,-0.0412/

      data tau1    /  0.3437, 0.3601, 0.3577/

      data tau2    /  0.2637, 0.3076, 0.3419/

      data sigma1  /  0.4458, 0.4816, 0.4581/

      data sigma2  /  0.3459, 0.3902, 0.4213/

      data sigma3  /  0.8000, 0.8000, 0.7504/
c  PERIOD       Spectral period (sec), enter -2 for peak velocity
c  M            Magnitude
c  F_RV         Enter 1 for reverse faulting, 0 otherwise default is SS)
c  F_NM         Enter 1 for normal faulting, 0 otherwise (default is SS)
c  Z_TOR        Depth to top of rupture (km)
c  DELTA        Fault dip (degrees)
c  V_S30        Average shear wave velocity (m/s)
c  F_Inferred   Enter 1 if Vs30 is inferred, 0 if measured
c  Z1           Depth to Vs = 1 km/sec (m), enter -999 for default
c  R_Rup        Distance to rupture plane (km)
c  R_JB         Joyner-Boore distance (km)
c  R_x          Normal distance to strike of fault (km), negative on footwall

        real cc, gamma, cosDELTA, psa, psa_ref, NL0, tau, sigma_NL0,
     1       total_app
        real pi, d2r, r1, r2, r3, r4, fw, hw, a, b, c, rkdepth
c      real	sigmaf,sqrt2
c      real	sigmaf
        integer iprd, sa


        pi = atan(1.0)*4.0
        d2r = pi/180.0

        cc = c5(iprd)* cosh(c6(iprd) * max(M-cHM(iprd),0.))
        gamma = cgamma1(iprd) +
     1          cgamma2(iprd)/cosh(max(M-cgamma3(iprd),0.))
        cosDELTA = cos(DELTA*d2r)
c Magnitude scaling
        r1 = c1(iprd) + c2(iprd) * (M-6.0) +
     1       (c2(iprd)-c3(iprd))/cn(iprd) *
     1             log(1.0 + exp(cn(iprd)*(cM(iprd)-M)))

c Near-field magnitude and distance scaling
        r2 = c4(iprd) * log(R_Rup + cc)
c Far-field distance scaling
        r3 = (c4a(iprd)-c4(iprd)) *
     1            log(sqrt(R_Rup*R_Rup+cRB(iprd)*cRB(iprd))) +
     1       R_Rup * gamma

c Scaling with other source variables (F_RV, F_NM, and Z_TOR)
        r4 = c1a(iprd)*F_RV +
     1       c1b(iprd)*F_NM +
     1       c7(iprd)*(Z_TOR - 4)
        if (M .lt. 6.0) r4 = r4 + c11(iprd)*min(DELTA-70,0.)

c HW effect
        if (R_x .lt. 0) then
          hw = 0.0
        else
          hw = c9(iprd) * tanh(R_x*cosDELTA**2/c9a(iprd)) *
     1        (1.0 - sqrt(R_JB**2+Z_TOR**2)/(R_Rup + 0.001))
        endif


        psa_ref = r1+r2+r3+r4+hw

c Soil effect: linear response
        a = phi1(iprd) * min(log(V_S30/1130.), 0.)
c Soil effect: nonlinear response
        b = phi2(iprd) *
     1(exp(phi3(iprd)*(min(V_S30,1130.)-360.))-exp(phi3(iprd)*(1130.-360.)))
        c = phi4(iprd)
c Corrections to ln(Vs30) scaling
c
c !!  NOTE: On the 3rd line, I cap max(0,Z1-15) at 300 to avoid
c !!          the numerical overflow of intrinsic function 'cosh'.
c
        rkdepth = phi5(iprd) *
     1        ( 1 - 1.0/cosh(phi6(iprd)*max(0.,Z1-phi7(iprd)))) +
     1        phi8(iprd)/cosh(0.15*min(max(0., Z1-15.),300.))


c....... Polulation mean of ln(psa) (eta=0)

c        psa = exp(psa_ref + (a + b * log((exp(psa_ref)+c)/c)) + rkdepth)
c dont exp the ground accel. keep logged. SH nov 14 2012.
      psa = psa_ref + (a + b * log((exp(psa_ref)+c)/c)) + rkdepth
        psa_ref = exp(psa_ref)
         gndout(1) = psa
         if(l_gnd_ep(ip))then
         gndout(2)= psa+gnd_ep(ide,ime,ip)
         gndout(3)= psa-gnd_ep(ide,ime,ip)
         endif
c....... Total variance of ln(psa) about the population mean:
c          The approximate method (Equation 21)

        NL0 = b * psa_ref/(psa_ref+c)

        tau = tau1(iprd) +
     1            (tau2(iprd)-tau1(iprd))/2*(min(max(M,5.),7.)-5.)
        sigma_NL0 = sigma1(iprd) +
     1              (sigma2(iprd)-sigma1(iprd))/2*(min(max(M,5.),7.)-5.)
        sigma_NL0 = sigma_NL0 *
     1        sqrt(0.7*F_Measured+F_Inferred*sigma3(iprd)+(1+NL0)**2)
      if(fix_sigma)then
      total_app=sigma_fx
      else
        total_app = sqrt((tau*(1+NL0))**2+sigma_NL0**2)
        endif
c            print *,gndout(1),total_app
        sigmaf=1./total_app/sqrt2
        return
        end subroutine CY2012_NGA

      subroutine gksa12(L,ip,Mw,x,gndout,sigmaf,Vs,mec,Bdepth,Q)
c Revised Graizer and Kalkan model with continuous response variation with 
c      basin depth. Also PGA has been identified with 0.01s SA due to basin
c      effect. 
c Input: L = period index in below per array.
c      x = Rcd (km).
c      Mw = moment magnitude.
c Bdepth = basin depth (km). new parm. 2012. Geotech: depth to Vs=1.5 km/s
c       Q quality factor 435 for California according to GK
c mec is sense-of-slip:
ccccc      Mec = 1 for strike slip    Normal Faults (3) are treated like strike slip
ccccc      Mec = 2 for thrust faults
ccccc      Mec = 4 combination of thrust and strike
c Output:
c gndout(j) = logged median motion (SA) (g).
c sigmaf = 1/sigma/sqrt2

      real gnd_ep(3,3,8),gndout(3),sqrt2,Mw
      parameter (nper=35,sqrt2=1.41421356)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
       common/fix_sigma/fix_sigma,sigma_fx
      logical l_gnd_ep(8),fix_sigma
      real sigma_fx

      real bv, Va, c1, c2, c3, c4, c5, c6, c7, c8, c9
      real R0, D0, R1, D1, Y1, A, Dsp, Mu
      real per(35), sigma(35)
      real e1, e2, e3, e4
      real a1, a2, a3
      real t1, t2, t3, t4
      real s1, s2, s3

      DATA PER/0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.12,0.14,
     &         0.16,0.18,0.20,0.22,0.24,0.27,0.30,0.33,
     &         0.36,0.4,0.46,0.5,0.6,0.75,0.85,1.0,1.5,
     &         2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/
      DATA sigma/0.542,0.543,0.544,0.547,0.549,0.549,0.553,0.558,0.568,
     &         0.577,0.584,0.591,0.597,0.602,0.610,0.616,0.622,
     &         0.627,0.634,0.642,0.647,0.658,0.667,0.672,0.681,0.703,
     &         0.718,0.74,0.756,0.768,0.775,0.78,0.78,0.78,0.78/
C***********************************************************
C*********Coefficients file, dec 2012 ***************
         c1 =  0.140
         c2 = -6.250
         c3 =  0.370
         c4 =  2.237
         c5 = -7.542
         c6 = -0.125
         c7 =  1.190
         c8 = -6.150
         c9 =  0.525
         bv = -0.240 
         Va = 484.5
         R1 = 100.0
         c11= 1.0
         D1 = 0.65
         sigpga = 0.550
         if(fix_sigma)sigpga=sigma_fx
         sigmaf=1./sigpga/sqrt2


c********* Distance Dependance ********************

c********* Amplification Factor depending upon style of faulting **********
ccccc      Mec = 1 for strike slip    Normal Faults are treated as strike
ccccc      Mec = 2 for thrust faults
ccccc      Mec = 4 combination of thrust and strike

      IF(Mec.EQ.1.or. Mec.eq.3) THEN
          Frv = 1.0
       ELSE IF(Mec.EQ.2) THEN
          Frv = 1.28
       ELSEIF(Mec.EQ.4) THEN
          Frv = 1.14
      ENDIF
      A = (c1*atan(Mw+c2)+c3)*Frv
      A = A/1.096
      Y1 = alog(A)

c********* Distance R0 Factor ********************
ccccc    Corner distance R0 depends upon magnitude

         R0 = c4 * Mw + c5

c********* Damping D0 Factor **************************
ccccc    Damping is magnitude dependent

         D0 = c6 * cos(c7*(Mw + c8)) + c9

c********* Basin Effect *******************************

c-------  Depth dependance ----------------------------

      Abdepth = 1.4 / sqrt((1.0-(1.5/(Bdepth+0.1))**2)**2

     &           + 1.96*(1.5/(Bdepth+0.1))**2)

c-------  Distance dependance -------------------------

      Abdist = 1. / sqrt((1.0-(40./(x+0.1))**2)**2

     &           + 1.96*(40./(x+0.1))**2)

c-------  Slope of SA at long-periods ----------------
      Slope = 1.763 - 0.25 * atan(1.4*(Bdepth-1.))
c*************** Calculating non-basin PGA *****************

c------- Vs30 site effect --------------------------
         Y1 = Y1 + bv * (alog(Vs/Va))
c-------- Core Filter and Anelastic ----------------

         x0 = x/R0
         x1 = sqrt(x/R1)
         Y = Y1 - 0.5 * alog((1-x0)**2 + 4*D0**2*x0)
     &       - (c11 * x/Q)
ccccccc End non-basin PGA Calculation cccccccccccccccccccccccc


c------SA calculations ---------------------------------------
c Treat PGA as 0.01 SA because of the basin effect below... SH. Dec 4 2012.
          e1=-0.0012
          e2=-0.4085
          e3= 0.0006
          e4= 3.63
          a1= 0.01686
          a2= 1.2695
          a3= 0.0001
          Dsp=0.75
          t1= 0.0022
          t2= 0.63
          t3=-0.0005
          t4=-2.1
          s1=0.001
          s2=0.077
          s3=0.3251
          Mu = e1*x+e2*Mw+e3*Vs+e4

          Amp = (a1*Mw+a2)*exp(a3*x)

          Si = s1*x-(s2*Mw+s3)

          Tspo = t1*x+t2*Mw+t3*Vs+t4
          Pern = (per(L)/Tspo)**Slope

          temp1 = (alog(per(L))+Mu)/Si
          temp2 = temp1**2
          SA1 = Amp*exp(-0.5*temp2)
          SA2 = 1./sqrt((1-Pern)**2 + 4.*Dsp**2
     &               *Pern)

          SA = SA1 + SA2
          SA1 = SA*exp(Y)
     &              * (1.+ABdist*ABdepth/1.3)
          SA = alog(SA1)
           gndout(1)=SA
         if(l_gnd_ep(ip))then
         gndout(2)= SA+gnd_ep(ide,ime,ip)
         gndout(3)= SA-gnd_ep(ide,ime,ip)
         endif
         if(fix_sigma)then
         sigmaf=1./sigma_fx/sqrt2
         else
           sigmaf=1./sigma(L)/sqrt2
           endif
           return
        END subroutine gksa12

      subroutine gksa(ipgk,ip,xmag,rcd,gndout,sigmaf,Vs,islip,IB)
c***  Graizer & Kalkan July 2007, Attenuation **************
c***  Graizer & Kalkan June 2009, SA from .01s to 10 s.**********
c *** G & K have a GMPE update according to Vladimir phone conversation Nov 2012.
c IB=1 if site over sed basin 0 otherwise
c ip = global period index, .ge. 1
c ipgk = local Graizer period index, 0 for PGA. Otherwise, see PER below.
c islip = 1 strike slip
c islip = 2 reverse
c islip = 3 normal
c islip = 4 mixed str-slip and reverse. (They may mean oblique slip. Or multiseg
c      with different modes on diff segs? Or whatever?)
c islip definition differs moderately from GK supplied software,
c  but corresponds with hazFX nomenclature
      real gnd_ep(3,3,8),gndout(3),sqrt2,Mw
      parameter (nper=35,sqrt2=1.41421356)
       common/fix_sigma/fix_sigma,sigma_fx
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      logical l_gnd_ep(8),fix_sigma
      Real X,x0,x1,Y, Z
      REAL Mu,Si,Amp,Tspo
      REAL SA,SA1,SA2
      Integer Mec, IB, ide,ime
c      INTEGER KEY(3)
      real bv, Va, c1, c2, c3, c4, c5, c6, c7, c8, c9
      real R0, D0, R1, D1, Y1, A, Dsp
      real per(35), sigma(35), Pern
      real e1, e2, e3, e4
      real a1, a2, a3
      real t1, t2, t3, t4
      real s1, s2, s3
      DATA PER/0.01,0.02,0.03,0.04,0.06,0.08,0.1,0.12,0.14,
     &         0.16,0.18,0.20,0.22,0.24,0.27,0.30,0.33,
     &         0.36,0.4,0.46,0.5,0.6,0.75,0.85,1.0,1.5,
     &         2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0/
      DATA sigma/0.542,0.543,0.544,0.547,0.549,0.549,0.553,0.558,0.568,
     &         0.577,0.584,0.591,0.597,0.602,0.610,0.616,0.622,
     &         0.627,0.634,0.642,0.647,0.658,0.667,0.672,0.681,0.703,
     &         0.718,0.74,0.756,0.768,0.775,0.78,0.78,0.78,0.78/
c      DATA KEY/0,0,1/
C***********************************************************
C     Mw MOMENT MAGNITUDE
C     Vs Vs30 in m/sec
C     Mec = 1 FOR STRIKE SLIP & NORMAL, 2 FOR DIP SLIP FAULT, 3 combination of dip and strike
C     IB = 0 No Basin Effect, 1 with Basin Effect
C     NPTS Number of data points in the input file
         Mw =xmag 
C*********Coefficients file, July_2006 ***************
         c1 =  0.140
         c2 = -6.250
         c3 =  0.370
         c4 =  2.237
         c5 = -7.542
         c6 = -0.125
         c7 =  1.190
         c8 = -6.150
         c9 =  0.525
         bv = -0.240
         Va = 484.5
         R1 = 100.0
         Frv= 1.28
         sigpga = 0.550
         sigmaf=1./sigpga/sqrt2
c********* Mechanism Dependance ********************
      IF(islip.EQ.1.or.islip.eq.3) THEN
          Frv = 1.0
c         A = c1*atan(Mw+c2)+c3
c        Write(12,'(/A)') '   Mechanism - Strike-slip or normal '
       ELSE IF(islip.EQ.2) THEN
          Frv = 1.28
c         A = (c1*atan(Mw+c2)+c3)*Frv
c        Write(12,'(/A)') '   Mechanism - Dip-slip   '
       ELSE
          Frv = 1.14
c         A = (c1*atan(Mw+c2)+c3)*1.14
      ENDIF
      A = (c1*atan(Mw+c2)+c3)*Frv
      A = A/1.12
      Y1 = alog(A)
c********* Distance Dependance ********************
c********* Distance R0 Factor ********************
ccccc    Corner distance R0 depends upon magnitude

         R0 = c4 * Mw + c5

c********* Damping D0 Factor **************************
ccccc    Damping is magnitude dependent

         D0 = c6 * cos(c7*(Mw + c8)) + c9

c********* Basin Effect *******************************
cccccc   If basin is deeper than 1 km IB=1, otherwise IB=0
      IF(IB.EQ.0) THEN
         D1 = 0.65
       ELSE
         D1 = 0.35
      ENDIF
c***************************************************
c******* Calculating linear Vs30 site effect **********
         Y1 = Y1 + bv * (alog(Vs/Va))

         x0 = rcd/R0
         x1 = sqrt(rcd/R1)

         Y = Y1 - 0.5 * alog((1-x0)**2 + 4*D0**2*x0) -

     &       0.5 * alog((1-x1)**2 + 4*D1**2*x1)
       pga=Y
ccc End PGA Calculation cccccccccccccccccccccccccccccccc
      if(ipgk.eq.0)then
           gndout(1)=pga
         if(l_gnd_ep(ip))then
         gndout(2)= pga+gnd_ep(ide,ime,ip)
         gndout(3)= pga-gnd_ep(ide,ime,ip)
         endif
         if(fix_sigma)then
         sigmaf=1./sigma_fx/sqrt2
         else
           sigmaf=1./sigpga/sqrt2
           endif
      return
      endif
c------SA --------------------------------------------
          e1=-0.0012
          e2=-0.4085
          e3= 0.0006
          e4= 3.63
          a1= 0.01686
          a2= 1.2695
          a3= 0.0001
          Dsp=0.75
          t1= 0.0022
          t2= 0.63
          t3=-0.0005
          t4=-2.1
          s1=0.001
          s2=0.077
          s3=0.3251

          Mu = e1*rcd+e2*Mw+e3*Vs+e4
c          Mu = e1*rcd+e2*Mw+e3

          Amp = (a1*Mw+a2)*exp(a3*rcd)

          Si = s1*rcd-(s2*Mw+s3)

          Tspo = t1*rcd+t2*Mw+t3*Vs+t4
c          Tspo = t1*rcd+t2*Mw+t3

          Pern = (per(ipgk)/Tspo)**1.5

          temp1 = (alog(per(ipgk))+Mu)/Si
          temp2 = temp1**2

          SA1 = Amp*exp(-0.5*temp2)

c         SA1 = Amp*exp(-0.5*((alog(per(ipgk))+Mu)/

c    &               Si))**2

          SA2 = 1./sqrt((1-Pern)**2 + 4.*Dsp**2

     &               *Pern)

          SA = SA1 + SA2
c Y is the logged PGA. 
          SA1 = SA*exp(Y)
           SA=alog(SA1)
c SA1, SA2 are the -sigma and +sigma ground motions. Not needed
c           SA1=SA/exp(sigma(ipgk))
c           SA2=SA*exp(sigma(ipgk))
           gndout(1)=SA
         if(l_gnd_ep(ip))then
         gndout(2)= SA+gnd_ep(ide,ime,ip)
         gndout(3)= SA-gnd_ep(ide,ime,ip)
         endif
         if(fix_sigma)then
         sigmaf=1./sigma_fx/sqrt2
         else
           sigmaf=1./sigma(ipgk)/sqrt2
           endif
      return
      end subroutine gksa
      ! --------------------------------------------------------- ngaw2_gm_sub4y
  ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
      function y_ngaw2_no_site_amps(m, rjb, mech,
     :                    e, mh, c, mref, rref, 
     :                    h,
     :                    r,
     :                    fm, fd, fs,
     :                    gspread)
     
! NOTE: y in g, unless pgv!!!!  

! ncoeff_s1, ncoeff_s2 are not included as arguments because it is
! assumed that they are 4 and 7, respectively.

! I should probably generalize thus function subprogram by using 
! funcs for stages 1 and 2.

! Assume no site amp

! mech = 1, 2, 3 for SS, NS, RS, respectively

! Dates: 12/14/11 - Modified to use calculations from ba08_gm_suby  
!        12/15/11 - Allow for possibility of M-dependent anelastic attenuation and 
!                   a quadratic M-dependence for M> Mhinge (the
!                   coefficients for these options were 0.0 in BA08)
!                 - Convert to natural logs
!        12/19/11 - There are now 7, not 8, e coefficients (because
!                   I no longer include a mechanism unspecified case).
!        01/06/12 - Added implicit none
!                 - Add fm, fd, fs to output
!                 - Remove ansig, sigt from argument list
!        01/10/12 - Change mech to be 1, 2, 3 for SS, NS, RS (instead of 0, 1, 2)
!                 - Add "r" to argument list.
!        01/11/12 - Add "gspread" to argument list.
!                 - Use "mech" rather than "mechindex"
!        11/10/12 - Renamed and revised version of y_ba08_no_site_amps
 
      IMPLICIT none
      
      real :: m, e(*), mh, c(*), mref, rref,
     :        h, alny, gspread
     
      integer :: mech
      
      real :: r, rjb, fm_mech, fm_mscale, fm, fd, fs
      
      real :: y_ngaw2_no_site_amps
      
       
      if (mech /= 1 .and. mech /= 2 .and. mech /= 3) then
        print *,' From y_ngaw2_no_site_amps, mech = ', 
     :             mech,
     :       ' not a valid value; QUITTING'
        stop
      end if
      
      if (e(mech) == 0.0) then   ! no coefficient
        print *,' From y_ngaw2_no_site_amps, '//
     :           'mech>= 0, but'//
     :       ' e(mech) == 0.0; QUITTING'
        stop
      else
        fm_mech = e(mech)
      end if
        

      if (m < mh ) then      
        fm_mscale = e(4)*(m-mh)
     :       + e(5)*(m-mh)**2 
      else
        fm_mscale = e(6)*(m-mh)
     :       + e(7)*(m-mh)**2 
      end if
 
      fm = fm_mech + fm_mscale
 
      r = sqrt(rjb**2+h**2)
      
      gspread = c(1) + c(2)*(m-mref)

      fd = gspread*alog(r/rref)
     :   + c(3)*(r-rref)
     :   + c(4)*(m-mref)*(r-rref)

      fs = 0.0  ! No amps
         
      alny = fm + fd + fs 
      
      y_ngaw2_no_site_amps = exp(alny)
      
      return
      end function y_ngaw2_no_site_amps
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      
      subroutine getIdriss2012(iper,ip,xmag,rrup,vs30,gndout,sigmaf)
c Dec 3 2012: for pga and many periods to 10s.
c ip is period index in calling program. 
c iper is period index in this subroutine. 
      parameter (pi=3.14159265,sqrt2=1.414213562,vref=760.,nper=22)
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      real gnd_ep(3,3,8),gndout(3)
      logical l_gnd_ep(8)
      common/mech/ss,rev,normal,obl
      real, dimension (nper):: a1, a2, a3, b1,b2, x,g,j, Period
c 22 periods in Idriss2012 NGA Dec 2012:
c----  This version assumes v30 in neighborhood of 760 m/sec (linear Vs30 dependency)
c----  uses ln coefficients
c Limits 2012:
c      5 <=M<=8.5
c      Vs30 > 450 m/s; for Vs30>1200, use 1200
c      Rrup <= 150 km (but we go to 200 km).
c
      real sig,xmag,T,vs30,vscap,sigma_fx
       common/fix_sigma/fix_sigma,sigma_fx
      logical fix_sigma,ss,rev,normal,obl
        Period  = (/0.01,0.02,0.03,0.04,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.,1.5,2.,3.,4.,5.,7.5,10./)
c coeffs. from dec 2012 powerpoint progress report. Two sets which to use depends on Mw (xmag)
      if(xmag.le.6.75)then
        a1 = (/4.0246,4.0496,4.1246,4.3982,3.6009,2.877,4.4729,5.0966,6.607,7.2428,7.9132,7.6416,
     & 7.6753,7.3511,6.2227,4.19,3.1218,0.1913,-2.2774,-4.3775,-7.3922,-8.8253/)
        a2 = (/0.2058,0.2058,0.2058,-0.0113,0.0625,0.1128,0.0848,0.1713,0.1041,0.0875,0.0003,0.0027,
     & 0.0399,0.0689,0.16,0.2429,0.3966,0.756,0.9283,1.1209,1.4016,1.5574/)
        a3 = (/0.0589,0.0589,0.0589,0.0265,0.0417,0.0527,0.0442,0.0329,0.0188,0.0095,-0.0039,
     & -0.0133,-0.0224,-0.0267,-0.0198,-0.0367,-0.0291,-0.0214,-0.024,-0.0202,-0.0219,-0.0035/)
        b1 = (/2.9827,2.9827,2.9827,2.9239,2.9629,3.0486,3.0899,2.8296,2.8792,2.8389,2.8378,
     & 2.8191,2.8196,2.7802,2.7707,2.7388,2.7343,2.7339,2.6448,2.6208,2.5649,2.5454/)
        b2 = (/-0.2287,-0.2287,-0.2287,-0.239,-0.2418,-0.2513,-0.2516,-0.2236,-0.2229,-0.22,
     &-0.2284,-0.2318,-0.2337,-0.2392,-0.2398,-0.2417,-0.245,-0.2389,-0.2514,-0.2541,-0.2593,-0.2586/)
        x = (/-0.36,-0.36,-0.36,-0.2,-0.15,-0.06,-0.23,-0.42,-0.55,-0.65,-0.7,-0.7,-0.77,-0.86,-0.83,
     & -0.69,-0.76,-0.73,-0.68,-0.62,-0.6,-0.65/)
      g = (/-0.0025,-0.0025,-0.0025,-0.004,-0.0038,-0.003,-0.0026,-0.0038,-0.0027,-0.0027,-0.0029,-0.0032,-0.0032,-0.0025,
     & -0.0026,-0.0022,-0.0021,-0.002,-0.0033,-0.0037,-0.0023,-0.002/)
        j = (/0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.06,0.04,0.02,0.02,0.,0.,0.,0./)
                else
        a1 = (/5.9497,5.9747,6.0497,6.3233,5.3208,4.7279,6.319,6.7018,8.4558,9.0253,9.8038,9.7171,9.9095,9.9059,9.0339,7.3956,
     & 6.6264,4.1307,1.9629,0.0772,-2.621,-3.9783/)
        a2 = (/-0.0794,-0.0794,-0.0794,-0.2965,-0.1923,-0.1614,-0.1887,-0.0665,-0.1698,-0.1766,-0.2798,-0.3048,-0.2911,
     & -0.3097,-0.2565,-0.232,-0.1226,0.1724,0.3001,0.4609,0.6948,0.8393/)
        a3 = (/0.0589,0.0589,0.0589,0.0265,0.0417,0.0527,0.0442,0.0329,0.0188,0.0095,-0.0039,-0.0133,-0.0224,-0.0267,-0.0198,
     & -0.0367,-0.0291,-0.0214,-0.024,-0.0202,-0.0219,-0.0035/)
        b1 = (/2.9827,2.9827,2.9827,2.9239,2.896,2.9223,2.884,2.4516,2.5119,2.3873,2.3727,2.3304,
     & 2.3113,2.23,2.1861,2.0996,2.0519,1.984,1.8429,1.777,1.6383,1.5727/)
        b2 = (/-0.2287,-0.2287,-0.2287,-0.239,-0.2319,-0.2326,-0.2211,-0.1676,-0.1685,-0.1531,-0.1595,-0.1594,-0.1584,-0.1577,
     & -0.1532,-0.147,-0.1439,-0.1278,-0.1326,-0.1291,-0.122,-0.1145/)
        x = (/-0.36,-0.36,-0.36,-0.2,-0.15,-0.06,-0.23,-0.42,-0.55,-0.65,-0.7,-0.7,-0.77,-0.86,-0.83,
     & -0.69,-0.76,-0.73,-0.68,-0.62,-0.6,-0.65/)
      g = (/-0.0025,-0.0025,-0.0025,-0.0040,-0.0038,-0.0030,-0.0026,-0.0038,-0.0027,-0.0027,-0.0029,-0.0032,-0.0032,
     & -0.0025,-0.0026,-0.0022,-0.0021,-0.0020,-0.0033,-0.0037,-0.0023,-0.0020/)
      j = (/0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.06,0.04,0.02,0.02,0.,0.,0.,0./)
                endif
c T is constrained spectral period see Idriss EqSpectra Mar 2008 for the details. use aleatory sigma from the article.
c This sigma may be revised.
      T = max(period(iper),0.05)
      T = min(T,3.0)
      if(fix_sigma)then
      sig=sigma_fx
      else
      sig = 1.28 + 0.05*alog(T) - 0.08 * xmag
      endif	!special study can fix all sigma at a specified value. jan 7 2012.
          sigmaf= 1./sig/sqrt2
          vscap=min(vs30,1200.)
c-- 
c 2nd line: extra mag scaling, anelastic attn scaling and site -term scaling
          gnd= a1(iper)+a2(iper)*xmag-(b1(iper)+b2(iper)*xmag)*alog(rrup+10.)
     & + a3(iper)*(8.5-xmag)**2+g(iper)*rrup+ x(iper)*alog(vscap)
c sense of slip sensitivity: oblique-reverse or reverse slip. Otherwise none.
          if(rev.or.obl)then
          gnd = gnd+ j(iper)
          endif
        gndout(1)=gnd
         if(l_gnd_ep(ip))then
         gndout(2)= gnd+gnd_ep(ide,ime,ip)
         gndout(3)= gnd-gnd_ep(ide,ime,ip)
         endif
        return
      end subroutine getIdriss2012

c------------------------------------------------------------------------------
      SUBROUTINE CB13_NGA_SPEC_IN
c read in the coeffs only. Harmsen Feb 25 2013. Additional coeffs such as phi_low
C     This program computes 5%-damped linear elastic response spectra from the 
C     2012 Campbell-Bozorgnia NGA-West2 Ground Motion Prediction Equation (GMPE)
C
C     It is Based on Campbell-Bozorgnia NGA-West2 Draft Model, November 12, 2012, Rev111
C
C     Written by Yousef Bozorgnia (November 25, 2012)

      parameter (nper=23)
c      parameter (nper=22)
c should this be 23 or 22 - declared with 22 in main routine
      common /cb13/Phi,Tau,c0,c1,c2low,c2,c3,c4,c5,c6,c7,c8,c9,c10,c10j,c10jlow,c11,c11j,
     + c12,c13low,c13hi,c14,c15,k1,k2,k3,a2,h1,h2,h3,h4,h5,h6,phi_lny,phi_low, phi_hi,
     + tau_lny,tau_lnyB,tau_low,tau_hi,phi_lnAF,rho
      common/cb13p/Per
      REAL, dimension(nper) ::  Phi, Tau, csoil,nsoil,tau_low,tau_hi
      REAL Per(nper), c0(nper), c1(nper), c2low(nper), c2(nper)
      REAL c3(nper), c4(nper)
      REAL c5(nper), c6(nper), c7(nper), c8(nper), c9(nper)
      REAL c10(nper), c10J(nper), c10Jlow(nper) 
      REAL c11(nper), c11J(nper) 
      REAL, dimension(nper) :: c12, c13low,c13hi, c14,c15,k1, k2, k3,a2,h1,h2,h3
      REAL, dimension(nper) :: h4,h5,h6, phi_low, phi_hi,phi_lny,phi_lnyB,tau_lny,tau_lnyB
      REAL phi_lnAF(nper), rho(nper)
      character*1 dum1
   
C.....
C.....READ MODEL COEFFICIENTS FROM AN EXCEL FILE (text FORMAT, csv does not work on sun workstation)
C.....
c      if(icase.ne.1) go to 21 
      call get_lun(m)
c updated input file name Feb 25 2013.
        open (m,file='GR/CB13_Coeffs_Atten113_RE_fix_all_Feb_15_2013.txt',status='old',err=2012)
      print *,'CB coefficient file opened OK unit ',m
      read (m,1) dum1
 1    format(a1)

      do 20 iper=1,nper
        read (m,*)Per(iper),c0(iper),c1(iper),c2low(iper),c2(iper),c3(iper),
     +c4(iper),
     +c5(iper),
     +c6(iper),
     +c7(iper),
     +c8(iper),
     +c9(iper),
     +c10(iper),
     +c10Jlow(iper),
     +c10J(iper),
     +c11(iper),
     +c11J(iper),
     +c12(iper),
     +c13low(iper),c13hi(iper),c14(iper),
     +c15(iper),a2(iper), h1(iper),h2(iper),h3(iper),h4(iper),h5(iper),h6(iper),
     +k1(iper), k2(iper), k3(iper),
     +csoil(iper), nsoil(iper),
     +phi_low(iper), phi_hi(iper),tau_low(iper), tau_hi(iper), phi_lnAF(iper), sigma_c, rho(iper)
c      print *,'period ',Per(iper)
c      print *,c14(iper),c15(iper),phi_lnAF(iper), rho(iper),'c14(iper),c15(iper),phi_lnAF(iper), rho(iper)'
 20   CONTINUE
c Although read in csoil and nsoil will be real variables in the actual CB_SPEC_IN
       close(m)
       return
2012      print *,'CB coeff file not found ; put in GR/ folder'
      stop
       end SUBROUTINE CB13_NGA_SPEC_IN

c add Abrahamson Silva code with no Ry_from_Rx (the simplified model).
c Example call from above:
       
      SUBROUTINE CB13_NGA_SPEC  (ip,iper, Mw,Rrup,Rx,Rjb, Frv,Fnm,Ztor,Hhyp,W,Z25,Vs30,Dip,SJ,gnd,sigmaf)
c     +(Mw,Rrup,Rjb,Rx, Frv,Fnm,Ztor,Hhyp,W,Dip,Vs30,Z25,SJ, 
c     &Per,Y,Phi,Tau,Sigmatot,icase)
c call this subroutine after calling CB13_NGA_SPEC_IN
c mods by  Harmsen Feb 25 2013. They added PGV index 23 this time.
c call this with PGA as first period to run for a given source. This is needed to get the A1100 term
c used with other periods. SH. Mod of Bozorgnia code which performs internal loop on period.
C.....Input parameters (from file CB12_NGA_infile.txt)
c      ip = global period index
c      iper = local period index
C     Frv    = 1 for Reverse-faulting, 0 otherwise
C     Fnm    = 1 for Normal-faulting, 0 otherwise
C     Ztor   = Depth to top of coseismic rupture (km)
C     Hhyp   = Hypocentral depth (km)
C     W      = Fault width (km)
C     Dip    = Dip of fault rupture plane (degrees)
C     Vs30   = Average velocity in top 30-m (m/s)
C     Z25    = Depth to 2.5 km/s velocity horizon (km)
C     SJ     = 1 for sites in Japan; =0 otherwise
C     iSpec  = 0 for generating Sa(g) [Not at this point: 1 for Sv(cm/s); 2 for Sd(cm)]
c output
      parameter (nper=23,sqrt2=1.4142136)
       real, dimension(3) ::gnd
       common/sdi/sdi,dy_sdi,fac_sde
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      integer ide,ime,ip
       common/fix_sigma/fix_sigma,sigma_fx
      common/Dipinf/dipang,cosDELTA,cdipsq,cyhwfac,cbhwfac
      common /cb13/Phi,Tau,c0,c1,c2low,c2,c3,c4,c5,c6,c7,c8,c9,c10,c10j,c10jlow,c11,c11j,
     + c12,c13low,c13hi,c14,c15,k1,k2,k3,a2,h1,h2,h3,h4,h5,h6,phi_lny,phi_low, phi_hi,
     + tau_lny,tau_lnyB,tau_low,tau_hi,phi_lnAF,rho
c       common/cb13p/Per
       real, dimension(8) :: fac_sde
      real, dimension(nper):: Phi,Tau,Per,c0,c1,c2low,c2,c3,c4,c5,c6,c7,c8,c9,c10,c10j,c10jlow,c11,
     + c11j,c12,c13low,c13hi,c14,c15,k1,k2,k3,a2,h1,h2,h3,h4,h5,h6,phi_lny,phi_low, phi_hi,
     + tau_lny,tau_lnyB,tau_low,tau_hi,phi_lnAF,rho
       real alpha,Mw,A1100,csoil,nsoil,Sigmatot,sigma_fx,phi_lnyB
       logical l_gnd_ep(8),fix_sigma,sdi
c If sdi is true, this subroutine will return inelastic spectral displ. (cm) instead
c of pSA. Tothong&Cornell approach
       real gnd_ep(3,3,8)
       save A1100
c-----Soil model constants (not using the constant vector from Bozorgnia code)
    
       nsoil = 1.18
      csoil = 1.88
        if(iper.eq.1)iper=22    

C*****Magnitude term
      IF (Mw .LE. 4.5) THEN
        F_mag = c0(iper) + c1(iper)*Mw
      ELSEIF (Mw .LE. 5.5) THEN
        F_mag = c0(iper) + c1(iper)*Mw + c2low(iper)*(Mw-4.5)
      ELSEIF (Mw .LE. 6.5) THEN
        F_mag = c0(iper) + c1(iper)*Mw + c2low(iper)*(Mw-4.5) + 
     &          c2(iper)*(Mw-5.5)
      ELSE
        F_mag = c0(iper) + c1(iper)*Mw + c2low(iper)*(Mw-4.5) + 
     &          c2(iper)*(Mw-5.5) + c3(iper)*(Mw-6.5)
      ENDIF
      
C*****Distance term
      R = SQRT(Rrup**2 + c6(iper)**2)
      F_dis = (c4(iper) + c5(iper)*Mw)*LOG(R)

C*****Style-of-fauting term
      f_flt_F = c7(iper)*Frv + c8(iper)*Fnm

C.....Note: Magnitude limits and equation have been changed
      IF (Mw.LE.4.5) THEN
        f_flt_M=0.
      ELSEIF (Mw.LE.5.5) THEN
        f_flt_M= Mw-4.5 
      ELSE
        f_flt_M= 1.
      ENDIF

      F_flt= f_flt_F * f_flt_M

C*****Hanging-wall term 
C     Jennifer Donahue's HW Model plus CB08 distance taper 
      R1= W * cos(Dip*PI/180.)
      R2= 62.*Mw - 350.

      f1_Rx= h1(iper) + h2(iper)*(Rx/R1) + h3(iper)*((Rx/R1)**2)
      f2_Rx= h4(iper) + h5(iper)*((Rx-R1)/(R2-R1)) + 
     &       h6(iper)*((Rx-R1)/(R2-R1))**2

C.... CB08 distance taper
      IF (Rrup.eq.0.0) THEN
        f_HW_Rrup= 1.0
      ELSE 
        f_HW_Rrup= (Rrup - Rjb)/Rrup
      ENDIF
C
C.....f_HW_R
      IF (Rx.lt.0.0) THEN
        f_HW_Rx= 0.0
      ELSEIF (Rx.lt.R1) THEN
        f_HW_Rx= f1_Rx 
      ELSE
        f_HW_Rx= (max(f2_Rx, 0.0))
      ENDIF

C.....f_HW_M
C.....Note: Equation for f_HW_M has been changed
      IF (Mw.le.5.5) THEN
        f_HW_M=0.
      ELSEIF (Mw.le.6.5) THEN 
        f_HW_M= (Mw-5.5)*(1+a2(iper)*(Mw-6.5))
      ELSE
        f_HW_M= 1. + a2(iper)*(Mw-6.5)
      ENDIF

C.....f_HW_Z
      IF(Ztor.le.16.66) THEN
        f_HW_Z=1.0-0.06* Ztor
      ELSE
        f_HW_Z=0.0
      ENDIF

C.....f_HW_dip
      f_HW_dip= (90.0 - DIP)/45.0

      F_HW= c9(iper)* f_HW_Rx * f_HW_Rrup * f_HW_M * f_HW_Z * f_HW_Dip

C*****Hypo depth term
      IF(Hhyp.le.7.) THEN
        f_Hhyp_H=0.
      ELSEIF (Hhyp.le.20.) THEN
        f_Hhyp_H= Hhyp-7.
      ELSE
        f_Hhyp_H=13.0
      ENDIF

      IF (Mw.le.5.5) THEN
        F_Hhyp_M= c13low(iper)
      ELSEIF (Mw.le.6.5) THEN        
        F_Hhyp_M= (c13low(iper)+ (c13hi(iper) - c13low(iper))*(Mw-5.5)) 
      ELSE
        F_Hhyp_M= c13hi (iper)
      ENDIF
      
      F_Hhyp= f_Hhyp_H * f_Hhyp_M  
      
C*****Dip term
C.....Dip term has been changed
      IF (Mw.le.4.5) THEN
        F_Dip= c14(iper)* DIP
      ELSEIF (Mw.le.5.5) THEN 
        F_Dip= c14(iper)* (5.5 - Mw)* Dip
      ELSE
        F_Dip= 0.
      ENDIF

C*****Shallow sediment depth and 3-D basin term
      IF (Z25 .LE. 1.0) THEN
        F_sed = (c11(iper) + c11J(iper)*SJ)*(Z25 - 1.0)
      ELSEIF (Z25 .LE. 3.0) THEN
        F_sed = 0.0
      ELSE
        F_sed = c12(iper)*k3(iper)*(EXP(-0.75))*
     +(1.0 - EXP(-0.25*(Z25 - 3.0)))
      ENDIF

C*****Anelastic attenuatin term
      F_atn=c15(iper)*Rrup

C*****For the first period (loop), computer A1100 *****************************
      IF (iperflag.eq.1)THEN
C........Shallow site conditions term for ROCK PGA (i.e., Vs30 = 1100 m/s)
C        Note csoil and nsoil are now arrays
         F_site_1100 = (c10(22) + k2(22)*nsoil)*LOG(1100.0/k1(22))

C........Rock PGA

         A1100 = EXP(F_mag + F_dis + F_flt + F_HW + 
     +               F_site_1100 + F_sed + F_Hhyp + F_Dip + F_atn)
         iperflag=0
         go to 1000
C     Note: Statment number 1000 is now at different place 
      ENDIF
c******************************************************************************

C*****Site term for other iper values: Note csoil and nsoil are now arrays 
       IF (Vs30 .LE. k1(iper)) THEN
         F_site = c10(iper)*LOG(Vs30/k1(iper))
     &             + k2(iper)*(LOG(A1100+csoil*
     &                ((Vs30/k1(iper))**nsoil)) 
     &             - LOG(A1100 + csoil))
       ELSE
         F_site = (c10(iper) + k2(iper)*nsoil)*LOG(Vs30/k1(iper))
       ENDIF

c*****Ground motion parameter, logged.

c       gnd(1) = EXP(F_mag + F_dis + F_flt + F_HW + 
c     &              F_site + F_sed + F_Hhyp + F_Dip + F_atn)
       gnd(1) = F_mag + F_dis + F_flt + F_HW + 
     &              F_site + F_sed + F_Hhyp + F_Dip + F_atn
C.....

C.....CALCULATE ALEATORY UNCERTAINTY
C     Note: This part has been changed

1000  IF (Vs30 .LT. k1(iper)) THEN
        alpha = 
     &    k2(iper)*A1100*(1.0/(A1100  
     &    +csoil*(Vs30/k1(iper))**nsoil) 
     &    -1.0/(A1100 + csoil))
      ELSE
        alpha = 0.0
      ENDIF

      If (Mw.le.4.5) then
         tau_lnyB(iper) =tau_low(iper)
      elseif (Mw.lt.5.5) then
         tau_lnyB(iper) =tau_hi(iper) + 
     &    (tau_low(iper) - tau_hi(iper))*(5.5-Mw)
      else
         tau_lnyB(iper) =tau_hi(iper) 
      endif

      tau_lnPGAB = tau_lnyB(22) 

      tau(iper) = SQRT(tau_lnyB(iper)**2 +  
     &           (alpha * tau_lnPGAB)**2 +
     &           2.0*alpha*rho(iper)*tau_lnyB(iper)*tau_lnPGAB)

      If (Mw.le.4.5) then
         phi_lny(iper) =phi_low(iper)
      elseif (Mw.lt.5.5) then
         phi_lny(iper) =phi_hi(iper) + 
     &    (phi_low(iper) - phi_hi(iper))*(5.5-Mw)
      else
         phi_lny(iper) =phi_hi(iper) 
      endif

      phi_lnyB = SQRT(phi_lny(iper)**2 - phi_lnAF(iper)**2)

      phi_lnPGAB = SQRT(phi_lny(22)**2 - phi_lnAF(22)**2)

      phi(iper) = SQRT(phi_lny(iper)**2 + 
     &           (alpha*phi_lnPGAB)**2 +
     &           2.0*alpha*rho(iper)*phi_lnyB*phi_lnPGAB)
      
      if(fix_sigma)then
      	Sigmatot=sigma_fx
      else      	
        Sigmatot = SQRT(phi(iper)**2 + Tau(iper)**2)
       endif
       if(sdi)then
       sde=gnd(1)+fac_sde(ip)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)      !10 is an upper bound for rhat.
       gnd(1) = sdi_ratio(per(iper),Mw,rhat,Sigmatot,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
         else
      sigmaf=1./sqrt2/Sigmatot
         endif      !if(sdi)
         if(l_gnd_ep(ip))then
         gnd(2)= gnd(1)+gnd_ep(ide,ime,ip)
         gnd(3)= gnd(1)-gnd_ep(ide,ime,ip)
         endif
      RETURN
      END SUBROUTINE CB13_NGA_SPEC

c------------------------------------------------------------------------------
      subroutine lin_interp(x, y, n, j, x_intrp, y_intrp)
      
* Computes linearly interpolated value of y

* Values out of range are assigned end values

* Dates: 03/16/05 - Written by D. Boore
*        07/24/05 - Added index j for end cases

      real x(*), y(*), slope, x_intrp, y_intrp
      integer j, n
      
      if (x_intrp .le. x(1)) then
        j = 1
        y_intrp = y(1)
        return
      end if
      
      if (x_intrp .ge. x(n)) then
        j = n
        y_intrp = y(n)
        return
      end if          
      
      call locate(x,n,x_intrp,j)
      
      slope = (y(j+1) - y(j))/(x(j+1)-x(j))
      y_intrp = y(j) + slope*(x_intrp - x(j))
      
      return
      end
      subroutine locate(x,n,y,j)
      real x(*)
      jp=2
      dowhile(y.gt.x(jp).and.jp.lt.n)
      jp=jp+1
      enddo
      j=jp-1
      return
      end subroutine locate
      


! --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

* Finds a logical unit number not in use; returns
* -1 if it cannot find one.

* Dates -- 05/19/98 - Written by D. Boore, following
*                     Larry Baker's suggestion

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
! --------------------------- END GET_LUN ----------------
     

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
3      continue
c      print *,'b1,b2,b3,b4,b5,b6,a,b:'
c      print *,b1,b2,b3,b4,b5,b6,a,b
       if(rhat.lt.0.2)then
        sdi_ratio=0.0
        sdi = sige
        elseif(rhat.lt.0.3)then
        sdi_ratio=g1(rhat,M,b1,b2,b3,b4,b5) - g1(0.2,M,b1,b2,b3,b4,b5)
c because in this interval g2 is zero we don't see it ablove
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
c-----------------------------------------------------------------------
c Below version of CY is March 5 2013...Includes some source directivity response
c
        subroutine CY2013_NGA(ip,iprd, M, R_Rup, R_JB, R_x, V_S30,
     1                   F_Measured, F_Inferred, deltaZ1, DELTA, Z_TOR,
     1                   F_RV, F_NM, cDPP, gndout, tau, sigma_NL0, sigmaf)

        implicit none
        integer nprd
        real pi,d2r,sqrt2
        parameter (nprd=24,pi=3.14159265,d2r=17.45329252e-3,sqrt2=1.414213562)
c Predictor variables
c cDPP = source directivity parameter. Need guidance on what to use Mar 6 2013.
c      Without such guidance, set direct=false, meaning do not use the cDPP model.
c ip = global period index
c iprd = local period index corresponding to below coeff. vectors.
      common/epistemic/l_gnd_ep,gnd_ep,ide,ime
      integer ide,ime,ip
       common/fix_sigma/fix_sigma,sigma_fx
       common/sdi/sdi,dy_sdi,fac_sde
      real gnd_ep(3,3,8),gndout(3), sigmaf,sigma_fx
      logical l_gnd_ep(8),fix_sigma,direct/.false./,sdi
        real PERIOD, M, Width, R_rup, R_JB, R_x, V_S30, F_Measured,
     1       F_Inferred, Z1, DELTA, Z_TOR, F_RV, F_NM, deltaZ1, cDPP
c     1       F_Inferred, Z1, DELTA, Z_TOR, F_RV, F_NM, deltaZ1, cDPP,sigmaf
      real fac_sde(8),dy_sdi,rhat,sdisd,sde,sdi_ratio
c Model cofficients
        real prd(nprd),
     1       c1(nprd), c1a(nprd), c1b(nprd), c2(nprd),
     1       c3(nprd), cn(nprd),  cM(nprd),  c4(nprd),
     1       c4a(nprd),cRB(nprd), c5(nprd),  c6(nprd),
     1       cHM(nprd),c7(nprd),  c8(nprd),  c8a(nprd), c8b(nprd),
     1       c9(nprd),  c9a(nprd), c11(nprd),
     1       cgamma1(nprd), cgamma2(nprd), cgamma3(nprd),
     1       phi1(nprd), phi2(nprd), phi3(nprd), phi4(nprd),
     1       phi5(nprd), phi6(nprd), phi7(nprd), phi8(nprd),
     1       tau1(nprd), tau2(nprd),
     1       sigma1(nprd), sigma2(nprd), sigma3(nprd)

      data prd     / -1.,     -2.,      0.0100, 0.0200, 0.0300,
     1                0.0400, 0.0500, 0.0750, 0.1000, 0.1500,
     1                0.2000, 0.2500, 0.3000, 0.4000, 0.5000,
     1                0.7500, 1.0000, 1.5000, 2.0000, 3.0000,
     1                4.0000, 5.0000, 7.5000,10.0000/

      data c1  / -1.336521 , 0.000000 ,-1.336521 ,-1.2782066,-1.1511862,
     1           -0.9720268,-0.8242828,-0.536952 ,-0.4395333,-0.5711291,
     1           -0.6313851,-0.780782 ,-0.8835749,-1.2348672,-1.5644566,
     1           -2.1300655,-2.5226237,-3.0518322,-3.2884456,-3.637457,
     1           -3.786648 ,-3.8581326,-3.9839868,-4.0284115/

      data c1a     /  0.1000, 0.1094, 0.1000, 0.1000, 0.1000,
     1                0.1000, 0.1000, 0.1000, 0.1000, 0.1000,
     1                0.1000, 0.1000, 0.0999, 0.0997, 0.0991,
     1                0.0936, 0.0766, 0.0022,-0.0591,-0.0931,
     1               -0.0982,-0.0994,-0.0999,-0.1000/

      data c1b     / -0.2550,-0.0626,-0.2550,-0.2550,-0.2550,
     1               -0.2550,-0.2550,-0.2540,-0.2530,-0.2500,
     1               -0.2449,-0.2382,-0.2313,-0.2146,-0.1972,
     1               -0.1620,-0.1400,-0.1184,-0.1100,-0.1040,
     1               -0.1020,-0.1010,-0.1010,-0.1000/

      data c2 /      1.06,1.06,1.06,1.06,1.06,
     1               1.06,1.06,1.06,1.06,1.06,
     1               1.06,1.06,1.06,1.06,1.06,
     1               1.06,1.06,1.06,1.06,1.06,
     1               1.06,1.06,1.06,1.06/

      data c3 /      1.965223,0.000000,1.965223,1.965223,1.965223,
     1               1.965223,1.965223,1.965223,1.965223,2.003098,
     1               2.125631,2.314981,2.415959,2.588119,2.688933,
     1               2.811177,2.864086,2.914551,2.928747,2.945157,
     1               2.952093,2.955774,2.960026,2.961823/


      data cn /      2.996,1.648,2.996,3.292,3.514,
     1               3.563,3.547,3.448,3.312,3.044,
     1               2.831,2.658,2.505,2.261,2.087,
     1               1.812,1.648,1.511,1.470,1.456,
     1               1.465,1.478,1.498,1.502/

      data cM /      5.593398,0.000000,5.593398,5.614333,5.621921,
     1               5.624729,5.623966,5.608298,5.573932,5.469153,
     1               5.418777,5.370959,5.349052,5.312861,5.301406,
     1               5.350243,5.440512,5.667317,5.792988,6.036039,
     1               6.228029,6.392756,6.718305,6.966668/



      data c4 /      -2.1,-2.1,-2.1,-2.1,-2.1,
     1               -2.1,-2.1,-2.1,-2.1,-2.1,
     1               -2.1,-2.1,-2.1,-2.1,-2.1,
     1               -2.1,-2.1,-2.1,-2.1,-2.1,
     1               -2.1,-2.1,-2.1,-2.1/

      data c4a/      -0.5,-0.5,-0.5,-0.5,-0.5,
     1               -0.5,-0.5,-0.5,-0.5,-0.5,
     1               -0.5,-0.5,-0.5,-0.5,-0.5,
     1               -0.5,-0.5,-0.5,-0.5,-0.5,
     1               -0.5,-0.5,-0.5,-0.5/

      data cRB/      50,50,50,50,50,
     1               50,50,50,50,50,
     1               50,50,50,50,50,
     1               50,50,50,50,50,
     1               50,50,50,50/

      data c5/       6.1600,5.1700,6.1600,6.1580,6.1550,
     1               6.1508,6.1441,6.1200,6.0850,5.9871,
     1               5.8699,5.7547,5.6527,5.4997,5.4029,
     1               5.2900,5.2480,5.2194,5.2099,5.2040,
     1               5.2020,5.2010,5.2000,5.2000/

      data c6/       0.4893,0.4407,0.4893,0.4892,0.4890,
     1               0.4888,0.4884,0.4872,0.4854,0.4808,
     1               0.4755,0.4706,0.4665,0.4607,0.4571,
     1               0.4531,0.4517,0.4507,0.4504,0.4501,
     1               0.4501,0.4500,0.4500,0.4500/

      data cHM/      3,3,3,3,3,
     1               3,3,3,3,3,
     1               3,3,3,3,3,
     1               3,3,3,3,3,
     1               3,3,3,3/

      data c7/       0.0848,0.0000,0.0848,0.0868,0.0920,
     1               0.0984,0.1015,0.0971,0.0875,0.0673,
     1               0.0565,0.0440,0.0376,0.0267,0.0206,
     1               0.0148,0.0124,0.0078,0.0056,0.0025,
     1               0.0011,0.0005,0.0000,0.0000/

      data c8/       0.2154,0.0000,0.2154,0.2154,0.2154,
     1               0.2154,0.2154,0.2154,0.2154,0.2154,
     1               0.2154,0.2154,0.2154,0.2154,0.2154,
     1               0.2154,0.2154,0.2154,0.2154,0.2154,
     1               0.2154,0.2154,0.2154,0.2154/

      data c8a/      0.2695,0.0000,0.2695,0.2695,0.2695,
     1               0.2695,0.2695,0.2695,0.2695,0.2695,
     1               0.2695,0.2695,0.2695,0.2695,0.2695,
     1               0.2695,0.2695,0.2695,0.2695,0.2695,
     1               0.2695,0.2695,0.2695,0.2695/

      data c8b/      4.610, 0.000, 4.610, 4.610, 4.610,
     1               4.610, 4.610, 4.610, 4.610, 4.610,
     1               4.610, 4.610, 4.610, 4.610, 4.610,
     1               5.038, 5.341, 5.837, 6.072, 6.500,
     1               6.803, 7.039, 7.467, 7.770/

      data c9/       0.7900,0.0000,0.7900,0.8129,0.8439,
     1               0.8740,0.8996,0.9442,0.9677,0.9660,
     1               0.9334,0.8946,0.8590,0.8019,0.7578,
     1               0.6788,0.6196,0.5101,0.3917,0.1244,
     1               0.0086,0.0000,0.0000,0.0000/

      data c9a/      15,15,15,15,15,
     1               15,15,15,15,15,
     1               15,15,15,15,15,
     1               15,15,15,15,15,
     1               15,15,15,15/

      data c11/      0.01289,0.00000,0.01289,0.01307,0.01365,
     1               0.01428,0.01463,0.01472,0.01361,0.01052,
     1               0.00916,0.00776,0.00708,0.00588,0.00507,
     1               0.00371,0.00278,0.00144,0.00093,0.00022,
     1               0.00000,0.00000,0.00000,0.00000/


      data cgamma1/  -0.007500, 0.000000,-0.007500,-0.007715,-0.007831,
     1               -0.007931,-0.008000,-0.009170,-0.010000,-0.009322,
     1               -0.009000,-0.008029,-0.007500,-0.006092,-0.005000,
     1               -0.003830,-0.003000,-0.002661,-0.002500,-0.002500,
     1               -0.002500,-0.002500,-0.002500,-0.002500/

      data cgamma2/  -0.006500, 0.000000,-0.006500,-0.006715,-0.006831,
     1               -0.006931,-0.007000,-0.005830,-0.005000,-0.004322,
     1               -0.004000,-0.003676,-0.003500,-0.002937,-0.002500,
     1               -0.002500,-0.002500,-0.002500,-0.002500,-0.002500,
     1               -0.002500,-0.002500,-0.002500,-0.002500/


      data cgamma3/  4,4,4,4,4,
     1               4,4,4,4,4,
     1               4,4,4,4,4,
     1               4,4,4,4,4,
     1               4,4,4,4/

      data phi1/     -0.5771, 0.0000,-0.5771,-0.5543,-0.5040,
     1               -0.4510,-0.4266,-0.4332,-0.4782,-0.6064,
     1               -0.6893,-0.7870,-0.8373,-0.9283,-0.9854,
     1               -1.0547,-1.0864,-1.1064,-1.1138,-1.1095,
     1               -1.0628,-1.0098,-0.8703,-0.7451/

      data phi2/     -0.1417,-0.0699,-0.1417,-0.1364,-0.1403,
     1               -0.1591,-0.1862,-0.2538,-0.2943,-0.3113,
     1               -0.2927,-0.2662,-0.2405,-0.1975,-0.1633,
     1               -0.1028,-0.0699,-0.0425,-0.0302,-0.0129,
     1               -0.0016, 0.0000, 0.0000, 0.0000/

      data phi3/     -0.007010,-0.008444,-0.007010,-0.007279,-0.007354,
     1               -0.006977,-0.006467,-0.005734,-0.005604,-0.005845,
     1               -0.006141,-0.006439,-0.006704,-0.007125,-0.007435,
     1               -0.008120,-0.008444,-0.007707,-0.004792,-0.001828,
     1               -0.001523,-0.001440,-0.001369,-0.001361/

      data phi4/     0.102151,5.410000,0.102151,0.108360,0.119888,
     1               0.133641,0.148927,0.190596,0.230662,0.266468,
     1               0.255253,0.231541,0.207277,0.165464,0.133828,
     1               0.085153,0.058595,0.031787,0.019716,0.009643,
     1               0.005379,0.003223,0.001134,0.000515/

      data phi5/       0.000,   0.000,   0.000,   0.000,   0.000,
     1                 0.000,   0.000,   0.000,   0.000,   0.000,
     1                 0.000,   0.000,   0.001,   0.004,   0.010,
     1                 0.034,   0.067,   0.155,   0.203,   0.277,
     1                 0.309,   0.321,   0.329,   0.330/


      data phi6/      300.00,  000.00,  300.00,  300.00,  300.00,
     1                300.00,  300.00,  300.00,  300.00,  300.00,
     1                300.00,  300.00,  300.00,  300.00,  300.00,
     1                300.00,  300.00,  300.00,  300.00,  300.00,
     1                300.00,  300.00,  300.00,  300.00/


      data phi7/       00.00,   00.00,   00.00,   00.00,   00.00,
     1                 00.00,   00.00,   00.00,   00.00,   00.00,
     1                 00.00,   00.00,   00.00,   00.00,   00.00,
     1                 00.00,   00.00,   00.00,   00.00,   00.00,
     1                 00.00,   00.00,   00.00,   00.00/

      data phi8/     4000.00, 4000.00, 4000.00, 4000.00, 4000.00,
     1               4000.00, 4000.00, 4000.00, 4000.00, 4000.00,
     1               4000.00, 4000.00, 4000.00, 4000.00, 4000.00,
     1               4000.00, 4000.00, 4000.00, 4000.00, 4000.00,
     1               4000.00, 4000.00, 4000.00, 4000.00/

      data tau1/     0.3437,0.2539,0.3437,0.3471,0.3603,
     1               0.3718,0.3848,0.3878,0.3835,0.3719,
     1               0.3601,0.3522,0.3438,0.3351,0.3353,
     1               0.3429,0.3577,0.3769,0.4023,0.4406,
     1               0.4784,0.5074,0.5328,0.5542/

      data tau2/     0.2637,0.2381,0.2637,0.2671,0.2803,
     1               0.2918,0.3048,0.3129,0.3152,0.3128,
     1               0.3076,0.3047,0.3005,0.2984,0.3036,
     1               0.3205,0.3419,0.3703,0.4023,0.4406,
     1               0.4784,0.5074,0.5328,0.5542/

      data sigma1/   0.4458,0.4496,0.4458,0.4458,0.4535,
     1               0.4589,0.4630,0.4702,0.4747,0.4798,
     1               0.4816,0.4815,0.4801,0.4758,0.4710,
     1               0.4621,0.4581,0.4493,0.4459,0.4433,
     1               0.4424,0.4420,0.4416,0.4414/

      data sigma2/   0.3459,0.3554,0.3459,0.3459,0.3537,
     1               0.3592,0.3635,0.3713,0.3769,0.3847,
     1               0.3902,0.3946,0.3981,0.4036,0.4079,
     1               0.4157,0.4213,0.4213,0.4213,0.4213,
     1               0.4213,0.4213,0.4213,0.4213/

      data sigma3/   0.8000,0.7504,0.8000,0.8000,0.8000,
     1               0.8000,0.8000,0.8000,0.8000,0.8000,
     1               0.8000,0.7999,0.7997,0.7988,0.7966,
     1               0.7792,0.7504,0.7136,0.7035,0.7006,
     1               0.7001,0.7000,0.7000,0.7000/

        real cc, gamma, cosDELTA, psa, psa_ref, NL0, tau, sigma_NL0,
     1       total_app
c        real pi, d2r, r1, r2, r3, r4, fw, hw, fd, a, b, c, rkdepth
        real r1, r2, r3, r4, fw, hw, fd, a, b, c, rkdepth
        integer iprd, sa

c        pi = atan(1.0)*4.0
c        d2r = pi/180.0
      fd=0.0	!forward directvity not computable without more knowledge
        cc = c5(iprd)* cosh(c6(iprd) * max(M-cHM(iprd),0.))
        gamma = cgamma1(iprd) +
     1          cgamma2(iprd)/cosh(max(M-cgamma3(iprd),0.))
        cosDELTA = cos(DELTA*d2r)

c Magnitude scaling
        r1 = c1(iprd) + c2(iprd) * (M-6.0) +
     1       (c2(iprd)-c3(iprd))/cn(iprd) *
     1             log(1.0 + exp(cn(iprd)*(cM(iprd)-M)))

c Near-field magnitude and distance scaling
        r2 = c4(iprd) * log(R_Rup + cc)

c Far-field distance scaling
        r3 = (c4a(iprd)-c4(iprd)) *
     1            log(sqrt(R_Rup*R_Rup+cRB(iprd)*cRB(iprd))) +
     1       R_Rup * gamma

c Scaling with other source variables (F_RV, F_NM, and Z_TOR)
        r4 = c1a(iprd)*F_RV +
     1       c1b(iprd)*F_NM +
     1       c7(iprd)*(Z_TOR - 4)
        if (M .lt. 6.0) r4 = r4 + c11(iprd)*min(DELTA-70.,0.)

c HW effect
        if (R_x .lt. 0) then
         hw = 0.0
        else
         hw = c9(iprd) * (cosDELTA**2) * (0.2+0.8*tanh(R_x/c9a(iprd))) *
     1        (1.0 - sqrt(R_JB**2+Z_TOR**2)/(R_Rup + 0.001))
        endif

c Directivity effect
        if(direct)fd = c8(iprd) * exp(-c8a(iprd) * (M-c8b(iprd))**2) *
     1       max(0., 1.-max(0.,R_Rup-40.)/30.) *
     1       min(max(0.,M-5.5)/0.8, 1.) * cDPP

        psa_ref = r1+r2+r3+r4+hw+fd
c        write (*,'(7f10.4)') r1, r2, r3, r4, hw, fd, psa_ref

c Soil effect: linear response
        a = phi1(iprd) * min(log(V_S30/1130), 0.0)

c Soil effect: nonlinear response
        b = phi2(iprd) *
     1(exp(phi3(iprd)*(min(V_S30,1130.)-360.))-exp(phi3(iprd)*(1130-360)))

        c = phi4(iprd)

c Corrections to ln(Vs30) scaling
c
c
        rkdepth = phi5(iprd) * ( 1.0 - exp(-deltaZ1/phi6(iprd) ) ) +
     1            phi7(iprd) * ( 1.0 - exp(-deltaZ1/phi8(iprd) ) )


c....... Polulation mean of ln(psa) (eta=0)
c keep the response logged.
c        psa = exp(psa_ref + (a + b * log((exp(psa_ref)+c)/c)) + rkdepth)
        psa = psa_ref + (a + b * log((exp(psa_ref)+c)/c)) + rkdepth
        psa_ref = exp(psa_ref)      ! using in NL0 computation below.
      gndout(1)=psa
c....... Total variance of ln(psa) about the population mean:
c          The approximate method (Equation 21)

        NL0 = b * psa_ref/(psa_ref+c)

        tau = tau1(iprd) +
     1            (tau2(iprd)-tau1(iprd))/2*(min(max(M,5.),7.)-5.)

        sigma_NL0 = sigma1(iprd) +
     1              (sigma2(iprd)-sigma1(iprd))/2*(min(max(M,5.),7.)-5.)

        sigma_NL0 = sigma_NL0 *
     1        sqrt(0.7*F_Measured+F_Inferred*sigma3(iprd)+(1.0+NL0)**2)

        total_app = sqrt((tau*(1.0+NL0))**2+sigma_NL0**2)
       if(sdi)then
       sde=gndout(1)+fac_sde(ip)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)      !10 is an upper bound for rhat.
       gndout(1) = sdi_ratio(prd(iprd),M,rhat,total_app,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
         else
      sigmaf=1./sqrt2/total_app
         endif      !if(sdi)
         if(l_gnd_ep(ip))then
         gndout(2)= gndout(1)+gnd_ep(ide,ime,ip)
         gndout(3)= gndout(1)-gnd_ep(ide,ime,ip)
         endif
        return

        end subroutine CY2013_NGA

      subroutine bssa2013drv(ip,indx_pga,indx_per,m,rjb,v30,z1,mech,iregion,lny,sigmaf)
      real sqrt2
      	parameter (npermx=23,sqrt2=1.414213562)
      real  c(3)
c this driver code runs the NGAW BSSA model for PGA and then for period with index indx_per
c the outputs are lny = logged SA, and sigmaf =1/sigma/sqrt2 
c ip is the global period index, used in SDI calcs for example.    
       common/sdi/sdi,dy_sdi,fac_sde
      real e(0:6)
      
      real, dimension(npermx):: T,clin, vclin, vref,
     :     c1,c2,c3,phi2, phi3,e0,e1,e2,e3,e4,e5,e6,
     :     amh,f1, f3, f4,h,Dc3CATW,Dc3China,Dc3Italy,
     :     lambda1, lambda2,Mref,Rref,delta_c3,
     :     tau1, tau2
     
      logical sdi
      real m, logy_pred, rjb, z1,dy_sdi
      real lny, gnd, sdisd,sig, rhat
        real, dimension(8) :: fac_sde
     
      integer mech, iregion,ip
 

      real r, r_pga4nl
      
      real :: sigt, 
     :  fs, fsb,
     :  fp, fpb,
     :  fe,
     :  amp_total, 
     :  amp_nl, 
     :  amp_lin, 
     :  fp_pga4nl, fpb_pga4nl,
     :  fe_pga4nl,
     :  gspread_pga4nl,
     :  pga4nl, 
     :  per1, 
     :  per2,
     :  per_max, 
     :  sigt_dummy, 
     :  v30, 
     :  slope_logy, 
     :  slope_phi, 
     :  slope_tau, 
     :  slope_sigma, 
     :  phi, tau, sigmaf,
     :  yg, 
     :  ycgs,
     :  pi, twopi,
     :  gspread,
     :  weight

      real :: 
     :  r_t1,
     :  y_t1,
     :  fsb_t1, 
     :  fs_t1, 
     :  fpb_t1,
     :  fp_t1,
     :  fe_t1,
     :  amp_total_t1, 
     :  amp_nl_t1, 
     :  amp_lin_t1, 
     :  phi_t1, tau_t1, sigma_t1,
     :  gspread_t1
     
      real :: 
     :  r_t2,
     :  y_t2,
     :  fsb_t2, 
     :  fs_t2, 
     :  fpb_t2,
     :  fp_t2,
     :  fe_t2,
     :  amp_total_t2, 
     :  amp_nl_t2, 
     :  amp_lin_t2, 
     :  phi_t2, tau_t2, sigma_t2,
     :  gspread_t2
     
      real ::
     :  slope_fe, 
     :  slope_fpb,
     :  slope_fp,
     :  slope_fs,
     :  slope_fsb,
     :  slope_logamp_lin,
     :  slope_logamp_nl,
     :  slope_logamp_total,
     :  slope_logr,
     :  slope_gspread
      T       =(/-1.000000, 0.000000, 0.010000, 0.020000, 0.030000, 0.040000, 0.050000,
     + 0.075000, 0.100000, 0.150000, 0.200000, 0.250000, 0.300000, 0.400000, 0.500000,
     + 0.750000, 1.000000, 1.500000, 2.000000, 3.000000, 4.000000, 5.000000,10.000000/)
      c1      =(/-1.235000,-1.128000,-1.128000,-1.133000,-1.133000,-1.129000,-1.099000,
     +-1.071000,-1.045000,-1.037000,-1.048000,-1.059000,-1.081000,-1.111000,-1.150000,
     +-1.176000,-1.197000,-1.205000,-1.225000,-1.196000,-1.210000,-1.222000,-1.328000/)
      c2      =(/ 0.141900, 0.183200, 0.183000, 0.181300, 0.179800, 0.180400, 0.179000,
     + 0.174800, 0.163100, 0.141900, 0.135200, 0.133100, 0.128800, 0.116500, 0.117200,
     + 0.109500, 0.099800, 0.089500, 0.093000, 0.090400, 0.099900, 0.107200, 0.149500/)
      c3      =(/-0.003440,-0.008088,-0.008088,-0.008074,-0.008336,-0.009030,-0.009819,
     +-0.010580,-0.010200,-0.008977,-0.007717,-0.006517,-0.005475,-0.004053,-0.003220,
     +-0.001931,-0.001210,-0.000365, 0.000000, 0.000000,-0.000052, 0.000000, 0.000000/)
      Mref    =(/ 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000,
     + 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000,
     + 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000, 4.500000/)
      Rref    =(/ 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
     + 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.00000,
     + 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000/)
      h       =(/ 5.500000, 4.900000, 4.900000, 4.900000, 4.700000, 4.700000, 4.300000,
     + 4.300000, 4.300000, 4.700000, 4.900000, 4.900000, 5.100000, 5.300000, 5.700000,
     + 5.700000, 5.900000, 6.300000, 6.900000, 6.700000, 7.500000, 7.700000, 9.500000/)
      Dc3CATW =(/ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     + 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     + 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000/)
      Dc3China=(/ 0.004345, 0.002858, 0.002840, 0.002750, 0.002710, 0.002830, 0.002980,
     + 0.002940, 0.002900, 0.002760, 0.002730, 0.002440, 0.002150, 0.002060, 0.002300,
     + 0.002660, 0.002970, 0.003010, 0.003070, 0.002430, 0.002680, 0.002640, 0.003030/)
      Dc3Italy=(/-0.000330,-0.002550,-0.002550,-0.002210,-0.002160,-0.002030,-0.001950,
     +-0.002120,-0.002440,-0.002700,-0.002880,-0.003220,-0.003340,-0.003230,-0.002960,
     +-0.002540,-0.002210,-0.001430,-0.001050,-0.001280,-0.001080,-0.000840, 0.001490/)
      e0      =(/ 5.059000, 0.466700, 0.472500, 0.497100, 0.560200, 0.673500, 0.722100,
     + 0.967300, 1.081000, 1.339000, 1.366000, 1.256000, 1.229000, 1.129000, 1.009000,
     + 0.681800, 0.422800,-0.139200,-0.545000,-1.229000,-1.656000,-1.982000,-3.022000/)
      e1      =(/ 5.100000, 0.509400, 0.515200, 0.539500, 0.605800, 0.722000, 0.772700,
     + 1.018000, 1.120000, 1.380000, 1.408000, 1.283000, 1.247000, 1.145000, 1.027000,
     + 0.709600, 0.446200,-0.112700,-0.508500,-1.188000,-1.587000,-1.904000,-2.901000/)
      e2      =(/ 4.872000, 0.256000, 0.261600, 0.296700, 0.381600, 0.522200, 0.562500,
     + 0.768000, 0.832600, 1.073000, 1.165000, 1.071000, 1.027000, 0.929500, 0.808400,
     + 0.495700, 0.239200,-0.308200,-0.687100,-1.231000,-1.740000,-2.020000,-3.355000/)
      e3      =(/ 5.053000, 0.468300, 0.474100, 0.495400, 0.543600, 0.640300, 0.688100,
     + 0.948500, 1.104000, 1.365000, 1.365000, 1.279000, 1.274000, 1.178000, 1.054000,
     + 0.702400, 0.451000,-0.122800,-0.558500,-1.307000,-1.754000,-2.117000,-3.124000/)
      e4      =(/ 1.125000, 1.501000, 1.491000, 1.505000, 1.505000, 1.481000, 1.477000,
     + 1.489000, 1.582000, 1.347000, 1.201000, 1.038000, 0.950800, 1.004000, 1.041000,
     + 1.300000, 1.495000, 1.749000, 1.883000, 2.181000, 2.155000, 2.238000, 1.774000/)
      e5      =(/-0.142300, 0.081320, 0.080300, 0.084200, 0.097690, 0.104200, 0.107500,
     + 0.111900, 0.119400,-0.029000,-0.094740,-0.157900,-0.195400,-0.222300,-0.238700,
     +-0.215400,-0.196400,-0.160600,-0.134900,-0.035170,-0.037120,-0.002483,-0.210900/)
      e6      =(/ 0.245600,-0.137900,-0.137300,-0.138400,-0.138000,-0.145400,-0.153400,
     +-0.171300,-0.152300,-0.145700,-0.137100,-0.102400,-0.096080, 0.002565, 0.050700,
     + 0.110500, 0.187500, 0.395500, 0.458400, 0.655600, 0.788600, 0.852700, 1.091000/)
        amh      =(/ 6.200000, 5.500000, 5.500000, 5.500000, 5.500000, 5.500000, 5.500000
     +, 5.500000, 5.500000, 5.740000, 5.920000, 6.050000, 6.160000, 6.200000, 6.200000
     +, 6.200000, 6.200000, 6.200000, 6.200000, 6.200000, 6.200000, 6.200000, 6.200000/)
        clin       =(/-0.805000,-0.515000,-0.515000,-0.547500,-0.550000,-0.524440,-0.4900,
     +-0.450000,-0.470000,-0.510000,-0.600000,-0.700000,-0.780000,-0.850000,-0.900,
     +-0.980000,-1.030000,-1.090000,-1.100000,-1.050000,-1.000000,-0.970000,-0.644000/)
      Vclin =(/   1130.0,    836.7,    836.7,    928.3,    939.6,    923.2,    908.8,
     +    931.4,   1000.0,   1100.0,   1250.0,   1300.0,   1350.0,   1400.0,   1370.0,
     +   1250.0,   1170.0,   1100.0,   1000.0,    900.0,    800.0,    750.0,    750.0/)
      Vref    =(/    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,
     +    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,
     +    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,    760.0,    760.0/)
      phi2    =(/-0.100000,-0.150000,-0.150000,-0.145000,-0.150000,-0.170000,-0.190000,
     +-0.230000,-0.250000,-0.260000,-0.250000,-0.235000,-0.220000,-0.195000,-0.175000,
     +-0.140000,-0.110000,-0.065000,-0.035000,-0.012900,-0.001600, 0.000000, 0.000000/)
      phi3    =(/-0.008440,-0.007010,-0.007010,-0.007280,-0.007350,-0.006980,-0.006470,
     +-0.005730,-0.005600,-0.005850,-0.006140,-0.006440,-0.006700,-0.007130,-0.007440,
     +-0.008120,-0.008440,-0.007710,-0.004790,-0.001830,-0.001520,-0.001440,-0.001360/)
      f1      =(/ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     + 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     + 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000/)
      f3      =(/ 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000,
     + 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000,
     + 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000, 0.100000/)
      f4      =(/ 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     + 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
     + 0.000000, 0.412000, 0.412000, 0.471000, 0.588000, 0.647000, 0.659000, 0.676000/)
      lambda1 =(/ 0.503500, 0.511000, 0.512000, 0.511600, 0.516800, 0.523000, 0.528800,
     + 0.546500, 0.567500, 0.580000, 0.570300, 0.571200, 0.569700, 0.550800, 0.570000,
     + 0.563000, 0.572400, 0.580200, 0.605200, 0.651800, 0.700100, 0.744300, 0.646500/)
      lambda2 =(/ 0.589700, 0.661000, 0.662000, 0.664200, 0.681300, 0.700700, 0.718800,
     + 0.726500, 0.712900, 0.703500, 0.698200, 0.681800, 0.667400, 0.649800, 0.638800,
     + 0.628700, 0.627400, 0.621100, 0.620000, 0.624200, 0.604500, 0.593300, 0.548000/)
      tau1    =(/ 0.386900, 0.432500, 0.433300, 0.437700, 0.453700, 0.477400, 0.501000,
     + 0.525100, 0.504700, 0.450500, 0.418800, 0.405100, 0.404400, 0.411600, 0.421700,
     + 0.465200, 0.512800, 0.549300, 0.563800, 0.552100, 0.519700, 0.513600, 0.514200/)
      tau2    =(/ 0.371400, 0.435000, 0.295000, 0.298800, 0.309000, 0.328400, 0.349500,
     + 0.405200, 0.410900, 0.400900, 0.362300, 0.324600, 0.309700, 0.291900, 0.287500,
     + 0.357800, 0.404900, 0.443400, 0.472200, 0.446200, 0.418700, 0.421800, 0.451900/)

c      pi = 4.0*atan(1.0)
c      twopi = 2.0 * pi
       c(1)=c1(indx_pga)
       c(2)=c2(indx_pga)
       c(3)=c3(indx_pga)
       e(0)=e0(indx_pga)
       e(1)=e1(indx_pga)
       e(2)=e2(indx_pga)
       e(3)=e3(indx_pga)
       e(4)=e4(indx_pga)
       e(5)=e5(indx_pga)
       e(6)=e6(indx_pga)
c       delta_c3 = dc3CATW      !Use CA-Taiwan version for now. Can add options for other regions. iregion=1
            r_pga4nl = sqrt(rjb**2+h(indx_pga)**2)
        call y_bssa13_no_site_amps( m, r_pga4nl, mech,  c, Mref(indx_pga), 
     :           Rref(indx_pga), dc3CATW(indx_pga),
     :           e, amh(indx_pga), pga4nl) 
c       if(abs(m-7.).lt.0.05.and.ip.eq.2)print *,rjb,pga4nl,r_pga4nl/Rref(indx_pga),m,amh(indx_pga),' pga4nl'    
c   pga4nl is logged in this code.
            r = sqrt(rjb**2+h(indx_per)**2)
        c(1)=c1(indx_per)
        c(2)=c2(indx_per)
        c(3)=c3(indx_per)
        e(0)=e0(indx_per)
        e(1)=e1(indx_per)
        e(2)=e2(indx_per)
        e(3)=e3(indx_per)
        e(4)=e4(indx_per)
        e(5)=e5(indx_per)
        e(6)=e6(indx_per)
     
            call bssa13_gm_sub4y(m, r, v30, mech,  z1, pga4nl,c, Mref(indx_per), 
     :        Rref(indx_per), dc3CATW(indx_per),
     :        e, amh(indx_per),
     :        clin(indx_per), vclin(indx_per), vref(indx_per),
     :        phi2(indx_per), phi3(indx_per),   
     :        f1(indx_per), f3(indx_per), f4(indx_per),
     :        lambda1(indx_per), lambda2(indx_per), 
     :        tau1(indx_per), tau2(indx_per),
     :        lny, sigmaf)  
c       if(abs(m-7.).lt.0.05.and.ip.eq.1)print *,rjb,exp(lny),r,' pgaw'   
       if(sdi)then
       sig=1./sigmaf/sqrt2
       sde=lny+fac_sde(ip)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)      !10 is an upper bound for rhat.
       lny = sdi_ratio(T(indx_per),m,rhat,sig,sdisd) + sde
       sigmaf=1./sdisd/sqrt2      !use sdi for all gnd_ep branches
       endif
        return
        end subroutine bssa2013drv

       subroutine y_bssa13_no_site_amps(m, r, mech, c, mref, rref, delta_c3,e, mh, y)     
     
! NOTE: y in g, unless pgv!!!!  

! ncoeff_s1, ncoeff_s2 are not included as arguments because it is
! assumed that they are 3 and 7, respectively.

! Assume no site amp

! mech = 0, 1, 2, 3 for unspecified, SS, NS, RS, respectively

! Dates: 02/27/13 - Modified from y_ngaw2_no_site_amps
 
      IMPLICIT none
       real, dimension(3):: c     
      real :: m, r, mref, rref, delta_c3, 
     :        e(0:6), mh, 
     :        alny, y, fpb, fp, fe, fs, fault_type_term, gspread     
     
      integer :: mech, iregion
      
        fault_type_term = e(mech)
         

      if (m < mh ) then      
        fe = e(4)*(m-mh)
     :       + e(5)*(m-mh)**2 
      else
        fe = e(6)*(m-mh)
      end if
 
      fe = fault_type_term + fe
 
      gspread = c(1) + c(2)*(m-mref)

      fpb = gspread*alog(r/rref)
     :   + c(3)*(r-rref)

      fp = fpb + delta_c3*(r-rref)
      
      fs = 0.0  ! No amps
         
      y = exp(fe + fp + fs) 
      return
      end subroutine y_bssa13_no_site_amps
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>      


! --------------------------------------------------------- bssa13_gm_sub4y
      subroutine bssa13_gm_sub4y( m, r, v30, mech,  z1,  pga4nl,c, amref, 
     :        rref, delta_c3,
     :        e, amh,
     :        clin, vclin, vref,
     :        phi2, phi3,   
     :        f1, f3, f4,
     :        lambda1, lambda2, 
     :        tau1, tau2,
     :        alny, sigmaf)     
c Note output log(SA) not SA in this version. No need for SA in downstream calcs. 
c removed iregion: taken care of before arriving here.    
c z1 is depth to 1 km/s Vs. Units km (?). Meters doesn't work. 
!Input arguments:

!          per, m, rjb, 
!          mech, v30, 
!          e, amh_gmpe, c_gmpe, amref_gmpe,
!          rref_gmpe, h_gmpe, 
!          v30ref, 
!          sigt_gmpe,
!          e_pga4nl, amh_pga4nl, c_pga4nl, amref_pga4nl,
!          rref_pga4nl, h_pga4nl,


!Output arguments:

!          alny, sigmaf,
!          


! Computes NGAW2 NGA motions for given input variables

! Note: For generality, include a magnitude-dependent anelastic
! term, although the coefficients are 0.0 for the 2007 BA NGA GMPEs

 
! Dates: 02/27/13 - Written by D.M. Boore, based on ngaw2_gm_sub4y

      implicit none
      real sqrt2
      parameter (sqrt2=1.414213562)
      real, dimension(3) :: c
      real m, r, v30, z1, pga4nl, amref, rref, delta_c3,e(0:6), amh,
     :        clin, vclin, vref,phi2, phi3,  f1, f2, f3, f4,lambda1, lambda2, 
     :        tau1, tau2,
     :        alny, fsb, fs, 
     :        phi, tau, sigmaf,    
     :        amp_lin, amp_nl, amp_total, bnl,
     :        fz1
     
      real y_xamp
     
      integer :: mech, iregion

        

!GET Y FOR GMPE, WITHOUT SITE AMPS:

      call y_bssa13_no_site_amps(
     :           m, r, mech, c, amref, 
     :           rref, delta_c3,
     :           e, amh,
     :           y_xamp)     
      
!Compute site amps

      if (v30 <= vclin) then
        amp_lin  =  (v30/vref)**clin
      else
        amp_lin  =  (vclin/vref)**clin
      end if
        
      f2 = phi2 * (exp(phi3*(amin1(v30,760.0)-360.0)) -
     :             exp(phi3*(760.0-360.0)           )   )
      
      bnl = f2
      
      amp_nl = exp(f1+f2*alog( (pga4nl+f3)/f3 ))

      amp_total = amp_lin * amp_nl
      
      fsb  = alog(amp_total)   !
      
      if (z1 <= 0.3) then
        fz1 = 0.0
      else if (z1>= 2.0) then
        fz1 = 1.7*f4
      else
        fz1 = f4*(z1-0.3)
      end if
      
      fs = fsb + fz1
      
!Compute phi, tau, and sigma   

      if (m <= 5.0) then
        tau = tau1
      else if (m >= 6.0) then
        tau = tau2
      else 
        tau = tau1+(tau2-tau1)*(m-5.0)
      end if
c alny=the fundamental output (logged Mar 11 2013. SH).     
      alny = fs + alog(y_xamp)
       
      if (v30 <= 250.0) then
        phi = lambda1
      else if (v30 >= 400.0) then
        phi = lambda2
      else
        phi = lambda1 +(lambda2-lambda1)*alog(v30/250.0)/0.47
      end if
      
      sigmaf = 1./sqrt2/sqrt(phi**2 + tau**2)         
       return
      end subroutine bssa13_gm_sub4y
! --------------------------------------------------------- bssa13_gm_sub4y
c This code is checked by Sanaz Rezaeian based on AS 2008 Earthquake Spectra paper
c and modified based on AS 2009 errata (3/11/2013)
c   (look for SR for modifications)
c
c ------------------------------------------------------------------            
C *** Abrahamson and Silva (07/2007 - Model) Horizontal ************
c ------------------------------------------------------------------            
      subroutine AS_072007 ( ip,iper,mag, dip, Fn,Frv, width, rRup, rjb,R_x,
     1                     vs30, hwflag, lnSa, sigma1, depthtop,  vs30_class,
     3                     depthvs10 )
c ip = global period index, used in SDI calcs, for example.
c specT = spectral period (s).   iper = index in T array of specT  
c no interpolation. SH. Added option to compute inelastic spectral displacement
c using the method of Tothong and Cornell March 2013.
c Output: lnSA = logged median gm (g)
c      	sigma1 =  SD(lnSA)
       common/sdi/sdi,dy_sdi,fac_sde
       real, dimension (108):: period
       real, dimension(8):: fac_sde
       logical sdi
      real mag, dip, fType, aspectratio, rRup, rjb, vs30, pgaRock,
     1       srcSiteA, lnSa, sigma, tau, period1, sigma1, depthtop, width,
     2       depthvs10, lnSaTD, lnSa1, lnSa2, lnSaRock,dy_sdi
      integer hwflag,  vs30_class,ip
      data period / 0.0, -1.0, -2.0, 0.01, 0.02, 0.022, 0.025, 0.029,
     1              0.03, 0.032, 0.035, 0.036, 0.04, 0.042, 0.044, 
     2              0.045, 0.046, 0.048, 0.05, 0.055, 0.06, 0.065, 
     3              0.067, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 
     4              0.11, 0.12, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 
     5              0.18, 0.19, 0.2, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 
     6              0.3, 0.32, 0.34, 0.35, 0.36, 0.38, 0.4, 0.42, 0.44,
     7              0.45, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 0.667, 0.7,
     8              0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 
     9              1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 
     1              3.0, 3.2, 3.4, 3.5, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 
     1              5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 /

c     VS30 class is to distinguish between the sigma if the VS30 is measured
c     vs the VS30 being estimated from surface geology.
c         Vs30_class = 1 for measured 
c         Vs30_class = 0 for estimated
      specT = period(iper)
c     For now, convert ftype to an equivalent rake      
C     fType     Mechanism                      Rake
C     ------------------------------------------------------
C      -1       Normal                   -120 < Rake < -60.0
C     1, 0.5    Reverse and Rev/Obl        30 < Rake < 150.0
C     0,-0.5    Strike-Slip and NMl/Obl        Otherwise
C 
c      if ( fType .eq. 1.0 ) then
c        Frv = 1.0
c        Fn = 0.0
c      elseif ( fType .eq. 0.5 ) then
c        Frv = 1.0
c        Fn = 0.0
c      elseif ( fType .eq. -1.0 ) then
c        Frv = 0.0
c        Fn = 1.0
c      else
c        Frv = 0.0
c        Fn = 0.0
c      endif
      specT = per
C HW factors......
c     For now, convert hwflag to an equivalent srcSiteA      
      if (  hwflag .eq. 1 ) then
        srcSiteA = 90.
      else
        srcSiteA = 0.
      endif
     
c     compute pga on rock
      period0 = 0.0
      pgaRock = 0.0
      vs30_rock = 1100.
      depthvs10rock = 0.006
      call AS_072007_model (iper, mag, dip, rake, width, rRup, rjb,R_x,
     1                     vs30_rock, pgaRock,  lnSa, sigma, tau,
     2                     period0, depthTop, Frv, Fn,  vs30_class,
     3                     depthvs10rock, hwflag )

      pgaRock = exp(lnSa)

c     Compute the spectral period for the constant displacement model
c     (Eq. 22)
      TD = 10**(-1.25+0.3*mag)
      if (specT .ge. TD) then
c         Compute Sa at TD spectral period for Vs30=1,100m/sec
          call AS_072007_model (iper, mag, dip, rake, width, rRup, rjb,R_x,
     1                     vs30_rock, pgaRock, lnSaTD, sigmaX, tauX, 
     2                     TD, depthtop, Frv, Fn,  vs30_class, 
     3                     depthvs10rock, hwflag )

          lnSaRock = lnSaTD + alog((TD*TD)/(specT*specT))

c         Compute Sa at spectral period for Vs30
          call AS_072007_model (iper, mag, dip, rake, width, rRup, rjb,R_x,
     1                     vs30, pgaRock,  lnSa1, sigma, tau, 
     2                     specT, depthtop, Frv, Fn,  vs30_class, 
     3                     depthvs10, hwflag )
c         Compute Sa at spectral period for Vs30=1,100m/sec
          call AS_072007_model (iper, mag, dip, rake, width, rRup, rjb,R_x,
     1                     vs30_rock, pgaRock,  lnSa2, sigmaX, tauX, 
     2                     specT, depthtop, Frv, Fn,  vs30_class, 
     3                     depthvs10rock, hwflag )
C         Compute soil amplification (i.e., 'lnSa1-lnSa2') and add to Constant displcaement rock spectrum
          lnSa = lnSaRock + (lnSa1-lnSa2)
 
      else
C         For cases where specT < TD compute regular ground motions. 
          call AS_072007_model (iper, mag, dip, rake, width, rRup, rjb,R_x,
     1                     vs30, pgaRock, lnSa, sigma, tau, 
     2                     specT, depthtop, Frv, Fn,  vs30_class, 
     3                     depthvs10, hwflag )
      endif

c     compute Sa (given the PGA rock value)
      sigma1 = sqrt( sigma**2 + tau**2 )
       if(sdi)then
       sde=lnSa+fac_sde(ip)      !fac_sde is log(T**2/(4pisq))
       rhat = min(10.,exp(sde)/dy_sdi)      !10 is an upper bound for rhat.
       lnSa = sdi_ratio(specT,mag,rhat,sigma1,sdisd) + sde
       sigma1=sdisd      !use sdi for all gnd_ep branches
       endif
c     Don't Convert units spectral acceleration to gal                                
c      lnSa = lnSa + 6.89                                                
      return
      end subroutine AS_072007 
c ----------------------------------------------------------------------

      subroutine AS_072007_model ( iper,mag, dip, rake, width, rRup, rjb,R_x,
     1                     vs30, pgaRock,  lnSa, sigma, tau, 
     2                     specT, depthtop, Frv, Fn,  vs30_class,
     3                     z10, hwflag )
      
      parameter (MAXPER=108)
      real a1(MAXPER), a2(MAXPER), a3(MAXPER), a4(MAXPER), a5(MAXPER),
     1     a6(MAXPER), a7(MAXPER), a8(MAXPER), a9(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER),
     1     c4(MAXPER), a16(MAXPER), a17(MAXPER), a18(MAXPER),
     2     s1e(MAXPER), s2e(MAXPER), s1m(MAXPER), s2m(MAXPER), s3(MAXPER), s4(MAXPER)
      real period(MAXPER), b_soil(MAXPER), vLin(MAXPER), M1
      real sigma, tau,lnSa, pgaRock, vs30, rjb, rRup, R_x, aspectratio, rake,
     1     dip, mag, sigcorr(MAXPER), taucorr(MAXPER) 
c      real sig1(MAXPER), sig2(MAXPER), tau1(MAXPER), tau2(MAXPER)
c      real tauCorr, sigCorr
      real taper1, taper2, taper3, taper4, taper5
      real a1T, a2T, a3T, a5T, a6T, a7T, a8T, a9T
      real a10T, a11T, a12T, a13T, a14T, a15T, a18T
      real c4T, M1T, a4T
      real s1eT, s2eT, s1mT, s2mT, s3T, s4T, s1T, s2T
      real sigcorrT, taucorrT, vLinT, b_soilT, sum, Frv, Fn
      real damp_dpga, sigamp, width, testv1
      real sigmanot, sigmanotPGA, taunot, taunotPGA, sigmaB, sigmaBPGA
      integer iper,count1, count2,  vs30_class, hwflag
      real n, c, z10, c2, e2, a21, a22, test, zhat

c      	Updated coefficients according to AS08 paper by SR from matlab code
      data period / 0.0, -1.0, -2.0, 0.01, 0.02, 0.022, 0.025, 0.029,
     1              0.03, 0.032, 0.035, 0.036, 0.04, 0.042, 0.044, 
     2              0.045, 0.046, 0.048, 0.05, 0.055, 0.06, 0.065, 
     3              0.067, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 
     4              0.11, 0.12, 0.13, 0.133, 0.14, 0.15, 0.16, 0.17, 
     5              0.18, 0.19, 0.2, 0.22, 0.24, 0.25, 0.26, 0.28, 0.29, 
     6              0.3, 0.32, 0.34, 0.35, 0.36, 0.38, 0.4, 0.42, 0.44,
     7              0.45, 0.46, 0.48, 0.5, 0.55, 0.6, 0.65, 0.667, 0.7,
     8              0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 
     9              1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.5, 2.6, 2.8, 
     1              3.0, 3.2, 3.4, 3.5, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 
     1              5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0 /
      data Vlin / 865.1, 400.0, 400.0, 865.1, 865.1, 865.1, 865.1, 898.6, 
     1            907.8, 926.4, 953.6, 962.3, 994.5, 1008.9, 1022.0, 
     2           1028.0, 1033.8, 1044.3, 1053.5, 1071.5, 1082.8, 1088.3, 
     3           1089.1, 1089, 1085.7, 1079.2, 1070.1, 1059.0, 1046.3, 
     4           1032.5, 1002.5, 970.9, 939.1, 929.7, 907.8, 877.6, 
     5            848.7, 821.2, 795.3, 771.0, 748.2, 706.8, 670.6, 
     6            654.3, 639.0, 611.4, 598.9, 587.1, 565.9, 547.1, 
     7            538.6, 530.6, 516.0, 503.0, 491.5, 481.2, 476.5, 
     8            472.1, 463.9, 456.6, 441.6, 430.3, 421.7, 419.3, 
     9            415.2, 410.5, 407.0, 404.6, 402.9, 401.8, 400.0, 400.0,
     1            400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 
     1            400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0,
     2            400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 
     3            400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0,
     4            400.0, 400.0, 400.0, 400.0 /
      data b_soil / -1.186, -1.955, 0.0, -1.186, -1.219, -1.232, -1.25,
     1              -1.269, -1.273, -1.281, -1.291, -1.295, -1.308,
     2              -1.315, -1.323, -1.326, -1.33, -1.338, -1.346, 
     3              -1.367, -1.391, -1.416, -1.426, -1.443, -1.471,
     4              -1.5, -1.53, -1.561, -1.592, -1.624, -1.687, -1.751,
     5              -1.813, -1.831, -1.873, -1.931, -1.988, -2.041,
     6              -2.093, -2.141, -2.188, -2.272, -2.347, -2.381,
     7              -2.412, -2.469, -2.494, -2.518, -2.558, -2.592,
     8              -2.607, -2.62, -2.641, -2.657, -2.668, -2.674,
     9              -2.676, -2.676, -2.675, -2.669, -2.643, -2.599,
     1              -2.543, -2.521, -2.476, -2.401, -2.319, -2.233,
     1              -2.142, -2.049, -1.955, -1.762, -1.57, -1.382, 
     2              -1.199, -1.025, -0.859, -0.703, -0.558, -0.423,
     3              -0.299, -0.086, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     4               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     5               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data c4 / 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     1          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     2          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     3          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     4          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     5          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 
     6          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     7          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     8          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5,
     9          4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5 /
      data a1 / 0.725, 5.772, 5.065, 0.731, 0.797, 0.826, 0.869, 0.917, 
     1          0.929, 0.954, 0.99, 0.998, 1.033, 1.054, 1.075, 1.086,
     2          1.096, 1.118, 1.138, 1.194, 1.25, 1.302, 1.322, 1.349, 
     3          1.393, 1.433, 1.469, 1.503, 1.538, 1.566, 1.608, 1.643,
     4          1.66, 1.663, 1.657, 1.667, 1.667, 1.667, 1.66, 1.653,
     5          1.64, 1.611, 1.574, 1.56, 1.544, 1.508, 1.487, 1.466,
     6          1.457, 1.433, 1.422, 1.411, 1.392, 1.372, 1.353, 1.335,
     7          1.326, 1.317, 1.298, 1.28, 1.239, 1.195, 1.149, 1.133, 
     8          1.1, 1.054, 1.002, 0.945, 0.893, 0.844, 0.816, 0.757,
     9          0.7, 0.636, 0.564, 0.496, 0.429, 0.366, 0.304, 0.244,
     1          0.187, 0.094, 0.007, -0.035, -0.074, -0.152, -0.225,
     1         -0.294, -0.361, -0.393, -0.425, -0.485, -0.544, -0.601,
     2         -0.654, -0.707, -0.766, -0.823, -0.959, -1.085, -1.202,
     3         -1.312, -1.416, -1.515, -1.609, -1.698, -1.783, -1.865 /
      data a2 / -0.968, -0.91, -0.88, -0.968, -0.9818, -0.9869, -0.9952, -1.0066, -1.0095, -1.0152, 
     1          -1.0235, -1.0262, -1.0367, -1.0416, -1.0463, -1.0486, -1.0509, -1.0552, -1.0593, 
     2          -1.0685, -1.0764, -1.0829, -1.0852, -1.0882, -1.0921, -1.0949, -1.0965, -1.097, 
     3          -1.0962, -1.0941, -1.0872, -1.0777, -1.0668, -1.0634, -1.0551, -1.0431, -1.0309, 
     4          -1.0188, -1.0069, -0.9952, -0.9837, -0.9616, -0.9406, -0.9305, -0.9207, -0.9016, 
     5          -0.8924, -0.8834, -0.8781, -0.8732, -0.8709, -0.8686, -0.8643, -0.8601, -0.8562, 
     6          -0.8524, -0.8506, -0.8488, -0.8454, -0.8421, -0.8344, -0.8273, -0.8209, -0.8188, 
     7          -0.8149, -0.8093, -0.8041, -0.7992, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946,
     8          -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, 
     9          -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, 
     1          -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, 
     1          -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946, -0.7946 /

      data a3 / 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     1          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     2          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     3          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     4          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     5          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     6          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     7          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     8          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     9          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     1          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
     1          0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265 /

      data a4 / -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     1     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     2     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     3     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     4     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     5     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     6     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     7     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     8     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     9     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     1     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     1     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     2     -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
     3     -0.231, -0.231, -0.231, -0.231 /

      data a5 / -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
     1     -0.398, -0.398, -0.398, -0.398 /

      data a6 / -1.44, 0.0, 0.0, -1.445, -1.4522, -1.4571, -1.4612, -1.4672, -1.4708, 
     1     -1.4774, -1.4862, -1.4891, -1.4964, -1.4998, -1.5043, -1.5079, 
     1     -1.5108, -1.517, -1.523, -1.5374, -1.5489, -1.5613, -1.5646, 
     1     -1.5681, -1.5776, -1.5806, -1.5863, -1.5854, -1.5778, -1.5814, 
     1     -1.5797, -1.5595, -1.5487, -1.5469, -1.5426, -1.5342, -1.5109, 
     1     -1.4782, -1.4512, -1.447, -1.4433, -1.4284, -1.4027, -1.3899, 
     1     -1.3823, -1.367, -1.3558, -1.3499, -1.3349, -1.2933, -1.2765, 
     1     -1.26, -1.2475, -1.2559, -1.2797, -1.293, -1.2965, -1.2991, 
     1     -1.2982, -1.2904, -1.3051, -1.316, -1.3201, -1.3282, -1.3393, 
     1     -1.3497, -1.3851, -1.4046, -1.4074, -1.4034, -1.3932, -1.3681, 
     1     -1.3589, -1.3236, -1.2717, -1.2353, -1.2152, -1.0434, -1.0328, 
     1     -1.0274, -1.0294, -0.9823, -0.9892, -0.9739, -0.9596, -0.951, 
     1     -0.9342, -0.8645, -0.8872, -0.8958, -0.8995, -0.8827, -0.8856,
     1     -1.0348, -1.0212, -0.9985, -0.9753, -0.9638, -1.0028, -0.995, 
     1     -0.9702, -0.9565, -0.8952, -0.8398, -0.7739, -0.7367, -0.696, -0.6738 /

      data a7 / 1.052, 0.0, 0.0, 1.1342, 1.108, 1.1018, 1.0817, 1.0532, 1.0538, 
     1     1.0534, 1.051, 1.0543, 1.0528, 1.0479, 1.0492, 1.0552, 1.0596, 
     1     1.0696, 1.0799, 1.101, 1.1179, 1.1465, 1.152, 1.1599, 1.1908, 
     1     1.1995, 1.2268, 1.2231, 1.1967, 1.2247, 1.2514, 1.1943, 1.1923, 
     1     1.2006, 1.2071, 1.2259, 1.1671, 1.0688, 0.9916, 1.0116, 1.0363,
     1     1.0702, 1.0293, 1.0022, 0.9936, 0.9791, 0.9577, 0.9576, 0.9298,
     1     0.7901, 0.7361, 0.6788, 0.6406, 0.6815, 0.7915, 0.8621, 0.882,
     1     0.8975, 0.8996, 0.872, 0.9525, 1.0075, 1.0243, 1.0537, 1.1042,
     1     1.1578, 1.3146, 1.4079, 1.4183, 1.3971, 1.3486, 1.239, 1.2047, 
     1     1.0786, 0.8978, 0.7434, 0.6499, 0.0272, -0.0402, -0.0778, 
     1     -0.0837, -0.2088, -0.1682, -0.2436, -0.3096, -0.1603, -0.2573, 
     1     -0.5093, -0.4603, -0.448, -0.455, -0.5916, -0.6181, 0.1125, 
     1     0.03, -0.0984, -0.2211, -0.2868, -0.0427, -0.1023, -0.2274, 
     1     -0.3748, -0.6794, -0.9725, -0.6304, -0.8524, -1.0953, -1.2714 /

      data a8 / 0.0, -0.1, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0008, 0.0038, 0.0052,
     1      0.0065, 0.0072, 0.0078, 0.009, 0.0102, 0.0129, 0.0154, 0.0177, 
     1      0.0186, 0.0198, 0.0218, 0.0236, 0.0254, 0.027, 0.0268, 0.0266, 
     1      0.0262, 0.0258, 0.0255, 0.0254, 0.0251, 0.0206, 0.0164, 0.0124, 
     1      0.0087, 0.0052, 0.0018, -0.0044, -0.0101, -0.0128, -0.0153, -0.0202, 
     1      -0.0225, -0.0247, -0.0289, -0.0329, -0.0348, -0.0366, -0.0402, -0.0435, 
     1      -0.0467, -0.0497, -0.0512, -0.0527, -0.0554, -0.0581, -0.0651, 
     1          -0.0715, -0.0774, -0.0793, -0.0829, -0.088, -0.0927, -0.0972, 
     1          -0.1014, -0.1054, -0.1092, -0.1162, -0.1226, -0.1285, -0.134, 
     1          -0.1391, -0.1438, -0.1483, -0.1525, -0.1565, -0.1603, -0.1681, 
     1          -0.1753, -0.1786, -0.1819, -0.188, -0.1936, -0.199, -0.2039, 
     1          -0.2063, -0.2086, -0.2131, -0.2173, -0.2213, -0.2252, -0.2288, 
     1          -0.2323, -0.2357, -0.2435, -0.2507, -0.2573, -0.2634, -0.2691, 
     1          -0.2744, -0.2794, -0.2841, -0.2885, -0.2927 /

      data a9 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     2          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     3          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     4          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     5          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     6          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     7          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     8          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     9          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

      data a10 / 0.9485, 1.587, -0.9, 0.9485, 0.9874, 1.0028, 1.024, 1.0464, 
     1         1.0511, 1.0606, 1.0724, 1.0771, 1.0924, 1.1007, 1.1101, 1.1137, 
     1         1.1184, 1.1278, 1.1373, 1.1621, 1.1904, 1.2199, 1.2317, 1.2517, 
     1         1.2848, 1.319, 1.3544, 1.391, 1.4276, 1.4653, 1.5397, 1.6152, 
     1         1.6883, 1.7069, 1.7504, 1.8107, 1.8703, 1.9257, 1.9803, 2.0306, 
     1         2.08, 2.1679, 2.2461, 2.2814, 2.3134, 2.3719, 2.3972, 2.4216, 
     1         2.4611, 2.4941, 2.5084, 2.5204, 2.5388, 2.5516, 2.5589, 2.5604,
     1          2.5602, 2.5576, 2.5514, 2.5395, 2.4975, 2.4354, 2.3598, 2.3308, 
     1          2.272, 2.1643, 2.0496, 1.9313, 1.808, 1.6833, 1.5581, 1.3038, 
     1          1.0531, 0.809, 0.5725, 0.348, 0.1341, -0.0668, -0.2538, -0.4281, 
     1          -0.5887, -0.84, -0.9415, -0.9415, -0.9415, -0.9415, -0.9415, 
     1          -0.9415, -0.9415, -0.9415, -0.9415, -0.9415, -0.9415, -0.9415, 
     1          -0.9415, -0.9415, -0.9263, -0.9118, -0.8779, -0.8469, -0.8184, 
     1          -0.792, -0.7675, -0.7445, -0.7229, -0.7026, -0.6833, -0.6651 /

      data a11 / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     2           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     3           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     4           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     5           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     6           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     7           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     8           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     9           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

      data a12 / -0.12, -0.05, -0.05, -0.12, -0.12, -0.12, -0.12, 
     1         -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, 
     1         -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, 
     1         -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.12, -0.1171, 
     1         -0.1145, -0.112, -0.1113, -0.1098, -0.1077, -0.1057, 
     1         -0.1039, -0.1021, -0.1005, -0.0989, -0.096, -0.0934, 
     1         -0.0921, -0.091, -0.0887, -0.0876, -0.0866, -0.0846, 
     1         -0.0828, -0.0819, -0.0811, -0.0794, -0.0779, -0.0764, 
     1         -0.075, -0.0743, -0.0736, -0.0723, -0.0711, -0.0682, 
     1         -0.0655, -0.0631, -0.0623, -0.0608, -0.0587, -0.0568, 
     1         -0.0549, -0.0532, -0.0516, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05 /

      data a13 / -0.05, -0.2, -0.2, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05,
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05,
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, 
     1         -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.06, -0.0632,
     1          -0.0692, -0.0778, -0.0858, -0.0934, -0.1005, -0.1073, 
     1          -0.1136, -0.1255, -0.1364, -0.1463, -0.1556, -0.1642, 
     1          -0.1722, -0.1798, -0.1869, -0.1936, -0.2, -0.2, -0.2, 
     1          -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 
     1          -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 
     1          -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2 /

      data a14 / 1.08, 0.7, 0.5, 1.08, 1.08, 1.0925, 1.1092, 1.1287, 
     1         1.1331, 1.1416, 1.1533, 1.157, 1.1708, 1.1772, 1.1833, 
     1         1.1862, 1.1891, 1.1947, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 
     1         1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.1886, 1.1854,
     1         1.1781, 1.1683, 1.1591, 1.1505, 1.1424, 1.1347, 1.1274, 
     1         1.1138, 1.1015, 1.0956, 1.0901, 1.0795, 1.0745, 1.0697,
     1         1.0605, 1.0519, 1.0478, 1.0438, 1.0361, 1.0288, 1.0219,
     1         1.0153, 1.0121, 1.009, 1.0029, 0.9971, 0.9835, 0.9712,
     1         0.9598, 0.9561, 0.9493, 0.9395, 0.9303, 0.9217, 0.9135,
     1         0.9058, 0.8985, 0.885, 0.8726, 0.8612, 0.8507, 0.8409,
     1         0.8317, 0.8231, 0.815, 0.8073, 0.8, 0.7246, 0.6558, 
     1         0.6235, 0.5925, 0.5339, 0.4793, 0.4283, 0.3804, 0.3574,
     1         0.3352, 0.2924, 0.2518, 0.2133, 0.1765, 0.1413, 0.1077,
     1         0.0754, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

      data a15 / -0.405, -0.284, 0, -0.4018, -0.4098, -0.4133, -0.4186, 
     1         -0.4277, -0.4295, -0.432, -0.4345, -0.4377, -0.4478, -0.451, 
     1         -0.4533, -0.4558, -0.4587, -0.4675, -0.4723, -0.4829, -0.4879,
     1         -0.4956, -0.5032, -0.5056, -0.5036, -0.5085, -0.5153, -0.5248,
     1         -0.5223, -0.5178, -0.5138, -0.5257, -0.537, -0.5364, -0.547, 
     1         -0.5227, -0.4967, -0.4784, -0.4685, -0.4587, -0.4407, -0.3986,
     1         -0.3678, -0.3646, -0.3631, -0.3803, -0.3942, -0.4061, -0.3844,
     1         -0.3674, -0.3601, -0.3543, -0.3453, -0.349, -0.3499, -0.3535,
     1         -0.3541, -0.3576, -0.361, -0.3614, -0.3573, -0.3552, -0.3612, 
     1         -0.345, -0.329, -0.3153, -0.3089, -0.2754, -0.2565, -0.2557, 
     1         -0.2591, -0.2289, -0.2198, -0.2159, -0.1606, -0.1626, -0.1738, 
     1         -0.0893, -0.0899, -0.112, -0.1236, -0.0906, -0.092, -0.0602, 
     1         -0.0434, -0.009, -0.0225, -0.086, -0.0818, -0.0936, -0.111, 
     1         -0.1482, -0.1531, 0.0879, 0.0715, 0.0409, -0.0034, -0.011, 
     1         -0.0711, -0.0503, -0.0237, 0.0069, -0.0099, -0.0227, -0.1238, 
     1         -0.2288, -0.2924, -0.3319 /

      data a16 / 0.65, 0.46, 0, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 
     1         0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 
     1         0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 
     1         0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 
     1         0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 
     1         0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 
     1         0.6202, 0.593, 0.568, 0.5599, 0.5448, 0.5233, 0.5031, 0.4841, 
     1         0.4663, 0.4494, 0.4333, 0.4035, 0.3763, 0.3513, 0.3282, 0.3066, 
     1         0.2864, 0.2675, 0.2496, 0.2327, 0.2167, 0.1869, 0.1597, 0.1469, 
     1         0.1347, 0.1115, 0.0899, 0.0698, 0.0508, 0.0417, 0.0329, 0.016, 
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1         0.0, 0.0, 0.0, 0.0, 0.0 /

      data a17 / 0.6, 0.35, 0.0, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 
     1         0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
     1         0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,
     1         0.6, 0.6, 0.6, 0.5789, 0.5596, 0.5506, 0.5419, 0.5255, 0.5177, 
     1         0.5102, 0.4959, 0.4824, 0.476, 0.4698, 0.4578, 0.4464, 0.4356, 
     1         0.4253, 0.4203, 0.4155, 0.406, 0.397, 0.3759, 0.3566, 0.3389, 
     1         0.3331, 0.3224, 0.3071, 0.2929, 0.2794, 0.2668, 0.2548, 0.2434, 
     1         0.2223, 0.203, 0.1853, 0.1689, 0.1536, 0.1393, 0.1258, 0.1132, 
     1         0.1012, 0.0898, 0.0687, 0.0494, 0.0404, 0.0317, 0.0153, 0.0, 0.0, 
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

      data a18 / -0.0067, 0.0, 0.0, -0.0067, -0.0067, -0.0067, -0.0067, -0.0067,
     1         -0.0067, -0.0067, -0.0067, -0.0067, -0.0067, -0.0069, -0.0071, 
     1         -0.0072, -0.0073, -0.0075, -0.0076, -0.008, -0.0084, -0.0087, 
     1         -0.0088, -0.009, -0.0093, -0.0093, -0.0093, -0.0093, -0.0093, 
     1         -0.0093, -0.0093, -0.0093, -0.0093, -0.0093, -0.0093, -0.0093, 
     1         -0.0093, -0.0093, -0.0089, -0.0086, -0.0083, -0.0077, -0.0071, 
     1         -0.0069, -0.0066, -0.0062, -0.006, -0.0057, -0.0053, -0.005, 
     1         -0.0048, -0.0046, -0.0043, -0.0039, -0.0036, -0.0033, -0.0032, 
     1         -0.0031, -0.0028, -0.0025, -0.0019, -0.0014, -0.0009, -0.0007, 
     1         -0.0004, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
     1         0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /

      data s1e / 0.576, 0.612, 0.733336393, 0.576, 0.575, 0.577, 0.579, 
     1           0.582, 0.582, 0.583, 0.585, 0.585, 0.587, 0.588, 0.589, 0.589,
     1           0.589, 0.59, 0.591, 0.592, 0.594, 0.595, 0.596, 0.596, 0.598, 
     1           0.599, 0.6, 0.601, 0.602, 0.603, 0.605, 0.607, 0.609, 0.61, 
     1           0.611, 0.613, 0.615, 0.616, 0.618, 0.619, 0.621, 0.623, 0.626,
     1           0.627, 0.628, 0.63, 0.631, 0.632, 0.633, 0.634, 0.635, 0.636,
     1           0.637, 0.638, 0.639, 0.64, 0.641, 0.641, 0.642, 0.643, 0.645,
     1           0.647, 0.649, 0.649, 0.65, 0.652, 0.654, 0.655, 0.657, 0.658,
     1           0.66, 0.662, 0.664, 0.666, 0.668, 0.67, 0.672, 0.673, 0.675, 
     1           0.676, 0.677, 0.678, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679,
     1           0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679,
     1           0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 0.679, 
     1           0.679, 0.679, 0.679 /

      data s2e / 0.456, 0.489, 0.739515254, 0.456, 0.455, 0.457, 0.459, 0.462,
     1           0.462, 0.464, 0.465, 0.466, 0.467, 0.468, 0.469, 0.469, 0.47, 0.47,
     1           0.471, 0.473, 0.474, 0.476, 0.476, 0.477, 0.478, 0.479, 0.481, 
     1           0.482, 0.483, 0.484, 0.486, 0.488, 0.49, 0.49, 0.492, 0.493, 
     1           0.495, 0.497, 0.498, 0.499, 0.501, 0.503, 0.505, 0.506, 0.507,
     1           0.509, 0.51, 0.511, 0.512, 0.514, 0.514, 0.515, 0.516, 0.517, 
     1           0.518, 0.519, 0.52, 0.52, 0.521, 0.522, 0.524, 0.526, 0.527, 
     1           0.528, 0.529, 0.53, 0.532, 0.533, 0.534, 0.535, 0.536, 0.538, 
     1           0.54, 0.541, 0.543, 0.544, 0.545, 0.547, 0.548, 0.549, 0.55, 
     1           0.552, 0.553, 0.554, 0.555, 0.556, 0.558, 0.559, 0.56, 0.561, 
     1           0.561, 0.563, 0.567, 0.572, 0.575, 0.579, 0.583, 0.586, 0.594, 
     1           0.601, 0.608, 0.614, 0.619, 0.625, 0.63, 0.634, 0.639, 0.643 /

      data s3 / 0.479, 0.416, 0.686401383, 0.479, 0.479, 0.481, 0.482, 0.484, 
     1          0.485, 0.486, 0.487, 0.488, 0.489, 0.49, 0.491, 0.491, 0.491, 0.492, 
     1          0.493, 0.494, 0.496, 0.497, 0.497, 0.498, 0.499, 0.5, 0.501, 0.502,
     1          0.503, 0.503, 0.505, 0.506, 0.507, 0.507, 0.508, 0.509, 0.51, 
     1          0.511, 0.512, 0.512, 0.513, 0.515, 0.516, 0.517, 0.517, 0.518, 
     1          0.519, 0.519, 0.52, 0.521, 0.521, 0.522, 0.522, 0.523, 0.524,
     1          0.524, 0.525, 0.525, 0.525, 0.526, 0.527, 0.528, 0.529, 0.529, 
     1          0.53, 0.531, 0.531, 0.532, 0.532, 0.533, 0.533, 0.534, 0.535, 
     1          0.535, 0.536, 0.536, 0.537, 0.537, 0.538, 0.538, 0.538, 0.539,
     1          0.539, 0.539, 0.54, 0.54, 0.541, 0.541, 0.541, 0.541, 0.541, 
     1          0.542, 0.542, 0.542, 0.543, 0.543, 0.543, 0.543, 0.544, 0.544,
     1          0.544, 0.545, 0.545, 0.545, 0.545, 0.546, 0.546, 0.546 /

      data s4 / 0.286, 0.262, 0.441385296, 0.286, 0.286, 0.288, 0.289, 0.291,
     1          0.292, 0.293, 0.294, 0.295, 0.296, 0.297, 0.298, 0.298, 0.299,
     1          0.299, 0.3, 0.302, 0.303, 0.304, 0.305, 0.306, 0.307, 0.308, 
     1          0.309, 0.31, 0.31, 0.311, 0.312, 0.314, 0.315, 0.315, 0.316, 
     1          0.316, 0.317, 0.317, 0.318, 0.318, 0.319, 0.32, 0.32, 0.321, 
     1          0.321, 0.322, 0.322, 0.322, 0.323, 0.323, 0.324, 0.324, 0.324,
     1          0.325, 0.325, 0.326, 0.326, 0.326, 0.326, 0.327, 0.327, 0.328,
     1          0.328, 0.329, 0.329, 0.329, 0.33, 0.33, 0.33, 0.331, 0.331, 
     1          0.331, 0.332, 0.332, 0.333, 0.333, 0.333, 0.334, 0.334, 0.334,
     1          0.334, 0.335, 0.335, 0.335, 0.335, 0.335, 0.336, 0.336, 0.336,
     1          0.336, 0.336, 0.336, 0.336, 0.336, 0.336, 0.336, 0.336, 0.336,
     1          0.337, 0.337, 0.337, 0.337, 0.337, 0.337, 0.337, 0.337, 0.337,
     1          0.337 /

      data s1m / 0.562, 0.578, 0.688678637, 0.562, 0.561, 0.563, 
     1         0.565, 0.568, 0.568, 0.57, 0.571, 0.572, 0.573, 0.574, 0.575, 
     1         0.575, 0.576, 0.576, 0.577, 0.579, 0.58, 0.582, 0.582, 0.583, 
     1         0.584, 0.585, 0.587, 0.588, 0.589, 0.59, 0.592, 0.594, 0.596, 
     1         0.597, 0.598, 0.599, 0.6, 0.601, 0.603, 0.604, 0.605, 0.606, 
     1         0.608, 0.609, 0.61, 0.611, 0.612, 0.613, 0.614, 0.614, 0.615, 
     1         0.615, 0.616, 0.617, 0.617, 0.618, 0.618, 0.619, 0.619, 0.62, 
     1         0.621, 0.622, 0.623, 0.624, 0.624, 0.624, 0.625, 0.625, 0.625, 
     1         0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 0.625, 
     1         0.624, 0.624, 0.624, 0.625, 0.625, 0.626, 0.626, 0.626, 0.626, 
     1         0.626, 0.626, 0.626, 0.626, 0.626, 0.626, 0.626, 0.626, 0.626, 
     1         0.627, 0.629, 0.633, 0.636, 0.639, 0.642, 0.644, 0.646, 0.648, 
     1         0.65, 0.651, 0.653 /

      data s2m / 0.438, 0.445, 0.695254494, 0.438, 0.437, 0.439, 
     1         0.441, 0.444, 0.445, 0.446, 0.448, 0.448, 0.45, 0.451, 0.452, 
     1         0.452, 0.452, 0.453, 0.454, 0.456, 0.457, 0.459, 0.459, 0.46, 
     1         0.461, 0.463, 0.464, 0.465, 0.466, 0.467, 0.469, 0.471, 0.473, 
     1         0.474, 0.475, 0.476, 0.477, 0.478, 0.479, 0.48, 0.48, 0.482, 
     1         0.483, 0.484, 0.485, 0.486, 0.487, 0.487, 0.488, 0.489, 0.489, 
     1         0.489, 0.49, 0.491, 0.491, 0.492, 0.492, 0.492, 0.493, 0.493, 
     1         0.494, 0.495, 0.496, 0.496, 0.497, 0.496, 0.496, 0.495, 0.495, 
     1         0.494, 0.494, 0.492, 0.491, 0.49, 0.489, 0.488, 0.487, 0.486, 
     1         0.485, 0.483, 0.482, 0.485, 0.487, 0.487, 0.488, 0.49, 0.491, 
     1         0.493, 0.494, 0.495, 0.495, 0.498, 0.503, 0.507, 0.511, 0.516, 
     1         0.522, 0.527, 0.541, 0.552, 0.563, 0.572, 0.581, 0.589, 0.596, 
     1         0.603, 0.61, 0.616 /
      data sigcorr / 1.000, 0.74, 0.47, 1.000, 1.000, 0.998, 0.996, 0.992,
     1               0.991, 0.989, 0.987, 0.986, 0.982, 0.98, 0.978, 0.978,
     2               0.977, 0.975, 0.973, 0.969, 0.965, 0.961, 0.96, 0.957, 
     3               0.952, 0.946, 0.942, 0.937, 0.933, 0.929, 0.921, 0.914,
     4               0.908, 0.906, 0.902, 0.896, 0.891, 0.886, 0.882, 0.878,
     5               0.874, 0.866, 0.859, 0.856, 0.853, 0.847, 0.844, 0.841,
     6               0.836, 0.831, 0.829, 0.827, 0.823, 0.818, 0.815, 0.811, 
     7               0.81, 0.804, 0.794, 0.783, 0.759, 0.737, 0.717, 0.71, 
     8               0.698, 0.68, 0.664, 0.648, 0.634, 0.62, 0.607, 0.583, 
     9               0.561, 0.541, 0.522, 0.504, 0.488, 0.472, 0.458, 0.444, 
     1               0.431, 0.407, 0.385, 0.374, 0.364, 0.346, 0.328, 0.312, 
     2               0.296, 0.289, 0.282, 0.268, 0.255, 0.243, 0.231, 0.22,
     3               0.209, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 /
      data taucorr / 1.000, 0.74, 0.47, 1.000, 1.000, 0.998, 0.996, 0.992,
     1               0.991, 0.989, 0.987, 0.986, 0.982, 0.98, 0.978, 0.978,
     2               0.977, 0.975, 0.973, 0.969, 0.965, 0.961, 0.96, 0.957, 
     3               0.952, 0.946, 0.942, 0.937, 0.933, 0.929, 0.921, 0.914,
     4               0.908, 0.906, 0.902, 0.896, 0.891, 0.886, 0.882, 0.878,
     5               0.874, 0.866, 0.859, 0.856, 0.853, 0.847, 0.844, 0.841,
     6               0.836, 0.831, 0.829, 0.827, 0.823, 0.818, 0.815, 0.811, 
     7               0.81, 0.804, 0.794, 0.783, 0.759, 0.737, 0.717, 0.71, 
     8               0.698, 0.68, 0.664, 0.648, 0.634, 0.62, 0.607, 0.583, 
     9               0.561, 0.541, 0.522, 0.504, 0.488, 0.472, 0.458, 0.444, 
     1               0.431, 0.407, 0.385, 0.374, 0.364, 0.346, 0.328, 0.312, 
     2               0.296, 0.289, 0.282, 0.268, 0.255, 0.243, 0.231, 0.22,
     3               0.209, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 /

      n = 1.18
      c = 1.88
      sigamp = 0.30

C Find the requested spectral period and corresponding coefficients
c      nPer = 3
      i1=iper
C  PGA, PGV, PGD cases  1,2,3 respectively
         period1 = period(i1)
         a1T = a1(i1)
         a2T = a2(i1)
         a3T = a3(i1)
         a4T = a4(i1)
         a5T = a5(i1)
         a6T = a6(i1)
         a7T = a7(i1)
         a8T = a8(i1)
         a9T = a9(i1)
         a10T = a10(i1)
         a11T = a11(i1)
         a12T = a12(i1)
         a13T = a13(i1)
         a14T = a14(i1)
         a15T = a15(i1)
         a16T = a16(i1)
         a17T = a17(i1)
         a18T = a18(i1)
         c4T = c4(i1)
         b_soilT = b_soil(i1)
         vLinT   = vLin(i1)
         if (vs30_class .eq. 0) then
            s1T = s1e(i1)
            s2T = s2e(i1)
         elseif (vs30_class .eq. 1) then
            s1T = s1m(i1)
            s2T = s2m(i1)
         endif
         s3T = s3(i1)
         s4T = s4(i1)
         taucorrT = taucorr(i1)
         sigcorrT = sigcorr(i1)

c     Set distance (Eq. 3)
      R = sqrt(rRup**2 + c4T**2)

C     No longer : Compute Rx value given source to site azimuth
c      if (srcsiteA .eq. 90.0) then         
c         Rx = Rrup/sin(dip*3.14159/180.0) - DepthTop/tan(dip*3.14159/180.0)
c      elseif (srcsiteA .gt. 0.0) then
c         Rx = Rjb*tan(srcSiteA*3.14159/180.0)
c      elseif (srcsiteA .lt. 0.0) then
c         Rx = 0.0
c      endif
       Rx=R_x
C     Set the 5 taper models (Eq. 8-12 SR)
C     Taper 1 (Eq. 8)
      if ( rjb .ge. 30. ) then
         taper1 = 0.
      else
         taper1 = (30. - rjb) / 30.
      endif

c     Taper 2 (Eq. 9)
      if ( dip .eq. 90.0 .or. Rx .gt. width*cos(dip*3.14159/180.0) ) then
         taper2 = 1.
       else
         taper2 = 0.5 + Rx/(2.0*width*cos(dip*3.14159/180.0))
       endif

c     Set Taper 3 (Eq. 10)
      if (Rx .ge. DepthTop) then
         taper3 = 1.0
      else
         taper3 = Rx/depthTop      
      endif

c     Set Taper 4 (Eq. 11)
      if ( mag .ge. 7. ) then
        taper4 = 1.0
      elseif ( mag .gt. 6.0) then
        taper4 = (mag-6.0)
      else
        taper4 = 0.0
      endif

c     Set Taper 5 (Eq. 12)
c     Modified by SR based on 2009 errata
      if ( dip .lt. 30. ) then
        taper5 = 1.
      else
        taper5 = 1. - (dip-30.)/60.
      endif
              
C     Base Model (Eq. 2)
      c1 = 6.75
       if ( mag .le. c1 ) then
        sum = a1T + a4T*(Mag-c1) + a8T*(8.5-Mag)**2 + (a2T + a3T*(Mag-c1))*alog(R)
       else
        sum = a1T + a5T*(Mag-c1) + a8T*(8.5-Mag)**2 + (a2T + a3T*(Mag-c1))*alog(R)
       endif

c     style of faulting
      sum = sum + a12T*Frv + a13T*Fn
      
c     Site response
c     Set Velocity for break in slope (Eq. 6)
      if ( specT .ge. 2.0) then
          v1 = 700.
      elseif ( specT .gt. 1.00 ) then
          v1 = exp(6.76-0.297*alog(specT))
      elseif ( specT .gt. 0.50 ) then
          v1 = exp(8.0-0.795*alog(specT/0.21))
      elseif (specT .eq. -1.0) then
c      	SR: this is different in paper, 826	
          v1 = 700.0
      elseif (specT .eq. -2.0) then
          v1 = 700.0
      else
            v1 = 1500.
      endif

C     Set the Vs30* (Eq. 5)
      if ( vs30 .lt. v1 ) then
          vs = vs30
      else
          vs = v1
      endif      	

c     Compute site amplification (Eq. 4)  
      if (vs30 .lt. vLinT) then
          soilamp = a10T*alog(vs/vLinT) - 1.0*b_soilT*alog(c+pgaRock) 
     1              + b_soilT*alog(pgaRock+c*((vs/vLinT)**(n)) )
      else
             soilamp = (a10T + b_soilT*n) * alog(vs/vLinT)
      endif
      sum = sum + soilamp


C     Soil Depth Model
C     (Eq. 17 SR)
      if (vs30 .lt. 180) then
         zhat = exp(6.745)
      elseif (vs30 .le. 500.0) then
         zhat = exp(6.745-1.35*alog(vs30/180.0) )
      else
         zhat = exp(5.394-4.48*alog(vs30/500.0) )
      endif
C     (Eq. 19 SR)
      if (vs30 .gt. 1000 .or. specT .lt. 0.35) then
         e2 = 0.0
      elseif (specT .ge. 0.35 .and. specT .le. 2.0) then
         e2 = -0.25*alog(vs30/1000)*alog(specT/0.35)
      else
         e2 = -0.25*alog(vs30/1000)*alog(2.00/0.35)
      endif
C     (Eq. 20 SR)
      if (specT .ge. 2.0) then
         a22 = 0.0625*(specT-2.00)
      else
         a22 = 0.0
      endif
C     Test Value for Eq. 18 SR
      c2 = 50.0
      if (v1 .lt. 1000.0) then
         testv1 = v1
      else
         testv1 = 1000.0
      endif
      
      test = (a10T+b_soilT*n)*alog(vs/(testv1)) + 
     1        e2*alog((z10*1000.0+c2)/(zhat+c2))

C     (Eq. 18 SR)
      if (vs30 .ge. 1000.0) then 
         a21 = 0.0
      elseif (test .lt. 0.0) then
         a21 = -1.0*(a10T+b_soilT*n)*alog(vs/testv1) /
     1           alog((z10*1000.0+c2)/(zhat+c2))
      else
         a21 = e2
      endif

C     (Eq. 16 SR)
      if (z10*1000.0 .ge. 200) then
         sum = sum + a21*alog((z10*1000.0+c2)/(zhat+c2)) + 
     1               a22*alog(z10*1000.0/200.0)
      else
         sum = sum+ a21*alog((z10*1000.0+c2)/(zhat+c2))
      endif

c     Hanging wall Model (Eq. 7 SR)
      if ( HWFlag .eq. 1 ) then
        sum = sum + a14T * taper1 * taper2 * taper3 * taper4 * taper5
      endif
      
c     Depth to top of Rupture Model
c        Modified by SR to follow Eq. 13 of AS08 paper
      z0 = 2.
      z1 = 5.
      z2 = 10.       
        if (depthTop .lt. z2) then
        	sum = sum + a16T*depthTop/10.0
        else
      	sum = sum + a16T
      endif

       
C     Large Distance Model
c        Modified by SR to follow Eq. 15 of AS08 paper
      if (mag .lt. 5.5) then
         T6 = 1.0
      elseif (mag .le. 6.5) then
         T6 = 0.5*(6.5 - mag)+0.5
      else
         T6 = 0.5
      endif

C     (Eq. 14 SR)
       if (Rrup .gt. 100.0)sum = sum + a18T*(Rrup - 100.0)*T6 

C     Set the Sigma Values

c     Compute parital derivative of alog(soil amp) w.r.t. alog(rock PGA)
c     (Eq. 26) SR: This is in agreement with Errata 2009
      if ( vs30 .ge. vLinT) then
        dAmp_dPGA = 0.
      else
        dAmp_dPGA = b_soilT*pgaRock * ( -1. / (pgaRock+c) 
     1              + 1./ (pgaRock + c*(vs30/vLinT)**(n)) )
      endif

c     (Eq. 27 SR)
      if (mag .lt. 5.0) then
         sigmanot = s1T
         if (vs30_class .eq. 0) then
            sigmanotPGA = s1e(1)
         else
            sigmanotPGA = s1m(1)
         endif
      elseif (mag .le. 7.0) then
         sigmanot = s1T + ((s2T-s1T)/2.0)*(mag-5.0)
         if (vs30_class .eq. 0) then
            sigmanotPGA = s1e(1) + ((s2e(1)-s1e(1))/2.0)*(mag-5.0)
         else
            sigmanotPGA = s1m(1) + ((s2m(1)-s1m(1))/2.0)*(mag-5.0)
         endif
      else
         sigmanot = s2T
         if (vs30_class .eq. 0) then
            sigmanotPGA = s2e(1)
         else
            sigmanotPGA = s2m(1)
         endif
      endif

c     (Eq. 28 SR)
      if (mag .lt. 5.0) then
         taunot = s3T
         taunotPGA = s3(1)
      elseif (mag .le. 7.0) then
         taunot = s3T + ((s4T-s3T)/2.0)*(mag-5.0)
         taunotPGA = s3(1) + ((s4(1)-s3(1))/2.0)*(mag-5.0)
      else
         taunot = s4T
         taunotPGA = s4(1)
      endif

C     (Eq. 23 SR)
      sigmaB = sqrt(sigmanot*sigmanot - sigamp*sigamp)
      sigmaBPGA = sqrt(sigmanotPGA*sigmanotPGA - sigamp*sigamp)

C     (Eq. 24) SR: checks with Errata 2009
      sigma = sqrt( sigmaB**2 + sigAmp**2 + (dAmp_dPGA * sigmaBPGA)**2
     1        + 2. * dAmp_dPGA * sigmaBPGA * sigmaB * sigCorrT )

C     (Eq. 25)
      tau = sqrt( taunot**2 + (dAmp_dPGA * taunotPGA)**2
     1        + 2. * dAmp_dPGA * taunotPGA * taunot * tauCorrT )

c     Set total to return
      lnSa = sum

      return
      end subroutine AS_072007_model
      			
      subroutine AS2013_v10_model ( iper,mag, dip, FltWidth, ZTOR, Frv, Fn, rRup, rjb, Rx, Ry0, 
     1                     vs30, Sa1100, Z1, Z1_ref, hwflag, vs30_class, lnSa, phi, tau,useRy0)
 
      implicit none
      integer MAXPER     
      parameter (MAXPER=22)
c z1_ref added as input mar 13 2013. But NA is using a diff. fcn compared to Chiou's latest...
      logical useRy0
      real c4(MAXPER), a1(MAXPER), a2(MAXPER), a3(MAXPER), a4(MAXPER), a5(MAXPER), a6(MAXPER),
     1     a7(MAXPER), a8(MAXPER), a9(MAXPER), a10(MAXPER), a11(MAXPER),
     1     a12(MAXPER), a13(MAXPER), a14(MAXPER), a15(MAXPER), a17(MAXPER),
     1     a43(MAXPER), a44(MAXPER), a45(MAXPER), a46(MAXPER),
     2     s1(MAXPER), s2(MAXPER), s3(MAXPER), s4(MAXPER)
      real period(MAXPER), b(MAXPER), vLin(MAXPER)
      real M1, M2
      real lnSa, SA1100, rjb, rRup, Rx, Ry0, dip, mag, vs30
      real HW_taper1, HW_taper2, HW_taper3, HW_taper4, HW_taper5
      real damp_dSA1100, sigAmp, fltWidth, testv1
      real f1, f4, f5, f6, f7, f8, f9, f10
      real Ry1, ZTOR, Frv, Fn, SpecT
      real phiA, phiB, tauA, tauB, phi, tau
c      integer iper, vs30_class, hwflag, iPer
      integer iper, vs30_class, hwflag
      real n, c, Z1, zhat, c4_mag
      real R, V1, Vs30Star, hw_a2, h1, h2, h3, R1, R2, z1_ref

      data period / 0.00, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 
     1              0.3, 0.4, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7.5, 10 /
      data Vlin/ 660, 680, 770, 800, 800, 800, 740, 590, 495, 430, 360, 340, 330, 
     1           330, 330, 330, 330, 330, 330, 330, 330, 330 /
      data b/ -1.37, -1.35, -1.26, -1.16, -1.15, -1.23, -1.53, -1.86, -2.16, -2.411,
     1        -2.785, -3.006, -3.113, -2.851, -1.891, -0.792, 0, 0, 0, 0, 0, 0 /
      data c4/ 6, 6, 6, 6, 6, 5.9, 5.8, 5.7, 5.6, 5.5, 5.2, 4.8, 4.4, 4, 3.75, 3.5,
     1         3.25, 3, 3, 3, 3, 3 /
      data s1/ 0.7541, 0.7593, 0.7742, 0.8179, 0.826, 0.805, 0.796, 0.7788, 
     1         0.75, 0.731, 0.6934, 0.6682, 0.6209, 0.5848, 0.5632, 0.5452, 
     1          0.5302, 0.5152, 0.4956, 0.4755, 0.44, 0.4365 /
      data s2/ 0.5439, 0.5466, 0.5488, 0.5587, 0.5837, 0.6035, 0.6069, 0.6112,
     1          0.6148, 0.6171, 0.6123, 0.6206, 0.6395, 0.6618, 0.6574, 0.6563, 
     1          0.6231, 0.5996, 0.5834, 0.5834, 0.56, 0.5273 /
      data s3/ 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47,
     1         0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47 /
      data s4/ 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38,
     1         0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 0.38, 00.38 /

      data a1/ 1.0889, 1.1027, 1.2548, 1.5266, 1.753, 1.79, 2.0379, 2.2006, 2.2707, 
     1         2.28, 2.2564, 2.152, 1.8695, 1.5505, 1.0223, 0.6529, 0.1946, -0.1555, 
     2        -0.4919, -0.9061, -1.414, -2.3488/
      data a2/ -1.03, -1.03, -1.07, -1.09, -1.09, -1.07, -1.05, -1.02, -1, -0.98, 
     1         -0.96, -0.94, -0.91, -0.88, -0.84, -0.82, -0.78, -0.75, -0.72, -0.66, -0.595, -0.48 /
      data a3/ 0.281, 0.281, 0.281, 0.281, 0.278, 0.27, 0.258, 0.25, 0.242, 0.239, 
     1         0.231, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23, 0.23 /
      data a4/ -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1,
     1         -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1/
      data a5/ -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49,
     1         -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49, -0.49/
      data a6/ 2.2882, 2.2863, 2.2497, 2.1728, 2.1246, 2.1728, 2.2013, 2.2591, 2.3286, 
     1         2.3237, 2.4414, 2.5174, 2.6745, 2.7543, 2.8332, 2.8674, 2.8927, 2.8194, 2.9062, 
     2         3.0193, 3.1474, 3.4479 /
      data a7/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a8/ 0, 0, 0, 0, 0, 0, -0.029, -0.05, -0.066, -0.079, -0.099, -0.115, 
     1         -0.144, -0.165, -0.194, -0.214, -0.243, -0.264, -0.27, -0.27, -0.27, -0.27 /
      data a9/ 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75,
     1         6.75, 6.75, 6.75, 6.75, 6.82, 6.92, 7, 7.06, 7.145, 7.25 /
      data a10/ 1.56, 1.54, 1.43, 1.35, 1.3, 1.31, 1.57, 1.99, 2.39, 2.73, 3.27, 
     1          3.58, 3.74, 3.35, 1.91, 0.26, -0.93, -0.93, -0.93, -0.93, -0.93, -0.93 /
      data a11/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a12/ -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, 
     1          -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2, -0.2 /
      data a13/ 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.58, 0.56, 0.53,
     1          0.5, 0.42, 0.35, 0.2, 0, 0, 0, 0, 0 /
      data a14/ -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.24, -0.19, -0.11,
     1          -0.04, 0.07, 0.15, 0.27, 0.35, 0.46, 0.54, 0.61, 0.65, 0.72, 0.8 /
      data a15/ 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.03, 0.92, 0.84, 0.68, 
     1          0.57, 0.42, 0.31, 0.16, 0.05, -0.04, -0.11, -0.19, -0.3 /
      data a17/ -0.0025, -0.0025, -0.0025, -0.003, -0.0037, -0.004, -0.0037, 
     1          -0.0032, -0.0028, -0.0023, -0.0013, -0.0007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /
      data a43/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.08, 0.1, 0.1, 0.1, 0.13, 0.24, 
     1          0.3, 0.36, 0.43, 0.51, 0.61, 0.55, 0.42/
      data a44/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.04, 0.06, 0.11, 
     1          0.15, 0.2, 0.21, 0.25, 0.29, 0.32, 0.33, 0.32, 0.265, 0.21 /
      data a45/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.02, 0.04, 0.07, 0.13, 0.17, 
     1          0.2, 0.22, 0.25, 0.23, 0.21, 0.18, 0.155, 0.13 /
      data a46/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.03, 0.07, 0.11, 0.13, 
     1         0.14, 0.15, 0.16, 0.14, 0.12, 0.11, 0.07, 0.07, 0.07 /
     
      n = 1.5
      c = 1.8
     
c     Set distance (eq 3)
      R = sqrt(rRup**2 + c4(iper)**2)
               
C     Base Model (eq 2)
      M1 = a9(iper)
      M2 = 5.0
      if ( mag .le. M2 ) then
        f1 = a1(iper) + a6(iper)*(Mag-M2) + a7(iper)*(Mag-M2)**2 + a4(iper)*(M2-M1) + a8(iper)*(8.5-M2)**2 +
     1                (a2(iper) + a3(iper)*(M2-M1)) * alog(R) + a17(iper)*R
      elseif ( mag .le. M1 ) then
        f1 = a1(iper) + a4(iper)*(Mag-M1) + a8(iper)*(8.5-Mag)**2 + (a2(iper) + a3(iper)*(Mag-M1)) * alog(R) + a17(iper)*R
      else
        f1 = a1(iper) + a5(iper)*(Mag-M1) + a8(iper)*(8.5-Mag)**2 + (a2(iper) + a3(iper)*(Mag-M1)) * alog(R) + a17(iper)*R
      endif

c     style of faulting (eq 15 an 16) 
      if ( mag .lt. 4. ) then
        f7 = 0
        f8 = 0
      elseif ( mag .le. 5. ) then
        f7 = Frv * a11(iper) * (mag-4.)
        f8 = Fn * a12(iper) * (mag-4.)
      else 
        f7 = Frv * a11(iper)
        f8 = Fn * a12(iper)
      endif

c     ZTOR (eq 14)
      if (ZTOR .le. 20.) then 
        f6 = a15(iper) * ZTOR/20.0
      else
        f6 = a15(iper)
      endif    

c     Set VS30_star (eq 5)
      if ( period(iper) .ge. 3.0 ) then
        V1 = 800.
      elseif ( period(iper) .gt. 0.5 ) then
        V1 = exp( -0.351 * alog(period(iper)/0.5)  + alog(1500.) )
      else
        V1=1500.
      endif
      if ( vs30 .lt. v1 ) then
         vs30Star = vs30
      else
      vs30Star = v1
      endif      	

c     Compute site amplification (Eq. 4)  
      n = 1.5
      c = 1.8
      if (vs30 .lt. vLin(iPer)) then
        f5 = a10(iper)*alog(vs30Star/vLin(iper)) - b(iper)*alog(c+Sa1100) 
     1              + b(iper)*alog(Sa1100+c*((vs30Star/vLin(iper))**(n)) )
      else
           f5 = (a10(iper) + b(iper)*n) * alog(vs30Star/vLin(iper))
      endif


C     Soil Depth Model (eq 17)
      if ( vs30 .lt. 200. ) then 
        f10 = a43(iper) * alog( (z1+0.01) /(z1_ref+0.01) )
      elseif ( vs30 .lt. 300. ) then
        f10 = a44(iper) * alog( (z1+0.01) /(z1_ref+0.01) )
      elseif ( vs30 .lt. 500. ) then
        f10 = a45(iper) * alog( (z1+0.01) /(z1_ref+0.01) )
      elseif ( vs30 .lt. 700. ) then
        f10 = a46(iper) * alog( (z1+0.01) /(z1_ref+0.01) )
      else
        f10 = 0.
      endif        

c     Compute HW taper1 (eq 9) 
      if ( dip .lt. 30. ) then
        HW_taper1 = 60./ 45.
      else
        HW_taper1 = (90.-dip)/45.
      endif

c     Compute HW taper2 (eq. 10)
      hw_a2 = 0.2
      if( mag .gt. 6.5 ) then
        HW_taper2 = 1. + hw_a2 * (mag-6.5) 
      elseif ( mag .gt. 5.5 ) then
c the - sign below is a correction made on Mar 14 2013. 
        HW_taper2 = 1. + HW_a2 * (mag-6.5) - (1-HW_a2)*(mag-6.5)**2
c        HW_taper2 = 1. + HW_a2 * (mag-6.5) + (1-HW_a2)*(mag-6.5)**2      !wrong NA
      else
        HW_taper2 = 0.
      endif

c     Compute HW taper 3 (eq. 11)
      h1 = 0.25
      h2 = 1.5
      h3 = -0.75
      R1 = fltWidth * cos(dip*3.1415926/180.)
      R2 = 3.*R1
      if ( Rx .le. R1 ) then
        HW_taper3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)**2
      elseif ( Rx . lt. R2 ) then
        HW_taper3 = 1. - (Rx-R1)/(R2-R1)
      else
        HW_taper3 = 1.
      endif 

c     Compute HW taper 4 (eq 12)
      if ( ZTOR .lt. 10. ) then
        HW_taper4 = 1. - (ZTOR**2) / 100.
      else
        HW_taper4 = 0.
      endif
      
c     Compute HW taper 5 (eq. 13)  **** Ry0 version, not used here ***     
c      Ry1 = Rx * tan(20.*3.1415926/180.)
c      if ( Ry0 .lt. Ry1 ) then
c        HW_taper5 = 1.
c      elseif ( Ry0-Ry1 .lt. 5. ) then
c        HW_taper5 = 1. - (Ry0-Ry1) / 5.
c      else
c        HW_taper5 = 0.
c      endif

c     Compute HW taper 5 (eq. 13a)  **** No Ry0 version ***     
      if (Rjb .eq. 0. ) then
        HW_taper5 = 1. 
      elseif ( Rjb .lt. 30. ) then
        HW_taper5 = 1 - Rjb/30.
      else
        HW_taper5 = 0.
      endif

c     Hanging wall Model (eq 8)
      if ( HWFlag .eq. 1 ) then
        f4 = a13(iper) * HW_taper1 * HW_taper2 * HW_taper3 * HW_taper4 * HW_taper5
      else
        f4 = 0.
      endif
      
C     Set the Sigma Values
      
c     Compute within-event term, phiA, at the surface for linear site response (eq 23)
      if (mag .lt. 3.5) then
         phiA = s1(iper)
      elseif (mag .le. 7.0) then
         phiA = s1(iper) + ((s2(iper)-s1(iper))/3.5)*(mag-3.5)
      else
         phiA = s2(iper)
      endif

c     Compute between-event term, tau (eq. 24)
      if (mag .lt. 5.0) then
         tauA = s3(iper)
      elseif (mag .le. 7.0) then
         tauA = s3(iper) + ((s4(iper)-s3(iper))/2.0)*(mag-5.0)
      else
         tauA = s4(iper)
      endif
      tauB = tauA

c     Compute phiB, within-event term with site amp variablity removed (eq. 19)
      sigAmp = 0.4
      phiB = sqrt( phiA**2 - sigAmp**2)
      
c     Compute parital derivative of alog(soil amp) w.r.t. alog(SA1100) (eq. 22)
      if ( vs30 .ge. vLin(iper)) then
        dAmp_dSA1100 = 0.
      else
        dAmp_dSA1100 = b(iper)*SA1100 * ( -1. / (SA1100+c) 
     1              + 1./ (SA1100 + c*(vs30/vLin(iper))**(n)) )
      endif

C     Compute phi, with non-linear effects (eq. 20)
      phi = sqrt( phiB**2 * (1. + dAmp_dSA1100)**2 + sigAmp**2 )

C     Compute tau, with non-linear effects (eq. 21)
      tau = tauB * (1. + dAmp_dSA1100)
      
c     Compute median ground motion (eq. 1)
      lnSa = f1 + f4 + f5 + f6 + f7 + f8 + f9 + f10  

      return
      end subroutine AS2013_v10_model
      
