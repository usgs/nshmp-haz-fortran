c Combine_cms.2013.f   Fortran 95 or gfortran code last revised Dec 18 2013. Steve Harmsen.
c Dec 17 2013: revise GMM set to NSHMP 2013 hazard model.
c Aug 5 2011: add lms0 logic check because when CEUS and WUS were combined some extraeous periods
c were showing up in the output. Can only show 7 periods in mean spectrum for combined CEUS, WUS because
c some CEUS GMPEs only have SA defined at those seven periods.
c June 2011: add ltmp array which tells if element j has valid data (ltmp(j)=.true.)
c   for spectra with variably mixed on/off periods.
c combine the binned cms at a site using epistemic or other branch weights read in here.
c For ease of IO, the temporary input files are binary or unformatted.
c The output file is ascii for delivery to requestor.
c Modification May 25 2011: The DCMS are now written to two files, the
c  former for just plotting, the latter for serving to users... The former is
c just the mean binned CMS from index 0 (mean of GMPEs) while the other one
c contains the mean and the individual cases. These will be the same if the
c user does not request individual cases.
c
c Modification Apr 12 2011: The CMS for individual GMPEs are ranked by importance where
c importance here means relative frequency of contributions to the exceedance or band gm.
c
c The one wrinkle in this combine as of jan 28 2011 is that the WUS crustral sources
c have rates that are concentrating all probability into the central GM-uncert branch.
c The cms has to be reduced for those crustal sources.  Stop gap kludge. In the future, we
c may wish to explicitly include the low-weight upper and lower gm-uncert branches, 
c each w/weight 0.185. If/when that is done, wfac should become 1.0 for each file (i.e., lose this factor)
c the storage for primed rates are primerate and prate, respectively.
c
c SUN Compile:
c   f95 combine_cms.2013.f -o ../bin/combine_cms.2013 -e -C
c PC Compile with gnu-fortran:
c	gfortran combine_cms.2013.f -O -o combine_cms.2013 -static -ffixed-line-length-none
c
c Usage:
c   combine_cms.2013 control out1 out2 out3
c
c where the input file name is control. Control is the "addhaz" file associated with this
c deagg request. The addhaz file was presumably written by program combine.2011 during execution
c of script "deagg.wus." S 
c Information in record 1 of the *.bin files is compared to test compatibility.
c For use 
c
	character*80 fin,control,fout
	parameter (imm=44,np=21,small=1.e-11)
	real, dimension(np) :: st,ms_p
	real cms(np,0:40),dms(18,imm,10,np,0:40)
	integer kms(0:40),immin/40/,immax/2/,irmax/1/
	real pack(15000)	!for storing sorted rates of contributor GMPEs
	integer,dimension(15000):: indx,dakine	!tag along
c cms stores the mean of mean cms for each atten model, 0 = mean of all
c dms stores binned mean cms. The M index accomodates M5 to M9.3. For CEUS may need to modify,
c because M<5 sources are present (due to mblg->M conversion).
	real, dimension(0:40):: rbar,mbar,ebar,primerate,rate,sumwt
	real, dimension(10) :: epsband
c epsband center of band of width 0.5 epsilon, min value of -1.25, max is 2.75. 3.25 is there, but
c should never be seen for e0 capped at 3.0.
c brbar = binned rbar, bmbar = binned mbar, etc
	real, dimension(18,imm,10,0:40):: brbar,bmbar,bebar,brate,prate
	real x,xm,rcd,eps00,prx,vs0,vs30
	integer ir,im,ie,imax,i,j,k,m,n
	logical rwa(0:40)	! true if data present for atten model k
	logical, dimension(18,imm,10,0:40):: rua
	logical ltmp(21),lms0(21)
	logical isok/.true./	!was file found?
c m_att: index of atten model 0 for mean. 21 = BA NGA if all is well; 22=CB NGA; 23=CY NGA;...
c Here are the expected models and their slots:
      character*24 att_type(0:40)
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
     +'Chiou-Youngs NGA2013','Abr-Silva-Kamai 2013',
     +'IMIdriss NGA-W 2013','GraizerKalkan 2013','Atten Model 39','Atten Model 40'/
c Srcdesc: the source-type descriptions behind index "isrc". These src descriptions might be an output
c index at some time, in which case this block should be uncommented out.
c Spectral periods, ms_p. 0.01 s and PGA have equivalent response. Some periods are not available for
c several GMPEs. The spectrum will be truncated at the largest available T.
       ms_p = (/0.01, 0.020, 0.030, 0.050, 0.075, 0.100,
     + 0.150, 0.200, 0.250, 0.300, 0.400, 0.500, 0.750, 1.000,
     + 1.500, 2.000, 3.000, 4.000, 5.0, 7.5, 10.0/)
       epsband = (/-1.25,-0.75,-0.25,0.25,0.75,1.25,1.75,2.25,2.75,3.25/)
	if(iargc() .lt.4)stop 'Usage:  combine_cms control out1 out2 out3'
	call getarg(1,control)
	open(1,file=control,status='old')
	read(1,*)nfi
	read(1,246)
246	format(/)	!skip 2 records
c	Initializations...
	cms=0.	!exp(cms) is 1
	rwa = .false.
	rua = .false.
	lms0 = .true.
	sumwt = 0.0
	iamx=1
	kms=np	!initialize array size to max number of elements, currently 21
c	print *,'SA  T(s)  LAT  LONG Vs30 INFI#'
	do i=1,nfi
	read(1,50)fin
	isok = .false.
50	format(a)
	j=index(fin,' ')-1
	fin=fin(1:j)//'.bin'	!put "bin" extension on name to retrieve CMS files
        read (1,*)wt,isrc
        wfac=1.
c        if(isrc.gt.3.and.isrc.ne.23.and.isrc.lt.31)wfac=0.63	!this is the central branch wt for WUS crustal sources. Not 
c for Cascadia or CEUS. Cascadia isrc is in the range 1 to 3. CEUS isrc >30(not yet ready Jan 2010).
c Deep intraplate, isrc=23. Intraplate also has full weight GMPEs (no additional epistemic uncert)
c
        inquire(file=fin(1:j+4),exist = isok)
        if(.not.isok)then
        print *,fin, 'not found'
        goto 2
		  else
c		  print *,'Working on file: ',fin
        endif
	open(33,file=fin,form='unformatted')
        read(33)safix,period,rlatd,rlond,vs30
        print *,safix,period,rlatd,rlond,vs30
	if(i.eq.1)then
	sa0=safix; per0=period
	rlat0=rlatd; rlon0=rlond
	vs0=vs30
      else	!compare ith data with 1st
c      print *,safix,period,rlatd,rlond,vs30,i
      if (abs(safix-sa0).gt.0.0002)stop'fixed gm do not agree'
      if (abs(period-per0).gt.0.01)stop'spectral periods do not agree'
      if(abs(rlatd-rlat0).gt.0.0001.or.abs(rlond-rlon0).gt.0.0001)stop'lat or long do not agree'
      if(vs30.ne.vs0)stop' Site Vs30 do not agree with first file'
      endif
c the do loop with no added paraphenalia.
c
      do	
      read(33, end=8)m_att,ir,im,ie,rcd,xm,prx,eps00
c	if(m_att.gt.0.and.m_att.lt.3)print *,m_att,ir,im,ie,rcd,xm,prx,eps00
      iamx=max(iamx,m_att)
c       if (m_att.eq.13) print *,rcd,xm,eps00,prx
       if(prx.gt.small)then
c nuisance quantities (rates) like 1e-15 can gum up the works.
       rwa(m_att)=.true.
       immin=min(immin,im)
       immax=max(immax,im)	!maximum im index for this run
       irmax=max(ir,irmax)
       endif
	if(0.lt.ir.and.ir.lt.18 .and. 0.lt.im.and.im.lt.45
     + .and.0.lt.ie.and.ie.lt.11.and.prx.gt.small)then
       rua(ir,im,ie,m_att)=.true.
c	else
c	print *,'File ',i,' omit ',ir,im,ie,prx,m_att
	endif
      read(33)kmsp,st,ltmp		!temp storage used to output CMS
       lms0 = lms0 .and. ltmp	!only report mean CMS if all GMPEs have a report at period (T)
	do mm=1,kmsp
c Test for the occasional NaN. Below will fix it: interim solution. Should find why
c NaN is occurring in the first place. Presumption is that the hazard is very low (underflow related?)
	if(st(mm).ne.st(mm))st(mm)=-1e-18
	enddo
	prx=prx*wt
	rcd=rcd*prx
	xm=xm*prx
	eps00=eps00*prx
c	print *,m_att,kmsp,kms(m_att)
	kms(m_att)=min(kms(m_att),kmsp)
c unbinned version of CMS, except for attn model.
	rbar(m_att)=rbar(m_att) + rcd
	mbar(m_att)=mbar(m_att) +  xm
	ebar(m_att)=ebar(m_att) + eps00
	rate(m_att)=rate(m_att) + prx
	primerate(m_att)=primerate(m_att) +prx*wfac
c binned version
	if(rua(ir,im,ie,m_att))then
	brbar(ir,im,ie,m_att)=brbar(ir,im,ie,m_att) + rcd
	bmbar(ir,im,ie,m_att)=bmbar(ir,im,ie,m_att) +  xm
	bebar(ir,im,ie,m_att)=bebar(ir,im,ie,m_att) + eps00
	brate(ir,im,ie,m_att)=brate(ir,im,ie,m_att) + prx
	prate(ir,im,ie,m_att)=prate(ir,im,ie,m_att) + prx*wfac
	endif
	do j=1,kms(m_att)
	if(ltmp(j)) then
	cms(j,m_att)=cms(j,m_att)+wt*st(j)	!st already has prx multiplier when input
	if(rua(ir,im,ie,m_att))dms(ir,im,ie,j,m_att) = dms(ir,im,ie,j,m_att)+wt*st(j)
	endif !computable period, in which case ltmp() is true
	enddo	!periods of the cms
	enddo	!dowhile
8	close(33)
2	continue
c Having  accumulated all available data, go on to next input file
	enddo	!i=1,nfi
c Now it is time to write something. Two write files: the short and long versions.
	close(1)
	call getarg(2,fout)
	open(2,file=fout,status='unknown')
        write(2,907)safix,period,control,rlatd,rlond,vs30
907   format('#combine_cms MEAN of MEANs @fixed GM=',f6.3,' g; T=',f6.3,' s, control ',a,
     + /,'#site lat long = ',f8.4,1x,f9.4,' Vs30 ',f6.1,' meters/s') 
     	i=0 
     	do n=1,iamx
     	if(rwa(n))then
     	i=i+1
     	pack(i)=rate(n)
     	indx(i)=n
     	endif
     	enddo
     	imax=i	
     	if(imax.eq.0)then
     	goto 77
     	elseif(imax.gt.1)then
     	call sort2(imax,pack,indx)
     	endif
77	pr=rate(0)
     	fac=1/pr
	write(2,255)rbar(0)*fac,mbar(0)*fac,pr,ebar(0)*fac,att_type(0)
	fac =1/primerate(0)	!primerate from kludge.
c add lms0 check Aug 5 2011. S Harmsen.
       do j=1,kms(0)
       if(cms(j,0).ne.0..and.lms0(j))write(2,257)ms_p(j),exp(cms(j,0)*fac)
257	format(f7.3,1x,1pe11.4)
       enddo	!periods
c output is now sorted by contribution from GMPE (k) 
	if(imax.eq.0)then
	iamx=0
	goto 80 
	endif   	    
	do m=imax,1,-1	!top contributor is last but should appear first
	n=indx(m)
	fac=1./rate(n)
	pr=rate(n)
	write(2,255)rbar(n)*fac,mbar(n)*fac,pr,ebar(n)*fac,att_type(n)
255      format(/,'#T(s)  Mean_SA(g) for Rcdbar(km), Mbar, haz, eps0, GMPE: ',
     + f6.1,1x,f5.2,1x,e11.5,1x,f6.2,1x,a)
256      format(/,'#T(s)  Mean_SA(g) for Rcdbar(km), Mbar, haz, eps*, GMPE, eps0: ',
     + f6.1,1x,f5.2,1x,e11.5,1x,f6.2,1x,a,1x,f6.2)
        fac=1./primerate(n)	!the kludge. crustal weights diff from subduction.
       do j=1,kms(n)
       if(cms(j,n).ne.0.)write(2,257)ms_p(j),exp(cms(j,n)*fac)
       enddo	!periods
       enddo	!attn models
c next R,M,e0 binned output rua?  ready unwilling able?
80	close(2)
	call getarg(3,fout)
	open(2,file= fout,status='unknown')
	call getarg(4,fout)
	open(3,file= fout,status='unknown')
        write(2,917)safix,period,control,rlatd,rlond,vs30
        write(3,917)safix,period,control,rlatd,rlond,vs30
c	print *,'just wrote headerr, safix is ',safix
917   format('#combine_cms R,M,E BINNED MEANs @fixed GM=',f6.3,' g; T=',f6.3,' s, control ',a,
     + /,'#site lat long = ',f8.4,1x,f9.4,' Vs30 ',f6.1,' m/s')      
	pack=0.0
	i=0	
	do ir=1,irmax
	do im=immin,immax
	do ie=1,10
	if(rua(ir,im,ie,0))then
	i=i+1
	pack(i)=brate(ir,im,ie,0)
	indx(i)=ie + 100*im + 10000*ir
	print *,i,pack(i),indx(i)
	endif
	enddo
	enddo
	enddo
	imax=i
	call sort2(imax,pack,indx)
	dakine=indx
c	print *,dakine(1:3)
c
c write top 16 binned contributors : all gmpes, then corresponding for indiv. GMPEs if they are there
c for example subduction bin (ir,im,ie) may be entirely different from crustal event (ir,im,ie)
c
	print *,'number of bins to potentially report DCMS ',imax
	imin=max(1,imax-16)
	do m=imax,imin,-1
c	print *,m,dakine(m)
	ir=dakine(m)/10000
	irfac=ir*10000
	im = (dakine(m)-irfac)/100
	imfac=irfac+im*100
	ie=dakine(m)-imfac
c	print *,ir,im,ie,m,dakine(m)
	pr=brate(ir,im,ie,0)
	fac=1./pr
c replace epilon0 by a band value near zero for low-value episilon0
	e0=	bebar(ir,im,ie,0)*fac
	if(e0.lt.-1.)then
	e0p=-0.3
	elseif(e0.lt.-0.5)then
	e0p=-0.2
	elseif(e0.lt.0.)then
	e0p=0.2
	else
	e0p=e0
	endif
	write(2,255)brbar(ir,im,ie,0)*fac,bmbar(ir,im,ie,0)*fac,pr,e0p,att_type(0)
	write(3,255)brbar(ir,im,ie,0)*fac,bmbar(ir,im,ie,0)*fac,pr,e0p,att_type(0)
	fac=1./prate(ir,im,ie,0)
       do j=1,kms(0)
       if(dms(ir,im,ie,j,0).ne.0.)then
       write(2,257)ms_p(j),exp(dms(ir,im,ie,j,0)*fac)
       write(3,257)ms_p(j),exp(dms(ir,im,ie,j,0)*fac)
       endif
       enddo	!periods
        i=0
        if(iamx.eq.0)goto 88
c may be performing no-attn model analysis. in that case this part of code is bypassed.
	do n=1,iamx
c sort atten models which may or may not be making an appearance given ir,im, and ie. Report if at least
c 10% of max in a bin... omit the mini-contributors for brevity and put big bad boys on top.
	if(rua(ir,im,ie,n))then
	pr=brate(ir,im,ie,n)
	i=i+1
	pack(i)=pr
	indx(i)=n
	endif	!rua?
	enddo
	imax=i
	if(imax.gt.1)then
	call sort2(imax,pack,indx)
	endif
	x=0.05*pack(imax)
	do k=imax,1,-1
	n=indx(k)
	pr=pack(k)
	if(pr.lt.x)goto 88	!too small to worry about
	fac=1./pr
c replace epilon0 by a band value for low-value episilon0
	e0=	bebar(ir,im,ie,n)*fac
	if(e0.lt.-1.)then
	e0p=-0.3
	elseif(e0.lt.-0.5)then
	e0p=-0.2
	elseif(e0.lt.0.)then
	e0p=0.2
	else
	e0p=e0
	endif
	write(3,255)brbar(ir,im,ie,n)*fac,bmbar(ir,im,ie,n)*fac,pr,e0p,
     + att_type(n)
	fac=1./prate(ir,im,ie,n)
       do j=1,kms(n)
       if(dms(ir,im,ie,j,n).ne.0.)write(3,257)ms_p(j),exp(dms(ir,im,ie,j,n)*fac)
       enddo	!periods
       enddo	!k loop attn models
 88	continue
        enddo	! m (ordered rates) loop imax,imax-9
      close(2)
      close(3)
		end

      SUBROUTINE SORT2(N,RA,IB)
c sort the real values in RA and do the same sort to tag-along integer values in IB.
c these integer values are typically indexes  for accessing array elements. From Num Recipes.    
      real RA(N)
      integer IB(N)
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          RRA=RA(L)
          IRB=IB(L)
        ELSE
          RRA=RA(IR)
          IRB=IB(IR)
          RA(IR)=RA(1)
          IB(IR)=IB(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            RA(1)=RRA
            IB(1)=IRB
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(RA(J).LT.RA(J+1))J=J+1
          ENDIF
          IF(RRA.LT.RA(J))THEN
            RA(I)=RA(J)
            IB(I)=IB(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        RA(I)=RRA
        IB(I)=IRB
      GO TO 10
      END

