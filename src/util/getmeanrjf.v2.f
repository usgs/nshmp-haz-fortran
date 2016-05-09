c getmeanrjf.v2.f: 
c revised feb 7, 2011 to include Geomatrix RL(M). Used for Samoa background.
c revised July 10, 2015 : l and l0 should be real variables. They were implicit integer
c in earlier versions. Effect on results was visible but quite small.
c Compute mean distance from fault with random strike measured from its midpoint
c Use Wells and C. (1) SRL or (2) Area for WUS, (3) Somerville A for CEUS or
c (4) Geomatrix subduction.
c V2 has Mmax of 8.55 and Nmag=26. Original version has Mmaxo of 7.55 and Nmag=16.
c V2 output requires a change to hazgridXnga5 to read in the larger array rjbmean.
c
c WC Area:
c West:		A=10**(0.91M-3.49)   km**2,  Zmax=15
c East: 	A= 10**(M-4.366)		Zmax=20 km
c Subduction:	RL= 10**((xmag-4.94)/1.39)  zmax=50 km
c	For M9, expected length is 833 km
c Solaris compile:
c
c	f95 getmeanrjf.v2.f -o getmeanrjf.v2 -e
c S Harmsen.
c
	logical ishere
	character*18 fn(4)/'rjbmean.bin.srl','rjbmean.bin.Awest','rjbmean.bin.Aeast',
     + 'rjbmean.bin.geoma'/
	character*18 fa(4)/'rjbmean.dat.srl','rjbmean.dat.Awest','rjbmean.dat.Aeast',
     + 'rjbmean.dat.geoma'/
        common/rjbmean/rjbmean(26,1001)
        real rjbb(26,1001)
	real mmin,mmax,dm
	mmin=6.05
	mmax=8.55
	dm=0.1
	rmax=1000.
	dr=1.0
1	write (6,50)'Enter 1 for WC SRL, 2 for WC WUS A, 3 for Somerville E',
     + ' or 4 Geomatrix subduction: '
	read *,itype
50	format(a,$)
	if(itype.lt.1.or.itype.gt.4)goto 1
	open(2,file=fa(itype),status='unknown')
	inquire(file=fn(itype),exist=ishere)
	open(3,file=fn(itype),form='unformatted')
	if(ishere)then
	read(3)rjbb
	endif
	nr =min(nint(Rmax/dr)+1,1001)
		nmag = nint((mmax-mmin)/dm)+1
	print *,nmag,' number of magnitude indexs'
	print *,nr,' number of distance indexs'
	call getmeanrjb(mmin,mmax,dm,rmax,dr,itype)
	xmag=mmin
	if(.not.ishere)rjbb=rjbmean
	do im=1,nmag
	r=0.001
	write(2,10)xmag
10	format(/,'#Mag ',f6.2,/,'#Dcenter(km) Dbar(km) from getmeanrjf')
	do ir=1,nr
	write(2,12)r,	rjbmean(im,ir)
	if(abs(rjbmean(im,ir)-rjbb(im,ir)).gt.1.e-4)then
	print *,rjbmean(im,ir)-rjbb(im,ir),' diff at im,ir',im,ir
	endif
12	format(f8.2,1x,f8.2)
	r=r+dr
	enddo
	xmag=xmag+dm
	enddo
	if(.not.ishere)write(3)rjbmean
	end
	
	subroutine getmeanrjb(mmin,mmax,dm,rmax,dr,itype)
c mmin = minimum M (moment), mmax= maximum M
c dm = discretization of M from mmin to mmax (typical 0.05 or 0.1)
c rmax = maximum epicentral distance (km) typical 200 to 1000 km
c dr = distance discretization typically 1 km
c itype one of 3 ways to compute lengt from M
c first dim: variation wrt Magnitude. flt length based on magnitude
c 2nd dim: variation wrt distance to center of fault
c from program getdist.f
c see Gradshetyn and Ryzhik, eqn. 2 at bottom of p 156. Then see Num. Recipes
c in Fortran, 1st edition, p 185
c Version 0.0 (orig)
c just computed mean distance when epicenter at midpt of fault. 
c You may add a
c do loop to integrate over portion of fault nearest the site & computes sd.
c Could write a file called mean.var.tmp that lists the mean distance and std. deviation
c as a function of l, where l varies from 0 to the input fault length.
c This file also reports theta, the angle where transition occurs from measuring
c distance from fault endpoint to dropping a perp. onto the fault.
c
c ev2 known bugs: 
c Some combinations may fail in ev2: for example, r=5 l=10. 
c Close neighbors work, e.g., 5.01 and 10, or 5 and 10.01. For most apps. a
c small dither of r or l may be acceptable.
c 2nd bad case: x and 4x for x reasonably small. 
c 3rd bad: r = l and r=6, 60, 8, 4.4, 44, 8.8, 88
c 4th bad case: r=100 l=101. Pattern aint exactly clear.
c Failure is a rare occurrence but needs to be patched for production work.
c Some of above potholes have been patched. Just need time to discover others.
c Faults should have positive length. Here l is forced to be >=0.001
c distance units are km or degrees depending on context. 
c Steve Harmsen, Apr 2008. Email harmsen@usgs.gov
	parameter(pi=3.141592654)
        common/rjbmean/rjbmean(26,1001)

	real mmin,mmax,dm,rmax,l,l0
c a reference point could be an epicenter or could be some arbitrary posn
c on a fault. 
	rjbmean=rmax+10.	!initialize all elements to a large distance
c	open(24,file='mean.var.tmp',status='unknown')
100	format('#l0  meandist std.dev. theta1, r= ',f6.1,' from getdist.f')	
	nok=0	!store number of ok estimates 
	nbad=0	!number of bad estimates (0 length fault is a bad case)
	nmag = nint((mmax-mmin)/dm)+1
	rad=pi/180.
	nr =min(nint(Rmax/dr)+1,1001)
	rmin=0.001
	r=rmin
	do ir=1,nr
	xmag=mmin
c mag loop
	do  m=1,nmag
	if(itype.eq.1)then
         arg= -3.22+0.69*xmag
	arg=10**arg
	elseif(itype.eq.2)then
	a=0.91*xmag -3.49
	a=10**a
	w=sqrt(a*0.5)
	w=min(15.,w)
	arg=a/w
	elseif(itype.eq.3)then
	a=xmag-4.366
	a=10**a
c golden ratio.
	w=min(20.,sqrt(a/1.61803))
	arg=a/w
	else
	arg = 10**((xmag-4.94)/1.39)
c this is Geomatrix subduction relation. 159 km for M=8.
	endif
c for the fault very close to the site, and gridded, set r to zero.
	if(r.lt.0.5)then
	rjbmean(m,ir)=0.0
	else	
	rsq=r*r
         l = arg
         l0=0.5*l	!Important: epicenter assumed at midpoint
	if(l0.ge.r)then
c the easy problem.
c expc = first moment = 2/pi Integral r sin(theta) dtheta
	expc = 2.*r/pi
c esq = second moment = 2/pi Integral r**2 sin**2(theta) dtheta 
	esq= 0.5*r**2
	else
c elliptic integral of the 2nd kind, use Num. Recipes		
	a=r**2+l0**2
	b=2.*r*l0
	costheta=l0/r
c 	
	theta=acos(costheta)
c	theta is the transition angle where distance ceases to be measured
c from fault endpoint. Beyond theta distance measured along the perp.
c	pi4=0.25*pi
	delta=sqrt((a+b)*(1-costheta)/2./(a-b*costheta))
	delta= min(0.999997,delta)
	delta=asin(delta)
	x=tan(delta)
	r2=2.*b/(a+b)
c s2 below is k_c**2 in eqn 6.7.8 of Num Recipes v1
	s2=1.0-r2
	if(s2.le.0.)then
c the protection here is needed because cant pass 0 to el2
c	write(6,*)x,s,' bad case. patch applied'
	s=1.e-6
	s2=1.e-12
	else
	s=sqrt(s2)
	endif
c In this instance, the function el2 returns a number
c related to the integral sqrt(a-bcos(theta)),
c theta from 0 to theta, for the above a b and theta (functions of r&l).
c The solution was checked by simple numerical integration commented out below.
c The second moment is an easy trig integral.
	e2= el2(x,s,1.,s2)
        fac1=2.*sqrt(a+b)
	term2=2.*b*sin(theta)/sqrt(a-b*costheta)
c	write(6,*)'fac1, e2, term2, a, b ',fac1,e2,term2,a,b
	ans=fac1*e2-term2
c	write(6,*)ans,' should be about the same as numerical soln below'
c try numerical integration. ( integrate only up to theta. not to pi/2)
c	dth=0.002*theta
c	th0=0.5*dth
c	sum=0.0
c	do i=1,500
c	sum=sum+dth*sqrt(a-b*cos(th0))
c	th0=th0+dth
c	enddo
c	write(6,*)' numerical sum and final theta ',sum,th0-dth
c exp below is the expected value of distance, r_jb
	expc=2.0/pi*(ans+r*costheta)
c
c esq=2nd moment should just be integral of (a-b*cos theta) from 0 to theta
c	plus r**2 sin**2(theta) effect beyond theta
c
	esq1=a*theta - b*sin(theta)
	esq2=(0.5*pi-theta+0.5*sin(2.0*theta))*rsq*0.5
	esq=2./pi*(esq1+esq2)
	endif
C	write(6,*)' Mean distance of randomly striking fault to site ',expc
	if(expc.ge.0.)then
	rjbmean(m,ir)=expc
c	rjbvar(m,ir)=esq-expc**2
	else
	nbad = nbad+1
	rjbmean(m,ir)=Rmax
c	rjbvar(m,ir)=10.	!dont know. check that this doesnt happen
	endif
	endif	!r < rmax
	xmag=xmag+dm
	enddo	!mag
	r=r+dr
	enddo	!distance index ir
	if(nbad.gt.0)write(6,*)'Houston, we have ',nbad,' problems'
24	format(f8.2,1x,f8.2,1x,i3,1x,i3,3(1x,f10.4))
6	format(a)
	return
	end subroutine getmeanrjb

      FUNCTION EL2(X,QQC,AA,BB)
      PARAMETER(PI=3.14159265, CA=.0003, CB=1.E-9)
      IF(X.EQ.0.)THEN
        EL2=0.
      ELSE IF(QQC.NE.0.)THEN
        QC=QQC
        A=AA
        B=BB
        C=X**2
        D=1.+C
        P=SQRT((1.+QC**2*C)/D)
        D=X/D
        C=D/(2.*P)
        Z=A-B
        EYE=A
        A=0.5*(B+A)
        Y=ABS(1./X)
        F=0.
        L=0
        EM=1.
        QC=ABS(QC)
1       B=EYE*QC+B
        E=EM*QC
        G=E/P
        D=F*G+D
        F=C
        EYE=A
        P=G+P
        C=0.5*(D/P+C)
        G=EM
        EM=QC+EM
        A=0.5*(B/EM+A)
        Y=-E/Y+Y
        IF(Y.EQ.0.)Y=SQRT(E)*CB
        IF(ABS(G-QC).GT.CA*G)THEN
          QC=SQRT(E)*2.
          L=L+L
          IF(Y.LT.0.)L=L+1
          GO TO 1
        ENDIF
        IF(Y.LT.0.)L=L+1
        E=(ATAN(EM/Y)+PI*L)*A/EM
        IF(X.LT.0.)E=-E
        EL2=E+C*Z
      ELSE
c        PAUSE 'failure in EL2'
	write(6,*)'failure in el2'
      ENDIF
      RETURN
      END function el2

