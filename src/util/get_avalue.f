c get_avalue gets the avalue based on conservation of moment rate
c and on conservation of event rate
c This version assumes rates are accumulated in an interval (M, M+dM], i.e., include
c upper end M but exclude lower end M. This can make a difference, compared to using [M,M+dm).
c Steve Harmsen
c SOLARIS compile:
c f95 get_avalue.f -e -o get_avalue
c 
c Or alternate, cumulative rate and use Herrmann's formula to get incremental
c Input a specified M, dM, rate. And a finer scale dM for output. Typically dMout=k dMin,
c where k is 5 or 10.
	parameter (pi=3.1415927)
	character*20 name1
	logical cumul
	real*8 emom,sum,bfac,a,a0,xsum,sm,dM
	xm(x)=10.**(x*1.5+16.)
c Th&Hanks
	r(f)=2.34*3.5/2./pi/f
	const=2.34*3.5*10**5/2./pi*(16./7.)**0.33333
 	print 5,'Enter b value: '
 	read *,b
 	print 5,'Enter final dM [0.1]: '
 	read *,dM
 	print 5,'Enter a magnitude and an original dM [0.5] > final dM: '
 	read *,xml,dMo
5	format(a,$)
      xmo= 1.5*xml + 9.05
      xmo= 10.**xmo
	print *,'Standard moment, n-m units: ',xmo
	print 5,'Enter a rate for that M: '
	read *,rate
	print 5,'Is that rate incremental (0) or cumulative (1): '
	read *,i
	cumul = i.ge.1
	if(cumul)then
	  print *,'The original dM will not be used in this case.'
	      dmag2= dM/2.
	      bval=b
     	   fac= 10.**(bval*dmag2)-10.**(-bval*dmag2)
	   rate=rate*fac
	   print *,rate,' is the incremental rate from Herrmann conversion'
	   print *,fac ,' is the factor to get there'
	   a=alog10(rate)+xml*bval
	   print *,'The value to enter into hazSUBXnga is ',a
	   a0=10**a
	   print *,'The rate of M0 is ',a0
	   else
	xmorate=rate*xmo
	print *,'Moment rate is ',xmorate
c starting Mag = sm
	n= nint(dMo/dM)-1
	sm= xml-dMo/2.+dm;sm0=sm
	print *,'starting mag ',sm
	bfac=10**(-b*sm)
	sum=0
	em = sm+ dM*n
	print *,'ending mag ',em
	do i=0,n
	print *,'mag and bfac now ',sm,bfac
	emom=10**(1.5*sm+9.05)
	print *,emom
	sum=sum+bfac*emom
	sm=sm+dM
	bfac=10**(-b*sm)
	enddo
	a=xmorate/sum
	a0=a
	a= dlog10(a)
	print *,'A  value based on conservation of moment rate:'
	print *,a,sum*a0
	sm=sm0
	xsum=0
	do i=0,n
	x=a0*10**(-b*sm)
	xsum=xsum+x
	print *,sm,x,xsum
	sm=sm+dM
	enddo
c next rate conservationn
	print *,'next a corresponds to rate-of-events conservation'
	sm=sm0
	sum=0.
	do i=0,n
	bfac=10**(-b*sm)
	print *,'mag and bfac now ',sm,bfac
	sum=sum+bfac
	sm=sm+dM
	enddo
	a=rate/sum
	a0=a
	a= dlog10(a)
	print *,'  **************WOW I have an estimate*************'
	print *,a,' a from event-rate conserve'
	endif
	ratex=0.
	sm=5.0
	do i=0,35

	ratex=ratex+a0*10**(-b*sm)
	print 20,sm,ratex
20	format(f6.2,1x,1pe11.5,1x,'M , accumulating rate')
	sm=sm+dM
	enddo
	
	end
	
