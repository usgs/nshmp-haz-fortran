c finds rate of eqs in a given M range.
c Input 10**a and b (slope of distribution or density function)
c Uses Herrmann formula to get cumu. rate, then brute forces the same calc
c in a loop.
c
	real mmin,mmax,m
	real fac,ac,a,rate_c,ratefac
	logical modrate/.false./
	character*1 ch1
c cumulative rate of eqs > M. Same form as interval rate...
	rate_c(a,b,xm)=a/10.**(b*xm)
1	continue
	print 5,'Enter 10**a '
	read *,a
	print 5,'Enter b '
5	format(a,$)
	read *,b
	afac=a
	print 5,'mmin,mmax,dm: '
	read *,mmin,mmax,dm
	print 5,'Modify rate factors for some intervals 1=yes 0=no: '
	read *,i
	if(i.eq.1)then
	modrate=.true.
	else
	ratefac=1.
	endif
	dm2=dm*0.5
	m=mmin
	im=nint((mmax-mmin)/dm) +1
	fac=10.**(b*dm2)
c ac = avalue for cumu distrib., from Herrmann eq notes jan-june 1977, eqn 5.
c  Herrmann notes use 2*dM as the interval,
c but a-value and our codes use 1*dM, i.e., M0+-dM/2
	ac=a/(fac-1./fac)
c	print *,fac,ac, ' term in den. and ac'
	print *,'Cumulative rate between Mmin and Mmax ',rate_c(ac,b,mmin-dm2)-rate_c(ac,b,mmax+dm2)
	ratelow=rate_c(ac,b,mmin-dm2)-rate_c(ac,b,mmin+dm2)
	print *,'Rate between Mmin-dM/2 and Mmin + dM/2',ratelow
	print 5,'Suppose you want to change sampling interval dm. Enter a new dm: '
	read *,dmn
	dmn=dmn/2.
	print *,'Then 10**a would be ',rate_c(ac,b,-dmn)-rate_c(ac,b,dmn)
	bfac=0.0
	print 24
24	format('Magnitude rate(Mmin <=M <=that Magnitude)')
	x=10**(-b*dm)
c untruncated
	crate=ratelow/(1.0-x)
	do i=1,im
	pow=10**(-b*m)
	bfac=bfac+pow
c	print *,bfac,' sum bfactors'
	if(mod(i,3).eq.1)print 25,m,afac*bfac
25	format(f5.2,1x,f11.8)
	if(modrate)then
	print 28,'enter factor for mag ',m
	read *,ratefac
	endif
28	format(a,f5.2,1x,': ',$)
	write(2,*)m,pow,pow*ratefac
	m=m+dm
	enddo
	rate = afac*bfac
	print 33,'Cumulative rate of eqs in ',mmin,mmax,' range: ',rate
33	format(a,f5.2,1x,f5.2,a,e10.5)
	print 35,'Mean recurrence time ',1/rate
35	format(a,f9.1,' years')
	print *,'Untruncated cumulative GR rate is ',crate
	print 35,'Mean recurrence time ',1/crate
	print 5,'to continue enter 1: '
	read 50,ch1
50	format(a)
	if(ch1.eq.'1')goto 1
	end
