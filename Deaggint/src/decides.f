c decides.f : gldjeff compile f77 decides.f -Bstatic -o ../decides
c decides returns a 1 if the data say iterate again, Sam, otherwise 0
c The decision is based on the maximum error, computed over the set of PSAs
c being considered.  But, never more than itmax iterations.
c Maximum number of periods  of PSA set is 10
	parameter (itmax=3)
	logical doover
	real rate(10)
	character*12 arg
c        character*30 ident / "@(#)PSHA Deagg Cntrl Program"/
	iper=iargc()-3
	if(iper.gt.10)stop 'sorry 10 PSA values maximum'
	call getarg(1,arg)
	read(arg,'(i1)')it
	if(it.gt.itmax)then
	i=0
	else
	call getarg(2,arg)
	read(arg,'(f5.3)')tol
	call getarg(3,arg)
	read(arg,'(f10.8)')ann
	doover=.false.
	do i=1,iper
	call getarg(i+3,arg)
	read(arg,'(e10.5)')rate(i)
	doover=doover.or.(abs(1.-rate(i)/ann).gt.tol)
	enddo
	i=0
	if(doover)i=1
	endif
	write(6,1)i
1	format(i2)
	end
