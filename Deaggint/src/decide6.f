c decide6.f : gldjeff compile f77 decide6.f -Bstatic -o ../decide6
c gldplone, use gfortran:
c   gfortran decide6.f -static -o decide6 -ffixed-line-length-none
c
c Decide6 returns a 1 if the data say iterate again, Sam, otherwise 0.
c currently the decision is based on 2nd rate, that is, the one associated
c with the 2% in 50 year pe. A collective decision, based on worst case or
c something like it, might be considered.
c That decision would be based on the maximum error, computed over the set of PGAs
c being considered.  But, never more than itmax iterations.
c Number of elements (levels) of PGA set is 6 . SHarmsen nov 10 2008
	parameter (itmax=3)
	logical doover
	real rate(6),ann(6)
	integer i(6)
	character*12 arg
	ann = (/0.000201, .000404, .001026, .002107, .0044629, .01386 /)
	call getarg(1,arg)
	read(arg,'(i1)')it
	if(it.gt.itmax)then
	i=0
	else
	call getarg(2,arg)
	read(arg,'(f5.3)')tol
	i=0
c Sometimes there is a stray: ( symbol creeping in the stream at arg j+2.
c However, this problem seems to be intermittant. Why?
	do j=1,6
	call getarg(j+2,arg)
	read(arg,'(e10.5)')rate(j)
	doover=abs(1.-rate(j)/ann(j)).gt.tol
	if(doover)i(j)=1
	enddo
	endif	!3 iterations or less?
	write(6,1)i(2)
1	format(i1)
	end
