c checksum6.f
c To compile: f95 checksum6.f -o checksum6
c 2008 version: 2 header records and then pure data
c 6 ground motions : for a given record, which one is indexed as field 5
c Use: checksum6 path control.file
c the path tells where to write some temp files to
c the control file tells which files to sum and what weights, w(i), to apply.
c sums w(i)*haz(i,j), i=1,nfil, j=1,...,6
c The jth sum is supposed to approximate a specified annual rate
c Input records have a distance, d, and magnitude, M, preceeding the hazard scalar,
c and the j for each record follows h (see read statement).
c
c Output: 6 computed mean annual exc_rates for the set of files input & their weights
c Steve Harmsen 11/07/2008. cosmetic mods 11/10/2008.
c
	parameter (nmax=100000)
c nmax=max number of records per file. Be careful here!
	real d,e,h,m,wt
	real, dimension(6) :: sumhaz
	character*80 loc
	character*40 name1,namein*40
	character*80 rec
	character*1 ch1
	integer iloc,i,j,k,nfil
	sumhaz=0.0
c	argc=iargc()
	if(iargc().ge.2)then
	loc='               '
	call getarg(1,loc)
	iloc=index(loc,' ')-1
	call getarg(2,namein)
	else
	stop 'checksum usage: checksum loc input.control.file'
	endif
        goto 8
1999	stop'checksum: control file not found'
8	open(4,file=loc(1:iloc)//namein,status='old',err=1999)
5	format(a)
50	format(a,$)
	read (4,*)nfil
	sumhaz=0.0
	do i=1,nfil
1	read (4,5)name1
	read (4,*)wt
	open(1,file=loc(1:iloc)//name1,status='old',err=2000)
	read (1,5,end=20)rec
	read (1,5,end=20)rec
c there are files that are empty: they correspond to 0-hazard cases.
c Non-empty files are expected to have a header record.
	wtsum=wtsum+wt
	do j=1,nmax
	read(1,*,err=20,end=20)d,m,h,e,k
	sumhaz(k)=sumhaz(k)+h*wt
	enddo
20	close(1)
2000	continue
        enddo
	write(6,26)sumhaz
26	format(6(1x,e10.5))
	close(4)
	end
