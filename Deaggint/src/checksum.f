c checksum.f
c To compile: f95 checksum.f -o checksum
c 2008 version: toss out records that begin with '#' character. can 
c occur just about anywhere inside the file. comment lines.
c Use: checksum path control.file
c the path tells where to write some temp files to
c the control file tells which files to sum and what weights, w(i), to apply.
c as of oct 3,2001 both input and output files should be in a "log dir" which
c is still read in as arg(1)
c sums w(i)*haz(i,j), i=1,nfil, j=1
c this sum is supposed to approximate a specified annual rate
c Input records have a distance, d, and magnitude, M, preceeding the hazard vector.
c
c Output: computed mean annual exc_rate for the set of files input & their weights
c Steve Harmsen 8/4/2008
c
	parameter (nmax=10000)
c nmax=max number of records per file. Be careful here!
	logical getwt
	real d,m
	character*80 loc
	character*40 name1,namein
	character*80 rec
	character*1 ch1
	integer iloc
	real h,sumhaz
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
	read (4,5)ch1
	getwt= (ch1.ne.'y') .and. (ch1.ne.'Y')
	wt=1./nfil
	read (4,5)ch1
	do i=1,nfil
1	read (4,5)name1
	if(getwt)read (4,*)wt
	open(1,file=loc(1:iloc)//name1,status='old',err=2000)
	read (1,5,end=20)rec
c there are files that are empty: they correspond to 0-hazard cases.
c Non-empty files are expected to have a header record.
	wtsum=wtsum+wt
	do j=1,nmax
	read(1,5,end=20)rec
	if(rec(1:1).ne.'#')then
	read(rec,*,err=20,end=20)d,m,h
	sumhaz=sumhaz+h*wt
	endif
	enddo
20	close(1)
2000	continue
        enddo
	write(6,26)sumhaz
26	format(e10.5)
	close(4)
	end
