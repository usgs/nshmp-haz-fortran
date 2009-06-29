c rescale.f. Compile with f77.
c    f77 rescale.f -o ~harmsen/bin/rescale
c or use the -f77 flag.
c mod dec 3,1999: Writes gmset.jpid, where pid makes filenames unique (pid=arg7)
c rescale returns a new PSA if the data say iterate again,  otherwise old PSA
c mod aug 7 2008: limit reduction to 40% of orig guess.
	real guess(12),rate(12)
c guess(j) will contain the jth PSA estimate and rate(j) will contain the
c corresponding (input) annual rate of exceedances. This set may be used to
c compute the next PSA estimate. 
c two or three updates should suffice.
c The update formula should be delta(PSA) = (delta rate)**(-1/S), where S is the
c slope of the hazard in log-log (e.g., McGuire, 1997 SSHAC report, App. G)
c Stephen Harmsen, US Geological Survey, (303) 273-8567
	character*12 arg
	character*72 loc
c arg(8) will be the location where file io goes (absolute path)
	character*1 n
c n is a PSA index (tested: 1 to 6) treated as a char*1, for file-name building
	call getarg(1,arg)
	read(arg,'(i1)')i01
	call getarg(2,arg)
	read(arg,'(i2)')it
	call getarg(3,arg)
	read(arg,'(f10.8)')ann
	call getarg(4,arg)
	read(arg,'(a1)')n
c n= variable number
	call getarg(5,arg)
	read(arg,'(f10.7)')psa
	if(i01.eq.1)then
c Flag says, "revise the psa estimate, based on available information:
c available information = set of previous estimates and their computed rates"
1	format(i2)
c The below "open" works if at first there is no gmset.n. Scripts should 
c sweep such files into the trash before running "rescale" 
	call getarg(7,arg)
	call getarg(8,loc)
	i=index(loc,' ')-1
	if(it.eq.1)then
	open(1,file=loc(1:i)//'gmset.'//n//arg,status='unknown')
	call getarg(9,arg)
	read(arg,'(f4.2)')SI
	else
	open(1,file=
     1loc(1:i)//'gmset.'//n//arg,access='append',status='old')
	endif
        guess(it)=psa
        call getarg(6,arg)
        read(arg,'(e10.5)')rate(it)
        write(1,44)guess(it),rate(it)
44	format(1x,f10.7,1x,e10.5)
	rewind(1)
	do i=1,it-1
	read(1,44)guess(i),rate(i)
	enddo
	close(1)
c File 1 stores the data. Contents are updated at each iteration.
c The above file manipulations may be unstable. Works on DecAlpha, Solaris systems OK.
	if(it.eq.1)then
c There is experience-based information to use: use default shape of haz curve.
c however limit reduction to some extent
	fac=max(0.4,(rate(1)/ANN)**SI)
	psa=psa*fac
	elseif(rate(it).ne.ANN.and.rate(it).ne.rate(it-1))then
	SI = alog(guess(it)/guess(it-1))/alog(rate(it-1)/rate(it))
c SI is now 1.0/(computed slope of haz curve) in appropriate region
	psa=guess(it)*(rate(it)/ANN)**SI
	endif
	endif
        write(6,60)psa
60	format(f10.7)
c	stop' '
c2000	stop'Rescale: please delete gmset.* before running me.'
	end
