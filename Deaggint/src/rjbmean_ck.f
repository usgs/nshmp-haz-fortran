      logical, dimension(128000):: craton,margin
c check the rjb mean values see if move to gldplone worked.
c
	open(2,file='../GR/margin',form='unformatted')
	open(1,file='../GR/craton',form='unformatted')
	read(2)margin
	read(1)craton
1	continue
5	format(a,$)
	print 5,'Enter CEUS lat and lon to check if craton: '
	read *,ylat,xlon
	indx = (50.-ylat)*501+(xlon-115.)*10+1
	print *,'Margin? ',margin(indx)
	print *,'Craton? ',craton(indx)
	print 5,'Enter a 1 to continue: '
	read *,i
	if(i.eq.1)goto 1
	end
