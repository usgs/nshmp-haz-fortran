c Check arrays craton and margin to determine if point is in which one.
c  Gfortran compile:
c		gfortran craton_ck.f -o craton_ck -static
      logical, dimension(128000):: craton,margin
c Seems to be valid when input specified at tenths of a degree. Dubious elsewhere.
c
        open(2,file='../CU/GR/margin',form='unformatted')
        open(1,file='../CU/GR/craton',form='unformatted')
        read(2)margin
        read(1)craton
	print *,'Enter location to nearest tenth of a degree.'
1       continue
5       format(a,$)
        print 5,'Enter CEUS lat and lon to check if craton: '
        read *,ylat,xlon
        indx = (50.-ylat)*5010+(xlon+115.)*10+1
        print *,'index for this loc is ',indx
	if(indx.lt.1)goto 1
        print *,'Margin? ',margin(indx)
        print *,'Craton? ',craton(indx)
        print 5,'Enter a 1 to continue: '
        read *,i
        if(i.eq.1)goto 1
        end
