c	program vancouver check for vancouver island and various other places
c or pacific nw
c or several USA states
	character*10 a
	character*1 c1(0:16)
	data c1/'0','1','2','3','4','5','6','7','8','9','C','T','M','O',
     +'o','N','S'/
	real lat,long
	call getarg(1,a)
	i=index(a,'.')
	j=index(a,' ')-1
	if(i.eq.0)a=a(1:j)//'.0'
	read(a,'(f6.3)')lat
	call getarg(2,a)
	i=index(a,'.')
	j=index(a,' ')-1
	if(i.eq.0)a=a(1:j)//'.0'
	read(a,'(f9.3)')long
	if(long.lt.-123.3.and.lat.gt.48.44)then
	i=1
	elseif(long.lt.-114.9.and.lat.gt.44.01)then
	i=2
c Pacific NW
	elseif(long.lt.-116.99.and.lat.gt.42.001)then
	i=14
c southern Oregon	
	elseif(long.lt.-114..and.long.gt.-120..and.lat.lt.42.
     +.and.lat.gt.38.8)then
	i=15
c northern nevada	
	elseif(long.lt.-113.99.and.lat.lt.35.03)then
	i=3
c southern Calif
	elseif(lat.lt.37.4.and.long.lt.-115.3-9./7.*(lat-35.5)) then
	i=4
	elseif(lat.le.42..and.long.lt.-120.)then
	i=5
	elseif(long.gt.-109.01.and.long.lt.-101..and.lat.gt.37.0
     +.and.lat.lt.41.01)then
c colorado
	i=10
	elseif(lat.gt.35.95.and.long.gt.-116..and.long.lt.-108.) then
	i=6
c isb
	elseif(lat.lt.43..and.long.gt.-120..and.long.lt.-113.5)then
	i=7
c Basin and Range +-
	elseif((lat.gt.45..and.long.lt.-71..and.long.gt.-77.).or.
     + (lat.gt.43.5.and.long.lt.-76.6.and.long.gt.-83.))then
	i=9
c southern Canada shouldnt do whatever
	elseif(abs(lat-33.5) .le.1.3 .and.abs(long+81.).lt.1.5)then
	i=16
	elseif((lat.lt.26.01.and.long.lt.-96.).or.(lat.lt.27..and.
     +long.lt.-100.))then
	i=12
c Mexico site	
	elseif(lat.lt.34.8.and.long.lt.-95.and.long.gt.-101.)then
c Texas
	i=11
	elseif(lat.lt.36.7.and.lat.gt.34.8.and.long.gt.-100.5.and.
     +long.lt.-95.)then
	i=13
c oklahoma	
	elseif(long.gt.-100.001)then
	i=8
	else
	i=0
	endif
	write(6,2)c1(i)
2	format(a)
	end
