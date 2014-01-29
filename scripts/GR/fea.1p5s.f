c interpolate 1.5s coeffs for Frankel and Atkinson tables CEUS. sh jan 17 2014.
c
c compile: f95 fea.1p5s.f -o fea.1p5s -e
c 
      character*30 hd
           real,dimension(20):: gma,gmb,gmc
c Arock frankel
	open(1,file='t1p0k006.tbl',status='old')
	open(2,file='t2p0k006.tbl',status='old')
	open(3,file='t1p5k006.tbl',status='unknown')
	read(1,5)hd
	read(2,5)hd
	write(3,5)hd
5	format(a)	!fix 
	frac=alog(1.5)/alog(2.0);frock=alog10(1.5)/alog10(2.0)
	fric=1.-frac
	print *,'fric frac',fric,frac,frock
	do i=1,21
	read(1,*)d,gma
	read(2,*)e,gmb
	if(d.eq.e)then
	gmc=fric*gma + frac*gmb
	else
	stop'mismatch distances'
	endif
	write(3,30)d,gmc
30	format(f6.2,(20(1x,f5.2)))
	enddo
	close(1);close(2);close(3)
c BC frankel
	open(1,file='t1p0k01l.tbl',status='old')
	open(2,file='t2p0k01l.tbl',status='old')
	open(3,file='t1p5k01l.tbl',status='unknown')
	read(1,5)hd
	read(2,5)hd
	write(3,5)hd
	frac=alog(1.5)/alog(2.0);frock=alog10(1.5)/alog10(2.0)
	fric=1.-frac
	print *,'fric frac',fric,frac,frock
	do i=1,21
	read(1,*)d,gma
	read(2,*)e,gmb
	if(d.eq.e)then
	gmc=fric*gma + frac*gmb
	write(3,30)d,gmc
	else
	stop'mismatch distances frankel bc'
	endif
	enddo
	close(1);close(2);close(3)
c BC Atkinson
	open(1,file='Abbc1p00.tbl',status='old')
	open(2,file='Abbc2p00.tbl',status='old')
	open(3,file='Abbc1p50.tbl',status='unknown')
	read(1,5)hd
	read(2,5)hd
	write(3,5)hd
	frac=alog(1.5)/alog(2.0);frock=alog10(1.5)/alog10(2.0)
	fric=1.-frac
	print *,'fric frac',fric,frac,frock
	do i=1,21
	read(1,*)d,gma
	read(2,*)e,gmb
	if(d.eq.e)then
	gmc=fric*gma + frac*gmb
	write(3,30)d,gmc
	else
	stop'mismatch distances atkinson bc'
	endif
	enddo
	close(1);close(2);close(3)
c Arock Atkinson
	open(1,file='ABHR1P00.TBL',status='old')
	open(2,file='ABHR2P00.TBL',status='old')
	open(3,file='ABHR1P50.TBL',status='unknown')
	read(1,5)hd
	read(2,5)hd
	write(3,5)hd
	frac=alog(1.5)/alog(2.0);frock=alog10(1.5)/alog10(2.0)
	fric=1.-frac
	print *,'fric frac',fric,frac,frock
	do i=1,21
	read(1,*)d,gma
	read(2,*)e,gmb
	if(d.eq.e)then
	gmc=fric*gma + frac*gmb
	write(3,30)d,gmc
	else
	stop'mismatch distances atkinson a'
	endif
	enddo
	close(1);close(2);close(3)
	print *,'Hand edit the header line 1 for period 1.5 s'
			end
			