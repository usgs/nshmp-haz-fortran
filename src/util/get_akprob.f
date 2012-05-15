c get_akprob.f
c convert binary rate files to ascii. fast retrieve of data
c files like eqrate.ak.50km.bin 5 fields/record work for input.
c from get_prob.f
	logical seven/.false./,six5/.false./,six/.false./
	real xmin,xmax,ymin,ymax,xmag,t
	character*80 inf,outf
	character*80 hd
	character*12 ax
	real rec(5)
	call getarg(1,inf)
	call getarg(2,outf)
	open(1,file=inf,status='old',form='unformatted')
	open(2,file=outf,status='unknown')
c	open(3,file='temp.diag',status='unknown')
	nrec=0
	call getarg(3,ax)
	read(ax,'(f8.3)')xmin
	call getarg(4,ax)
	read(ax,'(f8.3)')xmax
	call getarg(5,ax)
	read(ax,'(f8.3)')ymin
	call getarg(6,ax)
	read(ax,'(f8.3)')ymax
50	format(a)
C	Write(3,23)ymin,ymax
C23	format('ymin ymax in get+prob ',2f6.2)
	read(1)rec
	read(1)rec
	call getarg(7,ax)
	i=index(ax,'.')
	j=index(ax,' ')
	if(i.eq.0)ax=ax(1:j-1)//'.0'
	read(ax,'(f8.3)')xmag
	call getarg(8,ax)
c t=exposure time yrs.
	i=index(ax,'.')
	j=index(ax,' ')
	if(i.eq.0)ax=ax(1:j-1)//'.0'
	read(ax,'(f8.3)')t
	if(xmag.eq.7.0)then
	seven=.true.
	elseif(xmag.eq.6.5)then
	six5=.true.
	write(3,50) 'alaska not ready for 6.5 cutoff'
	goto 12
	elseif(xmag.eq.6.0)then
	six=.true.
	endif
	do i=1,1000000
	read(1,end=12)rec
	if(rec(1).ge.xmin-.05.and.rec(1).le.xmax+.05.and.
	1rec(2).ge.ymin.and.rec(2).le.ymax)then
	if(seven)then
	p=1.-exp(-rec(5)*t)
	write(2,10)rec(1),rec(2),p
10	format(f8.3,1x,f6.3,1x,e11.5)
c	elseif(six5)then
c	p=1.-exp((-rec(5)-rec(6))*t)
c	write(2,10)rec(1),rec(2),p
c	if(rec(1).gt.-100.)print *,rec(1),rec(2),p
	elseif(six)then
	p=1.-exp((-rec(4)-rec(5))*t)
	write(2,10)rec(1),rec(2),p
	else
	p=1.-exp((-rec(3)-rec(4)-rec(5))*t)
        write(2,10)rec(1),rec(2),p
	endif
c go inside the map somewhat to establish max prob
        if(rec(1).ge.xmin+.06.and.rec(1).le.xmax-.06.and.
	1rec(2).ge.ymin+.1.and.rec(2).le.ymax-.1)
	2pmax=max(p,pmax)
	elseif(rec(2).lt.ymin)then
c wus usually safe to stop. at crease must keep looking at records
	goto 12
	endif
	enddo
12	close(1)
	close(2)
c	write(3,*)pmax
c	close(3)
	if(pmax.gt.0.79)then
	outf= 'prob.map.cpt'
	elseif(pmax.gt.0.61)then
	outf = 'prob.0top8.cpt'
	elseif(pmax.gt.0.35)then
	outf= 'prob.0top5.cpt'
	elseif(pmax.gt.0.21)then
	outf= 'prob.0top3.cpt'
	elseif(pmax.gt.0.115)then
	outf= 'prob.0top2.cpt'
	else
	outf= 'prob.0top1.cpt'
	endif
	print 13,outf
13	format(a)
	end
