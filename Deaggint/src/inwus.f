c Is point (QLN,QLT) inside the wus
c Written by C. Mueller, USGS.
      parameter (NZ=52)
      dimension ZLN(NZ),ZLT(NZ)
      logical LXYIN, INWUS
	character*10 arg
	zln=(/
     + -126.,
     + -113.70,
     + -113.35,
     + -112.20,
     + -110.90,
     + -110.45,
     + -110.00,
     + -110.10,
     + -109.90,
     + -109.85,
     + -110.28,
     + -110.72,
     + -110.70,
     + -111.00,
     + -111.50,
     + -111.60,
     + -111.40,
     + -110.90,
     + -110.90,
     + -111.50,
     + -111.35,
     + -111.00,
     + -110.25,
     + -109.70,
     + -109.30,
     + -109.00,
     + -108.45,
     + -107.60,
     + -107.15,
     + -106.90,
     + -106.85,
     + -107.00,
     + -107.55,
     + -107.75,
     + -107.80,
     + -107.50,
     + -106.25,
     + -105.90,
     + -105.55,
     + -105.20,
     + -104.90,
     + -104.85,
     + -105.00,
     + -105.20,
     + -105.55,
     + -105.55,
     + -105.35,
     + -104.50,
     + -102.10,
     + -101.00,
     + -101.,
     + -126./)
	zlt=(/
     + 50.,
     + 50.00,
     + 49.00,
     + 47.55,
     + 46.45,
     + 45.75,
     + 45.50,
     + 44.90,
     + 44.55,
     + 44.15,
     + 43.17,
     + 42.18,
     + 39.55,
     + 39.40,
     + 38.50,
     + 37.95,
     + 37.45,
     + 37.05,
     + 36.80,
     + 36.20,
     + 35.70,
     + 35.25,
     + 34.80,
     + 34.55,
     + 34.55,
     + 34.70,
     + 35.00,
     + 35.50,
     + 35.90,
     + 36.25,
     + 36.65,
     + 37.00,
     + 37.30,
     + 37.55,
     + 37.85,
     + 38.40,
     + 38.55,
     + 38.30,
     + 38.05,
     + 37.80,
     + 37.45,
     + 37.00,
     + 36.00,
     + 35.00,
     + 34.35,
     + 33.85,
     + 33.25,
     + 31.70,
     + 29.60,
     + 29.20,
     + 24.6,
     + 24.6/)
	if(iargc().ge.2)then
	call getarg(1,arg)
	read(arg,'(f9.4)')QLN
	tmp=QLN
	call getarg(2,arg)
	read(arg,'(f8.4)')QLT
	if(tmp .gt.10.)then
	QLN=QLT
	QLT=tmp
	endif
	else
	print 5,'Enter site long and lat, dec deg: '
5	format(a,$)
	read *,QLN,QLT
	endif
      INWUS= LXYIN(QLN,QLT,ZLN,ZLT,NZ)
	print *,INWUS
      end

      logical function LXYIN (X,Y,PX,PY,N)
c Is point (X,Y) inside the (unclosed) polygon (PX,PY)?
c See "Application of the winding-number algorithm...",
c  by Godkin and Pulli, BSSA, pp. 1845-1848, 1984.
c LXYIN= .true. if (X,Y) is inside or on polygon, .false. otherwise.
c Written by C. Mueller, USGS.
      dimension PX(N),PY(N)
      LXYIN= .true.
      KSUM= 0
      do 1 I=1,N-1
        K= KPSCR(PX(I)-X,PY(I)-Y,PX(I+1)-X,PY(I+1)-Y)
        if (K.eq.4) return
        KSUM= KSUM+K
1       continue
      K= KPSCR(PX(N)-X,PY(N)-Y,PX(1)-X,PY(1)-Y)
      if (K.eq.4) return
      KSUM= KSUM+K
      if (KSUM.eq.0) LXYIN= .false.
      return
      end

      integer function KPSCR (X1,Y1,X2,Y2)
c Compute the signed crossing number of the segment from (X1,Y1) to (X2,Y2).
c See "Application of the winding-number algorithm...",
c  by Godkin and Pulli, BSSA, pp. 1845-1848, 1984.
c KPSCR= +4 if segment passes through the origin
c        +2 if segment crosses -x axis from below
c        +1 if segment ends on -x axis from below or starts up from -x axis
c         0 if no crossing
c        -1 if segment ends on -x axis from above or starts down from -x axis
c        -2 if segment crosses -x axis from above
c Written by C. Mueller, USGS.
      KPSCR= 0
      if (Y1*Y2.gt.0.) return
      if (X1*Y2.eq.X2*Y1.and.X1*X2.le.0.) then
        KPSCR= 4
      elseif (Y1*Y2.lt.0.) then
        if (Y1.lt.0..and.X1*Y2.lt.X2*Y1) KPSCR= +2
        if (Y1.gt.0..and.X1*Y2.gt.X2*Y1) KPSCR= -2
      elseif (Y1.eq.0..and.X1.lt.0.) then
        if (Y2.lt.0.) KPSCR= -1
        if (Y2.gt.0.) KPSCR= +1
      elseif (Y2.eq.0..and.X2.lt.0.) then
        if (Y1.lt.0.) KPSCR= +1
        if (Y1.gt.0.) KPSCR= -1
      endif
      return
      end

