c avg_dist find average distance of points of secondary flt to master flt.
c Motivation:
c hazFX needs a single fault width for all faults in the model. For secondary faults,
c this width is generally constrained by the master fault location and dip. 
c All points of master flt within latmax+0.05 degrees and latmin-0.05 degrees are considered
c where latmin is the minimum latitude of the secondary flt and latmax is the maximum
c lat of the secondary flt .
c Code run will write a file of characteristic events suitable for input to hazFXnga7c or similar
c It also writes a summary file called 'branches.out'.
c Assumed sense of slip on 2nd flt is normal or ss. Reverse is not considered here but
c could be with an additional instruction.
c
c compile instr solaris:
c f95 avg_dist.f -e -o ../bin/avg_dist
c compile instr gfortran:
c gfortran avg_dist.f -ffixed-line-length-none -o avg_dist
c
c Usage: avg_dist master.fault.file secondary.fault.file
c
c Each input file should consist of:
c line1: slip_rate(mm/yr) or event_rate(events/yr) preferred_dip(deg) dip_code(1 to 5) Name_of_flt (in quotes)
c line 2: an integer n 
c next lines: n records of x,y (lon,lat i.e., location of surface trace of the flt)
c n is the number of points in the fault trace.
c if slip_rate is greater than 0.005 it is assumed slip rate; otherwise assumed
c event-rate which will be fixed (new Apr 20 2012, to model eglington) 
c In distance calculations, master fault points are considered
c if they are reasonably close to secondary flt (0.05 d lat maximum difference).
c So 2nd flt can go past master to north or south and some of those points will be considered.
c Treatment when extending 2nd flt beyond master is to assume same separation as the last computed one.
c Although the 2nd fault width might jump down to base of seismogenic, we do not consider that. Instead, we
c assume master fault extends in the same direction even though it was not mapped. High level of uncert
c in this part of model.
c Steve Harmsen, Apr 12, 2012. harmsen@usgs.gov
	parameter (d2r = 0.0174533,pi=3.141592654)
	character(40) fil,name1,name2
	real, dimension(500):: solat,solon,twolat,twolon
	real maxd/0./,mind/100000./,maga
c dip1,wt1 are the dip and weight distribution to assign to fault 1 : Master
c dip2,wt2 are the dip and weight distribution to assign to fault 2 : Secondary
	real,dimension(5) :: ddip, dip1,dip2,wt1,wt2
	real wts(5,5)
	logical fixrate
c ddip is the range of dip models delta-dip (degrees) from the preferred value
   	amo(x)=10.0**(1.5*x+9.05)
	ddip = (/-15.,-10.,0.,10.,15./)
	wts(1:5,1) =(/0.05,0.15,0.6,0.15,0.05/)	!std symmetric
	wts(1:5,2) = (/0.1,0.25,0.5,0.1,0.05/)	!skew low
	wts(1:5,3) = (/0.05,0.1,0.5,0.25,0.1/)	!skew high
	wts(1:5,4) = (/0.1,0.2,0.4,0.2,0.1/)	!fat tails
	wts(1:5,5) = (/0.2,0.2,0.2,0.2,0.2/)	!uniform (know nothing)
	if(iargc().lt.2)stop'usage avg_dist master.flt second.flt'
        xmu= 3.0e10
	call getarg(1,fil)
	open(3,file='branches.out',status='unknown')
	write(3,5)'#avg dist.f find distances from master flt ',fil
	ymax=-90.
	ymin=90.
	open(1,file=fil,status='old',err=345)
	call getarg(2,fil)
	open(2,file=fil,status='old',err=345)
5	format(a,a30)
10	format(a)
50	format(a,$)
	print *,'Enter output file name for char input to hazFX (25 branches) : '
	read 10,name1
	open(4,file=name1,status='unknown')
	write(3,5)'#Secondary flt ',fil
	write(3,10)'#width(km)  bot_depth Area  Mag(Area) branch_wt dip1   dip2 Morate     RI(yrs)'
	print *,'Summary output in branches.out'
	print *,'Detailed output in ',name1
	read(2,*)slip2,dipp,ins2,name2
	if(slip2.lt.0.005)then
	fixrate=.true.
	print *,'A fixed event rate of ',slip2,' will be used in the output'
	rate=slip2
	else
	fixrate=.false.
	endif
	if(ins2.gt.5 .or. ins2.lt.1)stop'Dip uncertainty code must be in range 1 to 5 (2nd flt)'
c slip2 is vertical comp of slip for a dipping fault.
	wt2=wts(1:5,ins2)
	dip2=(dipp+ddip)*d2r
	read(2,*)ns2
	if(ns2.gt.500)stop'increase dim of twolat etc arrays'
	do i=1,ns2
	read(2,*)twolon(i),twolat(i)
	ymax=max(ymax,twolat(i)+0.05)
	ymin=min(ymin,twolat(i)-.05)
	enddo
c In this code fault length is based on distance from first to last point (like W&C)
	call delaz2(twolat(1),twolon(1),twolat(ns2),twolon(ns2),flen,azi)
	close(2)
	y1min=100.
	y1max=-100.0
	read(1,*)slip1,dipp,ins1,name1
	if(ins1.gt.5 .or. ins1.lt.1)stop'Dip uncertainty code must be in range 1 to 5 (master flt)'
	wt1=wts(1:5,ins1)
	dip1=(dipp+ddip)*d2r
	nomit=0
	read(1,*)ns1
	i1=0
	if(ns1.gt.500)stop'increase dim of solat etc arrays'
	do i=1,ns1
	read(1,*,err=456)x,y
	if(y.lt.ymax.and.y.gt.ymin)then
	i1=i1+1
	solat(i1)=y
	solon(i1)=x
	y1max=max(y1max,solat(i1))
	y1min=min(y1min,solat(i1))
	else
c	print *,'this point too far away from secondary flt',x,y
	nomit=nomit+1
	endif
	enddo
	if(i1.lt.2)stop'not enough overlap to consider antithetic'
	ns=i1
	avgdist=0.
	njump=0
	wid=5.
	ftop=0.
	arate=0.	!average event rate on 2nd flt.
	avgdip=0.0
	sumwt=0.0
	dipf=090.	!surface dist
	ymo=0.0		!avg moment rate from all branches
	disb=0.
	nc=0
	ndump=0
c next check verifies that first point projects onto master fault. (just an east-west projection)
c Should be more general.
	if(twolat(1).gt.y1max.or.twolat(1).lt.y1min)then
	i1=ns2; i2=1; i3=-1
	else
	i1=1; i2=ns2; i3=1
	endif
	do i=i1,i2,i3
	stlat=twolat(i); stlon=twolon(i)
	if(stlat.le.y1max.and.stlat.ge.y1min)then
	call mindist(dist,azm,disb,azb,rx,stlat,stlon,solat,solon,
     +                   ns,njump,dipf,wid,ftop,rb)
        avgdist=avgdist+disb
c        print *,stlat,stlon,disb,i,' 2nd flt lat/long dist count'
        maxd=max(maxd,disb)
        mind=min(mind,disb)
        nc=nc+1
        elseif(disb.gt.0.)then
c continue using disb as the distance
	avgdist=avgdist+disb
	nc=nc+1
	else
	ndump=ndump+1
	endif
        enddo
        avgdist=avgdist/nc
        print *,'avg distance is ',avgdist, ' km'
        print *,'max distance is ',maxd
        print *,'min distance is ',mind
        print *,'number of points on secondary fault omitted from dist calcs:',ndump
        
        do i=1,5	!master flt outer loop
        wtm=wt1(i)
        fac=avgdist*sin(dip1(i))
        do j=1,5	!secondary flt inner loop
        wtnow=wtm*wt2(j)
        ang=pi-dip1(i)-dip2(j)
        h=fac/sin(ang)
        z=h*sin(dip2(j))
        avgdip=avgdip+dip2(j)*wtnow
        slip=slip2*0.001/sin(dip2(j))	!m/year from mm
        sumwt=sumwt+wtnow
c B&R: max seismogenic depth 15 km
        if(z.gt.15.)then
        z=15.
        h=z/sin(dip2(j))	!h is slant distance or hypotenuse.
        endif
c Wells and Coppersmith M(A) where A is the area of the fault (km^2)
        area=flen*h
        maga= 3.98+1.02*alog10(area)
	if(.not.fixrate)then
        xmo=xmu*area*slip*1.e6	!moment rate (Nm/yr)
c slip-based case: event rate based on moment rate
        rate=xmo/amo(maga)
        else
c fixed-rate case: moment rate prop. to event rate
        xmo=amo(maga)*rate
        endif
        arate=arate+rate*wtnow
        dip=dip2(j)/d2r
        if(dip.gt.70.)then
        islip=1
        else
        islip=3
        endif
        write(3,20)h,z,area,maga,wtnow,dip1(i)/d2r,dip,xmo,nint(1./rate)
20	format(f6.2,5x,f6.2,4x,f8.2,1x,f6.2,1x,f8.5,2(1x,f6.1),1x,e10.5,1x,i6)
	write(4,30)islip,name2
30	format('1 ',i1,' 1 1  ',a)
	write(4,35)maga,rate,wtnow,xmo,dip1(i)/d2r,name1
35	format(f4.2,1x,1pe11.5,1x,0p,f7.5,1x,e10.5,1x,f4.1,1x,a)
	write(4,40)dip,h,0.,flen
40	format(f6.1,1x,f6.2,1x,f3.1,1x,f7.2)
	write(4,*)ns2
	do k=1,ns2
	write(4,45)twolat(k),twolon(k)
45	format(f11.6,1x,f12.6)
	enddo
	ymo=ymo+xmo*wtnow
	enddo
	enddo
	write(3,48)avgdist,nint(1./arate),fixrate
48	format(/,'#average flt separation ',f6.2,' km. Avg RI =',i6,' yrs. Fixed event rate? ',l1)
	avgdip=avgdip/d2r/sumwt
	print *,'Average dip, secondary flt: ',avgdip,' sumwt ',sumwt
	print *,'Average moment rate on secondary flt(Nm): ',ymo/sumwt
	print *,'Average recurrence interval (yrs) on 2nd flt: ',sumwt/arate
	print *,'Number of points on master flt that were omitted ',nomit
	close(4)
	stop 'standard exit'
345	print *,'An input file on command line not found .'
	stop 'please check that it is in the working directory or folder'
456	print *,x,y
	stop 'The master fault coordinates did not read in correctly'
        end
	
	
      subroutine mindist(dist,azm,disb,azb,rx,stlat,stlon,solat,solon,
     +                   ns,njump,dip,wid,ftop,rb)
c.....................................................................
      parameter (d2r = 0.0174533)
c  input:
c         stlat - station latitude
c         stlon - station longitude
c         solat - latitude of fault segment corners
c         solon - longitude of fault segment corners
c         ns    - number of fault segment corners
c         njump - the non-rupture segment order number, if no non-rupture
c                 segment, set njump=0
c         dip   - fault dip angle
c         wid   - fault width
c         ftop  - depth of fault top 
c  output:
c         dist  - minimum distance to rutpure plane
c         azm   - azimuth for rupture distance
c         disb  - JB distance to rutpure plane
c         azm   - azimuth for JB distance
c         rx    - Chiou and Yong's distance to the surface projection 
c                 of top edge of fault
c         rb    - Bartlett's distance to the surface projection of 
c                 top edge of fault
c                                              Yuehua Zeng, 2007
c.....................................................................
      dimension solat(*),solon(*)
c
      dp=dip*d2r
      cosdp=cos(dp)
      sindp=sin(dp)
      dist=1.e10
      disb=1.e10
      call delaz2(solat(1),solon(1),solat(ns),solon(ns),flen,str1)
      call delaz2(solat(1),solon(1),stlat,stlon,dista,az)
      sx=dista*cos(az-str1)
      y=dista*sin(az-str1)
      rb=sqrt(dista*dista+flen*(flen-sx-sx))
      if(sx.le.0.0.or.sx.ge.flen)then; rx=y; else; rx=sign(rb,y); endif
      do i=1,ns-1
         if(i.ne.njump)then
         call delaz2(solat(i),solon(i),solat(i+1),solon(i+1),flen,str)
         call delaz2(solat(i),solon(i),stlat,stlon,dista,az)
         sx=dista*cos(az-str)
         if(sx.gt.flen)then
           sx=sx-flen
         elseif(sx.gt.0.0)then
           sx=0.0
         endif
         y=dista*sin(az-str)
         r=sqrt(sx*sx+y*y)
         if(r.lt.abs(rx))rx=sign(r,y)
         if(r.lt.rb)rb=r
c
c  for JB distance
         sy=y
         if(sy.gt.wid*cosdp)then
           sy=sy-wid*cosdp
         elseif(sy.gt.0.0)then
           sy=0.0
         endif
         r=sqrt(sx*sx+sy*sy)
         if(r.lt.disb)then
           disb=r
           azb=az/d2r
         endif
c
c  for rupture distance
         sy=y*cosdp-ftop*sindp
         if(sy.gt.wid)then
           sy=sy-wid
         elseif(sy.gt.0.0)then
           sy=0.0
         endif
         sz=y*sindp+ftop*cosdp
         r=sqrt(sx*sx+sy*sy+sz*sz)
         if(r.lt.dist)then
           dist=r
           azm=az/d2r
         endif
         endif
      end do
c
      return
      end subroutine mindist
c
      subroutine delaz2(sorlat,sorlon,stnlat,stnlon,delta,az)
      parameter (coef= 3.1415926/180.)
      xlat= sorlat*coef
      xlon= sorlon*coef
      st0= cos(xlat)
      ct0= sin(xlat) 
      phi0= xlon
      xlat= stnlat*coef
      xlon= stnlon*coef
      ct1= sin(xlat)
      st1= cos(xlat)
      sdlon= sin(xlon-phi0)
      cdlon= cos(xlon-phi0)
      cdelt= st0*st1*cdlon+ct0*ct1
      x= st0*ct1-st1*ct0*cdlon
      y= st1*sdlon
      sdelt= sqrt(x*x+y*y)
      delta= atan2(sdelt,cdelt)
      delta= delta/coef
      az= atan2(y,x)
c     az= az/coef
      delta= delta*111.2
      return
      end    subroutine delaz2

      subroutine mindist1(dist,azm,disb,azb,rx,stlat,stlon,solat,solon,
     +                   ns,njump,dip,wid,ftop,rb)
c.....................................................................
c  For calculating JB or rupture distance to a curved fault defined 
c  by extending the fault top segment trace down its dip direction
c
c  input:
c         stlat - station latitude
c         stlon - station longitude
c         solat - latitude of fault segment corners
c         solon - longitude of fault segment corners
c         ns    - number of fault segment corners
c         njump - the non-rupture segment order number, if no non-rupture
c                 segment, set njump=0
c         dip   - fault dip angle
c         wid   - fault width along the dip direction. wid>0 is a good constraint
c         ftop  - depth of fault top 
c  output:
c         dist  - minimum distance to rutpure plane
c         azm   - azimuth for rupture distance
c         disb  - JB distance to rutpure plane
c         azb   - azimuth for JB distance
c         rx    - Chiou and Yong's distance to the surface projection 
c                 of top edge of fault
c         rb    - Bartlett's distance to the surface projection of 
c                 top edge of fault
c                                              Yuehua Zeng, 2007
c.....................................................................
      parameter (d2r = 3.1415926/180.)
      dimension solat(*),solon(*)
c
      dp=dip*d2r
      cosdp=cos(dp)
      sindp=sin(dp)
      dist=1.e10
      disb=1.e10
      call delaz2(solat(1),solon(1),solat(ns),solon(ns),flen,str1)
      call delaz2(solat(1),solon(1),stlat,stlon,dista,az)
      sx=dista*cos(az-str1)
      yy=dista*sin(az-str1)
      rb=sqrt(dista*dista+flen*(flen-sx-sx))
      if(sx.le.0.0.or.sx.ge.flen)then; rx=yy
      else; rx=sign(sqrt(dista*dista+flen*(flen-sx-sx)),yy); endif
      do i=1,ns-1
         if(i.ne.njump)then
         call delaz2(solat(i),solon(i),solat(i+1),solon(i+1),flen,str)
         call delaz2(solat(i),solon(i),stlat,stlon,dista,az)
c
         x=sindp
         y=cosdp*cos(str-str1)
         sindpp=x/sqrt(x*x+y*y)
         cosdpp=sqrt(1.0-sindpp*sindpp)
c
         sx=dista*cos(az-str)
         yy=dista*sin(az-str)
c
         if(sx.ge.0.0.and.sx.le.flen.and.abs(yy).lt.abs(rx))then;rx=yy
         elseif(dista.lt.abs(rx))then;rx=sign(dista,yy);endif
         if(sx.ge.0.0.and.sx.le.flen.and.abs(yy).lt.rb)then;rb=abs(yy)
         elseif(dista.lt.rb)then;rb=dista;endif
c
c  for JB distance
         sy=yy
         wid1=wid*cosdp
c
         cosfi=sin(str-str1)
         sinfi=sqrt(1.0-cosfi*cosfi)+1.0e-10
         syy=sy/sinfi
         sxx=sx-syy*cosfi
         r=-10.0
c
         if(sxx.lt.0.0.or.sxx.gt.flen)then
           y=sxx*cosfi+syy
           if(sxx.gt.flen)y=(sxx-flen)*cosfi+syy
           if(y.lt.0.0.or.y.gt.wid1)then
             if(y.gt.wid1)syy=syy-wid1
             x=syy*cosfi+sxx
             if(x.gt.flen)then
               sxx=sxx-flen
             elseif(x.ge.0.0)then
               r=abs(syy)*sinfi
               az0=90.0*sign(1.0,syy)
             endif
           else
             if(sxx.gt.flen)sxx=sxx-flen
             r=abs(sxx)*sinfi
             az0=(str1-str)/d2r-90.0+90.0*sign(1.0,sxx)
           endif
         elseif(syy.lt.0.0.or.syy.gt.wid1)then
           if(syy.gt.wid1)syy=syy-wid1
           x=syy*cosfi+sxx
           if(x.gt.flen)then
             sxx=sxx-flen
           elseif(x.ge.0.0)then
             r=abs(syy)*sinfi
             az0=90.0*sign(1.0,syy)
           endif
         else
           r=0.0
           az0=90.0
         endif
c
         if(r.lt.-1.0)then
           sxx=sxx+syy*cosfi
           syy=syy*sinfi
           r=sqrt(sxx*sxx+syy*syy)
           az0=atan2(syy,sxx)/d2r
         endif
         if(r.lt.disb)then
           disb=r
           azb=str/d2r+az0
         endif
c
c  for rupture distance
         sy=yy*cosdpp-ftop*sindpp
         sz=yy*sindpp+ftop*cosdpp
c
         cosfi=sin(str-str1)*cosdp
         sinfi=sqrt(1.0-cosfi*cosfi)+1.0e-10
         syy=sy/sinfi
         sxx=sx-syy*cosfi
         r=-10.0
c
         if(sxx.lt.0.0.or.sxx.gt.flen)then
           y=sxx*cosfi+syy
           if(sxx.gt.flen)y=(sxx-flen)*cosfi+syy
           if(y.lt.0.0.or.y.gt.wid)then
             if(y.gt.wid)syy=syy-wid
             x=syy*cosfi+sxx
             if(x.gt.flen)then
               sxx=sxx-flen
             elseif(x.ge.0.0)then
               r=abs(syy)*sinfi
               az0=90.0*sign(1.0,syy)
             endif
           else
             if(sxx.gt.flen)sxx=sxx-flen
             r=abs(sxx)*sinfi
             az0=(str1-str)/d2r-90.0+90.0*sign(1.0,sxx)
           endif
         elseif(syy.lt.0.0.or.syy.gt.wid)then
           if(syy.gt.wid)syy=syy-wid
           x=syy*cosfi+sxx
           if(x.gt.flen)then
             sxx=sxx-flen
           elseif(x.ge.0.0)then
             r=abs(syy)*sinfi
             az0=90.0*sign(1.0,syy)
           endif
         else
           r=0.0
           az0=90.0
         endif
c
         if(r.gt.-1.0)then
           r=sqrt(r*r+sz*sz)
         else
           sxx=sxx+syy*cosfi
           syy=syy*sinfi
           r=sqrt(sxx*sxx+syy*syy+sz*sz)
           az0=atan2(syy*cosdpp+sz*sindpp,sxx)/d2r
         endif
         if(r.lt.dist)then
           dist=r
           azm=str/d2r+az0
         endif
         endif
      end do
c
      return
      end subroutine mindist1
