
c *----------------- BEGIN SMSIM_TDWW2 -----------------------------
      subroutine smsim_tdww2(per, nper, targpsa,sa5,ts,
     :    prvsim_m1,pgasim, pgasim_std, avgpsa,
     :    avgscafac, nacc_save, acc_save)

c * Uses time-domain stochastic model to compute ground motion
c *  (SMSIM = Stochastic Model Simulation) 
c *** targpsa = sought-after psa or pga from PSHA. 

c * Dates:  05/10/95 - Written by D. M. Boore
c *         05/30/95 - Modified by Basil Margaris to allow choice of window
c *         06/02/95 - Further renamed (from psvsimnu) and modified.
c *         06/08/95 - Normalize by sqrt(avg square amp)
c *         06/09/95 - Get r, amag from common, not parameter list
c *         06/12/95 - Removed writing of acc, vel to a file; this should
c *                    be done in the front-end driver program.
c *         08/08/95 - Removed dc from noise segment before going to
c *                    the frequency domain.
c *         08/08/95 - Changed things so that remove dc from noise sample before
c *                    applying the window.
c *         08/18/95 - Added subroutine Write_Params
c *         10/17/95 - Added flag to Get_Params and Write_Params
c *         10/17/95 - Remove dc or not from random series, depending on a flag
c *                    (irmvdc) that is passed through the include statement
c *         11/14/95 - Added call to const_am0_gsprd
c *         11/14/95 - Broke out numerical recipes and subroutines used
c *                    by both rv and td into separate collections
c *                    of Fortran subroutines
c *         11/16/95 - Replaced 'shape = ...' statement with call to spect_amp
c *         12/06/95 - Major simplification in exponential window, removing
c *                    tapers (the window itself provides tapers) and correcting
c *                    a bug so that the window can be used for a duration
c *                    other than the duration used to determine the window 
c *                    parameters.
c *         12/14/95 - Assigned durex the source plus path duration
c *         12/22/95 - Removed rmean real part of work only
c *         01/22/96 - Use REALFT rather than FORK
c *         04/12/96 - Changed names psv, sd to prv, rd
c *         04/15/96 - Changed "rd_spect" to "rd_calc"
c *         08/01/98 - Calculate standard deviation peak values
c *         09/11/98 - Increased array dimensions from 16400 to 33000
c *         01/13/99 - Added computation of Arias intensity
c *         02/09/99 - Added computation of displacement; removed subtracting
c *                    linear trend from acc in acc2vd routine
c *         03/05/99 - Removed include rvtdsubs.for, recipes.for from the
c *                    end; these must now be included in the driver programs.
c *         07/02/99 - Added arguments to durpath to incorporate Atkinson and 
c *                    Silva (1999) magnitude-dependent modification to distance.
c *         08/03/99 - Added uniform deviate if desired
c *         02/10/99 - Transfer choice of npts from rvtdsubs.for to this
c *                    subroutine.  This also includes pad extension if a low-cut
c *                    filter is used.  
c *                    In all cases 
c *                    npts*dt >= max(tshift, tfpad) + dur_fctr * nmotion * dt 
c *                               + tfpad 
c *         03/12/99 - Suggest increasing dt if the time series is too long.
c *         01/07/00 - Upgraded to most recent version of acc2vd
c *         01/27/01 - Modified time series window computation (new parameters)
c *         01/28/01 - Increased dimensions to 66000 (should use allocatable 
c *                    arrays)
c *         02/12/01 - Combine smsim_td and get_acc, and use dynamically 
c *                    allocatable arrays.  Put calculation of ndimen into
c *                    a subroutine.  Rename "npts" to "npw2".
c *          04/27/01 - Harmsen: output up to nsavmx accelerograms scaled to PSHA
c 			input level for a specific period (Per(1)).
c ***        11/09/01 - Harmsen. try to approximate "approximate UHS" by also
c ***                     considering sa5, the ordinate where that spectrum is flat
c ***		4/2003- compute log averages, do not scale but do choose best fit
c * ALLOCATABLE: comment out following lines to disable dynamic allocation
c      real acc(:), vel(:), dis(:), a_sq_int(:), work(:)
c      allocatable :: acc, vel, dis, a_sq_int, work
c * Remove comment character in lines below to disable dynamic allocation
	parameter (nsavmx=60)
      real acc(66000)
c	real vel(66000), dis(66000), a_sq_int(66000),
      real work(66000),psasav(nsavmx)
	common/psasav/psasav
      real avgpsa, per(*), prvsim_m1(*), 
     :     acc_save(66000,nsavmx),targpsa(*),facrv(nsavmx)
c can save up to 10 seismograms for web run..
      logical rmv_trend

c      character prvcalc*1

c      include '\smsim\smsim.fi'
	include 'smsim.fi'
      
c * pi, twopi now computed in smsim.fi
c *      pi = 4.0 * atan(1.0)
c *      twopi = 2.0 * pi
c	print *,per(1),' =per(1) in smsim_tdwww'
      do i = 1, nper
c        prvcumsq(i) = 0.
      facrv(i)=1.0
c facrv is a www scaling factor. 1 for pga, othewise computed with prv below
c        prvsim_std(i) = 0.
      end do
      targprv=per(1)/twopi*targpsa(1)
c seismograms will be scaled so that at the period specified SA matches PSHA SA
      avgpga = 0.
c      avgpgv = 0.
c      avgpgd = 0.
c      avgasqint = 0.
c      pgacumsq = 0.
c      pgvcumsq = 0.
c      pgdcumsq = 0.
c      asqintcumsq = 0.
      pgasim_std = 0.
	avgpsa = 0.0
c      pgvsim_std = 0.
c      pgdsim_std = 0.
c      arias_std = 0.

      iseed = -abs(seed)

      call Get_Npts(nstart, nstop, npw2, te)

c * ALLOCATABLE: comment out following lines to disable dynamic allocation
c      ndimen = npw2
c      allocate ( acc(ndimen), vel(ndimen), dis(ndimen), 
c     :          a_sq_int(ndimen), work(ndimen) )
c * Remove comment character in lines below to disable dynamic allocation
      ndimen = 66000
      if (npw2 .gt. ndimen) then
        write(*,'(2x,a,i5,a,i5,a)') 
     :       ' *** FATAL ERROR: NPW2 (',
     :                                    npw2, 
     :       ') > NDIMEN (', 
     :                                     ndimen, 
     :       '); QUITTING!'
        write(*, '(2x,a)') ' TRY INCREASING DT'
        stop
      end if
c * ALLOCATABLE: end of group

c      write(*,*)
	avgscafac=0.0
      do isim = 1, nacc_save	!was nsims. here more are computed for per>0
c than for pga. so we change bound to nacc_save
c
c         write(*,'(3(a,i3),a)') 
c     :    '+ Patience! Computing accelerogram # ', isim,
c     :    ' of ', nsims,
c     :    ' and rs at ', nper,' periods.'

c * Initialize arrays:
         do i = 1, ndimen
           acc(i)=0.
           work(i) = 0.
         end do

c * Generate white noise with zero mean, unit variance, and proper distribution:

        if (iran_type .eq. 0) then
          do i = nstart, nstop
            work(i) = gasdev(iseed)
          end do
        else
          do i = nstart, nstop
            work(i) = ran1(iseed) - 0.5
          end do
        end if

        if (irmvdc .eq. 0) then
          rmean = 0.0
        else
          call mean(work, nstart, nstop, rmean)
        end if

c * Window the white noise:

         if(indxwind .eq. 0) then        ! BOX WINDOW
   
           do i = nstart, nstop
              work(i) = wind_box(i, nstart, nstop, ntaper) *
     :                  (work(i) - rmean)
           end do

         else                           ! EXPONENTIAL WINDOW

           do i = nstart, nstop
             t = float(i-nstart) * dt
             work(i) = wind_exp( t, te, eps_w, eta_w, new_mr) *
     :                 (work(i) - rmean)
           end do

         end if

c * Transform to frequency domain:
         call realft(work,npw2,+1)

c * Find sqrt of squared spectral amplitude for normalization:
         call avgsq_realft(work,npw2,avgsqamp) 
         sqrt_avgsqamp = sqrt(avgsqamp)

c * Get frequency independent parameter:
         call const_am0_gsprd()     ! returns freq_indep_factor through common

         iaorins = 1             ! ground motion
         idva = 2                ! acceleration

         df = 1.0/(float(npw2)*dt)
         do i = 2,npw2/2
             f = (i-1) * df
             sfact = spect_amp(f)/ sqrt_avgsqamp
             work(2*i-1) = sfact * work(2*i-1)
             work(2*i)   = sfact * work(2*i)
         end do
         work(1) = 0.0              ! OK for acceleration spectrum, but not for
                                 ! displacement spectrum
         fnyq = 1.0/(2.0*dt)
         work(2) = spect_amp(fnyq) * work(2)/ sqrt_avgsqamp

c * Transform back to time domain:
         call realft(work,npw2,-1)

c * Apply FFT normalization:
         afact = 2 * df
         do i = 1, npw2
            acc(i) = afact * work(i)
         end do

c *         rmv_trend = .true.
         rmv_trend = .false.
c         v0 = 0.0
c         d0 = 0.0
c         call acc2vd(acc, npw2, dt, rmv_trend, v0, d0, vel, dis) 

c         call accsqint(acc, npw2, dt, .false., a_sq_int)

         
c * Compute pga, pgv, pgd, asqint:

         call mnmax(acc, 1, npw2, 1, amin, amax)
          pga = amax1(abs(amin), abs(amax))

c         call mnmax(vel, 1, npw2, 1, vmin, vmax)
c         pgv = amax1(abs(vmin), abs(vmax))

c         call mnmax(dis, 1, npw2, 1, dmin, dmax)
c         pgd = amax1(abs(dmin), abs(dmax))

c         asqint = a_sq_int(npw2)
c assume lognormal distribution, find median estimate
         avgpga = avgpga + alog(pga)
c         avgpgv = avgpgv + pgv
c         avgpgd = avgpgd + pgd
c         avgasqint = avgasqint + asqint
c         pgacumsq = pgacumsq + pga*pga
c         pgvcumsq = pgvcumsq + pgv*pgv
c         pgdcumsq = pgdcumsq + pgd*pgd
c         asqintcumsq = asqintcumsq + asqint*asqint
 
         if (per(1).gt.0.) then
c * Compute prv & psa : nper is 1 for www runs
           do i = 1, nper
              omega = twopi/per(i)
              call rd_calc(acc, npw2, omega, damp, dt, rd)
              prv = omega * rd
		psa= omega * prv
c		print *,omega,prv,psa
c units of targprv cm/sec, units of prv cm/sec.
	      facrv(i)=targpsa(i)/psa
              avgpsa = avgpsa + alog(psa)
		if(i.eq.1)psasav(isim)=psa
          end do 
           else
           facrv(1)=targpsa(1)/pga                  ! iper=0
         end if
         if (isim .le. nacc_save) then
c         print *,isim,facrv(1)
	avgscafac=avgscafac+alog(facrv(1))
           do i = 1, npw2
c             acc_save(i,isim) = facrv(1)*acc(i)
		acc_save(i,isim) = acc(i)
c             vel_save(i) = vel(i)
c             dis_save(i) = dis(i)
           end do
           if(ts.lt.1.0)then
c ts should be the corner of the design spectrum if everything OK           
           t0=ts
           else
           t0=1./2.5
           endif
           th=0.2*t0
        if(per(1).gt.0.0)then
c       find l1 measure of uniformity in band 0.2 to 1 sec
	perx=0.1
        do iperx=1,9
        perx=perx+0.1
        wt=1.0
        if(per(1).eq.perx)wt=3.0
c triple the penalty for misfit at specified SA period        
              omega = twopi/perx
        call rd_calc(acc_save(1,isim),npw2,omega,damp,dt,rd)
        prv=omega*rd

c t0 is approx UH velocity spectrum long-period corner. ordinate units cm/sec
        if(perx.lt.t0)then
        ord=sa5*perx/twopi
        else
        ord=sa5*t0/twopi
        endif
        prvsim_m1(isim)=prvsim_m1(isim)+wt*abs(alog(ord/prv))
        enddo
c below ordinate is at the shortest period, th, where approx UH accel spectrum 
c is flat. This is an addition of Nov 14,2001.        
        omega=twopi/th
        call rd_calc(acc_save(1,isim),npw2,omega,damp,dt,rd)
        prv=omega*rd
	ord=sa5*th/twopi
        prvsim_m1(isim)=prvsim_m1(isim)+abs(alog(ord/prv))
	        
        endif	!if(per(1) > 0)
         end if

      end do                           ! isim

c * Compute averages
c change "nsims" to nacc_save in denom below
	avgscafac=exp(avgscafac/float(nacc_save))
	avgpsa=exp(avgpsa/float(nacc_save))
       pgasim = exp(avgpga/float(nacc_save))
c      pgvsim = avgpgv/float(nsims)
c      pgdsim = avgpgd/float(nsims)
c      asqint = avgasqint/float(nsims)

c      if (prvcalc .eq. 'y' .or. prvcalc .eq. 'Y') then
c        do i = 1, nper
c           prvsim(i) = avgprv(i)/float(nsims)
c        end do
c      end if

c * Compute standard deviations.  Note that this is a poor way to do so
c * (because of subtracting two possibly large numbers), but
c * it avoids having to save the values in an array.  This is particularly
c * important for the prv array, which could be quite large.  Hopefully
c * using double precision will minimize the roundoff problem.  Also note that
c * the standard deviation is not too meaningful; it does not capture
c * variability associated with variability of the input parameters.

c      if (nsims .gt. 1) then
c        factr = float(nsims)/float(nsims-1)
c        pgasim_std = dsqrt(factr * (pgacumsq/nsims - pgasim*pgasim))
c        pgvsim_std = dsqrt(factr * (pgvcumsq/nsims - pgvsim*pgvsim))
c        pgdsim_std = dsqrt(factr * (pgdcumsq/nsims - pgdsim*pgdsim))
c        asqint_std = dsqrt(factr * (asqintcumsq/nsims-asqint*asqint))
c        if (prvcalc .eq. 'y' .or. prvcalc .eq. 'Y') then
c          do i = 1, nper
cc             prvsim_std(i) = 
c     :        dsqrt(factr * (prvcumsq(i)/nsims - 
c     :                      prvsim(i)*prvsim(i)))
c          end do
c        end if
c      end if        
c
c      g = 980.0    ! acceleration of gravity, assuming acc units of cm/s^2

c      arias_fctr = pi/(2.0*g)
c      arias = arias_fctr * asqint
c      arias_std = arias_fctr * asqint_std 
      
c * ALLOCATABLE: comment out following line to disable dynamic allocation
c      deallocate (acc, vel, dis, a_sq_int, work)  
c * ALLOCATABLE: end of group

      return
      end
c *----------------- END SMSIM_TDWW2 -----------------------------

c * ------------------------------------------------------------------ Get_Npts
      subroutine Get_Npts(nstart, nstop, npw2, te)

c * Calculate various indices and te (=t_eta)

c * Dates: 02/12/01 - Extracted from smsim_td.

      include 'smsim.fi'
c      include '\smsim\smsim.fi'

c * Calculate the number of points in the noise sample:

      call spect_scale()  ! call this now because need corner frequencies
                          ! for window durations

      if(indxwind .eq. 0) then        ! BOX WINDOW

        durex = dursource(w_fa, fa, w_fb, fb)+
     :          durpath(r, nknots, rdur, dur, slast, numsource, amag)
        nmotion = durex/dt
        ntaper = ifix(taper*nmotion)  ! Increase duration of motion by tapers 
                                      ! (figure 1/2 front and back):
                                      ! taper is fractional, not percent
        nmotion = nmotion + ntaper  

      else                           ! EXPONENTIAL WINDOW

        durex = dursource(w_fa, fa, w_fb, fb)+
     :          durpath(r, nknots, rdur, dur, slast, numsource, amag)
        tb = durex
        te = f_tb2te * tb
        tw = f_te_xtnd * te
        nmotion = tw/dt   
      end if

c * Calculate nstart and nstop, depending on tshift and if a
c *   low-cut filter is used.

      if (fcut .eq. 0.0) then
        tfpad = 0.0
      else
        tfpad = 1.5 * (norder/2) / fcut ! (Converse, USGS OFR 92-296A, p. 2-3)
      end if

      if (tfpad .gt. tshift) tshift = tfpad

      nstart = tshift/dt + 1

      nstop = nstart + nmotion

c * Calculate npts, depending on tshift, dur_fctr, and if a
c *   low-cut filter is used.

c * compute smallest power of two for which the resulting duration is greater
c * than or equal to tsimdur:

      tsimdur = nstart*dt + dur_fctr * nmotion * dt + tfpad

      npw2 = 2.0**ifix(alog10(tsimdur/dt + 1.0)/alog10(2.0))
      if ( (npw2-1)*dt .lt. tsimdur) npw2 = 2 * npw2

      return
      end
c * ------------------------------------------------------------------ Get_Npts

c *----------------- BEGIN Acc2VD -----------------------------
      subroutine acc2vd(acc, npts, dt, rmv_trnd, v0, d0, vel, dis)


c * Compute velocity and displacement time series from acceleration,
c * assuming that the acceleration
c * is represented by straight lines connecting the digitized values.

c * Dates: 02/09/99 - Written by D.M. Boore
c *        01/07/00 - Added initial velocity and displacements (v0, d0).
c *                   This also requires the addition of a linear trend in
c *                   displacement.
c *                   Also, Bill Joyner, Chris Stephens, and I considered how to
c *                   handle the first point.  Assuming v0 = 0, BAP uses a 
c *                   trapezoidal rule such that the first vel(1) = 0.5*dt*a(1),
c *                   but in acc2vd, vel(1) = 0.  In effect, BAP starts the
c *                   integration at -dt rather than 0.0, with the result that
c *                   because vel(1) depends on dt, vel(1) will change if the
c *                   digitization interval is changed.  We agreed that this is 
c *                   not correct.

      real acc(*), vel(*), dis(*) 
      logical rmv_trnd
      double precision cumv, cumd, a1, a2, v1, 
     : ddt, ddt_2, ddtdt_6

      if (rmv_trnd) then      
c * remove trend first (straight line between first and last points)
c * Note: acc is replaced with detrended time series
c *        call dcdt(acc, dt, npts, 1, npts, .false., .true.)  ! old routine,
c *                                                         ! gives steps at ends

         call rmvtrend(acc, npts)
      end if

c * compute velocity and displacement, using analytical formulas based
c * on representing the acceleration as a series of straightline segments.

      ddt     = dble(dt)
      ddt_2   = dble(dt/2)
      ddtdt_6 = dble(dt**2/6)

      cumv = 0.0
      cumd = 0.0

      vel(1) = sngl(cumv) + v0
      dis(1) = sngl(cumd) + d0
      do j=2,npts
        a1 = acc(j-1)
        a2 = acc(j)
        v1 = vel(j-1)
        cumv = cumv + (a1 + a2)*ddt_2
        vel(j) = sngl(cumv) + v0
        cumd = cumd + v1*ddt + (2.0*a1 + a2)*ddtdt_6
        dis(j) = sngl(cumd) + d0  ! no linear trend neede; it's include in v1
      end do

      return
      end
c *----------------- END Acc2VD -----------------------------

c * ----------------------------- BEGIN RMVTREND ----------------
      subroutine rmvtrend(y, n)

c * Removes a straightline fit to first and last points, replacing
c * the input array with the detrended array

c * Dates: 02/09/99 - written by D. Boore


      real y(*)

      y1 = y(1)
      y2 = y(n)
      slope = (y2 - y1)/float(n-1)

      do i = 1, n
        y(i) = y(i) - (y1 + slope*float(i-1))
      end do

      return
      end
c * ----------------------------- END RMVTREND ----------------

c *----------------- BEGIN DCDT -----------------------------
c * Dates:  - Written by C.S. Mueller
      SUBROUTINE DCDT (Y,DT,NPTS,INDX1,INDX2,LDC,LDT)
c+
c  DCDT - Fits DC or trend between indices INDX1 and INDX2.
C         Then removes DC or detrends whole trace.
c         Y is real, DT = delta t.
c         If remove DC, LDC = .TRUE.
c         IF detrend, LDT = .TRUE.
c-
      real Y(1)
      logical LDC,LDT
c
c...Fit DC and trend between indices INDX1 and INDX2.
100   NSUM = INDX2-INDX1+1
      SUMX = 0.0
      SUMX2 = 0.0
      SUMY = 0.0
      SUMXY = 0.0
      DO 200 I=INDX1,INDX2
         XSUBI = (I-1)*DT
         SUMXY = SUMXY+XSUBI*Y(I)
         SUMX = SUMX+XSUBI
         SUMX2 = SUMX2+XSUBI*XSUBI
200      SUMY = SUMY+Y(I)
C
C... Remove DC.
300   IF (LDC) THEN
        AVY = SUMY/NSUM
        DO 360 I=1,NPTS
360        Y(I) = Y(I)-AVY
        RETURN
      END IF
C
C... Detrend. See Draper and Smith, p. 10.
400   IF (LDT) THEN
        BXY = (SUMXY-SUMX*SUMY/NSUM)/(SUMX2-SUMX*SUMX/NSUM)
        AXY = (SUMY-BXY*SUMX)/NSUM
        QXY = DT*BXY
        DO 450 I=1,NPTS
450        Y(I) = Y(I)-(AXY+(I-1)*QXY)
        RETURN
      END IF
C
      STOP
      END
c *----------------- END DCDT -----------------------------

c *----------------- BEGIN MNMAX -----------------------------
      subroutine mnmax(a,nstrt,nstop,ninc,amin,amax) 
c
c author: D. M. Boore
c last change: 9/7/84
c
      dimension a(1)
      amax = a( nstrt)
      amin=amax         
      do 10 i=nstrt,nstop,ninc
      if(a(i)-amax) 15,15,20  
20    amax=a(i)               
      go to 10                
15    if(a(i)-amin) 25,10,10  
25    amin=a(i)               
10    continue                
      return                  
      end                     
c *----------------- END MNMAX -----------------------------

c *----------------- BEGIN MEAN -----------------------------
      subroutine mean(work, nstart, nstop, rmean)
c *  Dates: 01/22/96 - Written by D. Boore
      real work(*)
c * find mean of the array:
      sum = 0.0
      do i = nstart, nstop
        sum = sum + work(i)
      end do
      rmean = sum/float(nstop-nstart+1)
      return
      end
c *----------------- END MEAN -----------------------------

c *------------------- BEGIN AVGSQ_REALFT -------------------------------
c * Dates: 01/22/96 - Written by D.M. Boore
      subroutine avgsq_realft( s, npts, avgsqamp)
      real s(*)
      sum=0.
      do j = 2, npts/2                   !don't include the dc or Nyquist values
        sum=sum + s(2*j-1)**2 + s(2*j)**2 ! odd, even = real, imag spect
      end do
      avgsqamp = sum/float(npts/2 - 1)
      return
      end
c *  ------------------- END AVGSQ_REALFT -------------------------------

c *  ------------------- BEGIN WIND_BOX -------------------------------
      function wind_box( i, nstart, nstop, ntaper)
c
c applies cosine tapered window.
c unit amplitude assumed
c
c written by D. M. Boore
c
c latest revision: 9/26/95
      real wind_box

      wind_box = 0.0
      if ( i .lt. nstart .or. i. gt. nstop) return
      wind_box = 1.0
      if ( i .ge. nstart+ntaper .and. i .le. nstop-ntaper ) return
c
      pi = 4.0 * atan(1.0)
c
      dum1 = (nstop+nstart)/2.0
      dum2 = (nstop-nstart-ntaper)/2.0
c
      wind_box = 0.5 * (1.0 - sin( pi*
     *  ( abs(float(i)-dum1) - dum2 ) /float(ntaper) ) )
      return
      end
c *  ------------------- END WIND_BOX -------------------------------
c *  ------------------- BEGIN WIND_EXP ---------------------------
      function  wind_exp( t, te, eps_w, eta_w, new_mr)
c
c     apply Sargoni and Hart (1974) window, with parameters
c     tw, eps (fraction of tw to reach peak), and
c     eta ( fraction of peak ampl. at tw).  See Boore (BSSA, 1983,
c *     p. 1869).  Note that t can be larger than tw.

c * Dates:
c *         05/30/95 - Initially written by Dave Boore, based on a routine
c *                    by G. Atkinson but with significant structural 
c *                    changes and the taper to be a raised cosine taper
c *         12/06/95 - Removed the taper and simplified the parameters in the
c *                    calling list.  Also added a switch so that the window
c *                    shape coefficients are only computed once.  This assumes
c *                    that the shape coefficients will be the same for
c *                    any run of the program... a good assumption.
c *         12/28/95 - Used new_mr, set from driver, to control computation
c *                    of b,c,a; before I used "firstcall", set in this 
c *                    subprogram, but this gave an error if the driver
c *                    looped over multiple values of m and r.
c *         01/27/01 - Modified time series window computation (changed variable
c *                    names, use normalized time)


      real wind_exp
      logical new_mr
      save a, b, c

      if (new_mr) then
        b = -eps_w * alog(eta_w)/
     :      (1. + eps_w*(alog(eps_w)-1.))
        c = b/eps_w
        a = (exp(1.0)/eps_w)**b
        new_mr = .false.
      end if

      wind_exp = 0.
      if( t .lt. 0.0) return
c
c     Apply Sargoni and Hart window.
c
      wind_exp = a*(t/te)**b * exp(-c*(t/te))

      return
      end
c *------------------- END WIND_EXP ---------------------------

c *----------------- BEGIN RD_CALC -----------------------------
      subroutine rd_calc(acc,na,omega,damp,dt,rd)
c * This is a modified version of "Quake.For", originally
c * written by J.M. Roesset in 1971 and modified by
c * Stavros A. Anagnostopoulos, Oct. 1986.  The formulation is that of
c * Nigam and Jennings (BSSA, v. 59, 909-922, 1969).  This modification 
c * eliminates the computation of the relative velocity and absolute 
c * acceleration; it returns only the relative displacement.  
c * Dates: 05/06/95 - Modified by David M. Boore
c *        04/15/96 - Changed name to RD_CALC and added comment lines
c *                   indicating changes needed for storing the oscillator 
c *                   time series and computing the relative velocity and 
c *                   absolute acceleration

c *   acc = acceleration time series
c *    na = length of time series
c * omega = 2*pi/per
c *  damp = fractional damping (e.g., 0.05)
c *    dt = time spacing of input
c *    rd = relative displacement of oscillator

      dimension acc(*)
      omt=omega*dt
      d2=1-damp*damp
      d2=sqrt(d2)
      bom=damp*omega
c *      d3 = 2.*bom                 ! for aa
      omd=omega*d2
      om2=omega*omega
      omdt=omd*dt
      c1=1./om2
      c2=2.*damp/(om2*omt)
      c3=c1+c2
      c4=1./(omega*omt)
      ss=sin(omdt)
      cc=cos(omdt)
      bomt=damp*omt
      ee=exp(-bomt)
      ss=ss*ee
      cc=cc*ee
      s1=ss/omd
      s2=s1*bom
      s3=s2+cc
      a11=s3
      a12=s1
      a21=-om2*s1
      a22=cc-s2
      s4=c4*(1.-s3)
      s5=s1*c4+c2
      b11=s3*c3-s5
      b12=-c2*s3+s5-c1
      b21=-s1+s4
      b22=-s4
      rd=0.
c *      rv = 0.                           ! for rv
c *      aa = 0.                           ! for aa
      n1=na-1
      y=0.
      ydot=0.
      do i=1,n1
        y1=a11*y+a12*ydot+b11*acc(i)+b12*acc(i+1)
        ydot=a21*y+a22*ydot+b21*acc(i)+b22*acc(i+1)
        y=y1    ! y is the oscillator output at time corresponding to index i
        z=abs(y)
        if (z.gt.rd) rd=z
c *        z1 = abs(ydot)                   ! for rv
c *        if (z1.gt.rv) rv = z1            ! for rv
c *        ra = -d3*ydot -om2*y1            ! for aa
c *        z2 = abs(ra)                     ! for aa
c *        if (z2.gt.aa) aa = z2            ! for aa
      end do
      return
      end
c *----------------- END RD_CALC -----------------------------

c *----------------- BEGIN AccSqInt -----------------------------
      subroutine accsqint(acc, npts, dt, rmv_trnd, a_sq_int)


c * Form integral of acceleration squared, assuming that the acceleration
c * is represented by straight lines connecting the digitized values.  This
c * routine can be used to compute Arias intensity, defined as

c *            Ixx = (pi/2g)*int(acc^2*dt), integrating from 0.0 to the total
c *  duration of the record.  The units of Ixx are 
c *  velocity [ l^(-1)t^2*(l/t^2)^2*t ] =  l^(-1+2)*t^(2-4+1) = l*t^(-1) = l/t

c * Be sure to use consistent units for the acceleration of gravity (g) and acc.
c * I am not sure what is conventionally used, but Wilson (USGS OFR 93-556) 
c * uses m/sec.

c * Dates: 01/13/99 - Written by D.M. Boore

      real a_sq_int(*), acc(*)
      logical rmv_trnd
      double precision cum, a1, a2, ddt_3

      if (rmv_trnd) then      
c * remove trend first
        call dcdt(acc, dt, npts, 1, npts, .false., .true.)
      end if

c * compute integral of squared acceleration (assume a_sq_int = 0 for first point)

      ddt_3 = dble(dt/3)

      cum = 0.0

      a_sq_int(1) = sngl(cum)
      do j=2,npts
        a1 = acc(j-1)
        a2 = acc(j)
        cum = cum + (a1**2+a1*a2+a2**2)*ddt_3
        a_sq_int(j) = sngl(cum)
      end do

c * high pass filter the velocity (may want to allow this in a future version;
c * as it is, the acceleration time series can be filtered, so there is no need
c * to do it again).

      return
      end
c *----------------- END AccSqInt -----------------------------

