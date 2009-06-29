c *----------------- BEGIN BANNER -----------------------------
      subroutine banner(nu)

c * Write version number to unit nu

      version = 2.20

      write(nu,*)
      write(nu,'(a,f6.3,a)') 
     :   '   **************** SMSIMWWW from SMSIM, Version ', 
     :    version,
     :   ' ************************'
      write(nu,*)

      return
      end
c *----------------- BEGIN BANNER -----------------------------

c *----------------- BEGIN GET_PARAMS -----------------------------
      subroutine get_params(f_params, tdflag)

c * Dates: 06/07/95 - Written by D.M. Boore
c *        10/17/95 - Added tdflag to control if more data need be read for
c *                   the time series processing; also added irmvdc as an input
c *                   parameter
c *        11/14/95 - Read eps_int from data file
c *        11/15/95 - changed input of spectral parameters
c *        11/16/95 - Read amp_cutoff from data file and determine fup
c *        11/27/95 - changed input of spectral parameters again
c *        12/06/95 - Get low-cut filter parameters before rv and td params
c *        12/07/95 - Get taper with window parameters
c *        12/19/95 - Added zup to input parameters
c *        02/28/97 - allow kappa to be magnitude dependent
c *                   ( kappa = akappa + dkappadmag*(M-amagkref) )
c *                   Use amagkref so that akappa is a reasonable value.  
c *                   One reason to do so is that in
c *                   this version of the program I determine fcut using
c *                   a magnitude-independent value for kappa (akappa), and
c *                   if the equation for kappa is given with amagkref = 0.0
c *                   the value for akappa could well be less than 0.0.
c *        01/18/99 - Added osc_crrctn to allow specification of which
c *                   correction for oscillator response to make (Boore and Joyner
c *                   or Liu and Pezeshk).
c *        01/21/99 - Skip over reading rv params if tdflag = .true.
c *        02/10/99 - Remove computation of npts for td sims; change name 
c *                   "tsimdur" to "dur_fctr"
c *        02/13/99 - Use get_lun, changed "fname" to "f", etc.
c *        08/03/99 - Added parameter (iran_type) to control type of random 
c *                   number (iran_type = 0 for normal, uniform otherwise)
c *        06/05/00 - Added velocity parameter c_q to Q parameters
c *        06/08/00 - Added r_ref to gsprd parameters and move definition of 
c *                   const to const_am0_gsprd
c *        01/27/01 - Added new time series window parameters
c *        01/27/02 - Allow geom. spreading to be magnitude dependent
c *        08/08/02 - Disabled use of remove-mean-from-noise parameter (irmvdc)

      character f_params*(*)
      logical tdflag

      include 'smsim.fi'

      pi = 4.0 * atan(1.0)

      call get_lun(nin)
      open(unit=nin, file=f_params, status='unknown')
      
c * title:
      read(nin, '(a)') title
c * DEBUG
      write(*, '(a)') title
c * DEBUG

c * rho, beta, prtitn, fs:
      call skip(nin, 1)
      read(nin, *) rho, beta, prtitn, rtp, fs
c * DEBUG
      write(*, '(a)') ' rho, beta, prtitn, rtp, fs'
      write(*, *)       rho, beta, prtitn, rtp, fs
c * DEBUG

c * source shape:
      call skip(nin, 4)
      read(nin, *) numsource, pf, pd

c * source scaling: (note: stress is now read in here; to write a driver
c * that loops over stress, just assign the desired values of stress
c * after calling get_params)
      call skip(nin, 4)
      read(nin, *) stressc, dlsdm, fbdfa, amagc

c * gsprd: 
      call skip(nin, 1)
      read(nin, *) r_ref
      read(nin, *) nsprd_segs
      do i = 1, nsprd_segs
        read(nin, *) rlow(i), a_s(i), b_s(i), m_s(i)
      end do
 
c * Q:
      call skip(nin, 1)
      read(nin, *) fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q

c * source duration:
      call skip(nin, 1)
      read(nin, *) w_fa, w_fb
c * DEBUG
      write(*, '(a)') ' w_fa, w_fb'
      write(*, *)       w_fa, w_fb
c * DEBUG

            

c * path duration:              
      call skip(nin, 1)
      read(nin, *) nknots
      do i = 1, nknots
        read(nin, *) rdur(i), dur(i)
      end do
      read(nin, *) slast

c * site amps:
      call skip(nin, 1)
      read(nin, *) namps
      do i = 1, namps
        read(nin, *) famp(i), amp(i)
      end do

c * site diminution:
      call skip(nin, 1)
      read(nin, *) fm, akappa, dkappadmag, amagkref

c * low-cut filter parameters:
      call skip(nin, 1)
      read(nin, *) fcut, norder

c * parameters for rv integration:
      if(.not. tdflag) then
        call skip(nin, 1)
        read(nin, *) zup, eps_int, amp_cutoff, osc_crrctn
        if( akappa .eq. 0.0) then
          fup = fm/amp_cutoff**0.25
        else  
          fup = 
     :     amin1(fm/amp_cutoff**0.25, -alog(amp_cutoff)/(pi*akappa))
        end if
        go to 999
      end if

c * Read more if time domain method:

c * First skip reading the rv params:
      call skip(nin, 2)

c * window parameters:
      call skip(nin, 1)
      read(nin, *) indxwind, taper, eps_w, eta_w, f_tb2te, f_te_xtnd

c * timing stuff:
      call skip(nin, 1)
      read(nin, *) dur_fctr, dt, tshift, seed, nsims, iran_type

c * Flag controlling the removal of the dc from the random series:
!      call skip(nin, 1)
!      read(nin, *) irmvdc
!
!      irmvdc = 0            ! force this value, so mean will not be removed
!

999   continue
      close(unit=nin)

      return
      end
c *----------------- END GET_PARAMS -----------------------------

c *----------------- BEGIN WRITE_PARAMS -----------------------------
      subroutine write_params(nout, tdflag)

c * Dates: 08/18/95 - Written by D. Boore
c *        10/17/95 - Added tdflag to control if more data need be read for
c *                   the time series processing; also added irmvdc as an input
c *                   parameter
c *        11/14/95 - Write eps_int
c *        11/16/95 - Write amp_cutoff, fup
c *        12/06/95 - Write low-cut filter params before rv and td params
c *        12/07/95 - Write taper with window parameters
c *        12/19/95 - Added zup to input parameters
c *        01/18/99 - Added osc_crrctn to allow specification of which
c *                   correction for oscillator response to make (Boore and Joyner
c *                   or Liu and Pezeshk).
c *        08/03/99 - Added parameter (iran_type) to control type of random 
c *                   number (iran_type = 0 for normal, uniform otherwise)
c *        06/05/00 - Added velocity parameter c_q to Q parameters
c *        01/27/01 - Added new time series window parameters
c *        01/27/02 - Allow geom. spreading to be magnitude dependent
c *        02/06/03 _ Use formatted write statements

      logical tdflag

      include 'smsim.fi'

      write( nout, '(a)') ' Title:'
      write( nout, '(4x,a)') title

      write( nout, '(a)') ' rho, beta, prtitn, rtp, fs:'
      write( nout, '(2x,f5.2, 1x,f5.2, 1x,f5.3, 1x,f4.2, 1x,f4.2)')   
     :       rho, beta, prtitn, rtp, fs

      write( nout, '(2a/2a/a/a)') 
     : ' spectral shape: source number (1=Single Corner;2=Joyner;',
     : '3=A93;4=custom)',
     : ' pf, pd (1-corner spectrum = ',
     : '1/(1+(f/fc)**pf)**pd; 0.0 otherwise)',
     : ' (usual model: pf=2.0,pd=1.0; Butterworth: pf=4.0,pd=0.5)',
     : ' (Note: power of high freq decay --> pf*pd)'
      write( nout, '(2x,i2, 1x,f4.2, 1x,f4.2)' ) numsource, pf, pd

      write( nout, '(a/a/a/a)') 
     : ' spectral scaling: stressc, dlsdm, fbdfa, amagc', 
     : ' (stress=stressc*10.0**(dlsdm*(amag-amagc))',
     : ' (fbdfa, amagc for Joyner model, usually 4.0, 7.0)',
     : ' (not used for srce 3, but placeholders still needed)'
      write( nout, '(2x,f7.2, 1x,1p,e10.3, 1x,0p,f5.2, 1x,f4.2)') 
     :     stressc, dlsdm, fbdfa, amagc

      write( nout, '(2a)') 
     :  ' gsprd: r_ref, nsegs, (rlow(i), a_s, b_s, m_s(i))',
     :                     '  (Usually set r_ref = 1.0 km)'
      write( nout, '(2x,f6.2)' ) r_ref
      write( nout, '(2x,i3)' ) nsprd_segs
      do i = 1, nsprd_segs
        write( nout, '(2x,f7.2,0p,1x,e10.3,1x,e10.3,0p,1x,f4.2)') 
     :    rlow(i), a_s(i), b_s(i), m_s(i)
      end do

      write( nout, '(a)') 
     :   ' q: fr1, Qr1, s1, ft1, ft2, fr2, qr2, s2, c_q'
      write( nout, '(2x,f7.3, 1x,f7.2, 1x,f6.3, 3(1x,f7.3),
     :     1x, f7.2, 1x, f6.3, 1x,f4.2)' ) 
     :       fr1, qr1, s1, ft1, ft2, fr2, qr2, s2, c_q

      write( nout, '(a)') ' source duration: weights of 1/fa, 1/fb'
      write( nout, '(2x,f4.2, 1x, f4.2)' ) w_fa, w_fb

      write( nout, '(2a)') ' path duration: nknots, (rdur(i), dur(i),',
     :                     ' slope of last segment'
      write( nout, '(2x,i3)' ) nknots
      do i = 1, nknots
        write( nout, '(2x,f6.1, 1x,f6.2)' ) rdur(i), dur(i)
      end do
      write( nout, '(2x,1p,e10.3)' ) slast

      write( nout, '(2a)') ' site amplification: namps, (famp(i), ',
     :                     'amp(i))'
      write( nout, '(2x,i3)' ) namps
      do i = 1, namps
        write( nout, '(2x,f7.3, 1x,f6.3)' ) famp(i), amp(i)
      end do

      write( nout, '(a)') 
     :  ' site diminution parameters: fm, akappa, dkappadmag, amagkref'
      write( nout, '(2x,f7.3, 1p, 2(1x,e10.3), 0p,1x,f4.2)' ) 
     :    fm, akappa, dkappadmag, amagkref

      write( nout, '(a)') ' low-cut filter parameters: fcut, norder'
      write( nout, '(2x,f7.3, 1x,i2)' ) fcut, norder

      if (.not. tdflag) then
        write( nout, '(a)') 
     :  ' parameters for rv calcs: zup, eps_int, amp_cutoff, osc_crrctn'
        write( nout, '(2x,f6.2, 1p, 2(1x,e10.3), 0p,1x, i2)' ) 
     :        zup, eps_int, amp_cutoff, osc_crrctn
        write(nout, '(a,1pe10.3)') ' calculated fup = ', fup
        goto 999
      end if

c * Write more if time domain method:

      write( nout, '(2a)') 
     :  ' window params: indxwind(0=box,1=exp), ',
     :  ' taper, eps_w, eta_w, f_tb2te, f_te_xtnd'
      write( nout, '(2x,i1, 1x,f4.2, 1x,f5.2, 1x,f6.3, 1x,f4.1, 
     :               1x, f4.1)' ) 
     :    indxwind, taper, eps_w, eta_w,f_tb2te,f_te_xtnd

      write( nout, '(a)') 
     :  ' timing stuff: dur_fctr, dt, tshift, seed, nsims, iran_type'
      write( nout, '(2x,f4.2, 1x,f6.4, 1x,f6.2, 1x,f9.1, 
     :               1x,i4, 1x,i1)' ) 
     :    dur_fctr, dt, tshift, seed, nsims, iran_type

!      write( nout, '(2a/2a)') 
!     :         ' parameter to control whether dc is',
!     :         ' removed from random series before',
!     :         ' transformation to frequency domain',
!     :         ' is no longer used'
c *      write( nout, *) irmvdc

999   continue
      return
      end
c *----------------- END WRITE_PARAMS -----------------------------
  
c *----------------- BEGIN SPECT_AMP -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
      function spect_amp(f)

      include 'smsim.fi'

      tempf = 1.0
      if ( idva .ne. 0 ) tempf = (twopi*f)**float(idva)
      spect_amp =   freq_indep_factor * ( tempf ) * 
     :      buttrlcf(f, fcut, norder) * 
     :      spect_shape(f, fa, fb, pf, pd, am0b_m0, numsource) * 
     :      site_amp_factor(f, namps, famp, amp) *
     :      dimin(f)
                                 ! Could save some multiplications
                                 ! by removing freq_indep_factor outside the
                                 ! function, but at the cost of possible
                                 ! confusion.

      go to ( 301, 302, 303 ), iaorins
c no instrument response.
301   h=1
      go to 300
c psv
302   v = twopi * fosc      ! converts from displacement to velocity response
      h = v * harmoscf( f, fosc, damp, idva )
      go to 300
c * for customized response (network instruments, etc)
303   h=1
      go to 300
300   continue
      spect_amp = h * spect_amp   ! spect_amp contains the spectral amplitude.

      return
      end
c *----------------- END SPECT_AMP -----------------------------

c *----------------- BEGIN CONST_AM0_GSPRD -----------------------------
c * Dates: 11/14/95 - Written by D.M. Boore
c *        07/02/99 - Added magnitude-dependent "depth" from Atkinson
c *                   and Silva, which required adding some parameters to
c *                   the passed arguments in gsprd
c *        06/08/00 - Moved computation of const into this subroutine

      subroutine const_am0_gsprd()

      include 'smsim.fi'

c * Define constant, for r=r_ref(km).Usually set r_ref = 1.0 km.  Be careful
c * if a different value or different units are used.  In particular, using 
c * different units will require a change in the scaling factor of 1.e-20 below
      const=prtitn*rtp*fs*(1.e-20)/(4.*pi*rho*beta**3*r_ref)

      freq_indep_factor = const*am0*
     :  gsprd(r, r_ref, nsprd_segs, rlow, a_s, b_s, m_s, 
     :        numsource, amag)
c *                         (const from Get_Params, am0 from Spect_Scale) 

      return
      end
c *----------------- END CONST_AM0_GSPRD -----------------------------

c *----------------- BEGIN GSPRD -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
c *        07/02/99 - Added magnitude-dependent "depth" from Atkinson
c *                   and Silva, which required adding some parameters to
c *                   the passed arguments
c *        06/05/00 - Added some explanation of r
c *        06/08/00 - Make gsprd nondimensional through the use of r_ref, which 
c *                   now appears in the definition of variable const
c *                   in const_am0_gsprd
c *        01/27/02 - Following Silva, parameters added to allow magnitude
c *                   dependent slope (to capture finite fault effects)

      function gsprd(r,r_ref,nsprd_segs,rlow,a_s,b_s,m_s,
     :               numsource,amag)
      real r_ref, rlow(*), a_s(*), b_s(*), m_s(*), geff(10)
c * Note that generally r = hypocentral distance.  For Atkinson and Silva 
c * (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
c * their paper; below their eq. 4), so that rmod is, in effect, accounting
c * source depth twice 
      rmod = r
      if (numsource .eq. 9) then            ! Atkinson and Silva (2000)
        deff = 10.0**(-0.05 + 0.15 * amag)
        rmod = sqrt(r**2 + deff**2)
      end if
      geff(1) = r_ref/rlow(1)  ! usually set r_ref = 1.0 km.  Be careful
                               ! if a different value or different units are
                               ! used.  In particular, using different units
                               ! will require a change in the scaling factor
                               ! of 1.e-20 used in the definition of const in
                               ! const_am0_gsprd

      do i = 2, nsprd_segs
        slope = a_s(i-1) + b_s(i-1)*(amag - m_s(i-1))
        geff(i) = geff(i-1)*(rlow(i)/rlow(i-1))**slope
      end do
      if (rmod .le. rlow(1)) then
        j = 1
      else if (rmod .ge. rlow(nsprd_segs)) then
        j = nsprd_segs
      else
        call locate(rlow, nsprd_segs, rmod, j)
      end if
      slope = a_s(j) + b_s(j)*(amag - m_s(j))

      gsprd = (geff(j)) * (rmod/rlow(j))**slope

      return
      end
c *----------------- END GSPRD -----------------------------

c *  ------------------- BEGIN BUTTRLCF -------------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
      function buttrlcf( f, fcut, norder)
c
c Computes the response of an norder, bidirectional
c * high-pass Butterworth filter.  This is the filter
c * used by the AGRAM processing (the equation was
c * provided by April Converse).

c * Modification: 3/27/89 - created by modifying HiPassF

      real buttrlcf
      buttrlcf = 1.0
      if ( fcut.eq.0.0 ) return

      buttrlcf = 0.0

      if ( f .eq. 0.0) return

      buttrlcf = 1.0/ (1.0+(fcut/f)**(2.0*norder))

c * Note: Because this filter is intended to simulate the
c * two-pass, zero-phase (acausal) Butterworth filter used in
c * accelerogram processing, buttrlcf = 1/2 when f = fcut, not 1/sqrt(2) as in
c * a single-pass, causal Butterworth filter.

      return
      end
c *  ------------------- END BUTTRLCF -------------------------------

c *----------------- BEGIN SPECT_SHAPE -----------------------------
c * Source Displacement Spectrum
c * Dates: 06/07/95 - Written by D.M. Boore
c *        11/15/95 - changed source types
c *        12/02/96 - added Atkinson and Silva (model 4) (I did this earlier)
c *                   and Haddon (model 5)
c *        02/28/97 - Added Atkinson's new version of the source spectra
c *                   derived from Atkinson and Silva (this will appear
c *                   in Atkinson and Boore... Evaluation of source models...).
c *                   (model 6)
c *        06/24/97 - added Boatwright and Choy scaling (model 7).
c *        07/21/97 - Added Joyner ENA model (model 8; the spectral shape is
c *                   the same as his model 2, but I because the corner frequency
c *                   relations are new I have to repeat the shape with the new 
c *                   model number).
c *        09/02/98 - Renamed model 6 to AB98-Ca to be consistent with usage
c *                   in Tables 3 and 4 in Atkinson and Boore, BSSA, v. 88, 
c *                   p. 917--934.
c *        02/16/99 - Added Atkinson and Silva, 1999, as model 9
c *        06/05/00 - Changed "AS99" to "AS2000" because the paper was published
c *                   in 2000 (BSSA 90, 255--274)

      function spect_shape(f, fa, fb, pf, pd, am0b_m0, numsource)
      real spect_shape
      goto (1, 2, 3, 4, 5, 6, 7, 8, 9, 10), numsource

      write(*, '(a, i5, a)') ' !!!!!! numsource = ',
     : numsource, ', .ne. legal value; quitting !!!!!!'
      stop
 
c * Single corner frequency:
1     sb = 1.0
      sa = 1.0/( 1.0 + (f/fa)**pf )**pd
      go to 900

c * Joyner model
2     sb = 1.0/ ( 1.0 + (f/fb)**2 )**0.25
      sa = 1.0/ ( 1.0 + (f/fa)**2 )**0.75
      go to 900

c * Atkinson 1993 model
3     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

c * Atkinson & Silva 1996 (same format as Atkinson 1993) model
4     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

c * Haddon (see 12/02/96 notes; approximate spectra in Fig. 10 of
c * Haddon's paper in BSSA, v. 86, p. 1312)
5     pda = 1.0/8.0
      pdb = 1.0/8.0
      pfa = 1/pda
      pfb = 1/pdb
      sa = 1.0/( 1.0 + (f/fa)**pfa )**pda
      sb = 1.0/( 1.0 + (f/fb)**pfb )**pdb
      go to 900

c * AB98-Ca (Atkinson & Boore 1998) (same format as Atkinson 1993) model
6     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

c * Boatwright and Choy (this is the functional form used by 
c *  Boore and Atkinson, BSSA 79, p. 1761)
7     sa = 1.0
      if (f .ge. fa) sa = fa/f
      sb = 1.0/sqrt( 1.0 + (f/fb)**2 )
      go to 900 

c * Joyner model (to be used with his ENA two-corner model)
8     sb = 1.0/ ( 1.0 + (f/fb)**2 )**0.25
      sa = 1.0/ ( 1.0 + (f/fa)**2 )**0.75
      go to 900

c * 
c * AS2000 (Atkinson and Silva, 2000, BSSA 90, 255--274) 
c * (same format as Atkinson 1993) model
9     sb = 1.0
      sa = (1.0 - am0b_m0)/( 1.0 + (f/fa)**2 )
     :      +      (am0b_m0)/( 1.0 + (f/fb)**2 )
      go to 900

c * For customized relation:
10    sb = 1.0
      sa = 1.0
      go to 900

 
900   continue
      spect_shape = sa*sb

      return
      end
c *----------------- END SPECT_SHAPE -----------------------------

c *----------------- BEGIN SPECT_SCALE -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
c *        06/05/96 - Added modified Atkinson & Silva scaling
c *        12/02/96 - added Haddon scaling (see spect_shape.for comments)
c *        02/28/97 - added Atkinson and Boore scaling 
c *                   (see spect_shape.for comments)
c *        06/24/97 - added Boatwright and Choy scaling
c *        07/10/97 - changed A93, AS96, AS97 scaling to be constant
c *                   stress for M < Mc, where Mc is the magnitude for which
c *                   am0b_m0 = 1.0.  The single corner model for smaller
c *                   magnitudes is determined so that the high frequency
c *                   level is matched at M = Mc.
c *        07/21/97 - Added Joyner 2-corner model for ENA, as specified 
c *                   in his notes prepared for the SSHAC workshop (published 
c *                   in SSHAC report, NUREG/CR-6372, p. B-303--B-305, 1997).
c *        08/06/97 - I discovered that Joyner erroneously fit vertical spectral
c *                   amplitudes.  He provided a revised model, fitting the
c *                   horizontal amplitudes.  I changed model 8 accordingly.
c *        09/02/98 - Renamed model 6 to AB98-Ca to be consistent with usage
c *                   in Tables 3 and 4 in Atkinson and Boore, BSSA, v. 88, 
c *                   p. 917--934.
c *        02/16/99 - Added Atkinson and Silva (1999)
c *        06/05/00 - Changed "AS99" to "AS2000" because the paper was published
c *                   in 2000 (BSSA 90, 255--274)

      subroutine spect_scale()

c *      include 'smsim.fi'
      include 'smsim.fi'

      am0 = 10.**(1.5*amag + 16.05)
      am0b_m0 = 0.0

      goto (1, 2, 3, 4, 5, 6, 7, 8, 9, 10), numsource

      write(*, '(a, i5, a)') ' !!!!!! numsource = ',
     : numsource, ', .ne. legal value; quitting !!!!!!'
      stop
 
c * Single corner frequency:
1     stress = stressc*10.0**(dlsdm*(amag-amagc))
      fa = (4.906e+06) * beta * (stress/am0)**(1.0/3.0)
      fb = fa
      return

c * Joyner scaling:
2     am0c = 10.0 ** ( 1.5*amagc + 16.05 )
      stress = stressc*10.0**(dlsdm*(amag-amagc))
      rat = stress/am0
      dum = 4.906e+06
      if ( am0 .gt. am0c ) rat = stress/am0c
      fb = ( dum*beta ) * ( fbdfa )**(3./4.) * ( rat )**(1./3.)
      fa = ( dum*beta ) * (stress)**(1./3.) * (am0c)**(1./6.)
     *     * ( fbdfa )**(-0.25) * ( am0 )**(-0.5)
      if ( am0 .lt. am0c ) fa = fb / fbdfa
      return

c * Atkinson 93 scaling:
3     if (amag .gt. 4.0) then
        fa = 10.0**(2.41 - 0.533 * amag)
        fb = 10.0**(1.43 - 0.188 * amag)      ! fa = fb for M = 2.84
        am0b_m0 = 10.0**(2.52 - 0.637 * amag)
      else
        fb = 10.0**(2.678 - 0.5 * amag)
        fa = fb
        am0b_m0 = 1.0
      endif
      return

c * Atkinson and Silva 96 scaling, with am0b_m0 modified by D. Boore on 6/04/96
4     if (amag .gt. 4.6) then
        fa = 10.0**(2.181 - 0.496 * amag)
        fb = 10.0**(1.778 - 0.302 * amag)   ! fa = fb for M = 2.08
        am0b_m0 = 10.0**(3.440 - 0.746 * amag)  ! DMB's fitting of spctrl ratios
c *        am0b_m0 = 10.0**(2.764 - 0.623 * amag) ! in Atkinson & Silva preprint
      else
        fb = 10.0**(2.689 - 0.5 * amag)
        fa = fb
        am0b_m0 = 1.0
      endif
      return

c * Haddon scaling:
5     fa = 10.0**(0.3 - (1.5/3)*(amag-4.0))
      fb = 10.0**(1.4 - (1.5/3)*(amag-4.0))  ! fa < fb for all M
      return

c * AB98-Ca (Atkinson and Boore 98 scaling, based on fitting a model to the 
c * Atkinson and Silva 1997 Fourier amplitude spectra for California; see
c * Atkinson and Boore, BSSA, v. 88, p. 917--934).
6     if (amag .gt. 4.8) then
        fa = 10.0**(2.181 - 0.496 * amag)
        fb = 10.0**(1.308 - 0.227 * amag)      ! fa=fb for M = 3.25
        am0b_m0 = 10.0**(3.223 - 0.670 * amag)
      else
        fb = 10.0**(2.617 - 0.5 * amag)
        fa = fb
        am0b_m0 = 1.0
      endif
      return

c * Boatwright and Choy (this is not from Boore and Atkinson, BSSA 79, p. 1761;
c *  it is based on new fits by DMB on 6/24/97 to data in Boat.&Choy, 1992 BSSA.
c *  See BC_BSSA.WQ1 in \haddon subdirectory and handwritten notes on 
c *  yellow sheets.
c *  except set fa=fb=constant stress scaling for M<=5.3)
7     fa = 10.0**(3.409 - 0.681 * amag)
      fb = 10.0**(1.495 - 0.319 * amag)
      if (amag .le. 5.3) then
        fa = 0.634*10.0**(0.5*(5.3 - amag)) ! 0.634 = 10^(logfa+logfb)/2 at M5.3
        fb = fa
      end if
      return

c * Joyner ENA scaling:
8     fa = 10.0**(2.312 - 0.5 * amag)
      fb = 10.0**(3.609 - 0.5 * amag)
      return

c * AS2000 -- Atkinson and Silva, 2000 scaling, based on fitting a point source
c * model to finite source calculations, with constraints on various modeling
c * studies, with modification for very small magnitude (when eps = 1) 
9     if (amag .gt. 2.4) then
        fa = 10.0**(2.181 - 0.496 * amag)
        fb = 10.0**(2.41  - 0.408 * amag)      ! fa=fb for M = -2.6
        am0b_m0 = 10.0**(0.605 - 0.255 * amag)
      else
        fb = 10.0**(1.431 - 0.5 * (amag - 2.4))
        fa = fb
        am0b_m0 = 1.0
      endif
      return

c * For customized scaling:
10    fa = 0.0
      fb = 0.0
      return

      end
c *----------------- END SPECT_SCALE -----------------------------

c *----------------- BEGIN SITE_AMP_FACTOR -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
      function site_amp_factor(f, namps, famp, amp)
      real famp(*), amp(*), site_amp_factor
      
      if ( f .le. famp(1) ) then
        site_amp_factor = amp(1)
      else if ( f .ge. famp(namps) ) then
        site_amp_factor = amp(namps)
      else
        call locate( famp, namps, f, j)
        site_amp_factor = amp(j)*10.0**(alog10(f/famp(j))
     :             *alog10(amp(j+1)/amp(j))
     :             /alog10(famp(j+1)/famp(j)))
       end if

       return
       end
c *----------------- END SITE_AMP_FACTOR -----------------------------

c *----------------- BEGIN DIMIN -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
c *        07/02/99 - Added modification to r required by Atkinson
c *                   and Silva (1999)
c *        06/05/00 - Substitute c_q for beta in akappaq and add comments
c *                   about r

      function dimin(f)
      real dimin, mag
      include 'smsim.fi'

c * Note that generally r = hypocentral distance.  For Atkinson and Silva 
c * (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
c * their paper; below their eq. 4), so that rmod is, in effect, accounting
c * source depth twice 
      rmod = r
      if (numsource .eq. 9) then            ! Atkinson and Silva (2000)
        deff = 10.0**(-0.05 + 0.15 * amag)  
        rmod = sqrt(r**2 + deff**2)
      end if

      akappaq = rmod/(c_q*q(f))

      mag = amag    
      dimin = exp( -pi*(kappa_f(mag) + akappaq) * f)/
     :   sqrt( 1. + (f/fm)**8.0 )
      

      return
      end
c *----------------- END DIMIN -----------------------------

c *----------------- BEGIN KAPPA_F -----------------------------
c * Dates: 02/28/97 - Written by D.M. Boore
      function kappa_f(mag)
      real mag
      include 'smsim.fi'
  
      kappa_f = akappa + dkappadmag*(mag-amagkref)

      return
      end
c *----------------- END KAPPA_F -----------------------------
        
c *----------------- BEGIN Q -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
c *        12/14/95 - Added check for equal transition frequencies
      function q(f) 
      logical firstcall / .true. /
      save qt1, qt2, st
      include 'smsim.fi'

      q = 9999.0
      if (f .eq. 0.0) return
        
      if (firstcall) then
        qt1 = qr1*(ft1/fr1)**s1
        qt2 = qr2*(ft2/fr2)**s2
        st = 0.0
        if (ft1 .ne. ft2) then
          st = alog10(qt2/qt1)/alog10(ft2/ft1)
        end if
        firstcall = .false.
      end if
      
      if ( f .le. ft1) then
        q = qr1*(f/fr1)**s1
      else if ( f .ge. ft2) then
        q = qr2*(f/fr2)**s2
      else
        q = qt1*(f/ft1)**st
      end if

      return
      end
c *----------------- END Q -----------------------------

c *----------------- BEGIN DURSOURCE -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
      function dursource(w_fa, fa, w_fb, fb)
      real dursource
      dursource = w_fa/fa + w_fb/fb
      return
      end
c *----------------- END DURSOURCE -----------------------------
      
c *----------------- BEGIN DURPATH -----------------------------
c * Dates: 06/07/95 - Written by D.M. Boore
c *        07/02/99 - Added magnitude-dependent "depth" from Atkinson
c *                   and Silva, which required adding some parameters to
c *                   the passed arguments of durpath
c *        06/05/00 - Changed reference for Atkinson and Silva (1999) to
c *                   2000, the year in which their paper was published
c *                   and added comments regarding r

      function durpath(r, nknots, rdur, dur, slast, numsource, amag)
      real rdur(*), dur(*), durpath


c * Note that generally r = hypocentral distance.  For Atkinson and Silva 
c * (BSSA 90, 255--274) r is the closest distance to the fault plane ("d" in 
c * their paper; below their eq. 4), so that rmod is, in effect, accounting
c * source depth twice 
      rmod = r
      if (numsource .eq. 9) then            ! Atkinson and Silva (2000)
        deff = 10.0**(-0.05 + 0.15 * amag)
        rmod = sqrt(r**2 + deff**2)
      end if

      if ( rmod .ge. rdur(nknots) ) then
        durpath = (rmod -rdur(nknots))*slast + dur(nknots)
      else
        call locate(rdur, nknots, rmod, j)
        durpath = (rmod - rdur(j))*( dur(j+1)-dur(j))
     :            /(rdur(j+1)-rdur(j)) + dur(j)
      end if

      return
      end
c *----------------- END DURPATH -----------------------------

c *----------------- BEGIN HARMOSCF -----------------------------
      function harmoscf( f, fosc, damp, idva)
c harmonic oscillator displacement response. 
c * idva = 0 for response to displacement
c * idva = 2 for response to acceleration
c * idva = 1 returns 0 for the response
c * The response is normalized to be unity in the flat portion.

c * Written by D. Boore, 12/01/83
c latest modification:  3/25/89
c *                       7/30/00 - Changed dum for both cases of idva
c *                                 (see notes from same date).

      pi = 4.0*atan(1.0)
      twopi = 2.0 * pi

      if (idva .ne. 0 .and. idva .ne. 2) dum = 0.0
      if (idva .eq. 0) dum = (twopi*f)**2/twopi**2
      if (idva .eq. 2) dum =          1.0/twopi**2

      harmoscf = dum/sqrt( ( fosc*fosc - f*f )**2
     * + ( 2.0*f*fosc*damp )**2 )
      return
      end
c *----------------- END HARMOSCF -----------------------------

c *----------------- BEGIN SMSIMFAS -----------------------------
      subroutineas(fas, freq, nfreq)
c * Dates: 12/16/95 - Written by D.M. Boore, based on SMSIM_RV
c *        03/05/99 - Removed include statements from the end (put them into the
c *                   driver programs)

      real fas(*), freq(*)
      include 'smsim.fi'

      pi = 4.0*atan(1.0)
      twopi = 2.0 * pi

c * Set spectral parameters:
      call spect_scale()   ! Be sure to call Get_Params in driver

c * Get frequency-independent factor:
      call const_am0_gsprd()

      do i = 1, nfreq
        fas(i) = spect_amp(freq(i))
      end do

      return
      end
c *----------------- END SMSIMFAS -----------------------------

c * --------------------------- BEGIN BJFV1V2F -------------------          
      function BJFV1V2F(per, V1, V2)

c * Returns the correction factor to apply to response spectra with V30 = V1
c * to obtain the response spectra for V30 = V2 (i.e., computes PRV(V2)/PRV(V1)).
c * The value for per = 0.1 is used for all smaller periods and the per=2.0 value 
c * is used for periods longer than 2.0 sec (in an earlier version I used
c * per = 1.0).  Note that the latter is
c * conservative; we expect the amplifications to reach unity for long enough
c * periods (at least for Fourier spectra, what about for response spectra?)

c * Dates: 07/17/96 - Written by D. Boore, based on BJFR2S_F.FOR
c *        07/18/96 - Changed endpoint period from 1.0 to 2.0
c * Dates: 07/24/96 - Added check of V1, V2 equal zero, in which case unit
c *                   amplification is returned.
c *                   Also, reversed the meaning of V1 and V2.
c *                   Added pga amps when per = 0.01
c *        10/08/00 - pga when per = 0.0, rather than 0.01

      if (v1 .eq. 0.0 .or. v2 .eq. 0.0) then
         bjfv1v2f = 1.0
         return
      end if

      velfact = alog10( V2/V1 )
      if(per .lt. 0.1) then
         if(per .eq. 0.0) then
           bjfv1v2f = 10.0**(velfact*(-0.371))
         else           
           bjfv1v2f = 10.0**(velfact*cubic(0.1))
         end if
      else if (per .gt. 2.0) then
         bjfv1v2f = 10.0**(velfact*cubic(2.0))
      else
         bjfv1v2f = 10.0**(velfact*cubic(per))
      end if

      return
      end

      function cubic(per)
        c0 = -0.21172
        c1 =  0.06619
        c2 = -1.35085
        c3 =  0.79809
        x = alog10(per/0.1)
        cubic = c0 + c1*x + c2*x**2 + c3*x**3
c *        a0 = 0.2102
c *        a1 = 0.0726
c *        a2 = -0.3142
c *        a3 = -0.2403
c *        x = alog10(per)
c *        cubic = a0 + a1*x + a2*x**2 + a3*x**3
      return
      end
c * --------------------------- END BJFV1V2F -------------------          

c *----------------- BEGIN SKIP -----------------------------
      subroutine SKIP(lunit, nlines)
        do i = 1, nlines
           read(lunit, *)
        end do
        return
      end
c *----------------- END SKIP -----------------------------

c * --------------------- BEGIN UPSTR ----------------------------------
      Subroutine UPSTR ( text )
c * Converts character string in TEXT to uppercase
c * Dates: 03/12/96 - Written by Larry Baker

C
      Implicit   None
C
      Character  text*(*)
C
      Integer    j
      Character  ch
C
      Do 1000 j = 1,LEN(text)
         ch = text(j:j)
         If ( LGE(ch,'a') .and. LLE(ch,'z') ) Then
            text(j:j) = CHAR ( ICHAR(ch) - ICHAR('a') + ICHAR('A') )
         End If
 1000    Continue
C
      Return
      End
c * --------------------- END UPSTR ----------------------------------

c * --------------------------- BEGIN GET_LUN ----------------
      subroutine get_lun(lun)

c * Finds a logical unit number not in use; returns
c * -1 if it cannot find one.

c * Dates -- 05/19/98 - Written by D. Boore, following
c *                     Larry Baker's suggestion

      logical isopen
      do i = 99,10,-1
        inquire (unit=i, opened=isopen)
        if(.not.isopen) then
          lun = i
          return
        end if
      end do
      lun = -1

      return
      end
c * --------------------------- END GET_LUN ----------------

c * --------------------------- BEGIN TRIM_C -----------------------
      subroutine trim_c(cstr, nchar)

c * strips leading and trailing blanks from cstr, returning the
c * result in cstr, which is now nchar characters in length

c * Here is a sample use in constructing a column header, filled out with 
c * periods:

c ** Read idtag:
c *        idtag = ' '
c *        read(nu_in, '(1x,a)') idtag
c *        call trim_c(idtag, nc_id)
c ** Set up the column headings:
c *        colhead = ' '
c *        colhead = idtag(1:nc_id)//'......' ! nc_id + 6 > length of colhead

c * Dates: 12/23/97 - written by D. Boore
c *        12/08/00 - pad with trailing blanks.  Otherwise some comparisons
c *                   of the trimmed character string with other strings
c *                   can be in error because the original characters are left
c *                   behind after shifting.  For example, here is a string
c *                   before and after shifting, using the old version:
c *                      col:12345
c *                           MTWH  before
c *                          MTWHH  after (but nc = 4).
c *        03/21/01 - Check for a zero length input string
c *        11/09/01 - Change check for zero length string to check for all blanks

      character cstr*(*)

      if(cstr .eq. ' ') then
        nchar = 0
        return
      end if

      nend = len(cstr)

c *      if(nend .eq. 0) then
c *        nchar = 0
c *        return
c *      end if

      do i = nend, 1, -1
        if (cstr(i:i) .ne. ' ') then
           nchar2 = i
           goto 10
        end if
      end do

10    continue

      do j = 1, nchar2
        if (cstr(j:j) .ne. ' ') then
          nchar1 = j
          goto 20
        end if
      end do

20    continue
   
      nchar = nchar2 - nchar1 + 1
      cstr(1:nchar) = cstr(nchar1: nchar2)
      if (nchar .lt. nend) then
        do i = nchar+1, nend
          cstr(i:i) = ' '
        end do
      end if

      return
      end
c * --------------------------- END TRIM_C -----------------------

c * ------------------------------------------------------------------ skipcmnt
      subroutine skipcmnt(nu, comment, ncomments)

c * Skip text comments in the file attached to unit nu, but save skipped 
c * comments in character array comment.  Skip at least one line, and more as 
c * long as the lines are preceded by "|" or "!".

c * Dates: 04/16/01 - Written by D. Boore
c *        12/07/01 - Added check for eof

      character comment(*)*(*), buf*80

      ncomments = 0
100   buf = ' '
      read (nu,'(a)',end=999) buf
      if (buf(1:1) .eq.'!' .or. buf(1:1) .eq.'|' .or. 
     :                     ncomments + 1 .eq. 1) then
        ncomments = ncomments + 1
        comment(ncomments) = buf
        goto 100
      else 
        backspace nu
      end if

999   continue

      return
      end
c * ------------------------------------------------------------------ skipcmnt


     
