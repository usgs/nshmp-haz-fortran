*--------------------- START GET_DATE  ---------------------
      subroutine get_date(date_c)

      character date_c*(*)
c      integer*2 iyr, imon, iday

* Microsoft compiler version:
*      call GETDAT (iyr, imon, iday)

*      date_c = ' '         
C	 character date_c*10

*      date_c(3:3) = '/'
*      date_c(6:6) = '/'

*      write(date_c(1:2),'(i2.2)') imon      
*      write(date_c(4:5),'(i2.2)') iday      
*      write(date_c(7:10),'(i4.4)') iyr      

* Lahey compiler version:
*      call DATE(date_c)   
C	 character date_c*8, but *10 is OK

*	UNIX SUNOS 1.?
	integer iarray(3)
	call idate(iarray)
      date_c = ' '         
      date_c(3:3) = '/'
      date_c(6:6) = '/'

      write(date_c(1:2),'(i2.2)') iarray(1)      
      write(date_c(4:5),'(i2.2)') iarray(2)      
      write(date_c(7:10),'(i4.4)') iarray(3)      


      return
      end
*--------------------- END GET_DATE  ---------------------

*--------------------- START GET_TIME  ---------------------
      subroutine get_time(time_c)

      character time_c*(*)
      integer*2 ihr, imin, isec, i100th

* Microsoft compiler version:
*      call GETTIM (ihr, imin, isec, i100th)

*      time_c = ' '             
C	 character time_c*11 for both compilers

*      time_c(3:3) = ':'
*      time_c(6:6) = ':'
*      time_c(9:9) = '.'

*      write(time_c(1:2),'(i2.2)') ihr      
*      write(time_c(4:5),'(i2.2)') imin      
*      write(time_c(7:8),'(i2.2)') isec      
*      write(time_c(10:11),'(i2.2)') i100th

* Lahey compiler version:
*      call TIME(time_c)
*	UNIX SUNOS 1.?
	integer iarray(3)
	call itime(iarray)
      time_c = ' '             
      time_c(3:3) = ':'
      time_c(6:6) = ':'
      time_c(9:9) = '.'

	i100th = 0
      write(time_c(1:2),'(i2.2)') iarray(1)      
      write(time_c(4:5),'(i2.2)') iarray(2)      
      write(time_c(7:8),'(i2.2)') iarray(3)      
      write(time_c(10:11),'(i2.2)') i100th

      return
      end
*--------------------- END GET_TIME  ---------------------
* --------------- BEGIN TIME_DIFF ---------------------------------
* Dates: 06/07/95 - Written by D.M. Boore
      subroutine time_diff(time_start, time_stop, time_elapsed)
      character time_start*(*), time_stop*(*)
      read(time_start(1:11),'(i2,1x,i2,1x,i2,1x,i2)') 
     :                       ihb, imb, isb, ihsb
      read(time_stop(1:11),'(i2,1x,i2,1x,i2,1x,i2)') 
     :                       ihe, ime, ise, ihse
      tbegin = 3600.0*ihb + 60.0*imb + isb + float(ihsb)/100.0
      tend   = 3600.0*ihe + 60.0*ime + ise + float(ihse)/100.0
      time_elapsed = tend - tbegin
      return
      end
* --------------- END TIME_DIFF ---------------------------------

