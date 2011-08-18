c----program hazinterp.nga.f
c--- interpolates hazard curve file to 0.5 original spacing in lat and lon
c writes new sampled dx,dy values into header record. See extra() in header.
c Can do this for several files. Put instructions in script and run as:
c      hazinterp.nga < script.of.instructions
c line 1 : number of files
c line 2 : vs30 (m/s). this was needed for frankel files but not for more recent
c line 3 : first infile
c line 4 : first outfile name
c line 5 : 2nd infile
c ...
c
c Compile Solaris:
c f95 hazinterp.nga.f iosubs.o -o hazinterp.nga -e
c -e to extend line length beyond the col. 72 default
c Gfortran:
c  this code has less testing on Windows platforms.
c Steve Harmsen May 16 2006. Minor revs to increase array sizes mar 15 2006.
c modified from A. Frankel program hazinterpnga.f 
c This program assumes that the geographic region information is contained in the
c binary header and that this information supercedes anything you may find elsewhere.
c      
      parameter (nin=560,nout=1102)
      type header 
        character*128 :: name(6)
        real*4 :: period
        integer*4 :: nlev
        real*4 :: xlev(20)
        real*4 :: extra(10)
      end type header
      type(header) :: headr,head
      integer*4 readn
c As usual, "prob" files are in general mean annual freq. of exceedance not probability
      dimension prob(1),xlev(20)
      character name*60,namout*60
      character namepar*30,name2*30,pgm*30
      real pout(nout,nout,20),mat(nin,nin,20)
      pointer (p1,prob)
      logical ishere,empty_hd
      call getarg(0,pgm)
      print 5,'Enter number of files to resample: '
      read *,nf
       print *,'enter vs30 m/s for the data or your best guess'
       read *,vs30
       head%extra(9)=vs30
      do ijk=1,nf
      write(6,5) "enter name of input hazard-curve file: "
      read(5,900,end=299) name
      write(6,900)name
5      format(a,$)
 900  format(a)
      write(6,5) "Enter name of output file: "
      read(5,900) namout
      write(6,900)namout
      inquire (file=namout,exist=ishere)
      if(ishere)print *,'Warning : this file already exists in WD'
      call openr(name)
c 30 character name 
c      ndata= 308
c 128 character name 
      ndata= 896
      call gethead(head,ndata,readn)
c------ get map parameters
c      write(6,*) "enter name for par file for new header"
c      read(5,900) name
         name=head%name(1)
         head%name(3)=pgm      !put program name into this spot for tracking
         print *,'levels: ',head%xlev
         empty_hd=head%extra(6).eq.0..and.head%extra(3).eq.0
         if(.not.empty_hd)then
c get information from header if at all possible.
         ymin=head%extra(5)
         ymax=head%extra(6)
         dy=head%extra(7)
         xmin=head%extra(2)
         xmax=head%extra(3)
         dx=head%extra(4)
         else
         open(unit=1,file=head%name(1),status='old')
         read(1,*) idum
         read(1,*) ymin,ymax,dy
        read(1,*) xmin,xmax,dx
       head%extra(2)=xmin
       head%extra(3)=xmax
       head%extra(5)=ymin
       head%extra(6)=ymax
        close(1)
        endif
        write(6,*) ymin,ymax,dy,xmin,xmax,dx,'  geographic sampling '
c make a  robust calculation of nx and ny using "nint" function. SH.
        nx= nint((xmax-xmin)/dx) +1
        ny= nint((ymax-ymin)/dy) +1
        if(nx.gt.nin.or.ny.gt.nin)then
        print *,'coarsemesh grid > nin. Increase mmat dims please'
        stop'hazinterp.nga.f could not continue'
        endif
        write(6,*) nx,ny
        nx2= 2*nx-1
        ny2= 2*ny-1
        if(nx2.gt.nout.or.ny2.gt.nout)then
        print *,'finemesh grid count > 1002. Please up dims of pout()'
        stop 'hazinterp.nga.f needs to be modified'
        else
        print *,nx2,ny2,' dims of fine output array'
        endif
        nrec= nx*ny
cccc
      nlev= head%nlev
      ndata=nlev*nrec
      write(6,*) "nlev=",nlev, " ndatain ",ndata
      write(6,*) 'resampling ...'
      ymin0=ymin
      p1 = malloc(ndata*4)
      do 300 i=1,nlev
300   xlev(i)= head%xlev(i)
      call getbuf2(prob,ndata,readn)
      if(ndata.ne.readn)write(6,*) ndata,readn,' mismatch. Bad.'
      izero=0
      do 10 iy=1,ny
      do 10 ix=1,nx
      j0= ix-1+ (iy-1)*nx
      j0= j0*nlev
      do 11 j=1,nlev
      pj=prob(j+j0)
      if(pj.lt.1e-15)izero=izero+1
      mat(ix,iy,j)=pj
  11  continue
  10  continue
        frac0=float(izero)/float(nlev+j0)
        
        print *,'fraction of zero data input ',frac0,' izero ',izero
c---- interpolate  matrix to half spacing
      do 105 ilev=1,nlev
c---- place existing values into new matrix
      do 1 iy=1,ny
      iy2= 2*iy -1
      do 1 ix=1,nx
      ix2= 2*ix -1
c      if(ilev.eq.1.and.mat(ix,iy,1).gt.1e-7)print *,mat(ix,iy,1)
   1  pout(ix2,iy2,ilev)= mat(ix,iy,ilev)
c-- fill in even x-values for odd rows
      do 2 iy=1,ny2,2
      do 2 ix=2,nx2-1,2
   2  pout(ix,iy,ilev)= 0.5*(pout(ix-1,iy,ilev)+pout(ix+1,iy,ilev))
c-- fill in even rows
      do 3 iy=2,ny2-1,2
      do 3 ix=1,nx2
  3   pout(ix,iy,ilev)= 0.5*(pout(ix,iy-1,ilev)+pout(ix,iy+1,ilev))
 105  continue
      call close(name)
      nrec= nx2*ny2
      call free(p1)
      nbyte=nrec*nlev*4
      p1=malloc(nbyte)
        head%extra(4)=dx*0.5      !save the change in sampling
        head%extra(7)=dy*0.5      !to header record
        head%extra(8)=float(nrec)
      do 50 i=1,nrec
      do 50 j=1,nlev
      index= (i-1)*nlev +j
      iy= (i-1)/nx2 +1
      ix= i-(iy-1)*nx2
      prob(index) = pout(ix,iy,j)
 50   continue
      write(6,*) namout
c----
      call openwx(ifp,namout)
c---- output probabilistic ground motions
c 30 character name 
c      ndata=308
c 128 character name 
      ndata= 896
      call puthead(ifp,head,ndata,readn)
      ndata= nrec*nlev
      call putbufx(ifp,prob,ndata,readn)
      if(ndata.ne.readn)write(6,*) ndata,readn

       enddo      !nfiles to be modfified
       stop 'normal exit'
 299      print *,'Your control file did not have enough file names ',nf,ijk
       stop 'hazinterp.nga: however, the ones you listed were resampled.'
      end
