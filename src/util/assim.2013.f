c assim.2013.f
c assimilates and distills information in hazFXnga13l deterministic output files (.DET)
c works for namx = number of atten. models =  3 to (5 nga-w models+GK) to 9
c CEUS mod, new Feb 19 2013: incl ab08' a06' p11 or whatever.
c Cascadia mod, new Mar 4 2013: 4 models incl BCHydro 
c f95 compile:
c f95 assim.2013.f -o ../bin/assim.2013 -e
c Steve Harmsen May 2013.
        parameter (namax=9)
        real x,y
        real, dimension(namax) :: wte
        integer, dimension(namax) :: indexe
        real, dimension(3) :: wtc
        real, dimension(4) :: wtc2
        real, dimension(6):: wtw
        integer, dimension(6):: indw
        integer, dimension(3) :: indexc
        integer, dimension(4) :: indc2
c wtc indexc are weights and indexes of Cascadia attenuation models. Code mod mar 4 2008
c
        character*36 fname(1000)
        logical written/.false./,usewt,usewtc,subduction,more,wus13,namef, ceus
c sa is the median sa for a given atten model and sd is the std dev. for that model. Assume namx models
        real, dimension (namax)  :: sa,samax,sd,sdmax
        character*88 rec,recold,fi*88,fi2
        indexe = (/6,2,20,7,10,19,25,26,27/)
        indexc = (/2,5,7/)
        indw = (/33,34,35,36,37,38/)
        wtw = (/0.18,0.18,0.18,0.18,0.18,0.1/)
        indc2 = (/5,7,20,21/)
        wtc2 = (/0.25,0.25,0.25,0.25/)
        wtc = (/0.25, 0.25, 0.5/)
c        wte = (/0.09,0.09,0.09,0.1,0.09,0.09,0.1,0.25,.1/)
c updated for current CEUS weighting scheme
        wte = (//0.06,0.11,0.06,0.10,0.11,0.11,0.08,0.22,.15/) !revised June 28 2013
        if(iargc().lt.2)stop'Usage: assim2013 infile outfile'
        call getarg(1,fi)
        open(1,file=fi,status='old',err=2013)
        ind=index(fi,'.DET')-4
        fi2=fi
        fi ='F'//fi(1:ind)//'in'
        print *,fi
        inquire(file=fi,exist=namef)
        if(namef)open(10,file=fi,status='old')
1001        format(i4,18x,a36)
3001        format(i4,t17,a36)
        irec=index(rec,' ')
        call getarg(2,fi)
        print *,'Enter a 0 if CEUS 2013 file;'
        print 50,'Enter 1 if this is a 2013 Cascadia SUBDUCTION-Source file;Enter 2 if this is a WUS 2013 source file: '
        read *,i
        print *,i
        ceus=i.eq.0
        subduction=i.eq.1
        usewt=.false.
        wus13=i.eq.2
        if(wus13)wtw=6.*wtw
1        read(1,5,end=24)rec
        dowhile(rec(1:1).eq.'#')
        read(1,5)rec
        enddo
        more = .true.
        if(namef)then
        do
        if(i.eq.0)then
        read(10,3001,end=12)ift,fname(ift)
        else
        read(10,1001,end=12)ift,fname(ift)
        endif	
        enddo
12        close(10)
        print 1001,ift,fname(ift)

        endif

32        print 50,'Enter number of attenuation models (1 to 9): '
        read *,namx
        ic=2*namx+1
        if(namx.eq.9)then
        print 50,'Enter 1 to use NSHMP CEUS model wts of 2014 when
     + averaging: '
         read *,i
         usewt=i.eq.1
         wte=wte*float(namx)        !average calc divides by namx later on
         elseif (subduction)then
         usewtc=.true.
         wtc=wtc*4.
         else
         usewtc=.false.
         usewt=.false.
         endif
        if(namx.lt.1.or.namx.gt.9)goto 32
50        format(a,$)
        open(2,file=fi,status='unknown')
        if(usewt)write(2,202)namx,fi2(1:30),usewt
        if(usewtc)write(2,204)namx,fi2(1:30),usewtc
202        format('#assim.2013.f Natten models: ',i1,' infile ',a,' CEUSwts?',
     + l4)
             if(wus13)write(2,2014)namx,fi2(1:30),wus13
204        format('#assim.2013.f Natten models: ',i1,' infile ',a,
     + ' Use Cascadia attn-model wts?',l4)
2014        format('#assim.2013.f Natten models: ',i1,' infile ',a,
     + ' Use WUS 2014 attn-model wts?',l4)
        do while (more)
        recold=rec
5        format(a)
        read(1,5,end=24)rec
        if(rec(1:1).eq.'!')goto 10	!keep reading
        read(recold,6)x,y
        backspace(1)	!data record for a site. have to take another look
        sbmax=0.
6        format(2x,f8.3,1x,f7.3)
4        sabar=0.0
        do j=1,namx
        read(1,5,end=24)rec
        if(rec(1:1).eq.'!')then
        recold=rec
        goto 2
        endif
        read(rec,*)sa(j),sd(j),iat,iflt,xmag,rjb,rrup
c        if(x.ge.-99.9 .and. y .eq. 34.8)print *,sa(j),sd(j),iat,iflt,xmag,rjb,rrup
        if(usewt)then
        j1=1
        dowhile(iat.ne.indexe(j1))
        j1=j1+1
        if(j1.gt.namx)stop'iat does not match any ceus index'
        enddo
        sabar=sabar+wte(j1)*sa(j)
        elseif(wus13)then
        j1=1
        dowhile(iat.ne.indw(j1))
        j1=j1+1
        if(j1.gt.namx)stop'iat does not match any WUS-2013 index'
        enddo
        sabar=sabar+wtw(j1)*sa(j)
        
        elseif(usewtc)then
        j1=1
        dowhile(iat.ne.indc2(j1))
        j1=j1+1
        if(j1.gt.4)stop'iat does not match valid Cascadia atten indexes'
        enddo
        sabar = sabar+wtc2(j1)*sa(j)
        else
        sabar=sabar+sa(j)
        endif
        enddo	!loop thru namx attn models
        if(sabar.gt.sbmax)then
        do j=1,namx
        samax(j)=sa(j)
        sdmax(j)=sd(j)
        xmagmax=xmag
        enddo
        rjbmax=rjb
        rrupmax=rrup
        ifltmx=iflt
        sbmax=sabar
        endif
        written=.false.
        goto 4
2        if(namef)then
           if(subduction) write(2,222)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),fname(ifltmx)
           if(wus13) write(2,224)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),fname(ifltmx)
           if(subduction) write(2,226)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),fname(ifltmx)
222        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,9(1x,e11.5),1x,a)
224        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,13(1x,e11.5),1x,a)
226        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,19(1x,e11.5),1x,a)
322        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,9(1x,e11.5),1x,f7.0)
324        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,13(1x,e11.5),1x,f7.0)
326        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,19(1x,e11.5),1x,f7.0)
        else
        write(2,22)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),float(ifltmx)
             endif	!write with or without name of fault?
        written=.true.
22        format(f8.2,1x,f6.2,1x,f7.2,1x,f7.2,1x,f7.2,20(1x,e11.5))
10        continue
        enddo	!dowhile more to do
24        if(.not.written)then
        if(namef)then
        if(subduction) write(2,222)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),fname(ifltmx)
        if(wus13) write(2,224)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),fname(ifltmx)
        if(ceus) write(2,226)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),fname(ifltmx)
        else
        write(2,22)x,y,xmagmax,rjbmax,rrupmax,sbmax/namx,(samax(k),k=1,namx),
     + (sdmax(k),k=1,namx),float(ifltmx)
             endif	!write with or without name of fault?
             endif	!if written
        close(2)
        print *,'output file is ',fi
        stop
2013        stop'input file not found'
        end
