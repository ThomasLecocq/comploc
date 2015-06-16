!-----------------------------------------------------------------------
!
!     comploc locates  earthquakes using Shrinking Box SSST method.
!
!     We provide Grid-search (nonlinear) method.
!
!     GRIDLOC only includes location part; 
!     the station terms calculation part is in the main program.
!     10/30/04  Guoqing Lin 
! 
!     For large regions, using all the events to calculate the ssst is 
!     unrealistic, so still apply starting and ending nmed
!     01/02/05   Guoqing Lin
!
!     Apply ROBUST MEAN to location algorithms
!     01/03/05/  Guoqing Lin
!
!     skip picks with tt picks le. zero .or. delta le. zero
!     03/01/05 Guoqing Lin
!
!     output location quality -1 --- bad 
!                              0 --- good
!     DO NOT output those events with bad locations. 
!                              06/23/05 Guoqing Lin
!
!     Allow to read in station terms
!                              08/22/05 Guoqing Lin 
!
!     In subroutine GET_STPPHASE, if the phase is not P or S , then skip.
!                                                09/21/2005 Guoqing Lin 
!
!     MEDIAN changes the order of an array.
!     residual index number is not consistent.  Fixed.
!                                                10/19/2005 Guoqing Lin
!             
!     subroutine STPPHASE format change due SCSN STP PHASE DATA format 
!     change                                     01/18/2006 Guoqing Lin
!
!-----------------------------------------------------------------------
!
      program comploc        
      implicit none
!
!     also in subroutines
      integer nq0               !maximum number of events
      parameter (nq0=500000)
      integer npick0            !maximum number of picks
      parameter (npick0=13000000)
      integer nsta0             !maximum number of stations
      parameter (nsta0=3000)
      integer mp                !maximum number of picks per event
      parameter (mp=1000) 
!     also in the subroutines

      integer nqtot          !total number of events from phase data file
      integer npickall       !total number of picks from phase data file
      integer maxev          !maximum number of events to locate given by user
      integer nq             !number of events to be located
      integer npick          !number of picks
      integer idstmax        !number of stations
      integer npick_min      !minimum number of picks per event
      integer nmed_ssst1     !starting number of nearby events for source specific station terms calculation
      integer nmed_ssst2     !ending number of nearby events for source specific station terms calculation
      integer nmed_ssst      !number of nearby events for each station terms calculation iteration
      integer ntermit1       !number of iterations for STATIC STATION TERMS CALCULATION
      integer ntermit2       !number of iterations for SOURCE SPECIFIC STATION TERMS CALCULATION
      integer ntermit        !number of totle iterations  
      integer ilocfix        !(0) normal  or  (1) fixed loc kluge
      integer idepfix        !(0) catalog depth  or  (1) qdepref given by user
      integer inorm          !(1) L1 norm, (2) L2 norm or (3) ROBUST MEAN
      integer iformat        !Output Location Format
!                            (1) SCSN Format
!                            (2) HYPOINVERSE Format
!                            (3) NCEDC READABLE Format
!                            (4) HYP71 Format
!                            (5) Our Format
      integer nit            !number of iterations for grid search
      integer iyr            !Year of event
      integer imon           !Month of event
      integer idy            !Day of event
      integer ihr            !Hour of event
      integer imn            !Minute of event
      integer ilat           !integer part of qlat
      integer ilon           !integer part of qlon 
      integer npspick        !total picks number of each event 
      integer iq
      integer is
      integer ipick
      integer ipha
      integer i4
      integer idphase
      integer i
      integer j
      integer i1
      integer i2
      integer itermit
      integer iter_ssst
      integer nrms
      integer ipick1(nq0)      !array with first pick index for each event
      integer ipick2(nq0)      !array with last pick index for each event 
      integer nppick(nq0)      !array with P picks number for each event
      integer nspick(nq0)      !array with S picks number for each event
      integer nqqual(npick0)   !array with earthquake location quality
!                               =-1, bad location, should not be used in station terms calculation
      integer*4 idcusp0(nq0)   !array with cuspid for each event
      integer*2 idph(npick0)   !array with phase id for each pick
      integer*2 idsta(npick0)  !array with station index recording each pick
      integer nstpick(nsta0,2) !array with P & S picks number recorded by each station

      real flat1,flat2,flon1,flon2     !location region window: minlat, maxlat, minlon, maxlon
      real delmax           !cutoff distance between stations and events
      real efrac            !random fraction of events to read from phase data file 
      real qdepref          !reference starting depth
      real dismax1          !starting cutoff distance for source specfic station terms calculation
      real dismax2          !ending cutoff distance for source specific station terms calculation
      real dismax           !cutoff distance for source specific station terms calculation
      real fracsk           !shrinking fraction for grid search (e.g. 0.67)
      real fsc              !seconds of event
      real elat             !decimal part of qlat
      real elon             !decimal part of qlon
      real maxazi           !maximum azimuth (unknown)
      real mindis           !minimum distance between stations and events (unknown) 
      real errx             !error in x location (unknown)
      real erry             !error in y location (unknown)
      real errz             !error in z location (unknown)
      real errt             !error in origin time (unknown) 
      real srms
      real s1
      real s2
      real s3
      real frac
      real slog
      real qlat(nq0)        !array with quake latitude
      real qlon(nq0)        !array with quake longitude
      real qdep(nq0)        !array with quake depth
      real qorg(nq0)        !array with quake origin time
      real qlat0(nq0)
      real qlon0(nq0)
      real qdep0(nq0)
      real qorg0(nq0)
      real qlat2(nq0)
      real qlon2(nq0)
      real qdep2(nq0)
      real qorg2(nq0)
      real qmag(nq0)        !array with quake magnitude
      real rms(nq0)         !array with RMS residuals for each event
      real rmed(nq0)        !array with Median Absolute Value of residuals for each event
      real*8 tsec0(nq0)     !event time (seconds since 1600)
      real stlat(nsta0)     !array with station latitude
      real stlon(nsta0)     !array with station longitude
      real stelev(nsta0)    !array with station elevation (not in use) 
      real pick(npick0)     !array with travel time picks for each pick from phase data file
      real resid(npick0)    !array with travel time residuals for each pick
      real term(npick0)     !station terms for each pick
      real term2(npick0) 
      real stterm(nsta0,2)  !station term at each station
      real stterm2(nsta0,2)
      real stterm0(2)
 
      character*2 cmon      !Month of event in character
      character*2 cdy       !Day of event in character 
      character*2 chr       !Hour of event in character
      character*2 cmn       !Minuter of event in character
      character*4 cyr       !Year of event in character
      character*2 pickph(npick0)
      character*12 stname
      character*12 sname0(npick0)
      character*100 infile        !input phase data file name
      character*100 outfile       !output location file name
      character*100 stlocfile     !station location file name
      character*100 tname(2)      !travel time table file name
!
!-----------------------------------------------------------------------
!
!     read input 
!
      print *,'Enter P-TT file name'
      read (*,'(a)') tname(1)
      print *,'Enter S-TT file name'
      read (*,'(a)') tname(2)
   
      print *,'Phase file format '
      print *,' (1) BED3 (recommended)'
      print *,' (2) SCSN (OLD STP)'
      print *,' (3) HYPOINVERSE'
      print *,' (4) SCSN (NEW STP)'
      read *,ipha

      print *,'Enter phase file name'
      read (*,'(a)') infile

      print *,'Enter stlist file name for station locations'
      read (*,'(a)') stlocfile

      print *,'Output Location format:'
      print *,'(1) SCSN format'
      print *,'(2) HYPOINVERSE (Y2K)'
      print *,'(3) NCEDC readable'
      print *,'(4) HYPO71 format'
      print *,'(5) Our Format '
      read *,iformat

      print *,'Enter output file name for locations'
      read (*,'(a)') outfile

      print *,'Enter delmax (km)'
      read *,delmax

      print *,'Enter npick_min'
      read *,npick_min

      print *,'Enter flat1,flat2,flon1,flon2 window for quakes'
      read *,flat1,flat2,flon1,flon2

      print *,'Enter:  (0) normal  or  (1) fixed loc kluge'
      read *,ilocfix

      print *,'Enter:  (0) catalog depth  or  (1) qdepref'
      read *,idepfix

      print *,'Enter the reference starting depth'
      read *,qdepref

      print *,'Enter max event number to read (999999 for all)'
      read *,maxev

      print *,'Enter random frac of events to get (2 for all)'
      read *,efrac
!
!     for single-event location, set ntermit1=0 and ntermit2=0
!
      print *,'Enter number of iterations for static station terms'
      read *,ntermit1

      print *,'Enter number of iterations for ssst'
      read *,ntermit2
      ntermit = 1 + ntermit1 + ntermit2

      print *,'Enter starting nmed, dmax for SSST'
      read *,nmed_ssst1,dismax1

      print *,'Enter ending nmed, dmax for SSST'
      read *,nmed_ssst2,dismax2

      print *,'Enter NORM (1) L1, (2) L2  or (3) ROBUST MEAN'
      read *,inorm

      print *,'Enter Shrinking iteration number '
      print *,'    and fraction for Grid-search'
      read *,nit,fracsk 
!
!-----------------------------------------------------------------------
!
!     read phase data
!
      print *,'Reading phase data...'
  
      if (ipha .eq. 1) then          !BED3 format 
      call GETBED3_R(infile,flat1,flat2,flon1,flon2,delmax,efrac,
     &            nqtot,npickall,idcusp0,tsec0,qlat0,qlon0,qdep0,
     &            qmag,ipick1,ipick2,sname0,pick,pickph)

      else if (ipha .eq. 2.or.ipha.eq.4) then     !SCSN STP_PHASE format
      call GET_STPPHASE(ipha,infile,flat1,flat2,flon1,flon2,
     &         delmax,efrac,
     &         nqtot,npickall,idcusp0,tsec0,qlat0,qlon0,qdep0,qmag,
     &         ipick1,ipick2,sname0,pick,pickph)

      else if (ipha .eq. 3) then     !HYPOINVERSE format
      call GET_HYPOINVERSE(infile,flat1,flat2,flon1,flon2,delmax,efrac,
     &         nqtot,npickall,idcusp0,tsec0,qlat0,qlon0,qdep0,qmag,
     &         ipick1,ipick2,sname0,pick,pickph)

      end if
      print *,'nqtot=',nqtot,'  npickall=',npickall

      print *,'Locations of first two events follow:'
      do i=1,2
         print *,i,idcusp0(i),tsec0(i),qlat0(i),qlon0(i),qdep0(i),
     &          qmag(i)
      enddo

      print *,'Locations of last two events follow:'
      do i=nqtot-1,nqtot
         print *,i,idcusp0(i),tsec0(i),qlat0(i),qlon0(i),qdep0(i),
     &          qmag(i)
      enddo

      npickall=0
      do 60 iq=1,nqtot
            qlat(iq)=qlat0(iq)
            qlon(iq)=qlon0(iq)
            qdep(iq)=qdep0(iq)
            qorg(iq)=qorg0(iq)
         do 50 ipick=ipick1(iq),ipick2(iq)
            stname=sname0(ipick)
            call GET_STANUM(stname,i4,idstmax,1)   !assign station number

            if (pickph(ipick)(1:1).eq.'P'.or.
     &          pickph(ipick)(1:1).eq.'p') then
               nppick(iq)=nppick(iq)+1
               npickall=npickall+1
               idph(ipick)=1
               idsta(ipick)=i4
               idphase=idph(ipick)   !convert to I*4
               nstpick(i4,idphase)=nstpick(i4,idphase)+1
            else if (pickph(ipick)(1:1).eq.'S'.or.
     &               pickph(ipick)(1:1).eq.'s') then
               nspick(iq)=nspick(iq)+1
               npickall=npickall+1
               idph(ipick)=2
               idsta(ipick)=i4
               idphase=idph(ipick)   !convert to I*4
               nstpick(i4,idphase)=nstpick(i4,idphase)+1
            end if

            if (idsta(ipick).eq.0) then
               print *,'***TROUBLE: ',iq,ipick,i4,idsta(ipick)
               print *,'   pickph(ipick) = ',pickph(ipick)
               print *,'Program requires all picks be P or S'
               stop
            end if

50       continue
60    continue
      print *,'nqtot=',nqtot,'   maxev=',maxev
      nq=min(nqtot,maxev)
      npick=npickall
      if (idstmax .gt. nsta0) then 
      print *,' idstmax > nsta0, increase nsta0 ...'
      print *,idstmax,nsta0
      stop
      end if
      print *, 'nq,npick,idstmax = ',nq,npick,idstmax

      do 8300 i=1,idstmax
         call GET_STANUM(stname,i,idstmax,2)   !get station name
         call GETSTAT_STLIST(stlocfile,stname,s1,s2,s3,stterm0) !get station location 
         stlat(i)=s1
         stlon(i)=s2
         stelev(i)=s3
         stterm(i,1)=stterm0(1)
         stterm(i,2)=stterm0(2)
         if (stlat(i) .ne. 999.) stelev(i)=stelev(i)/1000.    !in km
!         print *,i,s1,s2,s2,stterm(i,1),stterm(i,2)
8300  continue
!
!-----------------------------------------------------------------------
!
!     Location part
!
      print *,'Now locating quakes...'

      open(11,file=outfile)
      errx=0.0
      erry=0.0
      errz=0.0
      errt=0.0

      if (idepfix .eq. 1) then
      do iq=1,nq
         qdep(iq)=qdepref
      end do
      end if 

      do i=1,npick0
      resid(i)=0.0
      term(i)=stterm(idsta(i),idph(i))
      end do

      do 400 itermit=1,ntermit       !loop on station terms

      call GRIDLOC(tname,nq,idph,idsta,npick_min,
     &    ipick1,ipick2,pick,term,
     &    qlat,qlon,qdep,qorg,ilocfix,idepfix,inorm,nit,
     &    fracsk,stlat,stlon,qlat2,qlon2,
     &    qdep2,qorg2,nqqual,resid,rms,rmed)

!
      print *,'Finished location, iteration = ',itermit
      print *,'   '
 
      if (itermit .eq. ntermit) then !write new locations into outfile in different format
      do 299 iq=1,nq 
         i1=ipick1(iq)
         i2=ipick2(iq)
         if (i2-i1+1 .lt. npick_min) go to 299 !exclude events with too few picks
         if (nqqual(i1) .eq. -1) go to 299 !exclude events with bad locations 
         call DT_ADDTIME(1600,1,0,0,0,0.,
     &        iyr,imon,idy,ihr,imn,fsc,tsec0(iq)+dble(qorg2(iq)))

         if (iformat.eq.1) then   !SCSN format
         write(cyr,'(i4)') iyr
         write(cmon,'(i2)') imon
         if (cmon(1:1).eq.' ') cmon='0'//cmon(2:2) 
         write(cdy,'(i2)') idy
         if (cdy(1:1).eq.' ') cdy='0'//cdy(2:2) 
         write(chr,'(i2)') ihr
         if (chr(1:1).eq.' ') chr='0'//chr(2:2) 
         write(cmn,'(i2)') imn
         if (cmn(1:1).eq.' ') cmn='0'//cmn(2:2) 

         ilat=int(qlat2(iq))
         elat=(qlat2(iq)-float(ilat))*60.0

         ilon=int(qlon2(iq))
         elon=(float(ilon)-qlon2(iq))*60.0

         npspick=nppick(iq)+nspick(iq)

         write (11,161) cyr,cmon,cdy,chr,cmn,fsc,ilat,elat,
     &           ilon,elon,'Z',qmag(iq),qdep2(iq),
     &           npspick,rms(iq),idcusp0(iq)


         else if (iformat .eq. 2) then    ! HYPOINVERSE (Y2K)
         write(cyr,'(i4)') iyr
         write(cmon,'(i2)') imon
         if (cmon(1:1).eq.' ') cmon='0'//cmon(2:2)
         write(cdy,'(i2)') idy
         if (cdy(1:1).eq.' ') cdy='0'//cdy(2:2)
         write(chr,'(i2)') ihr
         if (chr(1:1).eq.' ') chr='0'//chr(2:2)
         write(cmn,'(i2)') imn
         if (cmn(1:1).eq.' ') cmn='0'//cmn(2:2)

         ilat=int(qlat2(iq))
         elat=(qlat2(iq)-float(ilat))*60.

         ilon=int(qlon2(iq))
         elon=(float(ilon)-qlon2(iq))*60.

         npspick=nppick(iq)+nspick(iq)
         write (11,162) cyr,cmon,cdy,chr,cmn,int(fsc*100.0),ilat,' ',
     &           int(elat*100.0),-ilon,' ',int(elon*100.0),
     &           int(qdep2(iq)*100.0),0,npspick,0,0,nint(rms(iq)*100.0),
     &           idcusp0(iq)
     
         else if (iformat .eq. 3) then    ! NCEDC readable 
         write(cyr,'(i4)') iyr
         write(cmon,'(i2)') imon
         if (cmon(1:1).eq.' ') cmon='0'//cmon(2:2)
         write(cdy,'(i2)') idy
        if (cdy(1:1).eq.' ') cdy='0'//cdy(2:2)
         write(chr,'(i2)') ihr
         if (chr(1:1).eq.' ') chr='0'//chr(2:2)
         write(cmn,'(i2)') imn
         if (cmn(1:1).eq.' ') cmn='0'//cmn(2:2)

         write (11,163) cyr,'/',cmon,'/',cdy,chr,':',cmn,':',fsc,
     &           qlat2(iq),qlon2(iq),qdep2(iq),qmag(iq),'ML',
     &           ipick2(iq)-ipick1(iq)+1,180,100,rms(iq),'SHLK',
     &           idcusp0(iq)


         else if (iformat.eq.4)  then      !HYPO71 format
         ilat=int(qlat2(iq))
         elat=(qlat2(iq)-float(ilat))*60.

         ilon=int(qlon2(iq))
         elon=(float(ilon)-qlon2(iq))*60.

         npspick=nppick(iq)+nspick(iq)

         maxazi=0.0
         mindis=0.0

         write (11,164) iyr,imon,idy,ihr,imn,fsc,float(ilat),' ',
     &           elat,float(-ilon),' ',elon,qdep2(iq),' ',qmag(iq),
     &           npspick,maxazi,mindis,rms(iq),errx,errz,' ',' ',' ',
     &           idcusp0(iq)

         else                              !Our own format
         write (11,165) iyr,imon,idy,ihr,imn,fsc,qlat2(iq),qlon2(iq),
     &           qdep2(iq),errx,erry,errz,errt,idcusp0(iq),
     &           qmag(iq),nppick(iq),nspick(iq),qorg2(iq),rms(iq),
     &           nqqual(i1)

         
161      format(a4,1x,a2,1x,a2,2x,a2,1x,a2,1x,f5.2,1x,i3,1x,f5.2,
     &          i4,1x,f5.2,1x,a1,1x,f3.1,4x,f6.2,i3,4x,f5.2,1x,i8)


162      format(a4,4a2,i4,i2,a1,i4,i3,a1,i4,i5,i3,i3,i3,i3,i4,
     &          84x,i10)

163      format(a4,1a,a2,1a,a2,1x,a2,a1,a2,a1,f5.2,f9.4,f10.4,
     &          f7.2,f6.2,3x,a2,i5,i4,i5,f5.2,1x,a4,1x,i10)

164      format(i4,2i2,1x,2i2,f6.2,f3.0,a1,f5.2,f4.0,a1,
     &          f5.2,f7.2,1x,a1,f5.2,i3,f4.0,f5.1,f5.2,
     &          f5.1,f5.1,a1,a1,a1,1x,i10,1x,a1,a3)

165      format (i4,4i3,f7.3,f8.4,f10.4,f6.2,4f6.2,i11,
     &           f4.1,2i4,2f8.2,2x,i3)
         end if

299      continue
         end if
!
!     static station terms calculation
!
      if (itermit.le.ntermit1.and.itermit.ne.ntermit
     &           .or.nmed_ssst1.eq.0) then
         print *,'Computing static station terms...'

         call STATERM(inorm,resid,idsta,idph,npick,idstmax,
     &                2,nstpick,nqqual,term2,stterm2)
         do i=1,npick
            term(i)=term(i)+term2(i)
         end do
         srms=0.
         nrms=0
         do is=1,idstmax
            do j=1,2
               stterm(is,j)=stterm(is,j)+stterm2(is,j)
               if (nstpick(is,j).ge.5) then
                  srms=srms+stterm2(is,j)**2
                  nrms=nrms+1
               end if
            end do
         end do 
         if (nrms .gt. 0) srms=sqrt(srms/float(nrms))
!
!     source-specific station terms calculation
!
      else if (itermit.ne.ntermit) then
         print *,'Computing SSST station terms...'

         iter_ssst = itermit - ntermit1
         frac = real(iter_ssst-1)/real(ntermit2-1)

         s1=alog10(real(nmed_ssst1))
         s2=alog10(real(nmed_ssst2))
         slog = s1+frac*(s2-s1)
         nmed_ssst = 10**slog

         s1=alog10(dismax1)
         s2=alog10(dismax2)
         slog = s1+frac*(s2-s1)
         dismax = 10**slog

         print *,'iter_ssst,frac,dismax = ',iter_ssst,frac,dismax

         call SSST_DISTF(inorm,resid,idsta,idph,npick,nstpick,nq,
     &         qlat2,qlon2,qdep2,nqqual,ipick1,ipick2,
     &         nmed_ssst,dismax,term2)
         srms=0.           
         do i=1,npick
            term(i)=term(i)+term2(i)
            srms=srms+term2(i)**2
         end do
         srms=sqrt(srms/float(npick))

      print *,'Finished station term calculation'
      print *,'rms of change in station terms = ',srms
      end if

      if (itermit.eq.ntermit) go to 801
 
400   end do
801   close(11)

      end
!
!
!-----------------------------------------------------------------------
!
!     subroutine(01)
!
!     GETBED3_R extracts quake info and picks from bed3 file
!     (which can be generated by phase2bed3.f program).
!
!     It's much faster to read phase data in BED3 format than
!     the direct phase data from the network for large data sets.
!  
!     Note:  only extracts events within qlatmin, etc. limits
!            and only gets picks within delmax distance
!            This version gets P or S picks only
!         
!            only takes first pick when more than one pick
!            from the same station is present, a common occurrence
!            with NCSN hypoinverse phase data
!
!            permits getting only random fraction of events
!
!
!  Input:  infile --- phase data file in BED3 format
!          qlatmin, qlatmax, qlonmin, qlonmax --- location limits window
!          delmax --- cutoff distance
!          frac --- fraction of total events to read 
!
!  Output:   nqtot --- number of earthquakes 
!            npickall --- number of picks
!            idcusp0 --- array with event identification number
!            tsec0 --- event time (seconds since 1600)
!            qlat0 --- array with event latitude
!            qlon0 --- array with event longitude
!            qdep0 --- array with event depth
!            qmag0 --- array with event magnitude
!            ipick1 --- array with the first pick index for each event
!            ipick2 --- array with the last pick index for each event
!            sname --- array with station id (12 characters)
!                      including network id, station name and channel.   
!            pick --- array with travel time pick for each pick 
!            pickph --- array with phase id for each pick 
!
!-----------------------------------------------------------------------
!
      subroutine GETBED3_R(infile,qlatmin,qlatmax,qlonmin,qlonmax,
     & delmax,frac,nqtot,npickall,idcusp0,tsec0,qlat0,qlon0,qdep0,
     & qmag0,ipick1,ipick2,sname,pick,pickph)
      implicit none

      integer nq0               !maximum number of events
      parameter (nq0=500000)
      integer npick0            !maximum number of picks
      parameter (npick0=13000000)
      integer mp                !maximum number of picks per event
      parameter (mp=1000) 

      integer i
      integer iend
      integer iev
      integer ipick
      integer iq
      integer j
      integer nn
      integer npickall
      integer nqmax
      integer nqtot
      integer*4 idcusp0(nq0)
      integer ipick1(nq0)
      integer ipick2(nq0)
!
!     variables for header
      integer*2 nevent
      integer*2 header
!
!     buffer variables
      integer*2 ibufp(18,mp)
      integer*2 ibufpick(18)
      integer*2 ibufq(40)
!
!     variables for quake info
      integer*4 cuspid
      integer*2 idip
      integer*2 imag
      integer*2 irake
      integer*2 istrike
      integer*2 npicktot
      integer*2 npickused
      integer*2 nppick
      integer*2 nspick
      integer*2 ntrace
!
      real delmax
      real qlatmin,qlatmax,qlonmin,qlonmax 
      real*4 qdep
      real*4 qdeperr
      real*4 qlat
      real*4 qlaterr
      real*4 qlon
      real*4 qlonerr
      real*4 rms
      real*8 t0
      real*4 t0err
!
!     variables for station pick info
      real*4 delta
      real frac
      real fran
      real RAND
      real*4 tt
      real pick(npick0)
      real qdep0(nq0)
      real qlat0(nq0)
      real qlon0(nq0)
      real qmag0(nq0)
      real*8 tsec0(nq0)
!
      character*1 errtype
      character*1 foctype
      character*1 mtype
      character*1 qtype
      character*1 qual
!
      character*1 acc
      character*3 comp
      character*4 junk
      character*1 onset
      character*6 phase
      character*1 phaseold
      character*1 pol
      character*12 stname
      character*7 stnameold
!
      character*100 infile
      character*2 pickph(npick0)
      character*12 sname(npick0)
!
!
      common/com1/t0,qlat,qlon,qdep,qlaterr,qlonerr,qdeperr,t0err,
     &            rms,cuspid,imag,istrike,idip,irake,ntrace,
     &            nppick,nspick,npicktot,npickused,qtype,mtype,qual,
     &            errtype,foctype

      common/com2/tt,delta,stname,comp,phase,onset,pol,acc,junk
!
      equivalence (ibufq(1),t0)
      equivalence (ibufpick(1),tt)
!
!-----------------------------------------------------------------------
!
      print *,'Reading BED3 file ...'

      open (11,file=infile,form='unformatted')

      call RDBLK(nevent,1,11,iend)
      call RDBLK(header,1,11,iend)

      nqmax=999999
      ipick=0
      iev=0
      do 50 iq=1,nqmax
10       call RDBLK(ibufq,40,11,iend)
         if (iend.eq.1) go to 150
         nn=18*npicktot
         call RDBLK(ibufp,nn,11,iend)

         if (qlat.lt.qlatmin.or.qlat.gt.qlatmax) go to 50
         if (qlon.lt.qlonmin.or.qlon.gt.qlonmax) go to 50

         fran=RAND(0)
         if (fran.gt.frac) go to 50

         iev=iev+1
         if (iev.gt.nq0) then
            print *,'***Too many events for GETBED3'
            print *,'***Truncated at nq = ',nq0
            iev=iev-1
            go to 150
         end if 

         idcusp0(iev)=cuspid
         tsec0(iev)=t0
         qlat0(iev)=qlat
         qlon0(iev)=qlon
         qdep0(iev)=qdep
         qmag0(iev)=float(imag)/1000.
         if (qmag0(iev).lt.0.) qmag0(iev)=0.
         ipick1(iev)=ipick+1

         stnameold='1234567'
         phaseold=' '
         do 20 i=1,npicktot
            do 15 j=1,18
               ibufpick(j)=ibufp(j,i)
15          continue
            if (delta.gt.delmax .or. delta .le. 0.0) go to 20
            if (tt.le.0.0) go to 20

!            if (onset.ne.'I'.and.
            if(stname(1:7).eq.stnameold .and.
     &         phase(1:1).eq.phaseold) go to 20

            if (phase(1:1).ne.'p'.and.phase(1:1).ne.'P'.and.
     &          phase(1:1).ne.'s'.and.phase(1:1).ne.'S') go to 20

            ipick=ipick+1
            sname(ipick)=stname
            pick(ipick)=tt
            pickph(ipick)=phase(1:2)

            stnameold=stname(1:7)
            phaseold=phase(1:1)
          
20       continue

         ipick2(iev)=ipick

50    continue
150   close (11)
      nqtot=iev
      npickall=ipick
      print *,'nqtot,npickall = ',nqtot,npickall

      print *,'Last quake follows: '
      print *,t0,qlat,qlon,qdep,cuspid,imag,npicktot
      do 200 i=1,npicktot
         do 195 j=1,18
            ibufpick(j)=ibufp(j,i)
195      continue
         print 197, stname,delta,comp,phase,tt,onset,pol,acc,junk
197      format (a12,x,f9.3,x,a3,x,a6,x,f8.3,x,a1,x,a1,x,a1,x,a4)
200   continue

      return
      end
!
!
!
!-----------------------------------------------------------------------
!
!  subroutine(02)
!
!  GET_HYPOINVERSE extracts quake info and picks from STP PHASE DATA file
!  Note:  only extracts events within qlatmin, etc. limits
!         and only gets picks within delmax distance
!         This version gets P or S picks only
!         
!         Only takes first pick when more than one pick
!         for the same station is present
!
!         permits getting only random fraction of events
!
!   inputs and outputs are the same as in GETBED3_R subroutine.
!
!-----------------------------------------------------------------------
!
      subroutine GET_HYPOINVERSE(infile,qlatmin,qlatmax,qlonmin,qlonmax,
     &    delmax,frac,nqtot,npickall,idcusp0,tsec0,qlat0,qlon0,qdep0,
     &    qmag0,ipick1,ipick2,stname0,pick,pickph)
      implicit none

      integer nq0                !maximum number of events
      parameter (nq0=500000)
      integer npick0             !maximum number of picks
      parameter (npick0=13000000)

      integer*4 cuspid
      integer i1
      integer i2
      integer i3
      integer i4
      integer i5
      integer i6
      integer iflag
      integer iline
      integer ipick
      integer iq
      integer isc
      integer iyr
      integer imon
      integer idy
      integer ihr
      integer imn
      integer isec
      integer idelkm
      integer npickall        !total number of picks (output)
      integer nqmax
      integer nqtot           !total number of events (output)
      integer qyr
      integer qmon
      integer qdy
      integer qhr
      integer qmn
      integer*4 idcusp0(nq0)  !array with event identification number (output)
      integer ipick1(nq0)     !array with first pick index for each event (output)
      integer ipick2(nq0)     !array with last pick index for each event (output)


      real qlatmin,qlatmax,qlonmin,qlonmax   !location region window (input)
      real delmax       !cutoff distance between stations and events (input)
!
!     variables for quake info
      real frac          !random fraction of events to read from phase data file (input)
      real fran          
      real qdep
      real qlat
      real qlon
      real qmag
      real qsc
      real qt0
      real RAND
      real*8 t0
      real pick(npick0)  !array with travel time pick for each pick (output)
      real qdep0(nq0)    !array with event depth (output)
      real qlat0(nq0)    !array with evetn latitude (output)
      real qlon0(nq0)    !array with event longitude (output)
      real qmag0(nq0)    !array with event magnitude (output)
      real*8 tsec0(nq0)  !event time (seconds since 1600) (output)
!
!     variables for station pick info
      real delkm
      real delta
      real sec
      real tt
!
      character acc,onset,pol
      character*1 chr1
      character*1 chr2
      character*3 comp
      character*100 infile         !input phase data file name (input)
      character*2 iph
      character*150 linebuf
      character phaseold
      character*11 phinfo
      character*12 stname
      character*12 stname2
      character*7 stnameold
      character*2 pickph(npick0)   !array with phase id for each pick (output)
      character*12 stname0(npick0)   !array with station id (12 characters)
!                                     including network id, station name and channel.  (output)  
!
!-----------------------------------------------------------------------
!     
      do iq=1,nq0
      ipick1(iq)=0
      ipick2(iq)=0
      enddo
      
      iflag=1
      nqmax=999999
      open(11,file=infile,status='old')
      iq=0
      ipick=0
      
      do 401 iline=1,nqmax
      read(11,'(a)',end=101) linebuf
      if (linebuf(1:2) .eq. '19' .or. linebuf(1:2) .eq. '20') then   !Quake Line
      iflag=1
      if (iq .ne. 0) ipick2(iq)=ipick
      
      read (linebuf,201) qyr,qmon,qdy,qhr,qmn,isc,i1,chr1,i2,
     &            i3,chr2,i4,i5,cuspid,i6
!      print 201, qyr,qmon,qdy,qhr,qmn,isc,i1,chr1,i2,
!     &            i3,chr2,i4,i5,cuspid,i6
201   format (i4,4i2,i4,i2,a1,i4,i3,a1,i4,i5,100x,i10,1x,i3)
      
         qsc=float(isc)/100.
         qlat=float(i1)+float(i2)/6000.
         if (chr1.eq.'s'.or.chr1.eq.'S') qlat=-qlat
         qlon=-(float(i3)+float(i4)/6000.)
         if (chr2.eq.'w'.or.chr2.eq.'W') qlon=-qlon
         qdep=float(i5)/100.
         qmag=float(i6)/100.

      if (qlat.lt.qlatmin.or.qlat.gt.qlatmax)  then
         iflag=0
         goto 401
      end if
      if (qlon.lt.qlonmin.or.qlon.gt.qlonmax)  then
         iflag=0
         goto 401
       end if
      
      fran=RAND(0)
      if (fran.gt.frac) then
         iflag=0
         go to 401
      end if

      iq=iq+1
      if (iq.gt.nq0) then
         print *, '***Too many events for HYPOINVERSE'
         print *, '***Truncated at nq = ',nq0 
         stop
      end if 

      qlat0(iq)=qlat
      qlon0(iq)=qlon
      qdep0(iq)=qdep
      qmag0(iq)=qmag
      qt0=0.0
      call DT_TIMEDIF(1600,1,0,0,0,0.,qyr,qmon,qdy,qhr,qmn,qsc,t0)
      tsec0(iq)=t0

      stnameold='1234567'
      phaseold=' '

      else if (linebuf(1:4) .ne. '    ') then  !Phase Line

         if (iq.eq.0 .or.
     &       linebuf(6:7).eq.'  ' .or.
     &       linebuf(31:31).eq.'.' .or.
     &       linebuf(32:32).eq.'.' .or.
     &       linebuf(33:33).eq.'.' .or.
     &       linebuf(34:34).eq.'.') go to 401

      if (iflag.eq.0) goto 401
         read (linebuf,301) stname,phinfo(1:4),iyr,imon,idy,ihr,imn,
     &            isec,idelkm
!         print 301, stname,phinfo(1:4),iyr,imon,idy,ihr,imn,
!     &            isec,idelkm
301      format (a12,1x,a4,i4,4i2,i5,40x,i4)

         sec=float(isec)/100.
         delkm=float(idelkm)/10.

         call DT_TIMEDIF(qyr,qmon,qdy,qhr,qmn,qsc,
     &                   iyr,imon,idy,ihr,imn,sec,t0)
         if (abs(t0).gt.1000.) then
            print *,'Timing problem, t0 = ',t0
            print *,qyr,qmon,qdy,qhr,qmn,qsc
            print *,iyr,imon,idy,ihr,imn,sec
            print *,linebuf
            stop
         endif

         tt=real(t0)-qt0

         delta=delkm
         stname2='            '
         stname2(1:2)=stname(6:7)
         stname2(4:7)=stname(1:4)
         stname2(9:12)=stname(9:12)
         comp=stname(10:12)
         iph(1:1)=phinfo(2:2)     !e.g., P or S
         onset=phinfo(1:1)          !e.g., I or E
         acc=phinfo(4:4)            !quality
         pol=phinfo(3:3)            !e.g., U or D


      if (delta .gt. delmax.or.delta.le.0.0) goto 401
      if (tt .le. 0.0) goto 401
      if (iph(1:1).ne.'P'.and.iph(1:1).ne.'p'.and.
     &    iph(1:1).ne.'S'.and.iph(1:1).ne.'s') goto 401
      
!      if (onset.ne.'I'.and.
      if(stname(1:7).eq.stnameold.and.
     &   phaseold.eq.iph(1:1)) goto 401

      ipick=ipick+1
      stname0(ipick)=stname2
      pickph(ipick)=iph
      pick(ipick)=tt
      
!      print *, ipick,stname2,iph,tt,iyr,imon,idy,ihr,imn,
!     &         isec,delta

      stnameold=stname
      phaseold=iph
      else if (linebuf(1:15) .eq. '               ') then  !End of file
      read (linebuf,*) idcusp0(iq)
      endif
401   continue
101   close(11)
      nqtot=iq
      npickall=ipick
      ipick2(nqtot)=npickall
      ipick1(1)=1
      do iq=2,nqtot
      ipick1(iq)=ipick2(iq-1)+1
      end do
      
      end
!
!
!-----------------------------------------------------------------------
!
!  subroutine(02)
!
!  subroutine GET_STANUM is utility routine to implement
!  an automatic numbering scheme for station names
!  we check only first 7 characters which should have
!  network id and station name, but not channel
!
!    menu = 1
!        Input:  stname  =  12 character station name
!        Return: ista    =  unique number for station name
!                istamax =  current maximum station number
!
!    menu = 2
!        Input:  ista    =  station number
!        Return: stname  =  12 character station name
!
!-----------------------------------------------------------------------
!
      subroutine GET_STANUM(stname,ista,istamax,menu)
      implicit none


      integer imax
      integer menu
      integer i
      integer ista
      integer istamax
      

      character*12 stname
      character*12 stlist(5000)


      save stlist,imax
!
!-----------------------------------------------------------------------
!
      if (menu.eq.1) then
         do 10 i=1,imax
            if (stlist(i)(1:7).eq.stname(1:7)) then
               ista=i
               return
            end if
10       continue
         imax=imax+1
         if (imax.lt.1) imax=1
         ista=imax
         istamax=imax
         stlist(ista)=stname
      else
         if (ista.lt.1.or.ista.gt.imax) then
            print *,'***ERROR in GET_STANUM: ',ista,imax,menu
            stop
         end if
         stname=stlist(ista)
      end if
      return
      end
!
!
!-----------------------------------------------------------------------
!
!  subroutine(03)
!
!  GETSTAT_STLIST gets station location from stlist file generated
!  from sactogfs_stp program.  It uses only first 7 characters of
!  12 character station name.
!
!  Input: stfile --- station list file (unsorted)
!         snam --- station code
!
!  Output: flat --- station latitude (999. if no match) 
!          flon --- station longitude (999. if no match)
!          felev --- station elevation (m) (999. if no match)
!
!-----------------------------------------------------------------------
!
      subroutine GETSTAT_STLIST(stfile,snam,flat,flon,felev,stterm)
      implicit none


      integer i
      integer nsta


      real felev
      real flat
      real flon
      real selev(5000)
      real slat(5000)
      real slon(5000)
      real stterm1(5000)
      real stterm2(5000)
      real stterm(2)


      character*12 snam
      character*100 stfile
      character*12 stname(5000)


      logical firstcall
      save firstcall,stname,slat,slon,selev,stterm1,stterm2,nsta
      data firstcall/.true./
!
!-----------------------------------------------------------------------
!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stfile,status='old')
         do i=1,5000
            read (19,11,end=12) stname(i),slat(i),slon(i),selev(i),
     &                          stterm1(i),stterm2(i)
11          format (a12,f10.5,f12.5,f10.1,2f8.2)
         enddo
         print *,'***Warning:  station file has more than 5000 lines'
         i=5001
12       nsta=i-1
         close (19)
         print *,'Station locations read.  nsta = ',nsta
      end if
      
      do 30 i=1,nsta
         if (snam(1:7) .eq. stname(i)(1:7)) then
            flat=slat(i)
            flon=slon(i)
            felev=selev(i)
            stterm(1)=stterm1(i)
            stterm(2)=stterm2(i)
            return
         endif
30    continue
      print *,'***station not found ',snam
      flat=999.
      flon=999.
      felev=999.
      stterm(1)=0.0
      stterm(2)=0.0
      return
      end
!
!
!-----------------------------------------------------------------------
!
!  subroutine(04)
!
!  GET_STPPHASE extracts quake info and picks from STP PHASE DATA file
!  Note:  only extracts events within qlatmin, etc. limits
!         and only gets picks within delmax distance
!         This version gets P or S picks only
!         
!         Only takes first pick when more than one pick
!         for the same station is present
!
!         permits getting only random fraction of events
!
!   inputs and outputs are the same as in GETBED3_R subroutine.
!
!-----------------------------------------------------------------------
!
      subroutine GET_STPPHASE(ipha,infile,
     &    qlatmin,qlatmax,qlonmin,qlonmax,
     &    delmax,frac,nqtot,npickall,idcusp0,tsec0,qlat0,qlon0,qdep0,
     &    qmag0,ipick1,ipick2,stname0,pick,pickph)
      implicit none

      integer nq0                !maximum number of events
      parameter (nq0=500000)
      integer npick0             !maximum number of picks
      parameter (npick0=13000000)

      integer ipha
      integer*4 cuspid
      integer i
      integer iline
      integer ipick
      integer iq
      integer npickall        !total number of picks (output)
      integer nqmax
      integer nqtot           !total number of events (output)
      integer qyr
      integer qmon
      integer qdy
      integer qhr
      integer qmn
      integer*4 idcusp0(nq0)  !array with event identification number (output)
      integer ipick1(nq0)     !array with first pick index for each event (output)
      integer ipick2(nq0)     !array with last pick index for each event (output)


      real qlatmin,qlatmax,qlonmin,qlonmax   !location region window (input)
      real delmax       !cutoff distance between stations and events (input)
!
!     variables for quake info
      real frac          !random fraction of events to read from phase data file (input)
      real fran          
      real qdep
      real qlat
      real qlon
      real qmag
      real qsc
      real RAND
      real*8 t0
      real pick(npick0)  !array with travel time pick for each pick (output)
      real qdep0(nq0)    !array with event depth (output)
      real qlat0(nq0)    !array with evetn latitude (output)
      real qlon0(nq0)    !array with event longitude (output)
      real qmag0(nq0)    !array with event magnitude (output)
      real*8 tsec0(nq0)  !event time (seconds since 1600) (output)
!
!     variables for station pick info
      real delta
      real tt
!
      character*1 evtype
      character*100 infile         !input phase data file name (input)
      character*2 iph
      character*100 linebuf        
      character*12 stname
      character*12 stname2
      character*2 pickph(npick0)   !array with phase id for each pick (output)
      character*12 stname0(npick0)   !array with station id (12 characters)
!                                     including network id, station name and channel.  (output)  
      character*12 stnameold
      character phaseold
!
!-----------------------------------------------------------------------
!     
      do iq=1,nq0
      ipick1(iq)=0
      ipick2(iq)=0
      enddo
      
      nqmax=999999
      open(11,file=infile,status='old')
      iq=0
      ipick=0
      
      do 401 iline=1,9999999
      read(11,'(a)',end=101) linebuf
      if (linebuf(1:1).eq.'#') go to 101

      if ((ipha.eq.2.and.linebuf(77:79).ne.'   ').or.             !OLD STP
     &    (ipha.eq.4.and.linebuf(1:1).eq.' ')) then               !NEW STP
!                                                                  QUAKE LINE

      if (iq .ne. 0) ipick2(iq)=ipick
      
      read (linebuf,201) cuspid,evtype,qyr,qmon,qdy,qhr,qmn,qsc,
     &                   qlat,qlon,qdep,qmag
!      print 201, cuspid,evtype,qyr,qmon,qdy,qhr,qmn,qsc,
!     &           qlat,qlon,qdep,qmag
201   format (i10,x,a1,2x,i4,x,i2,x,i2,x,i2,x,i2,x,f6.3,
     &        x,f9.4,x,f11.4,x,f6.2,x,f5.2)
      
      if (qlat.lt.qlatmin.or.qlat.gt.qlatmax)  goto 401
      if (qlon.lt.qlonmin.or.qlon.gt.qlonmax)  goto 401
      
      fran=RAND(0)
      if (fran.gt.frac) go to 401
      
      iq=iq+1
      idcusp0(iq)=cuspid
      qlat0(iq)=qlat
      qlon0(iq)=qlon
      qdep0(iq)=qdep
      qmag0(iq)=qmag
      call DT_TIMEDIF(1600,1,0,0,0,0.,qyr,qmon,qdy,qhr,qmn,qsc,t0)
      tsec0(iq)=t0
      
      stnameold='123456789012'
      phaseold=' '

      else                                      !PICK LINE

         if (ipha.eq.2) then !OLD
            read (linebuf,301) stname,iph,delta,tt
!            print 301, stname,iph,delta,tt
301         format(1x,a12,29x,a2,11x,f7.2,f8.3)

            if (delta .gt. delmax) goto 401
            if (tt.le.0.0) go to 401
            if (iph(1:1).ne.'P'.and.
     &          iph(1:1).ne.'p'.and.
     &          iph(1:1).ne.'S'.and.
     &          iph(1:1).ne.'s') goto 401

            if (stname(1:7).eq.stnameold(1:7) .and.
     &          iph(1:1).eq.phaseold) go to 401

               stname2=stname
               if(stname2(4:4).eq.' ') then
               do i=5,8
               stname2(i-1:i-1) = stname2(i:i)
               end do
               end if

         else if (ipha.eq.4.and.linebuf(44:44) .eq. ' ') then  !NEW
            read (linebuf,302)  stname(1:2),stname(3:12),iph,delta,tt
!            print 302, stname(1:2),stname(3:12),iph,delta,tt
302         format(1x,a2,2x,a10,28x,a2,11x,f7.2,f8.3)
            if (delta .gt. delmax) goto 401
            if (tt.le.0.0) go to 401
            if (iph(2:2).ne.'P'.and.
     &          iph(2:2).ne.'p'.and.
     &          iph(2:2).ne.'S'.and.
     &          iph(2:2).ne.'s') goto 401

            if (stname(1:7).eq.stnameold(1:7) .and.
     &          iph(2:2).eq.phaseold) go to 401

               stname2=stname
               if(stname2(4:4).eq.' ') then
               do i=5,8
               stname2(i-1:i-1) = stname2(i:i)
               end do
               end if

             iph(1:1)=iph(2:2)
             iph(2:2)=' '

        endif

      ipick=ipick+1
      stname0(ipick)=stname2
      pickph(ipick)=iph
      pick(ipick)=tt
      
      stnameold(1:7)=stname(1:7)
      phaseold=iph

      endif
401   continue
101   close(11)
      nqtot=iq
      npickall=ipick
      ipick2(nqtot)=npickall
      ipick1(1)=1
      do iq=2,nqtot
      ipick1(iq)=ipick2(iq-1)+1
      end do
      
      end
!
!
!-----------------------------------------------------------------------
!
!  subroutine(05)
!
!  GET_TTS_FAST obtains a travel time for a seismic phase
!  at a specified range and earthquake depth by interpolating
!  from a file containing a table of travel times. 
!
!  NOTE: Assumes x and d are evenly spaced
!
!  Input:  phase --- travel time table file
!          ip --- index number for phase (1 for P, 2 for S)
!          del --- source-receiver range
!          qdep --- event depth
!
!  Output:  tt --- travel time
!           iflag --- -1 if outside depth range
!                      0 for interpolation
!                      1 for extrapolation in range
!
!-----------------------------------------------------------------------
!
      subroutine GET_TTS_FAST(phase,ip,del,qdep,tt,iflag)
      implicit none


      integer nx0
      parameter (nx0=1001)
      integer nd0
      parameter (nd0=480)
      integer id
      integer id1
      integer id2
      integer iflag
      integer ip
      integer ix
      integer ix1
      integer ix2
      integer ixbest1
      integer ixbest2
      integer nd(2)
      integer nx(2)


      real del
      real dfrac
      real qdep
      real t1
      real t2
      real tt
      real tt1
      real tt2
      real xfrac
      real xfrac1
      real xfrac2
      real xoff
      real xoffmin1
      real xoffmin2
      real d(nd0,2)
      real dd(2)
      real dx(2)
      real t(nx0,nd0,2)
      real x(nx0,2)


      character*100 linebuf
      character*100 phase
      character*100 phaseold(2)

      save t,x,d,phaseold,nx,nd,dx,dd
!
!-----------------------------------------------------------------------
!
! read file if new phase file is specified
!
      if (phase.ne.phaseold(ip)) then
         print *,'reading phase file name: ',phase(1:40)
         open (3,file=phase,status='old',err=990)
         read (3,'(a40)') linebuf   !ignore first line
         read (3,*) nx(ip),nd(ip)
         if (nx(ip).gt.nx0) then
             print *,'***GET_TTS nx truncated ',nx(ip),nx0
             nx(ip)=nx0
         end if
         if (nd(ip).gt.nd0) then
             print *,'***GET_TTS nd truncated ',nd(id),nd0
             nd(ip)=nd0
         end if
         read (3,*) (d(id,ip),id=1,nd(ip))
         do 20 ix=1,nx(ip)
         read (3,*) x(ix,ip),(t(ix,id,ip),id=1,nd(ip))
20       continue
         close (3)
         dx(ip)=x(2,ip)-x(1,ip)
         dd(ip)=d(2,ip)-d(1,ip)
         print *,'nx,nd,dx,dd= ',nx(ip),nd(ip),dx(ip),dd(ip)
      end if
      phaseold(ip)=phase
!
! check if outside depth range
      if (qdep.lt.d(1,ip).or.qdep.gt.d(nd(ip),ip)) then
         iflag=-1
         tt=999
         return
      end if
!
! first check to see if interpolation alone will work
      id1=1+int((qdep-d(1,ip))/dd(ip))
      id2=id1+1
      ix1=1+int((del-x(1,ip))/dx(ip))
      ix2=ix1+1

      if (ix1.ge.nx0) then
         print *,'***Error in GET_TTS_FAST'
         print *,'del,qdep = ',del,qdep,ix1,nx0,x(1,ip),dx(ip)
         stop
      end if
      
37    if (t(ix1,id1,ip).eq.0.) go to 50
      if (t(ix1,id2,ip).eq.0.) go to 50
      if (t(ix2,id1,ip).eq.0.) go to 50
      if (t(ix2,id2,ip).eq.0.) go to 50
!
      if (t(ix1,id1,ip).eq.-9.9900.or.
     &    t(ix1,id2,ip).eq.-9.9900.or.
     &    t(ix2,id1,ip).eq.-9.9900.or.
     &    t(ix2,id2,ip).eq.-9.9900) then
          iflag=-1
          tt=999
          return
      end if
!
      if (x(ix2,ip).lt.del) go to 50
      iflag=0
      xfrac=(del-x(ix1,ip))/(x(ix2,ip)-x(ix1,ip))
      t1=t(ix1,id1,ip)+xfrac*(t(ix2,id1,ip)-t(ix1,id1,ip))
      t2=t(ix1,id2,ip)+xfrac*(t(ix2,id2,ip)-t(ix1,id2,ip))
      dfrac=(qdep-d(id1,ip))/(d(id2,ip)-d(id1,ip))
      tt=t1+dfrac*(t2-t1)
      return
!        
! extrapolate to get tt
50    iflag=1
      xoffmin1=999.
      xoffmin2=999.
      ixbest1=999
      ixbest2=999
      do 60 ix=2,nx(ip)
         if (t(ix-1,id1,ip).eq.0) go to 55
         if (t(ix,id1,ip).eq.0) go to 55
         xoff=abs((x(ix-1,ip)+x(ix,ip))/2.-del)
         if (xoff.lt.xoffmin1) then
            xoffmin1=xoff
            ixbest1=ix
         end if
55       if (t(ix-1,id2,ip).eq.0) go to 60
         if (t(ix,id2,ip).eq.0) go to 60
         xoff=abs((x(ix-1,ip)+x(ix,ip))/2.-del)
         if (xoff.lt.xoffmin2) then
            xoffmin2=xoff
            ixbest2=ix
         end if
60    continue
      if (ixbest1.eq.999.or.ixbest2.eq.999) then
         iflag=-1
         tt=999
         return
      end if
      
      xfrac1=(del-x(ixbest1-1,ip))/(x(ixbest1,ip)-x(ixbest1-1,ip))
      t1=t(ixbest1-1,id1,ip)
      t2=t(ixbest1,id1,ip)
      tt1=t1+xfrac1*(t2-t1)
      
      xfrac2=(del-x(ixbest2-1,ip))/(x(ixbest2,ip)-x(ixbest2-1,ip))
      t1=t(ixbest2-1,id2,ip)
      t2=t(ixbest2,id2,ip)
      tt2=t1+xfrac2*(t2-t1)
      
      dfrac=(qdep-d(id1,ip))/(d(id2,ip)-d(id1,ip))
      tt=tt1+dfrac*(tt2-tt1)
      
      go to 999
!      
990   print *,'*** phase file not found: ',phase
      stop
999   return
      end
!
!
!-----------------------------------------------------------------------
!
!   subroutine(06)
!
!   GRIDLOC calls LOCATE3DF to locate quakes using grid-search method.
! 
!   Input:  tname --- travel time table file
!           nq --- number of events
!           idph --- array with phase id for each pick
!           idsta --- array with station number for each pick
!           npick_min --- minimum number of picks per event
!           ipick1 --- array with the first pick index for each event
!           ipick2 --- array with the last pick index for each event
!           pick --- array with travel time picks
!           term --- array with station terms for each pick
!           qlat --- array with quake latitude
!           qlon --- array with quake longitude
!           qdep --- array with quake depth
!           qorg --- array with quake origin time
!           qdepref --- reference starting depth
!           ilocfix --- (0) normal or (2) fix loc kluge
!           idepfix --- (0) catalog depth  or  (1) qdepref given by user
!           inorm --- (1) L1 norm ,(2) L2 norm or (3) ROBUST MEAN
!           nit --- number of iterations for grid-search 
!           fracsk --- shrinking fraction for grid search (e.g. 0.67)
!           stlat --- array with station latitude
!           stlon --- array with station longitude
!      
!   Output:  qlat2 --- array with new quake latitude
!            qlon2 --- array with new quake longitude
!            qdep2 --- array with new quake depth
!            qorg2 --- array with new quake origin time
!            nqqual --- array with earthquake location quality
!                       if -1, bad location, should not be used in station terms calculation
!            resid --- array with travel time residuals for each pick
!            rms --- array with RMS residuals for each event
!            rmed --- array with Median Absolute Value of residuals for each event
! 
!-----------------------------------------------------------------------
!
      subroutine GRIDLOC(tname,nq,idph,idsta,npick_min,ipick1,ipick2,
     &    pick,term,qlat,qlon,qdep,qorg,ilocfix,idepfix,inorm,
     &    nit,fracsk,stlat,stlon,qlat2,qlon2,qdep2,qorg2,
     &    nqqual,resid,rms,rmed)
      implicit none

      integer nq0       !maximum number of events
      parameter (nq0=500000)
      integer npick0    !maximum number of picks
      parameter (npick0=13000000)
      integer nsta0     !maximum number of stations
      parameter (nsta0=3000)
      integer mp        !maximum number of picks per event
      parameter (mp=1000)

      integer i
      integer i1
      integer i2
      integer i4
      integer idepfix
      integer ilocfix
      integer inorm
      integer iq
      integer iter
      integer j
      integer n
      integer nit
      integer npick_min
      integer nq
      integer*2 idph(npick0)
      integer*2 idsta(npick0)
      integer ip(mp)
      integer ipick1(nq0)
      integer ipick2(nq0)
      integer nqqual(npick0)

      real boxwid
      real fracsk
      real pick(npick0)
      real qdep(nq0)
      real qlat(nq0)
      real qlon(nq0)
      real qorg(nq0)
      real qdep2(nq0)
      real qlat2(nq0)
      real qlon2(nq0)
      real qorg2(nq0)    
      real resid(npick0)
      real resol(nq0)
      real rmed(nq0)
      real rms(nq0)
      real slat(mp)
      real slon(mp)
      real stlat(nsta0)
      real stlon(nsta0)
      real term(npick0)
      
      character*100 tname(2)

      logical firstcall
      save firstcall,iter
      data firstcall/.true./
!
!-----------------------------------------------------------------------
!     
      if (firstcall) then
         iter=1
         boxwid=20.
         firstcall=.false.
      else
         iter=iter+1
!     allow to change boxwid
         boxwid=20.0
      end if
      print *,'iter=',iter,'  boxwid=',boxwid

      do 200 iq=1,nq
         
         i1=ipick1(iq)
         i2=ipick2(iq)
         n=i2-i1+1
         if (n.lt.npick_min) then
         qlat2(iq)=qlat(iq)
         qlon2(iq)=qlon(iq)
         qdep2(iq)=qdep(iq)
         qorg2(iq)=qorg(iq)
!         if (iter.eq.1) print *, 'event ',iq,' was skipped!'
         go to 200     !***skip if too few picks
         end if
         
         do 150 i=i1,i2
            j=i-i1+1
            i4=idsta(i)  !convert to I*4
            
            if (i4.eq.0) then
               print *,'***i4 = 0: ',iq,i,j,idsta(i)
               stop
            end if
            
            slat(j)=stlat(i4)
            slon(j)=stlon(i4)
            
            if (slat(j).eq.999.) go to 150
            
            ip(j)=idph(i)
150      continue


         if (qlon(iq).eq.0.) then
            print *,'***PROB1'
            print *,iq,qlat(iq),qlon(iq)
         end if
         
         call LOCATE3DF(qlat(iq),qlon(iq),qdep(iq),n,
     &      pick(i1),term(i1),ip,slat,slon,tname,boxwid,nit,ilocfix,
     &      inorm,fracsk,qlat2(iq),qlon2(iq),qdep2(iq),qorg2(iq),
     &      nqqual(i1),resid(i1),rms(iq),rmed(iq),resol(iq))
!
          if (iq .eq. 1 .or. mod(iq,500) .eq. 0) then
          print *,iq,n,rms(iq),rmed(iq)
          print *,qlat(iq),qlon(iq),qdep(iq),qorg(iq)
          print *,qlat2(iq),qlon2(iq),qdep2(iq),qorg2(iq)
          end if
200   continue
      
      return
      end
!
!
!-----------------------------------------------------------------------
!
!     subroutine(07)
!
!-----------------------------------------------------------------------
!
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      IMPLICIT NONE


      INTEGER J
      INTEGER L
      INTEGER I
      INTEGER IR
      INTEGER INDXT
      INTEGER N
      INTEGER INDX(N)


      REAL Q
      REAL ARRIN(N)
!
!-----------------------------------------------------------------------
!
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
!
!
!-----------------------------------------------------------------------
!
! subroutine(08)
!
! LOCATE3DF performs grid-search event location in three dimensions
!           This version (3DF) is for local event location
!           permits fixed locations to compute station terms
!
!  Input:  qlat0 --- initial event latitude
!          qlon0 --- initial event longitude
!          qdep0 --- initial event depth
!          npick --- number of picks
!          tt --- array with travel times (s)
!          term --- array with station terms
!          ip --- array with phase index numbers (1 for P or 2 for S)
!          slat --- array with station latitude (999. to skip)
!          slon --- array with station longitude (999. to skip)
!          tname --- phase data file
!          boxwid --- starting shrinking box width (km)
!          nit --- number of iterations to perform
!          ifix --- 0 for normal, 1 to fix locations
!          inorm --- 1 for L1, 2 for L2 and 3 for ROBUST MEAN
!          fracsk --- fraction of shrinking per iteration
!
!  Output:  qlat --- best-fitting latitude
!           qlon --- best-fitting longitude
!           qdep --- best-fitting depth
!           qorg --- best-fitting origin time
!           nqqual --- array with earthquake location quality
!                      -1, bad location, should not be used in station terms calculation
!           resid --- array with residuals
!           rms --- rms residual
!           rmed --- median absolute value residual
!           resol --- nominal resolution (m)
!
!-----------------------------------------------------------------------
!
      subroutine LOCATE3DF(qlat0,qlon0,qdep0,npick,tt,term,ip,slat,slon,
     &                  tname,boxwid,nit,ifix,inorm,fracsk,qlat,qlon,
     &                  qdep,qorg,nqqual,resid,rms,rmed,resol)
      implicit none


      integer i
      integer ifix
      integer iflag
      integer ii
      integer inorm
      integer it
      integer ix
      integer iy
      integer iz
      integer nit
      integer niter
      integer npick
      integer npick2
      integer ip(npick)
      integer nqqual(npick)


      real boxwid
      real cosqlat
      real ddep
      real degkm
      real degrad
      real delkm
      real dlat
      real dlon
      real dx
      real dy
      real fdep
      real fdep0
      real fdepbest
      real fit
      real fit2
      real fitbest
      real flat
      real flat0
      real flatbest
      real flon
      real flon0
      real flonbest
      real fracsk
      real qdep
      real qlat
      real qlon
      real qorg
      real qdep0
      real qlat0
      real qlon0
      real residmed
      real resol
      real rmed
      real rms
      real tbest
      real tsec
      real xgap
      real rabs(3000)
      real resid(npick)
      real resid2(npick)
      real slat(npick)
      real slon(npick)
      real term(npick)
      real tt(npick)


      character*100 tname(2)
!
!-----------------------------------------------------------------------
!
      degrad=180./3.1415927
      degkm=111.19493

      dlat=0.5*boxwid/degkm
      cosqlat=cos(qlat0/degrad)
      dlon=dlat/cosqlat
!
      ddep=0.5*boxwid
!
!      ddep=0.5*boxwid*0.5
!
!      ddep=0.5*boxwid*2.0
!
      flat0=qlat0
      flon0=qlon0
      fdep0=qdep0

      niter=nit
      if (ifix.eq.1) niter=1

      do 100 it=1,niter

      fitbest=9.e20
      do 60 iy=-1,1
         flat=flat0+dlat*float(iy)
         do 50 ix=-1,1
            flon=flon0+dlon*float(ix)
            do 40 iz=-1,1
               fdep=fdep0+ddep*float(iz)
               if (fdep.lt.0.) fdep=0.

               if (ifix.eq.1) then
                  flat=flat0
                  flon=flon0
                  fdep=fdep0
               endif

               ii=0                           
               do 20 i=1,npick
                  if (slat(i).eq.999.) go to 20
                  dy=flat-slat(i)
                  dx=(flon-slon(i))*cosqlat
                  delkm=sqrt(dx**2+dy**2)*degkm

                  if (delkm.ge.500.) then
!                     print *,'***PROBLEM in LOCATE3DF '
!                     print *,'ipick=',i
!                     print *,qlat0,qlon0,qdep0
!                     print *,slat(i),slon(i)
!                     stop
                      go to 20
                  end if

                  call GET_TTS_FAST(tname(ip(i)),ip(i),delkm,
     &                              fdep,tsec,iflag)
!                  if (iflag .eq. -1)  goto 20
                  if(iflag.eq.-1) then
                     resid(i)=-999.
                  else
                     ii=ii+1
                     resid2(ii)=tt(i)-term(i)-tsec
                     resid(i)=resid2(ii)
                  end if
20             continue
               npick2=ii

               if (inorm .eq. 1) then         !L1 NORM
                  call MEDIAN(resid2,npick2,residmed)
               else if (inorm .eq. 2) then    !L2 NORM
                  call MEAN(resid2,npick2,residmed)
               else if (inorm .eq. 3) then    !ROBUST MEAN
                  xgap=0.1
                  call ROBOMEAN2(resid2,npick2,xgap,10,residmed,fit2)
               end if

               fit=0.
               if (inorm .eq. 1) then
               do 30 i=1,npick2
                  fit=fit+abs(resid2(i)-residmed)
30             continue
               else if (inorm .eq. 2) then
               do i=1,npick2
                  fit=fit+(resid2(i)-residmed)**2
               end do
               else if (inorm .eq. 3) then
                  fit=fit2
               end if

               if (fit.lt.fitbest.and.npick2.gt.0) then
                  fitbest=fit
                  flatbest=flat
                  flonbest=flon
                  fdepbest=fdep
                  tbest=residmed
               end if

40          continue
50       continue
60    continue

      flat0=flatbest
      flon0=flonbest
      fdep0=fdepbest
      dlat=dlat*fracsk      !shrink box by fracsk each iteration
      dlon=dlon*fracsk
      ddep=ddep*fracsk

100   continue

      qlat=flat0
      qlon=flon0
      qdep=fdep0
      qorg=tbest
      resol=(dlat/fracsk)*degkm

      ii=0     
      rms=0.
      do 120 i=1,npick
         if (slat(i).eq.999.) go to 120
         dy=qlat-slat(i)
         dx=(qlon-slon(i))*cosqlat
         delkm=sqrt(dx**2+dy**2)*degkm
         if (delkm.ge.500.) goto 120
         call GET_TTS_FAST(tname(ip(i)),ip(i),delkm,
     &                     qdep,tsec,iflag)
!         if (iflag .eq. -1)  goto 120
         if(iflag.eq.-1) then
            resid(i)=-999.
         else
            ii=ii+1
            resid2(ii)=(tt(i)-term(i)-qorg)-tsec
            resid(i)=resid2(ii)
            rms=rms+resid2(ii)**2
            rabs(ii)=abs(resid2(ii))
         end if
120   continue 
      npick2=ii
      if (npick2 .ne. 0) rms=sqrt(rms/float(npick2))
      if (npick2 .eq. 0 .and. rms .eq. 0.0) then
         do i=1,npick
            nqqual(i)=-1
         end do   
      end if   
      if (inorm .eq. 1) then         !L1 NORM
         call MEDIAN(rabs,npick2,rmed)
      else if (inorm .eq. 2) then    !L2 NORM
         call MEAN(resid2,npick2,rmed)
      else if (inorm .eq. 3) then    !ROBUST MEAN
         xgap=0.1 
         call ROBOMEAN2(resid2,npick2,xgap,10,rmed,fit2)
      end if   

      return
      end
!
!     
!-----------------------------------------------------------------------
!
!     subroutine(09)
!
!-----------------------------------------------------------------------
!
      subroutine MEAN(x,n,xmean)
      implicit none

      integer i
      integer n

      real sum
      real xmean
      real x(n)
!
!
      if (n.eq.0) then
         xmean=0.
         return
      end if
      sum=0.
      do 10 i=1,n
         sum=sum+x(i)
10    continue
      xmean=sum/float(n)
      return
      end
!
!
!-----------------------------------------------------------------------
!
!     subroutine(10)
!
!-----------------------------------------------------------------------
!
      subroutine MEDIAN(x,n,xmed)
      implicit none

      integer n
      integer n2

      real xmed
      real x(n)
!
!
      if (n.eq.0) then
         xmed=0.
         return
      else if (n.eq.1) then
         xmed=x(1)
         return
      end if
      call SORT(n,x)
      n2=n/2
      if(2*n2.eq.n)then
        xmed=0.5*(x(n2)+x(n2+1))
      else
        xmed=x(n2+1)
      endif
      return
      end
!
!
!-----------------------------------------------------------------------
!
!     subroutine(11)
!
!-----------------------------------------------------------------------
!
      subroutine RDBLK(a,n,nun,iend)
      implicit none

      integer n
      integer nun
      integer iend
      integer*2 a(n)
!
      iend=0
      read (nun,end=99) a
      return
99    iend=1
      return
      end
!
!
!-----------------------------------------------------------------------
!
!     subroutine(12)
!
!-----------------------------------------------------------------------
!
      subroutine ROBOMEAN2(x,n,xgap,nit,xmean,fit2)
      parameter (nmax=10000000)
      real x(n),xw(nmax),fit2

      if (n.eq.0) then
         xmean=0.
         return
      else if (n.gt.nmax) then
         print *,'***ROBOMEAN error, n,nmax = ',n,nmax
         stop
      end if
      do i=1,n
         xw(i)=1.
      enddo
      do it=1,nit
         sum=0.
         sumw=0.
         do i=1,n
            sum=sum+x(i)*xw(i)
            sumw=sumw+xw(i)
         enddo
         xmean=sum/sumw
         if (it.ne.nit) then
            do i=1,n
               d=xgap+abs(x(i)-xmean)
               xw(i)=1./d
            enddo
         endif
      enddo

      fit2=0.0
      do i=1,n
         fit2=fit2+xw(i)*(x(i)-xmean)**2
      end do

      return
      end
!
!
!-----------------------------------------------------------------------
!
!     subroutine(13)
!
!-----------------------------------------------------------------------
!
      SUBROUTINE SORT(n,arr)
      implicit none


      INTEGER i
      INTEGER ir
      INTEGER j
      INTEGER k
      INTEGER l
      INTEGER M
      INTEGER n
      INTEGER NSTACK
      PARAMETER (M=7, NSTACK=50)
      INTEGER jstack
      INTEGER istack(NSTACK)


      REAL a
      REAL temp
      REAL arr(n)
!
!-----------------------------------------------------------------------
!      
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then  !Insertion sort when subarray small enough.
         do j=l+1,ir
               a=arr(j)
               do i=j-1,l,-1
                  if(arr(i).le.a)goto 2
                  arr(i+1)=arr(i)
               enddo
               i=l-1
2              arr(i+1)=a
         enddo
         if(jstack.eq.0)return
         ir=istack(jstack)   !Pop stack and begin a new round of partitioning.
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2   !Choose median of left, center, and right elements 
                      !as partitioning element a. Also rearrange 
                      !so that a(l) <= a(l+1) <= a(ir).
         temp=arr(k)
         arr(k)=arr(l+1)
         arr(l+1)=temp
         if(arr(l).gt.arr(ir))then
           temp=arr(l)
           arr(l)=arr(ir)
           arr(ir)=temp
         endif
         if(arr(l+1).gt.arr(ir))then
           temp=arr(l+1)
           arr(l+1)=arr(ir)
           arr(ir)=temp
         endif
         if(arr(l).gt.arr(l+1))then
           temp=arr(l)
           arr(l)=arr(l+1)
           arr(l+1)=temp
         endif      
         i=l+1      !Initialize pointers for partitioning.
         j=ir
         a=arr(l+1) !Partitioning element.
3        continue   !Beginning of innermost loop.
         i=i+1      !Scan up to find element > a.
         if(arr(i).lt.a)goto 3
4        continue   
         j=j-1      !Scan down to find element < a.
         if(arr(j).gt.a)goto 4
         if(j.lt.i)goto 5  !Pointers crossed. Exit with partitioning complete.
         temp=arr(i)       !Exchange elements.
         arr(i)=arr(j)
         arr(j)=temp       
         goto 3            !End of innermost loop.
5        arr(l+1)=arr(j)   !Insert partitioning element.
         arr(j)=a
         jstack=jstack+2
!        Push pointers to larger subarray on stack, 
!        process smaller subarray immediately.
         if(jstack.gt.NSTACK)pause  ! 'NSTACK too small in sort'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      END
!
!
!-----------------------------------------------------------------------
!
! subroutine(14)
! 
! SSST_DISTF computes station terms from residual list using closest events
!            this version is for local events (flat earth)
!
!  Input: inorm --- 1 for L1 norm, 2 for L2 norm or 3 for ROBUST MEAN
!         resid --- array with residuals
!         isin --- station number array
!         ipin --- array with phase numbers
!         nr --- number of residuals
!         nstpick --- array with number of picks per station+phase
!         nq --- number of events 
!         qlat --- array of event latitude
!         qlon --- array of event longitude
!         qdep --- array of event depth
!         nqqual --- array with earthquake location quality
!                   -1, bad location, should not be used in station terms calculation
!         ipick1 --- array with first pick index for each event
!         ipick2 --- array with last pick index for each event
!         nmed --- number of stations to compute median over
!         dismax --- cutoff distance 
!
!  Outout:  term --- array of station terms 
!
!-----------------------------------------------------------------------
!
      subroutine SSST_DISTF(inorm,resid,isin,ipin,nr,nstpick,nq,
     &              qlat,qlon,qdep,nqqual,ipick1,ipick2,
     &              nmed,dismax,term)
      implicit none


      integer nmax
      parameter (nmax=500000)
      integer nq0
      parameter (nq0=500000)
      integer nsta0
      parameter (nsta0=3000)
      integer i
      integer ii
      integer imed
      integer inorm
      integer ipick
      integer ip_targ
      integer iq
      integer iq2
      integer iscr1
      integer iscr2
      integer iscr3
      integer iscr4
      integer is_targ
      integer j
      integer jpick
      integer jq
      integer jq2
      integer njq
      integer nmed
      integer nq
      integer nr
      integer indexjq(nq0)
      integer indist(nq0)
      integer ipick1(nq)
      integer ipick2(nq)
      integer*2 ipin(nr)
      integer*2 isin(nr)
      integer nqqual(nr)
      integer nstpick(nsta0,2)


      real aspect
      real degoff
      real degrad
      real dismax
      real distavg
      real dx
      real dy
      real fit2
      real qlatmax
      real qlatmin
      real qlonmax
      real qlonmin
      real rmed
      real xgap
      real dist(nq0)
      real qdep(nq)
      real qlat(nq)
      real qlon(nq)
      real resid(nr)
      real rr(nmax)
      real term(nr)
!
!-----------------------------------------------------------------------
!
      degrad=180./3.1415927

      degoff=dismax/111.19     !this can greatly speed calculation

      do 30 i=1,nr
         term(i)=0.
30    continue

      do 200 iq=1,nq

         qlatmin=qlat(iq)-degoff
         qlatmax=qlat(iq)+degoff
         aspect=cos(qlat(iq)/degrad)
         qlonmin=qlon(iq)-degoff/aspect
         qlonmax=qlon(iq)+degoff/aspect

         distavg=0.
         jq=0
         do 60 j=1,nq
            if (qlat(j).lt.qlatmin) go to 60
            if (qlat(j).gt.qlatmax) go to 60
            if (qlon(j).lt.qlonmin) go to 60
            if (qlon(j).gt.qlonmax) go to 60
            if (j.eq.iq) go to 60
            jq=jq+1
            dy=(qlat(iq)-qlat(j))*111.19
            dx=(qlon(iq)-qlon(j))*aspect*111.19
            dist(jq)=sqrt(dx**2+dy**2+(qdep(iq)-qdep(j))**2)
            indexjq(jq)=j
60       continue
         njq=jq

         if (njq.ge.5) then
            call INDEXX(njq,dist,indist)  !most time consuming part (?)
         else
            go to 200
         end if

         if (iq.eq.0) then
            print *,'***TEST SECTION FOR SSST_DIST***'
            print *,'Target quake = ',qlat(iq),qlon(iq),qdep(iq)
            do 70 i=1,20
               ii=indist(i)
               print *,ii,i,qlat(ii),qlon(ii),qdep(ii),dist(ii)
70          continue  
         end if

         iscr1=ipick1(iq)
         iscr2=ipick2(iq)
         do 120 ipick=iscr1,iscr2
            imed=0
            is_targ=isin(ipick)
            ip_targ=ipin(ipick)
            if(resid(ipick).ne.-999.) then
               imed=imed+1
               rr(imed)=resid(ipick)
            end if
            
            do 100 jq=1,njq
               jq2=indist(jq)
               iq2=indexjq(jq2)
               iscr3=ipick1(iq2)
               iscr4=ipick2(iq2)
               do 80 jpick=iscr3,iscr4
                  if (isin(jpick).ne.is_targ) go to 80
                  if (ipin(jpick).ne.ip_targ) go to 80
                  if (nqqual(jpick) .eq. -1) go to 100
                  if(isin(jpick).eq.is_targ.and.
     &               ipin(jpick).eq.ip_targ.and.
     &               resid(jpick).eq.-999.) go to 100
                  imed=imed+1
                  rr(imed)=resid(jpick)
                  if (imed.eq.nmed) go to 105
                  go to 100     !there should only be one matching pick
80             continue
100         continue
            if (imed.lt.5) then
               term(ipick)=0.
               go to 120
            end if
105         if (inorm .eq. 1) then         !L1 NORM
               call MEDIAN(rr,imed,rmed)   
            else if (inorm .eq. 2) then    !L2 NORM
               call MEAN(rr,imed,rmed)
            else if (inorm .eq. 3) then    !ROBUST MEAN 
               xgap=0.1
               call ROBOMEAN2(rr,imed,xgap,10,rmed,fit2)
            end if
            term(ipick)=rmed
            
            if (iq.eq.0) then
               print *,ipick,resid(ipick),term(ipick),imed
               if (ipick.eq.iscr1) then
                  do 110 i=1,imed
                     print *,i,rr(i)
110               continue
               end if
            end if

120      continue !end loop on picks for target quake

200   continue               !end loop on quakes
      
      return

      end
!
!
!-----------------------------------------------------------------------
!
! subroutine(15)
!
! STATERM computes station terms from residual list
!
!  Input:  inorm --- 1 for L1 norm ,2 for L2 norm or 3 for ROBUST MEAN
!          resid --- array with residuals
!          isin --- array of station number
!          ipin --- array of phase numbers
!          nr --- number of residuals
!          ns --- number of stations
!          np --- number of phases
!          nstpick --- array with number of picks per station+phase
!          nqqual --- array with earthquake location quality
!                     -1, bad location, should not be used in station terms calculation
!          
!  Output:  term --- array with station terms (matches resid array)
!           stterm --- !array with terms at station+phase 
!
!-----------------------------------------------------------------------
!
      subroutine STATERM(inorm,resid,isin,ipin,nr,ns,np,nstpick,
     &                   nqqual,term,stterm)
      implicit none


      integer nmax
      parameter (nmax=500000)
      integer nsta0
      parameter (nsta0=3000)
      integer i
      integer inorm
      integer ip
      integer*2 ip2
      integer is
      integer*2 is2
      integer nbin
      integer np
      integer nr
      integer ns
      integer ntest
      integer*2 ipin(nr)
      integer*2 isin(nr)
      integer nqqual(nr)
      integer nstpick(nsta0,2)


      real fit2
      real rmed
      real xgap
      real resid(nr)
      real rr(nmax)
      real stterm(nsta0,2)
      real term(nr)
!
!-----------------------------------------------------------------------
!
      print *,'STATERM nr,ns,np = ',nr,ns,np
      ntest=0
      do i=1,nr
         if (resid(i).ne.0..and.
     &       resid(i).ne.-999.) ntest=ntest+1
      enddo 
      print *,'Number of not bad residuals = ',ntest
      
      do 30 i=1,nr
         term(i)=0.
30    continue
      
      do 100 ip=1,np
         ip2=ip
         do 90 is=1,ns
            stterm(is,ip)=0.
            if (nstpick(is,ip).lt.5) go to 90   !skip if less than 5 picks
            is2=is
            nbin=0
            do 50 i=1,nr
               if (isin(i).ne.is2) go to 50
               if (ipin(i).ne.ip2) go to 50 
               if (nqqual(i) .eq. -1) go to 50
               if (resid(i).eq.-999.) go to 50
               nbin=nbin+1
               rr(nbin)=resid(i)
50          continue
            if (nbin.eq.0) go to 90
            
            if (inorm .eq. 1) then    !L1 NORM
               call MEDIAN(rr,nbin,rmed)   
            else if (inorm .eq. 2) then    !L2 NORM
               call MEAN(rr,nbin,rmed)
            else if (inorm .eq. 3) then    !ROBUST MEAN
               xgap=0.1
               call ROBOMEAN2(rr,nbin,xgap,10,rmed,fit2)
            end if
            
            do 60 i=1,nr
               if (isin(i).ne.is2) go to 60
               if (ipin(i).ne.ip2) go to 60 
               if (nqqual(i) .eq. -1) go to 60
               term(i)=rmed
60          continue
            stterm(is,ip)=rmed
90       continue
100   continue
      
      return
      end
!
!
!-----------------------------------------------------------------------
