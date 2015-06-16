!-----------------------------------------------------------------------
!
!   getstlist is the program to convert station list from website into 
!   the station list format needed by our location program
!
!   Only keep the first one if a station listed more than once.
!
!   SCSN station list format changed corresponding to the format change
!   in the SCEDC station list format from web
!                                                 04/18/2007 Guoqing Lin
!-----------------------------------------------------------------------
!
      program getstlist       
      implicit none
    
      integer maxnst
      parameter (maxnst=30000)

      integer iformat
      integer is
      integer nsta

      real slat(maxnst)
      real slon(maxnst)
      real selev(maxnst)

      character*12  sname(maxnst)
      character*100 infile
      character*100 outfile
!
!-----------------------------------------------------------------------
!
      print *,'Enter the input station list format'
      print *,'(1) SCSN'
      print *,'(2) NCSN'
      read *,iformat

      print *,'Enter the input station list file name'
      read (*,'(a)') infile

      print *,'Enter the output station list file name'
      read (*,'(a)') outfile

      if (iformat .eq. 1) then
      call GETSTAT_SCSN(infile,nsta,sname,slat,slon,selev)

      else if (iformat .eq. 2) then
      call GETSTAT_NCSN(infile,nsta,sname,slat,slon,selev)
      end if
      print *,'# stations = ',nsta

      open(12, file=outfile)
      do 30 is=1,nsta 
         write(12,31) sname(is),slat(is),slon(is),selev(is)
31       format (a12,f10.5,f12.5,f10.1)
30    continue
      close(12)
      end  
!
!-----------------------------------------------------------------------
!
      subroutine GETSTAT_SCSN(infile,nsta,sname,slat,slon,selev)
      implicit none

      integer maxnst
      parameter (maxnst=30000)
      integer nsta
      integer i
      integer is
      integer ns

      real flat
      real flon
      real felev  
      real slat(maxnst)
      real slon(maxnst)
      real selev(maxnst)

      character*2  netname
      character*3  compn
      character*4  stname
      character*12 sname(maxnst)
      character*100 infile
!
!-----------------------------------------------------------------------
!
      open (19,file=infile,status='old')
      ns=0
      do 10 i=1,maxnst
         read (19,11,end=12) stname,compn,flat,flon,felev,netname
!11       format (a4,1x,a3,32x,f10.5,f11.5,f6.0,23x,a2)   !old format
11       format (a4,2x,a3,32x,f10.5,f11.5,f6.0,23x,a2)

         if (ns .gt. 0) then
         do is=1,ns
            if (stname .eq. sname(is)(4:7)) go to 10
         end do
         end if
            
         ns=ns+1
         sname(ns)=netname//' '//stname//' '//compn
         slat(ns)=flat
         slon(ns)=flon
         selev(ns)=felev

10    continue 

12    nsta=ns
      close (19)

      return
      end
!
!-----------------------------------------------------------------------
!
      subroutine GETSTAT_NCSN(infile,nsta,sname,slat,slon,selev)
      implicit none

      integer maxnst
      parameter (maxnst=30000)
      integer nsta
      integer i
      integer is 
      integer ns 
      integer slat1
      integer slon1

      real slat2
      real slon2
      real felev
      real slat(maxnst)
      real slon(maxnst)
      real selev(maxnst)

      character*2  netname
      character*3  compn
      character*4  stname
      character*12 snam
      character*12 sname(maxnst)
      character*100 infile
!
!-----------------------------------------------------------------------
!
      open (19,file=infile,status='old')
      ns=0
      do 10 i=1,maxnst
         read (19,11,end=12) snam,slat1,slat2,slon1,slon2,felev
11       format (a12,18x,i2,f8.4,i4,f8.4,f6.0)

         stname=snam(1:4)
         netname=snam(6:7)
         compn=snam(9:11)
         snam=netname//' '//stname//' '//compn

         if (ns .gt. 0) then
         do is=1,ns
            if (snam(4:7) .eq. sname(is)(4:7)) go to 10
         end do
         end if

         ns=ns+1
         sname(ns)=snam
         slat(ns) = real(slat1) + slat2/60.
         slon(ns) = -real(slon1) - slon2/60.
         selev(ns)=felev
10    continue 

12    nsta=ns
      close (19)

      return
      end
!
!-----------------------------------------------------------------------
