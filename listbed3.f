!-----------------------------------------------------------------------
!
!     program to list phase data in BED3 format
!
      program listbed3
      implicit none


      integer mp
      parameter (mp=3000)               !maximum number of picks per event
      integer menu
      integer nqmax
      integer iq
      integer ipick
      integer iend
      integer nn
      integer i
      integer j
      integer nqtot
      integer npickall
    

      character infile*100
!
! variables for header
      integer*4 nevent
      character*512 header
!
! variables for quake info
      integer iyr
      integer imon
      integer idy
      integer ihr 
      integer imn 
      integer*2 imag
      integer*2 istrike
      integer*2 idip
      integer*2 irake
      integer*2 ntrace
      integer*2 nppick
      integer*2 nspick
      integer*2 npicktot
      integer*2 npickused
      integer*4 cuspid

      real qsc
      real*4 qlat
      real*4 qlon
      real*4 qdep
      real*4 qlaterr
      real*4 qlonerr
      real*4 qdeperr
      real*4 t0err
      real*4 rms
      real*8 t0

      character*1 qtype
      character*1 mtype
      character*1 qual
      character*1 errtype
      character*1 foctype
!
! variables for station pick info
      real*4 tt
      real*4 delta

      character*1 onset
      character*1 pol
      character*1 acc
      character*3 comp
      character*4 junk
      character*6 phase
      character*12 stname
!
! buffer variables
      integer*2 ibufq(40)
      integer*2 ibufpick(18)
      integer*2 ibufp(18,mp)

      common/com1/t0,qlat,qlon,qdep,qlaterr,qlonerr,qdeperr,t0err,
     &            rms,cuspid,imag,istrike,idip,irake,ntrace,
     &            nppick,nspick,npicktot,npickused,qtype,mtype,qual,
     &            errtype,foctype

      common/com2/tt,delta,stname,comp,phase,onset,pol,acc,junk

      equivalence (ibufq(1),t0)
      equivalence (ibufpick(1),tt)
!
!-----------------------------------------------------------------------
!
      print *,'Enter input file name'
      read (*,'(a)') infile
      open (11,file=infile,form='unformatted')
      call RDBLK(nevent,2,11,iend)
      call RDBLK(header,256,11,iend)


      print *,'Print:  (1) all events,  (2) last only'
      read *,menu

      print *,'Enter maximum number of quakes to list'
      read *,nqmax

      print *,'Now reading quakes'

      ipick=0
      do 50 iq=1,nqmax

10       call RDBLK(ibufq,40,11,iend)
         if (iend.eq.1) go to 150

         if (menu.eq.1) then
            call DT_ADDTIME(1600,1,0,0,0,0.,
     &                   iyr,imon,idy,ihr,imn,qsc,t0)
            print *,t0,qlat,qlon,qdep,imag,cuspid,npicktot
            print *,iyr,imon,idy,ihr,imn,qsc
         end if

         nn=18*npicktot
         call RDBLK(ibufp,nn,11,iend)

         do 20 i=1,npicktot
            ipick=ipick+1
            do 15 j=1,18
               ibufpick(j)=ibufp(j,i)
15          continue
            if (menu.eq.1) then
               print 197, stname,delta,comp,phase,tt,
     &                    onset,pol,acc,junk
            end if
20       continue

50    continue

150   close (11)

      nqtot=iq-1
      npickall=ipick
      print *,'nqtot,npickall = ',nqtot,npickall

      call DT_ADDTIME(1600,1,0,0,0,0.,
     &                iyr,imon,idy,ihr,imn,qsc,t0)

      print *,'Last quake follows: '
      print *,t0,qlat,qlon,qdep,imag,cuspid,npicktot
      print *,iyr,imon,idy,ihr,imn,qsc
      do 200 i=1,npicktot
         do 195 j=1,18
            ibufpick(j)=ibufp(j,i)
195      continue
         print 197, stname,delta,comp,phase,tt,onset,pol,acc,junk
197      format (a12,x,f9.3,x,a3,x,a6,x,f8.3,x,a1,x,a1,x,a1,x,a4)
200   continue         

      stop
      end
!
!
!-----------------------------------------------------------------------
!
      subroutine RDBLK(a,n,nun,iend)
      integer*2 a(n)
      iend=0
      read (nun,end=99) a
      return
99    iend=1
      return
      end
!
!
!-----------------------------------------------------------------------
