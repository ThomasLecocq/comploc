!-----------------------------------------------------------------------
!
!  Program to convert phase output format to bed3 format
!
! 03/01/05  set tt for bad picks zero, not skip this pick as before
!           skip them in location progam. 
!
! 02/20/06 add the new SCSN STP phase data format           Guoqing Lin
!
!-----------------------------------------------------------------------
!
      program phase2bed3
      implicit none

      integer i
      integer j
      integer k
      integer npick
      integer npickall
      integer nn
      integer nq
      integer npha
      integer iquake
      integer iline
      integer iflag

      real delkm

      character*1 chr1
      character*1 chr2
      character*100 phasefile
      character*100 bedfile
      character*150 linebuf
!
!
! variables for BED3 header
      integer*4 nevent
      character*512 header
!
!
! variables for reading SCSN quake info
!
      integer qyr
      integer qmon
      integer qdy
      integer qhr
      integer qmn
      integer qid
      integer nph
      integer neq
      integer iyr
      integer imon
      integer idy
      integer ihr
      integer imn
      integer isec
      integer idelkm
      integer i1
      integer i2
      integer i3
      integer i4
      integer i5
      integer i6
      integer ierr1
      integer ierr2
      integer idcusp

      real qsc
      real qmag
      real qlat
      real qlon
      real qdep
      real qt0
      real elat
      real elon
      real sec
      real qrms
      real qlaterr
      real qlonerr
      real qdeperr
      real qt0err

      character*1 evtype
      character*1 magtype
!
! variables for reading SCSN pick info
!
      real pick
!      real stlat
!      real stlon
!      real stelev
      real phweight

      character*2 stnet
      character*4 sname
      character*3 stcomp
      character*9 phinfo_bad
      character*11 phinfo
      character*12 stname
!
! variables for BED3 quake info
!
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

      real*4 qlat2
      real*4 qlon2
      real*4 qdep2
      real*4 qlaterr2
      real*4 qlonerr2
      real*4 qdeperr2
      real*4 rms
      real*8 t0        !total seconds since 1600

      character*1 qtype
      character*1 mtype
      character*1 qual
      character*1 errtype
      character*1 foctype
!
! variables for BED3 station pick info
!
      real*4 tt
      real*4 delta

      character*1 onset
      character*1 pol
      character*1 acc
      character*3 comp
      character*4 junk
      character*6 phase
      character*12 stname2
!
! BED3 buffer variables
!
      integer mp
      parameter (mp=3000)               !maximum number of picks per event

      integer*2 ibufq(40)
      integer*2 ibufpick(18)
      integer*2 ibufp(18,mp)

      real t0err
      common/com1/t0,qlat2,qlon2,qdep2,qlaterr2,qlonerr2,qdeperr2,t0err,
     &            rms,cuspid,imag,istrike,idip,irake,ntrace,
     &            nppick,nspick,npicktot,npickused,qtype,mtype,qual,
     &            errtype,foctype

      common/com2/tt,delta,stname2,comp,phase,onset,pol,acc,junk

      equivalence (ibufq(1),t0)
      equivalence (ibufpick(1),tt)

      phase='      '

      i=0
      j=0     !counter for total number of picks
      nq=0
      npick=0
      iflag=1
      print *,'Enter output file name for bed3 file'
      read (*,'(a)') bedfile
      open (12,file=bedfile,form='unformatted')

      nevent=0      !unknown at this point
      header='Converted from SCSN ascii format using stpphase_to_bed3.f'
      call WRBLK(nevent,2,12)
      call WRBLK(header,256,12)

110   print *,'Enter input file name with phase data (or none)'
      read (*,'(a)') phasefile
      if (phasefile.eq.'none') go to 900

      print *,'Enter phase data format '
      print *,'(1) SCSN STP (OLD)'
      print *,'(2) SCSN STP (NEW)'
      print *,'(3) HYPOINVERSE' 
      read *,npha

      open (11,file=phasefile,status='old')

      do 100 iline=1,99999999 
      read (11,'(a150)',end=299) linebuf
      
      if (npha.ne.3.and.linebuf(1:1).eq.'#') go to 299

      if ((npha.eq.1.and.linebuf(77:79).ne.'   ').or.             !OLD STP
     &    (npha.eq.2.and.linebuf(1:1).eq.' ') .or.                !NEW STP
     &    (npha.eq.3.and.(linebuf(1:2).eq.'19' .or.               !HYPOINVERSE
     &     linebuf(1:2).eq.'20'))) then                !QUAKE LINE
          
         if (i.ge.1.and.iflag.eq.1) then
         call DT_TIMEDIF(1600,1,0,0,0,0.,qyr,qmon,qdy,qhr,qmn,qsc,t0)
         t0=t0+qt0
         qlat2=qlat
         qlon2=qlon
         qdep2=qdep
         qlaterr2=qlaterr
         qlonerr2=qlonerr
         qdeperr2=qdeperr
         t0err=qt0err
         rms=qrms
         cuspid=qid
         imag=nint(qmag*1000)
         npicktot=npick      !safer to use actual count
         npickused=neq

         call WRBLK(ibufq,40,12)
         nn=18*npicktot
         call WRBLK(ibufp,nn,12)
         
         end if

         i=i+1
         iflag=1
         if (mod(i,1000) .eq. 0 .or. i.eq. 1) then
         print *,'iq=',i
         end if

         iquake=1
         if (npha.eq.1.or.npha.eq.2) then     !SCSN STP
         read (linebuf,131) qid,evtype,qyr,qmon,qdy,qhr,qmn,qsc,
     &                      qlat,qlon,qdep,qmag
!         print 131, qid,evtype,qyr,qmon,qdy,qhr,qmn,qsc,
!     &              qlat,qlon,qdep,qmag
131      format (i10,x,a1,2x,i4,x,i2,x,i2,x,i2,x,i2,x,f6.3,
     &           x,f9.4,x,f11.4,x,f6.2,x,f5.2)

         qt0=0.
         qrms=0.
         qlaterr=0.
         qdeperr=0.
         qt0err=0.
         nph=0
         neq=0
         magtype=' '

         else if (npha.eq.3) then !HYPOINVERSE
         read (linebuf,132) qyr,qmon,qdy,qhr,qmn,isec,i1,chr1,i2,
     &            i3,chr2,i4,i5,idcusp,i6
!     &            i3,chr2,i4,i5,i6,ierr1,ierr2,idcusp
!         print 132, qyr,qmon,qdy,qhr,qmn,isec,i1,chr1,i2,
!     &              i3,chr2,i4,i5,i6,ierr1,ierr2,idcusp
132   format (i4,4i2,i4,i2,a1,i4,i3,a1,i4,i5,100x,i10,1x,i3)
!132      format (i4,4i2,i4,i2,a1,i4,i3,a1,i4,i5,i2,47x,2i4,43x,i10)

         qsc=float(isec)/100.
         elat=float(i1)+float(i2)/6000.
         if (chr1.eq.'s'.or.chr1.eq.'S') elat=-elat
         elon=-(float(i3)+float(i4)/6000.)
         if (chr1.eq.'w'.or.chr1.eq.'W') elon=-elon
         qlat=elat
         qlon=elon
         qdep=float(i5)/100.
         qmag=float(i6)/100.

         qid=idcusp
         qt0=0.
         qrms=0.
         qlaterr=float(ierr1)/100.
         qdeperr=float(ierr2)/100.
         qt0err=0.
         nph=0
         neq=0
         evtype=' '
         magtype=' '

         end if

         npick=0

      else if ((npha.ne.3).or.
     &        (npha.eq.3.and.linebuf(1:4).ne.'    '))  then   !PICK LINE


      if(npha.ne.3) then  !SCSN STP

         if (npha.eq.1) then !OLD
            read (linebuf,171) stnet,sname,stcomp,
     &       phinfo_bad,phweight,delkm,pick
!            print 171, stnet,sname,stcomp,
!     &       phinfo_bad,phweight,delkm,pick
171      format (1x,a2,1x,a4,1x,a3,30x,a8,f3.1,f9.2,f9.3)   !tab is first char

         else if (npha.eq.2.and.linebuf(44:44) .eq. ' ') then  !NEW
            read (linebuf,172) stnet,sname,stcomp,
     &       phinfo_bad,phweight,delkm,pick
!            print 172, stnet,sname,stcomp,
!     &       phinfo_bad,phweight,delkm,pick
172      format (1x,a2,3x,a4,1x,a3,30x,a8,f3.1,f9.2,f9.3)   !tab is first char
         endif

         j=j+1
         npick=npick+1
         
         tt=pick-qt0
         delta=delkm
         stname2='            '
         stname2(1:2)=stnet
         if (sname(1:1).ne.' ') then
            stname2(4:7)=sname(1:4)
         else
            stname2(4:6)=sname(2:4)
         end if
         stname2(9:11)=stcomp

         comp=stcomp
         phase(1:1)=phinfo_bad(1:1)     !e.g., P or S
         onset=phinfo_bad(6:6)          !e.g., I or E
         acc=' '                    
         pol=phinfo_bad(3:3)            !e.g., U or D

         do k=1,18
            ibufp(k,npick)=ibufpick(k)
         end do
       
 
      else if (npha.eq.3) then      !HYPOINVERSE


         if (iquake.eq.0 .or.
     &       linebuf(6:7).eq.'  ' .or.
     &       linebuf(31:31).eq.'.' .or. 
     &       linebuf(32:32).eq.'.' .or. 
     &       linebuf(33:33).eq.'.' .or. 
     &       linebuf(34:34).eq.'.') go to 100

         j=j+1
         npick=npick+1
         
         read (linebuf,173) stname,phinfo(1:4),iyr,imon,idy,ihr,imn,
     &            isec,idelkm
!         print 173, stname,phinfo(1:4),iyr,imon,idy,ihr,imn,
!     &            isec,idelkm
173      format (a12,1x,a4,i4,4i2,i5,40x,i4)

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
         pick=real(t0)

         tt=pick-qt0
         if (iquake.eq.0 .or.
     &       linebuf(6:7).eq.'  ' .or.
     &       linebuf(31:31).eq.'.' .or. 
     &       linebuf(32:32).eq.'.' .or. 
     &       linebuf(33:33).eq.'.' .or. 
     &       linebuf(34:34).eq.'.') tt=-99.

         delta=delkm
         stname2='            '
         stname2(1:2)=stname(6:7)
         stname2(4:7)=stname(1:4)
         stname2(9:12)=stname(9:12)
         comp=stname(10:12)
         phase(1:1)=phinfo(2:2)     !e.g., P or S
         onset=phinfo(1:1)          !e.g., I or E
         acc=phinfo(4:4)            !quality
         pol=phinfo(3:3)            !e.g., U or D

         do k=1,18
            ibufp(k,npick)=ibufpick(k)
         end do

      end if
      
      end if

100   continue
299   call DT_TIMEDIF(1600,1,0,0,0,0.,qyr,qmon,qdy,qhr,qmn,qsc,t0)
      t0=t0+qt0
      qlat2=qlat
      qlon2=qlon
      qdep2=qdep
      qlaterr2=qlaterr
      qlonerr2=qlonerr
      qdeperr2=qdeperr
      t0err=qt0err
      rms=qrms
      cuspid=qid
      imag=nint(qmag*1000)
      npicktot=npick                !safer to use actual count
      npickused=neq

      call WRBLK(ibufq,40,12)
      nn=18*npicktot
      call WRBLK(ibufp,nn,12)
      close(11)
      nq=i
      npickall=j
      iflag=0
      go to 110
900   print *,'nq,npickall = ',nq,npickall
      close(12)

      stop
      end
!
!
!-----------------------------------------------------------------------
!
      subroutine WRBLK(a,n,nun)
      implicit none
    
      integer n
      integer nun
      integer*2 a(n)

      write (nun) a
      return
      end
!
!-----------------------------------------------------------------------
