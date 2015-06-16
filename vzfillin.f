!-----------------------------------------------------------------------
!
!    VZFILLIN is designed to read (z, Vp, Vs) models and resamples 
!    them at any desired finer depth interval.
!
!    This is done with linear interpolation.
!
!-----------------------------------------------------------------------
!
      program vzfillin
      implicit none

      integer npts0
      parameter (npts0=5000)

      integer i
      integer j
      integer k
      integer kk
      integer npts
      integer nk


      real sfact
      real dz
      real z
      real z1
      real z2
      real fact
      real d(npts0)
      real v(npts0,2)
      real v0(2)
      real buf(npts0,3)


      character*100 infile
      character*100 outfile
!
!-----------------------------------------------------------------------
!
      print *,'Enter input model file name'
      read (*,'(a)') infile
      open (11,file=infile,status='old')

      print *,'Enter output file name'
      read (*,'(a)') outfile
      open (12,file=outfile)
      
      print *,'Enter Vp/Vs ratio if Vs=0 in input (e.g., 1.732)'
      read *,sfact

      do 20 i=1,100
         read (11,*,end=21) d(i),(v(i,j),j=1,2)
         if (v(i,1).eq.0.) v(i,1)=v(i,2)*sfact
         if (v(i,2).eq.0.) v(i,2)=v(i,1)/sfact
20    continue
21    npts=i-1
      close (11)

      print *,'Enter output dz spacing (e.g. 1)'
      read *,dz

      kk=0
      do 50 z=0.,d(npts),dz
         do 40 i=1,npts-1
            z1=d(i)
            z2=d(i+1)
            if (z1.le.z.and.z2.ge.z) then
               if (z1.eq.z2) then
                  fact=0.
               else
                  fact=(z-z1)/(z2-z1)
               end if
               do 30 j=1,2
                  v0(j)=v(i,j)+fact*(v(i+1,j)-v(i,j))
30             continue
               kk=kk+1
               buf(kk,1)=z
               buf(kk,2)=v0(1)
               buf(kk,3)=v0(2)               
33             format (f10.3,2f8.4)
               do 38 k=i+1,npts
                  if (d(k)-z.lt.dz) then
                     kk=kk+1
                     buf(kk,1)=d(k)
                     buf(kk,2)=v(k,1)
                     buf(kk,3)=v(k,2)
                  end if
38             continue
               go to 50
            end if
40       continue
50    continue

      nk=kk
      write (12,33) (buf(1,j),j=1,3)
      do k=2,nk
         if (buf(k,1).ne.buf(k-1,1).or.
     &       buf(k,2).ne.buf(k-1,2).or.
     &       buf(k,3).ne.buf(k-1,3)) then
           write (12,33) (buf(k,j),j=1,3)
         end if
      enddo
      close (12)

      stop
      end
!
!
!-----------------------------------------------------------------------
