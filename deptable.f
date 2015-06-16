!-----------------------------------------------------------------------
!
!    DEPTABLE creates tables of seismic travel times and ray takeoff
!    angles as a function of distance and source depth.  
!
!    These tables are designed to be read with the GET_TTS subroutine.
!
!-----------------------------------------------------------------------
!
      program deptable
      implicit none

      integer npts0   !npts0 = maximum number of lines in input velocity model
      parameter (npts0=1000)
      integer nz0     !nz0 = maximum number of depths in desired output tables
      parameter (nz0=100)
      integer nx0     !nx0 = maximum number of distances is desired output tables
      parameter (nx0=500)
      integer nray0   !nray0 = maximum number of rays during ray tracing
      parameter (nray0=40002)
      integer ncount0 !ncount0 = 2*nray0
      parameter (ncount0=80004)

      integer npts
      integer ndep
      integer nump
      integer np
      integer ncount
      integer ndel
      integer ideptype
      integer idep
      integer itype
      integer icount
      integer idel
      integer iw
      integer imth
      integer irtr
      integer i
      integer j
      integer i2
      integer ideprad


      real erad
      parameter (erad=6371.) 
      real pi
      parameter (pi=3.1415927)

      real ecircum
      real kmdeg
      real degrad
      real angle
      real p
      real pmin
      real pmax
      real pstep
      real plongcut
      real zmax
      real frac
      real xcore
      real tcore
      real h
      real dep
      real dep1
      real dep2
      real dep3
      real del
      real del1
      real del2
      real del3
      real deldel
      real scr1
      real scr2
      real scr3
      real scr4
      real xold
      real x
      real x1
      real x2
      real dx
      real t
      real t1
      real t2
      real dt
      real tbest
      real pbest
      real ubest
      real xsave(ncount0)
      real tsave(ncount0)
      real psave(ncount0)
      real usave(ncount0)
      real deptab(nz0)
      real ptab(nray0)
      real delttab(nray0)
      real z(npts0)
      real alpha(npts0)
      real beta(npts0)
      real z_s(npts0)
      real alpha_s(npts0)
      real beta_s(npts0)
      real slow(npts0,2)
      real deltab(nray0)
      real tttab(nray0)
      real angang(nx0,nz0)
      real tt(nx0,nz0)
      real rayray(nx0,nz0)
      real etaeta(nx0,nz0)
      real depxcor(nray0,nz0)
      real deptcor(nray0,nz0)
      real depucor(nray0,nz0)


      character*100 vmodel
      character*100 rfile
      character*100 ttfile
      character*100 angfile
      character*100 rayfile
      character*100 etafile
!
!-----------------------------------------------------------------------
!
      ecircum=2.*pi*erad
      kmdeg=ecircum/360.
      degrad=180./pi

      print *,'Enter input velocity model'
      read (*,'(a)') vmodel

      print *,'First column of input:  1=depth, 2=radius'
      read *,ideprad

!    read velocity model and transform to flat earth
!    value at center of earth is removed to
!    avoid singularity in transformation

      open (7,file=vmodel,status='old')
      do 10 i=1,npts0
         read (7,*,end=30) z_s(i),alpha_s(i),beta_s(i)
         if (ideprad.eq.2) z_s(i)=erad-z_s(i)
         if (z_s(i).eq.erad) go to 30
8        call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i))
         call FLATTEN(z_s(i),beta_s(i),z(i),beta(i))
10    continue
      print *,'***',npts0,' point maximum exceeded in model'
      stop
30    close (7)
      print *,'finished reading model'

38    z_s(i)=z_s(i-1)                  !set up dummy interface at bottom
      alpha_s(i)=alpha_s(i-1)
      beta_s(i)=beta_s(i-1)
      call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i))
      call FLATTEN(z_s(i),beta_s(i),z(i),beta(i))
      npts=i
      print *,'Depth points in model= ',npts

      do 40 i=1,npts
         slow(i,1)=1./alpha(i)
         if (beta(i).ne.0.) then
            slow(i,2)=1./beta(i)
         else
            slow(i,2)=1./alpha(i)              !fluid legs are always P!
         end if       
40    continue

      print *,'************************* Table of Model Interfaces *****
     &*****************'
      print *,' Depth  Top Velocities  Bot Velocities    -----Flat Earth
     & Slownesses-----'
      print *,'             vp1  vs1        vp2  vs2       p1      p2  
     &     s1      s2'

      do 45 i=2,npts
         if (i.eq.2.or.z(i).eq.z(i-1)) then
            scr1=1./alpha(i-1)
            scr2=1./alpha(i)
            scr3=999.
            if (beta(i-1).ne.0.) scr3=1./beta(i-1)
            scr4=999.
            if (beta(i).ne.0.) scr4=1./beta(i)
            print 42,z_s(i),i-1,alpha_s(i-1),beta_s(i-1),
     &              i,alpha_s(i),beta_s(i),
     &              scr1,scr2,scr3,scr4
42          format (f6.1,2(i5,f6.2,f5.2),2x,2f8.5,2x,2f8.5)
         end if
45    continue
!
!-----------------------------------------------------------------------
!
50    print *,'Enter maximum depth (9999 for none)'  
      read *,zmax                            !allow for reflected phases
      
      print *,'Source depths:  (1) Range, (2) Exact'
      read *,ideptype
      if (ideptype.eq.1) then
         print *,'Enter source dep1,dep2,dep3 (km)'
         read *,dep1,dep2,dep3
         dep2=dep2+dep3/20.
         idep=0.
         do 55 dep=dep1,dep2,dep3
            idep=idep+1
            deptab(idep)=dep
55       continue
         ndep=idep
      else if (ideptype.eq.2) then
         do 56 idep=1,80
            print *,'Enter source depth (km, -999. to stop)'
            read *,deptab(idep)
            if (deptab(idep).eq.-999.) go to 57
56       continue
         idep=51
57       ndep=idep-1
      end if
      
      print *,'(1) P-waves  or  (2) S-waves'
      read *,iw
      
      pmin=0.
      pmax=slow(1,iw)
      print *,'pmin, pmax = ', pmin, pmax
60    print *,'Enter number of rays to compute'
      read *, nump
      
      if (nump.gt.nray0) then
         print *,'nump must not exceed ',nray0
         go to 60
      end if

      pstep=(pmax-pmin)/float(nump)

      print *,'Enter min p at long range (.133 = no Pn, .238 = Sn)'
      read *,plongcut
!
! ------------------------------------------ ray tracing
!
      np=0
      do 200 p=pmin,pmax+pstep/2.,pstep
         np=np+1
         ptab(np)=p

         x=0.
         t=0.
         xcore=0.
         tcore=0.
         imth=3              !preferred v(z) interpolation method
         do 70 idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(np,idep)=0.
               deptcor(np,idep)=0.
               depucor(np,idep)=slow(1,iw)
            else
               depxcor(np,idep)=-999.
               deptcor(np,idep)=-999.
               depucor(np,idep)=-999.
            end if
70       continue

         do 100 i=1,npts-1
         if (z_s(i).ge.zmax) then                          !exceeds zmax
            deltab(np)=-999.
            tttab(np)=-999.
            go to 200
         end if

         h=z(i+1)-z(i)
         if (h.eq.0.) go to 100                       !skip if interface
         call LAYERTRACE(p,h,slow(i,iw),slow(i+1,iw),imth,dx,dt,irtr)
         x=x+dx
         t=t+dt

         if (irtr.eq.0.or.irtr.eq.2) go to 105           !ray has turned
         
         do 80 idep=1,ndep
            if (abs(z_s(i+1)-deptab(idep)).lt.0.1) then
               depxcor(np,idep)=x
               deptcor(np,idep)=t
               depucor(np,idep)=slow(i+1,iw)            
            end if
80       continue

!   
100      continue
105      x=2.*x
         t=2.*t

110      deltab(np)=x                   !stored in km
         tttab(np)=t                    !stored in seconds

200   continue             !---------------- end loop on ray parameter p
      print *,'Completed ray tracing loop'
      
!    special section to get (0,0) point
         np=np+1
         ptab(np)=slow(1,iw)
         deltab(np)=0.
         tttab(np)=0.
         do idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(np,idep)=0.
               deptcor(np,idep)=0.
               depucor(np,idep)=slow(1,iw)
            else
               depxcor(np,idep)=-999.
               deptcor(np,idep)=-999.
               depucor(np,idep)=-999.
            end if         
         enddo
    
      print *,'(1) del in km, t in sec  or  (2) del in deg., t in min'
      read *,itype
      
      print *,'Enter output file for ray table (or none)'
      read (*,'(a)') rfile
      if (rfile.ne.'none') then
         open (13,file=rfile) 
         do 207 i=1,np
            if (itype.eq.1) then
               write (13,*) i,ptab(i), deltab(i),tttab(i)
            else
               write (13,*) i,ptab(i),deltab(i)/kmdeg,tttab(i)/60.
            end if
207      continue
      end if

      print *,'Enter del1,del2,del3 (min, max, spacing)'
      read *,del1,del2,del3

      angle=-9.99
      do 250 idep=1,ndep
         icount=0
         xold=-999.
         if (deptab(idep).eq.0.) then
            i2=np
            go to 223
         end if
         do 220 i=1,np                         !upgoing rays from source
            x2=depxcor(i,idep)
            if (x2.eq.-999.) go to 221
            if (x2.le.xold) go to 221            !stop when heads inward
            t2=deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=-ptab(i)         !save as negative for upgoing from source
            usave(icount)=depucor(i,idep)
            xold=x2
220      continue
221      continue
         i2=i-1
223      do 225 i=i2,1,-1                    !downgoing rays from source
            if (depxcor(i,idep).eq.-999.) go to 225
            if (deltab(i).eq.-999.) go to 225
            x2=deltab(i)-depxcor(i,idep)
            t2=tttab(i)-deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=ptab(i)
            usave(icount)=depucor(i,idep)
            xold=x2
225      continue
226      ncount=icount
         idel=0
         do 240 deldel=del1,del2,del3
            del=deldel
            if (itype.eq.2) del=deldel*kmdeg
            idel=idel+1
            delttab(idel)=deldel
            tt(idel,idep)=99999.
            do 230 i=2,ncount
               x1=xsave(i-1)
               x2=xsave(i)
               if (x1.gt.del.or.x2.lt.del) go to 230
               if (psave(i).gt.0..and.psave(i).lt.plongcut) go to 230
               frac=(del-x1)/(x2-x1)
               tbest=tsave(i-1)+frac*(tsave(i)-tsave(i-1))
               if (psave(i-1).le.0..and.psave(i).le.0. .or.
     &                    psave(i-1).ge.0..and.psave(i).ge.0.) then
                  pbest=psave(i-1)+frac*(psave(i)-psave(i-1))
                  ubest=usave(i-1)+frac*(usave(i)-usave(i-1)) 
               else
                  if (frac.lt.0.5) then
                     pbest=psave(i-1)
                     ubest=usave(i-1)
                  else
                     pbest=psave(i)
                     ubest=usave(i)
                  end if
               end if
              
               if (tbest.lt.tt(idel,idep)) then
                  tt(idel,idep)=tbest
                  scr1=pbest/ubest
                  if (scr1.gt.1.) then
                     print *,'***Warning: p>u in angle calculation'
                     print *,'   Ray assumed horizontal'
                     print *,deptab(idep),del,tbest,pbest,ubest
                     scr1=1.
                  end if
                  angle=asin(scr1)*degrad
                  if (angle.lt.0.) then
                     angle=-angle
                  else
                     angle=180.-angle
                  end if
                  angang(idel,idep)=angle
                  rayray(idel,idep)=pbest
                  etaeta(idel,idep)=ubest*sqrt(1.-scr1**2)              
                  if (angang(idel,idep).lt.90.) then
                     etaeta(idel,idep)=-etaeta(idel,idep)
                  endif
               end if
230         continue
!            if (tt(idel,idep).eq.99999.) tt(idel,idep)=0.
            if (tt(idel,idep).eq.99999.) tt(idel,idep)=-9.99
            if (itype.eq.2) tt(idel,idep)=tt(idel,idep)/60.            

240      continue                                     !end loop on range
         ndel=idel
                  
250   continue                                        !end loop on depth

      if (delttab(1).eq.0.) then
         if (deptab(1).eq.0.) tt(1,1)=0.                 !set tt to zero at (0,0)
         do 255 idep=1,ndep
            angang(1,idep)=0.                 !straight up at zero range
            etaeta(1,idep)=-abs(etaeta(1,idep))
255      continue

         go to 280
         do i=2,ndel                          !now fix tt=0 at short distances                        
            if (tt(i,1).ne.0) then               !first non-zero time
               if (i.eq.2) go to 280             !no problem to fix
               do j=2,i-1
                  t1=tt(1,1)
                  t2=tt(i,1)
                  frac=(delttab(j)-delttab(1))/(delttab(i)-delttab(1))
                  tt(j,1)=t1+frac*(t2-t1)
                  rayray(j,1)=rayray(i,1)
                  etaeta(j,1)=etaeta(i,1)
                  angang(j,1)=angang(i,1)
               enddo
               go to 280
            end if
         enddo
280      continue

      end if

      print *,'Enter output file name for travel times'
      read (*,'(a)') ttfile
      open (11,file=ttfile)
      print *,'Enter output file name for source ray angles'
      read (*,'(a)') angfile
      open (12,file=angfile)
      print *,'Enter output file name for ray parameters'
      read (*,'(a)') rayfile
      open (13,file=rayfile)
      print *,'Enter output file name for vertical slowness at source'
      read (*,'(a)') etafile
      open (14,file=etafile)      
            
      write (11,407) vmodel,iw,pmin,pmax,np
      write (12,407) vmodel,iw,pmin,pmax,np
      write (13,407) vmodel,iw,pmin,pmax,np
      write (14,407) vmodel,iw,pmin,pmax,np      
407   format ('From deptable, file= ',a20,' iw =',i2,' pmin=',f8.5,
     &        ' pmax=',f8.5,' np=',i6)
      write (11,408) ndel,ndep
      write (12,408) ndel,ndep
      write (13,408) ndel,ndep
      write (14,408) ndel,ndep
      write (11,409) (deptab(j),j=1,ndep)
      write (12,409) (deptab(j),j=1,ndep)
      write (13,409) (deptab(j),j=1,ndep)
      write (14,409) (deptab(j),j=1,ndep)
      do 420 i=1,ndel
         if (itype.eq.1) then
            write (11,410) delttab(i),(tt(i,j),j=1,ndep)
         else
            write (11,413) delttab(i),(tt(i,j),j=1,ndep)
         end if
         write (12,411) delttab(i),(angang(i,j),j=1,ndep)
         write (13,415) delttab(i),(rayray(i,j),j=1,ndep)
         write (14,415) delttab(i),(etaeta(i,j),j=1,ndep)
408      format (2i5)
409      format (8x,100f8.1)
410      format (101f8.4)
411      format (f8.4,100f8.2)
413      format (f8.3,100f8.4)
415      format (f8.4,100f8.4)

420   continue
      close (11)
      close (12)
      close (13)
      close (14)

999   stop

      end
!
!-----------------------------------------------------------------------
! INTERP finds the y3 value between y1 and y2, using the
! position of x3 relative to x1 and x2.
      subroutine INTERP(x1,x2,x3,y1,y2,y3)
      implicit none

      real x1
      real x2
      real x3
      real y1
      real y2
      real y3
      real fract
  
      fract=(x3-x1)/(x2-x1)
      y3=y1+fract*(y2-y1)
      return
      end
!
!-----------------------------------------------------------------------
! FLATTEN calculates flat earth tranformation.
      subroutine FLATTEN(z_s,vel_s,z_f,vel_f)
      implicit none

      real erad
      real r
      real z_s
      real vel_s
      real z_f
      real vel_f
 
      erad=6371.
      r=erad-z_s
      z_f=-erad*alog(r/erad)
      vel_f=vel_s*(erad/r)
      return
      end

!-----------------------------------------------------------------------
! UNFLATTEN is inverse of FLATTEN.
      subroutine UNFLATTEN(z_f,vel_f,z_s,vel_s)
      implicit none
 
      real erad
      real r
      real z_f
      real vel_f
      real z_s
      real vel_s 

      erad=6371.
      r=erad*exp(-z_f/erad)
      z_s=erad-r
      vel_s=vel_f*(r/erad)
      return
      end
!
!-----------------------------------------------------------------------
!
!     THIS INTEGER*4 FUNCTION RETURNS INDEX OF LAST NONBLANK CHARACTER
!     IN <STRING>. 
!
      FUNCTION LENB (STRING)
      implicit none

      integer LENB
      integer N
      integer I

      CHARACTER*(*) STRING

      N=LEN(STRING)
      DO 10 I=N,1,-1
         IF (STRING(I:I).NE.' ') THEN
            LENB=I
            RETURN
         END IF
10    CONTINUE
11    LENB=0
      RETURN
      END
!
!-----------------------------------------------------------------------
!
! LAYERTRACE calculates the travel time and range offset
! for ray tracing through a single layer.
!
! Input:    p     =  horizontal slowness
!           h     =  layer thickness
!           utop  =  slowness at top of layer
!           ubot  =  slowness at bottom of layer
!           imth  =  interpolation method
!                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)     fastest to compute
!                         = 2,  v(z) = a - b*z               linear gradient
!                         = 3,  v(z) = a*exp(-b*z)           preferred when Earth Flattening is applied
!
! Returns:  dx    =  range offset
!           dt    =  travel time
!           irtr  =  return code
!                 = -1, zero thickness layer
!                 =  0,  ray turned above layer
!                 =  1,  ray passed through layer
!                 =  2,  ray turned within layer, 1 segment counted
!
! Note:  This version does calculation in double precision,
!        but all i/o is still single precision
!
      subroutine LAYERTRACE(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
      implicit none

      integer irtr
      integer imth

      real*4 p1
      real*4 h1
      real*4 utop1
      real*4 ubot1
      real*4 dx1
      real*4 dt1
      
      real*8 p
      real*8 h
      real*8 utop
      real*8 ubot 
      real*8 u
      real*8 y
      real*8 q
      real*8 qs
      real*8 qr
      real*8 b
      real*8 vtop
      real*8 vbot
      real*8 etau
      real*8 z
      real*8 dx
      real*8 dtau
      real*8 dt
      real*8 ex
!
!      implicit real*8 (a-h,o-z)
!
      p=dble(p1)
      h=dble(h1)
      utop=dble(utop1)
      ubot=dble(ubot1)
!
      if (h.eq.0.) then                  !check for zero thickness layer
         dx1=0.
         dt1=0.
         irtr=-1
         return         
      end if
!
      u=utop
      y=u-p
      if (y.le.0.) then                       !complex vertical slowness
         dx1=0.
         dt1=0.
         irtr=0
         return
      end if
!
      q=y*(u+p)
      qs=dsqrt(q)
!
! special function needed for integral at top of layer
      if (imth.eq.2) then
         y=u+qs
         if (p.ne.0.) y=y/p
         qr=dlog(y)
      else if (imth.eq.3) then
         qr=atan2(qs,p)
      end if      
!
      if (imth.eq.1) then
          b=-(utop**2-ubot**2)/(2.*h)
      else if (imth.eq.2) then
          vtop=1./utop
          vbot=1./ubot
          b=-(vtop-vbot)/h
      else
          b=-dlog(ubot/utop)/h
      end if  
!
      if (b.eq.0.) then                         !constant velocity layer
         b=1./h
         etau=qs
         ex=p/qs
         irtr=1
         go to 160
      end if
!
! integral at upper limit, 1/b factor omitted until end
      if (imth.eq.1) then
         etau=-q*qs/3.
         ex=-qs*p
      else if (imth.eq.2) then
         ex=qs/u                       !*** - in some versions (wrongly)
         etau=qr-ex
         if (p.ne.0.) ex=ex/p
      else
         etau=qs-p*qr
         ex=qr
      end if
!
! check lower limit to see if we have turning point
      u=ubot
      if (u.le.p) then                                !if turning point,
         irtr=2                                    !then no contribution
         go to 160                                    !from bottom point
      end if 
      irtr=1
      q=(u-p)*(u+p)
      qs=dsqrt(q)
!
      if (imth.eq.1) then
         etau=etau+q*qs/3.
         ex=ex+qs*p
      else if (imth.eq.2) then
         y=u+qs
         z=qs/u
         etau=etau+z
         if (p.ne.0.) then
            y=y/p
            z=z/p
         end if
         qr=dlog(y)
         etau=etau-qr
         ex=ex-z
      else
         qr=atan2(qs,p)
         etau=etau-qs+p*qr
         ex=ex-qr
      end if      
!
160   dx=ex/b
      dtau=etau/b
      dt=dtau+p*dx                                     !convert tau to t
!
      dx1=sngl(dx)
      dt1=sngl(dt)
      return
      end
!
!-----------------------------------------------------------------------
