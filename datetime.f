c DATETIME is a set of FORTRAN subroutines to perform date/time
c conversions for seismology applications.  They use 4 digit years 
c and are fully Y2K compatible.  Note that time differencing
c operations are in real*8 (double precision) seconds.  One
c limitation is that leap seconds are not implemented.
c
c                                   Last modified 11/11/99
c                                   Peter Shearer
c                                   pshearer@ucsd.edu


c DT_GET_JDAY gets Julian day from (year, month, day)
c   Example: 1999, 4, 15 gives 105
c
c Inputs:  iyear  =  year (4 digits)
c          imon   =  month (1 to 12)
c          iday   =  day
c Returns: jday   =  Julian day (1 to 366)
c
      subroutine DT_GET_JDAY(iyear,imon,iday,jday)
      integer jsum(12),jsum2(12)
      data jsum/ 0,31,59,90,120,151,181,212,243,273,304,334/
      data jsum2/0,31,60,91,121,152,182,213,244,274,305,335/
      if (mod(iyear,400).eq.0) then       !leap year
         jday=jsum2(imon)+iday
      else if (mod(iyear,100).eq.0) then  !not a leap year
         jday=jsum(imon)+iday
      else if (mod(iyear,4).eq.0) then    !leap year
         jday=jsum2(imon)+iday
      else                                !not a leap year
         jday=jsum(imon)+iday
      end if
      return
      end


c DT_GET_DAY gets standard (month, day) from (year, Julian day)
c   Example: 1999, 105 gives 4, 15
c
c Inputs:  iyear  =  year (4 digits)
c          jday   =  Julian day (1 to 366)
c Returns: imon   =  month (1 to 12)
c          iday   =  day
c
      subroutine DT_GET_DAY(iyear,jday,imon,iday)
      integer jsum(12),jsum2(12)
      logical leapyear
      data jsum/ 0,31,59,90,120,151,181,212,243,273,304,334/
      data jsum2/0,31,60,91,121,152,182,213,244,274,305,335/
      if (mod(iyear,400).eq.0) then       !leap year
         leapyear=.true.
      else if (mod(iyear,100).eq.0) then  !not a leap year
         leapyear=.false.
      else if (mod(iyear,4).eq.0) then    !leap year
         leapyear=.true.
      else                                !not a leap year
         leapyear=.false.
      end if
      if (leapyear) then
         do 10 i=2,12
            if (jsum2(i).ge.jday) then
               iday=jday-jsum2(i-1)
               imon=i-1
               return
            end if
10       continue
         iday=jday-jsum2(12)
         imon=12
      else
         do 20 i=2,12
            if (jsum(i).ge.jday) then
               iday=jday-jsum(i-1)
               imon=i-1
               return
            end if
20       continue
         iday=jday-jsum(12)
         imon=12
      end if
      return
      end


c DT_GET_TDAY gets total number of days since 0 Jan 1600
c   Example: 1999,3,13 gives 145803
c
c Inputs:  iyear  =  year (4 digits)
c          imon   =  month (1 to 12)
c          iday   =  day (1 to 31)
c Returns: itday  =  days from 0 Jan 1600
c
      subroutine DT_GET_TDAY(iyear,imon,iday,itday)
      call DT_GET_JDAY(iyear,imon,iday,jday)
      itday=(iyear-1600)*365+jday
      itday=itday+int((iyear-1601)/4)    !add mult 4 LY
      itday=itday-int((iyear-1601)/100)  !subtract mult 100 lack of LY
      itday=itday+int((iyear-1601)/400)  !add mult 400 LY
      return
      end


c DT_GET_YMD gets year,month,day from total days since 0 Jan 1600
c   Example: 145803 gives 1999,3,13
c
c Inputs:  itday  =  days from 0 Jan 1600
c Returns: iyear  =  year (4 digits)
c          imon   =  month (1 to 12)
c          iday   =  day (1 to 31)
c
      subroutine DT_GET_YMD(itday,iyear,imon,iday)
      iy1=1600+int(float(itday)/365.24)-2
      iy2=iy1+3
      do 10 iyr=iy1,iy2
         call DT_GET_TDAY(iyr,1,1,itday0)
         if (itday0.gt.itday) go to 20
10    continue
      print *,'***DT_GET_YMD problem.  Year not found',itday,iy1,iy2
      stop
20    iyear=iyr-1
      call DT_GET_TDAY(iyear,1,1,itday0)
      jday=itday-itday0+1
      call DT_GET_DAY(iyear,jday,imon,iday)
      return
      end


c DT_TIMEDIF finds time difference between two date/times
c
c Inputs:  iyr1,imon1,idy1,ihr1,imn1,sc1  =  1st year,mon,day,hour,min,sec
c          iyr2,imon2,idy2,ihr2,imn2,sc2  =  2nd year,mon,day,hour,min,sec 
c Returns: timdif  =  2nd - 1st time (real*8 seconds)
c
c Note:  timdif is real*8, sc1,sc2 are real*4, other I/O are integers
c
      subroutine DT_TIMEDIF(iyr1,imon1,idy1,ihr1,imn1,sc1,
     &                      iyr2,imon2,idy2,ihr2,imn2,sc2,timdif)
      real*8 timdif
      call DT_GET_TDAY(iyr1,imon1,idy1,itday1)
      call DT_GET_TDAY(iyr2,imon2,idy2,itday2)
      timdif =  dble(3600.)*dble(ihr2-ihr1)
     &         +dble(  60.)*dble(imn2-imn1)+dble(sc2-sc1)
      timdif=timdif+dble(86400.)*dble(itday2-itday1)
      return
      end


c DT_ADDTIME adds time offset (seconds) to yr,mn,dy,hr,mn,sc time
c
c Inputs:  iyr1,imon1,idy1,ihr1,imn1,sc1  =  year,mon,day,hour,min,sec
c          timdif                         =  seconds to add (real*8)
c Returns: iyr2,imon2,idy2,ihr2,imn2,sc2  =  year,mon,day,hour,min,sec
c
c Note:  timdif is real*8, sc1,sc2 are real*4, other I/O are integers
c
      subroutine DT_ADDTIME(iyr1,imon1,idy1,ihr1,imn1,sc1,
     &                      iyr2,imon2,idy2,ihr2,imn2,sc2,timdif)
      real*8 timdif,timdif2,sc8
      iyr2=iyr1
      imon2=imon1
      idays=int(timdif/dble(86400))
      timdif2=timdif-dble(idays)*dble(86400)
      idy2=idy1+idays
      ihr2=ihr1
      imn2=imn1
      sc8=dble(sc1)+timdif2
      call DT_CLEANTIME(iyr2,imon2,idy2,ihr2,imn2,sc8)
      sc2=real(sc8)
      return
      end


c DT_CLEANTIME cleans up oddball yr,mn,dy,hr,mn,sc times
c Note that sc is real*8
      subroutine DT_CLEANTIME(yr,mon,dy,hr,mn,sc)
      integer yr,mon,dy,hr,mn,dmn,dhr,ddy
      real*8 sc
c
      dmn=int(sc/60.)
      sc=sc-dble(60.*float(dmn))
      if (sc.lt.0.) then
         sc=sc+dble(60.)
         dmn=dmn-1
      end if
      mn=mn+dmn
c      
      dhr=mn/60
      mn=mn-60*dhr
      if (mn.lt.0) then
         mn=mn+60
         dhr=dhr-1
      end if
      hr=hr+dhr
c
      ddy=hr/24
      hr=hr-24*ddy
      if (hr.lt.0) then
         hr=hr+24
         ddy=ddy-1
      end if
c
      call DT_GET_TDAY(yr,mon,dy,itday)
      itday=itday+ddy
      call DT_GET_YMD(itday,yr,mon,dy)
c
      return
      end

