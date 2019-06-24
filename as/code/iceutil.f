! source file: /net/mare/home1/eby/as/ism/iceutil.F
      subroutine itor (iarr, arr, n)
      dimension iarr(n), arr(n)
      do i=1,n
        arr(i) = iarr(i)
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine rtoi (arr, iarr, n)
      dimension arr(n), iarr(n)
      do i=1,n
        iarr(i) = nint(arr(i))
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
      FUNCTION LENCHR(CH)
c     Returns length of character variable ch, less trailing blanks
      CHARACTER*(*) CH
      DO 10 I = LEN(CH), 1, -1
        LENCHR = I
        IF (CH(I:I).NE.' ' .AND. CH(I:I).NE.CHAR(0)) RETURN
   10 CONTINUE
      RETURN
      END
c
      SUBROUTINE ZERO (ARR, NAR)
c     Zeros nar words starting at arr(1)
      DIMENSION ARR(NAR)
      DO 10 J=1,NAR
        ARR(J) = 0.
   10 CONTINUE
      RETURN
      END
c
      SUBROUTINE IZERO (IARR, NAR)
c     Zeros nar words starting at iarr(1)
      DIMENSION IARR(NAR)
      DO 10 J=1,NAR
        IARR(J) = 0
   10 CONTINUE
      RETURN
      END
c
c-----------------------------------------------------------------------
c
      subroutine reseti (ia, n, ival)

c        Sets integer array ia(n) to ival

      dimension ia(*)

      do 10 i=1,n
        ia(i) = ival
   10 continue

      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine resetr (a, n, val)

c        Sets real array a(n) to val

      dimension a(*)

      do 10 i=1,n
        a(i) = val
   10 continue

      return
      end

c-----------------------------------------------------------------------

      subroutine scopy_i (n, iarr, inca, ibrr, incb)
c     Like Cray scopy, except for integers
      dimension iarr(*), ibrr(*)
      ib = 1
      do ia=1,n,inca
        ibrr(ib) = iarr(ia)
        ib = ib + incb
      enddo
      return
      end
c
c-----------------------------------------------------------------------
c
c#if defined SUN
cc     to enable nint() with large values (eg, showice) with -i4:
cc     also declare in comdriveice.h:
c      integer*8 function nint(x)
c      nint = x  + sign(0.5, x)
c      return
c      end
c#endif
c
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
c
      function cvmgt (x,y,l)
c     Duplicates Cray-Vector-Merge-GT function.
c     Only safe when first two arguments are reals (especially
c     if integers are 32 bits and reals are 64 bits).
      logical l
      if (l) then
        cvmgt = x
      else
        cvmgt = y
      endif
      return
      end

      subroutine scopy (n, a, inca, b, incb)
c     Duplicates Cray scopy
      dimension a(*), b(*)
      ib = 1
      do ia=1,n,inca
        b(ib) = a(ia)
        ib = ib + incb
      enddo
      return
      end

      integer function ishell (a)
c     Duplicates Cray ishell
      character*(*) a
      integer system
      ishell = system (a)
      return
      end
c
c-----------------------------------------------------------------------
c

c-----------------------------------------------------------------------

      subroutine insolhalf (iyear, alorb, fluxh, ecc,obliq,prec,iwhich)

c        Computes summer half-year insolation fluxh (W/m2)
c        at latitude alorb (radians), for calendar year iyear.
c        iyear will be interpreted by zenorb* as "BP", ie, relative
c        to 1950 AD.
c
c        If iwhich = 0, ecc,obliq,prec are passed.
c        If iwhich = 1, call zenorb to get ecc,obliq,prec.

c     parameter (nflux=365)
      parameter (nflux=180)
c     parameter (nflux=24)   ! old fur.f
      dimension flux(nflux)

      parameter (pi=3.14159265358979)

      parameter (solcon = 1365.)          ! W/m2

      dimension
     *  cosq(1), fraq(1), cosq24(1), fraq24(1), alatq(1)

      alatq(1) = alorb

c          Get orbital elements for zencal

      if (iwhich.eq.1) then
c       call zenorb78 (iyear, ecc, obliq, prec, vern)
        call zenorb91 (iyear, ecc, obliq, prec, vern)
      else
c       just calculate vern as in zenorb:
        zs = sin (0.5*prec) / sqrt ((1.+ecc)/(1.-ecc))
        zc = cos (0.5*prec)
        ze = 2. * atan2 (zs,zc)
        vern = ze - ecc * sin(ze)
      endif

c        Loop through year in equal time increments, storing
c        24-hour mean fluxes at latitude alorb in flux().
c        Then bubble sort flux and take the mean of the largest
c        half of its members, to get the half-year insolation
c        (fluxh) at alorb.

      isecdy = 86400/2
      dt  = 86400.

      do m=1,nflux
        isecyr = nint (86400.*365.*(m-.5)/nflux)
        call zencal (iyear, isecyr, isecdy, dt,
     *               ecc, obliq, prec, vern, dist, eccf,
     *               cosq, fraq, cosq24, fraq24, alatq, 1)
        flux(m) = cosq(1)*fraq(1)*eccf*solcon
      enddo

      call bubblesort (flux, nflux)

      fluxh = 0.
      nfluxh = nflux/2 + 1
      do m = nfluxh, nflux
        fluxh = fluxh + flux(m)
      enddo
      fluxh = fluxh / (nflux-nfluxh+1)

      return
      end

c----------------------------------------------------------------------

      subroutine bubblesort (a, n)

c        Bubble sorts a(n) into ascending order

      dimension a(n)

      do ibub = n-1,1,-1
        do j=ibub,n-1
          if (a(j).gt.a(j+1)) then
            z = a(j+1)
            a(j+1) = a(j)
            a(j) = z
          else
            go to 10
          endif
        enddo
   10   continue
      enddo

      return
      end

c----------------------------------------------------------------------

      subroutine bubblemean (iarr, a, w, zmean, ilon, ilat, nlon,nlat,n)

c        (1) Bubble sorts a(n) into ascending order,sorting weights
c            w(n) to match.
c        (2) Subtracts weighted mean from a, and normalize w.
c        (3) Packs a and normalized w values into iarr(ilon,ilat):
c            iarr = 10000*nint(a+10000) + min(9999,nint(10000*w))
c            (So a+10000 must be >= 0 and w must be between 0 and 1).

      dimension iarr(nlon,nlat,n), a(n), w(n)

c        Bubble sort a(n) into ascending order, with w(n) following.

      do ibub = n-1,1,-1
        do j=ibub,n-1
          if (a(j).gt.a(j+1)) then
            z = a(j+1)
            a(j+1) = a(j)
            a(j) = z

            z = w(j+1)
            w(j+1) = w(j)
            w(j) = z
          else
            go to 10
          endif
        enddo
   10   continue
      enddo

c        Subtract out the mean of a (weighted by w), and normalize w

      zmean = 0.
      zarea = 0.
      do i=1,n
        zmean = zmean + a(i)*w(i)
        zarea = zarea +      w(i)
      enddo
      zmean = zmean/zarea
      do i=1,n
        a(i) = a(i) - zmean
        w(i) = w(i) / zarea
      enddo

c        Pack values and weights into iarr

      do i=1,n
        iarr(ilon,ilat,i) =
     *    10000*nint(a(i)+10000.) + min (9999, nint(10000.*w(i)))
      enddo

      return
      end

c----------------------------------------------------------------------

      SUBROUTINE ZENORB78 (IYEAR, ECC, OBLIQ, PREC, VERN)
c
c        Calculates orbital parameters from trigonometric series in
c        Berger (1978) JAS,35,2362-2367.
c
c     IYEAR = calendar year (BP) (supplied)
c     ECC   = eccentricity (returned)
c     OBLIQ = obliquity in radians (returned)
c     PREC  = precession, ie, prograde angle from perihelion to n.h.
c             vernal equinox, 0 to 2*pi radians (returned)
c     VERN  = 2*pi *  time from perihelion to n.h. vernal equinox
c             / perihelion-to-perihelion year (returned)
c
      PARAMETER (PI = 3.14159265358979)
c
      DIMENSION OA(20), OF(20), OP(20)
      DIMENSION EA(19), EF(19), EP(19)
      DIMENSION PA(10), PF(10), PP(10)
      SAVE OA,OF,OP,EA,EF,EP,PA,PF,PP
c
      DATA OA/
     *-2462.22E0, -857.32E0, -629.32E0, -414.28E0, -311.76E0,
     *  308.94E0, -162.55E0, -116.11E0,  101.12E0,  -67.69E0,
     *   24.91E0,   22.58E0,  -21.16E0,  -15.65E0,   15.39E0,
     *   14.67E0,  -11.73E0,   10.27E0,    6.49E0,    5.85E0/
c
      DATA OF/
     *31.609970E0, 32.620499E0, 24.172195E0, 31.983780E0, 44.828339E0,
     *30.973251E0, 43.668243E0, 32.246689E0, 30.599442E0, 42.681320E0,
     *43.836456E0, 47.439438E0, 63.219955E0, 64.230484E0,  1.010530E0,
     * 7.437771E0, 55.782181E0,  0.373813E0, 13.218362E0, 62.583237E0/
c
      DATA OP/
     *251.90E0, 280.83E0, 128.30E0, 292.72E0,  15.37E0,
     *263.79E0, 308.42E0, 240.00E0, 222.97E0, 268.78E0,
     *316.79E0, 319.60E0, 143.80E0, 172.73E0,  28.93E0,
     *123.59E0,  20.20E0,  40.82E0, 123.47E0, 155.69E0/
c
c
      DATA EA/
     * .01860798E0, .01627522E0,-.01300660E0, .00988829E0,-.00336700E0,
     * .00333077E0,-.00235400E0, .00140015E0, .00100700E0, .00085700E0,
     * .00064990E0, .00059900E0, .00037800E0,-.00033700E0, .00027600E0,
     * .00018200E0,-.00017400E0,-.00012400E0, .00001250E0/
c
      DATA EF/
     *  4.207205E0,  7.346091E0, 17.857263E0, 17.220546E0, 16.846733E0,
     *  5.199079E0, 18.231076E0, 26.216758E0,  6.359169E0, 16.210016E0,
     *  3.065181E0, 16.583829E0, 18.493980E0,  6.190953E0, 18.867793E0,
     * 17.425567E0,  6.186001E0, 18.417441E0,  0.667863E0/
c
      DATA EP/
     *  28.62E0, 193.78E0, 308.30E0, 320.19E0, 279.37E0,
     *  87.19E0, 349.12E0, 128.44E0, 154.14E0, 291.26E0,
     * 114.86E0, 332.09E0, 296.41E0, 145.76E0, 337.23E0,
     * 152.09E0, 126.83E0, 210.66E0,  72.10E0/
c
c
      DATA PA/
     * 7391.02E0, 2555.15E0, 2022.76E0,-1973.65E0, 1240.23E0,
     *  953.87E0, -931.75E0,  872.38E0,  606.35E0, -496.03E0/
c
      DATA PF/
     * 31.609970E0, 32.620499E0, 24.172195E0,  0.636717E0, 31.983780E0,
     *  3.138886E0, 30.973251E0, 44.828339E0,  0.991874E0,  0.373813E0/
c
      DATA PP/
     * 251.90E0, 280.83E0, 128.30E0, 348.10E0, 292.72E0,
     * 165.16E0, 263.79E0,  15.37E0,  58.57E0,  40.82E0/
c
c
      TPI = 2.*PI
      RPD = PI/180.
      RPS = RPD/3600.
c
c        Time in Berger series is years from 1950.0 AD.
c        In gcm: If IYEAR is positive, it is in "AD". If zero,
c               it is in "BP", i.e., relative to 1950 AD.
c        Here, always BP.

c     IF (IYEAR.GE.0) THEN
c       TIME = IYEAR - 1950
c     ELSE
        TIME = IYEAR
c     ENDIF
c
c        Series for obliquity
c
      X = 0.
      DO 20 I=1,20
        X = X + RPS*OA(I) * COS (RPS*OF(I)*TIME + RPD*OP(I))
20    CONTINUE
      OBLIQ = RPD*23.320556 + X
c
c        Series for ECC*SIN(PIE), ECC*COS(PIE), where PIE is perihelion
c        position relative to fixed stars
c
      ESP = 0.
      ECP = 0.
      DO 30 I=1,19
        X = RPS*EF(I)*TIME + RPD*EP(I)
        ESP = ESP + EA(I) * SIN(X)
        ECP = ECP + EA(I) * COS(X)
30    CONTINUE
c
c        Separate ECC and PIE
c
      ECC = SQRT (ESP*ESP + ECP*ECP)
      PIE = ATAN2 (ESP,ECP)
c
c        Series and linear term for PSI, n.h. vernal equinox position
c        relative to fixed stars
c
      X = 0.
      DO 40 I=1,10
        X = X + RPS*PA(I) * SIN (RPS*PF(I)*TIME + RPD*PP(I))
40    CONTINUE
      PSI = RPS*50.439273*TIME + RPD*3.392506 + X
c
c        Set PREC, prograde angle from perihelion to n.h.vernal equinox.
c        Add pi radians to PIE + PSI since it has been subtracted by
c        convention in Berger series.
c
      PREC = MOD (PIE + PSI + PI, TPI)
      IF (PREC.LT.0.) PREC = PREC + TPI
      PREC = TPI - PREC
c
c        Calculate VERN, 2 * pi * time from perihelion to n.h. vernal
c        equinox / perihelion-to-perihelion year. Use eccentric
c        anomaly E (Danby, Fundamentals of Celestial Mechanics,
c        1962, Eqs. (6.3.12) and (6.3.19).)
c
      ZS = SIN (0.5*PREC) / SQRT ((1.+ECC)/(1.-ECC))
      ZC = COS (0.5*PREC)
      ZE = 2. * ATAN2 (ZS,ZC)
      VERN = ZE - ECC * SIN(ZE)
c
      RETURN
      END

c----------------------------------------------------------------------

      subroutine zenorb91 (iyear, ecc, obliq, prec, vern)

c        Reads orbital parameters from tables by Berger and Loutre,1991,
c        QSR (via NOAA/NGDC). Only has values from 1950 to 5 Ma BP.

c     iyear = calendar year (bp) (supplied)
c     ecc   = eccentricity (returned)
c     obliq = obliquity in radians (returned)
c     prec  = precession, ie, prograde angle from perihelion to n.h.
c             vernal equinox, 0 to 2*pi radians (returned)
c     vern  = 2*pi *  time from perihelion to n.h. vernal equinox
c             / perihelion-to-perihelion year (returned)

      parameter (pi = 3.14159265358979)
      parameter (norb=5000)
      dimension eccorb(0:norb), oblorb(0:norb), precorb(0:norb)
      logical first
      save eccorb, oblorb, precorb, first
      character*200 cfile
      data first /.true./

      iu = 103
      if (first) then
        open (iu, file = cfile, status='old')
        read(iu,'(//)')
        do m=0,norb
          read(iu,'(6x,f10.0,2f8.0)') eccorb(m), precorb(m),  oblorb(m)
          precorb(m) = 180.-precorb(m)
          if (precorb(m).lt.0.) precorb(m) = precorb(m) + 360.0
        enddo
        close (iu)
        first = .false.
      endif

c        Time in berger file is years from 1950.0 ad.
c        In gcm: if iyear is positive, it is in "ad". If zero,
c               it is in "bp", i.e., relative to 1950 ad.
c        Here, always bp.

c     if (iyear.ge.0) then
c       zyear = iyear - 1950
c     else
        zyear = iyear
c     endif

c        Data runs from 5 Ma BP to present (0 BP),
c        with data every 1k years

      zyear = max (0., min (5000000., -zyear)) / 1000.
      iom = min (int(zyear), norb-1)
      fop = zyear - iom
      iop = iom + 1

      ecc   = (1.-fop)*eccorb(iom) + fop*eccorb(iop)
      obliq = (1.-fop)*oblorb(iom) + fop*oblorb(iop)

      if (precorb(iop)-precorb(iom).gt.180.) then
        zpm = precorb(iom) + 360.
        zpp = precorb(iop)
      else if (precorb(iop)-precorb(iom).lt.-180.) then
        zpm = precorb(iom)
        zpp = precorb(iop) + 360.
      else
        zpm = precorb(iom)
        zpp = precorb(iop)
      endif
      if (abs(zpp-zpm).gt.180.) then
        write (6,500) zpm,zpp,iom,iop
 500    format(/'*** Error (zenorb91): zpm,zpp,iom,iop=',2f8.2,2i6)
        stop
      endif
      prec  = (1.-fop)*zpm + fop*zpp
      prec = mod (prec, 360.)

      obliq = obliq * pi/180.
      prec  = prec  * pi/180.

c        Calculate VERN, 2 * pi * time from perihelion to n.h. vernal
c        equinox / perihelion-to-perihelion year. Use eccentric
c        anomaly E (Danby, Fundamentals of Celestial Mechanics,
c        1962, eqs. (6.3.12) and (6.3.19).)

      zs = sin (0.5*prec) / sqrt ((1.+ecc)/(1.-ecc))
      zc = cos (0.5*prec)
      ze = 2. * atan2 (zs,zc)
      vern = ze - ecc * sin(ze)

      return
      end

c----------------------------------------------------------------------

      SUBROUTINE ZENCAL (IYEAR, ISECYR, ISECDY, DT,
     *                   ECC, OBLIQ, PREC,VERN,DIST, ECCF,
     *                   COSQ, FRAQ, COSQ24, FRAQ24, ALAT, NLAT)
c
c        Calculates COSQ(J) and FRAQ(J), cos (zenith angle),
c        and daylight fraction in interval DT, versus latitude.
c        Also calculates 24-hr means in COSQ24 and FRAQ24.
c        Also calculates DIST, current earth-sun distance relative to
c        the semi-major axis, and eccentricity factor ECCF = 1/DIST**2.
c
c     IYEAR = year number (not used) (supplied)
c     ISECYR= secs into current calendar year from 00:00 Jan 1st (supp)
c     ISECDY= secs into current calendar day  from 00:00 (supplied)
c     DT    = time interval (must be .le. 1 day) (seconds) (supplied)
c     ECC   = eccentricity (supplied)
c     OBLIQ = obliquity (radians) (supplied)
c     PREC  = precession (prograde angle between perihelion and n.h.
c             vernal equinox) (0 to 2*pi) (supplied)
c     VERN  = 2*pi *  time from perihelion to n.h. vernal equinox
c             / perihelion-to-perihelion year (supplied)
c     DIST  = current earth-sun distance / semi-major axis (returned)
c     ECCF  = eccentricity factor (1/DIST**2) (returned)
c     COSQ  = cos (zenith angle) averaged over DT vs lat (returned)
c     FRAQ  = fraction of DT with daylight vs lat (returned)
c     COSQ24= 24-hour mean of cos(zenith angle) vs lat (returned)
c     FRAQ24= 24-hr fraction with daylight vs lat (returned)
c     ALAT  = latitudes (radians, cannot be pi/2 or -pi/2) (supplied)
c     NLAT  = latitudinal grid size (supplied)
c
c     Local variables:
c     SECPD = day length (constant, approx midnight-to-midnight) (secs)
c     AYEAR = perihelion-to-perihelion year (anomalistic) (days)
c     TYEAR = equinox-to-equinox year (tropical) (days)
c
      DIMENSION COSQ(NLAT), FRAQ(NLAT), COSQ24(NLAT), FRAQ24(NLAT),
     *          ALAT(NLAT)
      DIMENSION NDAYPM(12)
      PARAMETER (PI = 3.14159265358979)
c
      DATA NDAYPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
c
      DATA SECPD, AYEAR, TYEAR / .864E5, 365.25964, 365.24219 /
c
      DATA EPSORB / 1.E-7 /    ! or 1.e-5 for real*4
c
      TPI = 2.*PI
      PIH = PI/2.
      RPD = PI/180.
c
c        Calculate time from present to n.h.vernal equinox, by computing
c        time T from Jan 1st to present, then subtracting time TVE
c        from Jan 1st to vernal equinox. TVE is incremented each year
c        from 1989 by the difference between the tropical year TYEAR
c        (equinox-to-equinox) and 365 days, and decremented by 1 day
c        following each leap year. We define the calendar by setting
c        N.H. vernal equinox at 12:00 GMT March 21, 79.5 days after
c        00:00 GMT Jan 1 (as recommended for 6K BP by PMIP)
c
      VERNCAL = 79.5
c
      T = ISECYR
      TVE = VERNCAL
c     However, model ignores leap years, so comment out the
c     adjustment to TVE from 1989.
c     TVE = TVE + (IYEAR-1989)*(TYEAR-365.) - (IYEAR-1989)/4
      T = T - SECPD*TVE
c
c        Change to "time" from perihelion, expressed as
c        2*pi * fraction of perihelion-to-perihelion year
c
      T = VERN + TPI*T/(AYEAR*SECPD)
      IF (T.LT.0.) T = T + TPI
c
c        Solve for THETA, angle between perihelion and current position,
c        using Newton's method for ecc anomaly E (Danby, Fundamentals
c        of Celestial Mechanics, 1962, Eqs. (6.3.12) and (6.3.19).)
c
      E = T
      DO 20 ILOOP=1,100
        DE = - (E-ECC*SIN(E)-T) / (1.-ECC*COS(E))
        IF (ABS(DE).LE.EPSORB) GOTO 22
        E = E + DE
   20 CONTINUE
   22 CONTINUE
c
      ZS = SIN (0.5*E) * SQRT((1.+ECC)/(1.-ECC))
      ZC = COS (0.5*E)
      THETA = 2. * ATAN2 (ZS,ZC)
c
c        Calculate DIST (earth-sun distance / semi-major axis) and ECCF
c
      DIST = (1.-ECC*ECC) /  (1. + ECC*COS(THETA))
      ECCF = 1./(DIST*DIST)
c
c        Calculate B (angle between Earth spin axis and Earth-Sun line)
c        and CB = COS(B), SB = SIN(B)
c
      CB = SIN(THETA-PREC) * COS(PIH-OBLIQ)
c     Need ACOS to return 0 to pi for SB to have the right sign
      B = ACOS(CB)
      SB = SIN(B)
c
c        Loop over latitudes (nb: ALAT cannot be pi/2 or -pi/2)
c
      DO 100 J=1,NLAT

        CL = COS (PIH-ALAT(J))
        SL = SIN (PIH-ALAT(J))
c
c          Set RISE and SET, hour angles (between 0 and 2*pi, from
c          local midnight) of sunrise and sunset
c
        IF (CL*CB .GT. SL*SB)  THEN
c         Polar day
          X = PI
        ELSE IF (CL*CB .LT. -SL*SB)  THEN
c         Polar night
          X = 0.
        ELSE
c         sunrise and sunset occur
c         (need ACOS to return 0 to pi)
          X = ACOS(-CL*CB/(SL*SB))
        ENDIF
        RISE = PI - X
        SET  = PI + X
c
c          Set 24-hr total cos(zen) TCOSQ24 and daylight time TIMLQ24
c
        TCOSQ24 = 2.* (CL*CB*X + SL*SB*SIN(X))
        TIMLQ24 = SET - RISE
c
c          For an arbitrary longitude, set hour angles (0 to 2*pi,
c          from local  midnight) of T1 and T2, the start and end of
c          interval DT.
c
        ALON = 0. ! (greenwich, arbitrarily)
        T1 = ALON + ISECDY*TPI/SECPD
        T2 = T1 + DT*TPI/SECPD
        T1 = MOD (T1+TPI, TPI)
        T2 = MOD (T2+TPI, TPI)
c
c          Flag if local midnight occurs during DT,
c          and if so switch T1 and T2. Then will do calcs below for
c          24 hrs - DT, and convert back at end. First test for daily
c          mean case (DT = 24 hrs) and handle similarly.
c
        IFLAG = 1
        X = ABS ( MOD (T2-T1, TPI) )
        IF (X.LT.EPSORB .OR. ABS(X-TPI).LT.EPSORB) THEN
c         daily mean (DT = 24hrs)
          IFLAG = -1
          T2 = T1
        ELSE IF (T1. GT. T2) THEN
c         midnight occurs during DT (and not daily mean)
          IFLAG = -1
          X = T1
          T1 = T2
          T2 = X
        ENDIF
c
c          Constrain integ limits A1, A2 to be between RISE and SET,
c          and integrate to get total cos(zen) TCOSQ and daylight time
c          TIMLQ. Subtract PI from A1,A2 to get in range -PI to PI.
c
        A1 = MAX (RISE, MIN (SET, T1)) - PI
        A2 = MAX (RISE, MIN (SET, T2)) - PI
        TCOSQ = CL*CB*(A2-A1) + SL*SB*(SIN(A2)-SIN(A1))
        TIMLQ = A2 - A1
c
c          If local midnight occurs within DT, subtract from 24-hr tots
c
        IF (IFLAG.EQ.-1) THEN
          TCOSQ = TCOSQ24 - TCOSQ
          TIMLQ = TIMLQ24 - TIMLQ
        ENDIF
c
c          Set avg daylight cos(zen) COSQ and daylight fraction FRAQ
c          for time interval DT
c
        IF (TIMLQ.GT.0.) THEN
          COSQ(J) = TCOSQ / TIMLQ
        ELSE
          COSQ(J) = 0.
        ENDIF
        FRAQ(J) = TIMLQ / (TPI*DT/SECPD)
c
c          Set 24-hr mean cos(zen) and 24-hr daylight fraction
c
        COSQ24(J) = TCOSQ24 / TPI
        FRAQ24(J) = TIMLQ24 / TPI
c
  100 CONTINUE
c
      RETURN
      END

