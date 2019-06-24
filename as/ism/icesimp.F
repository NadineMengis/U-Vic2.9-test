
c-----------------------------------------------------------------------

      subroutine climate_simp (budgsnow, budgrain, budgevap, budgmelt,
     *                         tsurf, tsurfi,
     *                         h, hs, maskwater, 
     *                         sealev, dtantann, dtantjan, dtseas,
     *                         rco2, ecc, obl, prec,
#if defined CO2INTER
     *                         tempmon, precmon,
#endif
     *                         timeice, weirun)

c     Sets surface climate from relatively simple parameterizations

      use comicephys
      use comicegrid

c     passed:
      dimension 
     *  budgsnow(nx,ny),     budgrain(nx,ny),
     *  budgevap(nx,ny),     budgmelt(nx,ny),
     *  tsurf(nx,ny),        tsurfi(nx,ny),
     *  h(nx,ny),            hs(nx,ny), 
     *  maskwater(nx,ny)
#if defined CO2INTER
     * ,tempmon(nx,ny,12),   precmon(nx,ny,12)
#endif

#if defined EISANTA &&  defined BASACCUM
      logical firstbas
      data firstbas /.true./
      dimension budgbas (nx,ny)
      save firstbas, budgbas
#endif

#if defined EISANTA
c    compute current ice-sheet vol and area for sfc feedback:
      ztotari = 0.
      ztotice = 0.
      do j=1,ny
        do i=1,nx
          if (h(i,j).gt.0.001) ztotari = ztotari + darea(i,j)
          ztotice = ztotice + h(i,j)*darea(i,j)
        enddo
      enddo
      ztotari = ztotari*1.e-6 ! m2 to km2
      ztotice = ztotice*1.e-9 ! m3 to km3
#if defined BASACCUM
      if (firstbas) then
        write (6,'(a)') 'Reading BAS modern surface accumulation data'  
        call readbas (alond, alatd, budgbas, nx, ny, iubas)
        firstbas = .false.
      endif
#endif
#endif

      call zero (budgsnow, nx*ny)
      call zero (budgrain, nx*ny)
      call zero (budgevap, nx*ny)
      call zero (budgmelt, nx*ny)

c================
      do j=1,ny
        do i=1,nx
c================

c-----------------------------------------------------
#if defined EIS1 || defined EIS2 || defined EIS3
c-----------------------------------------------------
          tsurf(i,j)  = tmelt
          budgsnow(i,j) = 0.

c-------------------
#elif defined EIS4
c-------------------
          tsurf(i,j)  = tmelt
          if (h(i,j).gt.0.) then
            budgsnow(i,j) = 0.2
          else
            budgsnow(i,j) = 0.
          endif

c-------------------
#elif defined EISA
c-------------------
          zr = .001 * max (abs(xh(i,j)), abs(yh(i,j)))
          tsurf(i,j)  = 239. + 8.e-8*(zr**3)
          budgsnow(i,j) = 0.3
#if defined EISVAR20
          tsurf(i,j)    = tsurf(i,j)    + 10.*sin(2.*pi*timeice/20000.)
          budgsnow(i,j) = budgsnow(i,j) + 0.2*sin(2.*pi*timeice/20000.)
#elif defined EISVAR40
          tsurf(i,j)    = tsurf(i,j)    + 10.*sin(2.*pi*timeice/40000.)
          budgsnow(i,j) = budgsnow(i,j) + 0.2*sin(2.*pi*timeice/40000.)
#endif
c         zero at edge (h=0) to avoid unif h first sparse call,if shelf:
          if ( (nx.gt.1 .and. (i.eq.1 .or. i.eq.nx)) .or.
     *         (ny.gt.1 .and. (j.eq.1 .or. j.eq.ny)) ) budgsnow(i,j) =0.

c-------------------
#elif defined EISB
c-------------------
          zr   = .001 * sqrt(xh(i,j)**2 + yh(i,j)**2)
          zrel = 450.
          tsurf(i,j) = 270. -.01*hs(i,j)
#if defined EISVAR20
          zrel       = zrel          + 100.*sin(2.*pi*timeice/20000.)
          tsurf(i,j) = tsurf(i,j)    + 10. *sin(2.*pi*timeice/20000.)
#elif defined EISVAR40
          zrel       = zrel          + 100.*sin(2.*pi*timeice/40000.)
          tsurf(i,j) = tsurf(i,j)    + 10. *sin(2.*pi*timeice/40000.)
#endif
          budgsnow(i,j) =  min (0.5, .01*(zrel-zr))

c-------------------------------------
#elif defined EISC || defined EISD
c-------------------------------------
          zr   = .001 * sqrt(xh(i,j)**2 + yh(i,j)**2)
          tsurf(i,j) = 238.15 + .0167*zr   ! expt A,G
c         tsurf(i,j) = 223.15 + .0167*zr   ! expt F
          budgsnow(i,j) =  min (0.5, .01*(450.-zr))

c----------------------
#elif defined EISLINE
c----------------------

#if defined LINEB
c         tsurf(i,j) = tmelt - 15. - 10.*weirun - .0051*hs(i,j)
c         budgsnow(i,j) = 0.15 - 0.10*weirun
          tsurf(i,j) = tmelt - 20.
          budgsnow(i,j) = 0.10
#elif defined LINEC
          tsurf(i,j) = tmelt - 15. - 10.*weirun - .0051*hs(i,j)
          budgsnow(i,j) = 0.15 * (2.**((tsurf(i,j)-tmelt+15.)/10.))
c         budgsnow(i,j) = 0.15*exp(-hs(i j)/3000.)
#elif defined LINED
          tsurf(i,j) = tmelt - 15. - 10.*weirun - .0051*hs(i,j)
          budgsnow(i,j) = 0.15 * (2.**((tsurf(i,j)-tmelt+15.)/10.))
          if (weirun.eq.0.) budgsnow(i,j) = -1.
#elif defined LINEE
          tsurf(i,j) = tmelt - 15. - 10.*weirun - .0051*hs(i,j)
          budgsnow(i,j) =  0.
#elif defined LINEF
          tsurf(i,j) = tmelt - 15. - 10.*weirun - .0051*hs(i,j)
c         budgsnow(i,j) = 0.30 ! schoof
c         budgsnow(i,j) = 0.20
          budgsnow(i,j) = 0.15
c         budgsnow(i,j) =  max(-1.,min(.15,.0005*(hs(i,j)-500.) ))
c         if (timeice.lt.50000.) then
c           budgsnow(i,j) = 0.10
c         else if (timeice.lt.100000.) then
c           budgsnow(i,j) = 0.15
c         endif
c         budgsnow(i,j) = 0.10 ! .15
c         if (timeice.lt.20000.) then
c           zr = xh(i,j)
c           if (zr.gt.500.e3) budgsnow(i,j) = -1.0
c         endif
          if (nx.gt.ny) then
            if (i.eq.nx) budgsnow(i,j) = 0.
          else
            if (j.eq.ny) budgsnow(i,j) = 0.
          endif
#elif defined LINEG
          tsurf(i,j) = tmelt - 20.
          budgsnow(i,j) =  0.20
c         budgsnow(i,j) =  0.   ! 999_1pt
#elif defined TRANSECTA
          tsurf(i,j) = tmelt - 20.
          if (timeice.lt.50000.) then
            budgsnow(i,j) =  0.10                           ! lgm
c           budgsnow(i,j)=  max(-1.,min(.2,.0005*(hs(i,j)-500.) ))!mis31
          else if (timeice.lt.100000.) then
            budgsnow(i,j) =  0.15                           ! mod
          else 
            budgsnow(i,j)=  max(-1.,min(.2,.0005*(hs(i,j)-500.) ))!mis31
c           budgsnow(i,j) =  0.15                           ! mod 
c           budgsnow(i,j) =  0.10                           ! lgm
          endif
#endif

c----------------------
#elif defined EISANTA
c----------------------

c        For most experiments, budgsurf represents whole surface budget,
c        and budgrain, budgmelt, etc, are zero. But for EISANTA, 
c        budgsnow is snowfall, budgrain is rain, budgmelt is surface 
c        melt. (Nb: budgrain and budgmelt contribute to baseperc in 
c        icetherm/vdif, and baseperc is included in budgall in icectl).
c        (budgevap is zero throughout).


#if defined FORCEEO
          tsurf(i,j) = tmelt         - .0050*hs(i,j)              ! 8888
     *                               - .35* (abs(alatd(i,j))-80.) ! 8888
     *                               - 20.                        ! 8888
c    *                               + dtantann
     *                               + dtantjan                   ! 8888
c                                      antarc warming per CO2 doubling:
     *                               + 10.*(log(rco2)/log(2.))    ! 8888
c                                      ice sheet albedo feedback:
c    *                               + 10.*(1.-ztotice/26.e6)
c    *                               + 10.*(1.-ztotari/13.e6)
     *                               + 10.                        ! 8888
#else
c         EISMINT II:
          tsurf(i,j) = tmelt + 34.46 - .00914*hs(i,j)
     *                               - .68775*abs(alatd(i,j))
     *                               + dtantann
     *                               + 10.*sealev/125. 
#endif

c         Ritz et al (2001):
c         if (hs(i,j).le.200.) then
c           tsurf(i,j) = tmelt + 49.642 - .943*abs(alatd(i,j))
c
c         else if (hs(i,j).gt.200. .and. hs(i,j).le.1500) then
c           tsurf(i,j) = tmelt + 36.689 - .725*abs(alatd(i,j))
c    *                                  - .005102*hs(i,j)
c
c         else if (hs(i,j).ge.1500.) then
c           tsurf(i,j) = tmelt +  7.405 - .180*abs(alatd(i,j))
c    *                                  - .014285*hs(i,j)
c         endif

#if defined BASACCUM
          zprecip = budgbas(i,j) * (2.**(dtantann/10.))
#else
          zprecip = 1.5 * (2.**((tsurf(i,j)-tmelt)/10.))
#endif

c           Compute ablation, based on pdd with annual seasonal cycle
c           Also modify precip for fraction falling as snow (zsn) 

          zann = tsurf(i,j)
#if defined FORCEEO
c         zamp = 11. +  3.*(dtseas-22.)/6. ! 8888
          zamp = 12.                       ! 8888
#else
          zamp = 0.5*dtseas
#endif

#if defined CO2INTER
          do m=1,12
            zt = tsurf(i,j) + zamp*cos(2.*pi*(m-1)/12.)
            zp = 1.5 * (2.**((zt-tmelt)/10.))
            tempmon(i,j,m) = zt - tmelt                  ! K to C
            precmon(i,j,m) = zp * 1000./365.             ! m/y to mm/day
          enddo
#endif

          if (tmelt-zann.le.-zamp) then
            pdd = (zann-tmelt)*365.
            zsn = 0.
          else if (tmelt-zann.ge. zamp) then
            pdd = 0.
            zsn = 1.
          else 
            ztim = acos ((tmelt-zann)/zamp)
            pdd  = ((zann-tmelt)*ztim + zamp*sin(ztim)) * (365./pi)
            zsn = 1. - ztim/pi
          endif
          zabl = .005 * pdd

          budgsnow(i,j) = zsn*zprecip 
          budgrain(i,j) = (1.-zsn)*zprecip 
          budgmelt(i,j) = zabl

c         ensure zero or slightly negative budget over open ocean:
          if (maskwater(i,j).eq.1) then
            zf = max(0., min (1.,  (h(i,j)-2.5)/2.5))              ! 777
            budgsnow(i,j) = zf * budgsnow(i,j) + (1.-zf)*(0   )    ! 777
            budgrain(i,j) = zf * budgrain(i,j) + (1.-zf)*(0.  )    ! 777
            budgmelt(i,j) = zf * budgmelt(i,j) + (1.-zf)*(0.02)    ! 777
          endif

c--------------------
#elif defined HEINO
c--------------------

          if (nx.eq.81) then
            nxmid = (nx+1)/2
          else
            nxmid = 1
          endif
          zd = 1.e-3 * sqrt (   (xh(i,j)-xh(nxmid,(ny+1)/2))**2  
     *                        + (yh(i,j)-yh(nxmid,(ny+1)/2))**2
     *                      )
#if defined HEINO_T1
          ztsmin = 223.15
#elif defined HEINO_T2
          ztsmin = 243.15
#else
          ztsmin = 233.15
#endif
          tsurf(i,j) = ztsmin + 2.5e-9*(zd**3)

#if defined HEINO_B1
          zbmin = 0.075
          zbmax = 0.15
#elif defined HEINO_B2
          zbmin = 0.30
          zbmax = 0.60
#else
          zbmin = 0.15
          zbmax = 0.30
#endif
          budgsnow(i,j) = zbmin + (zbmax-zbmin)*(zd/2000.)

c------------------
#elif defined HOM
c------------------

          tsurf(i,j) = tmelt
          budgsnow(i,j) = 0.

c----
#else
c----

          tsurf(i,j)  = tmelt

c         if (abs((j-0.5) - 0.5*ny).lt.10.) then
c           budgsnow(i,j) = 1.
c         else
c           budgsnow(i,j) = -10.
c         endif

          zr = sqrt (xh(i,j)**2 + yh(i,j)**2)
c         ze = dd0*(ny/4.)                           ! stable
          ze = 0.6 * max(abs(xh(1,1)),abs(yh(1,1)))  ! unstable
          if (zr.le.ze) then
c         zx = abs(xh(i,j))
c         zy = abs(yh(i,j))
c         if (zx.le.ze .and. zy.le.ze) then
            budgsnow(i,j) = .15
          else
            budgsnow(i,j) = -.3
          endif
c-----
#endif
c-----

c            For effective sfc air temp above open ocean, for possible
c            use in icetherm, rectify for seasonal cycle (assumed 
c            instantaneously insulated by sea ice when below freezing).
c            i.e., annual mean of max (tmelt, zann+zamp*cos(t)):
c
c         For annual mean air temp (zann), either use tsurf (set above),
c         or use linear fit vs lat, obs (Peixoto and Oort, tsurf.cgm).
c         zann = tsurf(i,j)
c         zann = -2. - 25.*(abs(alatd(i,j)) - 60.)/30. + tmelt
c         For seasonal half-amplitude (zamp), base on Peixoto and Oort.
c         zamp = 10. + 0. *(abs(alatd(i,j)) - 60.)/30.
c         if (tmelt-zann.le.-zamp) then
c           tsurfo(i,j) = zann
c         else if (tmelt-zann.ge. zamp) then
c           tsurfo(i,j) = tmelt 
c         else 
c           ztim = acos ((tmelt-zann)/zamp)
c           tsurfo(i,j) = (zann*ztim + zamp*sin(ztim) + (pi-ztim)*tmelt)
c    *                    / pi
c         endif

c            Similarly for effective surface temp of ice/snow (used in 
c            icetherm): ann mean of min (tmelt, zann+zamp*cos(t)):

c         tsurfi(i,j) = min (tmelt, tsurf(i,j))
          zann = tsurf(i,j)
          zamp = 10.
          if (tmelt-zann.le.-zamp) then
            tsurfi(i,j) = tmelt
          else if (tmelt-zann.ge. zamp) then
            tsurfi(i,j) = zann 
          else 
            ztim = acos ((tmelt-zann)/zamp)
            tsurfi(i,j) = (zann*(pi-ztim) - zamp*sin(ztim) + ztim*tmelt)
     *                    / pi
          endif

        enddo
      enddo

      return
      end
