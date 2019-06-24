! source file: /net/mare/home1/eby/as/ism/icesimp.F

c-----------------------------------------------------------------------

      subroutine climate_simp (budgsnow, budgrain, budgevap, budgmelt,
     *                         tsurf, tsurfi,
     *                         h, hs, maskwater,
     *                         sealev, dtantann, dtantjan, dtseas,
     *                         rco2, ecc, obl, prec,

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

      call zero (budgsnow, nx*ny)
      call zero (budgrain, nx*ny)
      call zero (budgevap, nx*ny)
      call zero (budgmelt, nx*ny)

c================
      do j=1,ny
        do i=1,nx
c================

c-----------------------------------------------------
c----------------------

c        For most experiments, budgsurf represents whole surface budget,
c        and budgrain, budgmelt, etc, are zero. But for 1,
c        budgsnow is snowfall, budgrain is rain, budgmelt is surface
c        melt. (Nb: budgrain and budgmelt contribute to baseperc in
c        icetherm/vdif, and baseperc is included in budgall in icectl).
c        (budgevap is zero throughout).

c         EISMINT II:
          tsurf(i,j) = tmelt + 34.46 - .00914*hs(i,j)
     *                               - .68775*abs(alatd(i,j))
     *                               + dtantann
     *                               + 10.*sealev/125.

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

          zprecip = 1.5 * (2.**((tsurf(i,j)-tmelt)/10.))

c           Compute ablation, based on pdd with annual seasonal cycle
c           Also modify precip for fraction falling as snow (zsn)

          zann = tsurf(i,j)
          zamp = 0.5*dtseas

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
