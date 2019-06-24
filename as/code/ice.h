! source file: /net/mare/home1/eby/as/updates/ice.h
!======================== include file "ice.h" =========================

!     arrays for the ice model
!     hice     = ice thickness (tau-1, tau, tau+1)
!     aice     = ice area (tau-1, tau, tau+1)
!     tice     = ice surface temperature
!     hsno     = effective thickness of snow
!     psno     = precipitation as snow (+ accumulation, - melt)
!     frzpt    = freezing temperature of sea water
!   cpts ice model ("ice_cpts")
!     A        = area
!     heff     = effective thickness of ice
!     Ts       = ice surface temperature
!     groice   = thermo ice thickness change from diagnostic
!     E        = enthalpy of each ice layer
!     strength = pressure (from Rothrock), must have ncat > 3
!     hseff    = effective thickness of snow
!   ice dynamics ("ice_evp")
!     uice     = u ice velocity
!     vice     = v ice velocity
!     pice     = pressure due to internal stress
!     xint     = x component of ice interaction
!     yint     = y component of ice interaction
!     del      = delta for ice dynamics models
!     eI       = divergence of ice velocity
!     nseg     = number of domains needed for ice calculations
!     jsi      = start of limited domain for ice calculations
!     jei      = end of limited domain for ice calculations
!   brine convection ("convect_brine")
!     cbf      = flux from rejected brine over ice categories
!     cba      = area of brine rejection for each categories
!     cba0     = area of the cell without brine rejection
!   land ice ("landice_data")
!     hicel    = land ice thickness
!     aicel    = land ice area

      integer nseg, jsi, jei

      real hice, aicelan, aiceocn, ticel, ticeo, hsnol, hsnoo
      real psnol, psnoo, frzpt
      real A, heff, Ts, groice, E, hseff, strength
      real uice, vice, pice, xint, yint, del, eI
      real cbf, cba, cba0, hicel, aicel

      common /ice_r/ hice(imt,jmt,2)
      common /ice_r/ aicelan(imt,jmt,2), aiceocn(imt,jmt,2)
      common /ice_r/ ticel(imt,jmt), ticeo(imt,jmt)
      common /ice_r/ hsnol(imt,jmt,2), hsnoo(imt,jmt,2)
      common /ice_r/ psnol(imt,jmt), psnoo(imt,jmt)
      common /ice_r/ frzpt(imt,jmt)
      integer ncat
      parameter (ncat=1)
      common /ice_r/ uice(imt,jmt), vice(imt,jmt), pice(imt,jmt)
      common /ice_r/ xint(imt,jmt), yint(imt,jmt)
      common /ice_r/ del(imt,jmt), eI(imt,jmt)
      common /ice_i/ nseg, jsi(jmt), jei(jmt)
!   time averaged arrays
!     ta_hsno     = time average snow thickness
!     ta_hice     = time average ice thickness
!     ta_aice     = time average ice area
!     ta_tice     = time average ice surface temperature
!   ice dynamics ("ice_evp")
!     ta_uice     = time average u ice velocity
!     ta_vice     = time average v ice velocity
!     ta_pice     = time average pressure
!     ta_xint     = time average x component of ice interaction
!     ta_yint     = time average y component of ice interaction
!   land ice ("landice_data")
!     ta_hicel    = time average land ice thickness
!     ta_aicel    = time average land ice area

      real ta_hice, ta_aice, ta_tice, ta_hsno, ta_heff, ta_A, ta_Ts
      real ta_hseff, ta_uice, ta_vice, ta_pice, ta_xint, ta_yint
      real ta_aicel, ta_hicel

      common /ice_r/ ta_hice(imt,jmt), ta_aice(imt,jmt)
      common /ice_r/ ta_tice(imt,jmt), ta_hsno(imt,jmt)
      common /ice_r/ ta_uice(imt,jmt), ta_vice(imt,jmt)
      common /ice_r/ ta_pice(imt,jmt), ta_xint(imt,jmt)
      common /ice_r/ ta_yint(imt,jmt)

