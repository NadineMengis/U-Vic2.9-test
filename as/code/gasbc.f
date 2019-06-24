! source file: /net/mare/home1/eby/as/updates/gasbc.F
      subroutine gasbc (is, ie, js, je)

!=======================================================================
!     calculate boundary conditions for the atmospheric model
!=======================================================================

      use subgrid

      implicit none

      integer ie, is, je, js, i, iem1, isp1, j, jem1, jsp1, k, n

      real sss, sst, xconv, t_in, s_in, dic_in, ta_in, co2_in, pt_in
      real sit_in, atmpres, pHlo, pHhi, pH, co2star,  dco2star, pCO2surf
      real dpco2, CO3surf, Omega_c, Omega_a, scco2, piston_vel, avgflxc
      real calday, f, sco2, o2sat, o2sato, o2surf, piston_o2, cfc11ccn
      real cfc12ccn, wt, sccfc, piston_cfc, sol_cfc, cfcsat, ao, tarea
      real tdc14ccn, h_r, d, f1, f2, f3, f4, f5, area, C2K, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "csbc.h"
      include "mw.h"
      include "ice.h"
      include "switch.h"
      include "tmngr.h"
      include "cembm.h"
      include "atm.h"
      include "insolation.h"
      include "calendar.h"
      include "grdvar.h"
      include "levind.h"
      include "solve.h"
      include "mtlm.h"
      include "ism.h"

      include "temp.h"

      real cosz(is:ie,js:je)
      real dmsk(is:ie,js:je)

      character(120) g_st
      integer ntrec
      data ntrec /0/
      save ntrec

      ntrec=ntrec+1

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

!     xconv is constant to convert piston_vel from cm/hr -> cm/s
!     here it is 100.*a*xconv (100 => m to cm, a=0.337, xconv=1/3.6e+05)
      xconv = 33.7/3.6e+05
      C2K = 273.15
      phlo = 6.
      phhi = 10.
      sit_in = 7.6875e-03 !mol/m^3
      pt_in = 0.5125e-3   !mol/m^3
      atmpres = 1.0       !atm

!-----------------------------------------------------------------------
!     zero totals for new accumulation
!-----------------------------------------------------------------------
      atatm = 0.
      fluxo(:,:,:) = 0.
      fluxl(:,:,:) = 0.
      sbc(:,:,ihflx) = 0.
      sbc(:,:,isflx) = 0.
      sbc(:,:,iro) = 0.
      sbc(:,:,iat) = 0.
      sbc(:,:,irh) = 0.
      sbc(:,:,ipr) = 0.
      sbc(:,:,ips) = 0.
      sbc(:,:,iaws) = 0.
      sbc(:,:,iswr) = 0.
      sbc(:,:,idicflx) = 0.
      sbc(:,:,ialkflx) = 0.
      sbc(:,:,io2flx) = 0.
      sbc(:,:,ipo4flx) = 0.
      sbc(:,:,iphytflx) = 0.
      sbc(:,:,izoopflx) = 0.
      sbc(:,:,idetrflx) = 0.
      sbc(:,:,ino3flx) = 0.
      sbc(:,:,idiazflx) = 0.

!-----------------------------------------------------------------------
!     set solar constant
!-----------------------------------------------------------------------
      call solardata

!-----------------------------------------------------------------------
!     update insolation for the current day
!-----------------------------------------------------------------------
!     subroutine decl is expecting a 365.25 day year
      calday = dayoyr*365.25/yrlen
      call decl (calday, eccen, obliq, mvelp, lambm0, sindec, eccf)
      i = (ie-is+1)*(je-js+1)
!     get average zenith angle
      call zenith (i, c0, daylen, daylen, tlat, tlon, sindec, cosz)
      solins(is:ie,js:je) = solarconst*eccf*cosz(is:ie,js:je)

!-----------------------------------------------------------------------
!     update any atmospheric data
!-----------------------------------------------------------------------
      call data (is, ie, js, je)

!-----------------------------------------------------------------------
!     calculate winds with new feedback
!-----------------------------------------------------------------------
      call add_awind (is, ie, js, je)

!-----------------------------------------------------------------------
!     calculate freezing point of sea water using UNESCO (1983)
!-----------------------------------------------------------------------

      do j=jsp1,jem1
        do i=isp1,iem1

          if (tmsk(i,j) .gt. 0.0) then

            sss = 1000.0*sbc(i,j,isss) + 35.0
            frzpt(i,j) = -.0575*sss + 1.71e-3*sss**1.5 - 2.155e-4*sss**2
            sst = sbc(i,j,isst)
            ao = 1. - aiceocn(i,j,2)

!-----------------------------------------------------------------------
!           calculate ocean carbon fluxes
!-----------------------------------------------------------------------
            t_in = sst
            s_in = sss
            dic_in = sbc(i,j,issdic)
            ta_in = sbc(i,j,issalk)
            co2_in = co2ccn
            call co2calc_SWS (t_in, s_in, dic_in, ta_in, co2_in, pt_in
     &,                       sit_in, atmpres, pHlo, pHhi, pH, co2star
     &,                       dco2star, pCO2surf, dpco2, CO3surf
     &,                       Omega_c, Omega_a)
!           Schmidt number for CO2
            scco2 = 2073.1 - 125.62*sst + 3.6276*sst**2
     &            - 0.043219*sst**3
            piston_vel = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *((scco2/660.)**(-0.5))
!           dic in umol cm-3 or (mol m-3) => flux in umol cm-2 s-1
            sbc(i,j,idicflx) = piston_vel*dco2star

!-----------------------------------------------------------------------
!           calculate ocean oxygen fluxes
!-----------------------------------------------------------------------
!           Schmidt number for O2
            sco2 = 1638.0 - 81.83*sst + 1.483*sst**2 - 0.008004*sst**3
!           piston velocity for O2
            piston_o2 = ao*xconv*((sbc(i,j,iws)*0.01)**2)
     &                  *(sco2/660.0)**(-0.5)
!           oxygen saturation concentration [mol/m^3]
            f1 = alog((298.15 - sst)/(C2K + sst))
            f2 = f1*f1
            f3 = f2*f1
            f4 = f3*f1
            f5 = f4*f1
            o2sat = exp (2.00907 + 3.22014*f1 + 4.05010*f2
     &             + 4.94457*f3 - 2.56847E-1*f4 + 3.88767*f5
     &             + sss*(-6.24523e-3 - 7.37614e-3*f1 - 1.03410e-2*f2
     &             - 8.17083E-3*f3) - 4.88682E-7*sss*sss)
!           Convert from ml/l to mol/m^3
            o2sat = o2sat/22391.6*1000.0
            sbc(i,j,io2flx) = piston_o2*(o2sat - sbc(i,j,isso2))

          else
!-----------------------------------------------------------------------
!           calculate land carbon fluxes
!-----------------------------------------------------------------------
!           convert from kg m-2 s-1 => umol cm-2 s-1
            sbc(i,j,idicflx) = (sbc(i,j,inpp) - sbc(i,j,isr)
     &                       - sbc(i,j,iburn))*0.1/12.e-6
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     set boundary conditions for carbon
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,idicflx), imt, jmt)
      dmsk(:,:) = 1.
      call areaavg (sbc(1,1,idicflx), dmsk, avgflxc)

!-----------------------------------------------------------------------
!     set boundary conditions for oxygen
!-----------------------------------------------------------------------
      call setbcx (sbc(1,1,io2flx), imt, jmt)

!-----------------------------------------------------------------------
!     calculate CO2 forcing
!-----------------------------------------------------------------------
      call co2forc

!-----------------------------------------------------------------------
!     set flags to calculate new coefficients
!-----------------------------------------------------------------------
      newcoef(:,:) = .true.

!-----------------------------------------------------------------------
!     zero time averages if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. timavgperts) call ta_embm_tavg (is, ie, js, je, 0)

!-----------------------------------------------------------------------
!     zero time step integrals if not in an averaging period
!-----------------------------------------------------------------------
      if (.not. tsiperts) call ta_embm_tsi (is, ie, js, je, 0)
      !calculate subgridded bias fields
      call calcbias

      return
      end
