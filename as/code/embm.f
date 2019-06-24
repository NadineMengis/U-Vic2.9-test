! source file: /net/mare/home1/eby/as/updates/embm.F
      subroutine embm (is, ie, js, je)

!=======================================================================

!     The atmospheric energy moisture balance model (EMBM)

!     Fanning, A.F. and A.J. Weaver, An atmospheric energy-moisture
!       balance model: climatology, interpentadal climate change,
!       and coupling to an ocean general circulation model,
!       J. Geophys. Res., 101, 15,111-15,128, 1996
!=======================================================================

      implicit none

      integer ie, is, je, js, n

      include "size.h"
      include "csbc.h"
      include "cembm.h"
      character(120) g_st
      include "atm.h"

!-----------------------------------------------------------------------
!     increment counter and set the time step type
!-----------------------------------------------------------------------

      nats = nats + 1
      if (nats .gt. namix) then
        lf = 2
        dts = dtatm
        nats = 1
      else
        lf = 1
        dts = 2.0*dtatm
      endif

      addflxa = .true.
      if (mod(nats,2) .ne. 0) addflxa = .false.

!-----------------------------------------------------------------------
!     calculate fluxes at tau
!-----------------------------------------------------------------------

      call fluxes (is, ie, js, je)

!-----------------------------------------------------------------------
!     compute ice fluxes at tau and ice thickness and area at tau+1
!-----------------------------------------------------------------------

      call ice (is, ie, js, je)

!-----------------------------------------------------------------------
!     compute atmospheric tracers at tau+1. start with humidity so that
!     the precipitation flux can be calculated for latent heat
!-----------------------------------------------------------------------

      call solve (2)
      call precipitate (is, ie, js, je)
      call solve (1)
      do n=3,nat
        call solve (n)
      enddo

!-----------------------------------------------------------------------
!     calculate the total atmospheric fluxes for coupling
!-----------------------------------------------------------------------

      if (addflxa) call sum_flux (is, ie, js, je)

!-----------------------------------------------------------------------
!     accumulate time averages
!-----------------------------------------------------------------------
      call ta_embm_tavg (is, ie, js, je, 1)

!-----------------------------------------------------------------------
!     accumulate time step integrals
!-----------------------------------------------------------------------
      call ta_embm_tsi (is, ie, js, je, 1)

      return
      end

      subroutine sum_flux (is, ie, js, je)

!=======================================================================
!     sum fluxes over atmospheric time steps
!=======================================================================

      implicit none

      integer i, ie, iem1, is, isp1, j, je, jem1, js, jsp1

      real fa, fb

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "atm.h"
      include "ice.h"
      include "csbc.h"

      isp1 = is + 1
      iem1 = ie - 1
      jsp1 = js + 1
      jem1 = je - 1

      do j=jsp1,jem1
        do i=isp1,iem1
          if (tmsk(i,j) .gt. 0.0) then
            fluxo(i,j,isat) = fluxo(i,j,isat) + dts*(
     &                	     dnswro(i,j)      *tmsk(i,j)
     &                     - uplwro(i,j)      *tmsk(i,j)
     &                     - upltnto(i,j)     *tmsk(i,j)
     &                     - upsenso(i,j)     *tmsk(i,j))
            fluxo(i,j,ishum) = fluxo(i,j,ishum) + dts*
     & 	                       ((precip(i,j)   *tmsk(i,j))
     &                      -  ( evapo(i,j)    *tmsk(i,j)))
          endif

	  if (tmsk(i,j) .lt. 1. .and. land_map(i,j) .ne. 0) then
            sbc(i,j,iat) = sbc(i,j,iat) + dts*(at(i,j,2,1)
     &                   - elev(i,j)*rlapse)
            sbc(i,j,irh) = sbc(i,j,irh) + dts*rh(i,j)
            sbc(i,j,iaws) = sbc(i,j,iaws) + dts*sbc(i,j,iws)
            sbc(i,j,iswr) = sbc(i,j,iswr) + dts*dnswrl(i,j)
            sbc(i,j,ipr) = sbc(i,j,ipr) + dts*precip(i,j)
            if (psnol(i,j) .ge. 0.) then
              sbc(i,j,ips) = sbc(i,j,ips) + dts*psnol(i,j)
              sbc(i,j,ipr) = sbc(i,j,ipr) - dts*psnol(i,j)
            endif
	  endif
          if (tmsk(i,j) .lt. 1..and. land_map(i,j) .eq. 0 ) then

            sbc(i,j,iro) = sbc(i,j,iro) + dts*runoff(i,j)
	  endif
          if (umsk(i,j) .ge. 0.5) then
            fluxo(i,j,nat+1) = fluxo(i,j,nat+1) + dts*sbc(i,j,itaux)
            fluxo(i,j,nat+2) = fluxo(i,j,nat+2) + dts*sbc(i,j,itauy)
            fluxo(i,j,nat+1) = fluxo(i,j,nat+1) + dts*xint(i,j)
            fluxo(i,j,nat+2) = fluxo(i,j,nat+2) + dts*yint(i,j)
          endif
          atbar(i,j) = atbar(i,j) + dts*at(i,j,2,isat)
        enddo
      enddo

      totaltime = totaltime + dts
      atatm = atatm + dts

      return
      end
