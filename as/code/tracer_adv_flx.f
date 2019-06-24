! source file: /net/mare/home1/eby/as/2.9.old3/source/mom/tracer_adv_flx.F
      subroutine adv_flux (joff, js, je, is, ie, n)

!=======================================================================
!     computes advective fluxes using a flux-corrected transport scheme

!        for reference see
!        Gerdes, R., C. Koeberle and J. Willebrandt, 1991
!        the influence of numerical advection schemes on the results of
!        ocean general circulation models. Clim Dynamics 5, 211-226
!        and
!        Zalesak, S. T., 1979: Fully multidimensional flux-corrected
!        transport algorithms for fluids. J. Comp. Phys. 31, 335-362.

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       jstrt = starting row in the MW for fct
!       jend  = ending row in the MW for fct
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       istrt = max(2,starting longitude index in the MW)
!       iend  = min(imt-1,ending longitude index in the MW)
!       n     = tracer index

!     output: ( via common mw in "mw.h" )
!       adv_fn = 2*advective flux across northern face of T-cell
!       adv_fe = 2*advective flux across eastern face of T-cell
!       adv_fb = 2*advective flux across bottom face of T-cell
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, istrt, is, iend, ie
      integer istrtm1, iendp1, jstrt, js, jend, je, joff, jlast, jp2
      integer jp1

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, aidif, totadv
      real fxa, fxb, pplus, pminus, qplus, qminus, reltim

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "grdvar.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "isopyc.h"

!---------------------------------------------------------------------
!     dimension local data
!---------------------------------------------------------------------

      real twodt(km), dcf(imt), Trmin(imt), Trmax(imt)
      real Cpos(imt), Cneg(imt), flxlft(imt), flxrgt(imt)
      real Rpl(imt,km), Rmn(imt,km), tmaski(imt,km,jmw)
     &,         t_lo(imt,km)
      include "fdift.h"

!-----------------------------------------------------------------------
!     limit the indices based on those from the argument list
!-----------------------------------------------------------------------

      istrt   = max(2,is)
      iend    = min(imt-1,ie)
      istrtm1 = istrt - 1
      iendp1  = iend + 1
      jstrt   = js
      jend    = min(je,jmt-1-joff)

!-----------------------------------------------------------------------
!     initialization when calculating jrow 2
!-----------------------------------------------------------------------

      if (joff + js .eq. 2) then
        jstrt = js - 1
        do k=1,km
          do i=istrt-1,iend
            adv_fn(i,k,1) = c0
            anti_fn(i,k,1,n) = c0
            R_plusY(i,k,1,n) = c0
            R_minusY(i,k,1,n) = c0
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     create an inverse land mask
!-----------------------------------------------------------------------

      do j=1,jmw
        do k=1,km
          do i=1,imt
            tmaski(i,k,j) = c1 - tmask(i,k,j)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculate 2*advective (low order scheme) flux across northern,
!     eastern and bottom faces of "T" cells
!-----------------------------------------------------------------------

      jlast = min(jend+1+joff,jmt-1) - joff
      do j=js-1,jlast
        do k=1,km
          do i=istrt,iend
            totadv = adv_vnt(i,k,j)
     &             + adv_vntiso(i,k,j)
            adv_fn(i,k,j) = totadv*
     &                       (t(i,k,j,n,taum1) + t(i,k,j+1,n,taum1))
     &                      + abs(totadv)*
     &                       (t(i,k,j,n,taum1) - t(i,k,j+1,n,taum1))
          enddo
        enddo
      enddo

      jlast = min(jend+1+joff,jmt-1) - joff
      do j=js,jlast
        do k=1,km
          do i=istrtm1,iend
            totadv = adv_vet(i,k,j)
     &             + adv_vetiso(i,k,j)
            adv_fe(i,k,j) = totadv*
     &                       (t(i,k,j,n,taum1) + t(i+1,k,j,n,taum1))
     &                      + abs(totadv)*
     &                       (t(i,k,j,n,taum1) - t(i+1,k,j,n,taum1))
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            totadv = adv_vbt(i,k,j)
     &             + adv_vbtiso(i,k,j)
            adv_fb(i,k,j) = totadv*
     &                      (t(i,k+1,j,n,taum1) + t(i,k,j,n,taum1))
     &                    + abs(totadv)*
     &                      (t(i,k+1,j,n,taum1) - t(i,k,j,n,taum1))
          enddo
        enddo
        do i=istrt,iend
          adv_fb(i,0,j) = adv_vbt(i,0,j)*c2*t(i,1,j,n,taum1)
          adv_fb(i,km,j) = c0
        enddo
      enddo

!-----------------------------------------------------------------------
!     main j loop
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow = (j+1) + joff
        jp2  = min(j+2+joff,jmt) - joff
        jp1  = min(j+1+joff,jmt-1) - joff

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer at center of "T" cells in row j+1
!       - low order solution -
!-----------------------------------------------------------------------

        do k=1,km
          twodt(k) = c2dtts*dtxcel(k)
          do i=istrt,iend
            t_lo(i,k) = (t(i,k,j+1,n,taum1) - twodt(k)
     &        *(ADV_Tx(i,k,jp1) + ADV_Ty(i,k,jp1,jrow,n) +
     &        ADV_Tz(i,k,jp1))*tmask(i,k,j+1))
          enddo
        enddo
        call setbcx (t_lo(1,1), imt, km)

!-----------------------------------------------------------------------
!       next calculate raw antidiffusive fluxes, that is high order
!       scheme flux (leap frog) minus the low order (upstream)
!-----------------------------------------------------------------------

        do k=1,km
          do i=istrtm1,iend
            totadv = adv_vet(i,k,jp1)
     &             + adv_vetiso(i,k,jp1)
            anti_fe(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i+1,k,j+1,n,tau)) - adv_fe(i,k,jp1)
          enddo
          do i=istrt,iend
            totadv = adv_vnt(i,k,jp1)
     &             + adv_vntiso(i,k,jp1)
            anti_fn(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i,k,jp2,n,tau)) - adv_fn(i,k,jp1)
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            totadv = adv_vbt(i,k,jp1)
     &             + adv_vbtiso(i,k,jp1)
            anti_fb(i,k,j+1,n) = totadv*(t(i,k,j+1,n,tau) +
     &                           t(i,k+1,j+1,n,tau)) - adv_fb(i,k,jp1)
     &                           *tmask(i,k,j+1)
          enddo
        enddo
        do i=istrt,iend
          anti_fb(i,0,j+1,n) = adv_vbt(i,0,j+1)*c2*t(i,1,j+1,n,taum1)
          anti_fb(i,km,j+1,n) = c0
        enddo

!-----------------------------------------------------------------------
!       now calculate and apply one-dimensional delimiters to these
!       raw antidiffusive fluxes

!       1) calculate T*, that are all halfway neighbors of T
!       2) calculate ratio R+- of Q+- to P+-, that is maximal/minimal
!          possible change of T if no limit would be active,
!          must be at least 1
!       3) choose correct ratio depending on direction of flow as a
!          delimiter
!       4) apply this delimiter to raw antidiffusive flux
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!       delimit x-direction
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function

!         running mean of two adjacent points

          do i=istrt,iendp1
            Trmax(i) = p5*(t(i-1,k,j+1,n,tau) + t(i,k,j+1,n,tau))
          enddo

!         extremum of low order solution central point and adjacent
!         halfway neighbours; check for land

          do i=istrt,iend
            fxa = tmask(i-1,k,j+1)*Trmax(i) +
     &            tmaski(i-1,k,j+1)*t_lo(i,k)
            fxb = tmask(i+1,k,j+1)*Trmax(i+1) +
     &            tmaski(i+1,k,j+1)*t_lo(i,k)
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
            dcf(i) = cstdxt2r(i,j+1)
            flxlft(i) = anti_fe(i-1,k,j+1,n)
            flxrgt(i) = anti_fe(i,k,j+1,n)
          enddo

!         calculate ratio R

          do i=istrt,iend
            Pplus  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
            Pminus = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
            Qplus  = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
            Rpl(i,k) = min(1.,tmask(i,k,j+1)*Qplus/(Pplus+epsln))
            Rmn(i,k) = min(1.,tmask(i,k,j+1)*Qminus/(Pminus+epsln))
          enddo
          call setbcx (Rpl, imt, km)
          call setbcx (Rmn, imt, km)

!         calculate delimiter using ratio at adjacent points

          do i=istrt,iendp1
            Cpos(i-1) = min(Rpl(i,k),Rmn(i-1,k))
            Cneg(i-1) = min(Rpl(i-1,k),Rmn(i,k))
          enddo

!         finally apply appropriate delimiter to flux

          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = p5*((Cpos(i) + Cneg(i))
     &                               *anti_fe(i,k,j+1,n) +
     &                               (Cpos(i) - Cneg(i))
     &                               *abs(anti_fe(i,k,j+1,n)))
          enddo
        enddo

!-----------------------------------------------------------------------
!       delimit y-direction
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function

          do i=istrt,iend
            fxa = p5*tmask(i,k,j)*(t(i,k,j,n,tau) +
     &            t(i,k,j+1,n,tau)) +
     &            tmaski(i,k,j)*t_lo(i,k)
            fxb = p5*tmask(i,k,jp2)*(t(i,k,j+1,n,tau) +
     &            t(i,k,jp2,n,tau)) +
     &            tmaski(i,k,jp2)*t_lo(i,k)
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
            dcf(i) = cstdyt2r(jrow)
            flxlft(i) = anti_fn(i,k,j,n)
            flxrgt(i) = anti_fn(i,k,j+1,n)
          enddo

!         calculate ratio R, related to a point

          do i=istrt,iend
            Pplus  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
            Pminus = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
            Qplus  = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
            R_plusY(i,k,j+1,n) =
     &        min(1.,tmask(i,k,j+1)*Qplus/(Pplus+epsln))
            R_minusY(i,k,j+1,n) =
     &        min(1.,tmask(i,k,j+1)*Qminus/(Pminus+epsln))
          enddo

!         calculate delimiter using ratio at adjacent points

          do i=istrt,iend
            Cpos(i) = min(R_plusY(i,k,j+1,n),R_minusY(i,k,j,n))
            Cneg(i) = min(R_plusY(i,k,j,n),R_minusY(i,k,j+1,n))
          enddo

!         finally get delimiter c dependent on direction of flux and
!         apply it to raw antidiffusive flux

          do i=istrt,iend
            anti_fn(i,k,j,n) = p5*((Cpos(i) + Cneg(i))
     &                             *anti_fn(i,k,j,n) +
     &                             (Cpos(i) - Cneg(i))
     &                             *abs(anti_fn(i,k,j,n)))
          enddo
        enddo

!-----------------------------------------------------------------------
!       delimit z-direction
!-----------------------------------------------------------------------

        do k=1,km

!         prepare some data for use in statement function

          do i=istrt,iend
            dcf(i) = dzt2r(k)
            flxlft(i) = anti_fb(i,k,j+1,n)
            flxrgt(i) = anti_fb(i,k-1,j+1,n)
            if (k .gt. 1)then
              fxa = p5*tmask(i,k-1,j+1)*
     &              (t(i,k-1,j+1,n,tau) + t(i,k,j+1,n,tau)) +
     &              tmaski(i,k-1,j+1)*t_lo(i,k)
            else
              fxa = t_lo(i,k)
            endif
            if (k .lt. km) then
              fxb = p5*tmask(i,k+1,j+1)*
     &              (t(i,k,j+1,n,tau)+t(i,k+1,j+1,n,tau)) +
     &              tmaski(i,k+1,j+1)*t_lo(i,k)
            else
              fxb = t_lo(i,k)
            endif
            Trmax(i) = max(fxa,fxb,t_lo(i,k))
            Trmin(i) = min(fxa,fxb,t_lo(i,k))
          enddo

!         calculate delimiter using ratio at adjacent points
!         this variable is related to an arc (between two points,
!         the same way as fluxes are defined.)

          do i=istrt,iend
            Pplus  = c2dtts*dcf(i)*(max(c0,flxlft(i))-min(c0,flxrgt(i)))
            Pminus = c2dtts*dcf(i)*(max(c0,flxrgt(i))-min(c0,flxlft(i)))
            Qplus = Trmax(i) - t_lo(i,k)
            Qminus = t_lo(i,k) - Trmin(i)
            Rpl(i,k) = min(1.,tmask(i,k,j+1)*Qplus/(Pplus+epsln))
            Rmn(i,k) = min(1.,tmask(i,k,j+1)*Qminus/(Pminus+epsln))
          enddo

        enddo

        do k=1,kmm1

!         calculate delimiter using ratio at adjacent points
!         this variable is related to an arc (between two points,
!         the same way as fluxes are defined.)

          do i=istrt,iend
            Cneg(i) = min(Rpl(i,k+1),Rmn(i,k))
            Cpos(i) = min(Rpl(i,k),Rmn(i,k+1))
          enddo

!         finally get delimiter c dependent on direction of flux and
!         apply it to raw antidiffusive flux

          do i=istrt,iend
            anti_fb(i,k,j+1,n) = p5*((Cpos(i)+Cneg(i))
     &                               *anti_fb(i,k,j+1,n) +
     &                               (Cpos(i)-Cneg(i))
     &                               *abs(anti_fb(i,k,j+1,n)))
          enddo
        enddo
        do i=istrt,iend
          anti_fb(i,0,j+1,n)  = c0
          anti_fb(i,km,j+1,n) = c0
        enddo

!-----------------------------------------------------------------------
!       complete advective fluxes by adding low order fluxes to
!       delimited antidiffusive fluxes
!-----------------------------------------------------------------------

        do k=1,km
          do i=istrtm1,iend
            anti_fe(i,k,j+1,n) = anti_fe(i,k,j+1,n) + adv_fe(i,k,jp1)
          enddo
          do i=istrt,iend
            anti_fn(i,k,j,n) = (anti_fn(i,k,j,n) + adv_fn(i,k,j))
     &                         *tmask(i,k,j)
            anti_fb(i,k,j+1,n) = (anti_fb(i,k,j+1,n) + adv_fb(i,k,jp1))
     &                           *tmask(i,k,j+1)
          enddo
        enddo

      enddo

!-----------------------------------------------------------------------
!     set 2*corrected advective fluxes across northern, eastern and
!     bottom faces of "T" cells
!-----------------------------------------------------------------------

      do j=js-1,jend
        do k=1,km
          do i=istrt,iend
            adv_fn(i,k,j) = anti_fn(i,k,j,n)
          enddo
        enddo
      enddo

      do j=js,jend
        do k=1,km
          do i=istrtm1,iend
            adv_fe(i,k,j) = anti_fe(i,k,j,n)
          enddo
        enddo

        do k=1,kmm1
          do i=istrt,iend
            adv_fb(i,k,j) = anti_fb(i,k,j,n)
          enddo
        enddo
      enddo

      return
      end
