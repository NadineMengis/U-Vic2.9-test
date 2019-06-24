! source file: /net/mare/home1/eby/as/updates/tracer.F
      subroutine tracer (joff, js, je, is, ie)

!=======================================================================
!     compute tracers at "tau+1" for rows js through je in the MW.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name

      integer istrt, iend, i, k, j, ip, kr, jq, n, jp, jrow, iou, js
      integer je, limit, joff, is, ie, kmx, m, kb, idiag, index
      integer it(10), iu(10), ib(10), ic(10), nfnpzd, mfnpzd, mxfnpzd
      integer id_time, id_xt, id_yt, id_zt

      logical inqvardef, exists

      real rctheta, declin, gl, impo, expo, npp, time
      real remi, excr, graz, morp, morpt, morz, temp, swr, dayfrac
      real phin, dz, prca, dprca, nud, bct, tap, fo2, so2, ai, hi, hs
      real npp_D, graz_D, morp_D, no3flag, deni, nfix
      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, daylen, yrlen
      real zmax, cont, drho, drhom1, wt, ahbi_cstr, ahbi_csu_dyur
      real gamma, rrstd, fy, fyz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
      include "emode.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "timeavgs.h"
      include "tmngr.h"
      include "vmixc.h"
      include "diaga.h"
      include "ice.h"
      include "atm.h"
      include "npzd.h"

      parameter (istrt=2, iend=imt-1)

      real twodt(km)
      real snpzd(ntnpzd), tnpzd(ntnpzd)
      real src(imt,km,jsmw:jemw,nsrc)
      include "isopyc.h"
      include "fdift.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude indices based on those from the argument list
!     Note: this is currently bypassed. istrt and iend are set as
!           parameters to optimize performance
!-----------------------------------------------------------------------

!      istrt = max(2,is)
!      iend  = min(imt-1,ie)

!-----------------------------------------------------------------------
!     build coefficients to minimize advection and diffusion computation
!-----------------------------------------------------------------------

      limit = min(je+1+joff,jmt) - joff
      do j=js,limit
        jrow = j + joff
        do i=istrt-1,iend
          cstdxtr(i,j)    = cstr(jrow)*dxtr(i)
          cstdxt2r(i,j)   = cstr(jrow)*dxtr(i)*p5
          cstdxur(i,j)    = cstr(jrow)*dxur(i)
          ah_cstdxur(i,j) = diff_cet*cstr(jrow)*dxur(i)
        enddo
      enddo

!-----------------------------------------------------------------------
!     calculation of biological interactions
!-----------------------------------------------------------------------

      declin = sin((mod(relyr,1.) - 0.22)*2.*pi)*0.4   ! declination

      do k=1,km
        twodt(k) = c2dtts*dtxcel(k)
        nbio(k) = twodt(k)/dtnpzd
        dtbio(k) = twodt(k)/nbio(k)
        rdtts(k) = 1./twodt(k)
        rnbio(k) = 1./nbio(k)
      enddo
      tap = 2.*alpha*par

      do j=js,je
        jrow = j + joff

        do i=is,ie

          if (kmt(i,jrow) .gt. 0) then

            ai = aiceocn(i,jrow,2)
            hi = hice(i,jrow,2)
            hs = hsnoo(i,jrow,2)
!           calculate day fraction and incoming solar
!           angle of incidence = lat - declin, refraction index = 1.33
            rctheta = max(-1.5, min(1.5, tlat(i,jrow)/radian - declin))
            rctheta = kw/sqrt(1. - (1. - cos(rctheta)**2.)/1.33**2.)
            dayfrac = min( 1., -tan(tlat(i,jrow)/radian)*tan(declin))
            dayfrac = max(1e-12, acos(max(-1., dayfrac))/pi)
            swr = dnswro(i,jrow)*1e-3*(1. + ai*(exp(-ki*(hi + hs))-1.))
            expo = 0.0
            impo = 0.0
            phin = 0.0 ! integrated phytoplankton
            prca = 0.0 ! integrated production of calcite

            kmx = min(kmt(i,jrow), kpzd)
            do k=1,kmx

!-----------------------------------------------------------------------
!             limit tracers to positive values
!-----------------------------------------------------------------------

              tnpzd(1) = max(t(i,k,j,ipo4,taum1), trcmin)
              tnpzd(2) = max(t(i,k,j,iphyt,taum1), trcmin)
              tnpzd(3) = max(t(i,k,j,izoop,taum1), trcmin)
              tnpzd(4) = max(t(i,k,j,idetr,taum1), trcmin)
              tnpzd(5) = max(t(i,k,j,ino3,taum1), trcmin)
              tnpzd(6) = max(t(i,k,j,idiaz,taum1), trcmin)

              swr = swr*exp(-kc*phin)
              phin = phin + tnpzd(2)*dzt(k)
              gl = tap*swr*exp(ztt(k)*rctheta)
              impo = expo*dztr(k)
              bct = bbio**(cbio*t(i,k,j,itemp,taum1))
!             decrease remineralisation rate in oxygen minimum zone
              nud = nud0*(0.65+0.35*tanh(t(i,k,j,io2,taum1)*1000.-6.))

!-----------------------------------------------------------------------
!             call the npzd model
!-----------------------------------------------------------------------

              call npzd_src (tnpzd, nbio(k), dtbio(k), gl, bct, impo
     &,                      dzt(k), dayfrac, wd(k), rkwz(k), nud
     &,                      snpzd, expo, graz, morp, morz
     &                      )

              snpzd(1:4) = snpzd(1:4)*rdtts(k)
              snpzd(5:6) = snpzd(5:6)*rdtts(k)
              expo = expo*rnbio(k)
!-----------------------------------------------------------------------
!             calculate detritus at the bottom and remineralize
!-----------------------------------------------------------------------

              if (k .eq. kmt(i,jrow)) then
                if (addflxo .and. eots)
     &            sbc(i,jrow,irorg) = sbc(i,jrow,irorg)
     &                              + expo*redctn*twodt(1)*dzt(k)
                snpzd(1) = snpzd(1) + redptn*expo
                snpzd(5) = snpzd(5) + expo
              endif

!-----------------------------------------------------------------------
!             set source terms
!-----------------------------------------------------------------------

              src(i,k,j,ispo4) = snpzd(1)
              src(i,k,j,isphyt) = snpzd(2)
              src(i,k,j,iszoop) = snpzd(3)
              src(i,k,j,isdetr) = snpzd(4)
              src(i,k,j,isno3) = snpzd(5)
              src(i,k,j,isdiaz) = snpzd(6)
!             production of calcite
              dprca = (morp+morz+graz*(1.-gamma1))*capr*redctn*rnbio(k)
              prca = prca + dprca*dzt(k)
              src(i,k,j,isdic) = (snpzd(1)*redctp - dprca)
              src(i,k,j,isalk) = (-snpzd(1)*redntp*1.e-3 - 2.*dprca)
!             calculate total export to get total import for next layer
              expo = expo*dzt(k)

            enddo

            kmx = kmt(i,jrow)
            do k=1,kmx
!             limit oxygen consumption below concentrations of
!             5umol/kg as recommended in OCMIP
              fo2 = 0.5*tanh(t(i,k,j,io2,taum1)*1000. - 5.)
!             sink of oxygen
              so2 = src(i,k,j,ispo4)*redotp
              src(i,k,j,iso2) = -so2*(0.5 + fo2)
!             add denitrification as source term for NO3
              no3flag = 0.5+sign(0.5,t(i,k,j,ino3,taum1)-trcmin)
!             800 = 0.8*1000 = (elec/mol O2)/(elec/mol NO3)*(mmol/mol)
              deni = 800.*no3flag*so2*(0.5 - fo2)
              src(i,k,j,isno3) = src(i,k,j,isno3) - deni
            enddo

!-----------------------------------------------------------------------
!           remineralize calcite
!-----------------------------------------------------------------------
            kmx = kmt(i,jrow)
            do k=1,kmx-1
              src(i,k,j,isdic) = src(i,k,j,isdic) + prca*rcak(k)
              src(i,k,j,isalk) = src(i,k,j,isalk) + 2.*prca*rcak(k)
            enddo
            if (addflxo .and. eots)
     &        sbc(i,jrow,ircal) = sbc(i,jrow,ircal) + (prca*rcab(kmx)
     &                          - prca*rcak(kmx))*twodt(1)*dzt(kmx)
            src(i,kmx,j,isdic) = src(i,kmx,j,isdic) + prca*rcab(kmx)
            src(i,kmx,j,isalk) = src(i,kmx,j,isalk) + 2.*prca*rcab(kmx)

          endif

        enddo

      enddo

!-----------------------------------------------------------------------
!     set source terms from sediment model
!-----------------------------------------------------------------------
      do j=js,je
        jrow = j + joff
        do i=is,ie
          if (kmt(i,jrow) .gt. 0) then
            k = kmt(i,jrow)
            src(i,k,j,isdic) = src(i,k,j,isdic) + sbc(i,j,ibdicfx)
            src(i,k,j,isalk) = src(i,k,j,isalk) + sbc(i,j,ibalkfx)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     solve for one tracer at a time
!     n = 1 => temperature
!     n = 2 => salinity
!     n > 2 => other tracers (if applicable)
!-----------------------------------------------------------------------

      do n=1,nt

!-----------------------------------------------------------------------
!       calculate advective tracer flux
!-----------------------------------------------------------------------

        call adv_flux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       calculate diffusive flux across eastern and northern faces
!       of "T" cells due to various parameterizations for diffusion.
!-----------------------------------------------------------------------

!       diffusive flux on eastern face of "T" cells

        do j=js,je
          do k=1,km
            do i=istrt-1,iend
              diff_fe(i,k,j) =
     &                         ah_cstdxur(i,j)*
     &                         (t(i+1,k,j,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo

!       diffusive flux on northern face of "T" cells
!       (background for isopycnal mixing)

        do j=js-1,je
          jrow = j + joff
          do k=1,km
            do i=istrt,iend
              diff_fn(i,k,j) =
     &                         diff_cnt*
     &           csu_dyur(jrow)*(t(i,k,j+1,n,taum1) - t(i,k,j,n,taum1))
            enddo
          enddo
        enddo
!-----------------------------------------------------------------------
!       calculate diffusive flux across bottom face of "T" cells
!-----------------------------------------------------------------------

        do j=js,je
          do k=1,km-1
            do i=istrt,iend
              diff_fb(i,k,j) = diff_cbt(i,k,j)*dzwr(k)*
     &                         (t(i,k,j,n,taum1) - t(i,k+1,j,n,taum1))
            enddo
          enddo
        enddo

!-----------------------------------------------------------------------
!       compute isopycnal diffusive flux through east, north,
!       and bottom faces of T cells.
!-----------------------------------------------------------------------

        call isoflux (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!       set surface and bottom vert b.c. on "T" cells for diffusion
!       and advection. for isopycnal diffusion, set adiabatic boundary
!       conditions.
!       note: the b.c. at adv_fb(i,k=bottom,j) is set by the above code.
!             However, it is not set when k=km so it is set below.
!             adv_fb(i,km,j) is always zero (to within roundoff).
!-----------------------------------------------------------------------

        do j=js,je
          jrow   = j + joff
          do i=istrt,iend
            kb              = kmt(i,jrow)
            diff_fb(i,0,j)  = stf(i,j,n)
            diff_fb(i,kb,j) = btf(i,j,n)
            adv_fb(i,0,j)   = adv_vbt(i,0,j)*(t(i,1,j,n,tau) +
     &                                        t(i,1,j,n,tau))
            adv_fb(i,km,j)  = adv_vbt(i,km,j)*t(i,km,j,n,tau)
          enddo
        enddo

!-----------------------------------------------------------------------
!       set source term for "T" cells
!-----------------------------------------------------------------------

        source(:,:,:) = c0

        if (itrc(n) .ne. 0) then
          do j=js,je
            do k=1,km
              do i=istrt,iend
                source(i,k,j) = src(i,k,j,itrc(n))
              enddo
            enddo
          enddo
        endif

!-----------------------------------------------------------------------
!       solve for "tau+1" tracer using statement functions to represent
!       each component of the calculation
!-----------------------------------------------------------------------

!       1st: solve using all components which are treated explicitly

        do j=js,je
          jrow   = j + joff
          do k=1,km
            twodt(k) = c2dtts*dtxcel(k)
            do i=istrt,iend
              t(i,k,j,n,taup1) = t(i,k,j,n,taum1) + twodt(k)*(
     &          DIFF_Tx(i,k,j) + DIFF_Ty(i,k,j,jrow,n) + DIFF_Tz(i,k,j)
     &          - ADV_Tx(i,k,j) -  ADV_Ty(i,k,j,jrow,n) -  ADV_Tz(i,k,j)
     &          + source(i,k,j)
     &          )*tmask(i,k,j)
            enddo
          enddo
        enddo

!       2nd: add in portion of vertical diffusion handled implicitly

        call ivdift (joff, js, je, istrt, iend, n, twodt)

        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo

!-----------------------------------------------------------------------
!       construct diagnostics associated with tracer "n"
!-----------------------------------------------------------------------

        call diagt1 (joff, js, je, istrt, iend, n, twodt)

!-----------------------------------------------------------------------
!       end of tracer component "n" loop
!-----------------------------------------------------------------------

      enddo

!-----------------------------------------------------------------------
!     explicit convection: adjust column if gravitationally unstable
!-----------------------------------------------------------------------

      call convct2 (t(1,1,1,1,taup1), joff, js, je, is, ie, kmt)
      do j=js,je
        do n=1,nt
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo
      if (timavgperts .and. eots) then
        if (joff .eq. 0) nta_conv = nta_conv + 1
        do j=js,je
          jrow = j + joff
          do i=istrt,iend
            ta_totalk(i,jrow) = ta_totalk(i,jrow) + totalk(i,j)
            ta_vdepth(i,jrow) = ta_vdepth(i,jrow) + vdepth(i,j)
            ta_pe(i,jrow) = ta_pe(i,jrow) + pe(i,j)
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     construct diagnostics after convection
!-----------------------------------------------------------------------

      idiag = 10
      call diagt2 (joff, js, je, istrt, iend, idiag)

!-----------------------------------------------------------------------
!     filter tracers at high latitudes
!-----------------------------------------------------------------------

      if (istrt .eq. 2 .and. iend .eq. imt-1) then
        call filt (joff, js, je)
      else
        write (stdout,'(a)')
     &  'Error: filtering requires is=2 and ie=imt-1 in tracer'
        stop '=>tracer'
      endif
      do n=1,nt
        do j=js,je
          call setbcx (t(1,1,j,n,taup1), imt, km)
        enddo
      enddo

!-----------------------------------------------------------------------
!     construct diagnostics after filtering (for total dT/dt)
!-----------------------------------------------------------------------

      idiag = 1
      call diagt2 (joff, js, je, istrt, iend, idiag)

!-----------------------------------------------------------------------
!     if needed, construct the Atmos S.B.C.(surface boundary conditions)
!     averaged over this segment
!     eg: SST and possibly SSS
!-----------------------------------------------------------------------

      call asbct (joff, js, je, istrt, iend, isst, itemp)
      call asbct (joff, js, je, istrt, iend, isss, isalt)
      call asbct (joff, js, je, istrt, iend, issdic, idic)
      call asbct (joff, js, je, istrt, iend, issalk, ialk)
      call asbct (joff, js, je, istrt, iend, isso2, io2)
      call asbct (joff, js, je, istrt, iend, isspo4, ipo4)
      call asbct (joff, js, je, istrt, iend, issphyt, iphyt)
      call asbct (joff, js, je, istrt, iend, isszoop, izoop)
      call asbct (joff, js, je, istrt, iend, issdetr, idetr)
      call asbct (joff, js, je, istrt, iend, issno3, ino3)
      call asbct (joff, js, je, istrt, iend, issdiaz, idiaz)
!-----------------------------------------------------------------------
!     accumulate bottom tracer values for the sediment model
!-----------------------------------------------------------------------

      if (addflxo .and. eots) then
        if (joff .eq. 0) atsed = atsed + twodt(1)
        do j=js,je
          jrow = j + joff
          do i=is,ie
            if (kmt(i,jrow) .gt. 0) then
              k = kmt(i,jrow)
              sbc(i,jrow,ibtemp) = sbc(i,jrow,ibtemp)
     &                           + t(i,k,j,itemp,taup1)*twodt(1)
              sbc(i,jrow,ibsalt) = sbc(i,jrow,ibsalt)
     &                           + t(i,k,j,isalt,taup1)*twodt(1)
              sbc(i,jrow,ibdic) = sbc(i,jrow,ibdic)
     &                          + t(i,k,j,idic,taup1)*twodt(1)
              sbc(i,jrow,ibalk) = sbc(i,jrow,ibalk)
     &                          + t(i,k,j,ialk,taup1)*twodt(1)
              sbc(i,jrow,ibo2) = sbc(i,jrow,ibo2)
     &                         + t(i,k,j,io2,taup1)*twodt(1)
            endif
          enddo
        enddo
      endif

      return
      end

      subroutine diagt1 (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     construct diagnostics associated with tracer component "n"

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = (1,2) = (u,v) velocity component
!       twodt = (2*dtts,dtts) on (leapfrog,mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer i, k, j, ip, kr, jq, n, jp, jrow, js, je, joff, is, ie
      integer mask, m

      real t_i, t_j, dz_t2r, dz_tr, dz_wtr, dx_t2r, dx_tr, dy_t2r
      real dy_tr, adv_tx, adv_ty, adv_tz, adv_txiso, adv_tyiso
      real adv_tziso, diff_tx, diff_ty, diff_tz, dtdx, dtdy, dtdz
      real r2dt, cosdyt, fx, darea, boxar, rtwodt, sumdx, delx
      real sumdxr, dxdy, dxdydz

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "accel.h"
      include "coord.h"
      include "cregin.h"
      include "csbc.h"
      include "ctavg.h"
      include "diag.h"
      include "diaga.h"
      include "emode.h"
      include "grdvar.h"
      include "hmixc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "vmixc.h"

      real temp1(imt,km), temp2(imt,km), temp3(imt,km)
      real twodt(km)

      include "isopyc.h"
      include "fdift.h"

!-----------------------------------------------------------------------
!     diagnostic: integrate |d(tracer)/dt|  and tracer variance on "tau"
!                 globally
!-----------------------------------------------------------------------

      if (tsiperts .and. eots) then
        do j=js,je
          jrow = j + joff
          r2dt    = c1/c2dtts
          cosdyt  = cst(jrow)*dyt(jrow)
          do k=1,km
            fx = r2dt/dtxcel(k)
            do i=is,ie
              darea      = dzt(k)*dxt(i)*cosdyt*tmask(i,k,j)
              temp3(i,k) = t(i,k,j,n,tau)*darea
              temp1(i,k) = t(i,k,j,n,tau)**2*darea
              temp2(i,k) = abs(t(i,k,j,n,taup1)-t(i,k,j,n,taum1))*
     &                     darea*fx
            enddo
            do i=is,ie
              tbar(k,n,jrow)   = tbar(k,n,jrow) + temp3(i,k)
              travar(k,n,jrow) = travar(k,n,jrow) + temp1(i,k)
              dtabs(k,n,jrow)  = dtabs(k,n,jrow) + temp2(i,k)
            enddo
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     diagnostic: accumulate tracers for averages under horizontal
!                 regions (use units of meters, rather than cm)
!-----------------------------------------------------------------------

      if (tavgts .and. eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            mask = mskhr(i,jrow)
            if (mask .ne. 0) then
              boxar = cst(jrow)*dxt(i)*dyt(jrow)*tmask(i,1,j)*0.0001
              sumbf(mask,n) = sumbf(mask,n) + stf(i,j,n)*boxar
              do k=1,km
                sumbk(mask,k,n) = sumbk(mask,k,n) + t(i,k,j,n,tau)
     &                             *boxar*dzt(k)*tmask(i,k,j)*0.01
             enddo
            endif
          enddo
        enddo
      endif

!-----------------------------------------------------------------------
!     diagnostic: compute the northward transport components of
!                 each tracer
!-----------------------------------------------------------------------

      if (gyrets .and. eots)  call gyre (joff, js, je, is, ie, n)

!-----------------------------------------------------------------------
!     diagnostic: integrate r.h.s. terms in the tracer equations
!                 over specified regional volumes.
!-----------------------------------------------------------------------

      if (trmbts .and. eots)  call ttb1 (joff, js, je, is, ie, n)

      return
      end

      subroutine diagt2 (joff, js, je, is, ie, idiag)

!-----------------------------------------------------------------------
!     construct d(tracer)/dt diagnostics

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       idiag = 1  => total tracer change
!       idiag = 10 => change of tracer due to filtering(also convection)
!-----------------------------------------------------------------------

      implicit none

      integer idiag, j, js, je, k, i, joff, iocv, jrow, is, ie

      real rdt, reltim, period

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "diaga.h"
      include "iounit.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "timeavgs.h"

!-----------------------------------------------------------------------
!     diagnostic: integrate d/dt(tracer) over specified regional volumes
!                  after convection and filtering
!-----------------------------------------------------------------------

      if (trmbts .and. eots) call ttb2 (joff, js, je, is, ie, idiag)

      return
      end

      subroutine asbct (joff, js, je, is, ie, isbc, itr)

!-----------------------------------------------------------------------
!     construct the Atmos S.B.C. (surface boundary conditions)

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       isbc  = index for sbc
!       itr   = index for tracer
!-----------------------------------------------------------------------

      implicit none

      integer isbc, itr, j, js, je, jrow, joff, i, is, ie

      real rts

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "levind.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

!     initialize the Atmos S.B.C. at the start of each ocean segment
!     (do not alter values in land)

      if (isbc .le. 0 .or. itr .le. 0) return

      if (eots .and. osegs) then
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0) sbc(i,jrow,isbc) = c0
          enddo
        enddo
      endif

!     accumulate surface tracers for the Atmos S.B.C. every time step

      if (eots) then
        do j=js,je
          jrow = j + joff
          do i=is,ie
            sbc(i,jrow,isbc) = sbc(i,jrow,isbc)+t(i,1,j,itr,taup1)
          enddo
        enddo
      endif

!     average the surface tracers for the Atmos S.B.C. at the end of
!     each ocean segment. (do not alter values in land)

      if (eots .and. osege) then
        rts = c1/ntspos
        do j=js,je
          jrow   = j + joff
          do i=is,ie
            if (kmt(i,jrow) .ne. 0)
     &        sbc(i,jrow,isbc) = rts*sbc(i,jrow,isbc)
          enddo
        enddo
      endif

      return
      end

      subroutine ivdift (joff, js, je, is, ie, n, twodt)

!-----------------------------------------------------------------------
!     solve vertical diffusion of tracers implicitly

!     input:
!       joff  = offset relating "j" in the MW to latitude "jrow"
!       js    = starting row in the MW
!       je    = ending row in the MW
!       is    = starting longitude index in the MW
!       ie    = ending longitude index in the MW
!       n     = tracer component
!       twodt = (2*dtts, dtts) on (leapfrog, mixing) time steps
!-----------------------------------------------------------------------

      implicit none

      integer j, js, je, k, i, is, ie, n, joff

      real rc2dt

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "levind.h"
      include "mw.h"
      include "switch.h"
      include "vmixc.h"

      real twodt(km)

!     store terms to compute implicit vertical mixing on
!     diagnostic time steps

      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            do i=is,ie
              zzi(i,k,j) = t(i,k,j,n,taup1)
            enddo
          enddo
        enddo
      endif

      call invtri (t(1,1,1,n,taup1), stf(1,1,n), btf(1,1,n)
     &, diff_cbt(1,1,jsmw), twodt, kmt, tmask(1,1,1), is, ie
     &, joff, js, je)

!     compute residual implicit vertical mixing

      if (trmbts .and. eots) then
        do j=js,je
          do k=1,km
            rc2dt = c1/twodt(k)
            do i=is,ie
              zzi(i,k,j) = rc2dt*(t(i,k,j,n,taup1) - zzi(i,k,j))
            enddo
          enddo
        enddo
      endif

      return
      end

      subroutine swflux0 (joff, js, je, is, ie, source)

      return
      end
