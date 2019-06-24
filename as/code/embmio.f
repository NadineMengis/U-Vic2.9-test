! source file: /net/mare/home1/eby/as/updates/embmio.F
       subroutine embmout (is, ie, js, je)

!=======================================================================
!     output routine for energy-moisture balance model

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!=======================================================================

      use subgrid

      implicit none

      character(120) :: fname, new_file_name
      character(32) :: nstamp
      character(3) :: a3

      integer i, ie, is, id_xt, id_yt, iou, j, je, js, n, L
      integer ndx, ntrec, it(10), ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      logical exists, inqvardef

      real avgper, ca, time, wt, wtp1
      real c100, c500, C2K, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "csbc.h"
      include "atm.h"
      include "solve.h"
      include "coord.h"
      include "grdvar.h"
      include "levind.h"
      include "ice.h"
      include "evp.h"
      include "ism.h"
      include "mtlm.h"
      include "cembm.h"
      include "iounit.h"
      include "scalar.h"
      include "switch.h"
      include "tmngr.h"
      include "riv.h"
      include "cregin.h"

      real tmp_at(imt,jmt,nat)
      real sm(imt,jmt), st(imt,jmt), hs(imt,jmt), ro(imt,jmt)
      real rntatsl, rntatil
      real tmp_dt(imt,jmt)

      real p_alb(is:ie,js:je), a_alb(is:ie,js:je), s_alb(is:ie,js:je)
      real sat(is:ie,js:je)

      real dmsk(imt,jmt)

      real elevi,totar
      integer sn,en

      c100 = 100.
      c500 = 500.
      C2K = 273.15

!-----------------------------------------------------------------------
!     write atmospheric diagnostics
!-----------------------------------------------------------------------

      if (tsits .and. ntatia .ne. 0) then

        call ta_embm_tsi (is, ie, js, je, 2)

        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)

        if (iotsi .eq. stdout .or. iotsi .lt. 0) then
          write (*,'(1x, a3, i7, 1x, a32, 3(a,1pe13.6))')
     &      'ts=',itt, nstamp, ' iterations =', tai_maxit
     &,     ' TAbar=', tai_sat, ' QAbar=', tai_shum
        endif

        if (iotsi .ne. stdout) then
          avgper = dtatm*ntatia/daylen
          time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
          rntatil = 0.
          if (ntatil .gt. 0) rntatil = 1./float(ntatil)
          tai_hsno = tai_hsno + tai_LYING_SNOW*1000.*rntatil/rhosno
          call def_tsi
          call def_tsi_embm (fname)
!         convert emissions from g cm-2 to kg s-1
          tmp = tai_co2emit*atmsa*1.e-3
          tai_catm = 0.
          call embm_tsi_out (fname, avgper, time, nstamp, tai_sat
     &,                      tai_shum, tai_precip, tai_evap, tai_ohice
     &,                      tai_oaice, tai_hsno, tai_lhice, tai_laice
     &,                      tai_co2ccn, tmp, tai_dc14ccn
     &,                      tai_cfc11ccn, tai_cfc12ccn, tai_maxit
     &,                      tai_nsat, tai_ssat, tai_nshum, tai_sshum
     &,                      tai_nprecip, tai_sprecip, tai_nevap
     &,                      tai_sevap, tai_nohice, tai_sohice
     &,                      tai_noaice, tai_soaice, tai_nhsno
     &,                      tai_shsno, tai_nlhice, tai_slhice
     &,                      tai_nlaice, tai_slaice, tai_lsat, tai_osat
     &,                      tai_lprecip, tai_oprecip, tai_levap
     &,                      tai_oevap, tai_solins, tai_outlwr
     &,                      tai_dnswr, tai_inswr, tai_netrad, tai_palb
     &,                      tai_aalb, tai_salb, tai_lsalb, tai_osalb
     &,                      tai_sst, tai_sss, tai_ssdic, tai_ssc14
     &,                      tai_ssalk, tai_sso2, tai_sspo4, tai_ssno3
     &,                      tai_sscfc11, tai_sscfc12, tai_sulph
     &,                      tai_volc, tai_agg,tai_catm
     &,			     isvol, nisvol, sisvol, nissmb,sissmb
     &,   		     nisiar,sisiar,nismext,sismext
     &,			     nisaccum,sisaccum,nismelt,sismelt
     &,			     zclim_out,sealev_out
     &,			     ntrec)
        endif

        call ta_embm_tsi (is, ie, js, je, 0)

      endif

      if (timavgts .and. ntatsa .ne. 0) then

!-----------------------------------------------------------------------
!       write atmospheric time averaged data
!-----------------------------------------------------------------------

!       calculate average values

        call ta_embm_tavg (is, ie, js, je, 2)

!       write time averaged data

        rntatsl = 0.
        if (ntatsl .gt. 0) rntatsl = 1./float(ntatsl)
        call unloadland (POINTS, TA_M, imt, jmt, land_map, sm)
        call unloadland (POINTS, TA_TSTAR_GB, imt, jmt, land_map, st)
        call unloadland (POINTS, TA_LYING_SNOW, imt, jmt, land_map, hs)
        call unloadland (POINTS, TA_SURF_ROFF, imt, jmt, land_map, ro)
        do j=js,je
          do i=is,ie
	    elevi=0.
	    sn=strtiE(i,j)
	    en=endiE(i,j)
	    if (sn .le. en) then
	      totar=sum(iariE(sn:en))
	      elevi=sum(eleviE(sn:en)*iariE(sn:en))/totar
	    endif
            sat(i,j) =  ta_at(i,j,isat)

     &	    -(elev(i,j)*rlapse*(1.-fraci(i,j)))
     &      -elevi*rlapse*fraci(i,j)

            if (ta_solins(i,j) .gt. 0) then
!             calculate incoming swr reaching the surface
              tmp = ta_solins(i,j) - ta_arswr(i,j) - ta_absin(i,j)
              s_alb(i,j) = 1. - ta_dnswr(i,j)/tmp
              a_alb(i,j) = ta_arswr(i,j)/ta_solins(i,j)
!             calculate outgoing swr leaving the atmosphere
              tmp = tmp - ta_dnswr(i,j) - ta_absout(i,j)
              p_alb(i,j) = (ta_arswr(i,j) + tmp)/ta_solins(i,j)
            else
              s_alb(i,j) = -1.
              a_alb(i,j) = -1.
              p_alb(i,j) = -1.
            endif
            if (land_map(i,j) .eq. 0) then
              sm(i,j) = ta_soilm(i,j)
              st(i,j) = ta_surf(i,j)
            else
!             convert from kg m-2 to cm (converted back later)
              sm(i,j) = sm(i,j)*rntatsl*0.1
!             convert from K to C (converted back later)
              st(i,j) = st(i,j)*rntatsl - C2K
!             convert from kg/m2/s to g/cm2/s (converted back later)
              ta_runoff(i,j) = ro(i,j)*rntatsl*0.1
      !       convert from kg/m2 to cm (converted back later)
              ta_hsno(i,j) = hs(i,j)*rntatsl*0.1/rhosno
            endif
            if (ta_hsno(i,j) .lt. 0) ta_hsno(i,j) = 0.
            if (ta_hice(i,j) .lt. 0) ta_hice(i,j) = 0.
            if (ta_aice(i,j) .lt. 0) ta_aice(i,j) = 0.
            if (ta_aice(i,j) .gt. 1) ta_aice(i,j) = 1.
          enddo
        enddo

        avgper = dtatm*ntatsa/daylen
        time = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
        nyear = year0 + accel_yr0 + (relyr - accel_yr0)*accel
        call mkstmp (nstamp, nyear, nmonth, nday, nhour, nmin, nsec)
        call def_tavg
        call def_tavg_embm (fname)
        call embm_tavg_out (fname, is, ie, js, je, imt, jmt, nat, ncat
     &,   xt, yt, xu, yu, dxt, dyt, dxu, dyu, avgper, time, nstamp
     &,   mapat, ta_at, sat, ta_precip, ta_evap, ta_disch, ta_vflux
     &,   ta_outlwr, ta_uplwr, ta_upsens, ta_dnswr, ta_upltnt
     &,   ta_solins, p_alb, a_alb, s_alb, elev, ta_psno, ta_ws
     &,   ta_runoff, ta_wx, ta_wy
     &,   ta_awx, ta_awy, ta_rtbar, ta_apress
     &,   sm, st
     &,   ta_tice, ta_hice, ta_aice, ta_hsno
     &,   ta_uice, ta_vice, ta_xint, ta_yint
     &,   tlat, tlon, ulat, ulon, tgarea, ugarea
     &,   tmsk, mskhr, nriv, ntrec)

        write (*,'(a,i5,a,a,a,a)') '=> Atm time means #'
     &,   ntrec, ' written to ',trim(fname),' on ', stamp

!       zero time average accumulators

        call ta_embm_tavg (is, ie, js, je, 0)

      endif

        ocnsa = 0.
        do j=2,jmtm1
          do i=2,imtm1
	    if (kmt(i,j) .gt. 0) ocnsa = ocnsa + dxt(i)*dyt(j)*cst(j)
          enddo
        enddo
        ocnsa=ocnsa/10000.
!-----------------------------------------------------------------------
!     do things that need to be done only once a year
!-----------------------------------------------------------------------

      if (eoyear) then

!       make new annual and running average
        call embmbc (atbar)
        if (totaltime .ge. yrlen*daylen - dtatm*0.5) then
          wt = 0.05
          do j=js,je
            do i=is,ie
              rtbar(i,j) = (1.-wt)*rtbar(i,j) + wt*atbar(i,j)/totaltime
              atbar(i,j) = 0.0
            enddo
          enddo
        endif
        totaltime = 0.0

        call calc_awind (is, ie, js, je)

        dmsk(:,:) = 1.
        tmp_dt(:,:) = rtbar(:,:) - tbar(:,:)
        call areaavg (tmp_dt, dmsk, dtbar)

      endif

!-----------------------------------------------------------------------
!     write atmospheric restart
!-----------------------------------------------------------------------

      if (restrt) then
        if (restts) then
          call def_rest (0)
          call def_rest_embm (0, fname)
          call embm_rest_out (fname, is, ie, js, je)
        endif
        if (eorun) then
          call def_rest (1)
          call def_rest_embm (1, fname)
          call embm_rest_out (fname, is, ie, js, je)
        endif
      endif

      return
      end

      subroutine ta_embm_tavg (is, ie, js, je, iflag)

!=======================================================================
!     atmospheric data time averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer is, ie, js, je, iflag, i, j, k, n, ndx

      real rntatsa

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "atm.h"
      include "riv.h"
      include "ice.h"
      include "csbc.h"

      real dmsk(imt,jmt)
      real evap(imt,jmt)

      evap(:,:) = evapo(:,:)*tmsk(:,:) + evapl(:,:)*(1.-tmsk(:,:))

!-----------------------------------------------------------------------
!     time averaged data
!-----------------------------------------------------------------------

      if (iflag .eq. 0.) then

!       zero
        ntatsa = 0
        ta_at(:,:,:) = 0.
        ta_solins(:,:) = 0.
        ta_arswr(:,:) = 0.
        ta_dnswr(:,:) = 0.
        ta_absin(:,:) = 0.
        ta_absout(:,:) = 0.
        ta_uplwr(:,:) = 0.
        ta_upsens(:,:) = 0.
        ta_upltnt(:,:) = 0.
        ta_outlwr(:,:) = 0.
        ta_precip(:,:) = 0.
        ta_psno(:,:) = 0.
        ta_evap(:,:) = 0.
        ta_disch(:,:) = 0.
        ta_vflux(:,:) = 0.
        ta_ws(:,:) = 0.
        ta_runoff(:,:) = 0.
        ta_wx(:,:,:) = 0.
        ta_wy(:,:,:) = 0.
        ta_awx(:,:) = 0.
        ta_awy(:,:) = 0.
        ta_rtbar(:,:) = 0.
        ta_apress(:,:) = 0.
        ta_soilm(:,:) = 0.
        ta_surf(:,:) = 0.
        ta_tice(:,:) = 0.
        ta_hice(:,:) = 0.
        ta_aice(:,:) = 0.
        ta_hsno(:,:) = 0.
        ta_uice(:,:) = 0.
        ta_vice(:,:) = 0.
        ta_pice(:,:) = 0.
        ta_xint(:,:) = 0.
        ta_yint(:,:) = 0.
        ta_psum(:) = 0.
      elseif (iflag .eq. 1) then

!       accumulate

        ntatsa = ntatsa + 1
        do n=1,nat
          ta_at(:,:,n) = ta_at(:,:,n) + at(:,:,2,n)
        enddo
        ta_solins(:,:) = ta_solins(:,:) + solins(:,:)
        ta_arswr(:,:) = ta_arswr(:,:) + solins(:,:)
     &                  *(1. - sbc(:,:,iaca))
        ta_dnswr(:,:) = ta_dnswr(:,:) + dnswro(:,:)*tmsk(:,:)
     &                  	      + dnswrl(:,:)*(1.-tmsk(:,:))
        ta_absin(:,:) = ta_absin(:,:) + solins(:,:)*sbc(:,:,iaca)
     &                  *scatter
        ta_absout(:,:) = ta_absout(:,:) + (solins(:,:)*sbc(:,:,iaca)
     &                   *pass - dnswro(:,:))*scatter
        ta_uplwr(:,:) = ta_uplwr(:,:) + uplwro(:,:)*tmsk(:,:)
     &                  	      + uplwrl(:,:)*(1.-tmsk(:,:))
        ta_upsens(:,:) = ta_upsens(:,:) + upsenso(:,:)*tmsk(:,:)
     &                  	        + upsensl(:,:)*(1.-tmsk(:,:))
        ta_upltnt(:,:) = ta_upltnt(:,:) + upltnto(:,:)*tmsk(:,:)
     &                  	        + upltntl(:,:)*(1.-tmsk(:,:))
        ta_outlwr(:,:) = ta_outlwr(:,:) + outlwr(:,:)
        ta_precip(:,:) = ta_precip(:,:) + precip(:,:)
        ta_psno(:,:) = ta_psno(:,:) + psnoo(:,:)*tmsk(:,:)
     &                  	    + psnol(:,:)*(1.-tmsk(:,:))
        ta_evap(:,:) = ta_evap(:,:) + evap(:,:)
        ta_disch(:,:) = ta_disch(:,:) + disch(:,:)
        ta_vflux(:,:) = ta_vflux(:,:) + vflux(:,:)
        ta_ws(:,:) = ta_ws(:,:) + sbc(:,:,iws)
        ta_runoff(:,:) = ta_runoff(:,:) + runoff(:,:)
        ta_wx(:,:,ishum) = ta_wx(:,:,ishum) + sbc(:,:,iwxq)
        ta_wy(:,:,ishum) = ta_wy(:,:,ishum) + sbc(:,:,iwyq)
        ta_wx(:,:,isat) = ta_wx(:,:,isat) + sbc(:,:,iwxt)
        ta_wy(:,:,isat) = ta_wy(:,:,isat) + sbc(:,:,iwyt)
        ta_awx(:,:) = ta_awx(:,:) + awx(:,:)
        ta_awy(:,:) = ta_awy(:,:) + awy(:,:)
        ta_rtbar(:,:) = ta_rtbar(:,:) + rtbar(:,:)
        ta_apress(:,:) = ta_apress(:,:) + apress(:,:)
        ta_soilm(:,:) = ta_soilm(:,:) + soilm(:,:,2)
        ta_surf(:,:) = ta_surf(:,:) + surf(:,:)
        ta_tice(:,:) = ta_tice(:,:) + ticeo(:,:)*tmsk(:,:)
     &                  	    + ticel(:,:)*(1.-tmsk(:,:))
        ta_hice(:,:) = ta_hice(:,:) + hice(:,:,2)
        ta_aice(:,:) = ta_aice(:,:) + aiceocn(:,:,2)*tmsk(:,:)
    ! &                  	    + aicelan(:,:,2)*(1.-tmsk(:,:))
        ta_hsno(:,:) = ta_hsno(:,:) + hsnoo(:,:,2)*tmsk(:,:)
     &                  	    + hsnol(:,:,2)*(1.-tmsk(:,:))
        ta_uice(:,:) = ta_uice(:,:) + uice(:,:)
        ta_vice(:,:) = ta_vice(:,:) + vice(:,:)
        ta_pice(:,:) = ta_pice(:,:) + pice(:,:)
        ta_xint(:,:) = ta_xint(:,:) + xint(:,:)
        ta_yint(:,:) = ta_yint(:,:) + yint(:,:)

        ta_psum(:) = ta_psum(:) + psum(:)

      elseif (iflag .eq. 2 .and. ntatsa .ne. 0) then

!       average
        rntatsa = 1./float(ntatsa)
        ta_at(:,:,:) = ta_at(:,:,:)*rntatsa
        ta_solins(:,:) = ta_solins(:,:)*rntatsa
        ta_arswr(:,:) = ta_arswr(:,:)*rntatsa
        ta_dnswr(:,:) = ta_dnswr(:,:)*rntatsa
        ta_absin(:,:) = ta_absin(:,:)*rntatsa
        ta_absout(:,:) = ta_absout(:,:)*rntatsa
        ta_uplwr(:,:) = ta_uplwr(:,:)*rntatsa
        ta_upsens(:,:) = ta_upsens(:,:)*rntatsa
        ta_upltnt(:,:) = ta_upltnt(:,:)*rntatsa
        ta_outlwr(:,:) = ta_outlwr(:,:)*rntatsa
        ta_precip(:,:) = ta_precip(:,:)*rntatsa
        ta_psno(:,:) = ta_psno(:,:)*rntatsa
        ta_evap(:,:) = ta_evap(:,:)*rntatsa
        ta_disch(:,:) = ta_disch(:,:)*rntatsa
        ta_vflux(:,:) = ta_vflux(:,:)*rntatsa
        ta_ws(:,:) = ta_ws(:,:)*rntatsa
        ta_runoff(:,:) = ta_runoff(:,:)*rntatsa
        ta_wx(:,:,:) = ta_wx(:,:,:)*rntatsa
        ta_wy(:,:,:) = ta_wy(:,:,:)*rntatsa
        ta_awx(:,:) = ta_awx(:,:)*rntatsa
        ta_awy(:,:) = ta_awy(:,:)*rntatsa
        ta_rtbar(:,:) = ta_rtbar(:,:)*rntatsa
        ta_apress(:,:) = ta_apress(:,:)*rntatsa
        ta_soilm(:,:) = ta_soilm(:,:)*rntatsa
        ta_surf(:,:) = ta_surf(:,:)*rntatsa
        ta_hsno(:,:) = ta_hsno(:,:)*rntatsa
        ta_tice(:,:) = ta_tice(:,:)*rntatsa
        ta_hice(:,:) = ta_hice(:,:)*rntatsa
        ta_aice(:,:) = ta_aice(:,:)*rntatsa
        ta_uice(:,:) = ta_uice(:,:)*rntatsa
        ta_vice(:,:) = ta_vice(:,:)*rntatsa
        ta_pice(:,:) = ta_pice(:,:)*rntatsa
        ta_xint(:,:) = ta_xint(:,:)*rntatsa
        ta_yint(:,:) = ta_yint(:,:)*rntatsa
        ta_psum(:) = ta_psum(:)*rntatsa

      endif

      return
      end

      subroutine ta_embm_tsi (is, ie, js, je, iflag)

!=======================================================================
!     atmospheric data time integral averaging

!     input:
!       is, ie, js, je = starting and ending indicies for i and j
!       iflag = flag (0 = zero, 1 = accumulate, 2 = write)
!=======================================================================

      implicit none

      integer i, is, ie, j, js, je, iflag, maxit, n

      real rntatia, sum, tmp

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "grdvar.h"
      include "cembm.h"
      include "atm.h"
      include "ice.h"
      include "solve.h"

      real sat(imt,jmt), dmsk(imt,jmt)
      real ts_salb(imt,jmt), ts_aalb(imt,jmt), ts_palb(imt,jmt)
      real ts_netrad(imt,jmt), ts_inswr(imt,jmt)
      real evap(imt,jmt)

!-----------------------------------------------------------------------
!     time averaged integrated data
!-----------------------------------------------------------------------

      evap(:,:) = evapo(:,:)*tmsk(:,:) + evapl(:,:)*(1.-tmsk(:,:))

      if (iflag .eq. 0.) then

!       zero
        ntatia = 0
        tai_maxit = 0.
        tai_sat  = 0.
        tai_shum = 0.
        tai_precip = 0.
        tai_evap = 0.
        tai_co2ccn = 0.
        tai_co2emit = 0.
        tai_dc14ccn = 0.
        tai_cfc11ccn = 0.
        tai_cfc12ccn = 0.
        tai_osat = 0.
        tai_oprecip = 0.
        tai_oevap = 0.
        tai_ohice = 0.
        tai_oaice = 0.
        tai_hsno = 0.
        tai_lsat = 0.
        tai_lprecip = 0.
        tai_levap = 0.
        tai_lhice = 0.
        tai_laice = 0.
        tai_nsat = 0.
        tai_nshum = 0.
        tai_nprecip = 0.
        tai_nevap = 0.
        tai_nohice = 0.
        tai_noaice = 0.
        tai_nhsno = 0.
        tai_nlhice = 0.
        tai_nlaice = 0.
        tai_ssat = 0.
        tai_sshum = 0.
        tai_sprecip = 0.
        tai_sevap = 0.
        tai_sohice = 0.
        tai_soaice = 0.
        tai_shsno = 0.
        tai_slhice = 0.
        tai_slaice = 0.
	tai_solins = 0.
	tai_dnswr = 0.
	tai_outlwr = 0.
	tai_inswr = 0.
	tai_netrad = 0.
	tai_palb = 0.
	tai_aalb = 0.
	tai_salb = 0.
	tai_lsalb = 0.
	tai_osalb = 0.
        tai_sst = 0.
        tai_sss = 0.
        tai_ssdic = 0.
        tai_ssc14 = 0.
        tai_ssalk = 0.
        tai_sso2 = 0.
        tai_sspo4 = 0.
        tai_ssno3 = 0.
        tai_sscfc11 = 0.
        tai_sscfc12 = 0.
        tai_sulph = 0.
        tai_volc = 0.
        tai_agg = 0.

      elseif (iflag .eq. 1) then

!       accumulate
        ntatia = ntatia + 1
        maxit = 0
        do n=1,nat
          if (itout(n) .gt. maxit) maxit = itout(n)
        enddo
        tai_maxit = tai_maxit + maxit
        sat(:,:) = at(:,:,2,isat) - elev(:,:)*rlapse
!-----------------------------------------------------------------------
!       set data mask for global
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.

        call areaavg (sat, dmsk, tmp)
        tai_sat  = tai_sat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_shum = tai_shum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_precip = tai_precip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_evap = tai_evap + tmp
        tai_co2ccn = tai_co2ccn + co2ccn
        tai_co2emit = tai_co2emit + co2emit

!-----------------------------------------------------------------------
!       set data mask for global ocean
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_osat  = tai_osat + tmp
        call areaavg (precip, dmsk, tmp)
        tai_oprecip  = tai_oprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_oevap  = tai_oevap + tmp
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_ohice = tai_ohice + tmp
        call areatot (aiceocn(1,1,2), dmsk, tmp)
        tai_oaice = tai_oaice + tmp
        call areatot (hsnoo(1,1,2), dmsk, tmp)
        tai_hsno = tai_hsno + tmp

!-----------------------------------------------------------------------
!       set data mask for global land
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_lsat  = tai_lsat + tmp
        call areaavg (precip, dmsk, tmp)
        tai_lprecip  = tai_lprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_levap  = tai_levap + tmp
        call areatot (hsnoo(1,1,2), dmsk, tmp)
        tai_hsno = tai_hsno + tmp

!-----------------------------------------------------------------------
!       set data mask for northern hemisphere
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0.) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_nsat  = tai_nsat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_nshum  = tai_nshum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_nprecip  = tai_nprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_nevap  = tai_nevap + tmp

!-----------------------------------------------------------------------
!       set data mask for northern hemisphere ocean
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0. .or. tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_nohice = tai_nohice + tmp
        call areatot (aiceocn(1,1,2), dmsk, tmp)
        tai_noaice = tai_noaice + tmp
        call areatot (hsnoo(1,1,2), dmsk, tmp)
        tai_nhsno = tai_nhsno + tmp

!-----------------------------------------------------------------------
!       set data mask for northern hemisphere land
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .lt. 0. .or. tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.

        call areatot (hsnoo(1,1,2), dmsk, tmp)
        tai_nhsno = tai_nhsno + tmp

!-----------------------------------------------------------------------
!       set data mask for southern hemisphere
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .ge. 0.) dmsk(:,:) = 0.

        call areaavg (sat, dmsk, tmp)
        tai_ssat  = tai_ssat + tmp
        call areaavg (at(1,1,2,ishum), dmsk, tmp)
        tai_sshum  = tai_sshum + tmp
        call areaavg (precip, dmsk, tmp)
        tai_sprecip  = tai_sprecip + tmp
        call areaavg (evap, dmsk, tmp)
        tai_sevap  = tai_sevap + tmp

!-----------------------------------------------------------------------
!       set data mask for southern hemisphere ocean
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .ge. 0. .or. tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areatot (hice(1,1,2), dmsk, tmp)
        tai_sohice = tai_sohice + tmp
        call areatot (aiceocn(1,1,2), dmsk, tmp)
        tai_soaice = tai_soaice + tmp
        call areatot (hsnoo(1,1,2), dmsk, tmp)
        tai_shsno = tai_shsno + tmp

!-----------------------------------------------------------------------
!       set data mask for southern hemisphere land
!-----------------------------------------------------------------------
        dmsk(:,:) = 1.
        where (tlat(:,:) .ge. 0. .or. tmsk(:,:) .ge. 0.5) dmsk(:,:) = 0.

        call areatot (hsnoo(1,1,2), dmsk, tmp)
        tai_shsno = tai_shsno + tmp
        do j=1,jmt
          do i=1,imt
            if (solins(i,j) .gt. 0) then
    !         calculate incoming swr reaching the surface
              tmp = solins(i,j) - solins(i,j)*(1. - sbc(i,j,iaca))
     &                 - solins(i,j)*sbc(i,j,iaca)*scatter
              ts_salb(i,j) = 1. - dnswro(i,j)/tmp
              ts_aalb(i,j) = solins(i,j)*(1.-sbc(i,j,iaca))/solins(i,j)
    !         calculate outgoing swr leaving the atmosphere
              tmp = tmp - dnswro(i,j) - (solins(i,j)*sbc(i,j,iaca)
     &                   *pass - dnswro(i,j))*scatter
              ts_palb(i,j) = (solins(i,j)*(1. - sbc(i,j,iaca)) + tmp)
     &                   /solins(i,j)
              ts_inswr(i,j) = solins(i,j)*(1. - ts_palb(i,j))
              ts_netrad(i,j) = ts_inswr(i,j) - outlwr(i,j)
            else
              ts_salb(i,j) = -1.
              ts_aalb(i,j) = -1.
              ts_palb(i,j) = -1.
	      ts_inswr(i,j) = 0.
	      ts_netrad(i,j) = 0. - outlwr(i,j)
            endif
          enddo
        enddo

        dmsk(:,:) = 1.
        call areaavg (solins, dmsk, tmp)
	tai_solins = tai_solins + tmp
        call areaavg (outlwr, dmsk, tmp)
	tai_outlwr = tai_outlwr + tmp
	call areaavg (dnswro, dmsk, tmp)
	tai_dnswr = tai_dnswr + tmp
        call areaavg (ts_inswr, dmsk, tmp)
	tai_inswr = tai_inswr + tmp
        call areaavg (ts_netrad, dmsk, tmp)
	tai_netrad = tai_netrad + tmp
        dmsk(:,:) = 1.
        where (ts_palb(:,:) .lt. 0.) dmsk(:,:) = 0.
        call areaavg (ts_palb, dmsk, tmp)
	tai_palb = tai_palb + tmp
        dmsk(:,:) = 1.
        where (ts_aalb(:,:) .lt. 0.) dmsk(:,:) = 0.
	call areaavg (ts_aalb, dmsk, tmp)
	tai_aalb = tai_aalb + tmp
        dmsk(:,:) = 1.
        where (ts_salb(:,:) .lt. 0.) dmsk(:,:) = 0.
        call areaavg (ts_salb, dmsk, tmp)
	tai_salb = tai_salb + tmp
        dmsk(:,:) = 1.
        where (tmsk(:,:) .ge. 0.5 .or. ts_salb(:,:) .lt. 0.)
     &    dmsk(:,:) = 0.
	call areaavg (ts_salb, dmsk, tmp)
	tai_lsalb = tai_lsalb + tmp
        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5 .or. ts_salb(:,:) .lt. 0.)
     &    dmsk(:,:) = 0.
	call areaavg (ts_salb, dmsk, tmp)
	tai_osalb = tai_osalb + tmp
        dmsk(:,:) = 1.
        where (tmsk(:,:) .lt. 0.5) dmsk(:,:) = 0.
        call areaavg (sbc(1,1,isst), dmsk, tmp)
	tai_sst = tai_sst + tmp
        call areaavg (sbc(1,1,isss), dmsk, tmp)
	tai_sss = tai_sss + tmp
        call areaavg (sbc(1,1,issdic), dmsk, tmp)
	tai_ssdic = tai_ssdic + tmp
        call areaavg (sbc(1,1,issalk), dmsk, tmp)
	tai_ssalk = tai_ssalk + tmp
        call areaavg (sbc(1,1,isso2), dmsk, tmp)
	tai_sso2 = tai_sso2 + tmp
        call areaavg (sbc(1,1,isspo4), dmsk, tmp)
	tai_sspo4 = tai_sspo4 + tmp
        call areaavg (sbc(1,1,issno3), dmsk, tmp)
	tai_ssno3 = tai_ssno3 + tmp

      elseif (iflag .eq. 2 .and. ntatia .ne. 0) then

!       average
        rntatia = 0.
        if (ntatia .gt. 0.) rntatia = 1./float(ntatia)
        tai_maxit = tai_maxit*rntatia
        tai_sat  = tai_sat*rntatia
        tai_shum = tai_shum*rntatia
        tai_precip = tai_precip*rntatia
        tai_evap = tai_evap*rntatia
        tai_co2ccn = tai_co2ccn*rntatia
        tai_co2emit = tai_co2emit*rntatia
        tai_dc14ccn = tai_dc14ccn*rntatia
        tai_cfc11ccn = tai_cfc11ccn*rntatia
        tai_cfc12ccn = tai_cfc12ccn*rntatia
        tai_osat = tai_osat*rntatia
        tai_oprecip = tai_oprecip*rntatia
        tai_oevap = tai_oevap*rntatia
        tai_ohice = tai_ohice*rntatia
        tai_oaice = tai_oaice*rntatia
        tai_hsno = tai_hsno*rntatia
        tai_lsat = tai_lsat*rntatia
        tai_lprecip = tai_lprecip*rntatia
        tai_levap = tai_levap*rntatia
        tai_lhice = tai_lhice*rntatia
        tai_laice = tai_laice*rntatia
        tai_nsat = tai_nsat*rntatia
        tai_nshum = tai_nshum*rntatia
        tai_nprecip = tai_nprecip*rntatia
        tai_nevap = tai_nevap*rntatia
        tai_nohice = tai_nohice*rntatia
        tai_noaice = tai_noaice*rntatia
        tai_nhsno = tai_nhsno*rntatia
        tai_nlhice = tai_nlhice*rntatia
        tai_nlaice = tai_nlaice*rntatia
        tai_ssat = tai_ssat*rntatia
        tai_sshum = tai_sshum*rntatia
        tai_sprecip = tai_sprecip*rntatia
        tai_sevap = tai_sevap*rntatia
        tai_sohice = tai_sohice*rntatia
        tai_soaice = tai_soaice*rntatia
        tai_shsno = tai_shsno*rntatia
        tai_slhice = tai_slhice*rntatia
        tai_slaice = tai_slaice*rntatia
	tai_solins = tai_solins*rntatia
	tai_outlwr = tai_outlwr*rntatia
	tai_dnswr = tai_dnswr*rntatia
	tai_inswr = tai_inswr*rntatia
	tai_netrad = tai_netrad*rntatia
	tai_palb = tai_palb*rntatia
	tai_aalb = tai_aalb*rntatia
	tai_salb = tai_salb*rntatia
	tai_lsalb = tai_lsalb*rntatia
	tai_osalb = tai_osalb*rntatia
        tai_sst = tai_sst*rntatia
        tai_sss = tai_sss*rntatia
        tai_ssdic = tai_ssdic*rntatia
        tai_ssc14 = tai_ssc14*rntatia
        tai_ssalk = tai_ssalk*rntatia
        tai_sso2 = tai_sso2*rntatia
        tai_sspo4 = tai_sspo4*rntatia
        tai_ssno3 = tai_ssno3*rntatia
        tai_sscfc11 = tai_sscfc11*rntatia
        tai_sscfc12 = tai_sscfc12*rntatia
        tai_sulph = tai_sulph*rntatia
        tai_volc = tai_volc*rntatia
        tai_agg = tai_agg*rntatia

      endif

      return
      end