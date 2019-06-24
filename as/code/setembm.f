! source file: /net/mare/home1/eby/as/updates/setembm.F
      subroutine setembm (is, ie, js, je)

!=======================================================================
!     initialize the energy-moisture balance model
!=======================================================================

      implicit none

      character(120) :: fname, new_file_name, text
      character(3) :: a3

      integer i, ie, ii, iou, is, j, je, jj, js, jz, k, m, n, nsolve
      integer nu, nsum, ib(10), ic(10)

      logical exists, inqvardef
      real dlam, dphi, dtatms, dte, dyz, eccice, grarea, saltmax
      real si, ssh, t1, tair, yz_max, yz_min, wz, calday, tmp
      real zrel, c100, C2K

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "calendar.h"
      include "solve.h"
      include "switch.h"
      include "coord.h"
      include "grdvar.h"
      include "cembm.h"
      include "atm.h"
      include "insolation.h"
      include "ice.h"
      include "evp.h"
      include "riv.h"
      include "tmngr.h"
      include "levind.h"
      include "csbc.h"
      include "scalar.h"
      include "veg.h"
      include "ism.h"
      real rveg(imt,jmt)
      real dmsk(imt,jmt), tmpij(imtm2,jmtm2)
      real tmp_dt(imt,jmt)

      c100      = 100.
      C2K       = 273.15

      cdatm     = 1.e-3
      cpatm     = 1.004e7
      sht       = 8.4e5
      shq       = 1.8e5
      shc       = 8.049e5
      rhoatm    = 1.250e-3
      esatm     = 4.6e-05
      pcfactor  = 0.
      cssh      = 3.8011e-3
      cfc11ccnn = 0.
      cfc11ccns = 0.
      cfc12ccnn = 0.
      cfc12ccns = 0.
      dc14ccnn  = 0.
      dc14ccne  = 0.
      dc14ccns  = 0.

      rhoocn    = 1.035
      esocn     = 5.4e-5
      vlocn     = 2.501e10

      cdice     = 5.5e-3
      rhoice    = 0.910 !to match ice sheet model value
      rhosno    = 0.330
      esice     = 5.347e-5
      slice     = 2.835e10
      flice     = 3.34e9
      condice   = 2.1656e5

      soilmax   = 15.
      eslnd     = 5.347e-5

      nivc      = 1
      dtatms    = 1800.
      ns        = 30

      dalt_v    = 3.3e-3
      dalt_o    = 1.4e-3
      dalt_i    = 1.4e-3

!     ensure pass is between zero and one.
      pass =  min(max((1. - scatter), 0.), 1.)

!     gtoppm is used in converting g carbon cm-2 => ppmv CO2
!     4.138e-7 => 12e-6 g/umol carbon / 29 g/mol air
      gtoppm = 1./(4.138e-7*rhoatm*shc)

!     calculate atmospheric surface area
      atmsa = 0.
      do j=2,jmtm1
        do i=2,imtm1
          atmsa = atmsa + dxt(i)*dyt(j)*cst(j)
        enddo
      enddo

!     sort out forcing years
      if (pyear  .gt. 1.e20) pyear = 1850.
      if (co2_yr .gt. 1.e20) co2_yr = pyear
      if (c14_yr .gt. 1.e20) c14_yr = pyear
      if (orbit_yr .gt. 1.e20) orbit_yr = pyear
      if (crops_yr .gt. 1.e20) crops_yr = pyear
      if (ice_yr .gt. 1.e20) ice_yr = pyear
      if (solar_yr .gt. 1.e20) solar_yr = pyear
      if (sealev_yr .gt. 1.e20) sealev_yr = pyear

      if (mod(timavgint, segtim) .gt. 1.e-6 .and. timavgint .gt. 0.)
     &  then
        t1 = nint(timavgint/segtim)*segtim
        if (t1 .lt. segtim) t1 = segtim
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: "timavgint" does not contain an integral number'
     &,   '              of coupling time steps "segtim".              '
        write (stdout,*) '              (changed "timavgint" from '
     &, timavgint,' days to ', t1,' days to insure this condition)'
        timavgint = t1
      endif
      if (timavgint .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging interval "timavgint" = 0. implies no '
     &,   '             averaging when "time_averages" is enabled      '
      endif
      if (timavgint .gt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the interval "timavgint" exceeds the averaging '
     &,   '             period "timavgper" for option "time_averages"  '
      endif
      if (timavgint .lt. timavgper) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: averaging period "timavgper" exceeds interval  '
     &,   '             "timavgint". Setting timavgper = timavgint     '
        timavgper = timavgint
      endif
      if (timavgper .eq. 0.) then
        write (stdout,'(/,(1x,a))')
     &    '==> Warning: the averaging period "timavgper" is zero. The  '
     &,   '             average will be over only one time step!       '
      endif
      write (stdout,'(/,1x,a,f10.2,a,/,1x,a,f10.2,a)')
     &  '==> Time averages will be written every ', timavgint, ' days, '
     &, '    with an averaging period of         ', timavgper, ' days. '

!-----------------------------------------------------------------------
!     calculate the relative CO2 forcing term
!-----------------------------------------------------------------------

      call co2forc

      write (stdout,*)
      write (stdout,*) 'CO2 ratio (reference = 280 ppmv) =',co2ccn/280.
      write (stdout,*) 'Yields radiative forcing (W/m2) = ',anthro*1.e-3

!-----------------------------------------------------------------------
!     calculate the expansion coefficients for Berger's solution for
!     the year of the initial conditions
!-----------------------------------------------------------------------

      call orbit (orbit_yr, eccen, obliq, mvelp, lambm0)
      write (stdout,*)
      write (stdout,*) 'Initial Orbital Parameters:'
      write (stdout,*) '  Orbital Year:', orbit_yr
      write (stdout,*) '  Eccentricity:', eccen
      write (stdout,*) '  Obliquity:   ', obliq
      write (stdout,*) '  Longitude of Perihelion:', mvelp+180.

!-----------------------------------------------------------------------
!     calculate annual average insolation and Coriolis factor
!-----------------------------------------------------------------------

      radian = 360./(2.*pi)
      do j=1,jmt
        do i=1,imt
!         calculate coriolis parameter
          fcor(i,j) = 2.*omega*sin(ulat(i,j)/radian)
        enddo
      enddo

!-----------------------------------------------------------------------
!     read diffusion
!-----------------------------------------------------------------------

      dn(:,:,:) = 5.e9
      de(:,:,:) = 5.e9
      fname = new_file_name ("diff.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        call openfile (fname, iou)
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        do n=1,nat
          if (n .lt. 1000) write(a3,'(i3)') n
          if (n .lt. 100) write(a3,'(i2)') n
          if (n .lt. 10) write(a3,'(i1)') n
          call getvara ('dn_'//trim(a3), iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          dn(2:imtm1,2:jmtm1,n) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (dn(:,:,n))
          call getvara ('de_'//trim(a3), iou, imtm2*jmtm2, ib, ic
     &,     tmpij, c1, c0)
          de(2:imtm1,2:jmtm1,n) = tmpij(1:imtm2,1:jmtm2)
          call embmbc (de(:,:,n))
        enddo
      endif

!-----------------------------------------------------------------------
!     set solver parameters
!-----------------------------------------------------------------------

      nsolve = 0
      itin(:)  = 500        ! max solver iterations
      epsin(:) = 5.e-7
      epsin(ishum) = 1.e-5
      epsin(isat) = 1.e-3
      nsolve = nsolve + 1
      levelin = 20              ! max coarse grid level
      if (nsolve .ne. 1) then
        write(*,*) '==> Error: more or less than one solver defined.'
        write(*,*) '           Use only one of embm_adi, embm_mgrid,'
     &,   ' embm_slap, embm_essl, embm_sparskit or embm_explicit'
        stop '=>setembm'
      endif

!-----------------------------------------------------------------------
!     check latent heats will sum to zero
!-----------------------------------------------------------------------

      if (slice .ne. vlocn + flice) write (stdout,'(/,a)')
     &   '==> Warning: changing latent heat of fusion to conserve heat'
        flice = slice - vlocn

!-----------------------------------------------------------------------
!     calculate grid terms for the atmospheric solver
!-----------------------------------------------------------------------

      do j=2,jmtm1
        dsgrd(j) = csu(j-1)/(dyu(j-1)*cst(j)*dyt(j))
        dngrd(j) = csu(j)/(dyu(j)*cst(j)*dyt(j))
        asgrd(j) = csu(j-1)/(2.*cst(j)*dyt(j))
        angrd(j) = csu(j)/(2.*cst(j)*dyt(j))
      enddo
      do i=2,imtm1
        dwgrd(i) = 1./(dxu(i-1)*dxt(i))
        degrd(i) = 1./(dxu(i)*dxt(i))
        azgrd(i) = 1./(2.*dxt(i))
      enddo

!-----------------------------------------------------------------------
!     set initial conditions or read a restart
!-----------------------------------------------------------------------

      newcoef(:,:) = .true.

      nats = namix
      dayoyr = 1.
      itt = 0
      irstdy = 0
      msrsdy = 0
      totaltime = 0.
      atbar(:,:) = 0.
      rtbar(:,:) = 0.
      at(:,:,:,:) = 0.
      tair = 13.
      at(:,:,:,isat) = tair
      ssh = cssh*exp(17.67*tair/(tair + 243.5))
      rh(:,:) = rhmax
      at(:,:,:,ishum) = rhmax*ssh
      at(:,:,:,ishum) = 0. !jer
      precip(:,:) = 0.
      awx(:,:) = 0.
      awy(:,:) = 0.
      soilm(:,:,:) = 0.
      surf(:,:) = 0.
      hice(:,:,:) = 0.
      aiceocn(:,:,:) = 0.
      aicelan(:,:,:) = 0.
      ticel(:,:) = 0.
      ticeo(:,:) = 0.
      hsnol(:,:,:) = 0.
      hsnoo(:,:,:) = 0.
      nmbal(:,:,:) = 0.
      dicevol(:)   = 0.
      diceheat(:)  = 0.
      isflxm(:)    = 0.
      isflxh(:)    = 0.
      uice(:,:) = 0.
      vice(:,:) = 0.
      sbc(:,:,isu) = 0.
      sbc(:,:,isv) = 0.
      sbc(:,:,igu) = 0.
      sbc(:,:,igv) = 0.
      sig11n(:,:) = 0.
      sig11e(:,:) = 0.
      sig11s(:,:) = 0.
      sig11w(:,:) = 0.
      sig22n(:,:) = 0.
      sig22e(:,:) = 0.
      sig22s(:,:) = 0.
      sig22w(:,:) = 0.
      sig12n(:,:) = 0.
      sig12e(:,:) = 0.
      sig12s(:,:) = 0.
      sig12w(:,:) = 0.
      bv(:) = 0.
      xv(:) = 0.

!-----------------------------------------------------------------------
!     set land ice data and tracer grid ocean mask
!-----------------------------------------------------------------------
      !set basic land mask, to be redone during ice sheet initialization or
      !via reading of ice-sheet-enabled restart with tmsk included
      tmsk(:,:) = 0.0
      do j=1,jmtm1
        do i=1,imtm1
	  if (kmt(i,j) .gt. 0) tmsk(i,j) = 1.
	enddo
      enddo
      dsealev = sealev

      if (.not. init) then
        fname = new_file_name ("restart_embm.nc")
        inquire (file=trim(fname), exist=exists)
        if (exists) call embm_rest_in (fname, is, ie, js, je)
      endif
      ocnsa = 0.
      do j=2,jmtm1
        do i=2,imtm1
          if (kmt(i,j).gt.0.)
     &	  ocnsa = ocnsa + dxt(i)*dyt(j)*cst(j)
        enddo
      enddo
      ocnsa = ocnsa/100./100./1000./1000. !km^2
      print*, 'Ocean area = ', ocnsa, 'km^2'

!-----------------------------------------------------------------------
!     read average air temperature
!-----------------------------------------------------------------------

      tbar(:,:) = 0.
      fname = new_file_name ("tbar.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Error => ", trim(fname), " does not exist."
        stop 'tbar in setembm.f'
      endif
      ib(:) = 1
      ic(:) = 1
      ic(1) = imtm2
      ic(2) = jmtm2
      call openfile (fname, iou)
      exists = inqvardef('A_slat', iou)
      if (.not. exists) then
        print*, "Error => A_slat does not exist."
        stop 'A_slat in setembm.f'
      endif
      call getvara ('A_slat', iou, imtm2*jmtm2, ib, ic, tmpij, c1, c0)
      tbar(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
      text = "C"
      call getatttext (iou, 'A_slat', 'units', text)
!     convert to model units (C)
      if (trim(text) .eq. "K")
     &  where (tbar(:,:) .lt. 1.e30) tbar(:,:) = tbar(:,:) - C2K
      call embmbc (tbar)

!# if defined 1
      !if ice sheet model defined, and initial conditions, then obtain initial elevation
      !from data, which will then be adjusted during ice sheet initialization.
      !If reading from a restart, the restart will include the adjusted elevation so not
      !necessary to re-read it here.
      !if (init) then
!# endif
!-----------------------------------------------------------------------
!     read land elevations
!-----------------------------------------------------------------------

      elev(:,:) = 0.
      fname = new_file_name ("elev.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        call openfile (fname, iou)
        call getvara ('elev', iou, imtm2*jmtm2, ib, ic, tmpij, c100, c0)
        elev(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (elev)
      endif
!     check for negative elevations
      where (elev(:,:) .lt. 0.) elev(:,:) = 0.

!# if defined 1
      !endif
!# endif

!-----------------------------------------------------------------------
!     initialize running annual averages
!-----------------------------------------------------------------------

      fname = new_file_name ("restart_embm.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (fname, iou)
        exists = inqvardef('rtbar', iou)
      endif
      if (.not. exists .or. init) then
        totaltime = 0.
        atbar(:,:) = 0.
        rtbar(:,:) = tbar(:,:)
      endif

      dmsk(:,:) = 1.
      tmp_dt(:,:) = rtbar(:,:) - tbar(:,:)
      call areaavg (tmp_dt, dmsk, dtbar)

!-----------------------------------------------------------------------
!     set velocity grid ocean mask
!-----------------------------------------------------------------------

      call recalc_umsk

      do j=1,jmt
        do i=1,imt
          if (tmsk(i,j) .gt. 0.0) then !jer
            if (hice(i,j,1) .le. 0.) aiceocn(i,j,1) = 0.
            if (hice(i,j,2) .le. 0.) aiceocn(i,j,2) = 0.
          endif
        enddo
      enddo
      call embmbc (aiceocn(1,1,1))
      call embmbc (aiceocn(1,1,2))
      call embmbc (aicelan(1,1,1))
      call embmbc (aicelan(1,1,2))
!-----------------------------------------------------------------------
!     set the river model
!-----------------------------------------------------------------------

      call rivinit

!-----------------------------------------------------------------------
!     set ocean coalbedo
!-----------------------------------------------------------------------

      do j=1,jmt
        do i=1,imt
          if (kmt(i,j) .gt. 0) then
!           varies from 0.895 at the equator to 0.815 at the pole
            sbc(i,j,iscao) = 0.855 + 0.04*cos(abs(tlat(i,j))*2./radian)
          endif
        enddo
      enddo

!-----------------------------------------------------------------------
!     read vegetation class
!-----------------------------------------------------------------------

      rveg(:,:) = 0.
      fname = new_file_name ("veg_class.nc")
      inquire (file=trim(fname), exist=exists)
      if (.not. exists) then
        print*, "Warning => ", trim(fname), " does not exist."
      else
        ib(:) = 1
        ic(:) = 1
        ic(1) = imtm2
        ic(2) = jmtm2
        call openfile (fname, iou)
        call getvara ('veg', iou, imtm2*jmtm2, ib, ic, tmpij, c1, c0)
        rveg(2:imtm1,2:jmtm1) = tmpij(1:imtm2,1:jmtm2)
        call embmbc (rveg)
      endif
      do j=1,jmt
        do i=1,imt
          if (kmt(i,j) .gt. 0) then

            iveg(i,j) = iice
          elseif (rveg(i,j) .gt. 0.6 .and. rveg(i,j) .lt. 7.4) then
            iveg(i,j) = nint(rveg(i,j))
          else
            print*, "==> Warning: invalid vegetation type, ", rveg(i,j)
     &,             " at point (i,j), ", i, j
            print*, "    Setting point to vegetation type ice"
            iveg(i,j) = iice
          endif
        enddo
      enddo

      do j=1,jmt
        do i=1,imt
	  !reset underlying surface albedo to tundra if default albedo set to ice
          if (kmt(i,j) .le . 0 .and. iveg(i,j) == iice) then
	    iveg(i,j) = 5 !tundra
	  endif
	enddo
      enddo

      call gvsbc

!----------------------------------------------------------------------
!     initialize elastic viscous plastic variables
!-----------------------------------------------------------------------

      dlam = dxu(int(imt/2))/100.
      dphi = dyu(int(jmt/2))/100.
      diff1 = 0.004
      diff1 = diff1*dlam
      diff2 = diff1*dlam**2
      eccice = 2.
      ecc2 = 1./(eccice**2)
      ecc2m = 2.*(1.-ecc2)
      ecc2p = (1.+ecc2)
      zetamin = 4.e11
      eyc = 0.25
      dte = dtatm/float(ndte)
      dtei = 1./dte
      floor = 1.e-11
      do j=2,jmtm1
        do i=2,imtm1
           xyminevp = (min(cst(j)*dxt(i),dyt(j)))**2
        enddo
      enddo

!-----------------------------------------------------------------------
!     check ice velocity calculation
!-----------------------------------------------------------------------

      if (nivts .gt. nint(segtim*daylen/dtatm)) then
        write(*,*) '==> Warning: ice velocities will be calculated'
        write(*,*) '             every coupling time.'
        nivts =  nint(segtim*daylen/dtatm)
      endif
!-----------------------------------------------------------------------
!     zero time average accumulators
!-----------------------------------------------------------------------

      call ta_embm_tavg (is, ie, js, je, 0)

!-----------------------------------------------------------------------
!     zero integrated time average accumulators
!-----------------------------------------------------------------------

      call ta_embm_tsi (is, ie, js, je, 0)
      return
      end

      subroutine recalc_umsk

      implicit none

      include "size.h"        !imt,jmt
      include "param.h"       !imtm1,jmtm1
      include "atm.h"         !tmsk,umsk
      include "levind.h"      !kmt

      real imsk(imt,jmt)

      integer i,j

      character(120) g_st

      !reset velocity grid ocean mask
      call embmbc(tmsk)

      imsk(:,:) = 0.
      where (tmsk .gt. 0.) imsk = 1.  !anywhere there is partial/full ocean, set imsk = 1

      umsk(:,:) = 0.

      !sweep from SE, northwards.  If tmsk point, or leading tmsk points, is 0, umsk =0
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (imsk(i,j), imsk(i+1,j), imsk(i,j+1)
     &, 		   imsk(i+1,j+1))
        enddo
      enddo
      call embmbc (umsk)
      !remove isolated bays

      !sweep from SE, northwards.  If umsk = 1
      do j=2,jmtm1
        do i=2,imtm1
          imsk(i,j) = max (umsk(i,j), umsk(i-1,j), umsk(i,j-1)
     &, 		 umsk(i-1,j-1))
        enddo
      enddo
      call embmbc (imsk)
      do j=2,jmtm1
        do i=2,imtm1
          umsk(i,j) = min (imsk(i,j), imsk(i+1,j), imsk(i,j+1)
     &, 		   imsk(i+1,j+1))
        enddo
      enddo

      call embmbc(umsk)

      !where (imsk .eq. 0.) tmsk = 0.
      !call embmbc(tmsk)

      return
      end subroutine
