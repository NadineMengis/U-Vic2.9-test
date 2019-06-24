! source file: /net/mare/home1/eby/as/2.9.old3/source/embm/embm_tavg.F
      subroutine embm_tavg_def (fname, imt, jmt, nat, ncat, xt, yt
     &,                         calendar, expnam, runstamp, mapat)

!=======================================================================
!     definition routine for atmospheric time averages

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     calendar     = calendar
!     expnam       = experiment name
!     runstamp     = run stamp
!     mapat        = tracer map
!=======================================================================

      implicit none

      integer iou, j, n, imt, jmt, nat, ncat, igs, ige, ig, jgs, jge
      integer jg, it(10), iu(10), id_time, id_xt, id_xu, id_yt, id_yu
      integer id_cat, id_cat_e, id_xt_e, id_xu_e, id_yt_e, id_yu_e

      character(*) :: fname, calendar, expnam, runstamp
      character(3) :: a3
      character(10) :: mapat(nat)

      real xt(imt), yt(jmt)
      real c0, c1, c100, c500, c1e3, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     start definitions
!-----------------------------------------------------------------------
      call redef (iou)

!-----------------------------------------------------------------------
!     write global attributes
!-----------------------------------------------------------------------
      call putatttext (iou, 'global', 'Conventions', 'CF-1.0')
      call putatttext (iou, 'global', 'experiment_name', expnam)
      call putatttext (iou, 'global', 'run_stamp', runstamp)

!-----------------------------------------------------------------------
!     define dimensions
!-----------------------------------------------------------------------
      call defdim ('time', iou, 0, id_time)
      call defdim ('xt', iou, ig, id_xt)
      call defdim ('yt', iou, jg, id_yt)
      call defdim ('xu', iou, ig, id_xu)
      call defdim ('yu', iou, jg, id_yu)
      call defdim ('xt_edges', iou, ig+1, id_xt_e)
      call defdim ('yt_edges', iou, jg+1, id_yt_e)
      call defdim ('xu_edges', iou, ig+1, id_xu_e)
      call defdim ('yu_edges', iou, jg+1, id_yu_e)

!-----------------------------------------------------------------------
!       define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time', 'time', 'year since 0000-01-01 00:00:00')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_year', iou, 1, it, c0, c0, ' ', 'F'
     &, 'T_year', ' ','year')
      call defvar ('T_month', iou, 1, it, c0, c0, ' ', 'F'
     &, 'T_month', ' ','month')
      call defvar ('T_day', iou, 1, it, c0, c0, ' ', 'F'
     &, 'T_day', ' ','day')
      call defvar ('T_hour', iou, 1, it, c0, c0, ' ', 'F'
     &, 'T_hour', ' ','hour')
      call defvar ('T_minute', iou, 1, it, c0, c0, ' ', 'F'
     &, 'T_minute', ' ','minute')
      call defvar ('T_second', iou, 1, it, c0, c0, ' ', 'F'
     &, 'T_second', ' ','second')
      call defvar ('T_avgper', iou, 1, it, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','day')

!-----------------------------------------------------------------------
!     define 1d data (x, y or z)
!-----------------------------------------------------------------------
      it(1) = id_xt
      call defvar ('xt', iou, 1, it, c0, c0, 'X', 'D'
     &, 't grid longitude', 'grid_longitude', 'degree')
      call defvar ('G_dxt', iou, 1, it, c0, c0, ' ', 'D'
     &, 't grid width', ' ', 'm')
      it(1) = id_yt
      call defvar ('yt', iou, 1, it, c0, c0, 'Y', 'D'
     &, 't grid latitude', 'grid_latitude', 'degree')
      call defvar ('G_dyt', iou, 1, it, c0, c0, ' ', 'D'
     &, 't grid height', ' ', 'm')
      it(1) = id_xu
      call defvar ('xu', iou, 1, it, c0, c0, 'X', 'D'
     &, 'u grid longitude', 'grid_longitude', 'degree')
      call defvar ('G_dxu', iou, 1, it, c0, c0, ' ', 'D'
     &, 'u grid width', ' ', 'm')
      it(1) = id_yu
      call defvar ('yu', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'u grid latitude', 'grid_latitude', 'degree')
      call defvar ('G_dyu', iou, 1, it, c0, c0, ' ', 'D'
     &, 'u grid height', ' ', 'm')
      it(1) = id_xt_e
      call defvar ('xt_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude of t grid edges', ' ', 'degrees')
      it(1) = id_yt_e
      call defvar ('yt_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude of t grid edges', ' ', 'degrees')
      it(1) = id_xu_e
      call defvar ('xu_edges', iou, 1, it, c0, c0, 'X', 'D'
     &, 'longitude of u grid edges', ' ', 'degrees')
      it(1) = id_yu_e
      call defvar ('yu_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &, 'latitude of u grid edges', ' ', 'degrees')

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      call defvar ('L_elev', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'land elevation and ocean depth', 'surface_altitude', 'm')
      call defvar ('G_mskt', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'ocean mask', ' ' ,'1')
      call defvar ('G_mskhr', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'horizontal region mask', ' ' ,'1')
      call defvar ('G_mskriv', iou, 2, it, c0, c1e3, ' ', 'I'
     &, 'river mask', ' ' ,'1')
      call defvar ('G_latT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid latitude', 'latitude', 'degree_north')
      call defvar ('G_lonT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid longitude', 'longitude', 'degree_east')
      call defvar ('G_latU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid latitude', 'latitude', 'degree_north')
      call defvar ('G_lonU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid longitude', 'longitude', 'degree_east')
      call defvar ('G_areaT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid area', ' ', 'm2')
      call defvar ('G_areaU', iou, 2, iu, -c1e6, c1e6, ' ', 'F'
     &, 'velocity grid area', ' ', 'm2')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      iu(1) = id_xu
      it(2) = id_yt
      iu(2) = id_yu
      it(3) = id_time
      iu(3) = id_time
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          call defvar ('A_slat', iou, 3, it, -c100, c500, ' ', 'F'
     &,     'sea level atmospheric temperature'
     &,     'air_temperature', 'K')
        elseif (trim(mapat(n)) .eq. 'shum') then
          call defvar ('A_shum', iou, 3, it, -c100, c100, ' ', 'F'
     &,     'atmospheric surface specific humidity'
     &,     'specific_humidity', '1')
        elseif (trim(mapat(n)) .eq. 'co2') then
          call defvar ('A_co2', iou, 3, it, c0, c1e6, ' ', 'F'
     &,     'atmospheric co2', 'atmospheric_co2', 'ppm')
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          call defvar ('A_tracer'//trim(a3), iou ,3, it, -c1e6, c1e6
     &,     'F', 'tracer '//trim(a3)
     &,     ' ', 'tracer_'//trim(a3), 'unknown')
        endif
      enddo
      call defvar ('A_sat', iou, 3, it, -c100, c500, ' ', 'F'
     &, 'atmospheric surface temperature', 'air_temperature', 'K')
      call defvar ('F_precip', iou, 3, it, c0, c1, ' ', 'F'
     &, 'precipitation (includes snow in water equivalent)'
     &, 'precipitation_flux', 'kg m-2 s-1')
      call defvar ('F_snow', iou, 3, it, c0, c1, ' ', 'F'
     &, 'precipitation as snow'
     &, 'snowfall_flux', 'kg m-2 s-1')
      call defvar ('F_evap', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'upward evaporation plus sublimation'
     &, 'water_evaporation_flux', 'kg m-2 s-1')
      call defvar ('F_rivdis', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'river discharge', 'river_discharge_flux', 'kg m-2 s-1')
      call defvar ('F_virtual', iou, 3, it, -c1, c1, ' ', 'F'
     &, 'normalized virtual flux', 'virtual_flux', 'm-2 s-1')
      call defvar ('F_outlwr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'outgoing longwave from top of atmosphere'
     &, 'toa_outgoing_longwave_flux', 'W m-2')
      call defvar ('F_uplwr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward longwave'
     &, 'surface_net_upward_longwave_flux', 'W m-2')
      call defvar ('F_upsens', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward sensible heat'
     &, 'surface_upward_sensible_heat_flux', 'W m-2')
      call defvar ('F_dnswr', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'net surface downward short wave (absorbed)'
     &, 'surface_net_downward_shortwave_flux', 'W/m^2')
      call defvar ('F_upltnt', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'surface upward latent heat'
     &, 'surface_upward_latent_heat_flux', 'W m-2')
      call defvar ('F_solins', iou, 3, it, -c1e3, c1e3, ' ', 'F'
     &, 'incoming solar insolation'
     &, 'toa_incoming_shortwave_flux', 'W m-2')
      call defvar ('A_albplt', iou, 3, it, c0, c1, ' ', 'F'
     &, 'planetary albedo', 'planetary_albedo', '1')
      call defvar ('A_albatm', iou, 3, it, c0, c1, ' ', 'F'
     &, 'atmospheric albedo', ' ', '1')
      call defvar ('A_albsur', iou, 3, it, c0, c1, ' ', 'F'
     &, 'surface albedo', 'surface_albedo', '1')
      call defvar ('A_awindX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward wind  anomaly', 'eastward_wind_anomaly', 'm s-1')
      call defvar ('A_awindY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward wind anomaly', 'northward_wind_anomaly', 'm s-1')
      call defvar ('A_avgslat', iou, 3, iu, -c100, c500, ' ', 'F'
     &, 'running average sea level atmospheric temperature'
     &, 'air_temperature', 'K')
      call defvar ('A_apress', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'anomalous surface pressure', 'pressure', 'Pa')
      call defvar ('L_soilmois', iou, 3, it, c0, c1e3, ' ', 'F'
     &, 'soil moisture', 'soil_moisture_content', 'kg m-2')
      call defvar ('L_tempsur', iou, 3, it, -c100, c500, ' ', 'F'
     &, 'land surface temperature', 'surface_temperature', 'K')
      call defvar ('O_icetemp', iou, 3, it, -c100, c500, ' ', 'F'
     &, 'surface ice temperature', 'surface_temperature'
     &, 'K')
      call defvar ('O_icethk', iou, 3, it, c0, c1e6, ' ', 'F'
     &, 'ice thickness', 'sea_ice_thickness', 'm')
      call defvar ('O_icefra', iou, 3, it, c0, c100, ' ', 'F'
     &, 'ice area fraction (includes land ice area fraction)'
     &, 'sea_ice_area_fraction', '1')
      call defvar ('L_snothk', iou, 3, it, -0.01, c1e6, ' ', 'F'
     &, 'surface snow thickness', 'surface_snow_thickness', 'm')
      call defvar ('O_icevelX', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'eastward ice velocity', 'eastward_sea_ice_velocity'
     &, 'm s-1')
      call defvar ('O_icevelY', iou, 3, iu, -c1e3, c1e3, ' ', 'F'
     &, 'northward ice velocity', 'northward_sea_ice_velocity'
     &, 'm s-1')
      call defvar ('O_iceintX', iou, 3, iu, -c1e20, c1e20, ' ', 'F'
     &, 'eastward ice interaction'
     &, 'downward_eastward_stress_at_sea_ice_base', 'Pa')
      call defvar ('O_iceintY', iou, 3, iu, -c1e20, c1e20, ' ', 'F'
     &, 'northward ice interaction'
     &, 'downward_northward_stress_at_sea_ice_base', 'Pa')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine embm_tavg_out (fname, ids, ide, jds, jde, imt, jmt, nat
     &,                         ncat, xt, yt, xu, yu, dxt, dyt, dxu, dyu
     &,                         avgper, time, stamp, mapat, at, sat
     &,                         precip, evap, disch, vflux, outlwr
     &,                         uplwr, upsens, dnswr, upltnt, solins
     &,                         p_alb, a_alb, s_alb, elev, psno, ws
     &,                         runoff, wx, wy
     &,                         awx, awy, rtbar, apress
     &,                         soilm, surf
     &,                         tice, hice, aice, hsno
     &,                         uice, vice, xint, yint
     &,                         tlat, tlon, ulat, ulon, tgarea, ugarea
     &,                         tmsk, mskhr, nriv, ntrec)
!=======================================================================
!     output routine for atmospheric time averages

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname        = file name
!     ids, ide ... = start and end index for data domain
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     dxt, dyt ... = grid widths

!     time         = time in years
!     stamp        = time stamp
!     at, ...      = data to be written

!   outputs:
!     ntrec        = number of time record in file
!=======================================================================

      implicit none

      integer iou, j, k, ln, n, ntrec, imt, jmt, nat, ncat
      integer ids, ide, jds, jde, igs, ige, ig, jgs, jge, jg
      integer ils, ile, jls, jle, ib(10), ic(10)
      integer mskhr(ids:ide,jds:jde), nriv(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      character(*) :: fname, stamp
      character(3) :: a3
      character(10) :: mapat(nat)

      real xt(imt), xu(imt), yt(jmt), yu(jmt)
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt)
      real avgper, at(ids:ide,jds:jde,nat)
      real sat(ids:ide,jds:jde), precip(ids:ide,jds:jde)
      real evap(ids:ide,jds:jde), disch(ids:ide,jds:jde)
      real vflux(ids:ide,jds:jde), outlwr(ids:ide,jds:jde)
      real uplwr(ids:ide,jds:jde), upsens(ids:ide,jds:jde)
      real dnswr(ids:ide,jds:jde), upltnt(ids:ide,jds:jde)
      real solins(ids:ide,jds:jde), p_alb(ids:ide,jds:jde)
      real a_alb(ids:ide,jds:jde), s_alb(ids:ide,jds:jde)
      real elev(ids:ide,jds:jde), psno(ids:ide,jds:jde)
      real ws(ids:ide,jds:jde), runoff(ids:ide,jds:jde)
      real wx(ids:ide,jds:jde,nat), wy(ids:ide,jds:jde,nat)
      real awx(ids:ide,jds:jde), awy(ids:ide,jds:jde)
      real rtbar(ids:ide,jds:jde), apress(ids:ide,jds:jde)
      real soilm(ids:ide,jds:jde), surf(ids:ide,jds:jde)
      real tice(ids:ide,jds:jde), hice(ids:ide,jds:jde)
      real aice(ids:ide,jds:jde), hsno(ids:ide,jds:jde)
      real uice(ids:ide,jds:jde), vice(ids:ide,jds:jde)
      real xint(ids:ide,jds:jde), yint(ids:ide,jds:jde)
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real ulat(ids:ide,jds:jde), ulon(ids:ide,jds:jde)
      real tgarea(ids:ide,jds:jde), ugarea(ids:ide,jds:jde)
      real tmsk(ids:ide,jds:jde)
      real time, tmp, xt_e(imt+1), xu_e(imt+1), yt_e(jmt+1)
      real yu_e(jmt+1), tmpmask(ids:ide,jds:jde)
      real c0, c1, c10, c100, c1e3, c1e4, p1, C2K, cal2J
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      p1 = 0.1
      C2K = 273.15
      cal2J = 2.389e-05

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     set global write domain size (may be less than global domain)
!-----------------------------------------------------------------------
      igs = 1
      ige = imt
      if (xt(1) + 360. lt. xt(imt)) then
!       assume cyclic boundary
        igs = 2
        ige = imt-1
      endif
      ig  = ige-igs+1
      jgs = 1
      jge = jmt
      do j=2,jmt
        if (yt(j-1) .lt. -90. .and. yt(j) .gt. -90.) jgs = j
        if (yt(j-1) .lt.  90. .and. yt(j) .gt. 90.) jge = j-1
      enddo
      jg  = jge-jgs+1

!-----------------------------------------------------------------------
!     local domain size (minimum of data domain and global write domain)
!-----------------------------------------------------------------------
      ils = max(ids,igs)
      ile = min(ide,ige)
      jls = max(jds,jgs)
      jle = min(jde,jge)

      allocate ( tmpij(ils:ile,jls:jle) )
      allocate ( tmpijm(ils:ile,jls:jle) )

!-----------------------------------------------------------------------
!     write 1d data (t)
!-----------------------------------------------------------------------
      call putvars ('time', iou, ntrec, time, c1, c0)
      call rdstmp (stamp, nyear, nmonth, nday, nhour, nmin, nsec)
      tmp = nyear
      call putvars ('T_year', iou, ntrec, tmp, c1, c0)
      tmp = nmonth
      call putvars ('T_month', iou, ntrec, tmp, c1, c0)
      tmp = nday
      call putvars ('T_day', iou, ntrec, tmp, c1, c0)
      tmp = nhour
      call putvars ('T_hour', iou, ntrec, tmp, c1, c0)
      tmp = nmin
      call putvars ('T_minute', iou, ntrec, tmp, c1, c0)
      tmp = nsec
      call putvars ('T_second', iou, ntrec, tmp, c1, c0)
      call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)

      if (ntrec .eq. 1) then

!-----------------------------------------------------------------------
!       write 1d data (x, y or z)
!-----------------------------------------------------------------------
        allocate ( tmpi(igs:ige) )
        allocate ( tmpj(jgs:jge) )
        allocate ( tmpie(igs:ige+1) )
        allocate ( tmpje(jgs:jge+1) )

        ib(1) = 1
        ic(1) = ig
        tmpi(igs:ige) = xt(igs:ige)
        call putvara ('xt', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = dxt(igs:ige)
        call putvara ('G_dxt', iou, ig, ib, ic, tmpi, c100, c0)
        tmpi(igs:ige) = xu(igs:ige)
        call putvara ('xu', iou, ig, ib, ic, tmpi, c1, c0)
        tmpi(igs:ige) = dxu(igs:ige)
        call putvara ('G_dxu', iou, ig, ib, ic, tmpi, c100, c0)

        ic(1) = jg
        tmpj(jgs:jge) = yt(jgs:jge)
        call putvara ('yt', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyt(jgs:jge)
        call putvara ('G_dyt', iou, jg, ib, ic, tmpj, c100, c0)
        tmpj(jgs:jge) = yu(jgs:jge)
        call putvara ('yu', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyu(jgs:jge)
        call putvara ('G_dyu', iou, jg, ib, ic, tmpj, c100, c0)

        ic(1) = ig + 1
        call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xt_e(igs:ige+1)
        call putvara ('xt_edges', iou, ig+1, ib, ic, tmpie, c1, c0)
        call edge_maker (2, xu_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xu_e(igs:ige+1)
        call putvara ('xu_edges', iou, ig+1, ib, ic, tmpie, c1, c0)

        ic(1) = jg + 1
        call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
        call putvara ('yt_edges', iou, jg+1, ib, ic, tmpje, c1, c0)
        call edge_maker (2, yu_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yu_e(jgs:jge+1)
        call putvara ('yu_edges', iou, jg+1, ib, ic, tmpje, c1, c0)

        deallocate ( tmpi )
        deallocate ( tmpj )
        deallocate ( tmpie )
        deallocate ( tmpje )

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
        tmpij(ils:ile,jls:jle) = elev(ils:ile,jls:jle)
        call putvara ('L_elev', iou, ln, ib, ic, tmpij, c100, c0)
        tmpij(ils:ile,jls:jle) = tmsk(ils:ile,jls:jle)
        call putvara ('G_mskt', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = mskhr(ils:ile,jls:jle)
        call putvara ('G_mskhr', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = nriv(ils:ile,jls:jle)
        call putvara ('G_mskriv', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlat(ils:ile,jls:jle)
        call putvara ('G_latT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlon(ils:ile,jls:jle)
        call putvara ('G_lonT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulat(ils:ile,jls:jle)
        call putvara ('G_latU', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = ulon(ils:ile,jls:jle)
        call putvara ('G_lonU', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tgarea(ils:ile,jls:jle)
        call putvara ('G_areaT', iou, ln, ib, ic, tmpij, c1e4, c0)
        tmpij(ils:ile,jls:jle) = ugarea(ils:ile,jls:jle)
        call putvara ('G_areaU', iou, ln, ib, ic, tmpij, c1e4, c0)

      endif

!-----------------------------------------------------------------------
!     write 3d data (x,y,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ib(3) = ntrec
      ic(3) = 1
      ln = ic(1)*ic(2)*ic(3)
      do n=1,nat
        if (trim(mapat(n)) .eq. 'sat') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_slat', iou, ln, ib, ic, tmpij, c1, -C2K)
        elseif (trim(mapat(n)) .eq. 'shum') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_shum', iou, ln, ib, ic, tmpij, c1, c0)
        elseif (trim(mapat(n)) .eq. 'co2') then
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_co2', iou, ln, ib, ic, tmpij, c1, c0)
        else
          if (n .lt. 1000) write(a3, '(i3)') n
          if (n .lt. 100) write(a3, '(i2)') n
          if (n .lt. 10) write(a3, '(i1)') n
          tmpij(ils:ile,jls:jle) = at(ils:ile,jls:jle,n)
          call putvara ('A_tracer'//trim(a3), iou, ln, ib, ic, tmpij
     &,     c1, c0)
        endif
      enddo
      tmpij(ils:ile,jls:jle) = sat(ils:ile,jls:jle)
      call putvara('A_sat', iou, ln, ib, ic, tmpij, c1, -C2K)
      tmpij(ils:ile,jls:jle) = precip(ils:ile,jls:jle)
      call putvara ('F_precip', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = psno(ils:ile,jls:jle)
      call putvara ('F_snow', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = evap(ils:ile,jls:jle)
      call putvara ('F_evap', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = disch(ils:ile,jls:jle)
      call putvara ('F_rivdis', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = vflux(ils:ile,jls:jle)
      call putvara ('F_virtual', iou, ln, ib, ic, tmpij, p1, c0)
      tmpij(ils:ile,jls:jle) = outlwr(ils:ile,jls:jle)
      call putvara ('F_outlwr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = uplwr(ils:ile,jls:jle)
      call putvara ('F_uplwr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = upsens(ils:ile,jls:jle)
      call putvara ('F_upsens', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = dnswr(ils:ile,jls:jle)
      call putvara ('F_dnswr', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = upltnt(ils:ile,jls:jle)
      call putvara ('F_upltnt', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpij(ils:ile,jls:jle) = solins(ils:ile,jls:jle)
      call putvara ('F_solins', iou, ln, ib, ic, tmpij, c1e3, c0)
      tmpmask(:,:) = 1.
      where (p_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = p_alb(ils:ile,jls:jle)
      call putvaramsk ('A_albplt', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpmask(:,:) = 1.
      where (a_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = a_alb(ils:ile,jls:jle)
      call putvaramsk ('A_albatm', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpmask(:,:) = 1.
      where (s_alb(:,:) .lt. 0.) tmpmask(:,:) = 0.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = s_alb(ils:ile,jls:jle)
      call putvaramsk ('A_albsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = awx(ils:ile,jls:jle)
      call putvara ('A_awindX', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) =  awy(ils:ile,jls:jle)
      call putvara ('A_awindY', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = rtbar(ils:ile,jls:jle)
      call putvara ('A_avgslat', iou, ln, ib, ic, tmpij, c1, c0)
      call putvara ('A_avgslat', iou, ln, ib, ic, tmpij, c1, -C2K)
      tmpij(ils:ile,jls:jle) =  apress(ils:ile,jls:jle)
      call putvara ('A_apress', iou, ln, ib, ic, tmpij, c10, c0)
      tmpmask(:,:) = 0.
      where (tmsk(:,:) .lt. 0.5) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = soilm(ils:ile,jls:jle)
      call putvaramsk ('L_soilmois', iou, ln, ib, ic, tmpij, tmpijm
     &, p1, c0)
      tmpij(ils:ile,jls:jle) = surf(ils:ile,jls:jle)
      call putvaramsk ('L_tempsur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, -C2K)
      tmpmask(:,:) = 0.
      where (hice(:,:) .gt. 0) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      tmpij(ils:ile,jls:jle) = tice(ils:ile,jls:jle)
      call putvaramsk ('O_icetemp', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, -C2K)
      tmpij(ils:ile,jls:jle) = hice(ils:ile,jls:jle)
      call putvara ('O_icethk', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = aice(ils:ile,jls:jle)
      call putvara ('O_icefra', iou, ln, ib, ic, tmpij, c1, c0)
      tmpij(ils:ile,jls:jle) = hsno(ils:ile,jls:jle)
      call putvara ('L_snothk', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = uice(ils:ile,jls:jle)
      call putvara ('O_icevelX', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = vice(ils:ile,jls:jle)
      call putvara ('O_icevelY', iou, ln, ib, ic, tmpij, c100, c0)
      tmpij(ils:ile,jls:jle) = xint(ils:ile,jls:jle)
      call putvara ('O_iceintX', iou, ln, ib, ic, tmpij, c10, c0)
      tmpij(ils:ile,jls:jle) = yint(ils:ile,jls:jle)
      call putvara ('O_iceintY', iou, ln, ib, ic, tmpij, c10, c0)

      deallocate ( tmpij )
      deallocate ( tmpijm )

      return
      end
