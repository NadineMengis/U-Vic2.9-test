      subroutine sed_tavg_def (fname, imt, jmt, xt, yt, calendar, expnam
     &,                        runstamp)

#if defined O_sed
!=======================================================================
!     definition routine for ocean time averages

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     calendar     = calendar
!     expnam       = experiment name
!     runstamp     = run stamp
!     mapt         = tracer map
!=======================================================================

      implicit none

      integer iou, j, imt, jmt, igs, ige, ig, jgs, jge, jg
      integer it(10), id_time, id_xt, id_yt, id_xt_e, id_yt_e
      character(*) :: fname, calendar, expnam, runstamp

      real xt(imt), yt(jmt)
      real c0, c1, c100, c500, c1e3, c1e4, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c500 = 500.
      c1e3 = 1.e3
      c1e4 = 1.e4
      c1e6 = 1.e6
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
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
!     write global atributes
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
      call defdim ('xt_edges', iou, ig+1, id_xt_e)
      call defdim ('yt_edges', iou, jg+1, id_yt_e)

!-----------------------------------------------------------------------
!     define 1d data (t)
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
!      call defvar ('T_avgper', iou, 1, it, c0, c0, ' ', 'F'
!     &, 'averaging period', ' ','day')

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
      it(1) = id_xt_e
      call defvar ('xt_edges', iou, 1, it, c0, c0, 'X', 'D'
     &,   'longitude of t grid edges', ' ', 'degrees')
      it(1) = id_yt_e
      call defvar ('yt_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &,   'latitude of t grid edges', ' ', 'degrees')

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      call defvar ('G_kmt', iou, 2, it, c0, c1e6, ' ', 'I'
     &, 'ocean grid depth level', 'model_level_number' ,'1')
      call defvar ('G_latT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid latitude', 'latitude', 'degrees')
      call defvar ('G_lonT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid longitude', 'longitude', 'degrees')
      call defvar ('G_areaT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid area', ' ', 'm2')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_time

      call defvar ('F_caco3dis', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'sediment surface upward flux of calcite', ' ', 'mol m-2 s-1')
      call defvar ('F_caco3rai', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'sediment surface downward flux of calcite', ' ', 'mol m-2 s-1')
      call defvar ('S_caco3per', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'pore layer CaCO3 portion', ' ', 'percent')
      call defvar ('S_caco3mas', iou, 3, it, -c1e20, c1e20, ' ', 'F'
     &, 'pore layer CaCO3 mass', ' ', 'kg C m-2')
      call defvar ('S_caco3bur', iou, 3, it, -c1e20, c1e20, ' ', 'F'
     &, 'buried CaCO3 mass', ' ', 'kg C m-2')
      call defvar ('S_co3', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'sediment surface CO3 concentration', ' ', 'mol m-3')
      call defvar ('S_co3sat', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'sediment surface CO3 saturation', ' ', 'mol m-3')
      call defvar ('S_rainrat', iou, 3, it, -c1e6, c1e6, ' ', 'F'
     &, 'rain ratio', ' ', '1')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine sed_tavg_out (fname, ids, ide, jds, jde, imt, jmt, xt
     &,                        yt, xu, yu, dxt, dyt, dxu, dyu, avgper
     &,                        time, stamp, ttrcal, rain_cal, cal
     &,                        calmass, calmass_bur, co3, co3sat, rainr
     &,                        map_sed, kmt, tlat, tlon, tgarea, ntrec)
!=======================================================================
!     output routine for ocean time averages

!     data may be sized differently in x and y from the global fields.
!     fields may be written with or without a time dimension. data
!     should be defined with the routine defvar and written with putvar.
!     if  no time dimension, then data is only written once per file.
!     make sure the it, iu, ib, and ic arrays and are defining the
!     correct dimensions. ln may also need to be recalculated.

!   inputs:
!     fname        = file name
!     ids, ide ... = start and end index for data domain
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     dxt, dyt ... = grid widths
!     avgper       = length of averaging period
!     time         = time in years
!     t, ...       = data to be written

!   outputs:
!     ntrec        = number of time record in file
!=======================================================================

      implicit none

      integer iou, j, ln, ntrec, imt, jmt, ids, ide, jds, jde, igs
      integer ige, ig, jgs, jge, jg, ils, ile, jls, jle, ib(10), ic(10)
      integer kmt(ids:ide,jds:jde), map_sed(ids:ide,jds:jde)
      integer nyear, nmonth, nday, nhour, nmin, nsec

      character(*) :: fname, stamp

      real xt(imt), xu(imt), yt(jmt), yu(jmt)
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt), avgper
      real ttrcal(ids:ide,jds:jde), rain_cal(ids:ide,jds:jde)
      real cal(ids:ide,jds:jde), calmass(ids:ide,jds:jde)
      real calmass_bur(ids:ide,jds:jde), co3(ids:ide,jds:jde)
      real co3sat(ids:ide,jds:jde), rainr(ids:ide,jds:jde)

      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real tgarea(ids:ide,jds:jde),  xt_e(imt+1), yt_e(jmt+1)
      real time, tmp, c0, c1, c100, c1e3, c1e4, p01, p0001

      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpi(:), tmpj(:), tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c1e3 = 1.e3
      c1e4 = 1.e4
      p01 = 0.01
      p0001 = 0.0001

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

!-----------------------------------------------------------------------
!     global write domain size (may be less than global domain)
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
!      if (avgper .gt. 1.e-6) then
!        call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)
!      endif

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
        ic(1) = jg
        tmpj(jgs:jge) = yt(jgs:jge)
        call putvara ('yt', iou, jg, ib, ic, tmpj, c1, c0)
        tmpj(jgs:jge) = dyt(jgs:jge)
        call putvara ('G_dyt', iou, jg, ib, ic, tmpj, c100, c0)
        ic(1) = ig + 1
        call edge_maker (1, xt_e, xt, dxt, xu, dxu, imt)
        tmpie(igs:ige+1) = xt_e(igs:ige+1)
        call putvara ('xt_edges', iou, ig+1, ib, ic, tmpie, c1, c0)
        ic(1) = jg + 1
        call edge_maker (1, yt_e, yt, dyt, yu, dyu, jmt)
        tmpje(jgs:jge+1) = yt_e(jgs:jge+1)
        call putvara ('yt_edges', iou, jg+1, ib, ic, tmpje, c1, c0)

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
        tmpij(ils:ile,jls:jle) = kmt(ils:ile,jls:jle)
        call putvara ('G_kmt', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlat(ils:ile,jls:jle)
        call putvara ('G_latT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tlon(ils:ile,jls:jle)
        call putvara ('G_lonT', iou, ln, ib, ic, tmpij, c1, c0)
        tmpij(ils:ile,jls:jle) = tgarea(ils:ile,jls:jle)
        call putvara ('G_areaT', iou, ln, ib, ic, tmpij, c1e4, c0)

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
      tmpijm(:,:) = 0.
      where (map_sed(ils:ile,jls:jle) .gt. 0)
     &  tmpijm(ils:ile,jls:jle) = 1.
      tmpij(ils:ile,jls:jle) = ttrcal(ils:ile,jls:jle)
      call putvaramsk ('F_caco3dis', iou, ln, ib, ic, tmpij, tmpijm
     &, p0001, c0)
      tmpij(ils:ile,jls:jle) = rain_cal(ils:ile,jls:jle)
      call putvaramsk ('F_caco3rai', iou, ln, ib, ic, tmpij, tmpijm
     &, p0001, c0)
      tmpij(ils:ile,jls:jle) = cal(ils:ile,jls:jle)
      call putvaramsk ('S_caco3per', iou, ln, ib, ic, tmpij, tmpijm
     &, p01, c0)
      tmpij(ils:ile,jls:jle) = calmass(ils:ile,jls:jle)
      call putvaramsk ('S_caco3mas', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = calmass_bur(ils:ile,jls:jle)
      call putvaramsk ('S_caco3bur', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      tmpij(ils:ile,jls:jle) = co3(ils:ile,jls:jle)
      call putvaramsk ('S_co3', iou, ln, ib, ic, tmpij, tmpijm
     &, c1e3, c0)
      tmpij(ils:ile,jls:jle) = co3sat(ils:ile,jls:jle)
      call putvaramsk ('S_co3sat', iou, ln, ib, ic, tmpij, tmpijm
     &, c1e3, c0)
      tmpij(ils:ile,jls:jle) = rainr(ils:ile,jls:jle)
      call putvaramsk ('S_rainrat', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)

      deallocate (tmpij)
      deallocate (tmpijm)
#endif

      return
      end
