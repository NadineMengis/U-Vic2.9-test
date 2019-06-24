! source file: /net/mare/home1/eby/as/2.9.old3/source/mtlm/mtlm_tavg.F
      subroutine mtlm_tavg_def (fname, imt, jmt, NPFT, NTYPE, xt, yt
     &,                         calendar, expnam, runstamp)

!=======================================================================
!     definition routine for land time averages

!   inputs:
!     fname        = file name
!     imt, jmt ... = global array dimensions
!     xt, yt ...   = global axes
!     calendar     = calendar
!     expnam       = experiment name
!     runstamp     = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer iou, j, imt, jmt, igs, ige, ig, jgs, jge, jg
      integer id_time, id_xt, id_yt, id_pft, id_type
      integer id_xt_e, id_yt_e, id_pft_e, id_type_e
      integer it(10)

      integer NPFT, NTYPE

      real xt(imt), yt(jmt)

      real c0, c1, c1e6, c1e20

      c0 = 0.
      c1 = 1.
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
      call defdim ('pft', iou, NPFT, id_pft)
      call defdim ('type', iou, NTYPE, id_type)
      call defdim ('xt_edges', iou, ig+1, id_xt_e)
      call defdim ('yt_edges', iou, jg+1, id_yt_e)
      call defdim ('pft_edges', iou, NPFT+1, id_pft_e)
      call defdim ('type_edges', iou, NTYPE+1, id_type_e)

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      it(1) = id_time
      call defvar ('time', iou, 1, it, c0, c0, 'T', 'D'
     &, 'time', 'time', 'year since 0000-01-01 00:00:00')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_year', iou, 1, it, c0, c0, ' ', 'F'
     &,   'T_year', ' ','year')
      call defvar ('T_month', iou, 1, it, c0, c0, ' ', 'F'
     &,   'T_month', ' ','month')
      call defvar ('T_day', iou, 1, it, c0, c0, ' ', 'F'
     &,   'T_day', ' ','day')
      call defvar ('T_hour', iou, 1, it, c0, c0, ' ', 'F'
     &,   'T_hour', ' ','hour')
      call defvar ('T_minute', iou, 1, it, c0, c0, ' ', 'F'
     &,   'T_minute', ' ','minute')
      call defvar ('T_second', iou, 1, it, c0, c0, ' ', 'F'
     &,   'T_second', ' ','second')
      call defvar ('T_avgper', iou, 1, it, c0, c0, ' ', 'F'
     &,   'averaging period', ' ','day')

!-----------------------------------------------------------------------
!     define 1d data (x, y or n)
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
      it(1) = id_pft
      call defvar ('pft', iou, 1, it, c0, c0, 'Z', 'F'
     &,   'plant functional type', ' ', ' ')
      it(1) = id_type
      call defvar ('type', iou, 1, it, c0, c0, 'Z', 'F'
     &,   'land type', ' ', ' ')
      it(1) = id_xt_e
      call defvar ('xt_edges', iou, 1, it, c0, c0, 'X', 'D'
     &,   'longitude of t grid edges', ' ', 'degrees')
      it(1) = id_yt_e
      call defvar ('yt_edges', iou, 1, it, c0, c0, 'Y', 'D'
     &,   'latitude of t grid edges', ' ', 'degrees')
      it(1) = id_pft_e
      call defvar ('pft_edges', iou, 1, it, c0, c0, 'Z', 'F'
     &,   'plant functional type', ' ', ' ')
      it(1) = id_type_e
      call defvar ('type_edges', iou, 1, it, c0, c0, 'Z'
     &,   'F', 'land type', ' ', ' ')

!-----------------------------------------------------------------------
!     define 2d data (x,y)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      call defvar ('G_latT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &,   'tracer grid latitude', 'latitude', 'degrees')
       call defvar ('G_lonT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &,   'tracer grid longitude', 'longitude', 'degrees')
      call defvar ('G_areaT', iou, 2, it, -c1e6, c1e6, ' ', 'F'
     &, 'tracer grid area', ' ', 'm2')

!-----------------------------------------------------------------------
!     define 3d data (x,y,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_time

      call defvar ('L_soiltemp', iou , 3, it, -c1e20, c1e20, ' ', 'F'

     &,   'soil temperature', 'soil_temperature', 'K')

      call defvar ('L_soilcarb', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'soil carbon', 'soil_carbon_content', 'kg m-2 ')
      call defvar ('L_soilresp', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'soil respiration', 'soil_respiration_carbon_flux'
     &,   'kg m-2 s-1')
      call defvar ('L_veglit', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'total leaf litter', 'leaf_litter_carbon_flux'
     &,   'kg m-2 s-1')
      call defvar ('L_vegburn', iou , 3, it, -c1e20, c1e20, ' ', 'F'
     &,   'total vegetation burning', 'burning_carbon_flux'
     &,   'kg m-2 s-1')

!-----------------------------------------------------------------------
!     define 4d data (x,y,pft,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_pft
      it(4) = id_time
      call defvar ('L_veggpp', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'gross primary productivity'
     &,   'gross_primary_productivity_of_carbon', 'kg m-2 s-1')
      call defvar ('L_vegnpp', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'net primary productivity'
     &,   'net_primary_productivity_of_carbon', 'kg m-2 s-1')
      call defvar ('L_veghgt', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'canopy height', 'canopy_height', 'L_soilmois')
      call defvar ('L_veglai', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'leaf area index ', 'leaf_area_index', '1')
      call defvar ('L_vegcarb', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &,   'vegetation carbon', 'vegetation_carbon_content', 'kg m-2 ')

!-----------------------------------------------------------------------
!     define 4d data (x,y,type,t)
!-----------------------------------------------------------------------
      it(1) = id_xt
      it(2) = id_yt
      it(3) = id_type
      it(4) = id_time
       call defvar ('L_vegfra', iou , 4, it, -c1e20, c1e20, ' ', 'F'
     &, 'areal coverage', 'vegetation_area_fraction', '1')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mtlm_tavg_out (fname, ids, ide, jds, jde, imt, jmt
     &,                         POINTS, NPFT, NTYPE, xt, yt, xu, yu
     &,                         dxt, dyt, dxu, dyu, avgper, time
     &,                         stamp, land_map, TS1,  CS, RESP_S
     &,                         LIT_C_T, BURN, FRAC, GPP, NPP, HT, LAI
     &,                         C_VEG
     &,                         tlat, tlon, tgarea, ntrec)
!=======================================================================
!     output routine for land time averages

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
!     avgper       = length of averaging period
!     time         = time in years
!     stamp        = time stamp
!     land_map     = land map
!     TS1, ...     = data to be written

!   outputs:
!     ntrec        = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, j, ln, n, ntrec, imt, jmt, ids, ide, jds, jde, igs
      integer ige, ig, jgs, jge, jg, ils, ile, jls, jle, ib(10), ic(10)
      integer nyear, nmonth, nday, nhour, nmin, nsec
      integer POINTS, NPFT, NTYPE, land_map(imt,jmt)

      real xt(imt), xu(imt), yt(jmt), yu(jmt)
      real dxt(imt), dxu(imt), dyt(jmt), dyu(jmt)
      real tmpmask(imt,jmt), data(imt,jmt), pft(NPFT), type(NTYPE)
      real xt_e(imt+1), yt_e(jmt+1), pft_e(NPFT+1), type_e(NTYPE+1)
      real avgper, time, tmp, c0, c1, c100, c1e4, C2K

      real TS1(POINTS), CS(POINTS), RESP_S(POINTS), LIT_C_T(POINTS)
      real BURN(POINTS), FRAC(POINTS,NTYPE), GPP(POINTS,NPFT)
      real NPP(POINTS,NPFT), HT(POINTS,NPFT), LAI(POINTS,NPFT)
      real C_VEG(POINTS,NPFT)
      real tlat(ids:ide,jds:jde), tlon(ids:ide,jds:jde)
      real tgarea(ids:ide,jds:jde)
      real, allocatable :: tmpij(:,:), tmpijm(:,:)
      real, allocatable :: tmpi(:), tmpj(:)
      real, allocatable :: tmpie(:), tmpje(:)

      c0 = 0.
      c1 = 1.
      c100 = 100.
      c1e4 = 1.e4
      C2K = 273.15

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
      if (avgper .gt. 1.e-6) then
        call putvars ('T_avgper', iou, ntrec, avgper, c1, c0)
      endif

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

        do n=1, NPFT
          pft(n) = n
          pft_e(n) = pft(n) - 0.5
        enddo
        pft_e(NPFT+1) = pft(NPFT) + 0.5
        ic(1) = NPFT
        call putvara ('pft', iou, NPFT, ib, ic, pft, c1, c0)
        ic(1) = NPFT + 1
        call putvara ('pft_edges', iou, NPFT+1, ib, ic, pft_e, c1, c0)
        do n=1, NTYPE
          type(n) = float(n)
          type_e(n) = type(n) - 0.5
        enddo
        type_e(NTYPE+1) = type(NTYPE) + 0.5
        ic(1) = NTYPE
        call putvara ('type', iou, NTYPE, ib, ic, type, c1, c0)
        ic(1) = NTYPE + 1
        call putvara ('type_edges', iou, NTYPE+1, ib, ic, type_e
     &,   c1, c0)

!-----------------------------------------------------------------------
!       write 2d data (x,y)
!-----------------------------------------------------------------------
        ib(1) = ils-igs+1
        ic(1) = ile-ils+1
        ib(2) = jls-jgs+1
        ic(2) = jle-jls+1
        ln = ic(1)*ic(2)
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
      tmpmask(:,:) = 0.
      where (land_map(:,:) .gt. 0) tmpmask(:,:) = 1.
      tmpijm(ils:ile,jls:jle) = tmpmask(ils:ile,jls:jle)
      call unloadland (POINTS, TS1, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('L_soiltemp', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      call unloadland (POINTS, CS, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('L_soilcarb', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      call unloadland (POINTS, RESP_S, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('L_soilresp', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      call unloadland (POINTS, LIT_C_T, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('L_veglit', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)
      call unloadland (POINTS, BURN, imt, jmt, land_map, data)
      tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
      call putvaramsk ('L_vegburn', iou, ln, ib, ic, tmpij, tmpijm
     &, c1, c0)

!-----------------------------------------------------------------------
!     write 4d data (x,y,pft,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ic(3) = 1
      ib(4) = ntrec
      ic(4) = 1
      ln = ic(1)*ic(2)*ic(3)*ic(4)
      do n=1,npft
        ib(3) = n
        call unloadland (POINTS, GPP(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('L_veggpp', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
        call unloadland (POINTS, NPP(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('L_vegnpp', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
        call unloadland (POINTS, HT(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('L_veghgt', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
        call unloadland (POINTS, LAI(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('L_veglai', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
        call unloadland (POINTS, C_VEG(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        call putvaramsk ('L_vegcarb', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
      enddo

!-----------------------------------------------------------------------
!     write 4d data (x,y,type,t)
!-----------------------------------------------------------------------
      ib(1) = ils-igs+1
      ic(1) = ile-ils+1
      ib(2) = jls-jgs+1
      ic(2) = jle-jls+1
      ic(3) = 1
      ib(4) = ntrec
      ic(4) = 1
      do n=1,ntype
        ib(3) = n
        call unloadland (POINTS, FRAC(1,n), imt, jmt, land_map, data)
        tmpij(ils:ile,jls:jle) = data(ils:ile,jls:jle)
        where (tmpij(ils:ile,jls:jle) .lt. 0.)
     &    tmpij(ils:ile,jls:jle) = c0
        call putvaramsk ('L_vegfra', iou, ln, ib, ic, tmpij, tmpijm
     &,   c1, c0)
      enddo

      deallocate ( tmpij )
      deallocate ( tmpijm )

      return
      end
