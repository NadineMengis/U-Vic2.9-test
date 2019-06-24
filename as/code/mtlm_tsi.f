! source file: /net/mare/home1/eby/as/2.9.old3/source/mtlm/mtlm_tsi.F
      subroutine mtlm_tsi_def (fname, calendar, expnam, runstamp)

!=======================================================================
!     output routine for land time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou, ntrec

      real c0, c1, c1e20

      c0 = 0.
      c1 = 1.
      c1e20 = 1.e20

!-----------------------------------------------------------------------
!     open file
!-----------------------------------------------------------------------
      call openfile (fname, iou)

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
      id(1) = id_time

!-----------------------------------------------------------------------
!     define 1d data (t)
!-----------------------------------------------------------------------
      call defvar ('time', iou, 1, id, c0, c0, 'T', 'D'
     &, 'time', 'time', 'year since 0000-01-01 00:00:00')
      call putatttext (iou, 'time', 'calendar', calendar)
      call defvar ('T_year', iou, 1, id, c0, c0, ' ', 'F'
     &, 'T_year', ' ','year')
      call defvar ('T_month', iou, 1, id, c0, c0, ' ', 'F'
     &, 'T_month', ' ','month')
      call defvar ('T_day', iou, 1, id, c0, c0, ' ', 'F'
     &, 'T_day', ' ','day')
      call defvar ('T_hour', iou, 1, id, c0, c0, ' ', 'F'
     &, 'T_hour', ' ','hour')
      call defvar ('T_minute', iou, 1, id, c0, c0, ' ', 'F'
     &, 'T_minute', ' ','minute')
      call defvar ('T_second', iou, 1, id, c0, c0, ' ', 'F'
     &, 'T_second', ' ','second')
      call defvar ('T_avgper', iou, 1, id, c0, c0, ' ', 'F'
     &, 'averaging period', ' ','day')
      call defvar ('L_soiltemp', iou, 1, id, -c1e20, c1e20, ' '

     &, 'F', 'global average soil temperature', ' ', 'K')

      call defvar ('L_soilcarb', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total soil carbon', ' ', 'kg')
      call defvar ('L_soilresp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total soil respiration flux', ' ', 'kg s-1')
      call defvar ('L_veglit', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total leaf litter flux', ' ', 'kg s-1')
      call defvar ('L_vegburn', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total vegetation burning flux', ' ', 'kg s-1')
      call defvar ('L_vegcarb', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'global total vegetation carbon', ' ', 'kg')
      call defvar ('L_vegnpp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total net primary productivity', ' ', 'kg s-1')
      call defvar ('L_veggpp', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global total gross primary productivity', ' ', 'kg s-1')
      call defvar ('L_veghgt', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average vegetation height', ' ', 'L_soilmois')
      call defvar ('L_veglai', iou, 1, id, -c1e20, c1e20, ' '
     &, 'F', 'global average leaf area index', ' ', '1')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mtlm_tsi_out (fname, avgper, time, stamp, CS, RESP_S
     &,                        LIT_C_T, BURN, CV, NPP, GPP, HT, LAI
     &,                        LYING_SNOW, TSOIL, TSTAR, M, ET, clnd
     &,                        cfa2l, ntrec)
!=======================================================================
!     output routine for land time step integrals

!   inputs:
!     fname   = file name
!     avgper  = length of averaging period
!     time    = time in years
!     stamp   = time stamp
!     CS, ... = data to be written

!   outputs:
!     ntrec   = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real CS, RESP_S, LIT_C_T, BURN, CV, NPP, GPP, HT, LAI, LYING_SNOW
      real TSOIL, TSTAR, M, ET, clnd, cfa2l, avgper, time, tmp
      real c0, c1, C2K

      c0 = 0.
      c1 = 1.
      C2K = 273.15

!-----------------------------------------------------------------------
!     open file and get latest record number
!-----------------------------------------------------------------------
      call opennext (fname, time, ntrec, iou)
      if (ntrec .le. 0) ntrec = 1

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
      call putvars ('L_soiltemp', iou, ntrec, TSOIL, c1, c0)
      call putvars ('L_soilcarb', iou, ntrec, CS, c1, c0)
      call putvars ('L_soilresp', iou, ntrec, RESP_S, c1, c0)
      call putvars ('L_veglit', iou, ntrec, LIT_C_T, c1, c0)
      call putvars ('L_vegburn', iou, ntrec, BURN, c1, c0)
      call putvars ('L_vegcarb', iou, ntrec, CV, c1, c0)
      call putvars ('L_vegnpp', iou, ntrec, NPP, c1, c0)
      call putvars ('L_veggpp', iou, ntrec, GPP, c1, c0)
      call putvars ('L_veghgt', iou, ntrec, HT, c1, c0)
      call putvars ('L_veglai', iou, ntrec, LAI, c1, c0)

      return
      end
