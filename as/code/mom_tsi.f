! source file: /net/mare/home1/eby/as/2.9.old3/source/mom/mom_tsi.F
      subroutine mom_tsi_def (fname, calendar, expnam, runstamp)

!=======================================================================
!     output routine for ocean time step integrals

!   inputs:
!     fname      = file name
!     calendar   = calendar
!     expnam     = experiment name
!     runstamp   = run stamp
!=======================================================================

      implicit none

      character(*) :: fname, calendar, expnam, runstamp

      integer id(1), id_time, iou, ntrec

      real c0, c1, c10, c100, c500, c1e3, c1e6, c1e20

      c0 = 0.
      c1 = 1.
      c10 = 10.
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
      id = id_time

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
      call defvar ('O_ke', iou, 1, id, c0, c1e6, ' ', 'F'
     &, 'kinetic energy per unit volume', ' ', 'J m-3')
      call defvar ('O_temp', iou, 1, id, -c100, c500, ' ', 'F'
     &, 'global average ocean temperature', ' ', 'K')
      call defvar ('O_sal', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean salinity', ' ', '1e-3')
      call defvar ('O_tempvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'variance of ocean temperature', ' ', 'K2')
      call defvar ('O_salvar', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'variance of ocean salinity', ' ', '1e-6')
      call defvar ('O_absdtemp', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'rate of change of ocean temperature', ' ', 'K s-1')
      call defvar ('O_absdsal', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'rate of change of ocean salinity', ' ', '1e-3 s-1')
      call defvar ('O_avgiter', iou, 1, id, c0, c1e6, ' ', 'F'
     &, 'ocean solver iterations', ' ', '1')
      call defvar ('F_heat', iou, 1, id, -c1e20, c1e20, ' '
     &,   'F', 'global average ocean heat flux', ' ','W m-2')
      call defvar ('F_salt', iou, 1, id, -c1e20, c1e20, ' '
     &,   'F', 'global average ocean salt flux', ' ','kg m-2 s-1')
      call defvar ('O_dic', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean dic', ' ', 'mol m-3')
      call defvar ('F_dic', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean dic flux', ' ', 'mol m-2 s-1')
      call defvar ('O_alk', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean alkalinity', ' ', 'mol m-3')
      call defvar ('O_o2', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean oxygen', ' ', 'mol m-3')
      call defvar ('F_o2', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean oxygen flux', ' ', 'mol m-2 s-1')
      call defvar ('O_po4', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean phosphate', ' ', 'mol m-3')
      call defvar ('O_phyt', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean phytoplankton', ' ', 'mol N m-3')
      call defvar ('O_zoop', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean zooplankton', ' ', 'mol N m-3')
      call defvar ('O_detr', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean detritus', ' ', 'mol N m-3')
      call defvar ('O_no3', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean nitrate', ' ', 'mol m-3')
      call defvar ('O_diaz', iou, 1, id, c0, c100, ' ', 'F'
     &, 'global average ocean diazotrophs', ' ', 'mol N m-3')
      call defvar ('O_motmax', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'maximum meridional overturning streamfunction', ' '
     &, 'm3 s-1')
      call defvar ('O_motmin', iou, 1, id, c0, c1e20, ' '
     &, 'F', 'minimum meridional overturning streamfunction', ' '
     &, 'm3 s-1')
      call defvar ('O_dsealev', iou, 1, id, c0, c1e20, ' ', 'F'
     &, 'relative sea level height', ' ', 'm')

!-----------------------------------------------------------------------
!     end definitions
!-----------------------------------------------------------------------
      call enddef (iou)

      return
      end

      subroutine mom_tsi_out (fname, avgper, time, stamp, ektot, tbar1
     &,                       tbar2, travar1, travar2, dtabs1, dtabs2
     &,                       scan, hflx, sflx, tbar_dic, dicflx
     &,                       tbar_alk, tbar_o2, o2flx, tbar_po4, tbar_p
     &,                       tbar_z, tbar_d, tbar_no3, tbar_di
     &,                       tbar_c14, tbar_dc14, c14flx, tbar_cfc11
     &,                       cfc11flx, tbar_cfc12, cfc12flx, otmax
     &,                       otmin, slh, sspH, ssCO3, ssOc, ssOa
     &,                       sspCO2, cocn, cfa2o, ntrec)
!=======================================================================
!     output routine for ocean time step integrals

!   inputs:
!     fname      = file name
!     avgper     = length of averaging period
!     time       = time in years
!     stamp      = time stamp
!     ektot, ... = data to be written

!   outputs:
!     ntrec      = number of time record in file
!=======================================================================

      implicit none

      character(*) :: fname, stamp

      integer iou, ntrec, nyear, nmonth, nday, nhour, nmin, nsec

      real ektot, tbar1, tbar2, travar1, travar2, dtabs1, dtabs2, scan
      real hflx, sflx, tbar_dic, dicflx, tbar_alk, tbar_o2, o2flx
      real tbar_po4, tbar_p, tbar_z, tbar_d, tbar_no3, tbar_di, tbar_c14
      real tbar_dc14, c14flx, tbar_cfc11, cfc11flx, tbar_cfc12, cfc12flx
      real otmax, otmin, slh, sspH, ssCO3, ssOc, ssOa, sspCO2, cocn
      real cfa2o, avgper, time, tmp, c0, c1, c10, c100, c1e3, c1e6
      real p001, p035, p1, C2K, cal2J

      c0 = 0.
      c1 = 1.
      c10 = 10.
      c100 = 100.
      c1e3 = 1.e3
      c1e6 = 1.e6
      p001 = 0.001
      p035 = 0.035
      p1 = 0.1
      C2K = 273.15
      cal2J = 2.389e-05

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
      call putvars ('O_ke', iou, ntrec, ektot, c10, c0)
      call putvars ('O_temp', iou, ntrec, tbar1, c1, -C2K)
      call putvars ('O_sal', iou, ntrec, tbar2, p001, -p035)
      tmp = 0.0
      if (travar1 .gt. 0.0) tmp = (sqrt(travar1) + C2K)**2
      call putvars ('O_tempvar', iou, ntrec, tmp, c1, c0)
      call putvars ('O_salvar', iou, ntrec, tmp, c1, c0)
      call putvars ('O_absdtemp', iou, ntrec, dtabs1, c1, c0)
      call putvars ('O_absdsal', iou, ntrec, dtabs2, p001, c0)
      call putvars ('O_avgiter', iou, ntrec, scan, c1, c0)
      call putvars ('F_heat', iou, ntrec, hflx, cal2J, c0)
      call putvars ('F_salt', iou, ntrec, sflx, p1, c0)
      call putvars ('O_dic', iou, ntrec, tbar_dic, c1, c0)
      call putvars ('F_dic', iou, ntrec, dicflx, c100, c0)
      call putvars ('O_alk', iou, ntrec, tbar_alk, c1, c0)
      call putvars ('O_o2', iou, ntrec, tbar_o2, c1, c0)
      call putvars ('F_o2', iou, ntrec, o2flx, c100, c0)
      call putvars ('O_po4', iou, ntrec, tbar_po4, c1e3, c0)
      call putvars ('O_phyt', iou, ntrec, tbar_p, c1e3, c0)
      call putvars ('O_zoop', iou, ntrec, tbar_z, c1e3, c0)
      call putvars ('O_detr', iou, ntrec, tbar_d, c1e3, c0)
      call putvars ('O_no3', iou, ntrec, tbar_no3, c1e3, c0)
      call putvars ('O_diaz', iou, ntrec, tbar_di, c1e3, c0)
      call putvars ('O_motmax', iou, ntrec, otmax, c1e6, c0)
      call putvars ('O_motmin', iou, ntrec, otmin, c1e6, c0)
      call putvars ('O_dsealev', iou, ntrec, slh, c100, c0)

      return
      end
