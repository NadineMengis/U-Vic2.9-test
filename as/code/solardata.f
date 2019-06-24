! source file: /net/mare/home1/eby/as/2.9.old3/source/common/solardata.F
      subroutine solardata

!=======================================================================
!     routine to read and interpolate solar forcing data
!=======================================================================

      implicit none

      character(120) :: fname, name, new_file_name

      integer iou, n, ln, ib(10), ic(10)

      logical exists

      real dat(3), data_time, tim(3)
      real c1e3, wt1, wt3

      real, allocatable :: data(:), time(:)

      save dat, data, ln, tim, time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "cembm.h"
      include "tmngr.h"

      c1e3 = 1.e3
      if (.not. allocated (time)) then
        fname = "solar.nc"
        name = new_file_name (fname)
        inquire (file=trim(name), exist=exists)
        if (.not. exists) then
          print*, "==> Warning: ", trim(name), " does not exist."
          ln = 3
          allocate ( time(ln) )
          allocate ( data(ln) )
          time(:) = year0
          data(:) = solarconst
        else
          call openfile (name, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          allocate ( data(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          call getvara ('solar_constant', iou, ln, ib, ic, data
     &,     c1e3, c0)
        endif
        tim(:) = time(1)
        dat(:) = data(1)
      endif

      tim(2) = min(time(ln), max(time(1), solar_yr))

      if (tim(2) .le. time(1)) then
        dat(2) = data(1)
      elseif (tim(2) .ge. time(ln)) then
        dat(2) = data(ln)
      else
        if (tim(2) .gt. tim(3)) then
          do n=2,ln
            if (time(n-1) .le. tim(2) .and. time(n) .ge. tim(2)) then
              tim(1) = time(n-1)
              dat(1) = data(n-1)
              tim(3) = time(n)
              dat(3) = data(n)
            endif
          enddo
        endif
        wt1 = 1.
        if (tim(3) .ne. tim(1)) wt1 = (tim(3)-tim(2))/(tim(3)-tim(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        dat(2) = dat(1)*wt1 + dat(3)*wt3
      endif

      solarconst = dat(2)

      return
      end
