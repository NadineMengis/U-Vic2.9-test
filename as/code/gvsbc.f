! source file: /net/mare/home1/eby/as/updates/gvsbc.F
      subroutine gvsbc

!=======================================================================
!     calculates albedo and dalton numbers over vegetation
!     may read and interpolate cropland data
!=======================================================================

      implicit none

      integer i, iou, j, n, ln, ib(10), ic(10)

      logical first_time, intrp

      real data_time, wt3, wt1, z0, yrv(3), iyr(3)

      real, allocatable :: time(:)

      save time, ln, yrv, first_time

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "csbc.h"
      include "cembm.h"
      include "atm.h"
      include "veg.h"
      include "levind.h"
      include "tmngr.h"

      first_time = .true.
      intrp = .false.

!-----------------------------------------------------------------------
!     calculate surface coalbedo and Dalton number for land
!-----------------------------------------------------------------------

      if (intrp .or. first_time) call recalc_land_albedo

      return
      end
