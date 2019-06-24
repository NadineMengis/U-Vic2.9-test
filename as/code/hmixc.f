! source file: /net/mare/home1/eby/as/2.9.old3/source/mom/hmixc.F
      subroutine hmixc (joff, js, je, is, ie)

!=======================================================================
!     set horizontal mixing coeffs on north and east face of "t" and
!     "u" cells.

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer js ,je ,jrowstart ,joff ,jrowend ,jrow ,jm1 ,jp1 ,k ,ie
      integer is

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "grdvar.h"
      include "hmixc.h"
      include "mw.h"
      include "scalar.h"
      include "switch.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     set all horizontal mixing coefficients
!-----------------------------------------------------------------------

      jrowstart = js + joff
      jrowend   = min(je + joff + 2,jmt)

!     for momentum... set coefficients for all latitudes

      if (first) then
        visc_cnu  = am
        visc_ceu  = am
        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          amc_north(jrow) = visc_cnu*cst(jp1)*dytr(jp1)
     &                              *csur(jrow)*dyur(jrow)
          amc_south(jrow) = visc_cnu*cst(jrow)*dytr(jrow)
     &                              *csur(jrow)*dyur(jrow)
        enddo
      endif

!     for tracers... set coefficients for all latitudes

      if (first) then
        diff_cnt  = ah
        diff_cet  = ah

        do jrow=jrowstart,jrowend
          jm1 = max(1,jrow-1)
          jp1 = min(jmt,jrow+1)
          ahc_north(jrow) = diff_cnt*csu(jrow)*dyur(jrow)*cstr(jrow)
     &                              *dytr(jrow)
          ahc_south(jrow) = diff_cnt*csu(jm1)*dyur(jm1)*cstr(jrow)
     &                              *dytr(jrow)
        enddo
      endif

      return
      end
