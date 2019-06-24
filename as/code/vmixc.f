! source file: /net/mare/home1/eby/as/2.9.old3/source/mom/vmixc.F
      subroutine vmixc (joff, js, je, is, ie)

!=======================================================================
!     set viscosity coefficient on bottom face of "u" cells
!     set diffusion coefficient on bottom face of "t" cells

!     input:
!       joff = offset relating "j" in the MW to latitude "jrow"
!       js   = starting row in the MW
!       je   = ending row in the MW
!       is   = starting longitude index in the MW
!       ie   = ending longitude index in the MW
!=======================================================================

      implicit none

      integer i, k, j, ip, kr, jq, js, je, istrt, is, iend, ie, jstrt
      integer jend, jrow, joff, ks

      real zn2, hab, zkappa, edr

      include "size.h"
      include "param.h"
      include "pconst.h"
      include "stdunits.h"
      include "coord.h"
      include "mw.h"
      include "switch.h"
      include "vmixc.h"
      include "isopyc.h"

!-----------------------------------------------------------------------
!     bail out if starting row exceeds ending row
!-----------------------------------------------------------------------

      if (js .gt. je) return

!-----------------------------------------------------------------------
!     limit the longitude and latitude indices
!-----------------------------------------------------------------------

      istrt = max(2,is)
      iend  = min(imt-1,ie)
      jstrt = max(2,js-1)
      jend  = je-1

!-----------------------------------------------------------------------
!     constant vertical mixing coefficients
!-----------------------------------------------------------------------

      do j=jstrt,jend
        jrow = j + joff
        do i=istrt,iend
          do k=1,km
            visc_cbu(i,k,j) = kappa_m
            diff_cbt(i,k,j) = Ahv(k)
          enddo
        enddo
      enddo

!-----------------------------------------------------------------------
!     Add K33 component to vertical diffusion coefficient
!-----------------------------------------------------------------------

      do j=jstrt,jend
        do i=istrt,iend
          do k=1,km
            diff_cbt(i,k,j) = diff_cbt(i,k,j) + K33(i,k,j)
          enddo
        enddo
      enddo

      return
      end
