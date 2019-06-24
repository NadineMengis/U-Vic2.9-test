! source file: /net/mare/home1/eby/as/updates/riv.h
!======================== include file "riv.h" =========================

!     parameters for use in the river model

!     maxnb   = maximum number of basins
!     nb      = number of basins
!     nbp     = number of basin points per basin
!     ndp     = number of discharge points per basin
!     ndis    = map of discharge points
!     nrfill  = map of filled river numbers (extrapolated over water)
!     nriv    = map of river (basin) numbers
!     psum    = total discharge for a basin
!     hsum    = total heat discharge for a basin
!     wdar    = discharge weights over discharge area
!     ta_psum = time average total discharge for a basin

      integer maxnb
      parameter (maxnb=200)

      integer nb, nbp, ndp, ndis, nrfill, nriv
      common /river_i/ nb, nbp(maxnb), ndp(maxnb), ndis(imt,jmt)
      common /river_i/ nrfill(imt,jmt), nriv(imt,jmt)

      real psum, hsum, wdar, ta_psum
      common /river_r/ psum(maxnb), hsum(maxnb),wdar(imt,jmt)

      common /river_r/ ta_psum(maxnb)