! source file: /net/mare/home1/eby/as/2.9.old3/source/common/cpolar.h
!====================== include file "cpolar.h" =========================

!     polar transform coefficients used to transform velocities near
!     poles before filtering

      real spsin, spcos
      common /cpolar_r/ spsin(imt), spcos(imt)
